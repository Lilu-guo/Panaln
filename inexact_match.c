#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "bwt.h"
#include "FMapi.h"
#include "align.h"
#include "inexact_match.h"
#include "exact_match.h"
#include <omp.h>
long nbt=0; 
void calculate_d(bwt_t* BWT, char* read, const int readLen, diff_lower_bound_t* D, aln_params_t* params);
static inline int aln_score(const int m, const int o, const int e, const aln_params_t *p) {             
	return m*p->mm_score + o*p->gapo_score + e*p->gape_score;
}
int align_reads_inexact_fm(FILE* samFile, AUX* aux, struct FMA *fm, reads_t* reads, sa_intv_list_t* precalc_sa_intervals_table, aln_params_t* params, char* alnFname) { 
	for (int i = 0; i < aux->nchrom; i++) {          
		fprintf(samFile, "@SQ\tSN:%s\tLN:%d\n", aux->annotations->seq_anns[i].name, (int) (aux->annotations->seq_anns[i].end_index - aux->annotations->seq_anns[i].start_index)); 
	}
	fprintf(samFile, "@PG\tID:panaln\tPN:panaln\tVN:0.1\n");                
	int myseed = 25;
	srand48(myseed);                                                                                     
	jj("gll_max_len:%d \n", reads->max_len);
	diff_lower_bound_t* D = (diff_lower_bound_t*) calloc(reads->max_len+1, sizeof(diff_lower_bound_t));  
	diff_lower_bound_t* D2 = (diff_lower_bound_t*) calloc(reads->max_len+1, sizeof(diff_lower_bound_t)); 
	int cds=0, cde=0;
	sa_intv_list_t* pD2;
	sa_intv_list_t* pD;
	int num_processed = 0;
	while(num_processed < reads->count) {
		clock_t t = clock();
		int batch_size = ((reads->count - num_processed) > READ_BATCH_SIZE ) ? READ_BATCH_SIZE : (reads->count - num_processed);
		for(int i = num_processed; i < num_processed + batch_size; i++) {
			cds =cde =0;
			read_t* read = &reads->reads[i];                                                             
			printf("%d ...start %s   %d \n", i, read->name, read->len);                                  
			pD2 =calculate_d_fm2(precalc_sa_intervals_table, fm, read->rc, read->len, D2, params); 
			jj("D2-bw... num_diff: %d \n", D2[read->len-1].num_diff);                                    
			int sz2=D2[read->len-1].num_diff - D2[0].num_diff;
			if(sz2 !=0){
				free_sa_interval_list(pD2);                                                              
				int p2[sz2];                                                                             
				int b2=0, last2=D2[read->len-1].num_diff, curr2=0;
				for(int i=read->len-1; i>=0; i--){                                                       
					curr2 =D2[i].num_diff;                                                               
					if(curr2!=last2){                                                                    
						p2[b2] =read->len-i-1;                                                           
						b2++;
					}
					last2 =curr2;
					jj("%d ",D2[i].num_diff);                                                            
				}
				jj("\n");
				pD =calculate_d_fm2(precalc_sa_intervals_table, fm, read->seq, read->len, D, params);
				free_sa_interval_list(pD);
				jj("D-fw... num_diff: %d \n", D[read->len-1].num_diff);                                  
				int sz=D[read->len-1].num_diff - D[0].num_diff;
				int p[sz];                                          
				int b=0, last=0, curr=0;
				for(int i=0; i<read->len; i++){                                                          
					curr=D[i].num_diff;                                                                  
					if(curr!=last){                                                                      
						p[b]=i;                                                                          
						b++;
					}
					last=curr;
					jj("%d ",D[i].num_diff);                                                             
				}                                                                                 
				jj("\n");
				int C[sz2+sz+2];                        
				C[0]=0; 
				int _i=0, _j=0, _k=1;
				while(_i!=sz2 && _j!=sz){                                                                 
					if(p2[_i]<=p[_j]){
						C[_k]=p2[_i];
						_i++;
					}else{
						C[_k]=p[_j];
						_j++;
					}
					_k++;
				}
				if(_i==sz2){
					for(; _j<sz; _j++){
						C[_k]=p[_j];
						_k++;
					}
				}else{
					for(; _i<sz2; _i++){
						C[_k]=p2[_i];
						_k++;
					}
				}
				C[sz2+sz+1]=read->len-1;
				int len=0, lC=0, gap=0, lgap=0, s=0, e=0;
				for(int k=1; k<sz2+sz+2; k++){
					gap=C[k]-lC;
					if(gap>len){                                                                         
						len=gap;                                                                         
						s=lC;
						e=C[k];                                                                          
					}
					lC=C[k];
					lgap=gap;
				}
				jj("s:%d e:%d \n",s,e);
				cds =s;                                                                                  
				cde =e-1;
				seed_and_extend(precalc_sa_intervals_table, samFile, cds, cde, aux, fm, read, params);
			}else{ 
				cds=0; 
				bwtint_t fm_n;
				Fm_GetN(fm, &fm_n);
				bwtint_t woff;
				int readLen =read->len;
				char* read2=read->rc;
				char* qual2=read->qual;
				uint32_t tidx =0, toff =0, qlen =readLen;
				uint8_t rd[qlen+1], qu[qlen+1];
				fasta_annotations_t* annotations =aux->annotations;
				uint8_t iFW=0;
				sa_intv_t* intv =pD2->first_intv;                                                        
				int sz=0, z=0;
				for(int i=0; i<pD2->size; i++){
					sa_intv_t* tmp =intv;
					sz =sz+(tmp->U-tmp->L+1);                                                            
					intv =intv->next_intv;
				}
				int rr =lrand48()%(pD2->size);
				jj("itvSz: %d   SZ: %d \n", pD2->size, sz);
				intv =pD2->first_intv;
				for(int i=0; i<rr; i++){                                                                 
					intv =intv->next_intv;
				}
				for(int i=0; i<(pD2->size -rr); i++){
					sa_intv_t* tmp =intv;
					if(z>0){break;}
					for(bwtint_t j=tmp->L; j<=tmp->U; j++){
						if(z>0){break;} 
						int r =lrand48()%(tmp->U - tmp->L+1);              
						woff =Fm_Lookup(fm, j+r);
						iFW =0;
						if(woff >(fm_n/2)){                                
							iFW =1;
							woff =fm_n-woff-(readLen-cds);
						}else{
							woff =woff-(cds);                               
						}
						for(int k=0; k<annotations->num_seq; k++){         
							if((woff >=annotations->seq_anns[k].start_index) && (woff <=annotations->seq_anns[k].end_index)){
								tidx =k; break;
							}
						}
						toff =woff-annotations->seq_anns[tidx].start_index+1;       
						if(iFW){
							for(int i=0; i<qlen; i++){
								rd[qlen-i-1] =RCP2asc[read2[i]];      
								qu[i] =qual2[i];
							}
						}else{
							for(int i=0; i<qlen; i++){
								rd[i] =pan2asc[read2[i]];
								qu[qlen-i-1] =qual2[i];                               
							}
						}
						rd[qlen]=0;                                  
						qu[qlen]=0;
						fprintf(samFile, "%s\t%d\t%s\t%u\t", read->name, iFW?0:16, annotations->seq_anns[tidx].name, toff);
						fprintf(samFile, "%d\t", (sz>1)?0:((tidx<22)?37:3));       
						fprintf(samFile, "%dM\t", readLen);          
						fprintf(samFile, "*\t0\t0\t");               
						fprintf(samFile, "%s\t", rd);                
						fprintf(samFile, "%s\t", qu);                
						fprintf(samFile, "AS:i:%d\n", 0);            
						z++;
					}
					intv =intv->next_intv;
				}
				free_sa_interval_list(pD2);
			}
		}
		printf("Processed %d reads. Inexact matching time: %.2f sec. \n", num_processed+batch_size, (float)(clock() - t) / CLOCKS_PER_SEC);
		num_processed += batch_size;
	}
	free(D);
	free(D2);
	fclose(samFile);
	return 0;
}
void calculate_d(bwt_t *BWT, char* read, int readLen, diff_lower_bound_t* D, aln_params_t* params) {      
	int z = 0;
	bwtint_t L = 0;
	bwtint_t U = BWT->length-1;
	if(!params->is_multiref) {                                                                            
		for (int i = readLen-1; i >= 0; i--) {                                                            
			unsigned char c = nt4_gray[(int) read[i]];
			if(c == 10) {                                                                                 
				L = 0;
				U = BWT->length-1;
				z++;
			} else {
				bwtint_t occL, occU;
				if((L-1) == U) {                                                                          
					occL = O(BWT, c, L-1);
					occU = occL;
				} else {
					occL = O(BWT, c, L-1);
					occU = O(BWT, c, U);
				}
				L = BWT->C[c] + occL + 1;
				U = BWT->C[c] + occU;
				if(L > U) {                                                                               
					L = 0;
					U = BWT->length-1;                                                                    
					z++;
				}
			}
			D[readLen-1-i].num_diff = z;                                                                  
		}
		D[readLen].num_diff = ++z;                                                                        
		return;
	}
	sa_intv_list_t* intv_list_curr = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));
	sa_intv_list_t* intv_list_next = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));
	add_sa_interval(intv_list_curr, L, U);
	for (int i = readLen-1; i >= 0; i--) {
		unsigned char c = read[i];
		int num_matches = 0;
		if(c > 3) {
			clear_sa_interval_list(intv_list_curr); 
		} else {
			sa_intv_t* intv = intv_list_curr->first_intv;
			for(int s = 0; s < intv_list_curr->size; s++) {
				for(int b = 0; b < BASES_PER_NUCLEOTIDE; b++) {
					unsigned char base = nucl_bases_table[c][b];                                         
					if(base == 10) continue; 
					bwtint_t L = BWT->C[base] + O(BWT, base, intv->L-1) + 1;
					bwtint_t U = BWT->C[base] + O(BWT, base, intv->U);
					if (L <= U) {
						num_matches += U - L + 1;
						add_sa_interval(intv_list_next, L, U);                                          
					}
				}
				intv = intv->next_intv;
			}
		}
		sa_intv_list_t* tmp = intv_list_curr;
		intv_list_curr = intv_list_next;
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);
		if(intv_list_curr->size == 0) {                                                                  
			add_sa_interval(intv_list_curr, 0, BWT->length-1);
			z++;
			num_matches = U - L + 1;
		}
		D[readLen-1-i].num_diff = z;                                                                     
	}
	D[readLen].num_diff = ++z;
	free_sa_interval_list(intv_list_curr);
	free_sa_interval_list(intv_list_next);
}
void calculate_d_fm(struct FMA *fm, char* read, int readLen, diff_lower_bound_t* D, aln_params_t* params) {
	int z = 0;                                                                                             
	bwtint_t L = 0, U, fm_n;
	Fm_GetN(fm, &fm_n);
	U=fm_n;                                                                                               
	if(!params->is_multiref) {                                                               
	}
	sa_intv_list_t* intv_list_curr = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));                 
	sa_intv_list_t* intv_list_next = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));                 
	add_sa_interval(intv_list_curr, L, U);                                                                
	for (int i = readLen-1; i >= 0; i--) {
		unsigned char c = read[i]; 
		int num_matches = 0;
		if(c > 3) {
			clear_sa_interval_list(intv_list_curr); 
		} else {
			sa_intv_t* intv = intv_list_curr->first_intv;                                                 
			for(int s = 0; s < intv_list_curr->size; s++) {
				bwtint_t L[BASES_PER_NUCLEOTIDE] = { 0 }; 
				bwtint_t U[BASES_PER_NUCLEOTIDE] = { 0 };
				Occ8p(c,intv->L-1,intv->U,L,U,fm);                                                        
				for(int i=0; i<BASES_PER_NUCLEOTIDE; i++){
					if (L[i] <= U[i]) {
						num_matches += U[i] - L[i] + 1;
						add_sa_interval(intv_list_next, L[i], U[i]);                                      
					}
				}
				intv = intv->next_intv;                                                                   
			}
		}
		sa_intv_list_t* tmp = intv_list_curr;
		intv_list_curr = intv_list_next;                                                                  
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);                                                           
		if(intv_list_curr->size == 0) {                                                                   
			add_sa_interval(intv_list_curr, 0, fm_n);
			z++;                                                                                          
			num_matches = U - L + 1;                                                                      
		}
		D[readLen-1-i].num_diff = z;                                                                      
	}
	free_sa_interval_list(intv_list_curr);                                                                
	free_sa_interval_list(intv_list_next);
}
sa_intv_list_t* calculate_d_fm2(sa_intv_list_t* precalc_sa_intervals_table, struct FMA *fm, char* read, int readLen, diff_lower_bound_t* D, aln_params_t* params) {
	int z=0, sPos=readLen-1;                                                                              
	bwtint_t L = 0, U, fm_n;
	Fm_GetN(fm, &fm_n);
	jj("fm_n: %lld \n",fm_n);
	U=fm_n;                                                                                               
	sa_intv_list_t* intv_list_curr = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));                 
	sa_intv_list_t* intv_list_next = (sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));                 
	if(params->use_precalc) {
		int read_index = read2index(read, readLen);                                              
		if(read_index >= 0) {
			copy_sa_interval_list(intv_list_curr, &(precalc_sa_intervals_table[read_index]));
		}
	}
	if(intv_list_curr->size){                                                                             
		for(int i=sPos; i>sPos-PRECALC_INTERVAL_LENGTH; i--){                                             
			D[readLen-1-i].num_diff = z;
		}
		sPos =sPos-PRECALC_INTERVAL_LENGTH;
	}else{
		add_sa_interval(intv_list_curr, L, U);                                                            
	}                                       
	for (int i = sPos; i >= 0; i--) {
		unsigned char c = read[i];                                                                        
		int num_matches = 0;
		if(c > 3) {
			clear_sa_interval_list(intv_list_curr); 
		} else {
			sa_intv_t* intv = intv_list_curr->first_intv;                                                 
			for(int s = 0; s < intv_list_curr->size; s++) {
				bwtint_t L[BASES_PER_NUCLEOTIDE] = { 0 };                                                 
				bwtint_t U[BASES_PER_NUCLEOTIDE] = { 0 };
				Occ8p(c,intv->L-1,intv->U,L,U,fm);                                                        
				for(int i=0; i<BASES_PER_NUCLEOTIDE; i++){
					if (L[i] <= U[i]) {
						num_matches += U[i] - L[i] + 1;
						add_sa_interval(intv_list_next, L[i], U[i]);                                      
					}
				}
				intv = intv->next_intv;                                                                   
			}
		}
		sa_intv_list_t* tmp = intv_list_curr;
		intv_list_curr = intv_list_next;                                                                  
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);                                                           
		if(intv_list_curr->size == 0) {                                                                   
			z++;                                                                                          
			if(params->use_precalc && i>PRECALC_INTERVAL_LENGTH) {                                        
				int read_index = read2index(read, i);                                              
				if(read_index >= 0) {
					copy_sa_interval_list(intv_list_curr, &(precalc_sa_intervals_table[read_index]));
				}
			}
			if(intv_list_curr->size){                                                                     
				for(int j=i; j>i-PRECALC_INTERVAL_LENGTH; j--){
					D[readLen-1-j].num_diff = z;
				}
				i =i-PRECALC_INTERVAL_LENGTH;
			}else{
				add_sa_interval(intv_list_curr, 0, fm_n);                                                 
			}                               
		}
		D[readLen-1-i].num_diff = z;                                                                      
	}
	free_sa_interval_list(intv_list_next);
	return intv_list_curr;                                                                                   
}
void inexact_match(bwt_t * BWT, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals,
		const aln_params_t *params, diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns) {
}
int seed_and_extend(sa_intv_list_t* precalc_sa_intervals_table, FILE* samFile, int cs, int ce, AUX *aux, struct FMA *fm, read_t* rPt, const aln_params_t *params){
	jj("into...seed_and_extend cems-Len: %d\n",ce-cs);
	int RUN=1000; 
	double H=0.01;    
	int pads=15;
	int readLen=rPt->len;
	char* read=rPt->rc;   
	char* qual=rPt->qual;
	uint32_t tidx =0, toff =0, qlen =readLen;                              
	uint8_t  flank[qlen+2*pads], rd[qlen+1], qu[qlen+1];
	if((ce-cs) <(H*readLen)){                                              
		jj("filted. CEMS: %d \n",ce-cs);
		for(int i=0; i<qlen; i++){
			rd[qlen-i-1] =RCP2asc[read[i]];                                
			qu[i] =qual[i];                           
		}
		rd[qlen]='\0';
		qu[qlen]='\0';
		fprintf(samFile, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", rPt->name, 4);
		fprintf(samFile, "%s\t", rd);                
		fprintf(samFile, "%s\t", qu);                
		fprintf(samFile, "AS:i:%d\n", -999999);      
		return;
	}
	uint8_t iFW=0;                                   
	int xmS=0;                                       
	bwtint_t fm_n;
	Fm_GetN(fm, &fm_n);
	jj("################################################################################################################ \n");
	sa_intv_list_t* sa_intervals = NULL;
	if(params->use_precalc) {
		int read_index = read2index(read, ce);       
		if(read_index >= 0) {
			sa_intervals =(sa_intv_list_t*) calloc(1, sizeof(sa_intv_list_t));
			copy_sa_interval_list(sa_intervals, &(precalc_sa_intervals_table[read_index]));
		}
	}
	if(sa_intervals!=NULL && sa_intervals->size){
		exact_match_bounded_fm2(fm, read, (ce-cs-2)-PRECALC_INTERVAL_LENGTH, 0, fm_n, (ce-1)-PRECALC_INTERVAL_LENGTH, &sa_intervals, params); 
	}else{
		exact_match_bounded_fm(fm, read, ce-cs-1, 0, fm_n, ce-1, &sa_intervals, params); 
	}
	jj("guo...1 sa_intervals:%d \n",sa_intervals->size);
	int sz=0, z=0, zz=0; 
	sa_intv_t* intv =sa_intervals->first_intv;
	for(int i=0; i<sa_intervals->size; i++){
		sa_intv_t* tmp =intv;
		sz =sz+(tmp->U-tmp->L+1); 
		intv =intv->next_intv;
	}
	jj("Read-Len: %u   CEMS-LEN: %d   Rd-Off: %d, %d   intvSz: %d   SZ: %d   FIRST-CEMS-SA: %u, %u \n",
			readLen, ce-cs, cs, ce, sa_intervals->size, sz, sa_intervals->first_intv->L, sa_intervals->first_intv->U);
	bwtint_t woff;            
	fasta_annotations_t* annotations =aux->annotations;
	wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
	attributes.distance_metric = gap_affine;
	attributes.affine_penalties.match = 0;
	attributes.affine_penalties.mismatch = 4;
	attributes.affine_penalties.gap_opening = 6;
	attributes.affine_penalties.gap_extension = 2;
	attributes.alignment_form.span = alignment_endsfree;
	attributes.alignment_form.pattern_begin_free = 0;
	attributes.alignment_form.pattern_end_free = 0;
	attributes.alignment_form.text_begin_free = readLen; 
	attributes.alignment_form.text_end_free = readLen;
  	wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
	Aln * aln = (Aln*) calloc(sz, sizeof(Aln));
	int bSC=-999999;         
	int sSC=-999999;         
	int run=0;               
	intv =sa_intervals->first_intv;
	for(int i=0; i<sa_intervals->size; i++){ 
		sa_intv_t* tmp =intv;
		for(bwtint_t j=tmp->L; j<=tmp->U; j++){
			if(bSC==0 || run>RUN){ break;}
			xmS =0;
			woff =Fm_Lookup(fm, j);               
			zz++;
			jj("   i:%d /%d (z:%d CEMS-LEN:%d)   WholePos: %ld \n", zz, sz, z, ce-cs, woff);
			iFW =0;
			if(woff >(fm_n/2)){                   
				iFW =1;
				woff =fm_n-woff-(readLen);     
				jj("cov2fw...woff: %ld \n", woff);
			}
			if(woff > aux->pchrom){
				for(int k=0; k<annotations->num_seq; k++){              
					if((woff >=annotations->seq_anns[k].start_index) && (woff <=annotations->seq_anns[k].end_index)){
						tidx =k; break;                                 
					}
				}
				bubble_t bu =aux->bubble[tidx- aux->nchrom];
				woff =annotations->seq_anns[bu.ann].start_index + bu.A + (woff-annotations->seq_anns[tidx].start_index+1); 
			}
			if(iFW){                   
				woff =woff+cs; 
			}else{
				woff =woff-cs;
			}
			unsigned char* flank = Fm_Extract(fm, woff+1-pads, readLen+2*pads); 
			int strimI=0, etrimI=0, nD=0, nI=0, nX=0, nos=cs, noe=ce, wnD=0, wnI=0;
			for(int i=0; i<readLen+2*pads; i++){ 
				jj("%c",iupacChar[flank[i]]); 
				flank[i]=grayVal[flank[i]];                                          
				if(flank[i]&(flank[i]-1)){ 
					nX=(nX %readLen) +readLen;
				}
			}
			jj("\n");
			if(iFW){
				for(int i=0; i<qlen; i++){
					jj("%c",pan2asc[nt4_complement[read[readLen-i-1]]]); 
					rd[i] =nt4_gray_val[nt4_complement[read[readLen-i-1]]];          
					qu[i] =qual[i];
				}	
			}else{
				for(int i=0; i<qlen; i++){
					jj("%c",pan2asc[read[i]]); 
					rd[i] =nt4_gray_val[read[i]]; 
					qu[readLen-i-1] =qual[i];
				}
			}
			jj("\n");
			rd[qlen]='\0';
			qu[qlen]='\0';
  			wavefront_align(wf_aligner,rd,readLen,flank,readLen+2*pads); 
			int score =wf_aligner->cigar->score;
			jj("SCORE: %d \n", score);
			if(score>=bSC) {
				sSC =bSC;                                     
				bSC =score;
				run =0;	                                      
			}else{
				run++;
				continue;                                     
			}
			cigar_t* const cigar =wf_aligner->cigar;
			char* const operations = cigar->operations;
  			const int begin_offset = cigar->begin_offset;
  			const int end_offset = cigar->end_offset;
			if(operations[begin_offset]=='I'){ 
				for(int i=begin_offset; i<end_offset; ++i) {
					if(operations[i]=='I'){
						++strimI;
					}else{
						break; 
					}
				}
			}
			if(operations[end_offset-1]=='I'){ 
				for(int i=end_offset-1; i>=begin_offset; --i){
					if(operations[i]=='I'){
						++etrimI;
					}else{
						break; 
					}
				}
			}
			int clen =end_offset -begin_offset -strimI -etrimI;                
			char rawCigar[clen];
			for (int i=begin_offset+strimI,j=0; i<end_offset-etrimI; ++i) {
				if(operations[i]=='I'){ 
					rawCigar[j]='D'; 
					wnD++;
				}else if(operations[i]=='D'){ 
					rawCigar[j]='I'; 
					wnI++;
				}else{
					rawCigar[j] =operations[i]; 
					if(rawCigar[j]=='X'){                                      
						nX++;
						if(xmS<qu[j+wnD]){ xmS =qu[j+wnD]; }                   
					}
				}
				jj("%c", rawCigar[j]);                                         
				++j;
			}
			jj("\n");
			char* cigarStr =CigarFormat3(rawCigar, readLen+wnD);               
			jj("CIGAR: %s \n", cigarStr);
			if(iFW){
				noe =readLen-cs;                                             
				nos =readLen-ce;
			}
			for(int i=0; i<noe; i++){
				switch (rawCigar[i]){
					case 'D': nD++; break;             
					case 'I': nI++; break;
				}
			}
			for(int k=0; k<annotations->num_seq; k++){                  
				if((woff >=annotations->seq_anns[k].start_index) && (woff <=annotations->seq_anns[k].end_index)){
					tidx =k; break;                                     
				}
			}
			toff =woff - annotations->seq_anns[tidx].start_index+1;     
			int tidx2 =tidx;
			if(tidx >1){ 
				jj("bubble.tidx: %d   toff: %u \n", tidx, toff);
				bubble_t bu =aux->bubble[tidx-1];   
				tidx =bu.ann;                                           
				if(toff >= 1 && toff <= bu.B_minus_A)                   
				{   
					toff =toff + bu.A - 1;
				}else if(toff >=bu.B_minus_A + bu.alt_len + 1 && toff <=bu.B_minus_A + bu.alt_len + bu.D_minus_C + 1)
				{   
					toff =toff + bu.C - (bu.B_minus_A + bu.alt_len + 1);
				}else
				{   
					toff =bu.B_minus_A + bu.A + (bu.ref_len / 2) - 1;
				}
				jj("tidx...bu2chr: Chr%d    toff: %d \n", tidx, toff);
				jj("### A:%ld  B-A:%ld  C:%ld  D-C:%ld  ref-len:%d  alt-len:%d \n", bu.A, bu.B_minus_A, bu.C, bu.D_minus_C, bu.ref_len, bu.alt_len); 
				jj("### start:%ld  end:%ld \n", annotations->seq_anns[tidx2].start_index, annotations->seq_anns[tidx2].end_index);
				bSC =sSC;        
				score =score-33; 
			}
			jj("CEMS-Chr: %s   Pos: %u   iFW: %d \n", annotations->seq_anns[tidx].name, toff, iFW);
			uint32_t _toff =toff;                                       
			toff =toff-pads;
			jj("nos:%d   noe:%d   nwI:%d   nwD:%d   nX:%d   xmS:%d \n", nos,noe,wnI,wnD,nX,xmS);
			jj("Final-Chr: %s   Pos: %u \n", annotations->seq_anns[tidx].name, _toff+nI-nD);    
			aln[z].chr =tidx;                                                                   
			aln[z].pos =_toff+nI-nD;
			aln[z].cigar =cigarStr;
			aln[z].score =score;
			aln[z].ifw =iFW;
			aln[z].nx =nX;   
			aln[z].xmS =xmS; 
			z++;                                     
		}
		intv =intv->next_intv;                       
	}
	free_sa_interval_list(sa_intervals);
  	wavefront_aligner_delete(wf_aligner);
	int ii=0, done=0, nTop1=0, nTop2=0;
	for(int i=0; i<z; i++){                          
		if(!done && aln[i].score==bSC){              
			ii =i;
			done =1;
		}
		if(aln[i].score==bSC){
			nTop1++;                                 
		}else{
			nTop2++;                                 
		}
	}
	int mapq=0;
	if(z){                                               
		if(aln[ii].chr >21){                             
			mapq =3;
		}else{
			if(nTop1 >1){                                
				mapq =0;
			}else if(nTop2 ==0){                         
				if( (aln[ii].xmS > 62) || (aln[ii].nx >= readLen) ){
					mapq =37;                         
					if(aln[ii].score < -32){
						mapq =23;
					}
				}else{
					mapq =23;                            
				}
			}else{
				if( (aln[ii].xmS > 62) || (aln[ii].nx >= readLen) ){
					mapq =37;                         
					if(aln[ii].score < -32){
						mapq =23;
					}
				}else{
					mapq =23;                            
				}
			}
		}
	}             
	iFW=aln[ii].ifw;                                 
	if(iFW){
		for(int i=0; i<readLen; i++){
			rd[i] =iupacChar[nt4_gray[nt4_complement[read[readLen-i-1]]]]; 
			qu[i] =qual[i];
		}    
	}else{
		for(int i=0; i<readLen; i++){
			rd[i] =iupacChar[nt4_gray[read[i]]]; 
			qu[readLen-i-1] =qual[i];    
		}                           
	}
	rd[qlen]=0;
	qu[qlen]=0;
	if(z) {jj("---bestAln (ii:%d sz:%d z:%d)  Chr: %s  Pos: %u  MapQ: %d  Score: %d  CIGAR: %s \n", ii, sz, z, annotations->seq_anns[aln[ii].chr].name, aln[ii].pos, mapq, aln[ii].score, aln[ii].cigar);}
	if(z) {
		fprintf(samFile, "%s\t%d\t%s\t%u\t", rPt->name, (aln[ii].ifw)?0:16, annotations->seq_anns[aln[ii].chr].name, aln[ii].pos); 
		fprintf(samFile, "%d\t", mapq);              
		fprintf(samFile, "%s\t", aln[ii].cigar);     
		fprintf(samFile, "*\t0\t0\t");               
		fprintf(samFile, "%s\t", rd);                
		fprintf(samFile, "%s\t", qu);                
		fprintf(samFile, "AS:i:%d\t", aln[ii].score);                
		fprintf(samFile, "XM:i:%d\n", aln[ii].nx %readLen);
	}else{
		fprintf(samFile, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", rPt->name, 4);
		fprintf(samFile, "%s\t", rd);                
		fprintf(samFile, "%s\t", qu);                
		fprintf(samFile, "AS:i:%d\n", -999999);                
	}
	jj("sz: %d   z: %d \n", sz, z);
	for(int i=0; i<z; i++){
		free(aln[i].cigar);                                 
	}
	free(aln);
}
inline char* CigarFormat3(char rawCigar[], int clen){
	char* res =(char *)calloc(3000, 1);          
	char last=0;
	int n=0, k=0;
	for(int i=0; i<clen; i++){
		char c=rawCigar[i];
		if(c=='X') {c='M';} 
		if(c!=last && i!=0){
			char str[10];                       
			sprintf(str, "%d", n);               
			for(int j=0; j<strlen(str); j++){        
				res[k++] =str[j];
			}
			res[k++] =last; 
			n=0;                                 
		}
		n++;
		if(i==clen-1 && n!=0){                   
			char str[10];
			sprintf(str, "%d", n);
			for(int j=0; j<strlen(str); j++){        
				res[k++] =str[j];
			}
			res[k++] =last;
		}
		last=c;
	}
	return res; 
}
int CigarFormat2( 
    char* const buffer,
    cigar_t* const cigar,
    const bool print_matches,
	int strimI, int etrimI) {
  if (cigar_is_null(cigar)) {
    buffer[0] = '\0';
    return 0;
  }
  const char* const operations = cigar->operations;
  const int begin_offset = cigar->begin_offset;
  const int end_offset = cigar->end_offset;
  char last_op = operations[begin_offset+strimI]; 
  int last_op_length = 1;
  int i, cursor = 0;
  for (i=begin_offset+1+strimI; i<end_offset-etrimI; ++i) {
    if (operations[i]==last_op) {
      ++last_op_length;
    } else {
      if (print_matches || last_op != 'M') {
        cursor += sprintf(buffer+cursor,"%d%c",last_op_length,last_op);
      }
      last_op = operations[i];
      last_op_length = 1;
    }
  }
  if (print_matches || last_op != 'M') {
    cursor += sprintf(buffer+cursor,"%d%c",last_op_length,last_op);
  }
  buffer[cursor] = '\0';
  return cursor;
}
void JoinedToTextOff(uint32_t* tidx, uint32_t* toff, uint32_t* off, uint32_t* qlen){
	char var[]="/home/lab/gll/hisat2-2.2.1/grch38_snp/genome_snp.1.ht2";
	FILE* varFile = (FILE*) fopen(var, "r");
	if(varFile == NULL){
		printf("BitPairReference: Cannot open var file: %s \n", var);
		exit(1);
	}
	uint64_t bytesRead =0; 
	uint32_t one, idxVer, len, gbwtLen, numNodes, lineRate, linesPerSide, offRate, ftabChars, eftabLen, flags, _nPat, _nFrag;
	fread(&one, 4, 1, varFile); bytesRead +=4;
	fread(&idxVer, 4, 1, varFile); bytesRead +=4; 
	fread(&len, 4, 1, varFile); bytesRead +=4;
	fread(&gbwtLen, 4, 1, varFile); bytesRead +=4;
	fread(&numNodes, 4, 1, varFile); bytesRead +=4;
	fread(&lineRate, 4, 1, varFile); bytesRead +=4;
	fread(&linesPerSide, 4, 1, varFile); bytesRead +=4;
	fread(&offRate, 4, 1, varFile); bytesRead +=4;
	fread(&ftabChars, 4, 1, varFile); bytesRead +=4;
	fread(&eftabLen, 4, 1, varFile); bytesRead +=4;
	fread(&flags, 4, 1, varFile); bytesRead +=4;
	fread(&_nPat, 4, 1, varFile); bytesRead +=4;
	jj("one:%u idxVer:%u len:%u gbwtLen:%u numNodes:%u lineRate:%u linesPerSide:%u offRate:%u ftabChars:%u eftabLen:%u flags:%u _nPat:%u \n", 
		one, idxVer, len, gbwtLen, numNodes, lineRate, linesPerSide, offRate, ftabChars, eftabLen, flags, _nPat);
	uint32_t* plen = (uint32_t*) calloc(_nPat, sizeof(uint32_t));
	for(uint32_t i=0; i<_nPat; i++){
		fread(&(plen[i]), 4, 1, varFile); bytesRead +=4;
	}
	fread(&_nFrag, 4, 1, varFile); bytesRead +=4;
	jj("_nFrag:%u \n", _nFrag); 
	uint32_t* rstarts = (uint32_t*) calloc(_nFrag*3, sizeof(uint32_t));
	for(uint32_t i=0; i<_nFrag*3; i+=3){ 
		fread(&(rstarts[i]), 4, 1, varFile); bytesRead +=4; 
		fread(&(rstarts[i+1]), 4, 1, varFile); bytesRead +=4; 
		fread(&(rstarts[i+2]), 4, 1, varFile); bytesRead +=4; 
	}
	uint32_t _gbwtLen =(gbwtLen == 0 ? len + 1 : gbwtLen);
	uint32_t _gbwtSz =_gbwtLen/2 + 1;                      
	uint32_t _sideSz =(1<<lineRate);
	uint32_t _sideGbwtSz =_sideSz - (4 * 6); 
	uint32_t _numSides =(_gbwtSz+(_sideGbwtSz)-1)/(_sideGbwtSz);
	uint32_t _gbwtTotLen =_numSides * _sideSz;     fseek(varFile, _gbwtTotLen, SEEK_CUR);    
	uint32_t numZOffs;                             fread(&numZOffs, 4, 1, varFile);
	fseek(varFile, numZOffs*4, SEEK_CUR); 
	fseek(varFile, 5*4, SEEK_CUR);        
	uint32_t _ftabLen =(1<<(ftabChars*2))+1;       fseek(varFile, _ftabLen*4, SEEK_CUR);     
	uint32_t _eftabLen =eftabLen;                  fseek(varFile, _eftabLen*4, SEEK_CUR);    
	char chrName[256];
	int i=0, nChr=0;          
	while(1){
		char c ='\0';
		c= getc(varFile);
		if(c==EOF) break;
		if(c=='\0'){          
			break;
		}else if(c=='\n'){    
			chrName[i] ='\0'; 
			i =0;
			nChr++;
		}else{                
			chrName[i] =c;
			i++;
		}
	}
	uint32_t top =0; 
	uint32_t bot =_nFrag; 
	uint32_t elt =0xffffffff;
	while(1){
		elt =top+((bot-top)>>1); 
		uint32_t lower =rstarts[elt*3]; 
		uint32_t upper; 
		if(elt==_nFrag-1){
			upper =len; 
		}else{
			upper =rstarts[(elt+1)*3]; 
		}
		uint32_t fraglen =upper-lower; 
		if(lower <= (*off)){
			if(upper > (*off)){
				if((*off)+(*qlen) >upper){ 
					printf("exit: Straddled! \n");
					return;
				}
				(*tidx) =rstarts[(elt*3)+1]; 
				uint32_t fragoff =(*off)-rstarts[elt*3]; 
				if(1){
					fragoff =fraglen-fragoff-1;
					fragoff -=((*qlen)-1);
				}
				(*toff) =fragoff+rstarts[(elt*3)+2]; 
				break;
			}else{
				top =elt; 
			}
		}else{
			bot =elt; 
		}
	}
}
