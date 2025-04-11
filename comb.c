#include "comb.h"
long long get_min(long long a, long long b)
{
  	if(a > b) return b; else return a;
}
long long get_max(long long a, long long b)
{
  	if(a > b) return a; else return b;
}
int inSet(char test, char test_base) 
{
  	int i, j;
  	for(i = 0; i < BASE_SIZE; i++) 
  	{
		if(test_base == base[i]) 
		{
			for(j = 0; j < BASE_SET_SIZE; j++) 
			{
				if(test == base_set[i][j] || test == base_set[i][j] - 'A' + 'a')
				{
					return 1; 
				}
			}
			return 0;
		}
  	}
} 
void print_multigenome(char* multifasta_filename, char* bubble_filename, char *chromosome, long long start, char* header, long long flag, long long * total_snp_number, long long * low_end_snp_number, long long * high_end_snp_number, pars_t* pars)
{
  	char chr[20]; int j=0;
	memset(chr, '\0', strlen(chr));
	for(int i=1; i<strlen(header); i++){ 
		if(header[i] !=' '){
			chr[j++] =header[i];
		}else{
			break;
		}
	}
  	char extract_filename[strlen(chr)+12];
	sprintf(extract_filename, "dirPan/%s.SNP", chr);
	FILE *ext = (FILE*) fopen(extract_filename, "r");
	FILE *multifasta, *bubble;
  	if(flag == 1)
  	{
		multifasta = (FILE*) fopen(multifasta_filename, "w");
		bubble = (FILE*) fopen(bubble_filename, "w");
  	}
  	else
  	{
		multifasta = (FILE*) fopen(multifasta_filename, "a+");
		bubble = (FILE*) fopen(bubble_filename, "a+");
  	}
	fprintf(multifasta, "%s", header);
	fprintf(bubble, "%s", header);
	long long i;
	if(ext !=NULL)
	{
		long long pos, occ;
  		char ref[20], alt[20], c_pos[20], c_occ[20];
	  	int bit[BASE_SIZE], total;
		char line[1024];
		while(fgets(line, sizeof(line), ext) !=NULL){
			int k=0, j=0;
			char numS[20];
			int lineLen =strlen(line);
			memset(c_occ, '\0', strlen(c_occ));
			memset(c_pos, '\0', strlen(c_pos));
			memset(ref, '\0', strlen(ref));
			memset(alt, '\0', strlen(alt)); 
			memset(numS, '\0', strlen(numS));
			for(int i=0; i<lineLen; i++){
				char tmp =line[i];
				if(tmp=='\t' || i==lineLen-1){ 
					k++;
					if(k==1){strncpy(c_pos, numS, j);}
					else if(k==2){strncpy(ref, numS, j);}
					else if(k==3){strncpy(alt, numS, j);}
					else if(k==4){strncpy(c_occ, numS, j);}
					else{}
					memset(numS, '\0', strlen(numS)); 
					j =0; 
				}else{
					numS[j++] =tmp;
				}
			}
			pos = atoll(c_pos);
			occ = atoll(c_occ);
			pos = pos + 1; 
			if(pars->min_occ_specified && occ < pars->min_occ){ 
				(*low_end_snp_number)++;
				continue;
			}
			if(pars->max_occ_specified && occ > pars->max_occ){
				(*high_end_snp_number)++;
				chromosome[pos] = alt[0];
				continue;
			}
			(*total_snp_number)++;
			for(i = 0; i < BASE_SIZE; i++){
				if(inSet(chromosome[pos], base[i]) || inSet(ref[0], base[i]) || inSet(alt[0], base[i]))
				{
					bit[i] = 1; 
				}
				else
				{
					bit[i] = 0;
				}
			}
			total = 0; 
			for(i = 0; i < BASE_SIZE; i++){
				total += (bit[i] * base_order[i]); 
			} 
			for(i = 0; i < ALPHABET_SIZE; i++){ 
				if(gray_code[i] == total) 
				{
					chromosome[pos] = abbr[i];
				} 
			}
	  	}
		fclose(ext);
	}
  	for(i = 1; i < start; i++){
		fputc(chromosome[i], multifasta);
		fputc(chromosome[i], bubble);
		if(!(i % 60))
		{
			fputc('\n', multifasta);
			fputc('\n', bubble);
		}
  	}
  	if((start - 1) % 60){
		fputc('\n', multifasta);
		fputc('\n', bubble);
  	}
  	fclose(multifasta);
  	fclose(bubble);
}
void insert_SNP(char* fasta_filename, char* multifasta_filename, char* bubble_filename, pars_t* pars) 
{
	FILE *fasta = (FILE*) fopen(fasta_filename, "r");
	char line[1024]; 
  	char header[1024]; 
  	long long i, start, g = 0; 
	long long total_snp_number = 0, low_end_snp_number = 0, high_end_snp_number = 0; 
  	char *chromosome; 
  	chromosome = (char*)malloc(MAX_CHR_LENGTH * sizeof(char)); 
	memset(line, '\0', strlen(line)); 
  	while(fgets(line, sizeof(line), fasta) !=NULL) 
  	{
		if(strlen(line) >= 1 && line[0] == '>') 
		{
			if(g) 
			{
				print_multigenome(multifasta_filename, bubble_filename, chromosome, start, header, g, &total_snp_number, &low_end_snp_number, &high_end_snp_number, pars); 
			}
			g++; 
			start = 1; 
			memset(header, '\0', strlen(header)); 
			strcpy(header, line);
			memset(line, '\0', strlen(line)); 
		}
		else 
		{
			for(i = 0; i < strlen(line)-1; i++) 
			{
				chromosome[start] = line[i], start++; 
			}
			memset(line, '\0', strlen(line));
		}
  	}
  	print_multigenome(multifasta_filename, bubble_filename, chromosome, start, header, g, &total_snp_number, &low_end_snp_number, &high_end_snp_number, pars); 
	printf("total snp number is %lld\n", total_snp_number);
  	fclose(fasta);
}
void print_bubble(int chr, char* schr, char* bubble_filename, char* data_filename, char* chromosome, long long* indel_count, long long start, int *flag, long long *total_indel_number, long long *low_end_indel_number, pars_t* pars)
{
  	long long pos, occ;
  	int i;
  	char indel_ref[COLUMN_LEN], indel_alt[COLUMN_LEN];
	long long A, B_minus_A, C, D_minus_C, ref_len, alt_len;
  	FILE *ext;
	FILE *bubble, *data;                                                                     
  	bubble = (FILE*) fopen(bubble_filename, "a+"); 
	if(!(*flag))
	{
		data = (FILE*) fopen(data_filename, "w");
		(*flag) = 1;
	}
	else
	{
		data = (FILE*) fopen(data_filename, "a+");
	}
	char extract_filename[strlen(schr)+14];
	sprintf(extract_filename, "dirPan/%s.INDEL", schr);     
  	ext = (FILE*) fopen(extract_filename, "r");
	if(ext !=NULL)
	{
		char indel_ref[COLUMN_LEN], indel_alt[COLUMN_LEN], c_pos[COLUMN_LEN], c_occ[COLUMN_LEN];
		char line[COLUMN_LEN * 10];
	  	while(fgets(line, sizeof(line), ext) !=NULL)                               
  		{
			int k=0, j=0;
			char numS[COLUMN_LEN];
			int lineLen =strlen(line);
			memset(c_occ, '\0', strlen(c_occ));
			memset(c_pos, '\0', strlen(c_pos));
			memset(indel_ref, '\0', strlen(indel_ref));
			memset(indel_alt, '\0', strlen(indel_alt)); 
			memset(numS, '\0', strlen(numS));
			for(int i=0; i<lineLen; i++){
				char tmp =line[i];
				if(tmp=='\t' || i==lineLen-1){ 
					k++;
					if(k==1){strncpy(c_pos, numS, j);}
					else if(k==2){strncpy(indel_ref, numS, j);}
					else if(k==3){strncpy(indel_alt, numS, j);}
					else if(k==4){strncpy(c_occ, numS, j);}
					else{}
					memset(numS, '\0', strlen(numS)); 
					j =0; 
				}else{
					numS[j++] =tmp;
				}
			}
			pos = atoll(c_pos);
			occ = atoll(c_occ);
			pos = pos + 1;                                                               
			(*total_indel_number)++;
			fprintf(bubble, ">bubble%lld %s %lld\n", (*indel_count), schr, get_max(pos - pars->window_size, 1));
			A = get_max(pos - pars->window_size, 1);                                                                             
			B_minus_A = get_min(pars->window_size, pos - 1);                                                                     
			C = pos + strlen(indel_ref);                                                                                        
			D_minus_C = get_min(pars->window_size, start - pos - strlen(indel_ref)) - 1;                                        
			if(indel_ref[0] != '.') ref_len = strlen(indel_ref); else ref_len = 0;
			if(indel_alt[0] != '.') alt_len = strlen(indel_alt); else alt_len = 0;
			fprintf(data, "%d\n", chr);
			fprintf(data, "%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n", A, B_minus_A, C, D_minus_C, ref_len, alt_len);
			for(i = get_min(pars->window_size, pos - 1); i > 0; i--)
			{
				fputc(chromosome[pos - i], bubble);
			}
			if(indel_alt[0] != '.')
			{
				fprintf(bubble, "%s", indel_alt);
			}
			for(i = 0; i < get_min(pars->window_size, start - pos - strlen(indel_ref)); i++)
			{
				fputc(chromosome[pos + strlen(indel_ref) + i], bubble);
			}
			fputc('\n', bubble);
			(*indel_count)++;
	  	}
  		fclose(ext);
	}
  	fclose(bubble);
	fclose(data);
}
void comp_bubble(char* multifasta_filename, char* bubble_filename, char* data_filename, pars_t* pars) 
{
	FILE *multifasta = (FILE*) fopen(multifasta_filename, "r");
  	char line[COLUMN_LEN * 10], schr[COLUMN_LEN];
  	int g = 0, flag = 0;
  	long long indel_count = 0, start;
	long long total_indel_number = 0, low_end_indel_number = 0;
  	char *chromosome;
	int xu_chr = 0; 
  	chromosome = (char*)malloc(MAX_CHR_LENGTH * sizeof(char));
	memset(line, '\0', strlen(line)); 
  	while(fgets(line, sizeof(line), multifasta) !=NULL) 
  	{
		if(strlen(line) >= 1 && line[0] == '>') 
		{
			if(!g) 
			{
				g++;
			}
			else 
			{   
				print_bubble(xu_chr-1, schr, bubble_filename, data_filename, chromosome, &indel_count, start, &flag, &total_indel_number, &low_end_indel_number, pars); 
			} 
			start = 1; 
			memset(schr, '\0', strlen(schr)); 
			int j = 0;
			for(int i=1; i<strlen(line); i++){
				if(line[i] !=' '){
					schr[j++] =line[i];
				}else{
					break;
				}
			} 
			xu_chr++; 
		}
		else
		{ 
			for(int i = 0; i < strlen(line)-1; i++) 
			{
				chromosome[start] = line[i], start++; 
			}
		}
		memset(line, '\0', strlen(line));
  	}
	print_bubble(xu_chr-1, schr, bubble_filename, data_filename, chromosome, &indel_count, start, &flag, &total_indel_number, &low_end_indel_number, pars);  
	printf("total indel number is %lld\n", total_indel_number);
  	fclose(multifasta);
}
static int usage()
{
  	fprintf(stderr, "\n");
  	fprintf(stderr, "Usage: comb <input.fasta> <output.fasta> <output_bubble.fasta> <bubble.data>\n");
	fprintf(stderr, "Option:  -w INT  window size [default: 124]\n");
	fprintf(stderr, "         -i INT  minimum occurrence\n");
	fprintf(stderr, "         -a INT  maximum occurrence\n");
  	fprintf(stderr, "\n");
  	return 1;
}
void set_default_pars(pars_t* pars)                                                
{
	pars->window_size = 150; 
	pars->min_occ_specified = 0;
	pars->max_occ_specified = 0;
}
int main1(int argc, char* argv[])
{
	pars_t* pars = (pars_t*) calloc(1, sizeof(pars_t));
	set_default_pars(pars);
	int c;
	while ((c = getopt(argc, argv, "w:i:a:")) >= 0) 
	{
		switch (c) {
			case 'w': 	pars->window_size = atoi(optarg);                     
					if(pars->window_size < 0)
					{
						fprintf(stderr, "window size shouldn't be negative.\n");
						return 1;
					}
					break;
			case 'i':       pars->min_occ_specified = 1;                      
					pars->min_occ = atoi(optarg);
					break;
			case 'a':       pars->max_occ_specified = 1;                      
					pars->max_occ = atoi(optarg);
					break;
			case '?': 	usage(); return 1;
			default: 	return 1;
		}
	}
	char *fasta_filename = argv[optind]; 
  	char *basename = argv[optind+1];
	char multifasta_filename[strlen(basename)+5], bubble_filename[strlen(basename)+5], data_filename[strlen(basename)+5];
	sprintf(multifasta_filename, "%s.tmp", basename);
	sprintf(bubble_filename, "%s.pan", basename);
	sprintf(data_filename, "%s.aux", basename);
  	insert_SNP(fasta_filename, multifasta_filename, bubble_filename, pars); 
  	comp_bubble(multifasta_filename, bubble_filename, data_filename, pars); 
	return 0;
}
