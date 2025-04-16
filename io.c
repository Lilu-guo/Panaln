#ifndef BWTSNP_IO_H
#define BWTSNP_IO_H
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "io.h"
void fasta_error(const char* fastaFname) {
	printf("Error: File %s does not comply with the FASTA file format \n", fastaFname);
	exit(1);
}
void fastq_error(const char* fastqFname) {
	printf("Error: File %s does not comply with the FASTQ file format \n", fastqFname);
	exit(1);
}
void fasta2pac(const char *fastaFname, const char* pacFname, const char* annFname) {
	FILE * fastaFile = (FILE*) fopen(fastaFname, "r");
	if (fastaFile == NULL) {
		printf("fasta2pac: Cannot open FASTA file: %s!\n", fastaFname);
		exit(1);
	}
	FILE * pacFile = (FILE*) fopen(pacFname, "wb");
	if (pacFile == NULL) {
		printf("fasta2pac: Cannot open PAC file: %s!\n", pacFname);
		exit(1);
	}
	FILE * annFile = (FILE*) fopen(annFname, "wb");
	if (annFile == NULL) {
		printf("fasta2pac: Cannot open ANN file: %s!\n", annFname);
		exit(1);
	}
	unsigned char seqBuffer[SEQ_BATCH_ALLOC_LEN];
	unsigned char packedSeqBuffer[SEQ_BATCH_ALLOC_LEN/CHARS_PER_BYTE];
	bwtint_t allocatedSeqLen = SEQ_BATCH_ALLOC_LEN;
	char* seq = (char*) malloc(allocatedSeqLen * sizeof(char));
	bwtint_t totalSeqLen = 0;
	bwtint_t seqBufferCount = 0;
	fasta_annotations_t* annotations = (fasta_annotations_t*) calloc(1, sizeof(fasta_annotations_t));
	bwtint_t allocatedAnnNum = ANN_ALLOC_LEN;
	annotations->seq_anns = (seq_annotation_t*) malloc(allocatedAnnNum * sizeof(seq_annotation_t));
	char c = (char) getc(fastaFile);
	if(c != '>') fasta_error(fastaFname);
	while(!feof(fastaFile)) {
		if(allocatedAnnNum == annotations->num_seq) {
			allocatedAnnNum <<= 1;
			annotations->seq_anns = (seq_annotation_t*)realloc(annotations->seq_anns, allocatedAnnNum * sizeof(seq_annotation_t));
		}
		seq_annotation_t* seqAnnotation = &(annotations->seq_anns[annotations->num_seq]);
		c = (char) getc(fastaFile);
		bwtint_t seqNameLen = 0;
		while(c != '\n' && seqNameLen < MAX_SEQ_NAME_LEN && !feof(fastaFile)){
			seqAnnotation->name[seqNameLen] = c;
			seqNameLen++;
			c = (char) getc(fastaFile);
		}
		seqAnnotation->name[seqNameLen]='\0';
		while(c != '\n' && !feof(fastaFile)){
			c = (char) getc(fastaFile);
		}
		if(feof(fastaFile)) fasta_error(fastaFname);
		bwtint_t seqLen = 0;
		while(c != '>' && !feof(fastaFile)){
			if (c != '\n'){
				if (c >= 'a' && c <= 'z'){
					c += 'A'-'a';
				}
				if (seqLen >= allocatedSeqLen) {
					allocatedSeqLen <<= 1;
					seq = (char*)realloc(seq, sizeof(char)*allocatedSeqLen);
				}
				*(seq + seqLen) = c;
				seqLen++;
			}
			c = (char) getc(fastaFile);
		}
		*(seq + seqLen) = '$';
		seqLen++;
		printf("Done reading a sequence of size %" PRIbwtint_t " from FASTA\n", seqLen);
		for(bwtint_t i = 0; i < seqLen; i++) {
			if (seqBufferCount >= SEQ_BATCH_ALLOC_LEN) {
				pack_byte(seqBuffer, packedSeqBuffer, SEQ_BATCH_ALLOC_LEN);
				fwrite(packedSeqBuffer, 1, SEQ_BATCH_ALLOC_LEN / CHARS_PER_BYTE, pacFile);
				seqBufferCount = 0;
			}
			seqBuffer[seqBufferCount] = nt16_table[(unsigned int) seq[i]];
			seqBufferCount++;
		}
		seqAnnotation->start_index = totalSeqLen;
		seqAnnotation->end_index = totalSeqLen + seqLen - 1;
		totalSeqLen += seqLen;
		annotations->num_seq++;
	}
	if (seqBufferCount > 0) {
		pack_byte(seqBuffer, packedSeqBuffer, seqBufferCount);
		fwrite(packedSeqBuffer, 1, ceil((float)seqBufferCount/ CHARS_PER_BYTE), pacFile);
		seqBufferCount = 0;
	}
	c = (char)(totalSeqLen % CHARS_PER_BYTE);
	fwrite(&c, 1, 1, pacFile);
	fprintf(annFile, "%" PRIbwtint_t "\t%d\n", totalSeqLen, annotations->num_seq);
	for(int i = 0; i < annotations->num_seq; i++) {
		seq_annotation_t seqAnnotation = annotations->seq_anns[i];
		fprintf(annFile, "%s\t%" PRIbwtint_t "\t%" PRIbwtint_t "\n", seqAnnotation.name, seqAnnotation.start_index, seqAnnotation.end_index);
	}
	printf("Done reading FASTA file. Total sequence length read = %" PRIbwtint_t "\n", totalSeqLen);
	fclose(fastaFile);
	fclose(pacFile);
	fclose(annFile);
	free(annotations->seq_anns);
	free(annotations);
	free(seq);
}
void ref2seq(const char* refFname, unsigned char** seq, bwtint_t* totalSeqLen) {
	FILE* refFile = (FILE*) fopen(refFname, "r");
	if (refFile == NULL) {
		printf("ref2seq: Cannot open .ref file: %s!\n", refFname);
		exit(1);
	}
	bwtint_t allocatedSeqLen = 50*1024*1024;
	*seq = (unsigned char*) malloc(allocatedSeqLen * sizeof(unsigned char));
	bwtint_t seqLen = 0;
	unsigned char c = (unsigned char) getc(refFile);
	while(!feof(refFile)) {
		if (seqLen >= allocatedSeqLen) {
			allocatedSeqLen <<= 1;
			*seq = (unsigned char*)realloc(*seq, sizeof(unsigned char)*allocatedSeqLen);
			if(*seq == NULL) {
				printf("ref2seq: Could not allocate memory for the input sequence, size alloc'd = %" PRIbwtint_t "\n", allocatedSeqLen);
				exit(1);
			}
		}
		(*seq)[seqLen] = c;
		seqLen++;
		c = (unsigned char) getc(refFile);
	}
	printf("Done reading FASTA file. Total sequence length read = %" PRIbwtint_t "\n", seqLen);
	*totalSeqLen = seqLen;
}
void fasta2ref(const char *fastaFname, const char* refFname, const char* annFname, unsigned char** seq, bwtint_t *totalSeqLen) { 
	int myseed = 11;
	srand48(myseed);                                                                                        
	FILE* fastaFile = (FILE*) fopen(fastaFname, "r");                                                       
	if (fastaFile == NULL) {
		printf("fasta2ref: Cannot open FASTA file: %s!\n", fastaFname);
		exit(1);
	}
	FILE * annFile = (FILE*) fopen(annFname, "wb");                                                         
	if (annFile == NULL) {
		printf("fasta2ref: Cannot open .ann file: %s!\n", annFname);
		exit(1);
	}
	fasta_annotations_t* annotations = (fasta_annotations_t*) calloc(1, sizeof(fasta_annotations_t));
	bwtint_t allocatedAnnNum = ANN_ALLOC_LEN;
	annotations->seq_anns = (seq_annotation_t*) malloc(allocatedAnnNum * sizeof(seq_annotation_t));
	FILE* refFile;
	int save_ref = 0;
	if(refFname != NULL) {
		save_ref = 1;                                                                                      
		refFile = (FILE*) fopen(refFname, "wt");                                                           
		if (refFile == NULL) {
			printf("fasta2ref: Cannot open .ref file: %s!\n", refFname);
			exit(1);
		}
	}
	bwtint_t allocatedSeqLen = 50*1024*1024;
	*seq = (unsigned char*) malloc(allocatedSeqLen * sizeof(unsigned char));
	bwtint_t seqLen = 0;                                                                         
	char c = (char) getc(fastaFile);
	if(c != '>') fasta_error(fastaFname);
	while(!feof(fastaFile)) {                                                                    
		if(allocatedAnnNum == annotations->num_seq) {
			allocatedAnnNum <<= 1;
			annotations->seq_anns = (seq_annotation_t*)realloc(annotations->seq_anns, allocatedAnnNum * sizeof(seq_annotation_t));  
		}
		seq_annotation_t* seqAnnotation = &(annotations->seq_anns[annotations->num_seq]);
		c = (char) getc(fastaFile);
		bwtint_t seqNameLen = 0;
		while(c != '\n' && seqNameLen < MAX_SEQ_NAME_LEN && !feof(fastaFile)){
			if(c==' ') break;                                                              
			seqAnnotation->name[seqNameLen] = c;
			seqNameLen++;
			c = (char) getc(fastaFile);
		}
		seqAnnotation->name[seqNameLen]='\0';
		while(c != '\n' && !feof(fastaFile)){ 
			c = (char) getc(fastaFile); 
		}
		if(feof(fastaFile)) fasta_error(fastaFname);
		bwtint_t subseqLen = 0;                                                            
		while(c != '>' && !feof(fastaFile)){                                               
			if (c != '\n'){
				if (c >= 'a' && c <= 'z'){                                                 
					c += 'A'-'a';
				}
				if (seqLen >= allocatedSeqLen) {
					allocatedSeqLen <<= 1;                                                 
					*seq = (unsigned char*)realloc(*seq, sizeof(unsigned char)*allocatedSeqLen);
					if(*seq == NULL) {
						printf("fasta2ref: Could not allocate memory for the input sequence, size alloc'd = %" PRIbwtint_t "\n", allocatedSeqLen);
						exit(1);
					}
				}
				if(c == 'N') c = iupacChar[nt4_gray[lrand48()&3]];                        
				(*seq)[seqLen] = nt16_table[(unsigned int) c];                            
				seqLen++;
				subseqLen++;
				if(save_ref) {
					putc((int) (*seq)[seqLen-1], refFile);                                
				}
			}
			c = (char) getc(fastaFile); 
		}
		(*seq)[seqLen] = nt16_table[(unsigned int)'$'];          
		seqLen++;                                                
		subseqLen++;                                             
		if(save_ref) {
			putc((*seq)[seqLen-1], refFile);
		}
		seqAnnotation->start_index = seqLen - subseqLen;         
		seqAnnotation->end_index = seqLen - 1;
		annotations->num_seq++;                                  
	}
	printf("Done reading FASTA file. Total sequence length read = %" PRIbwtint_t "\n", seqLen);
	fprintf(annFile, "%" PRIbwtint_t "\t%d\n", seqLen, annotations->num_seq);
	for(int i = 0; i < annotations->num_seq; i++) {
		seq_annotation_t seqAnnotation = annotations->seq_anns[i];
		fprintf(annFile, "%s\t%" PRIbwtint_t "\t%" PRIbwtint_t "\n", seqAnnotation.name, seqAnnotation.start_index, seqAnnotation.end_index);
	}
	*seq = (unsigned char*)realloc(*seq, sizeof(unsigned char)*2*seqLen);
	if(*seq == NULL) {
		printf("fasta2ref: Could not allocate memory for the reference sequence (including complement), of length = %" PRIbwtint_t "\n", 2*seqLen);
		exit(1);
	}
	for(bwtint_t i = 0; i < seqLen; i++) {                       
		(*seq)[2*seqLen-i-1] = iupacCompl[(int) (*seq)[i]];      
	}
	(*totalSeqLen) = 2*seqLen;
	if(save_ref) {
		for(bwtint_t i = 0; i < seqLen; i++) {
			putc((*seq)[seqLen+i], refFile);
		}
		printf("Wrote %" PRIbwtint_t " chars to ref file\n", 2*seqLen);
	}
	fclose(fastaFile);
	fclose(annFile);
	if(save_ref) fclose(refFile);
	free(annotations->seq_anns);
	free(annotations);
}
bubble_t* bubf2bub(const char *bub){
	FILE *bubFile = (FILE*) fopen(bub, "r");
	char buf[1024];                                
	int bufLen=0;
	int z=0;
	int nline=0; 
	while(fgets(buf, sizeof(buf), bubFile)!=NULL){
		nline++;
	}
	fclose(bubFile);
	bubble_t* bubble =(bubble_t*) calloc(nline/2, sizeof(bubble_t)); 
	bubFile = (FILE*) fopen(bub, "r"); 
	while(fgets(buf, sizeof(buf), bubFile)!=NULL){ 
		bubble[z].ann =atoi(buf); 
		fgets(buf, sizeof(buf), bubFile);          
		bufLen =strlen(buf);
		uint32_t AB[6];
		char numS[12];                             
		int j=0, k=0;
		for(int i=0; i<bufLen; i++){
			char tmp =buf[i];
			if(tmp=='\t' || i==bufLen-1){          
				uint32_t num =atol(numS);
				AB[k++] =num;
				memset(numS, '\0', 12);            
				j =0;                              
			}else{
				numS[j++] =tmp;
			}
		}
		bubble[z].A =AB[0];
		bubble[z].B_minus_A =AB[1];
		bubble[z].C =AB[2];
		bubble[z].D_minus_C =AB[3];
		bubble[z].ref_len =AB[4];
		bubble[z].alt_len =AB[5];
		z++;
	}
	fclose(bubFile);
	return bubble;
}
fasta_annotations_t* annf2ann(const char *annFname) { 
	FILE * annFile = (FILE*) fopen(annFname, "r");
	if (annFile == NULL) {
		printf("annf2ann: Cannot open ANN file: %s!\n", annFname);
		perror(annFname);
		exit(1);
	}
	fasta_annotations_t* annotations = (fasta_annotations_t*) calloc(1, sizeof(fasta_annotations_t));
	bwtint_t totalSeqLen;
	if(fscanf(annFile, "%" SCNbwtint_t "\t%d\n", &totalSeqLen, &(annotations->num_seq)) != 2) { 
		printf("annf2ann: Could not parse ANN file: %s!\n", annFname);
                exit(1);
	}
	annotations->seq_anns = (seq_annotation_t*) malloc(annotations->num_seq * sizeof(seq_annotation_t));
	for(int i = 0; i < annotations->num_seq; i++) {
		seq_annotation_t* seqAnnotation = &(annotations->seq_anns[i]);
		if(fscanf(annFile, "%[^\n\t]\t%" SCNbwtint_t "\t%" SCNbwtint_t "\n",(char*) &(seqAnnotation->name), &(seqAnnotation->start_index), &(seqAnnotation->end_index)) < 3) {
			printf("annf2ann: Could not parse ANN file: %s!\n", annFname);
                	exit(1);
		}
	}
	fclose(annFile);
	return annotations;
}
fasta_annotations_t* annf2Ann(const char *annFname, uint32_t *nchrom, uint64_t *pchrom) { 
	FILE * annFile = (FILE*) fopen(annFname, "r");
	if (annFile == NULL) {
		printf("annf2ann: Cannot open ANN file: %s!\n", annFname);
		perror(annFname);
		exit(1);
	}
	fasta_annotations_t* annotations = (fasta_annotations_t*) calloc(1, sizeof(fasta_annotations_t));
	bwtint_t totalSeqLen;
	if(fscanf(annFile, "%" SCNbwtint_t "\t%d\n", &totalSeqLen, &(annotations->num_seq)) != 2) { 
		printf("annf2ann: Could not parse ANN file: %s!\n", annFname);
                exit(1);
	}
	annotations->seq_anns = (seq_annotation_t*) malloc(annotations->num_seq * sizeof(seq_annotation_t));
	for(int i = 0; i < annotations->num_seq; i++) {
		seq_annotation_t* seqAnnotation = &(annotations->seq_anns[i]);
		if(fscanf(annFile, "%[^\n\t]\t%" SCNbwtint_t "\t%" SCNbwtint_t "\n",(char*) &(seqAnnotation->name), &(seqAnnotation->start_index), &(seqAnnotation->end_index)) < 3) {
			printf("annf2ann: Could not parse ANN file: %s!\n", annFname);
                	exit(1);
		}
		if(!(seqAnnotation->name[0]=='b')) { 
			(*pchrom) = seqAnnotation->end_index; 
			(*nchrom)++; 
		}
	}
	fclose(annFile);
	return annotations;
}
void free_ann(fasta_annotations_t* annotations) {
	free(annotations->seq_anns);
	free(annotations);
}
void pac2seq(const char *pacFname, unsigned char** seq, bwtint_t *totalSeqLen) {
	FILE* pacFile = (FILE*) fopen(pacFname, "rb");
	if (pacFile == NULL) {
		printf("pac2seq: Cannot open PAC file %s!\n", pacFname);
		exit(1);
	}
	fseek(pacFile, -1, SEEK_END);
	bwtint_t pacFileLen = ftell(pacFile);
	if (pacFileLen < 0) {
		printf("pac2seq: Cannot determine the length of the PAC file!\n");
		exit(1);
	}
	unsigned char endByte;
	if(fread(&endByte, sizeof(unsigned char), 1, pacFile) < 1) {
		printf("pac2seq: Cannot read the PAC file!\n");
                exit(1);
	}
	bwtint_t seqLength = pacFileLen*CHARS_PER_BYTE - endByte;
	fseek(pacFile, 0, SEEK_SET);
	unsigned char* packedSeq = (unsigned char*) malloc(sizeof(char) * pacFileLen);
	if(fread(packedSeq, 1, pacFileLen, pacFile) < pacFileLen) {
		printf("pac2seq: Could not read the expected length of the PAC file!\n");
                exit(1);
	}
	fclose(pacFile);
	*seq = (unsigned char*) malloc(sizeof(char) * (2*seqLength + 1));
	unpack_byte(packedSeq, *seq, seqLength);
	for(bwtint_t i = 0; i < seqLength; i++) {
		(*seq)[2*seqLength-i-1] = iupacCompl[(int) (*seq)[i]];
	}
	(*totalSeqLen) = 2*seqLength;
	free(packedSeq);
}
void seq2rev_compl(unsigned char* seq, bwtint_t seqLen, unsigned char** rcSeq) {
	*rcSeq = (unsigned char*) malloc(sizeof(char) * (seqLen + 1));
	for(bwtint_t i = 0; i < seqLen; i++) {
		(*rcSeq)[seqLen-i-1] = iupacCompl[seq[i]];
	}
}
reads_t* fastq2reads(const char *readsFname) {
	FILE *readsFile = (FILE*) fopen(readsFname, "r");
	if (readsFile == NULL) {
		printf("load_reads_fastq: Cannot open reads file: %s !\n", readsFname);
		exit(1);
	}
	reads_t *reads = (reads_t*) calloc(1, sizeof(reads_t));                                                   
	int allocatedReads = NUM_READS_ALLOC;
	reads->reads = (read_t*) malloc(allocatedReads*sizeof(read_t));
	reads->count = 0;
	char c;
	while(!feof(readsFile)) {                                                                                
		if (reads->count >= allocatedReads) {
			allocatedReads += NUM_READS_ALLOC;
			reads->reads = (read_t*) realloc(reads->reads, allocatedReads*sizeof(read_t));
		}
		read_t* read = &(reads->reads[reads->count]);                                                  
		c = (char) getc(readsFile);
		while(c != '@' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) break;
		int seqNameLen = 0;
		c = (char) getc(readsFile); 
		while(c != '\n' && seqNameLen < MAX_SEQ_NAME_LEN && !feof(readsFile)){                     
			read->name[seqNameLen] = c;
			seqNameLen++;
			c = (char) getc(readsFile);                                                               
		}
		read->name[seqNameLen] = '\0';
		if(feof(readsFile)) fastq_error(readsFname);
		while (c != '\n' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) fastq_error(readsFname);
		int readLen = 0;
		int allocatedReadLen = READ_LENGTH_ALLOC;
		read->seq = (char*) malloc(allocatedReadLen*sizeof(char));
		read->rc = (char*) malloc(allocatedReadLen*sizeof(char));
		read->qual = (char*) malloc(allocatedReadLen*sizeof(char));
		c = (char) getc(readsFile);
		while (c != '\n' && !feof(readsFile)) {
			if (readLen >= allocatedReadLen) {
				allocatedReadLen += READ_LENGTH_ALLOC;
				read->seq = (char*) realloc(read->seq, allocatedReadLen*sizeof(char));                  
				read->rc = (char*) realloc(read->rc, allocatedReadLen*sizeof(char));
				read->qual = (char*) realloc(read->qual, allocatedReadLen*sizeof(char));
			}
			read->seq[readLen] = nt4_table[(int) c];                                                    
			c = (char) getc(readsFile);
			readLen++;
		}
		read->len = readLen;
		if(feof(readsFile)) fastq_error(readsFname);
		while (c != '+' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) fastq_error(readsFname);
		while(c != '\n' && !feof(readsFile)){
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) fastq_error(readsFname);
		int qualLen = 0;
		c = (char) getc(readsFile);
		while(c != '\n' && !feof(readsFile)) {
			if(qualLen <= readLen) {
				read->qual[qualLen] = c;                                                               
			}
			qualLen++;
			c = (char) getc(readsFile);
		}
		if(qualLen != readLen) {
			printf("Error: The number of quality score symbols does not match the length of the read sequence.\n");
			exit(1);
		}
		read->qual[qualLen] = '\0';
		for(int i = 0; i < read->len; i++) {
			read->rc[read->len-1-i] = nt4_complement[(int)read->seq[i]];                               
		}
		if(read->len > reads->max_len) {                                                               
			reads->max_len = read->len;
		}
		reads->count++;
	}
	printf("Loaded %d reads from %s.\n", reads->count, readsFname);
	fclose(readsFile);
	return reads;
}
char *my_strdup(const char *str) {
    int n = strlen(str) + 1;
    char *dup = malloc(n * sizeof(char));
    if(dup)
    {
        strcpy(dup, str);
    }
    return dup;
}
void parse_read_mapping(read_t* read) {
	int _count = 0;
	int idx = 0;
	char c = read->name[idx];
	while (c != '\0') {
		if(c == '_') {
			_count++;
		}
		idx++;
		c = read->name[idx];
	}
	read->num_mref_pos = _count - 3;
	read->mref_pos = (bwtint_t*) malloc(read->num_mref_pos*sizeof(bwtint_t));
	int token_index = 0;
	const char delimiters[] = "_";
	char* read_name = my_strdup(read->name);
	char* token = strtok(read_name, delimiters);
	while (token != NULL) {
		if(token_index == 1) {
			sscanf(token, "%" SCNbwtint_t "", &read->ref_pos_l);
		} else if(token_index == 2) {
			sscanf(token, "%" SCNbwtint_t "", &read->ref_pos_r);
		} else if(token_index == 3) {
			read->strand = (strcmp(token, "nm") == 0) ? 0 : 1;
		} else if(token_index > 3){
            		sscanf(token, "%" SCNbwtint_t "", &read->mref_pos[token_index-4]);
		}
		token = strtok(NULL, delimiters);
		token_index++;
	}
	free(read_name);
}
void free_reads(reads_t* reads) {
	for(int i = 0; i < reads->count; i++) {
		if(reads->reads[i].seq) free(reads->reads[i].seq);
		if(reads->reads[i].rc) free(reads->reads[i].rc);
		if(reads->reads[i].qual) free(reads->reads[i].qual);
		if(reads->reads[i].aln_path) free(reads->reads[i].aln_path);
		if(reads->reads[i].mref_pos) free(reads->reads[i].mref_pos);
	}
	free(reads->reads);
	free(reads);
}
void print_read(read_t* read) {
	printf("READ %s \n", read->name);
	printf("FWD: ");
	for(int i = 0; i < read->len; i++) {
		printf("%c", iupacChar[nt4_gray[(int) read->seq[i]]]);
	} printf("\n");
	printf("RC: ");
	for(int i = 0; i < read->len; i++) {
		printf("%c", iupacChar[nt4_gray[(int) read->rc[i]]]);
	} printf("\n");
}
void pack_word(const unsigned char *input, uint32_t *output, const bwtint_t length) {
	uint32_t c;
	bwtint_t i, j, k;
	for (i = 0; i< length / CHARS_PER_WORD; i++) {
		c = 0;
		j = i * CHARS_PER_WORD;
		for (k = 0; k < CHARS_PER_WORD; k++) {
			c = c | (input[j+k] << (BITS_IN_WORD - (k+1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
	if (i * CHARS_PER_WORD < length) {
		c = 0;
		j = i * CHARS_PER_WORD;
		for (k = 0; j + k < length; k++) {
			c = c | (input[j+k] << (BITS_IN_WORD - (k+1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
}
void unpack_word(const unsigned char *input, unsigned char *output, const bwtint_t length) {
	unsigned char c;
	bwtint_t i, j, k;
	for (i = 0; i < length / CHARS_PER_WORD; i++) {
		c = input[i];
		j = i * CHARS_PER_WORD;
		for (k = 0; k < CHARS_PER_WORD; k++) {
			output[j + k] = c >> (BITS_IN_WORD - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
	if (i * CHARS_PER_WORD < length) {
		c = input[i];
		j = i * CHARS_PER_WORD;
		for (k = 0; j + k < length; k++) {
			output[j + k] = c >> (BITS_IN_WORD - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
}
void pack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length) {
	unsigned char c;
	bwtint_t i, j, k;
	for (i = 0; i < length/CHARS_PER_BYTE; i++) {
		c = 0;
		j = i * CHARS_PER_BYTE;
		for (k = 0; k < CHARS_PER_BYTE; k++) {
			c = c | (unsigned char)(input[j+k] << (BITS_IN_BYTE - (k+1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
	if (i * CHARS_PER_BYTE < length) {
		c = 0;
		j = i * CHARS_PER_BYTE;
		for (k=0; j+k < length; k++) {
			c = c | (unsigned char)(input[j+k] << (BITS_IN_BYTE - (k+1) * BITS_PER_CHAR));
		}
		output[i] = c;
	}
}
void unpack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length) {
	unsigned char c;
	bwtint_t i, j, k;
	for (i = 0; i < length / CHARS_PER_BYTE; i++) {
		c = input[i];
		j = i * CHARS_PER_BYTE;
		for (k = 0; k < CHARS_PER_BYTE; k++) {
			output[j + k] = c >> (BITS_IN_BYTE - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
	if (i * CHARS_PER_BYTE < length) {
		c = input[i];
		j = i * CHARS_PER_BYTE;
		for (k = 0; j + k < length; k++) {
			output[j + k] = c >> (BITS_IN_BYTE - BITS_PER_CHAR);
			c <<= BITS_PER_CHAR;
		}
	}
}
#endif
