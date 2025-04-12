#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bwt.h"
#include "align.h"
#include <ctype.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>  
#include <unistd.h>
#include "data_prep.h"
#include "comb.h"
static int usage() {
	printf("Usage:   panaln [combine /index /align] <options> \n");
	fprintf(stderr, "--------------\n");
	printf("Function: combine  generate pangenome file\n");
	printf("          index    construct a pangenome index\n");
	printf("          align    perform read alignment\n");
	fprintf(stderr, "--------------\n");
	fprintf(stderr, "Feedback Email: <guolilu@stu.xidian.edu.cn>\n");
	printf("\n");
	return 0;
}
static int combine_usage() {
	printf("Usage: panaln combine -s <input.fasta> -v <input.vcf> -b <basename> \n");
	fprintf(stderr, "--------------\n");
    fprintf(stderr, "Specific:  -s STRING [required] reference genome\n");
	fprintf(stderr, "           -v STRING [required] vcf file\n");
	fprintf(stderr, "           -b STRING [required] basename\n");
	fprintf(stderr, "           -c STRING [optional] context size (default:150)\n");
    fprintf(stderr, "Please use ABSOLUTE PATHs when specifying files.\n");
    fprintf(stderr, "--------------\n");
    fprintf(stderr, "Feedback Email: <guolilu@stu.xidian.edu.cn>\n");
	printf("\n");
	return 1;
}
static int index_usage() {
	printf("Usage: panaln index -p <input.pan> \n");
	fprintf(stderr, "--------------\n");
    fprintf(stderr, "Specific:  -p STRING [required] pangenome file\n");
    fprintf(stderr, "Please use ABSOLUTE PATHs when specifying files.\n");
    fprintf(stderr, "--------------\n");
    fprintf(stderr, "Feedback Email: <guolilu@stu.xidian.edu.cn>\n");
	printf("\n");
	return 1;
}
static int align_usage() {
	printf("Usage: panaln align -x <index_basename> -f <input.fastq> -s <output.sam> \n");
	fprintf(stderr, "--------------\n");
    fprintf(stderr, "Specific:  -x STRING [required] basename\n");
	fprintf(stderr, "           -f STRING [required] fastq file\n");
	fprintf(stderr, "           -s STRING [required] smm file\n");
    fprintf(stderr, "Please use ABSOLUTE PATHs when specifying files.\n");
    fprintf(stderr, "--------------\n");
    fprintf(stderr, "Feedback Email: <guolilu@stu.xidian.edu.cn>\n");
	printf("\n");
	return 1;
}
int main(int argc, char *argv[]) {
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0) {                                          
		char *npan = NULL, * use_esa_file = NULL;  
		int c;
		while ((c = getopt(argc-1, argv+1, "p:")) >= 0) {
			switch (c) {
				case 'p': npan = optarg; break;                           
				case 'h': index_usage(); break;
				default: index_usage(); return 1;
			}
		}
		if(npan == NULL){index_usage(); return 1;}
		index_bwt(npan, use_esa_file);  
	} else if (strcmp(argv[1], "align") == 0) {                                   
		aln_params_t* params = (aln_params_t*) calloc(1, sizeof(aln_params_t));
		set_default_aln_params(params);                                           
		char *nidx = NULL, *nfastq = NULL, *nsam = NULL;  
		int c;
		while ((c = getopt(argc-1, argv+1, "x:f:s:h:")) >= 0) {
			switch (c) {
				case 'x': nidx = optarg; break;     
				case 'f': nfastq = optarg; break;   
				case 's': nsam = optarg; break;     
				case 'h': align_usage(); break;     
				default: align_usage(); return 1;
			}
		}
		if(nidx==NULL || nfastq==NULL || nsam==NULL){align_usage(); return 1;}
		align_reads(nidx, nfastq, nsam, params);
		free(params);
	} else if (strcmp(argv[1], "combine") == 0) {                          
		printf("**** Generate Pangenome **** \n");
		char* path = "dirPan/";
		int par_clear =1; 
		pars_t* pars = (pars_t*) calloc(1, sizeof(pars_t));
		char* nfasta, *nbase, *nvcf;
		set_default_pars(pars);
		int c;
		while ((c = getopt(argc-1, argv+1, "c:s:v:b:h:")) >= 0) 
		{
			switch (c) {
				case 'c': pars->window_size = atoi(optarg);                     
					if(pars->window_size < 0){
						fprintf(stderr, "Context size shouldn't be negative.\n");
						return 1;}
					break;
				case 's': nfasta = optarg; break;
				case 'v': nvcf = optarg; break;
				case 'b': nbase = optarg; break;
				case 'h': combine_usage(); break;
				default: combine_usage(); return 1;
			}
		}
		if(nfasta==NULL || nvcf==NULL || nbase==NULL){combine_usage(); return 1;}
		char **chr_list, **new_chr_list; 
		chr_list = (char **)malloc(MAX_CHROM * sizeof(char *));
		for (int i = 0; i < MAX_CHROM; i++) {
			chr_list[i] = (char *)malloc(COLUMN_LEN * sizeof(char));
		}
		int n=0; 
		new_chr_list = vcf_extract(nvcf, par_clear, chr_list); 
		int row = sizeof(new_chr_list) / sizeof(new_chr_list[0]);
		for(int i=0; i<row; i++){
			if(strlen(new_chr_list[i]) >0){ 
				strcpy(chr_list[n], new_chr_list[i]);
				n++;
			}
		}
		for (int i = 0; i < row; i++) { 
			free(new_chr_list[i]);
		}
		free(new_chr_list); 
		row = sizeof(chr_list) / sizeof(chr_list[0]);
		for(int i=0; i<row; i++){
      		printf("%s ", chr_list[i]); 
		}
		printf("\n");
		char multifasta[strlen(path)+strlen(nbase)+5];
		char bubblefile[strlen(path)+strlen(nbase)+5];
		char datafile[strlen(path)+strlen(nbase)+5];
		sprintf(multifasta, "%s%s.tmp", path, nbase);
		sprintf(bubblefile, "%s%s.pan", path, nbase);
		sprintf(datafile, "%s%s.aux", path, nbase);  
		insert_SNP(nfasta, multifasta, bubblefile, pars); 
		comp_bubble(multifasta, bubblefile, datafile, pars); 
		struct stat info;
		if (stat(multifasta, &info) == 0){ 
			if (unlink(multifasta) == 0){ 
			} 
		}
		for(int i=0; i<n; i++){ 
			char snpfile[strlen(path)+strlen(chr_list[i])+5];
			char indelfile[strlen(path)+strlen(chr_list[i])+7];
			sprintf(snpfile, "%s%s.SNP", path, chr_list[i]);
			sprintf(indelfile, "%s%s.INDEL", path, chr_list[i]);
			if (stat(snpfile, &info) == 0){ 
				unlink(snpfile); 
			}
			if (stat(indelfile, &info) == 0){ 
				unlink(indelfile); 
			}
		}
		printf("The generated file is in: ./dirPan\n");
	} else {
		printf("Error: Unknown command '%s'\n", argv[1]);
		usage();
	}
	return 0;
}
