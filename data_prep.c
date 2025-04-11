#include "data_prep.h"
int has(char** chr_list, char* chr)
{
	int row = sizeof(chr_list) / sizeof(chr_list[0]);
	for(int i=0; i<row; i++){
		if(strcpy(chr_list[i], chr) ==0){
			return 1;
		}
	}
  	return 0; 
}
char** vcf_extract(char* input_filename, int par_clear, char** chr_list)
{
  	int n=0; 
  	char line[COLUMN_LEN * 10]; 
	char chr[COLUMN_LEN], chr_ori[COLUMN_LEN], pos[COLUMN_LEN], ref[COLUMN_LEN], alt[COLUMN_LEN]; 
	char* basename = "dirPan/";
	char **new_chr_list = (char **)malloc(MAX_CHROM * sizeof(char *));
	for (int i = 0; i < MAX_CHROM; i++) {
        new_chr_list[i] = (char *)malloc(COLUMN_LEN * sizeof(char));
    }
	FILE *vcf = (FILE*) fopen(input_filename, "r");
	FILE *SNP=NULL, *INDEL=NULL;
	while(fgets(line, sizeof(line), vcf) !=NULL) 
  	{
		if(line[0] == '#') { 
			continue;
		}
		char* allele_freq = "0";
		int k=0, j=0; 
		char numS[COLUMN_LEN]; 
		int lineLen =strlen(line);
		memset(chr, '\0', strlen(chr));
		memset(pos, '\0', strlen(pos));
		memset(ref, '\0', strlen(ref));
		memset(alt, '\0', strlen(alt)); 
		memset(numS, '\0', strlen(numS));
		for(int i=0; i<lineLen; i++){
			char tmp =line[i];
			if(tmp=='\t' || i==lineLen-1){ 
				k++;
				if(k==1){strncpy(chr, numS, j);}
				else if(k==2){strncpy(pos, numS, j);}
				else if(k==4){strncpy(ref, numS, j);}
				else if(k==5){strncpy(alt, numS, j);}
				else{}
				memset(numS, '\0', strlen(numS)); 
				j =0; 
			}else{
				numS[j++] =tmp;
			}
		}
		if(strcmp(chr, chr_ori) != 0) 
		{
			char SNP_filename[strlen(basename)+strlen(chr)+5]; 
			sprintf(SNP_filename, "%s%s.SNP", basename, chr); 
			char INDEL_filename[strlen(basename)+strlen(chr)+7]; 
			sprintf(INDEL_filename, "%s%s.INDEL", basename, chr);
			struct stat info;
			if (stat(basename, &info) != 0){ 
				if (mkdir(basename, 0777) != 0){ 
					printf("Can't create file: %s\n", basename);
					exit(1); 
				} 
			}
			if(SNP !=NULL){fclose(SNP);}
			if(INDEL !=NULL){fclose(INDEL);}
			if(par_clear && !has(chr_list, chr) && !has(new_chr_list, chr)) 
			{
				SNP = (FILE*) fopen(SNP_filename, "w");
				INDEL = (FILE*) fopen(INDEL_filename, "w");
				strcpy(new_chr_list[n], chr);
				n++;
			}
			else
			{
				SNP = (FILE*) fopen(SNP_filename, "a+"); 
				INDEL = (FILE*) fopen(INDEL_filename, "a+");
			}
			strcpy(chr_ori, chr); 
		}
		if(strlen(ref) == 1 && strlen(alt) == 1){ 
			fprintf(SNP, "%s\t%s\t%s\t%s\n", pos, ref, alt, allele_freq); 
		}
		else if(strlen(ref) != strlen(alt)){
			fprintf(INDEL, "%s\t%s\t%s\t%s\n", pos, ref, alt, allele_freq);
		}
  	}
  	if(vcf !=NULL){fclose(vcf);}
  	if(SNP !=NULL){fclose(SNP);}
  	if(INDEL !=NULL){fclose(INDEL);}
  	return new_chr_list;
}
static int usage()
{
  	fprintf(stderr, "\n");
  	fprintf(stderr, "Usage:   data_prep [option] <input1.vcf> <input2.vcf> ... \n");
  	fprintf(stderr, "Option:  -c  clear all SNP.extract.chrxx.data and INDEL.extract.chrxx.data files before usage\n");
  	fprintf(stderr, "Example: data_prep -c ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf\n");
  	fprintf(stderr, "\n");
  	return 1;
}
int main0(int argc, char* argv[])
{
  	if(argc < 2) return usage();
  	int par_clear;
  	if(!strcmp(argv[1], "-c")) 
  	{
		par_clear = 1; 
  	}
  	else
  	{
  	  	par_clear = 0;
  	}
  	int i;
  	char **chr_list, **new_chr_list; 
  	char* input_filename;
	chr_list = (char **)malloc(MAX_CHROM * sizeof(char *));
	for (int i = 0; i < MAX_CHROM; i++) {
        chr_list[i] = (char *)malloc(COLUMN_LEN * sizeof(char));
    }
	int n=0; 
  	for(i = 1 + par_clear; i < argc; i++) 
  	{
		input_filename = argv[i]; 
		printf("vcf file: %s\n", input_filename);
		new_chr_list = vcf_extract(input_filename, par_clear, chr_list); 
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
  	}
  	return 0;
}
