#include <sys/stat.h>  
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define CHROM_POS 0
#define POS_POS 1
#define REF_POS 3
#define ALT_POS 4
#define FILTER_POS 6
#define INFO_POS 7
#define FORMAT_POS 8
#define COLUMN_LEN 1024
#define MAX_CHROM 30
int has(char** chr_list, char* chr);
char** vcf_extract(char* input_filename, int par_clear, char** chr_list);
