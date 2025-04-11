#include <getopt.h>
#include <sys/stat.h>  
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define MAX_CHR_LENGTH 1000000000 
#define ALPHABET_SIZE 16
#define BASE_SIZE 4
#define BASE_SET_SIZE 8
#define COLUMN_LEN 1024
#define MAX_CHROM 30
static const int gray_code[] =      { 0,   1,   3,   2,   6,   7,   5,   4,  12,  13,  15,  14,  10,  11,   9,   8};
static const char abbr[] = {'$', 'T', 'K', 'G', 'S', 'B', 'Y', 'C', 'M', 'H', 'N', 'V', 'R', 'D', 'W', 'A'};
static const int base_order[] = {8 , 4 , 2 , 1 };
static const char base[] = {'A', 'C', 'G', 'T'};
static const char base_set[BASE_SIZE][BASE_SET_SIZE] = {{'A', 'N', 'M', 'H', 'V', 'R', 'D', 'W'}, 
                                {'C', 'N', 'S', 'B', 'Y', 'M', 'H', 'V'}, 
                                {'G', 'N', 'K', 'S', 'B', 'V', 'R', 'D'}, 
                                {'T', 'N', 'K', 'B', 'Y', 'H', 'D', 'W'}}; 
typedef struct
{
	int window_size;
	int min_occ;
	int max_occ;
	int min_occ_specified;
        int max_occ_specified;
} pars_t;
long long get_min(long long a, long long b);
long long get_max(long long a, long long b);
int inSet(char test, char test_base);
void print_multigenome(char* multifasta_filename, char* bubble_filename, char *chromosome, long long start, char* header, long long flag, long long * total_snp_number, long long * low_end_snp_number, long long * high_end_snp_number, pars_t* pars);
void insert_SNP(char* fasta_filename, char* multifasta_filename, char* bubble_filename, pars_t* pars); 
void print_bubble(int chr, char* schr, char* bubble_filename, char* data_filename, char* chromosome, long long* indel_count, long long start, int *flag, long long * total_indel_number, long long * low_end_indel_number, pars_t* pars);
void comp_bubble(char* multifasta_filename, char* bubble_filename, char* data_filename, pars_t* pars); 
void set_default_pars(pars_t* pars);