#ifndef BWTSNP_ALIGN_H
#define BWTSNP_ALIGN_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "io.h"
#define READ_BATCH_SIZE 0x40000
#define STATE_M 0
#define STATE_I 1
#define STATE_D 2
#define MAX_DIFF 5
#define ALN_PATH_ALLOC 15000               
#define MAX_SNPS 5
#define ALN_NOMATCH 0
#define ALN_UNIQUE  1
#define ALN_REPEAT  2
#define MAPQ_CONFIDENT 10
#define NUM_PRECALC 16777216            
#define PRECALC_INTERVAL_LENGTH 12      
#define NBT 200000 
typedef struct sa_intv_ {
	bwtint_t L; 
	bwtint_t U; 
	struct sa_intv_* next_intv; 
} sa_intv_t;
typedef struct {
	int size;                  
	sa_intv_t* first_intv;     
	sa_intv_t* last_intv;      
} sa_intv_list_t;
typedef struct {
	int mm_score; 
	int gapo_score; 
	int gape_score; 
	int use_precalc;
	int is_multiref;
	int n_threads;
} aln_params_t;
typedef struct {
	bwtint_t L;                                                               
	bwtint_t U;                                                               
	bwtint_t l;                                                               
	bwtint_t u;                                                               
	uint32_t num_mm:8, num_gapo:8, num_gape:8, num_snps:8;                    
	uint32_t score, i, state:8, aln_length;                               
	char aln_path[ALN_PATH_ALLOC];
} aln_entry_t;
typedef struct {
	fasta_annotations_t* annotations;
	bubble_t* bubble;
	uint32_t nchrom; 
	uint64_t pchrom; 
} AUX;
int align_reads(char* fastaFname, char* readsFname, char* alnsFname, aln_params_t* params);
void add_sa_interval(sa_intv_list_t* intv_list, bwtint_t L, bwtint_t U);
void copy_sa_interval_list(sa_intv_list_t* dest, sa_intv_list_t* source); 
void clear_sa_interval_list(sa_intv_list_t* intv_list);
void free_sa_interval_list(sa_intv_list_t* intv_list);
void print_sa_interval_list(sa_intv_list_t* intv_list);
void store_sa_interval_list(sa_intv_list_t* intv_list, FILE* saFile);
void load_sa_interval_list(sa_intv_list_t* intv_list, FILE* saFile);
int read2index(char* read, int readLen);
void set_default_aln_params(aln_params_t* params);
alns_t* init_alignments();
alns_t* sa_intervals2alns(sa_intv_list_t* intv_list, int aln_length);
void free_alignments(alns_t* alns);
void reset_alignments(alns_t* alns);
void add_alignment(aln_entry_t* e, bwtint_t L, bwtint_t U, int score, alns_t* alns, const aln_params_t* params);
void print_alignments(alns_t* alns);
void alns2alnf(alns_t* alns, FILE* alnFile);
void alns2alnf_bin(alns_t* alns, FILE* alnFile);
alns_t* alnsf2alns(int* num_alns, char *alnFname);
alns_t* alnsf2alns_bin(int* n_alns, char *alnFname);
void eval_alns(char *fastaFname, char *readsFname, char *alnFname, int is_multiref, int max_diff);
void alns2sam(char *fastaFname, char *readsFname, char *alnsFname, char* samFname, int is_multiref, int max_diff);
void alns2sam_fm(char *fastaFname, char *readsFname, char *alnsFname, char* samFname, int is_multiref, int max_diff);
uint8_t code[256];
#endif
