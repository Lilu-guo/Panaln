#ifndef myPrint
#define jj 
#else
#define jj printf
#endif
#ifndef BWTSNP_INEXACT_MATCH_H
#define BWTSNP_INEXACT_MATCH_H
#define EMBED_PAD 4
#define NUM_STR 1
#define NUM_CHAR 5
#define MAX_ELEN 150
#define RBITS_PER_STRING (MAX_ELEN * NUM_CHAR) 
#define BITPOS(STR_ID, OFFSET, CHAR_ID) (STR_ID * RBITS_PER_STRING +\
        OFFSET * NUM_CHAR + CHAR_ID)
#include <limits.h>          
#include <stdint.h>
#include "semiWFA/utils/commons.h"
#include "semiWFA/wavefront/wavefront_align.h"
typedef struct {
	uint32_t chr;
	uint32_t pos;
	char* cigar;             
	int score;
	uint8_t ifw;             
	int nx;                  
	int xmS;                 
} Aln;                       
typedef struct {
	int num_diff;
} diff_lower_bound_t;
typedef struct {
	int num_entries;
	int max_entries;
	aln_entry_t* entries;
} heap_bucket_t;
typedef struct {
	int best_score;
	int num_buckets;
	int num_entries;
	heap_bucket_t* buckets;
} priority_heap_t;
static const unsigned char num2asc[5] =  {'A', 'C', 'G', 'T', 'N'};  
static const unsigned char pan2asc[5] =  {'A', 'G', 'C', 'T', 'N'};
static const unsigned char RCP2asc[5] =  {'T', 'C', 'G', 'A', 'N'};
void JoinedToTextOff(uint32_t* tidx, uint32_t* toff, uint32_t* off, uint32_t* qlen); 
char* CigarFormat3(char rawCigar[], int clen);
sa_intv_list_t* calculate_d_fm2(sa_intv_list_t* precalc_sa_intervals_table, struct FMA *fm, char* read, int readLen, diff_lower_bound_t* D, aln_params_t* params);
int align_reads_inexact_fm(FILE* samFile, AUX* aux, struct FMA *fm, reads_t* reads, sa_intv_list_t* precalc_sa_intervals_table, aln_params_t* params, char* alnFname); 
int align_reads_inexact(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals, aln_params_t* params, char* alnsFname);
int align_reads_inexact_parallel(bwt_t *BWT, reads_t* reads, sa_intv_list_t* precalc_sa_intervals_table, aln_params_t* params, char* alnFname);
void inexact_match(bwt_t * BWT, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals, const aln_params_t *params,
		diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns);
void inexact_match_fm(AUX* aux, struct FMA *fm, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals, aln_params_t *params,
		diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns);
int inexact_match_fm_mid(struct FMA *fm, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals, const aln_params_t *params,
		diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns); 
int inexact_match_fm_sid(AUX *aux, struct FMA *fm, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals, const aln_params_t *params,
		diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns); 
int seed_and_extend(sa_intv_list_t* precalc_sa_intervals_table, FILE* samFile, int cds, int cde, AUX *aux, struct FMA *fm, read_t* read, const aln_params_t *params); 
int inexact_match_fm_sid2(struct FMA *fm, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals, const aln_params_t *params,
		diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns); 
void inexact_match_fm_raw(struct FMA *fm, char *read, int readLen, priority_heap_t* heap, sa_intv_list_t* precalc_sa_intervals, const aln_params_t *params,
		diff_lower_bound_t* D, diff_lower_bound_t* D_seed, alns_t* alns); 
priority_heap_t *heap_init(const aln_params_t *p);
void heap_free(priority_heap_t *heap);
void heap_reset(priority_heap_t *heap);
void heap_push(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U, const bwtint_t l, const bwtint_t u, const int num_mm, const int num_gapo, const int num_gape,
		const int state, const int is_diff, const int aln_length, const char* aln_path, const aln_params_t *params);
void heap_push_s(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U,  const int num_mm, const int num_gapo, const int num_gape,
		const int state, const int is_diff, const int aln_length, const char* aln_path, const aln_params_t *params);
void heap_pop(priority_heap_t *heap, aln_entry_t* e);
int cmp(const void *a, const void *b);
#endif
