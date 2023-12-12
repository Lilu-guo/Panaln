
#ifndef myPrint
#define jj // printf
#else
#define jj printf
#endif

#ifndef BWTSNP_INEXACT_MATCH_H
#define BWTSNP_INEXACT_MATCH_H

#include <stdint.h>

typedef struct
{
	int num_diff;
	int sa_intv_width;
	sa_intv_list_t *sa_list;
} diff_lower_bound_t;

typedef struct
{
	int num_entries;
	int max_entries;
	aln_entry_t *entries;
} heap_bucket_t;

typedef struct
{
	int best_score;
	int num_buckets;
	int num_entries;
	heap_bucket_t *buckets;
} priority_heap_t;

void inexact_match_fm(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals, aln_params_t *params,
					  diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns);
int inexact_match_fm_mid(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals, const aln_params_t *params,
						 diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns);
int inexact_match_fm_sid(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals, const aln_params_t *params,
						 diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns);
int inexact_match_fm_sid2(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals, const aln_params_t *params,
						  diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns);
void inexact_match_fm_raw(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals, const aln_params_t *params,
						  diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns);

priority_heap_t *heap_init(const aln_params_t *p);
void heap_free(priority_heap_t *heap);
void heap_reset(priority_heap_t *heap);
void heap_push(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U, const bwtint_t l, const bwtint_t u, const int num_mm, const int num_gapo, const int num_gape,
			   const int state, const int is_diff, const int aln_length, const char *aln_path, const aln_params_t *params);
void heap_push_s(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U, /* const bwtint_t l, const bwtint_t u, */ const int num_mm, const int num_gapo, const int num_gape,
				 const int state, const int is_diff, const int aln_length, const char *aln_path, const aln_params_t *params);
void heap_pop(priority_heap_t *heap, aln_entry_t *e);

#endif
