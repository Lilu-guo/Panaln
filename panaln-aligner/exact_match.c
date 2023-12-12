
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "align.h"
#include "FMapi.h"
#include "exact_match.h"
#include "io.h"

int align_reads_exact_fm(struct FMA *fm, reads_t *reads, sa_intv_list_t *precalc_sa_intervals, const aln_params_t *params, char *alnFname)
{
	sa_intv_list_t **sa_intervals = (sa_intv_list_t **)malloc(reads->count * sizeof(sa_intv_list_t *));
	clock_t t = clock();
	for (int i = 0; i < reads->count; i++)
	{
		if (params->use_precalc)
		{
		}
		else
		{
			exact_match_fm(fm, &reads->reads[i], &sa_intervals[i], params);
		}
	}
	printf("Exact matching time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	t = clock();
	alns_t **alns = (alns_t **)malloc(reads->count * sizeof(alns_t *));
	for (int i = 0; i < reads->count; i++)
	{
		reads->reads[i].alns = sa_intervals2alns(sa_intervals[i], reads->reads[i].len);
		free_sa_interval_list(sa_intervals[i]);
	}
	FILE *alnFile = (FILE *)fopen(alnFname, "a+");
	if (alnFile == NULL)
	{
		printf("align_reads_exact: Cannot open ALN file: %s!\n", alnFname);
		perror(alnFname);
		exit(1);
	}
	for (int i = 0; i < reads->count; i++)
	{
		read_t *read = &reads->reads[i];
		alns2alnf(read->alns, alnFile);
		free_alignments(read->alns);
	}
	fclose(alnFile);
	printf("Storing results time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	free(sa_intervals);
	free(alns);
	return 0;
}

int exact_match_fm(struct FMA *fm, read_t *read, sa_intv_list_t **sa_intervals, const aln_params_t *params)
{
	bwtint_t fm_n;
	Fm_GetN(fm, &fm_n);
	return exact_match_bounded_fm(fm, read->seq, read->len, 0, fm_n, read->len - 1, sa_intervals, params);
}

int exact_match_1to1_bounded_fm(struct FMA *fm, const char *read, const int readLen, const bwtint_t l, const bwtint_t u, const int i, bwtint_t *sa_begin, bwtint_t *sa_end);

int exact_match_bounded_fm(struct FMA *fm, char *read, int readLen, bwtint_t l, bwtint_t u, int i, sa_intv_list_t **sa_intervals, const aln_params_t *params)
{
	sa_intv_list_t *intv_list_curr = (sa_intv_list_t *)calloc(1, sizeof(sa_intv_list_t));
	if (!params->is_multiref)
	{
		bwtint_t L, U;
		int matched = exact_match_1to1_bounded_fm(fm, read, readLen, l, u, i, &L, &U);
		if (matched > 0)
		{
			add_sa_interval(intv_list_curr, L, U);
		}
		*sa_intervals = intv_list_curr;
		return (intv_list_curr->size != 0);
	}
	sa_intv_list_t *intv_list_next = (sa_intv_list_t *)calloc(1, sizeof(sa_intv_list_t));
	add_sa_interval(intv_list_curr, l, u);
	for (int r = i; r >= 0; r--)
	{
		unsigned char c = read[r];
		if (c == nt4_table[(int)'N'])
		{
			clear_sa_interval_list(intv_list_curr);
			break;
		}
		sa_intv_t *intv = intv_list_curr->first_intv;
		for (int s = 0; s < intv_list_curr->size; s++)
		{
			if (params->is_multiref)
			{
				bwtint_t L[BASES_PER_NUCLEOTIDE] = {0};
				bwtint_t U[BASES_PER_NUCLEOTIDE] = {0};
				Occ8p(c, intv->L - 1, intv->U, L, U, fm);
				for (int i = 0; i < BASES_PER_NUCLEOTIDE; i++)
				{
					if (L[i] <= U[i])
					{
						add_sa_interval(intv_list_next, L[i], U[i]);
					}
				}
			}
			else
			{
			}
			intv = intv->next_intv;
		}
		sa_intv_list_t *tmp = intv_list_curr;
		intv_list_curr = intv_list_next;
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);
		if (intv_list_curr->size == 0)
			break;
	}
	*sa_intervals = intv_list_curr;
	free(intv_list_next);
	return (intv_list_curr->size != 0);
}

int exact_match_1to1_bounded_fm(struct FMA *fm, const char *read, const int readLen, const bwtint_t l, const bwtint_t u, const int i, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t L = l;
	bwtint_t U = u;
	for (int j = i; j >= 0; j--)
	{
		if (read[j] > 3)
		{
			return 0;
		}
		const unsigned char c = nt4_gray[(int)read[j]];
		bwtint_t occL, occU;
		if ((L - 1) == U)
		{
			Fm_occ(L - 1, U, iupacChar[c], &occL, &occU, fm);
			occU = occL;
		}
		else
		{
			Fm_occ(L, U, iupacChar[c], &occL, &occU, fm);
		}
		L = occL;
		U = occU;
		if (L > U)
		{
			return 0;
		}
	}
	(*sa_begin) = L;
	(*sa_end) = U;
	return U - L + 1;
}