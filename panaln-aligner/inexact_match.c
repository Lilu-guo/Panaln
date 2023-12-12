
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "FMapi.h"
#include "align.h"
#include "inexact_match.h"
#include "exact_match.h"
#include <omp.h>

long nbt = 0;

static inline int aln_score(const int m, const int o, const int e, const aln_params_t *p)
{
	return m * p->mm_score + o * p->gapo_score + e * p->gape_score;
}

int align_reads_inexact_fm(struct FMA *fm, reads_t *reads, sa_intv_list_t *precalc_sa_intervals_table, aln_params_t *params, char *alnFname)
{
	FILE *alnFile = (FILE *)fopen(alnFname, "a+b");
	if (alnFile == NULL)
	{
		printf("align_reads_inexact: Cannot open ALN file: %s!\n", alnFname);
		perror(alnFname);
		exit(1);
	}
	diff_lower_bound_t *D = (diff_lower_bound_t *)calloc(reads->max_len + 1, sizeof(diff_lower_bound_t));
	diff_lower_bound_t *D_seed = (diff_lower_bound_t *)calloc(params->seed_length + 1, sizeof(diff_lower_bound_t));
	diff_lower_bound_t *D2 = (diff_lower_bound_t *)calloc(reads->max_len + 1, sizeof(diff_lower_bound_t));
	diff_lower_bound_t *D2_seed = (diff_lower_bound_t *)calloc(params->seed_length + 1, sizeof(diff_lower_bound_t));
	int imax = 0, max = 0, smax = 0, emax = 0, imax2 = 0, max2 = 0, smax2 = 0, emax2 = 0, cds = 0, cde = 0;
	priority_heap_t *heap = heap_init(params);
	int num_processed = 0;
	while (num_processed < reads->count)
	{
		clock_t t = clock();
		int batch_size = ((reads->count - num_processed) > READ_BATCH_SIZE) ? READ_BATCH_SIZE : (reads->count - num_processed);
		int n_sid2 = 0, n_sid = 0, n_mid = 0, n_raw = 0, n_exact = 0;
		for (int i = num_processed; i < num_processed + batch_size; i++)
		{
			imax = max = smax = emax = imax2 = max2 = smax2 = emax2 = cds = cde = 0;
			read_t *read = &reads->reads[i];
			read->alns = init_alignments();
			sa_intv_list_t *precalc_sa_intervals = NULL;
			if (params->use_precalc)
			{
				int read_index = read2index(read->rc, read->len);
				if (read_index < 0)
				{
					continue;
				}
				precalc_sa_intervals = &(precalc_sa_intervals_table[read_index]);
			}
			calculate_d_fm(fm, read->rc, read->len, D2, params);
			int p2[D2[read->len - 1].num_diff];
			int b2 = 0, last2 = D2[read->len - 1].num_diff, curr2 = 0;
			for (int i = read->len - 1; i >= 0; i--)
			{
				curr2 = D2[i].num_diff;
				if (curr2 != last2)
				{
					p2[b2] = read->len - i - 1;
					b2++;
				}
				last2 = curr2;
			}

			if (D2[read->len - 1].num_diff > params->max_diff)
			{
				continue;
			}

			int dv = D2[read->len - 1].num_diff;
			if (dv != 0)
			{
				calculate_d_fm(fm, read->seq, read->len, D, params);
				int p[D[read->len - 1].num_diff];
				int b = 0, last = 0, curr = 0;
				for (int i = 0; i < read->len; i++)
				{
					curr = D[i].num_diff;
					if (curr != last)
					{
						p[b] = i;
						b++;
					}
					last = curr;
					jj("%d", D[i].num_diff);
				}

				int imax = D2[read->len - 1].num_diff;
				int jmax = D[read->len - 1].num_diff;
				int C[imax + jmax + 2];
				for (int i = 0; i < imax; i++)
				{
					C[i] = p2[i];
				}
				for (int j = 0; j < jmax; j++)
				{
					C[imax + j] = p[j];
				}
				C[imax + jmax] = 0;
				C[imax + jmax + 1] = read->len - 1;
				for (int i = 0, t = 0; i < imax + jmax + 1; i++)
				{
					for (int j = 0; j < imax + jmax - i + 1; j++)
					{
						if (C[j] > C[j + 1])
						{
							t = C[j];
							C[j] = C[j + 1];
							C[j + 1] = t;
						}
					}
				}

				int len = 0, lC = 0, gap = 0, lgap = 0, s = 0, e = 0;
				for (int k = 1; k < imax + jmax + 2; k++)
				{
					gap = C[k] - lC;
					if (gap > len)
					{
						len = gap;
						s = lC;
						e = C[k];
					}
					lC = C[k];
					lgap = gap;
				}
				cds = s;
				cde = e - 1;
			}
			else
			{
				n_exact++;
				cds = 0;
				cde = 0;
			}

			params->no_s = cds;
			params->no_e = cde;

			int wj = 0;
			if (cds == 0)
			{
				n_sid2++;
				wj = inexact_match_fm_sid2(fm, read->seq, read->len, heap, precalc_sa_intervals, params, D2, D2_seed, read->alns);
			}
			else if (cde == read->len - 2)
			{
				n_sid++;
				wj = inexact_match_fm_sid(fm, read->rc, read->len, heap, precalc_sa_intervals, params, D, D_seed, read->alns);
			}
			else
			{
				n_mid++;
				wj = inexact_match_fm_mid(fm, read->rc, read->len, heap, precalc_sa_intervals, params, D, D_seed, read->alns);
			}

			if (read->alns->num_entries == 0 && wj == 0)
			{
				n_raw++;
				calculate_d_fm(fm, read->seq, params->seed_length, D_seed, params);
				inexact_match_fm_raw(fm, read->rc, read->len, heap, precalc_sa_intervals, params, D, D_seed, read->alns);
			}
		}
		printf("read_n: sid2:%d sid:%d mid:%d raw:%d exact:%d nbt:%ld \n", n_sid2, n_sid, n_mid, n_raw, n_exact, nbt);
		printf("Processed %d reads. Inexact matching time: %.2f sec.", num_processed + batch_size, (float)(clock() - t) / CLOCKS_PER_SEC);

		clock_t ts = clock();
		for (int i = num_processed; i < num_processed + batch_size; i++)
		{
			read_t *read = &reads->reads[i];
			alns2alnf(read->alns, alnFile);
			free_alignments(read->alns);
			free(read->seq);
			free(read->rc);
			free(read->qual);
			read->seq = read->rc = read->qual = NULL;
		}
		printf("Storing results time: %.2f sec\n", (float)(clock() - ts) / CLOCKS_PER_SEC);
		num_processed += batch_size;
	}

	free(D);
	free(D_seed);
	heap_free(heap);
	fclose(alnFile);
	return 0;
}

void calculate_d_fm(struct FMA *fm, char *read, int readLen, diff_lower_bound_t *D, aln_params_t *params)
{
	int z = 0;
	bwtint_t L = 0, U, fm_n;
	Fm_GetN(fm, &fm_n);
	U = fm_n;

	sa_intv_list_t *intv_list_curr = (sa_intv_list_t *)calloc(1, sizeof(sa_intv_list_t));
	sa_intv_list_t *intv_list_next = (sa_intv_list_t *)calloc(1, sizeof(sa_intv_list_t));
	add_sa_interval(intv_list_curr, L, U);

	for (int i = readLen - 1; i >= 0; i--)
	{
		unsigned char c = read[i];
		int num_matches = 0;
		if (c > 3)
		{
			clear_sa_interval_list(intv_list_curr);
		}
		else
		{
			sa_intv_t *intv = intv_list_curr->first_intv;
			for (int s = 0; s < intv_list_curr->size; s++)
			{
				bwtint_t L[BASES_PER_NUCLEOTIDE] = {0};
				bwtint_t U[BASES_PER_NUCLEOTIDE] = {0};
				Occ8p(c, intv->L - 1, intv->U, L, U, fm);
				for (int i = 0; i < BASES_PER_NUCLEOTIDE; i++)
				{
					if (L[i] <= U[i])
					{
						num_matches += U[i] - L[i] + 1;
						add_sa_interval(intv_list_next, L[i], U[i]);
					}
				}
				intv = intv->next_intv;
			}
		}
		sa_intv_list_t *tmp = intv_list_curr;
		intv_list_curr = intv_list_next;
		intv_list_next = tmp;
		clear_sa_interval_list(intv_list_next);

		if (intv_list_curr->size == 0)
		{
			add_sa_interval(intv_list_curr, 0, fm_n);
			z++;
			num_matches = U - L + 1;
		}
		D[readLen - 1 - i].num_diff = z;
		D[readLen - 1 - i].sa_intv_width = num_matches;
	}

	D[readLen - 1].sa_list = intv_list_curr;
	D[readLen].sa_intv_width = 0;
	D[readLen].num_diff = ++z;
	// free_sa_interval_list(intv_list_curr);
	free_sa_interval_list(intv_list_next);
}

void inexact_match_fm_raw(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals,
						  const aln_params_t *params, diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns)
{
	bwtint_t fm_n;
	Fm_GetN(fm, &fm_n);

	int countN = 0;
	for (int i = 0; i < readLen; i++)
	{
		if (read[i] > 3)
			countN++;
	}
	if (countN > params->max_diff)
	{
		return;
	}

	heap_reset(heap);
	if (precalc_sa_intervals != NULL)
	{
		if (precalc_sa_intervals->size != 0)
		{
			char aln_path[ALN_PATH_ALLOC] = {0};
			sa_intv_t *intv = precalc_sa_intervals->first_intv;
			for (int i = 0; i < precalc_sa_intervals->size; i++)
			{
				heap_push_s(heap, readLen - PRECALC_INTERVAL_LENGTH, intv->L, intv->U, 0, 0, 0, 0, 0, PRECALC_INTERVAL_LENGTH - 1, aln_path, params);
				intv = intv->next_intv;
			}
		}
		else
		{
			return;
		}
	}
	else
	{
		heap_push_s(heap, readLen, 0, fm_n, 0, 0, 0, 0, 0, 0, 0, params);
	}

	int best_score = aln_score(params->max_diff + 1, params->max_gapo + 1, params->max_gape + 1, params);
	int best_diff = params->max_diff + 1;
	int max_diff = params->max_diff;
	int num_best = 0;
	int max_entries = 0;

	int total_entries = 0;
	int last_num_entries = 1;
	int tid = 0;

	int nbt0 = nbt;
	while (heap->num_entries != 0)
	{
		total_entries += heap->num_entries - last_num_entries + 1;
		last_num_entries = heap->num_entries;
		if (heap->num_entries > max_entries)
			max_entries = heap->num_entries;
		if (heap->num_entries > params->max_entries)
		{
			break;
		}
		aln_entry_t e_;
		heap_pop(heap, &e_);
		aln_entry_t *e = &e_;
		nbt++;

		if (e->score > (best_score + params->mm_score))
		{
			break;
		}
		int diff_left = max_diff - e->num_mm - e->num_gapo - e->num_gape;
		if (diff_left < 0)
		{
			continue;
		}
		if ((e->i > 0) && (diff_left < D[e->i - 1].num_diff))
		{
			continue;
		}
		int diff_left_seed = params->max_diff_seed - e->num_mm - e->num_gapo - e->num_gape;
		int seed_index = e->i - (readLen - params->seed_length);
		if ((seed_index > 0) && (diff_left_seed < D_seed[seed_index - 1].num_diff))
		{
			continue;
		}

		if (e->i == 0)
		{
			int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
			if (alns->num_entries == 0)
			{
				best_score = score;
				best_diff = e->num_mm + e->num_gapo + e->num_gape;
				max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
			}
			if (score == best_score)
			{
				num_best += e->U - e->L + 1;
			}
			else if (num_best > params->max_best)
			{
				break;
			}
			add_alignment(e, e->L, e->U, score, alns, params);
			continue;
		}
		else if (diff_left == 0)
		{
			sa_intv_list_t *sa_intervals;
			if (exact_match_bounded_fm(fm, read, readLen, e->L, e->U, e->i - 1, &sa_intervals, params) > 0)
			{
				int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
				if (alns->num_entries == 0)
				{
					best_score = score;
					best_diff = e->num_mm + e->num_gapo + e->num_gape;
					max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
				}
				if (score == best_score)
				{
					sa_intv_t *intv = sa_intervals->first_intv;
					for (int k = 0; k < sa_intervals->size; k++)
					{
						num_best += intv->U - intv->L + 1;
						intv = intv->next_intv;
					}
				}
				else if (num_best > params->max_best)
				{
					break;
				}

				e->aln_length += e->i;
				sa_intv_t *intv = sa_intervals->first_intv;
				for (int k = 0; k < sa_intervals->size; k++)
				{
					add_alignment(e, intv->L, intv->U, score, alns, params);
					intv = intv->next_intv;
				}
			}

			free_sa_interval_list(sa_intervals);
			continue;
		}

		long L[ALPHABET_SIZE] = {0};
		long U[ALPHABET_SIZE] = {0};
		int alphabet_size = ALPHABET_SIZE;
		int is_multiref = params->is_multiref;
		if (params->is_multiref)
		{
			Occ16p(e->L - 1, e->U, L, U, fm);
			for (int i = 1; i < ALPHABET_SIZE; i++)
			{
				if (U[i] - L[i] + 1 == 0)
				{
					continue;
				}
			}
		}
		else
		{
		}

		int allow_diff = 1;
		int allow_indels = 1;
		int allow_mm = 1;
		int allow_open = 1;
		int allow_extend = 1;

		if (e->i - 1 > 0)
		{
			if ((diff_left - 1) < D[e->i - 2].num_diff)
			{
				allow_diff = 0;
			}
			else if (((D[e->i - 1].num_diff == diff_left - 1) && (D[e->i - 2].num_diff == diff_left - 1)) && (D[e->i - 1].sa_intv_width == D[e->i - 2].sa_intv_width))
			{
				allow_mm = 0;
			}
		}

		if (seed_index - 1 > 0)
		{
			if ((diff_left_seed - 1) < D_seed[seed_index - 2].num_diff)
			{
				allow_diff = 0;
			}
			else if ((D_seed[seed_index - 1].num_diff == diff_left_seed - 1) && (D_seed[seed_index - 2].num_diff == diff_left_seed - 1) && (D_seed[seed_index - 1].sa_intv_width == D_seed[seed_index - 2].sa_intv_width))
			{
				allow_mm = 0;
			}
		}

		int tmp = e->num_gapo + e->num_gape;
		if ((e->i - 1 < (params->no_indel_length + tmp)) || ((readLen - (e->i - 1)) < (params->no_indel_length + tmp)))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo) && (e->num_gape >= params->max_gape))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo))
		{
			allow_open = 0;
		}
		if ((e->num_gape >= params->max_gape))
		{
			allow_extend = 0;
		}

		if (allow_diff && allow_indels)
		{
			if (e->state == STATE_I)
			{
				if (allow_extend)
				{
					heap_push_s(heap, e->i - 1, e->L, e->U, e->num_mm, e->num_gapo, e->num_gape + 1, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
				}
			}
			else
			{
				if (allow_open && (e->state == STATE_M))
				{
					heap_push_s(heap, e->i - 1, e->L, e->U, e->num_mm, e->num_gapo + 1, e->num_gape, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
				}
				for (int j = 1; j < alphabet_size; j++)
				{
					if (L[j] <= U[j])
					{
						if (e->state == STATE_M)
						{
							if (allow_open)
							{
								heap_push_s(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo + 1, e->num_gape, STATE_D,
											e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
						else
						{
							if (allow_extend)
							{
								heap_push_s(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape + 1, STATE_D,
											e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
					}
				}
			}
		}

		int c = read[e->i - 1];
		if (allow_diff && allow_mm)
		{
			for (int j = 1; j < alphabet_size; j++)
			{
				if (L[j] <= U[j])
				{
					int is_mm = 0;
					if (is_multiref)
					{
						if ((c > 3) || (j == 10) /*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0))
						{
							is_mm = 1;
						}
					}
					else
					{
						if ((c > 3) || (c != (j - 1)))
						{
							is_mm = 1;
						}
					}
					if (!is_mm)
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
				}
			}
		}
		else if (c < 4)
		{
			if (is_multiref)
			{
				for (int b = 0; b < BASES_PER_NUCLEOTIDE; b++)
				{
					unsigned char base = nucl_bases_table[c][b];
					if (L[base] <= U[base])
					{
						heap_push_s(heap, e->i - 1, L[base], U[base], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
					}
				}
			}
			else
			{
				if (L[c + 1] <= U[c + 1])
				{
					heap_push_s(heap, e->i - 1, L[c + 1], U[c + 1], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps, e->aln_length, e->aln_path, params);
				}
			}
		}
	}
}

int inexact_match_fm_sid(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals,
						 const aln_params_t *params, diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns)
{

	int cdlen = params->no_e - params->no_s + 1;

	bwtint_t fm_n;
	Fm_GetN(fm, &fm_n);

	int countN = 0;
	for (int i = 0; i < readLen; i++)
	{
		if (read[i] > 3)
			countN++;
	}
	if (countN > params->max_diff)
	{
		return;
	}

	heap_reset(heap);
	if (precalc_sa_intervals != NULL)
	{
		if (precalc_sa_intervals->size != 0)
		{
			char aln_path[ALN_PATH_ALLOC] = {0};
			sa_intv_t *intv = precalc_sa_intervals->first_intv;
			for (int i = 0; i < precalc_sa_intervals->size; i++)
			{
				heap_push_s(heap, readLen - PRECALC_INTERVAL_LENGTH, intv->L, intv->U, 0, 0, 0, 0, 0, PRECALC_INTERVAL_LENGTH - 1, aln_path, params);
				intv = intv->next_intv;
			}
		}
		else
		{
			return;
		}
	}
	else
	{
		heap_push_s(heap, readLen, 0, fm_n, 0, 0, 0, 0, 0, 0, 0, params);
	}

	int best_score = aln_score(params->max_diff + 1, params->max_gapo + 1, params->max_gape + 1, params);
	int best_diff = params->max_diff + 1;
	int max_diff = params->max_diff;
	int num_best = 0;
	int max_entries = 0;

	int total_entries = 0;
	int last_num_entries = 1;
	int tid = 0;

	long nbt0 = nbt;
	while (heap->num_entries != 0)
	{
		total_entries += heap->num_entries - last_num_entries + 1;
		last_num_entries = heap->num_entries;

		if (heap->num_entries > max_entries)
			max_entries = heap->num_entries;
		if (heap->num_entries > params->max_entries)
		{
			break;
		}

		aln_entry_t e_;
		heap_pop(heap, &e_);
		aln_entry_t *e = &e_;
		nbt++;
		if (e->i == readLen - cdlen)
		{
			heap_reset(heap);
		}

		if (e->score > (best_score + params->mm_score))
		{
			break;
		}
		int diff_left = max_diff - e->num_mm - e->num_gapo - e->num_gape;
		if (diff_left < 0)
		{
			continue;
		}
		if ((e->i > 0) && (diff_left < D[e->i - 1].num_diff))
		{
			continue;
		}

		int diff_left_seed = params->max_diff_seed - e->num_mm - e->num_gapo - e->num_gape;
		int seed_index = e->i - (readLen - params->seed_length);

		if (e->i == 0)
		{
			int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
			if (alns->num_entries == 0)
			{
				best_score = score;
				best_diff = e->num_mm + e->num_gapo + e->num_gape;
				max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
			}
			if (score == best_score)
			{
				num_best += e->U - e->L + 1;
			}
			else if (num_best > params->max_best)
			{
				break;
			}
			add_alignment(e, e->L, e->U, score, alns, params);
			continue;
		}
		else if (diff_left == 0)
		{
			sa_intv_list_t *sa_intervals;
			if (exact_match_bounded_fm(fm, read, readLen, e->L, e->U, e->i - 1, &sa_intervals, params) > 0)
			{
				int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
				if (alns->num_entries == 0)
				{
					best_score = score;
					best_diff = e->num_mm + e->num_gapo + e->num_gape;
					max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
				}
				if (score == best_score)
				{
					sa_intv_t *intv = sa_intervals->first_intv;
					for (int k = 0; k < sa_intervals->size; k++)
					{
						num_best += intv->U - intv->L + 1;
						intv = intv->next_intv;
					}
				}
				else if (num_best > params->max_best)
				{
					break;
				}

				e->aln_length += e->i;
				sa_intv_t *intv = sa_intervals->first_intv;
				for (int k = 0; k < sa_intervals->size; k++)
				{
					add_alignment(e, intv->L, intv->U, score, alns, params);
					intv = intv->next_intv;
				}
			}

			free_sa_interval_list(sa_intervals);
			continue;
		}

		long L[ALPHABET_SIZE] = {0};
		long U[ALPHABET_SIZE] = {0};
		int alphabet_size = ALPHABET_SIZE;
		int is_multiref = params->is_multiref;
		if (params->is_multiref)
		{
			Occ16p(e->L - 1, e->U, L, U, fm);
			for (int i = 1; i < ALPHABET_SIZE; i++)
			{
				if (U[i] - L[i] + 1 == 0)
				{
					continue;
				}
			}
		}
		else
		{
		}

		int allow_diff = 1;
		int allow_indels = 1;
		int allow_mm = 1;
		int allow_open = 1;
		int allow_extend = 1;
		if ((readLen - e->i) < cdlen)
		{
			allow_diff = 0;
			allow_indels = 0;
			allow_mm = 0;
			allow_open = 0;
			allow_extend = 0;
		}
		if (e->i - 1 > 0)
		{
			if ((diff_left - 1) < D[e->i - 2].num_diff)
			{
				allow_diff = 0;
			}
			else if (((D[e->i - 1].num_diff == diff_left - 1) && (D[e->i - 2].num_diff == diff_left - 1)) && (D[e->i - 1].sa_intv_width == D[e->i - 2].sa_intv_width))
			{
				allow_mm = 0;
			}
		}
		int tmp = e->num_gapo + e->num_gape;
		if ((e->i - 1 < (params->no_indel_length + tmp)) || ((readLen - (e->i - 1)) < (params->no_indel_length + tmp)))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo) && (e->num_gape >= params->max_gape))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo))
		{
			allow_open = 0;
		}
		if ((e->num_gape >= params->max_gape))
		{
			allow_extend = 0;
		}

		if (allow_diff && allow_indels)
		{
			if (e->state == STATE_I)
			{
				if (allow_extend)
				{
					heap_push_s(heap, e->i - 1, e->L, e->U, e->num_mm, e->num_gapo, e->num_gape + 1, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
				}
			}
			else
			{
				if (allow_open && (e->state == STATE_M))
				{
					heap_push_s(heap, e->i - 1, e->L, e->U, e->num_mm, e->num_gapo + 1, e->num_gape, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
				}
				for (int j = 1; j < alphabet_size; j++)
				{
					if (L[j] <= U[j])
					{
						if (e->state == STATE_M)
						{
							if (allow_open)
							{
								heap_push_s(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo + 1, e->num_gape, STATE_D,
											e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
						else
						{
							if (allow_extend)
							{
								heap_push_s(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape + 1, STATE_D,
											e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
					}
				}
			}
		}

		int c = read[e->i - 1];
		if (allow_diff && allow_mm)
		{
			for (int j = 1; j < alphabet_size; j++)
			{
				if (j == 1 || j == 3 || j == 7 || j == 15)
				{
					continue;
				}
				if (L[j] <= U[j])
				{
					int is_mm = 0;
					if (is_multiref)
					{
						if ((c > 3) || (j == 10) /*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0))
						{
							is_mm = 1;
						}
					}
					else
					{
						if ((c > 3) || (c != (j - 1)))
						{
							is_mm = 1;
						}
					}
					if (!is_mm)
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
				}
			}
			for (int i = 1; i < 5; i++)
			{
				int j = pow(2, i) - 1;
				if (L[j] <= U[j])
				{
					int is_mm = 0;
					if (is_multiref)
					{
						if ((c > 3) || (j == 10) /*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0))
						{
							is_mm = 1;
						}
					}
					else
					{
						if ((c > 3) || (c != (j - 1)))
						{
							is_mm = 1;
						}
					}
					if (!is_mm)
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
				}
			}
		}
		else if (c < 4)
		{
			if (is_multiref)
			{
				for (int b = 0; b < BASES_PER_NUCLEOTIDE; b++)
				{
					unsigned char base = nucl_bases_table[c][b];
					if (base == 1 || base == 3 || base == 7 || base == 15)
					{
						continue;
					}
					if (L[base] <= U[base])
					{
						heap_push_s(heap, e->i - 1, L[base], U[base], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
					}
				}
				unsigned char base = nt4_gray[c];
				if (L[base] <= U[base])
				{
					heap_push_s(heap, e->i - 1, L[base], U[base], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
				}
			}
			else
			{
				if (L[c + 1] <= U[c + 1])
				{
					heap_push_s(heap, e->i - 1, L[c + 1], U[c + 1], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps, e->aln_length, e->aln_path, params);
				}
			}
		}
	}
	return 0;
}

int inexact_match_fm_sid2(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals,
						  const aln_params_t *params, diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns)
{

	if (D[readLen - 1].num_diff == 0)
	{
		sa_intv_t *intv = D[readLen - 1].sa_list->first_intv;
		aln_entry_t a_;
		memset(a_.aln_path, STATE_M, readLen);
		a_.num_mm = 0;
		a_.num_gapo = 0;
		a_.num_gape = 0;
		a_.num_snps = 0;
		a_.L = -1;
		a_.U = -1;
		a_.score = 0;
		a_.aln_length = readLen;
		for (int k = 0; k < D[readLen - 1].sa_list->size; k++)
		{
			add_alignment(&a_, intv->L, intv->U, 0, alns, params);
			intv = intv->next_intv;
		}
		free_sa_interval_list(D[readLen - 1].sa_list);
		return;
	}

	int cdlen = params->no_e - params->no_s + 1;

	bwtint_t fm_n;
	Fm_GetN(fm, &fm_n);

	int countN = 0;
	for (int i = 0; i < readLen; i++)
	{
		if (read[i] > 3)
			countN++;
	}
	if (countN > params->max_diff)
	{
		return;
	}

	heap_reset(heap);
	if (precalc_sa_intervals != NULL)
	{
		if (precalc_sa_intervals->size != 0)
		{
			char aln_path[ALN_PATH_ALLOC] = {0};
			sa_intv_t *intv = precalc_sa_intervals->first_intv;
			for (int i = 0; i < precalc_sa_intervals->size; i++)
			{
				heap_push_s(heap, readLen - PRECALC_INTERVAL_LENGTH, intv->L, intv->U, 0, 0, 0, 0, 0, PRECALC_INTERVAL_LENGTH - 1, aln_path, params);
				intv = intv->next_intv;
			}
		}
		else
		{
			return;
		}
	}
	else
	{
		heap_push_s(heap, readLen, 0, fm_n, 0, 0, 0, 0, 0, 0, 0, params);
	}

	int best_score = aln_score(params->max_diff + 1, params->max_gapo + 1, params->max_gape + 1, params);
	int best_diff = params->max_diff + 1;
	int max_diff = params->max_diff;
	int num_best = 0;
	int max_entries = 0;

	int total_entries = 0;
	int last_num_entries = 1;
	int tid = 0;

	long nbt0 = nbt;
	while (heap->num_entries != 0)
	{
		total_entries += heap->num_entries - last_num_entries + 1;
		last_num_entries = heap->num_entries;

		if (heap->num_entries > max_entries)
			max_entries = heap->num_entries;
		if (heap->num_entries > params->max_entries)
		{
			break;
		}

		aln_entry_t e_;
		heap_pop(heap, &e_);
		aln_entry_t *e = &e_;
		nbt++;
		if (e->i == readLen - cdlen)
		{
			heap_reset(heap);
		}

		if (e->score > (best_score + params->mm_score))
		{
			break;
		}
		int diff_left = max_diff - e->num_mm - e->num_gapo - e->num_gape;
		if (diff_left < 0)
		{
			continue;
		}
		if ((e->i > 0) && (diff_left < D[e->i - 1].num_diff))
		{
			continue;
		}
		int diff_left_seed = params->max_diff_seed - e->num_mm - e->num_gapo - e->num_gape;
		int seed_index = e->i - (readLen - params->seed_length);

		if (e->i == 0)
		{
			int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
			if (alns->num_entries == 0)
			{
				best_score = score;
				best_diff = e->num_mm + e->num_gapo + e->num_gape;
				max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
			}
			if (score == best_score)
			{
				num_best += e->U - e->L + 1;
			}
			else if (num_best > params->max_best)
			{
				break;
			}
			add_alignment(e, e->L + 10000000000, e->U + 10000000000, score, alns, params);
			continue;
		}
		else if (diff_left == 0)
		{
			sa_intv_list_t *sa_intervals;
			if (exact_match_bounded_fm(fm, read, readLen, e->L, e->U, e->i - 1, &sa_intervals, params) > 0)
			{
				int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
				if (alns->num_entries == 0)
				{
					best_score = score;
					best_diff = e->num_mm + e->num_gapo + e->num_gape;
					max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
				}
				if (score == best_score)
				{
					sa_intv_t *intv = sa_intervals->first_intv;
					for (int k = 0; k < sa_intervals->size; k++)
					{
						num_best += intv->U - intv->L + 1;
						intv = intv->next_intv;
					}
				}
				else if (num_best > params->max_best)
				{
					break;
				}

				e->aln_length += e->i;
				sa_intv_t *intv = sa_intervals->first_intv;
				for (int k = 0; k < sa_intervals->size; k++)
				{
					add_alignment(e, intv->L + 10000000000, intv->U + 10000000000, score, alns, params);
					intv = intv->next_intv;
				}
			}

			free_sa_interval_list(sa_intervals);
			continue;
		}

		long L[ALPHABET_SIZE] = {0};
		long U[ALPHABET_SIZE] = {0};
		int alphabet_size = ALPHABET_SIZE;
		int is_multiref = params->is_multiref;
		if (params->is_multiref)
		{
			Occ16p(e->L - 1, e->U, L, U, fm);
			for (int i = 1; i < ALPHABET_SIZE; i++)
			{
				if (U[i] - L[i] + 1 == 0)
				{
					continue;
				}
			}
		}
		else
		{
		}

		int allow_diff = 1;
		int allow_indels = 1;
		int allow_mm = 1;
		int allow_open = 1;
		int allow_extend = 1;
		if ((readLen - e->i) < cdlen)
		{
			allow_diff = 0;
			allow_indels = 0;
			allow_mm = 0;
			allow_open = 0;
			allow_extend = 0;
		}
		if (e->i - 1 > 0)
		{
			if ((diff_left - 1) < D[e->i - 2].num_diff)
			{
				allow_diff = 0;
			}
			else if (((D[e->i - 1].num_diff == diff_left - 1) && (D[e->i - 2].num_diff == diff_left - 1)) && (D[e->i - 1].sa_intv_width == D[e->i - 2].sa_intv_width))
			{
				allow_mm = 0;
			}
		}
		int tmp = e->num_gapo + e->num_gape;
		if ((e->i - 1 < (params->no_indel_length + tmp)) || ((readLen - (e->i - 1)) < (params->no_indel_length + tmp)))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo) && (e->num_gape >= params->max_gape))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo))
		{
			allow_open = 0;
		}
		if ((e->num_gape >= params->max_gape))
		{
			allow_extend = 0;
		}

		if (allow_diff && allow_indels)
		{
			if (e->state == STATE_I)
			{
				if (allow_extend)
				{
					heap_push_s(heap, e->i - 1, e->L, e->U, e->num_mm, e->num_gapo, e->num_gape + 1, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
				}
			}
			else
			{
				if (allow_open && (e->state == STATE_M))
				{
					heap_push_s(heap, e->i - 1, e->L, e->U, e->num_mm, e->num_gapo + 1, e->num_gape, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
				}
				for (int j = 1; j < alphabet_size; j++)
				{
					if (L[j] <= U[j])
					{
						if (e->state == STATE_M)
						{
							if (allow_open)
							{
								heap_push_s(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo + 1, e->num_gape, STATE_D,
											e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
						else
						{
							if (allow_extend)
							{
								heap_push_s(heap, e->i, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape + 1, STATE_D,
											e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
					}
				}
			}
		}

		int c = read[e->i - 1];
		if (allow_diff && allow_mm)
		{
			for (int j = 1; j < alphabet_size; j++)
			{
				if (j == 1 || j == 3 || j == 7 || j == 15)
				{
					continue;
				}
				if (L[j] <= U[j])
				{
					int is_mm = 0;
					if (is_multiref)
					{
						if ((c > 3) || (j == 10) /*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0))
						{
							is_mm = 1;
						}
					}
					else
					{
						if ((c > 3) || (c != (j - 1)))
						{
							is_mm = 1;
						}
					}
					if (!is_mm)
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
				}
			}
			for (int i = 1; i < 5; i++)
			{
				int j = pow(2, i) - 1;
				if (L[j] <= U[j])
				{
					int is_mm = 0;
					if (is_multiref)
					{
						if ((c > 3) || (j == 10) /*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0))
						{
							is_mm = 1;
						}
					}
					else
					{
						if ((c > 3) || (c != (j - 1)))
						{
							is_mm = 1;
						}
					}
					if (!is_mm)
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push_s(heap, e->i - 1, L[j], U[j], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
					}
				}
			}
		}
		else if (c < 4)
		{
			if (is_multiref)
			{
				for (int b = 0; b < BASES_PER_NUCLEOTIDE; b++)
				{
					unsigned char base = nucl_bases_table[c][b];
					if (base == 1 || base == 3 || base == 7 || base == 15)
					{
						continue;
					}
					if (L[base] <= U[base])
					{
						heap_push_s(heap, e->i - 1, L[base], U[base], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
					}
				}
				unsigned char base = nt4_gray[c];
				if (L[base] <= U[base])
				{
					heap_push_s(heap, e->i - 1, L[base], U[base], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
				}
			}
			else
			{
				if (L[c + 1] <= U[c + 1])
				{
					heap_push_s(heap, e->i - 1, L[c + 1], U[c + 1], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								e->num_snps, e->aln_length, e->aln_path, params);
				}
			}
		}
	}
	return 0;
}

int inexact_match_fm_mid(struct FMA *fm, char *read, int readLen, priority_heap_t *heap, sa_intv_list_t *precalc_sa_intervals,
						 const aln_params_t *params, diff_lower_bound_t *D, diff_lower_bound_t *D_seed, alns_t *alns)
{

	bwtint_t fm_n;
	Fm_GetN(fm, &fm_n);

	int countN = 0;
	for (int i = 0; i < readLen; i++)
	{
		if (read[i] > 3)
			countN++;
	}
	if (countN > params->max_diff)
	{
		return;
	}

	heap_reset(heap);
	if (precalc_sa_intervals != NULL)
	{
		if (precalc_sa_intervals->size != 0)
		{
		}
		else
		{
			return;
		}
	}
	else
	{
		heap_push(heap, params->no_e - 1, 0, fm_n, 0, fm_n, 0, 0, 0, 0, 0, 0, 0, params);
	}

	int best_score = aln_score(params->max_diff + 1, params->max_gapo + 1, params->max_gape + 1, params);
	int best_diff = params->max_diff + 1;
	int max_diff = params->max_diff;
	int num_best = 0;
	int max_entries = 0;

	int total_entries = 0;
	int last_num_entries = 1;
	int tid = 0;

	long nbt0 = nbt;
	while (heap->num_entries != 0)
	{
		total_entries += heap->num_entries - last_num_entries + 1;
		last_num_entries = heap->num_entries;

		if (heap->num_entries > max_entries)
			max_entries = heap->num_entries;
		if (heap->num_entries > params->max_entries)
		{
			break;
		}

		aln_entry_t e_;
		heap_pop(heap, &e_);
		aln_entry_t *e = &e_;
		nbt++;
		if (e->i == params->no_s)
		{
			heap_reset(heap);
		}

		if (e->score > (best_score + params->mm_score))
		{
			break;
		}
		int diff_left = max_diff - e->num_mm - e->num_gapo - e->num_gape;
		if (diff_left < 0)
		{
			continue;
		}

		if (e->i == readLen)
		{
			int score = aln_score(e->num_mm, e->num_gapo, e->num_gape, params);
			if (alns->num_entries == 0)
			{
				best_score = score;
				best_diff = e->num_mm + e->num_gapo + e->num_gape;
				max_diff = (best_diff + 1 > params->max_diff) ? params->max_diff : best_diff + 1;
			}
			if (score == best_score)
			{
				num_best += e->U - e->L + 1;
			}
			else if (num_best > params->max_best)
			{
				break;
			}
			add_alignment(e, e->L, e->U - 1, score, alns, params);
			continue;
		}

		long L[ALPHABET_SIZE] = {0};
		long U[ALPHABET_SIZE] = {0};
		bwtint_t l[ALPHABET_SIZE] = {0};
		bwtint_t u[ALPHABET_SIZE] = {0};
		int alphabet_size = ALPHABET_SIZE;
		int is_multiref = params->is_multiref;
		if (params->is_multiref)
		{
			if (e->i == params->no_e)
			{
				Occ16p(e->l, e->u, l, u, fm);
				for (int j = 0; j < ALPHABET_SIZE; j++)
				{
					L[j] = (j == 0) ? (e->L) : (U[j - 1]);
					U[j] = L[j] + u[iupacCompl[j]] - l[iupacCompl[j]] + 1;
				}
			}
			else if (e->i < params->no_e)
			{
				Occ16p(e->L - 1, e->U, L, U, fm);
				for (int j = 0; j < ALPHABET_SIZE; j++)
				{
					l[j] = (j == 0) ? (e->l) : (u[j - 1]);
					u[j] = l[j] + U[iupacCompl[j]] - L[iupacCompl[j]] + 1;
				}
			}
			else
			{
				Occ16p(e->l - 1, e->u, l, u, fm);
				for (int j = 0; j < ALPHABET_SIZE; j++)
				{
					L[j] = (j == 0) ? (e->L) : (U[j - 1]);
					U[j] = L[j] + u[iupacCompl[j]] - l[iupacCompl[j]] + 1;
				}
			}
		}
		else
		{
		}

		int allow_diff = 1;
		int allow_indels = 1;
		int allow_mm = 1;
		int allow_open = 1;
		int allow_extend = 1;
		if ((e->i > params->no_s) && (e->i < params->no_e))
		{
			allow_diff = 0;
			allow_indels = 0;
			allow_mm = 0;
			allow_open = 0;
			allow_extend = 0;
		}
		if (diff_left == 0)
		{
			allow_diff = 0;
		}

		int tmp = e->num_gapo + e->num_gape;
		if ((e->i - 1 < (params->no_indel_length + tmp)) || ((readLen - (e->i - 1)) < (params->no_indel_length + tmp)))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo) && (e->num_gape >= params->max_gape))
		{
			allow_indels = 0;
		}
		if ((e->num_gapo >= params->max_gapo))
		{
			allow_open = 0;
		}
		if ((e->num_gape >= params->max_gape))
		{
			allow_extend = 0;
		}

		if (allow_diff && allow_indels)
		{
			if (e->state == STATE_I)
			{
				if (allow_extend)
				{
					if (e->i == 0)
					{
						heap_push(heap, params->no_e, e->L, e->U, e->l, e->u, e->num_mm, e->num_gapo, e->num_gape + 1, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
					}
					else if (e->i < params->no_e)
					{
						heap_push(heap, e->i - 1, e->L, e->U, e->l, e->u, e->num_mm, e->num_gapo, e->num_gape + 1, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push(heap, e->i + 1, e->L, e->U, e->l, e->u, e->num_mm, e->num_gapo, e->num_gape + 1, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
					}
				}
			}
			else
			{
				if (allow_open && (e->state == STATE_M))
				{
					if (e->i == 0)
					{
						heap_push(heap, params->no_e, e->L, e->U, e->l, e->u, e->num_mm, e->num_gapo + 1, e->num_gape, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
					}
					else if (e->i < params->no_e)
					{
						heap_push(heap, e->i - 1, e->L, e->U, e->l, e->u, e->num_mm, e->num_gapo + 1, e->num_gape, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push(heap, e->i + 1, e->L, e->U, e->l, e->u, e->num_mm, e->num_gapo + 1, e->num_gape, STATE_I, e->num_snps, e->aln_length, e->aln_path, params);
					}
				}
				for (int j = 1; j < alphabet_size; j++)
				{
					if (((e->i > 0) && (e->i < params->no_e)) ? (L[j] <= U[j]) : (l[iupacCompl[j]] <= u[iupacCompl[j]]))
					{
						if (e->state == STATE_M)
						{
							if (allow_open)
							{
								heap_push(heap, e->i, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo + 1, e->num_gape, STATE_D,
										  e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
						else
						{
							if (allow_extend)
							{
								heap_push(heap, e->i, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo, e->num_gape + 1, STATE_D,
										  e->num_snps, e->aln_length, e->aln_path, params);
							}
						}
					}
				}
			}
		}

		int c = (e->i == 0) ? read[0] : read[e->i];
		if (allow_diff && allow_mm)
		{
			for (int j = 1; j < alphabet_size; j++)
			{
				if (j == 1 || j == 3 || j == 7 || j == 15)
				{
					continue;
				}
				if (((e->i > 0) && (e->i < params->no_e)) ? (L[j] <= U[j]) : (l[iupacCompl[j]] <= u[iupacCompl[j]]))
				{
					int is_mm = 0;
					if (is_multiref)
					{
						if ((c > 3) || (j == 10) /*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0))
						{
							is_mm = 1;
						}
					}
					else
					{
						if ((c > 3) || (c != (j - 1)))
						{
							is_mm = 1;
						}
					}
					if (!is_mm)
					{
						if (e->i == 0)
						{
							heap_push(heap, params->no_e, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else if (e->i < params->no_e)
						{
							heap_push(heap, e->i - 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else
						{
							heap_push(heap, e->i + 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
					}
					else
					{
						if (e->i == 0)
						{
							heap_push(heap, params->no_e, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else if (e->i < params->no_e)
						{
							heap_push(heap, e->i - 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else
						{
							heap_push(heap, e->i + 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
					}
				}
			}

			for (int i = 1; i < 5; i++)
			{
				int j = pow(2, i) - 1;
				if (((e->i > 0) && (e->i < params->no_e)) ? (L[j] <= U[j]) : (l[iupacCompl[j]] <= u[iupacCompl[j]]))
				{
					int is_mm = 0;
					if (is_multiref)
					{
						if ((c > 3) || (j == 10) /*N*/ || ((nt4_gray_val[c] & grayVal[j]) == 0))
						{
							is_mm = 1;
						}
					}
					else
					{
						if ((c > 3) || (c != (j - 1)))
						{
							is_mm = 1;
						}
					}
					if (!is_mm)
					{
						if (e->i == 0)
						{
							heap_push(heap, params->no_e, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else if (e->i < params->no_e)
						{
							heap_push(heap, e->i - 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else
						{
							heap_push(heap, e->i + 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
					}
					else
					{
						if (e->i == 0)
						{
							heap_push(heap, params->no_e, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else if (e->i < params->no_e)
						{
							heap_push(heap, e->i - 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
						else
						{
							heap_push(heap, e->i + 1, L[j], U[j], l[iupacCompl[j]], u[iupacCompl[j]], e->num_mm + 1, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_multiref && is_snp[j]), e->aln_length, e->aln_path, params);
						}
					}
				}
			}
		}
		else if (c < 4)
		{
			if (is_multiref)
			{
				for (int b = 0; b < BASES_PER_NUCLEOTIDE; b++)
				{
					unsigned char base = nucl_bases_table[c][b];
					if (base == 1 || base == 3 || base == 7 || base == 15)
					{
						continue;
					}
					if (((e->i > 0) && (e->i < params->no_e)) ? (L[base] <= U[base]) : (l[iupacCompl[base]] <= u[iupacCompl[base]]))
					{
						if (e->i == 0)
						{
							heap_push(heap, params->no_e, L[base], U[base], l[iupacCompl[base]], u[iupacCompl[base]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
						}
						else if (e->i < params->no_e)
						{
							heap_push(heap, e->i - 1, L[base], U[base], l[iupacCompl[base]], u[iupacCompl[base]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
						}
						else
						{
							heap_push(heap, e->i + 1, L[base], U[base], l[iupacCompl[base]], u[iupacCompl[base]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
									  e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
						}
					}
				}
				unsigned char base = nt4_gray[c];
				if (((e->i > 0) && (e->i < params->no_e)) ? (L[base] <= U[base]) : (l[iupacCompl[base]] <= u[iupacCompl[base]]))
				{
					if (e->i == 0)
					{
						heap_push(heap, params->no_e, L[base], U[base], l[iupacCompl[base]], u[iupacCompl[base]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								  e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
					}
					else if (e->i < params->no_e)
					{
						heap_push(heap, e->i - 1, L[base], U[base], l[iupacCompl[base]], u[iupacCompl[base]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								  e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
					}
					else
					{
						heap_push(heap, e->i + 1, L[base], U[base], l[iupacCompl[base]], u[iupacCompl[base]], e->num_mm, e->num_gapo, e->num_gape, STATE_M,
								  e->num_snps + (is_snp[base]), e->aln_length, e->aln_path, params);
					}
				}
			}
			else
			{
			}
		}
	}
	return 0;
}

priority_heap_t *heap_init(const aln_params_t *p)
{
	priority_heap_t *heap = (priority_heap_t *)calloc(1, sizeof(priority_heap_t));
	heap->num_buckets = aln_score(p->max_diff + 1, p->max_gapo + 1, p->max_gape + 1, p);
	heap->buckets = (heap_bucket_t *)calloc(heap->num_buckets, sizeof(heap_bucket_t));
	if (heap == NULL || heap->buckets == NULL)
	{
		printf("Could not allocate memory for the heap \n");
		exit(1);
	}
	for (int i = 0; i < heap->num_buckets; i++)
	{
		heap_bucket_t *hb = &(heap->buckets[i]);
		hb->max_entries = 4;
		hb->entries = (aln_entry_t *)calloc(hb->max_entries, sizeof(aln_entry_t));
		if (hb->entries == NULL)
		{
			printf("Could not allocate memory for the heap \n");
			exit(1);
		}
	}
	heap->best_score = heap->num_buckets;
	return heap;
}

void heap_free(priority_heap_t *heap)
{
	for (int i = 0; i < heap->num_buckets; i++)
	{
		free(heap->buckets[i].entries);
	}
	free(heap->buckets);
	free(heap);
}

void heap_reset(priority_heap_t *heap)
{
	for (int i = 0; i < heap->num_buckets; i++)
	{
		heap->buckets[i].num_entries = 0;
	}
	heap->best_score = heap->num_buckets;
	heap->num_entries = 0;
}

void heap_push_s(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U, /*const bwtint_t l, const bwtint_t u, gll*/
				 const int num_mm, const int num_gapo, const int num_gape,
				 const int state, const int num_snps, const int aln_length, const char *aln_path,
				 const aln_params_t *params)
{

	int score = aln_score(num_mm, num_gapo, num_gape, params);
	heap_bucket_t *hb = &(heap->buckets[score]);
	if (hb->num_entries == hb->max_entries)
	{
		hb->max_entries <<= 1;
		hb->entries = (aln_entry_t *)realloc(hb->entries, sizeof(aln_entry_t) * hb->max_entries);
		if (hb->entries == NULL)
		{
			printf("Could not reallocate memory for the heap bucket! \n");
			exit(1);
		}
	}
	aln_entry_t *p = &(hb->entries[hb->num_entries]);
	p->i = i;
	p->score = score;
	p->L = L;
	p->U = U;
	// p->l = l;
	// p->u = u;
	p->num_mm = num_mm;
	p->num_gapo = num_gapo;
	p->num_gape = num_gape;
	p->state = state;
	p->num_snps = num_snps;
	p->aln_length = 0;
	if (aln_path != NULL)
	{
		memset(&(p->aln_path), 0, ALN_PATH_ALLOC * sizeof(char));
		memcpy(&(p->aln_path), aln_path, aln_length * sizeof(char));
		p->aln_path[aln_length] = state;
		p->aln_length = aln_length + 1;
	}

	hb->num_entries++;
	heap->num_entries++;

	if (heap->best_score > score)
	{
		heap->best_score = score;
	}
}

void heap_push(priority_heap_t *heap, const int i, const bwtint_t L, const bwtint_t U, const bwtint_t l, const bwtint_t u,
			   const int num_mm, const int num_gapo, const int num_gape,
			   const int state, const int num_snps, const int aln_length, const char *aln_path,
			   const aln_params_t *params)
{

	int score = aln_score(num_mm, num_gapo, num_gape, params);
	heap_bucket_t *hb = &(heap->buckets[score]);
	if (hb->num_entries == hb->max_entries)
	{
		hb->max_entries <<= 1;
		hb->entries = (aln_entry_t *)realloc(hb->entries, sizeof(aln_entry_t) * hb->max_entries);
		if (hb->entries == NULL)
		{
			printf("Could not reallocate memory for the heap bucket! \n");
			exit(1);
		}
	}
	aln_entry_t *p = &(hb->entries[hb->num_entries]);
	p->i = i;
	p->score = score;
	p->L = L;
	p->U = U;
	p->l = l;
	p->u = u;
	p->num_mm = num_mm;
	p->num_gapo = num_gapo;
	p->num_gape = num_gape;
	p->state = state;
	p->num_snps = num_snps;
	p->aln_length = 0;
	if (aln_path != NULL)
	{
		memset(&(p->aln_path), 0, ALN_PATH_ALLOC * sizeof(char));
		memcpy(&(p->aln_path), aln_path, aln_length * sizeof(char));
		p->aln_path[aln_length] = state;
		p->aln_length = aln_length + 1;
	}

	hb->num_entries++;
	heap->num_entries++;

	if (heap->best_score > score)
	{
		heap->best_score = score;
	}
}

void heap_pop(priority_heap_t *heap, aln_entry_t *e)
{
	heap_bucket_t *hb = &(heap->buckets[heap->best_score]);
	aln_entry_t *et = &(hb->entries[hb->num_entries - 1]);
	hb->num_entries--;
	heap->num_entries--;

	if ((hb->num_entries == 0) && heap->num_entries)
	{
		int i;
		for (i = heap->best_score + 1; i < heap->num_buckets; i++)
		{
			if (heap->buckets[i].num_entries != 0)
				break;
		}
		heap->best_score = i;
	}
	else if (heap->num_entries == 0)
	{
		heap->best_score = heap->num_buckets;
	}
	memcpy(e, et, sizeof(aln_entry_t));
}
