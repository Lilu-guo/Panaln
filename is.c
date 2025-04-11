#include <stdlib.h>
#include "bwt.h"
#include "io.h"
typedef unsigned char ubyte_t;
typedef long int seqint_t;
#define chr(i) (cs == sizeof(seqint_t) ? ((const seqint_t *)T)[i]:((const unsigned char *)T)[i])
static void getCounts(const unsigned char *T, seqint_t *C, seqint_t n, seqint_t k, int cs)
{
	seqint_t i;
	for (i = 0; i < k; ++i) C[i] = 0;
	for (i = 0; i < n; ++i) ++C[chr(i)];
}
static void getBuckets(const seqint_t *C, seqint_t *B, seqint_t k, seqint_t end)
{
	seqint_t i, sum = 0;
	if (end) {
		for (i = 0; i < k; ++i) {
			sum += C[i];
			B[i] = sum;
		}
	} else {
		for (i = 0; i < k; ++i) {
			sum += C[i];
			B[i] = sum - C[i];
		}
	}
}
static void induceSA(const unsigned char *T, seqint_t *SA, seqint_t *C, seqint_t *B, seqint_t n, seqint_t k, int cs)
{
	seqint_t *b;
	seqint_t  c0, c1;
	seqint_t i;
	seqint_t j;
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 0);	
	j = n - 1;
        b = SA + B[c1 = chr(j)];
	*b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
	for (i = 0; i < n; ++i) {
		j = SA[i], SA[i] = ~j;
		if (0 < j) {
			--j;
			if ((c0 = chr(j)) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}
			*b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
		}
	}
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	
	for (i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
		if (0 < (j = SA[i])) {
			--j;
			if ((c0 = chr(j)) != c1) {
				B[c1] = b - SA;
				b = SA + B[c1 = c0];
			}
			*--b = ((j == 0) || (chr(j - 1) > c1)) ? ~j : j;
		} else SA[i] = ~j;
	}
}
static int sais_main(const unsigned char *T, seqint_t *SA, seqint_t fs, seqint_t n, seqint_t k, int cs)
{
	seqint_t *C, *B, *RA;
	seqint_t  j, c, m, p, q, plen, qlen, name;
	long int i;
	seqint_t  c0, c1;
	seqint_t  diff;
	if (k <= fs) {
		C = SA + n;
		B = (k <= (fs - k)) ? C + k : C;
	} else if ((C = B = (seqint_t *) malloc(k * sizeof(seqint_t))) == NULL) return -2;
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	
	for (i = 0; i < n; ++i) SA[i] = 0;
	for (i = n - 2, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr(i)) < (c1 + c)) c = 1;
		else if (c != 0) SA[--B[c1]] = i + 1, c = 0;
	}
	induceSA(T, SA, C, B, n, k, cs);
	if (fs < k) free(C);
	for (i = 0, m = 0; i < n; ++i) {
		p = SA[i];
		if ((0 < p) && (chr(p - 1) > (c0 = chr(p)))) {
			for (j = p + 1; (j < n) && (c0 == (c1 = chr(j))); ++j);
			if ((j < n) && (c0 < c1)) SA[m++] = p;
		}
	}
	for (i = m; i < n; ++i) SA[i] = 0;	
	for (i = n - 2, j = n, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr(i)) < (c1 + c)) c = 1;
		else if (c != 0) {
			SA[m + ((i + 1) >> 1)] = j - i - 1;
			j = i + 1;
			c = 0;
		}
	}
	for (i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
		p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
		if (plen == qlen) {
			for (j = 0; (j < plen) && (chr(p + j) == chr(q + j)); j++);
			if (j == plen) diff = 0;
		}
		if (diff != 0) ++name, q = p, qlen = plen;
		SA[m + (p >> 1)] = name;
	}
	if (name < m) {
		RA = SA + n + fs - m;
		for (i = n - 1, j = m - 1; m <= i; --i) {
			if (SA[i] != 0) RA[j--] = SA[i] - 1;
		}
		if (sais_main((unsigned char *) RA, SA, fs + n - m * 2, m, name, sizeof(seqint_t)) != 0) return -2;
		for (i = n - 2, j = m - 1, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
			if ((c0 = chr(i)) < (c1 + c)) c = 1;
			else if (c != 0) RA[j--] = i + 1, c = 0; 
		}
		for (i = 0; i < m; ++i) SA[i] = RA[SA[i]]; 
	}
	if (k <= fs) {
		C = SA + n;
		B = (k <= (fs - k)) ? C + k : C;
	} else if ((C = B = (seqint_t *) malloc(k * sizeof(seqint_t))) == NULL) return -2;
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	
	for (i = m; i < n; ++i) SA[i] = 0; 
	for (i = m - 1; 0 <= i; --i) {
		j = SA[i], SA[i] = 0;
		SA[--B[chr(j)]] = j;
	}
	induceSA(T, SA, C, B, n, k, cs);
	if (fs < k) free(C);
	return 0;
}
int is_sa(const ubyte_t *T, seqint_t *SA, seqint_t n)
{
	if ((T == NULL) || (SA == NULL) || (n < 0)) return -1;
	SA[0] = n;
	if (n <= 1) {
		if (n == 1) SA[1] = 0;
		return 0;
	}
	return sais_main(T, SA+1, 0, n, 256, 1);
}
bwtint_t is_bwt(unsigned char *T, bwtint_t n, bwtint_t* bwtSA) {
	seqint_t *SA = (seqint_t*) calloc(n+1, sizeof(seqint_t));
	if(SA == 0) {
		printf("is_bwt: Could not allocate memory for the BWT index construction (in SAIS), alloc'ing %" PRIbwtint_t " Mb\n", (n+1)*sizeof(seqint_t)/1024/1024);
		exit(1);
	}
	is_sa(T, SA, n);                                              
	bwtint_t sa0_index;
	for (bwtint_t i = 0; i <= n; i++) {
		if(i % SA_INTERVAL == 0) {
			bwtSA[i/SA_INTERVAL] = SA[i];
		}
		if (SA[i] == 0) {                                         
			sa0_index = i;
			SA[i] = nt16_table['$'];                              
		}
		else {
			SA[i] = T[SA[i] - 1];                                 
		}
	}
	for (bwtint_t i = 0; i <= n; i++) {
		T[i] = SA[i];
	}
	free(SA);
	return sa0_index;
}
