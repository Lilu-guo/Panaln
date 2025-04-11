#ifndef BWTSNP_BWT_H
#define BWTSNP_BWT_H
#include <stdint.h>
#include "common.h"
#include "io.h"
#include <sys/stat.h>  
#define OCC_INTERVAL 128
#define SA_INTERVAL 32
typedef struct {
	bwtint_t length;
	bwtint_t num_words;
	uint32_t *bwt;
	bwtint_t C[ALPHABET_SIZE + 1];
	bwtint_t* O;
	bwtint_t num_occ;
	uint8_t occ_count_table[(1 << (BITS_IN_BYTE << 1))];
	bwtint_t *SA;
	bwtint_t num_sa;
	bwtint_t sa0_index;
} bwt_t;
int index_bwt(const char *fastaFname, const char* extSAFname);
bwt_t* construct_bwt(unsigned char *seq, const bwtint_t length, const char* extSAFname);
void free_bwt(bwt_t* BWT);
void store_bwt(const bwt_t* BWT, const char* bwtFname);
bwt_t* load_bwt(const char* bwtFname, const int loadSA);
void print_bwt(const bwt_t *BWT);
bwtint_t is_bwt(unsigned char *T, bwtint_t n, bwtint_t* SA); 
unsigned char B(const bwt_t* BWT, const bwtint_t i);
bwtint_t O(const bwt_t* BWT, const unsigned char c, bwtint_t i);
void O_alphabet(const bwt_t* BWT, const bwtint_t i, int alphabet_size, bwtint_t* occ, int L);
void O_actg_alphabet(const bwt_t* BWT, const bwtint_t i, bwtint_t* occ, int L);
void O_LU(const bwt_t* BWT, const unsigned char c, const bwtint_t L, const bwtint_t U, bwtint_t* occL, bwtint_t* occU);
bwtint_t C(const bwt_t* BWT, const unsigned char c);
bwtint_t SA(const bwt_t* BWT, const bwtint_t i);
bwtint_t invPsi(const bwt_t* BWT, const bwtint_t i);
#endif
