#ifndef IO_H_
#define IO_H_
#ifndef myPrint
#define jj 
#else
#define jj printf
#endif
#include "common.h"
#define SEQ_BATCH_ALLOC_LEN			262144
#define READ_LENGTH_ALLOC           15000 
#define ANN_ALLOC_LEN				256
#define MAX_SEQ_NAME_LEN 			256
#define NUM_READS_ALLOC 			1000000
#define BITS_PER_CHAR               4 
#define CHARS_PER_BYTE              (BITS_IN_BYTE / BITS_PER_CHAR)
#define CHARS_PER_WORD              8
#define CHARS_PER_128BITS           32
#define BYTES_IN_WORD 				4
#define BITS_IN_WORD 				32
#define BITS_IN_BYTE 				8
#define CHAR_MASK					15 
#define CHAR_COUNT_2BIT_MASK		3
#define LS_2BYTES_MASK				65535
#define MS_2BYTES_MASK				(LS_2BYTES_MASK << (2*BITS_IN_BYTE)) 
#define ALPHABET_SIZE 16
static const unsigned char iupacChar[16] =  {'$', 'T', 'K', 'G', 'S', 'B', 'Y', 'C', 'M', 'H', 'N', 'V', 'R', 'D', 'W', 'A'};  
static const unsigned char grayVal[16] =    { 0,   1,   3,   2,   6,   7,   5,   4,   12,  13,  15,  14,  10,  11,  9,   8 };  
static const unsigned char iupacGrayOrd[16]={ 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15};  
static const unsigned char grayIupac[16] =  { 15,  0,   2,   1,   6,   5,   3,   4,   14,  13,  11,  12,  7,   8,   10,  9 };  
static const unsigned char iupacCompl[16] = { 0,   15,  8,   7,   4,   11,  12,  3,   2,   13,  10,  5,   6,   9,   14,  1};   
static const unsigned char is_snp[16] = 	{ 0,   0,   1,   0,   1,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1,   0};   
static const unsigned char is_snp_$[16] = 	{ 1,   0,   1,   0,   1,   1,   1,   0,   1,   1,   1,   1,   1,   1,   1,   0};   
static const uint32_t char_masks[16] = { 0x0, 0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555, 0x66666666,
											  0x77777777, 0x88888888, 0x99999999, 0xAAAAAAAA, 0xBBBBBBBB, 0xCCCCCCCC,
											  0xDDDDDDDD, 0xEEEEEEEE, 0xFFFFFFFF};
static const uint32_t char_masks128[16][4] __attribute__((aligned(16))) =
													  { {0x0, 		0x0, 		0x0, 		0x0},
														{0x11111111, 0x11111111, 0x11111111, 0x11111111},
														{0x22222222, 0x22222222, 0x22222222, 0x22222222},
														{0x33333333, 0x33333333, 0x33333333, 0x33333333},
														{0x44444444, 0x44444444, 0x44444444, 0x44444444},
														{0x55555555, 0x55555555, 0x55555555, 0x55555555},
														{0x66666666, 0x66666666, 0x66666666, 0x66666666},
														{0x77777777, 0x77777777, 0x77777777, 0x77777777},
														{0x88888888, 0x88888888, 0x88888888, 0x88888888},
														{0x99999999, 0x99999999, 0x99999999, 0x99999999},
														{0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA},
														{0xBBBBBBBB, 0xBBBBBBBB, 0xBBBBBBBB, 0xBBBBBBBB},
														{0xCCCCCCCC, 0xCCCCCCCC, 0xCCCCCCCC, 0xCCCCCCCC},
														{0xDDDDDDDD, 0xDDDDDDDD, 0xDDDDDDDD, 0xDDDDDDDD},
														{0xEEEEEEEE, 0xEEEEEEEE, 0xEEEEEEEE, 0xEEEEEEEE},
														{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}};
#define CHAR_MASK_A 0xFFFFFFFF
#define CHAR_MASK_G 0x33333333
#define CHAR_MASK_C 0x77777777
#define CHAR_MASK_T 0x11111111
static const uint32_t word_fractions_masks128[32][4] __attribute__((aligned(16))) =
													 { {0x0FFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x000FFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x0000FFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000FFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x000000FF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x0000000F, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x0FFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00FFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x000FFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x0000FFFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000FFF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x000000FF, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x0000000F, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x0FFFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00FFFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x000FFFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x0000FFFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00000FFF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x000000FF, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x0000000F, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x0FFFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x00FFFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x000FFFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x0000FFFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x00000FFF},
													   {0x00000000, 0x00000000, 0x00000000, 0x000000FF},
													   {0x00000000, 0x00000000, 0x00000000, 0x0000000F},
													   {0x00000000, 0x00000000, 0x00000000, 0x00000000},
													};
static const uint32_t word_fractions_masks[8] = { 0x0FFFFFFF, 0x00FFFFFF, 0x000FFFFF, 0x0000FFFF,
                                                  0x00000FFF, 0x000000FF, 0x0000000F, 0x0};
#define NUM_NUCLEOTIDES 4
#define BASES_PER_NUCLEOTIDE 7 
static const unsigned char nucl_bases_table[NUM_NUCLEOTIDES][BASES_PER_NUCLEOTIDE] = {
			{8,  9, 11, 12, 13, 14, 15} ,                                        
			{2,  3,  4,  5, 11, 12, 13} ,
			{4,  5,  6,  7, 8,  9,  11} ,
			{1,  2,  5,  6, 9,  13, 14} };
static const unsigned char nt4_gray[5] = {15, 3, 7, 1, 10};      
static const unsigned char nt4_gray_val[5] = {8, 2, 4, 1, 15};   
static const unsigned char nt4_complement[5] = {3, 2, 1, 0, 4};  
static const unsigned char nt4_table[256] = {                              
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4,0,4,2,4,4,4,1,4,4,4,4,4,4,4,4,
	4, 4, 4, 4,  3,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4,0,4,2,4,4,4,1,4,4,4,4,4,4,4,4,
	4, 4, 4, 4,  3,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
static const unsigned char nt16_table[256] = {                             
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    0,10,10,10, 10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 15,5,7,  13,10, 10,3,      9,10,10,2,  10,8,10,10,
	10, 10, 12,4,     1,10,11,14,  10,6,10,10,      10, 10, 10, 10,
	10, 15,5,7,  13,10, 10,3,      9,10,10,2,  10,8,10,10,
	10, 10, 12,4,     1,10,11,14,  10,6,10,10,      10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
};
typedef struct {                                                              
	int score;
	bwtint_t L, U;
	int num_mm;
	int num_gapo;
	int num_gape;
	int num_snps;
	int aln_length;
	char* aln_path;
} aln_t;
typedef struct {
	int num_entries;
	int max_entries;
	aln_t* entries;
} alns_t;
typedef struct {                                           
	int len;
	char name[MAX_SEQ_NAME_LEN+1];
	char* seq;                                             
	char* rc;
	char* qual;
	int strand;
	bwtint_t ref_pos_l;
	bwtint_t ref_pos_r;
	bwtint_t* mref_pos;
	int num_mref_pos;
	alns_t* alns; 
	int aln_type;
	int aln_top1_count;
	int aln_top2_count;
	int mapQ; 
	int num_mm;
	int num_gapo;
	int num_gape;
	int aln_strand; 
	int aln_score; 
	bwtint_t aln_pos; 
	bwtint_t aln_sa;
	int aln_length;
	char* aln_path;
} read_t;
typedef struct {
	unsigned int count;                                               
	unsigned int max_len;
	read_t* reads;                                                    
} reads_t;                                                            
typedef struct {                                                   
	char name[MAX_SEQ_NAME_LEN];
	bwtint_t start_index;                                             
	bwtint_t end_index;
} seq_annotation_t;
typedef struct {                                                   
	int num_seq;
	seq_annotation_t* seq_anns;                                       
} fasta_annotations_t;
typedef struct
{
	int ann;
	long long A;
	long long B_minus_A;
	long long C;
	long long D_minus_C;
	int ref_len;
	int alt_len;
} bubble_t;
void fasta2ref(const char *fastaFname, const char* refFname, const char* annFname, unsigned char** seq, bwtint_t *totalSeqLen);
void ref2seq(const char* refFname, unsigned char** seq, bwtint_t* seqLen);
void fasta2pac(const char *fastaFname, const char* pacFname, const char* annFname);
void pac2seq(const char *pacFname, unsigned char** seq, bwtint_t *seqLength);
reads_t* fastq2reads(const char *readsFname);
void seq2rev_compl(unsigned char* seq, bwtint_t seqLen, unsigned char** rcSeq);
void parse_read_mapping(read_t* read);
void print_read(read_t* read);
void free_reads(reads_t* reads);
bubble_t* bubf2bub(const char *bub);      
fasta_annotations_t* annf2ann(const char *annFname);                             
fasta_annotations_t* annf2Ann(const char *annFname, uint32_t* nchrom, uint64_t* pchrom); 
void free_ann(fasta_annotations_t* annotations);
void pack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length);
void unpack_byte(const unsigned char *input, unsigned char *output, const bwtint_t length);
void pack_word(const unsigned char *input, unsigned int *output, const bwtint_t length);
void unpack_word(const unsigned char *input, unsigned char *output, const bwtint_t length);
#endif
