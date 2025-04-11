/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: WFA Sample-Code
 */

#include "utils/commons.h"
#include "wavefront/wavefront_align.h"
static const unsigned char iupacChar[16] =  {'$', 'T', 'K', 'G', 'S', 'B', 'Y', 'C', 'M', 'H', 'N', 'V', 'R', 'D', 'W', 'A'};  //iupac，转asc码
static const unsigned char grayVal[16] =    { 0,   1,   3,   2,   6,   7,   5,   4,   12,  13,  15,  14,  10,  11,  9,   8 };  //iupac/下标，转gray
static const unsigned char nt16_table[256] = {                                                                                 //将明文的ASC码（16字符），转为Fig3的iupac码（0~15, A:15
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    0/*$*/,10,10,10, 10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,    10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 15/*A*/,5/*B*/,7/*C*/,  13/*D*/,10, 10,3/*G*/,      9/*H*/,10,10,2/*K*/,  10,8/*M*/,10/*N*/,10,
	10, 10, 12/*R*/,4/*S*/,     1/*T*/,10,11/*V*/,14/*W*/,  10,6/*Y*/,10,10,      10, 10, 10, 10,
	10, 15/*a*/,5/*b*/,7/*c*/,  13/*d*/,10, 10,3/*g*/,      9/*h*/,10,10,2/*k*/,  10,8/*m*/,10/*n*/,10,
	10, 10, 12/*r*/,4/*s*/,     1/*t*/,10,11/*v*/,14/*w*/,  10,6/*y*/,10,10,      10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
	10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,  10, 10, 10, 10,
};

int main(int argc,char* argv[]) {
  // Patter & Text
  // char* pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
  // char* text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
  // char* pattern    =         "TTTGAGACACTCTTTGTATAGCATGTGGAAATGGATATTTGGAGCGCTTTGAGGCCCATGCTGAAGAAGGAAATATCTTCCCAAAAAAACTAGACGAAAGCATTCTCGGAATCTTGTTTGCCATGTGTGTACTCAACTAACAGAGTTGAACCTATCTTTTGACAGAGCAGTTTTGAAACACTCTTTTTGTGGAATCTGCAAGTGGATATTTGGATAGCTTCGAGGATTTCGTTGGAAACGGGAATATCCT";
  // char* text = "ATTCATACAGCAGGTTTGAGACACTCTTTGTATAGCATGCGGAAATGGATATTTGGAGCGCTTTGAGGCCTATGGTGAAGAAGGAAATATCTTCCCAAAAAAACTAGACGAAAGCATTCTCGGAATCTTGTTTGCCATGTGTGTACTCAACTAACAGAGTTGAACCTATCTTTTGACAGAGCAGTTTTGAAACACTCTTTTTGTGGAATCTGCAAGTGGATATTTGGATAGCTTCGAGGATTTCGTTGGAAACGGGAATATCCTCATTTAAAATCTAGAT";
  char* pattern0    =    "ACGCT"; //中间的'G'
  // char* text0 =       "TTCACACTAGT"; 
  char* text0 =       "TTCACKCTAGT"; //第一处mism考虑用iupac, G->K(T,G)
  
  char* pattern =(char*)malloc(strlen(pattern0));
  char* text =(char*)malloc(strlen(text0));
  for(int i=0; i<strlen(pattern0); i++){
      pattern[i] =grayVal[nt16_table[pattern0[i]]]; //ASCII -> iupac/下标 -> gray
  }
  for(int i=0; i<strlen(text0); i++){
      text[i] =grayVal[nt16_table[text0[i]]];
  }

  // Configure alignment attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = 0;
  attributes.affine_penalties.mismatch = 4;
  attributes.affine_penalties.gap_opening = 6;
  attributes.affine_penalties.gap_extension = 2;

  attributes.alignment_form.span = alignment_endsfree;
  attributes.alignment_form.pattern_begin_free = 0;
  attributes.alignment_form.pattern_end_free = 0;
  attributes.alignment_form.text_begin_free = strlen(text); //开头或结尾允许free的最大gap长度
  attributes.alignment_form.text_end_free = strlen(text);

  // Initialize Wavefront Aligner
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  // Align
  wavefront_align(wf_aligner,pattern,strlen(pattern),text,strlen(text));
  fprintf(stderr,"WFA-Alignment returns score %d\n",wf_aligner->cigar->score); //正确的score
  // Display alignment
  fprintf(stderr,"  PATTERN  %s\n",pattern);
  fprintf(stderr,"  TEXT     %s\n",text);
  fprintf(stderr,"  SCORE (RE)COMPUTED %d\n",
      cigar_score_gap_affine(wf_aligner->cigar,&attributes.affine_penalties)); //按end-to-end得到的score
  cigar_print_pretty(stderr,wf_aligner->cigar,pattern,strlen(pattern),text,strlen(text));
  // Free
  wavefront_aligner_delete(wf_aligner);
}
