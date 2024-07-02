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
 * DESCRIPTION: WaveFront aligner external C-ABI functions with accompanying data structures
 */

#pragma once

#include "../system/mm_allocator.h"
#include "wavefront_slab.h"
#include "wavefront_penalties.h"
#include "wavefront_attributes.h"
#include "wavefront_components.h"
#include "wavefront_sequences.h"
#include "wavefront_bialigner.h"

/*
 * Error codes & messages
 */
// [OK]
#define WF_STATUS_ALG_COMPLETED            0  // Success (Complete alignment found)
#define WF_STATUS_ALG_PARTIAL              1  // Success (Partial alignment found)
// [FAIL]
#define WF_STATUS_MAX_STEPS_REACHED     -100  // Maximum number of WFA-steps reached
#define WF_STATUS_OOM                   -200  // Maximum memory limit reached
#define WF_STATUS_UNATTAINABLE          -300  // Alignment unattainable under configured heuristics
// [INTERNAL]
#define WF_STATUS_OK                      -1  // Computing alignment (in progress)
#define WF_STATUS_END_REACHED             -2  // Alignment end reached
#define WF_STATUS_END_UNREACHABLE         -3  // Alignment end unreachable under current configuration (eg Z-drop)
// Error messages
char* wavefront_align_strerror(const int error_code);
char* wavefront_align_strerror_short(const int error_code);

/*
 * Alignment status
 */
typedef struct _wavefront_aligner_t wavefront_aligner_t;
typedef struct {
  // Status
  int status;                                                     // Status code
  int score;                                                      // Current WF-alignment score
  bool dropped;                                                   // Heuristically dropped
  int num_null_steps;                                             // Total contiguous null-steps performed
  uint64_t memory_used;                                           // Total memory used
  // Wavefront alignment functions
  void (*wf_align_compute)(wavefront_aligner_t* const,const int); // WF Compute function
  int (*wf_align_extend)(wavefront_aligner_t* const,const int);   // WF Extend function
} wavefront_align_status_t;

/*
 * Alignment type
 */
typedef enum {
  wf_align_regular = 0,
  wf_align_biwfa = 1,
  wf_align_biwfa_breakpoint_forward = 2,
  wf_align_biwfa_breakpoint_reverse = 3,
  wf_align_biwfa_subsidiary = 4
} wavefront_align_mode_t;

/*
 * Wavefront Aligner
 */
typedef struct _wavefront_aligner_t {
  // Mode and status
  wavefront_align_mode_t align_mode;          // WFA alignment mode
  char* align_mode_tag;                       // WFA mode tag
  wavefront_align_status_t align_status;      // Current alignment status
  // Sequences
  wavefront_sequences_t sequences;            // Input sequences
  // Alignment Attributes
  alignment_scope_t alignment_scope;          // Alignment scope (score only or full-CIGAR)
  alignment_form_t alignment_form;            // Alignment form (end-to-end/ends-free)
  wavefront_penalties_t penalties;            // Alignment penalties
  wavefront_heuristic_t heuristic;            // Heuristic's parameters
  wavefront_memory_t memory_mode;             // Wavefront memory strategy (modular wavefronts and piggyback)
  // Wavefront components
  wavefront_components_t wf_components;       // Wavefront components
  affine2p_matrix_type component_begin;       // Alignment begin component
  affine2p_matrix_type component_end;         // Alignment end component
  wavefront_pos_t alignment_end_pos;          // Alignment end position
  // Bidirectional Alignment
  wavefront_bialigner_t* bialigner;           // BiWFA aligner
  // CIGAR
  cigar_t* cigar;                             // Alignment CIGAR
  // MM
  bool mm_allocator_own;                      // Ownership of MM-Allocator
  mm_allocator_t* mm_allocator;               // MM-Allocator
  wavefront_slab_t* wavefront_slab;           // MM-Wavefront-Slab (Allocates/Reuses the individual wavefronts)
  // Display
  wavefront_plot_t* plot;                     // Wavefront plot
  // System
  alignment_system_t system;                  // System related parameters
} wavefront_aligner_t;

/*
 * Setup
 */
wavefront_aligner_t* wavefront_aligner_new(
    wavefront_aligner_attr_t* attributes);
void wavefront_aligner_reap(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_delete(
    wavefront_aligner_t* const wf_aligner);

/*
 * Span configuration
 */
void wavefront_aligner_set_alignment_end_to_end(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_set_alignment_free_ends(
    wavefront_aligner_t* const wf_aligner,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free);
void wavefront_aligner_set_alignment_extension(
    wavefront_aligner_t* const wf_aligner);

/*
 * Heuristic configuration
 */
void wavefront_aligner_set_heuristic_none(
    wavefront_aligner_t* const wf_aligner);
void wavefront_aligner_set_heuristic_wfadaptive(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int score_steps);
void wavefront_aligner_set_heuristic_wfmash(
    wavefront_aligner_t* const wf_aligner,
    const int min_wavefront_length,
    const int max_distance_threshold,
    const int score_steps);
void wavefront_aligner_set_heuristic_xdrop(
    wavefront_aligner_t* const wf_aligner,
    const int xdrop,
    const int score_steps);
void wavefront_aligner_set_heuristic_zdrop(
    wavefront_aligner_t* const wf_aligner,
    const int ydrop,
    const int score_steps);
void wavefront_aligner_set_heuristic_banded_static(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k);
void wavefront_aligner_set_heuristic_banded_adaptive(
    wavefront_aligner_t* const wf_aligner,
    const int band_min_k,
    const int band_max_k,
    const int score_steps);

/*
 * System configuration
 */
void wavefront_aligner_set_max_alignment_steps(
    wavefront_aligner_t* const wf_aligner,
    const int max_alignment_steps);
void wavefront_aligner_set_max_memory(
    wavefront_aligner_t* const wf_aligner,
    const uint64_t max_memory_resident,
    const uint64_t max_memory_abort);
void wavefront_aligner_set_max_num_threads(
    wavefront_aligner_t* const wf_aligner,
    const int max_num_threads);
void wavefront_aligner_set_min_offsets_per_thread(
    wavefront_aligner_t* const wf_aligner,
    const int min_offsets_per_thread);

/*
 * Wavefront Align
 */
int wavefront_align(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
int wavefront_align_lambda(
    wavefront_aligner_t* const wf_aligner,
    alignment_match_funct_t const match_funct,
    void* match_funct_arguments,
    const int pattern_length,
    const int text_length);
int wavefront_align_packed2bits(
    wavefront_aligner_t* const wf_aligner,
    const uint8_t* const pattern,
    const int pattern_length,
    const uint8_t* const text,
    const int text_length);
