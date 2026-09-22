/*
** autocal.h -- automatic tau selection for StarcodeV, WITHOUT ground truth.
**
** C mirror of autocal.py (the validated reference; see
** CORRECTNESS_REPORT_GENERALIZATION.md: 100% purity on 11/11 synthetic
** regimes, passes the 3-permutation determinism test). See the long
** comment block in autocal.py for the full story of the problem (a fixed
** tau overfits) and the solution (the bimodal separability gap in the
** consensus-to-consensus distance histogram).
**
** DETERMINISM (identical to the Python version): no RNG, no
** floating-point reduction anywhere on the decision path (the only double
** is MARGIN, used right before it gets rounded to an int), no dependency
** on input order. The gap is chosen as the LONGEST empty gap, with the
** smallest starting index as the tie-breaker.
*/
#ifndef _AUTOCAL_HEADER
#define _AUTOCAL_HEADER

#include <stdint.h>
#include "starcodev.h"

#define CAL_MARGIN 0.25   /* safety margin below the ceiling */

/* Pairwise edit-distance histogram over candidate consensus sequences
** (integer, deterministic). hist_out must be allocated by the caller with
** capacity >= dmax+2 -- call this once with hist_out=NULL to get dmax_out,
** then allocate, or use cal_estimate_dmax() to get dmax up front.
** Returns 0 when ncons < 2 (no histogram, matching Python's None). */
int32_t cal_estimate_dmax(const int32_t *cons_lens, int32_t ncons, double frac_dmax);

int cal_distance_histogram(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                            int32_t n_iter, int32_t sig_abs, int32_t dmax,
                            int64_t *hist_out /* size dmax+2 */);

/* Returns 1 if a gap was found (gap_start_out/gap_len_out are filled), 0 otherwise. */
int cal_longest_gap(const int64_t *hist, int32_t dmax,
                     int32_t *gap_start_out, int32_t *gap_len_out);

typedef struct {
    int64_t *hist;      /* allocated; the caller frees it */
    int32_t  dmax;
    int32_t  gap_start;
    int32_t  gap_len;
    int32_t  ceiling;
    int      has_gap;  /* 0 means "no gap found" */
} cal_diagnostics_t;

/* Pick tau for S3/S4 from the data alone. Returns 1 on success (tau_out
** and ceiling_out are filled), 0 if ncons < 2 or no gap was found
** (tau_out/ceiling_out are set to -1 in that case, matching Python's
** (None, None)). diag_out is optional (may be NULL); if it is provided,
** the caller MUST free(diag_out->hist). */
int cal_auto_tau(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                  double margin, int32_t n_iter,
                  int32_t *tau_out, int32_t *ceiling_out, cal_diagnostics_t *diag_out);

typedef struct {
    int32_t core_count;
    int32_t orphan_count;
    double  frac_reads_in_cores;
} cal_coverage_t;

/* Pre-run diagnostic: fraction of reads that ended up inside a core
** (as opposed to an orphan) after S1. */
cal_coverage_t cal_core_coverage(const int32_t *labels_s1, int32_t n);

#endif /* _AUTOCAL_HEADER */
