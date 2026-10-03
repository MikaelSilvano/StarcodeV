/*
** autocal.h -- automatic tau selection for StarcodeV, without ground truth.
**
** A fixed tau is tuned to one dataset. The consensus-to-consensus distance
** histogram is bimodal (pairs of the same oligo close together, pairs of
** different oligos far apart), and the empty gap between the two modes gives
** a safe ceiling for tau.
**
** Gap rule (cal_ceiling_below_peak): find the histogram PEAK, which lies in
** the different-oligo mode because unrelated consensus pairs far outnumber
** same-oligo pairs, and walk LEFT from it to the first empty run of at least
** CAL_MIN_GAP bins. The ceiling is the first occupied bin to the right of
** that run. If the walk reaches d = 1 without meeting such a run, the data
** show no same-oligo mode and the ceiling is the smallest observed distance.
** Looking only to the left of the peak matters on real data: the sparse
** right tail of the Microsoft histogram has an empty run at d = 43..47 that
** is longer than the true gap at 21..23, and a single empty bin inside a
** mode is not a gap.
**
** No RNG and no floating-point accumulation on the decision path (the only
** double is the margin, applied once right before rounding to an int), no
** dependency on input order.
*/
#ifndef _AUTOCAL_HEADER
#define _AUTOCAL_HEADER

#include <stdint.h>
#include "starcodev.h"

#define CAL_MARGIN  0.25  /* safety margin below the ceiling: tau = floor((ceiling-1)*(1-margin)) */
#define CAL_MIN_GAP 2     /* an empty run must span >= this many bins to count as a separability gap */

/* Pairwise edit-distance histogram over candidate consensus sequences
** (integer, deterministic). hist_out must be allocated by the caller with
** capacity >= dmax+2 -- call this once with hist_out=NULL to get dmax_out,
** then allocate, or use cal_estimate_dmax() to get dmax up front.
** Returns 0 when ncons < 2 (no histogram). */
int32_t cal_estimate_dmax(const int32_t *cons_lens, int32_t ncons, double frac_dmax);

int cal_distance_histogram(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                            int32_t n_iter, int32_t sig_abs, int32_t dmax,
                            int64_t *hist_out /* size dmax+2 */);


/* The gap rule described at the top. Finds the peak bin (largest count, smallest d on ties),
** then walks left to the first empty run of length >= min_gap. Fills
** peak_out, gap_start_out/gap_len_out (0/0 when the leading-edge case
** applies) and ceiling_out. kind_out: 1 = interior gap below the peak,
** 2 = no such gap (ceiling = smallest observed distance). Returns 0 only if
** the histogram is empty. */
int cal_ceiling_below_peak(const int64_t *hist, int32_t dmax, int32_t min_gap,
                            int32_t *peak_out, int32_t *gap_start_out, int32_t *gap_len_out,
                            int32_t *ceiling_out, int *kind_out);

typedef struct {
    int64_t *hist;      /* allocated; the caller frees it */
    int32_t  dmax;
    int32_t  gap_start;
    int32_t  gap_len;
    int32_t  ceiling;
    int32_t  peak;      /* histogram peak */
    int      kind;      /* 1 = interior gap below peak, 2 = leading edge (no same-oligo mode) */
    int64_t  n_pairs;   /* candidate pairs counted in the histogram */
    int      has_gap;  /* 0 means "no gap found" */
} cal_diagnostics_t;

/* Pick tau for S3/S4 from the data alone, using a histogram of filtered
** candidate pairs (--approx mode; the default mode builds the exact
** histogram in exactcons.c instead). Returns 1 on success (tau_out and
** ceiling_out are filled), 0 if ncons < 2 or no gap was found (tau_out and
** ceiling_out are set to -1). diag_out is optional (may be NULL); if it is
** provided, the caller MUST free(diag_out->hist). */
int cal_auto_tau(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                  double margin, int32_t n_iter, int32_t sig_abs,
                  int32_t *tau_out, int32_t *ceiling_out, cal_diagnostics_t *diag_out);

typedef struct {
    int32_t core_count;
    int32_t orphan_count;
    double  frac_reads_in_cores;
} cal_coverage_t;

/* Diagnostic after S1, no ground truth needed: fraction of DISTINCT reads
** that ended up inside a core (component with >= 2 distinct reads) rather
** than as an orphan. Printed by starcodev_seq in verbose mode. */
cal_coverage_t cal_core_coverage(const int32_t *labels_s1, int32_t n);

#endif /* _AUTOCAL_HEADER */
