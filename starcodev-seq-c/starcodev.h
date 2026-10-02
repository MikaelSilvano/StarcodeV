/*
** starcodev.h -- StarcodeV Stage 1-4, union-find, canonical partition.
**
** Structurally equivalent to src/starcode.c/.h from the original Starcode
** (gui11aume/starcode v1.4) -- see CONVENTIONS_STARCODEV.md. Candidate-pair generation is delegated
** to pairgen.c/h (our replacement for trie.c/h).
**
** The restrictions from §4.2 of the StarcodeV design doc apply across this
** whole module: no runtime RNG, no floating-point reduction on any
** decision path, no order-dependent accumulators, no dependency on input
** ordering. The canonical tie-breakers K1-K5 are documented next to each
** function that uses them.
*/
#ifndef _STARCODEV_HEADER
#define _STARCODEV_HEADER

#include <stdint.h>
#include <stddef.h>

/* ---------------------------------------------------------------- constants
** Values were calibrated on the Microsoft clustered-nanopore data (see
** VALIDATION_v2.md and AUTOTAU_REVISION.md). */
#define SV_TAU_CORE   8     /* Starcode's own semantics, fixed (trie.h:52) */
#define SV_TAU_ABS    16    /* orphan-absorption threshold (S3), default; tau_auto() can override */
#define SV_TAU_CONS   20    /* consensus-merge threshold (S4), default; tau_auto() can override */
#define SV_SIG_ABS    96    /* Hamming pre-filter threshold for S3/S4 */
#define SV_THETA_LOW  40
#define SV_THETA_HIGH 60    /* Hamming pre-filter threshold for S1 */
#define SV_N_ITER_S1  64    /* calibrated: 100% edge recall on the 300-block subsample */
#define SV_N_ITER_S3  96
#define SV_N_ITER_S4  96
#define SV_M_MAX      64    /* medoid cost cap (§3.2) */

/* ---------------------------------------------------------------- union-find */
typedef struct {
    int32_t *p;
    int32_t  n;
} sv_uf_t;

void sv_uf_init(sv_uf_t *uf, int32_t n);
int32_t sv_uf_find(sv_uf_t *uf, int32_t x);
/* Returns 1 if the two components were actually merged (previously distinct). */
int sv_uf_union(sv_uf_t *uf, int32_t a, int32_t b);
/* out must be allocated by the caller, sized uf->n. */
void sv_uf_labels(sv_uf_t *uf, int32_t *out);
void sv_uf_free(sv_uf_t *uf);

/* ---------------------------------------------------------------- K5: canonical partition
** Final cluster ID = lexicographic rank of the smallest member string in
** each cluster. Independent of union-find merge order. out must be
** allocated by the caller (size n). */
void sv_canonical_partition(const int32_t *labels, const char *const *seqs,
                             int32_t n, int32_t *out);

/* Number of canonical labels that differ between two partitions. WARNING:
** this is NOT a measure of error magnitude -- a single missing edge can
** cascade into thousands of "different label" positions purely from an ID
** shift. Use it ONLY as a binary gate (==0). */
int64_t sv_partitions_equal(const int32_t *la, const int32_t *lb,
                             const char *const *seqs, int32_t n);

/* ---------------------------------------------------------------- work-item counters */
typedef struct {
    int64_t cand;
    int64_t sig_pass;
    int64_t edit;
    int64_t uni;      /* number of successful unions/absorptions */
} sv_work_t;

/* ---------------------------------------------------------------- STAGE 1: exact cores
** lab_out: size n, raw (NOT canonical) union-find label per read.
** anchor_occ: keys per read per iteration (1 = v1 behaviour, see pairgen.h).
** iters_done_out: number of iterations actually run. */
void sv_s1_cores(const char *const *seqs, const int32_t *lens, int32_t n,
                  int32_t tau_core, int32_t n_iter, int32_t theta, int32_t anchor_occ,
                  int32_t *lab_out, sv_work_t *work_out, int32_t *iters_done_out);

/* List of cores (components with >1 member) and orphans (singletons), both
** sorted ascending. */
typedef struct {
    int32_t **members;   /* array of arrays */
    int32_t  *sizes;
    int32_t   count;
} sv_core_list_t;

void sv_split_cores_orphans(const int32_t *lab, int32_t n,
                             sv_core_list_t *cores_out,
                             int32_t **orphans_out, int32_t *n_orph_out);
void sv_core_list_free(sv_core_list_t *cl);

/* ---------------------------------------------------------------- STAGE 2: consensus
** cons_out: array of char* (allocated; caller frees each element plus the array). */
typedef struct {
    char   **cons;
    int32_t *sizes;
    int32_t *anchors;
    int32_t  count;
    int32_t  tie_medoid;
    int32_t  tie_column;
} sv_s2_result_t;

void sv_s2_consensus(const sv_core_list_t *cores, const char *const *seqs,
                      sv_s2_result_t *out);
void sv_s2_result_free(sv_s2_result_t *r);

/* ---------------------------------------------------------------- STAGE 3: orphan absorption
** bd_out/bi_out/absorbed_out: size norph. bi_out holds the (local, 0-based)
** consensus index, or -1 if not absorbed. */
void sv_s3_absorb(const char *const *orph_seqs, const int32_t *orph_lens, int32_t norph,
                   const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                   int32_t tau_abs, int32_t sig_abs, int32_t n_iter, int32_t anchor_occ,
                   int32_t *bd_out, int32_t *bi_out, uint8_t *absorbed_out,
                   sv_work_t *work_out);

/* ---------------------------------------------------------------- STAGE 4: consensus merge
** lab_out: size ncons, raw (not yet canonical) union-find label. */
void sv_s4_merge_cons(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                       int32_t tau_cons, int32_t sig_abs, int32_t n_iter,
                       int32_t *lab_out, sv_work_t *work_out);

/* ---------------------------------------------------------------- final assembly
** out: size n_reads, final per-read label (not yet canonical -- call
** sv_canonical_partition() if you need the canonical form). */
void sv_assemble(const sv_core_list_t *cores,
                  const int32_t *orph_idx, int32_t norph,
                  const int32_t *bi3, const uint8_t *absorbed,
                  const int32_t *lab_s4, int32_t n_reads,
                  int32_t *out);

/* ---------------------------------------------------------------- Levenshtein helper
** Also used by autocal.c. score_cutoff exists purely as a banded early-exit
** optimization: the return value is ALWAYS the true edit distance when it's
** <= cutoff, and may be any value > cutoff once the cutoff is exceeded (the
** caller may only test "<= cutoff", the usual contract of a
** score cutoff). */
int32_t sv_levenshtein(const char *a, int32_t la, const char *b, int32_t lb,
                        int32_t score_cutoff);

#endif /* _STARCODEV_HEADER */
