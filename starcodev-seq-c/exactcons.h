/*
** exactcons.h -- exact work on the consensus set.
**
** S1 has to use a filtered candidate search because the read set is large.
** Everything after S2 works on one consensus per core, and the number of
** cores is close to the number of oligos (10,066 on the Microsoft data,
** against 267,558 distinct reads). That set is small enough to handle
** exactly, so every decision that could join two oligos uses exact
** distances:
**
**   H   all consensus pairs -> exact distance histogram -> gap -> ceiling C
**   S3  union-find over every pair with d <= tau_cons (no filter, no misses)
**   S4  an orphan joins a component only if that component holds its
**       nearest consensus, checked exactly when it can't be proven cheaply
**
** When something can't be settled, the read stays on its own. A split costs
** one extra cluster; a wrong merge corrupts a cluster.
**
** All decisions are integer min/sum reductions, so the order in which pairs
** are visited never changes the result, which also holds for a parallel
** implementation where that order is not fixed.
*/
#ifndef _EXACTCONS_HEADER
#define _EXACTCONS_HEADER

#include <stdint.h>

/* ---------------------------------------------------------------- distance
** Bit-parallel Levenshtein (Myers 1999, block version). Exact global
** distance, no cutoff. The pattern side is pre-processed once into match
** masks so it can be reused against many texts. */
typedef struct {
    int32_t   len;     /* pattern length */
    int32_t   words;   /* 64-bit blocks */
    uint64_t *eq;      /* 4 * words masks, one row per base A,C,G,T */
    const char *seq;   /* kept for the rare non-ACGT text character */
} ec_peq_t;

void    ec_peq_build(const char *s, int32_t len, ec_peq_t *p);
void    ec_peq_free(ec_peq_t *p);
int32_t ec_dist(const ec_peq_t *p, const char *t, int32_t n);

/* ---------------------------------------------------------------- all pairs
** hist[d] for d in [0, dmax]: number of consensus pairs at exact distance d.
** Pairs with d <= keep are also stored as edges (i < j). */
typedef struct {
    int32_t *a, *b, *d;
    int64_t  count, cap;
} ec_edges_t;

void ec_edges_free(ec_edges_t *e);

void ec_all_pairs(const ec_peq_t *peq, const char *const *cons, const int32_t *lens,
                  int32_t k, int32_t dmax, int32_t keep,
                  int64_t *hist /* dmax+1 */, ec_edges_t *edges /* may be NULL */);

/* ---------------------------------------------------------------- S3
** Union-find over the stored edges with d <= tau. lab_out: size k.
** Returns the number of successful unions. c_post_out gets the smallest
** distance between two consensus left in different components, or keep+1
** when no such pair was stored (a lower bound, which is all S4 needs). */
int64_t ec_s3_merge(const ec_edges_t *e, int32_t k, int32_t tau, int32_t keep,
                 int32_t *lab_out, int32_t *c_post_out);

/* ---------------------------------------------------------------- S4
** Input: the filtered S4 result (bd/bi per orphan, bi = -1 when nothing
** was found within tau_abs). Output: final consensus index per orphan, or
** -1 when it stays a singleton.
**
** Rule per orphan, with a = distance to the filtered candidate:
**   2a < c_post      accept. Any consensus c' in another component has
**                    d(r,c') >= d(c,c') - a >= c_post - a > a, so the
**                    candidate's component holds the nearest consensus
**                    even if the filter never proposed it.
**   otherwise        find the nearest consensus exactly. With a candidate,
**                    only its neighbours within a + tau_abs can be closer
**                    than tau_abs (triangle inequality), so the edge list
**                    is enough. Without a candidate, scan all of them.
**                    Accept only if d <= tau_abs and no consensus from a
**                    different component is equally close. */
typedef struct {
    int64_t fast;          /* accepted by the 2a < c_post rule */
    int64_t local;         /* checked against the candidate's neighbours */
    int64_t full;          /* no candidate: full scan */
    int64_t moved;         /* exact nearest is in a different component than the filter's pick */
    int64_t added;         /* absorbed now, filter had found nothing */
    int64_t tie_rejected;  /* left alone because two components were equally close */
    int64_t dist_calls;
} ec_s4_stats_t;

void ec_s4_resolve(const char *const *orph, const int32_t *orph_lens, int32_t norph,
                   const ec_peq_t *peq, const char *const *cons, const int32_t *cons_lens, int32_t k,
                   const int32_t *lab3, const ec_edges_t *edges, int32_t keep,
                   int32_t c_post, int32_t tau_abs,
                   const int32_t *bd_in, const int32_t *bi_in,
                   int32_t *bi_out, uint8_t *absorbed_out, ec_s4_stats_t *st);

/* ---------------------------------------------------------------- S1 audit
** Exact check of S1 on a deterministic sample of reads (every n/sample-th
** read). For each sampled read, find every read within tau_core with an
** exact search (pigeonhole filter, then exact distance) and count those
** that S1 put in a different component. No ground truth involved. */
typedef struct {
    int32_t sampled;
    int64_t neighbours;       /* true d <= tau_core pairs found from the sample */
    int64_t cross;            /* of those, endpoints in different S1 components */
    int32_t reads_with_cross; /* sampled reads with at least one such pair */
} ec_audit_t;

void ec_s1_audit(const char *const *seqs, const int32_t *lens, int32_t n,
                 const int32_t *lab1, int32_t tau_core, int32_t sample, ec_audit_t *out);

#endif
