/*
** exactcons.c -- see exactcons.h for what this module is for and why.
*/
#include "exactcons.h"
#include "starcodev.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline int base_code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

/* ================================================================ distance */

void ec_peq_build(const char *s, int32_t len, ec_peq_t *p) {
    p->len = len;
    p->words = len > 0 ? (len + 63) / 64 : 1;
    p->seq = s;
    p->eq = (uint64_t *)calloc((size_t)4 * p->words, sizeof(uint64_t));
    for (int32_t i = 0; i < len; i++) {
        int b = base_code(s[i]);
        if (b >= 0) p->eq[b * p->words + i / 64] |= 1ULL << (i % 64);
    }
}

void ec_peq_free(ec_peq_t *p) {
    free(p->eq);
    p->eq = NULL;
}

/* Myers' advance_block, run over every block for each text character.
** hin/hout are the horizontal deltas entering/leaving a block; the top row
** of the DP is 0,1,2,... so block 0 always gets +1. The score tracks the
** bottom-right cell, starting at m (empty text against the full pattern). */
int32_t ec_dist(const ec_peq_t *p, const char *t, int32_t n) {
    int32_t m = p->len;
    if (m == 0) return n;
    if (n == 0) return m;
    int32_t W = p->words;

    uint64_t pv_stack[16], mv_stack[16], tmp_stack[16];
    uint64_t *Pv = pv_stack, *Mv = mv_stack, *tmp = tmp_stack;
    if (W > 16) {
        Pv = (uint64_t *)malloc((size_t)W * sizeof(uint64_t));
        Mv = (uint64_t *)malloc((size_t)W * sizeof(uint64_t));
        tmp = (uint64_t *)malloc((size_t)W * sizeof(uint64_t));
    }
    for (int32_t w = 0; w < W; w++) { Pv[w] = ~0ULL; Mv[w] = 0; }
    const uint64_t last = 1ULL << ((m - 1) % 64);
    int32_t score = m;

    for (int32_t j = 0; j < n; j++) {
        const uint64_t *eq;
        int b = base_code(t[j]);
        if (b >= 0) {
            eq = p->eq + (size_t)b * W;
        } else {
            /* not A/C/G/T: build the mask by hand (rare, keeps the result exact) */
            memset(tmp, 0, (size_t)W * sizeof(uint64_t));
            for (int32_t i = 0; i < m; i++)
                if (p->seq[i] == t[j]) tmp[i / 64] |= 1ULL << (i % 64);
            eq = tmp;
        }
        int hin = 1;
        for (int32_t w = 0; w < W; w++) {
            uint64_t Eq = eq[w], P = Pv[w], M = Mv[w];
            uint64_t Xv = Eq | M;
            if (hin < 0) Eq |= 1;
            uint64_t Xh = (((Eq & P) + P) ^ P) | Eq;
            uint64_t Ph = M | ~(Xh | P);
            uint64_t Mh = P & Xh;
            uint64_t high = (w == W - 1) ? last : (1ULL << 63);
            int hout = 0;
            if (Ph & high) hout = 1;
            else if (Mh & high) hout = -1;
            Ph <<= 1;
            Mh <<= 1;
            if (hin < 0) Mh |= 1;
            else if (hin > 0) Ph |= 1;
            Pv[w] = Mh | ~(Xv | Ph);
            Mv[w] = Ph & Xv;
            hin = hout;
        }
        score += hin;
    }
    if (W > 16) { free(Pv); free(Mv); free(tmp); }
    return score;
}

/* ================================================================ edges */

static void edges_push(ec_edges_t *e, int32_t a, int32_t b, int32_t d) {
    if (e->count == e->cap) {
        e->cap = e->cap ? e->cap * 2 : 1024;
        e->a = (int32_t *)realloc(e->a, (size_t)e->cap * sizeof(int32_t));
        e->b = (int32_t *)realloc(e->b, (size_t)e->cap * sizeof(int32_t));
        e->d = (int32_t *)realloc(e->d, (size_t)e->cap * sizeof(int32_t));
    }
    e->a[e->count] = a; e->b[e->count] = b; e->d[e->count] = d;
    e->count++;
}

void ec_edges_free(ec_edges_t *e) {
    free(e->a); free(e->b); free(e->d);
    e->a = e->b = e->d = NULL;
    e->count = e->cap = 0;
}

void ec_all_pairs(const ec_peq_t *peq, const char *const *cons, const int32_t *lens,
                  int32_t k, int32_t dmax, int32_t keep,
                  int64_t *hist, ec_edges_t *edges) {
    memset(hist, 0, (size_t)(dmax + 1) * sizeof(int64_t));
    for (int32_t i = 0; i < k; i++) {
        for (int32_t j = i + 1; j < k; j++) {
            int32_t dl = lens[i] - lens[j];
            if (dl < 0) dl = -dl;
            if (dl > dmax) continue;           /* distance >= length gap */
            int32_t d = ec_dist(&peq[i], cons[j], lens[j]);
            if (d > dmax) continue;
            hist[d]++;
            if (edges && d <= keep) edges_push(edges, i, j, d);
        }
    }
}

/* ================================================================ S3 */

int64_t ec_s3_merge(const ec_edges_t *e, int32_t k, int32_t tau, int32_t keep,
                 int32_t *lab_out, int32_t *c_post_out) {
    sv_uf_t uf;
    sv_uf_init(&uf, k);
    int64_t unions = 0;
    for (int64_t x = 0; x < e->count; x++)
        if (e->d[x] <= tau && sv_uf_union(&uf, e->a[x], e->b[x])) unions++;

    int32_t c_post = keep + 1;
    for (int64_t x = 0; x < e->count; x++)
        if (e->d[x] < c_post && sv_uf_find(&uf, e->a[x]) != sv_uf_find(&uf, e->b[x]))
            c_post = e->d[x];

    sv_uf_labels(&uf, lab_out);
    sv_uf_free(&uf);
    *c_post_out = c_post;
    return unions;
}

/* ================================================================ S4 */

/* neighbour lists in CSR form, built from the undirected edge list */
typedef struct { int64_t *start; int32_t *nb; int32_t *d; } csr_t;

static void csr_build(const ec_edges_t *e, int32_t k, csr_t *g) {
    g->start = (int64_t *)calloc((size_t)k + 1, sizeof(int64_t));
    for (int64_t x = 0; x < e->count; x++) { g->start[e->a[x] + 1]++; g->start[e->b[x] + 1]++; }
    for (int32_t i = 0; i < k; i++) g->start[i + 1] += g->start[i];
    int64_t tot = g->start[k];
    g->nb = (int32_t *)malloc((size_t)(tot > 0 ? tot : 1) * sizeof(int32_t));
    g->d  = (int32_t *)malloc((size_t)(tot > 0 ? tot : 1) * sizeof(int32_t));
    int64_t *fill = (int64_t *)malloc((size_t)k * sizeof(int64_t));
    memcpy(fill, g->start, (size_t)k * sizeof(int64_t));
    for (int64_t x = 0; x < e->count; x++) {
        int32_t a = e->a[x], b = e->b[x];
        g->nb[fill[a]] = b; g->d[fill[a]++] = e->d[x];
        g->nb[fill[b]] = a; g->d[fill[b]++] = e->d[x];
    }
    free(fill);
}

static void csr_free(csr_t *g) { free(g->start); free(g->nb); free(g->d); }

void ec_s4_resolve(const char *const *orph, const int32_t *orph_lens, int32_t norph,
                   const ec_peq_t *peq, const char *const *cons, const int32_t *cons_lens, int32_t k,
                   const int32_t *lab3, const ec_edges_t *edges, int32_t keep,
                   int32_t c_post, int32_t tau_abs,
                   const int32_t *bd_in, const int32_t *bi_in,
                   int32_t *bi_out, uint8_t *absorbed_out, ec_s4_stats_t *st) {
    (void)cons;
    memset(st, 0, sizeof(*st));
    csr_t g;
    csr_build(edges, k, &g);
    int32_t *cand = (int32_t *)malloc((size_t)(k > 0 ? k : 1) * sizeof(int32_t));
    int32_t *cd   = (int32_t *)malloc((size_t)(k > 0 ? k : 1) * sizeof(int32_t));

    for (int32_t o = 0; o < norph; o++) {
        int32_t bi = bi_in[o], a = bd_in[o];
        bi_out[o] = -1;
        absorbed_out[o] = 0;

        if (bi >= 0 && 2 * a < c_post) {
            bi_out[o] = bi;
            absorbed_out[o] = 1;
            st->fast++;
            continue;
        }

        /* collect the candidates whose exact distance we need */
        int32_t nc = 0;
        int32_t len = orph_lens[o];
        if (bi >= 0 && a + tau_abs <= keep) {
            st->local++;
            cand[nc] = bi; cd[nc] = a; nc++;
            for (int64_t x = g.start[bi]; x < g.start[bi + 1]; x++) {
                if (g.d[x] > a + tau_abs) continue;
                int32_t c = g.nb[x];
                int32_t dl = cons_lens[c] - len; if (dl < 0) dl = -dl;
                if (dl > tau_abs) continue;
                cand[nc] = c; cd[nc] = ec_dist(&peq[c], orph[o], len); nc++;
                st->dist_calls++;
            }
        } else {
            st->full++;
            for (int32_t c = 0; c < k; c++) {
                int32_t dl = cons_lens[c] - len; if (dl < 0) dl = -dl;
                if (dl > tau_abs) continue;
                cand[nc] = c; cd[nc] = ec_dist(&peq[c], orph[o], len); nc++;
                st->dist_calls++;
            }
        }

        /* nearest consensus, ties to the smaller index */
        int32_t best = -1, best_d = INT_MAX;
        for (int32_t x = 0; x < nc; x++)
            if (cd[x] < best_d || (cd[x] == best_d && cand[x] < best)) { best_d = cd[x]; best = cand[x]; }
        if (best < 0 || best_d > tau_abs) continue;

        /* someone from another component just as close: don't pick a side */
        int tie = 0;
        for (int32_t x = 0; x < nc; x++)
            if (cd[x] == best_d && lab3[cand[x]] != lab3[best]) { tie = 1; break; }
        if (tie) { st->tie_rejected++; continue; }

        bi_out[o] = best;
        absorbed_out[o] = 1;
        if (bi < 0) st->added++;
        else if (lab3[bi] != lab3[best]) st->moved++;
    }
    free(cand); free(cd);
    csr_free(&g);
}

/* ================================================================ S1 audit */

/* With at most tau edits, cutting a into tau+1 pieces leaves at least one
** piece untouched, and that piece shows up in b shifted by at most tau.
** So if no piece matches, d(a,b) > tau for sure. */
static int pigeonhole_pass(const char *a, int32_t la, const char *b, int32_t lb, int32_t tau) {
    int32_t pieces = tau + 1;
    if (la < pieces) return 1;
    for (int32_t t = 0; t < pieces; t++) {
        int32_t s = (int32_t)((int64_t)t * la / pieces);
        int32_t e = (int32_t)((int64_t)(t + 1) * la / pieces);
        int32_t L = e - s;
        for (int32_t delta = -tau; delta <= tau; delta++) {
            int32_t q = s + delta;
            if (q < 0 || q + L > lb) continue;
            if (a[s] == b[q] && memcmp(a + s, b + q, (size_t)L) == 0) return 1;
        }
    }
    return 0;
}

void ec_s1_audit(const char *const *seqs, const int32_t *lens, int32_t n,
                 const int32_t *lab1, int32_t tau_core, int32_t sample, ec_audit_t *out) {
    memset(out, 0, sizeof(*out));
    if (sample <= 0 || n < 2) return;
    if (sample > n) sample = n;
    out->sampled = sample;
    for (int32_t s = 0; s < sample; s++) {
        int32_t i = (int32_t)((int64_t)s * n / sample);
        ec_peq_t p;
        ec_peq_build(seqs[i], lens[i], &p);
        int flagged = 0;
        for (int32_t j = 0; j < n; j++) {
            if (j == i) continue;
            int32_t dl = lens[i] - lens[j]; if (dl < 0) dl = -dl;
            if (dl > tau_core) continue;
            if (!pigeonhole_pass(seqs[i], lens[i], seqs[j], lens[j], tau_core)) continue;
            if (ec_dist(&p, seqs[j], lens[j]) > tau_core) continue;
            out->neighbours++;
            if (lab1[i] != lab1[j]) { out->cross++; flagged = 1; }
        }
        if (flagged) out->reads_with_cross++;
        ec_peq_free(&p);
    }
}
