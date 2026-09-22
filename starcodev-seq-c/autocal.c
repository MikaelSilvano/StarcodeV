/*
** autocal.c -- see autocal.h and the long comment block in autocal.py
** (the reference implementation) for the full problem/solution story.
*/
#include "autocal.h"
#include "pairgen.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* A small open-addressing hash set for deduplicating pairs (i<<32)|j -- a
** copy independent from the one in starcodev.c (that one is static and
** not exported via a header); kept separate so autocal.c doesn't reach
** into starcodev.c's internals. Same semantics: it only ever skips
** redundant work, never affects ordering/decisions (an integer histogram
** is commutative with respect to insertion order). */
typedef struct { uint64_t *keys; uint8_t *used; size_t cap; size_t count; } cal_u64set_t;

static void cal_set_init(cal_u64set_t *s, size_t cap_hint) {
    size_t cap = 16;
    while (cap < cap_hint * 2) cap <<= 1;
    s->keys = (uint64_t *)calloc(cap, sizeof(uint64_t));
    s->used = (uint8_t *)calloc(cap, 1);
    s->cap = cap; s->count = 0;
}
static void cal_set_free(cal_u64set_t *s) {
    free(s->keys); free(s->used); s->keys = NULL; s->used = NULL; s->cap = 0; s->count = 0;
}
static inline size_t cal_hash(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return (size_t)x;
}
static void cal_set_grow(cal_u64set_t *s) {
    size_t newcap = s->cap * 2;
    uint64_t *nk = (uint64_t *)calloc(newcap, sizeof(uint64_t));
    uint8_t *nu = (uint8_t *)calloc(newcap, 1);
    for (size_t i = 0; i < s->cap; i++) {
        if (!s->used[i]) continue;
        size_t h = cal_hash(s->keys[i]) & (newcap - 1);
        while (nu[h]) h = (h + 1) & (newcap - 1);
        nu[h] = 1; nk[h] = s->keys[i];
    }
    free(s->keys); free(s->used);
    s->keys = nk; s->used = nu; s->cap = newcap;
}
static int cal_set_add(cal_u64set_t *s, uint64_t key) {
    if (s->count * 2 >= s->cap) cal_set_grow(s);
    size_t h = cal_hash(key) & (s->cap - 1);
    while (s->used[h]) {
        if (s->keys[h] == key) return 0;
        h = (h + 1) & (s->cap - 1);
    }
    s->used[h] = 1; s->keys[h] = key; s->count++;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Mirrors np.median: odd n -> the middle element; even n -> the average of
** the two middle elements. */
static double median_len(const int32_t *lens, int32_t n) {
    int32_t *tmp = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    memcpy(tmp, lens, (size_t)n * sizeof(int32_t));
    /* insertion sort -- n is usually small (number of cores/consensus sequences) */
    for (int32_t i = 1; i < n; i++) {
        int32_t key = tmp[i], j = i - 1;
        while (j >= 0 && tmp[j] > key) { tmp[j + 1] = tmp[j]; j--; }
        tmp[j + 1] = key;
    }
    int32_t mid = n / 2;
    double result = (n % 2 == 1) ? (double)tmp[mid]
                                  : ((double)tmp[mid - 1] + (double)tmp[mid]) / 2.0;
    free(tmp);
    return result;
}

int32_t cal_estimate_dmax(const int32_t *cons_lens, int32_t ncons, double frac_dmax) {
    if (ncons == 0) return 4;
    double med = median_len(cons_lens, ncons);
    int32_t v = (int32_t)(frac_dmax * med);
    return v > 4 ? v : 4;
}

int cal_distance_histogram(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                            int32_t n_iter, int32_t sig_abs, int32_t dmax,
                            int64_t *hist_out) {
    if (ncons < 2) return 0;
    memset(hist_out, 0, (size_t)(dmax + 2) * sizeof(int64_t));

    pg_encoded_t enc;
    pg_encode(cons, cons_lens, ncons, &enc);
    pg_signatures_t sig;
    pg_signatures(&enc, &sig);

    cal_u64set_t seen;
    cal_set_init(&seen, 4096);
    uint64_t *keys = (uint64_t *)malloc((size_t)ncons * sizeof(uint64_t));

    for (int32_t it = 0; it < n_iter; it++) {
        pg_anchor_bucket(&enc, it, keys);
        pg_pairlist_t pl;
        pg_candidate_pairs(keys, ncons, PG_BUCKET_CAP, &pl);
        for (int64_t k = 0; k < pl.count; k++) {
            int64_t i = pl.a[k], j = pl.b[k];
            int dh = pg_sig_hamming(&sig, i, j);
            if (dh > sig_abs) continue;
            uint64_t key = ((uint64_t)i << 32) | (uint64_t)(uint32_t)j;
            if (!cal_set_add(&seen, key)) continue;
            int32_t d = sv_levenshtein(cons[i], cons_lens[i], cons[j], cons_lens[j], dmax);
            if (d <= dmax) hist_out[d]++;
        }
        pg_pairlist_free(&pl);
    }
    free(keys);
    cal_set_free(&seen);
    pg_free_signatures(&sig);
    pg_free_encoded(&enc);
    return 1;
}

int cal_longest_gap(const int64_t *hist, int32_t dmax,
                     int32_t *gap_start_out, int32_t *gap_len_out) {
    int32_t best_len = 0, best_start = -1;
    int32_t d = 1;
    while (d <= dmax) {
        if (hist[d] == 0) {
            int32_t st = d;
            while (d <= dmax && hist[d] == 0) d++;
            if (d - st > best_len) { best_len = d - st; best_start = st; }
        } else {
            d++;
        }
    }
    if (best_start < 0) return 0;
    *gap_start_out = best_start;
    *gap_len_out = best_len;
    return 1;
}

int cal_auto_tau(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                  double margin, int32_t n_iter,
                  int32_t *tau_out, int32_t *ceiling_out, cal_diagnostics_t *diag_out) {
    if (diag_out) memset(diag_out, 0, sizeof(*diag_out));
    if (ncons < 2) { *tau_out = -1; *ceiling_out = -1; return 0; }

    int32_t dmax = cal_estimate_dmax(cons_lens, ncons, 0.5);
    int64_t *hist = (int64_t *)malloc((size_t)(dmax + 2) * sizeof(int64_t));
    cal_distance_histogram(cons, cons_lens, ncons, n_iter, SV_SIG_ABS, dmax, hist);

    int32_t gap_start, gap_len;
    int found = cal_longest_gap(hist, dmax, &gap_start, &gap_len);
    if (!found) {
        *tau_out = -1; *ceiling_out = -1;
        if (diag_out) { diag_out->hist = hist; diag_out->dmax = dmax; diag_out->has_gap = 0; }
        else free(hist);
        return 0;
    }
    int32_t ceiling = gap_start + gap_len;
    int32_t tau = (int32_t)((double)(ceiling - 1) * (1.0 - margin));
    if (tau < 2) tau = 2;

    *tau_out = tau; *ceiling_out = ceiling;
    if (diag_out) {
        diag_out->hist = hist; diag_out->dmax = dmax;
        diag_out->gap_start = gap_start; diag_out->gap_len = gap_len;
        diag_out->ceiling = ceiling; diag_out->has_gap = 1;
    } else {
        free(hist);
    }
    return 1;
}

cal_coverage_t cal_core_coverage(const int32_t *labels_s1, int32_t n) {
    sv_core_list_t cores;
    int32_t *orph; int32_t norph;
    sv_split_cores_orphans(labels_s1, n, &cores, &orph, &norph);
    int32_t n_cores = cores.count;
    int64_t ntotal_in_cores = 0;
    for (int32_t c = 0; c < n_cores; c++) ntotal_in_cores += cores.sizes[c];
    int64_t total = ntotal_in_cores + norph;
    cal_coverage_t r;
    r.core_count = n_cores;
    r.orphan_count = norph;
    r.frac_reads_in_cores = total > 0 ? (double)ntotal_in_cores / (double)total : 0.0;
    free(orph);
    sv_core_list_free(&cores);
    return r;
}
