/*
** pairgen.c -- see pairgen.h and CONVENTIONS_STARCODEV.md before making
** changes here.
**
** This file produces the candidate pairs that replace Starcode's trie search.
** It must stay deterministic: no RNG, no dependency on input order, and a
** total order on every sort, so that the same input always gives the same
** pair list.
*/
#include "pairgen.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define PG_BASES "ACGT"

/* ---------------------------------------------------------------- anchor table */
static int32_t g_anchors[1 << (2 * PG_W_ANCHOR)];   /* 4^PG_W_ANCHOR */
static int32_t g_anchors_len = 0;
static int     g_anchors_ready = 0;

void pg_init_anchors(void) {
    if (g_anchors_ready) return;
    const int32_t nmax = 1 << (2 * PG_W_ANCHOR);   /* 4^W_ANCHOR */
    uint8_t *seen = (uint8_t *)calloc((size_t)nmax, 1);
    uint32_t x = 1;
    int32_t count = 0;
    while (count < nmax) {
        x = (uint32_t)(1664525u * x + 1013904223u);   /* fixed LCG constants (Numerical Recipes), so the anchor order never changes */
        int32_t v = (int32_t)(x % (uint32_t)nmax);
        if (!seen[v]) {
            seen[v] = 1;
            g_anchors[count++] = v;
        }
    }
    free(seen);
    g_anchors_len = count;
    g_anchors_ready = 1;
}

const int32_t *pg_anchor_table(int32_t *len_out) {
    pg_init_anchors();
    if (len_out) *len_out = g_anchors_len;
    return g_anchors;
}

/* ---------------------------------------------------------------- popcount */
static uint8_t g_popc[256];
static int     g_popc_ready = 0;

static void ensure_popc(void) {
    if (g_popc_ready) return;
    for (int i = 0; i < 256; i++) {
        int c = 0, v = i;
        while (v) { c += v & 1; v >>= 1; }
        g_popc[i] = (uint8_t)c;
    }
    g_popc_ready = 1;
}

/* ---------------------------------------------------------------- encode */
static inline int base_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 255;
    }
}

void pg_encode(const char *const *seqs, const int32_t *lens, int32_t n,
               pg_encoded_t *out) {
    pg_init_anchors();
    int32_t maxlen = 0;
    for (int32_t i = 0; i < n; i++) if (lens[i] > maxlen) maxlen = lens[i];
    int32_t glen = maxlen - PG_W_ANCHOR + 1;
    if (glen < 0) glen = 0;

    uint8_t *T = (uint8_t *)malloc((size_t)n * (size_t)maxlen);
    int32_t *L = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t *G = (int32_t *)malloc((size_t)n * (size_t)glen * sizeof(int32_t));
    memset(T, 255, (size_t)n * (size_t)maxlen);

    for (int32_t i = 0; i < n; i++) {
        L[i] = lens[i];
        const char *s = seqs[i];
        uint8_t *row = T + (size_t)i * maxlen;
        for (int32_t k = 0; k < lens[i]; k++) row[k] = (uint8_t)base_index(s[k]);
    }

    for (int32_t i = 0; i < n; i++) {
        const uint8_t *row = T + (size_t)i * maxlen;
        int32_t *grow = G + (size_t)i * glen;
        for (int32_t k = 0; k < glen; k++) {
            int32_t code = 0;
            int ok = 1;
            for (int j = 0; j < PG_W_ANCHOR; j++) {
                uint8_t b = row[k + j];
                if (b == 255) { ok = 0; break; }
                code = code * 4 + b;
            }
            grow[k] = ok ? code : -1;
        }
    }

    out->n = n;
    out->maxlen = maxlen;
    out->glen = glen;
    out->T = T;
    out->L = L;
    out->G = G;
}

void pg_free_encoded(pg_encoded_t *enc) {
    if (!enc) return;
    free(enc->T); free(enc->L); free(enc->G);
    enc->T = NULL; enc->L = NULL; enc->G = NULL;
}

/* ---------------------------------------------------------------- signatures */
void pg_signatures(const pg_encoded_t *enc, pg_signatures_t *out) {
    int32_t n = enc->n, maxlen = enc->maxlen;
    uint8_t *S = (uint8_t *)calloc((size_t)n * PG_SIG_BYTES, 1);
    const int32_t nqg = 1 << (2 * PG_Q_SIG);             /* 4^Q_SIG */
    const int32_t max_blkid = (PG_SIG_BITS / nqg) - 1;   /* assumes nqg <= SIG_BITS */
    int32_t qlen = maxlen - PG_Q_SIG + 1;
    if (qlen < 0) qlen = 0;

    for (int32_t i = 0; i < n; i++) {
        const uint8_t *row = enc->T + (size_t)i * maxlen;
        uint8_t *srow = S + (size_t)i * PG_SIG_BYTES;
        for (int32_t k = 0; k < qlen; k++) {
            int32_t code = 0;
            int ok = 1;
            for (int j = 0; j < PG_Q_SIG; j++) {
                uint8_t b = row[k + j];
                if (b == 255) { ok = 0; break; }
                code = code * 4 + b;
            }
            if (!ok) continue;
            int32_t blkid = k / PG_SIG_BLOCK;
            if (blkid > max_blkid) blkid = max_blkid;
            int32_t bit = (blkid * nqg + code) % PG_SIG_BITS;
            srow[bit >> 3] |= (uint8_t)(1u << (bit & 7));
        }
    }
    out->n = n;
    out->S = S;
}

void pg_free_signatures(pg_signatures_t *sig) {
    if (!sig) return;
    free(sig->S); sig->S = NULL;
}

int pg_sig_hamming(const pg_signatures_t *sig, int64_t a, int64_t b) {
    ensure_popc();
    const uint8_t *sa = sig->S + (size_t)a * PG_SIG_BYTES;
    const uint8_t *sb = sig->S + (size_t)b * PG_SIG_BYTES;
    int d = 0;
    for (int k = 0; k < PG_SIG_BYTES; k++) d += g_popc[sa[k] ^ sb[k]];
    return d;
}

/* ---------------------------------------------------------------- anchor_bucket */
void pg_anchor_bucket(const pg_encoded_t *enc, int32_t it, uint64_t *keys_out) {
    int32_t n = enc->n, maxlen = enc->maxlen, glen = enc->glen;
    int32_t alen; const int32_t *anchors = pg_anchor_table(&alen);
    int32_t a = anchors[it % alen];

    for (int32_t i = 0; i < n; i++) {
        const int32_t *grow = enc->G + (size_t)i * glen;
        int32_t pos = 0;
        int has = 0;
        for (int32_t k = 0; k < glen; k++) {
            if (grow[k] == a) { pos = k; has = 1; break; }
        }
        if (!has) {
            keys_out[i] = ((uint64_t)1 << 63) | (uint64_t)(uint32_t)i;
            continue;
        }
        int32_t start = pos + PG_W_ANCHOR;
        uint64_t key = 0;
        const uint8_t *row = enc->T + (size_t)i * maxlen;
        for (int j = 0; j < PG_L_HASH; j++) {
            int32_t p = start + j;
            int64_t ch;
            if (p >= enc->L[i]) {
                ch = 4;
            } else {
                int32_t idx = p; if (idx < 0) idx = 0; if (idx >= maxlen) idx = maxlen - 1;
                ch = (int64_t)row[idx];
            }
            key = key * 5u + (uint64_t)ch;
        }
        keys_out[i] = key;
    }
}

/* ---------------------------------------------------------------- candidate_pairs */
typedef struct { uint64_t key; int32_t idx; } pg_keyed_t;

static int pg_keyed_cmp(const void *pa, const void *pb) {
    const pg_keyed_t *a = (const pg_keyed_t *)pa;
    const pg_keyed_t *b = (const pg_keyed_t *)pb;
    if (a->key < b->key) return -1;
    if (a->key > b->key) return 1;
    return 0;
}

void pg_pairlist_init(pg_pairlist_t *pl) {
    pl->a = NULL; pl->b = NULL; pl->count = 0; pl->cap = 0;
}

void pg_pairlist_push(pg_pairlist_t *pl, int64_t a, int64_t b) {
    if (pl->count >= pl->cap) {
        int64_t newcap = pl->cap ? pl->cap * 2 : 1024;
        pl->a = (int64_t *)realloc(pl->a, (size_t)newcap * sizeof(int64_t));
        pl->b = (int64_t *)realloc(pl->b, (size_t)newcap * sizeof(int64_t));
        pl->cap = newcap;
    }
    pl->a[pl->count] = a < b ? a : b;
    pl->b[pl->count] = a < b ? b : a;
    pl->count++;
}

void pg_pairlist_free(pg_pairlist_t *pl) {
    free(pl->a); free(pl->b);
    pl->a = NULL; pl->b = NULL; pl->count = 0; pl->cap = 0;
}

void pg_candidate_pairs(const uint64_t *keys, int32_t n, int32_t cap,
                         pg_pairlist_t *out) {
    pg_pairlist_init(out);
    if (n <= 0) return;
    pg_keyed_t *ks = (pg_keyed_t *)malloc((size_t)n * sizeof(pg_keyed_t));
    for (int32_t i = 0; i < n; i++) { ks[i].key = keys[i]; ks[i].idx = i; }
    qsort(ks, (size_t)n, sizeof(pg_keyed_t), pg_keyed_cmp);

    int32_t st = 0;
    while (st < n) {
        int32_t en = st + 1;
        while (en < n && ks[en].key == ks[st].key) en++;
        int32_t m = en - st;
        if (m >= 2 && m <= cap) {
            for (int32_t x = st; x < en; x++)
                for (int32_t y = x + 1; y < en; y++)
                    pg_pairlist_push(out, ks[x].idx, ks[y].idx);
        }
        st = en;
    }
    free(ks);
}

/* ---------------------------------------------------------------- multi-occurrence keys */
int64_t pg_anchor_keys_multi(const pg_encoded_t *enc, int32_t it, int32_t max_occ,
                             uint64_t *keys_out, int32_t *owner_out) {
    int32_t n = enc->n, maxlen = enc->maxlen, glen = enc->glen;
    int32_t alen; const int32_t *anchors = pg_anchor_table(&alen);
    int32_t a = anchors[it % alen];
    int64_t nk = 0;

    for (int32_t i = 0; i < n; i++) {
        const int32_t *grow = enc->G + (size_t)i * glen;
        const uint8_t *row = enc->T + (size_t)i * maxlen;
        int32_t found = 0;
        for (int32_t k = 0; k < glen && found < max_occ; k++) {
            if (grow[k] != a) continue;
            found++;
            /* same key as pg_anchor_bucket for this position */
            int32_t start = k + PG_W_ANCHOR;
            uint64_t key = 0;
            for (int j = 0; j < PG_L_HASH; j++) {
                int32_t p = start + j;
                int64_t ch;
                if (p >= enc->L[i]) {
                    ch = 4;
                } else {
                    int32_t idx = p; if (idx < 0) idx = 0; if (idx >= maxlen) idx = maxlen - 1;
                    ch = (int64_t)row[idx];
                }
                key = key * 5u + (uint64_t)ch;
            }
            keys_out[nk] = key;
            owner_out[nk] = i;
            nk++;
        }
    }
    return nk;
}

typedef struct { uint64_t key; int32_t owner; } pg_owned_t;

static int pg_owned_cmp(const void *pa, const void *pb) {
    const pg_owned_t *a = (const pg_owned_t *)pa, *b = (const pg_owned_t *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return (a->owner > b->owner) - (a->owner < b->owner);   /* full order, so the sort is deterministic */
}

void pg_candidate_pairs_owned(const uint64_t *keys, const int32_t *owner, int64_t count,
                              int32_t cap, pg_pairlist_t *out) {
    pg_pairlist_init(out);
    if (count <= 0) return;
    pg_owned_t *ks = (pg_owned_t *)malloc((size_t)count * sizeof(pg_owned_t));
    for (int64_t i = 0; i < count; i++) { ks[i].key = keys[i]; ks[i].owner = owner[i]; }
    qsort(ks, (size_t)count, sizeof(pg_owned_t), pg_owned_cmp);

    int64_t st = 0;
    while (st < count) {
        int64_t en = st + 1;
        while (en < count && ks[en].key == ks[st].key) en++;
        int64_t m = en - st;
        if (m >= 2 && m <= cap) {
            for (int64_t x = st; x < en; x++)
                for (int64_t y = x + 1; y < en; y++)
                    if (ks[x].owner != ks[y].owner)
                        pg_pairlist_push(out, ks[x].owner, ks[y].owner);
        }
        st = en;
    }
    free(ks);
}
