/*
** starcodev.c -- see starcodev.h and CONVENTIONS_STARCODEV.md.
**
** Union-find, canonical partition, S1, S2 and the filtered searches of S3
** and S4. Wherever two choices tie (medoid, vote, absorption target) the
** rule is fixed and documented next to the code (K1-K5 in
** CONVENTIONS_STARCODEV.md), so the result does not depend on input order
** or on the order in which pairs are visited.
*/
#include "starcodev.h"
#include "pairgen.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define BASES "ACGT"
static inline int base_idx(char c) {
    switch (c) {
        case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3;
        default:  return -1;
    }
}

/* ================================================================== u64set
** A small open-addressing hash set for deduplicating pairs (i<<32)|j,
** used by S1, S3 and S4 to skip pairs that were already checked. Deterministic:
** it's only ever used to SKIP redundant work, never to influence the
** actual decision order (unions/absorptions still go through the combined
** K1/K5 keys). */
typedef struct {
    uint64_t *keys;
    uint8_t  *used;
    size_t    cap;
    size_t    count;
} u64set_t;

static void u64set_init(u64set_t *s, size_t cap_hint) {
    size_t cap = 16;
    while (cap < cap_hint * 2) cap <<= 1;
    s->keys = (uint64_t *)calloc(cap, sizeof(uint64_t));
    s->used = (uint8_t *)calloc(cap, 1);
    s->cap = cap;
    s->count = 0;
}

static void u64set_free(u64set_t *s) {
    free(s->keys); free(s->used);
    s->keys = NULL; s->used = NULL; s->cap = 0; s->count = 0;
}

static inline size_t u64hash(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return (size_t)x;
}

static void u64set_grow(u64set_t *s) {
    size_t newcap = s->cap * 2;
    uint64_t *nk = (uint64_t *)calloc(newcap, sizeof(uint64_t));
    uint8_t  *nu = (uint8_t *)calloc(newcap, 1);
    for (size_t i = 0; i < s->cap; i++) {
        if (!s->used[i]) continue;
        size_t h = u64hash(s->keys[i]) & (newcap - 1);
        while (nu[h]) h = (h + 1) & (newcap - 1);
        nu[h] = 1; nk[h] = s->keys[i];
    }
    free(s->keys); free(s->used);
    s->keys = nk; s->used = nu; s->cap = newcap;
}

/* Returns 1 if the key was NEWLY inserted (not seen before), 0 if it was
** already present. */
static int u64set_add(u64set_t *s, uint64_t key) {
    if (s->count * 2 >= s->cap) u64set_grow(s);
    size_t h = u64hash(key) & (s->cap - 1);
    while (s->used[h]) {
        if (s->keys[h] == key) return 0;
        h = (h + 1) & (s->cap - 1);
    }
    s->used[h] = 1; s->keys[h] = key; s->count++;
    return 1;
}

/* ================================================================== union-find */
void sv_uf_init(sv_uf_t *uf, int32_t n) {
    uf->p = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    for (int32_t i = 0; i < n; i++) uf->p[i] = i;
    uf->n = n;
}

int32_t sv_uf_find(sv_uf_t *uf, int32_t x) {
    int32_t *p = uf->p;
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

int sv_uf_union(sv_uf_t *uf, int32_t a, int32_t b) {
    int32_t ra = sv_uf_find(uf, a), rb = sv_uf_find(uf, b);
    if (ra == rb) return 0;
    /* K5: the root is always the smaller index -> the resulting tree does
    ** not depend on call order. */
    if (rb < ra) { int32_t t = ra; ra = rb; rb = t; }
    uf->p[rb] = ra;
    return 1;
}

void sv_uf_labels(sv_uf_t *uf, int32_t *out) {
    for (int32_t i = 0; i < uf->n; i++) out[i] = sv_uf_find(uf, i);
}

void sv_uf_free(sv_uf_t *uf) {
    free(uf->p); uf->p = NULL; uf->n = 0;
}

/* ================================================================== K5: canonical partition */
void sv_canonical_partition(const int32_t *labels, const char *const *seqs,
                             int32_t n, int32_t *out) {
    /* rep[l] = representative string (lexicographically smallest member)
    ** for raw label l. Mapped through a plain array of size n -- valid
    ** because raw union-find labels are always < n (a root is always some
    ** index). */
    const char **rep = (const char **)calloc((size_t)n, sizeof(char *));
    int32_t *rep_len = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    uint8_t *has = (uint8_t *)calloc((size_t)n, 1);
    for (int32_t i = 0; i < n; i++) {
        int32_t l = labels[i];
        const char *s = seqs[i];
        if (!has[l] || strcmp(s, rep[l]) < 0) {
            rep[l] = s; has[l] = 1;
        }
    }
    /* Collect the unique labels that have a representative, sorted by
    ** their representative string. */
    int32_t *uniq = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t nu = 0;
    for (int32_t l = 0; l < n; l++) if (has[l]) uniq[nu++] = l;

    /* strcmp-based insertion sort -- nu is usually much smaller than n */
    for (int32_t i = 1; i < nu; i++) {
        int32_t key = uniq[i];
        int32_t j = i - 1;
        while (j >= 0 && strcmp(rep[uniq[j]], rep[key]) > 0) {
            uniq[j + 1] = uniq[j]; j--;
        }
        uniq[j + 1] = key;
    }
    int32_t *remap = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < nu; k++) remap[uniq[k]] = k;
    for (int32_t i = 0; i < n; i++) out[i] = remap[labels[i]];

    free((void *)rep); free(rep_len); free(has); free(uniq); free(remap);
}

int64_t sv_partitions_equal(const int32_t *la, const int32_t *lb,
                             const char *const *seqs, int32_t n) {
    int32_t *ca = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t *cb = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    sv_canonical_partition(la, seqs, n, ca);
    sv_canonical_partition(lb, seqs, n, cb);
    int64_t diff = 0;
    for (int32_t i = 0; i < n; i++) if (ca[i] != cb[i]) diff++;
    free(ca); free(cb);
    return diff;
}

/* ================================================================== Levenshtein
** Row-scrolling Wagner-Fischer DP (O(min(la,lb)) memory). The return value
** is ALWAYS the true edit distance (not strictly clamped to the cutoff) --
** score_cutoff is only used for a banded early exit once a whole DP row
** has already gone past the cutoff on every reachable cell. */
int32_t sv_levenshtein(const char *a, int32_t la, const char *b, int32_t lb,
                        int32_t score_cutoff) {
    if (la == 0) return lb;
    if (lb == 0) return la;
    /* make sure b is the shorter string, to keep the scrolling row minimal */
    if (la < lb) { const char *t = a; a = b; b = t; int32_t tl = la; la = lb; lb = tl; }
    int32_t *prev = (int32_t *)malloc((size_t)(lb + 1) * sizeof(int32_t));
    int32_t *cur  = (int32_t *)malloc((size_t)(lb + 1) * sizeof(int32_t));
    for (int32_t j = 0; j <= lb; j++) prev[j] = j;
    for (int32_t i = 1; i <= la; i++) {
        cur[0] = i;
        int32_t rowmin = cur[0];
        char ca = a[i - 1];
        for (int32_t j = 1; j <= lb; j++) {
            int32_t cost_sub = prev[j - 1] + (ca != b[j - 1] ? 1 : 0);
            int32_t cost_del = prev[j] + 1;
            int32_t cost_ins = cur[j - 1] + 1;
            int32_t v = cost_sub < cost_del ? cost_sub : cost_del;
            if (cost_ins < v) v = cost_ins;
            cur[j] = v;
            if (v < rowmin) rowmin = v;
        }
        if (score_cutoff >= 0 && rowmin > score_cutoff) {
            /* the whole row is already past the cutoff: the final result is
            ** guaranteed to be > cutoff. Return a sentinel > cutoff (the
            ** contract is: the caller only ever tests "<= cutoff" or not). */
            free(prev); free(cur);
            return score_cutoff + 1;
        }
        int32_t *t = prev; prev = cur; cur = t;
    }
    int32_t result = prev[lb];
    free(prev); free(cur);
    return result;
}

/* ================================================================== STAGE 1 */
void sv_s1_cores(const char *const *seqs, const int32_t *lens, int32_t n,
                  int32_t tau_core, int32_t n_iter, int32_t theta, int32_t anchor_occ,
                  int32_t *lab_out, sv_work_t *work_out, int32_t *iters_done_out) {
    pg_encoded_t enc;
    pg_encode(seqs, lens, n, &enc);
    pg_signatures_t sig;
    pg_signatures(&enc, &sig);

    sv_uf_t uf;
    sv_uf_init(&uf, n);

    u64set_t seen;
    u64set_init(&seen, 4096);

    sv_work_t work = {0, 0, 0, 0};
    if (anchor_occ < 1) anchor_occ = 1;
    uint64_t *keys = (uint64_t *)malloc((size_t)n * anchor_occ * sizeof(uint64_t));
    int32_t *owner = anchor_occ > 1 ? (int32_t *)malloc((size_t)n * anchor_occ * sizeof(int32_t)) : NULL;
    int32_t it;
    for (it = 0; it < n_iter; it++) {
        pg_pairlist_t pl;
        if (anchor_occ == 1) {          /* one key per read: plain bucket sort is enough */
            pg_anchor_bucket(&enc, it, keys);
            pg_candidate_pairs(keys, n, PG_BUCKET_CAP, &pl);
        } else {
            int64_t nk = pg_anchor_keys_multi(&enc, it, anchor_occ, keys, owner);
            pg_candidate_pairs_owned(keys, owner, nk, PG_BUCKET_CAP, &pl);
        }
        work.cand += pl.count;
        int64_t merged = 0;
        for (int64_t k = 0; k < pl.count; k++) {
            int64_t i = pl.a[k], j = pl.b[k];
            int dh = pg_sig_hamming(&sig, i, j);
            if (dh > theta) continue;
            work.sig_pass++;
            uint64_t key = ((uint64_t)i << 32) | (uint64_t)(uint32_t)j;
            if (!u64set_add(&seen, key)) continue;
            work.edit++;
            int32_t d = sv_levenshtein(seqs[i], lens[i], seqs[j], lens[j], tau_core);
            if (d <= tau_core) {
                if (sv_uf_union(&uf, (int32_t)i, (int32_t)j)) merged++;
            }
        }
        work.uni += merged;
        pg_pairlist_free(&pl);
    }
    free(keys);
    free(owner);
    sv_uf_labels(&uf, lab_out);
    if (work_out) *work_out = work;
    if (iters_done_out) *iters_done_out = it; /* == n_iter (there is no early stop) */

    sv_uf_free(&uf);
    pg_free_signatures(&sig);
    pg_free_encoded(&enc);
}

void sv_split_cores_orphans(const int32_t *lab, int32_t n,
                             sv_core_list_t *cores_out,
                             int32_t **orphans_out, int32_t *n_orph_out) {
    /* Count the size of each raw label (label < n, since union-find roots
    ** are always indices). */
    int32_t *cnt = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    for (int32_t i = 0; i < n; i++) cnt[lab[i]]++;

    /* Collect labels with >1 member, sorted ascending by label value (this
    ** gives a fixed component order, since raw union-find labels are
    ** already directly sortable integers). */
    int32_t *core_labels = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t ncore = 0;
    for (int32_t l = 0; l < n; l++) if (cnt[l] > 1) core_labels[ncore++] = l;

    int32_t **members = (int32_t **)malloc((size_t)ncore * sizeof(int32_t *));
    int32_t *sizes = (int32_t *)malloc((size_t)ncore * sizeof(int32_t));
    int32_t *fillpos = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    int32_t *label_to_coreidx = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    for (int32_t l = 0; l < n; l++) label_to_coreidx[l] = -1;
    for (int32_t c = 0; c < ncore; c++) {
        members[c] = (int32_t *)malloc((size_t)cnt[core_labels[c]] * sizeof(int32_t));
        sizes[c] = cnt[core_labels[c]];
        label_to_coreidx[core_labels[c]] = c;
    }
    for (int32_t i = 0; i < n; i++) {
        int32_t c = label_to_coreidx[lab[i]];
        if (c >= 0) {
            members[c][fillpos[lab[i]]++] = i;  /* i is ascending because the loop is ascending -> already sorted */
        }
    }
    free(fillpos); free(label_to_coreidx); free(core_labels);

    cores_out->members = members;
    cores_out->sizes = sizes;
    cores_out->count = ncore;

    int32_t norph = n;
    for (int32_t c = 0; c < ncore; c++) norph -= sizes[c];
    int32_t *orph = (int32_t *)malloc((size_t)norph * sizeof(int32_t));
    int32_t op = 0;
    for (int32_t i = 0; i < n; i++) if (cnt[lab[i]] <= 1) orph[op++] = i;
    free(cnt);
    *orphans_out = orph;
    *n_orph_out = norph;
}

void sv_core_list_free(sv_core_list_t *cl) {
    for (int32_t c = 0; c < cl->count; c++) free(cl->members[c]);
    free(cl->members); free(cl->sizes);
    cl->members = NULL; cl->sizes = NULL; cl->count = 0;
}

/* ================================================================== STAGE 2 */

/* K2 -- medoid selection with a canonical tie-breaker. Among the members
** tied for minimum summed distance, pick the lexicographically smallest;
** if the strings themselves are identical, pick the smallest global
** index. K4 subsampling: when |members| > M_MAX, take the M_MAX
** lexicographically smallest members (not the first M_MAX in input order). */
static int cmp_seq_then_idx(const void *pa, const void *pb, const char *const *seqs) {
    int32_t ia = *(const int32_t *)pa, ib = *(const int32_t *)pb;
    int c = strcmp(seqs[ia], seqs[ib]);
    if (c != 0) return c;
    return ia < ib ? -1 : (ia > ib ? 1 : 0);
}
/* qsort_r isn't portable (BSD and GNU disagree on argument order) -- we
** use a small temporary global for the comparison context instead. This is
** safe because everything here is single-threaded and calls are never
** re-entered while one is still in flight. */
static const char *const *g_cmp_seqs;
static int cmp_seq_then_idx_g(const void *pa, const void *pb) {
    return cmp_seq_then_idx(pa, pb, g_cmp_seqs);
}

static int32_t medoid_pick(const int32_t *members, int32_t nmem,
                            const char *const *seqs, int32_t *n_tie_out) {
    int32_t nsub = nmem > SV_M_MAX ? SV_M_MAX : nmem;
    int32_t *sub = (int32_t *)malloc((size_t)nmem * sizeof(int32_t));
    memcpy(sub, members, (size_t)nmem * sizeof(int32_t));
    g_cmp_seqs = seqs;
    qsort(sub, (size_t)nmem, sizeof(int32_t), cmp_seq_then_idx_g);
    /* sub is now sorted ascending by (seq, idx); take the first nsub = K4 */

    const char **ss = (const char **)malloc((size_t)nsub * sizeof(char *));
    int32_t *slen = (int32_t *)malloc((size_t)nsub * sizeof(int32_t));
    for (int32_t k = 0; k < nsub; k++) { ss[k] = seqs[sub[k]]; slen[k] = (int32_t)strlen(ss[k]); }

    int64_t *tot = (int64_t *)calloc((size_t)nsub, sizeof(int64_t));
    for (int32_t x = 0; x < nsub; x++) {
        for (int32_t y = x + 1; y < nsub; y++) {
            int32_t d = sv_levenshtein(ss[x], slen[x], ss[y], slen[y], -1);
            tot[x] += d; tot[y] += d;
        }
    }
    int64_t best = tot[0];
    for (int32_t k = 1; k < nsub; k++) if (tot[k] < best) best = tot[k];
    int32_t n_tie = 0, best_k = -1;
    for (int32_t k = 0; k < nsub; k++) {
        if (tot[k] == best) {
            n_tie++;
            /* the subsample is already sorted ascending by (seq, idx), so
            ** the first k that reaches best is automatically the minimal
            ** (ss[k], sub[k]) among the ties. */
            if (best_k < 0) best_k = k;
        }
    }
    int32_t result = sub[best_k];
    free(ss); free(slen); free(tot); free(sub);
    if (n_tie_out) *n_tie_out = n_tie;
    return result;
}

/* Full Levenshtein backtrace (Wagner-Fischer, O(la*lb) memory) producing a
** position-annotated edit script, used for the column consensus. When
** several paths cost the same, the backtrace prefers delete, then diagonal
** (match/substitute), then insert; see the note inside the function. */
typedef enum { OP_EQ, OP_SUB, OP_DEL, OP_INS } op_tag_t;
typedef struct { op_tag_t tag; int32_t pos; char ch; } elem_op_t;

static int32_t align_ops(const char *a, int32_t la, const char *b, int32_t lb,
                          elem_op_t *ops_out /* capacity >= la+lb */) {
    int32_t *D = (int32_t *)malloc((size_t)(la + 1) * (size_t)(lb + 1) * sizeof(int32_t));
    #define DIDX(i, j) ((size_t)(i) * (size_t)(lb + 1) + (size_t)(j))
    for (int32_t j = 0; j <= lb; j++) D[DIDX(0, j)] = j;
    for (int32_t i = 0; i <= la; i++) D[DIDX(i, 0)] = i;
    for (int32_t i = 1; i <= la; i++) {
        for (int32_t j = 1; j <= lb; j++) {
            int32_t sub = D[DIDX(i - 1, j - 1)] + (a[i - 1] != b[j - 1] ? 1 : 0);
            int32_t del = D[DIDX(i - 1, j)] + 1;
            int32_t ins = D[DIDX(i, j - 1)] + 1;
            int32_t v = sub < del ? sub : del;
            if (ins < v) v = ins;
            D[DIDX(i, j)] = v;
        }
    }
    int32_t i = la, j = lb, k = 0;
    elem_op_t rev[/* VLA not used */1];
    (void)rev;
    /* Write directly into ops_out in reverse order, then flip it at the end.
    **
    ** TIE-BREAK PRIORITY ORDER: when several edit paths cost the same,
    ** delete is preferred over diagonal (match/replace), which is preferred
    ** over insert. The order is arbitrary but fixed, and it matters: in a
    ** homopolymer run ("XAAAAY" -> "XAAAY") every position of the run is an
    ** equally good place for the deletion, and this order always picks the
    ** first one, so all members of a core vote on the same column. */
    while (i > 0 || j > 0) {
        if (i > 0 && D[DIDX(i, j)] == D[DIDX(i - 1, j)] + 1) {
            ops_out[k].tag = OP_DEL;
            ops_out[k].pos = i - 1;
            ops_out[k].ch = 0;
            k++; i--;
        } else if (i > 0 && j > 0 && D[DIDX(i, j)] == D[DIDX(i - 1, j - 1)] + (a[i - 1] != b[j - 1] ? 1 : 0)) {
            ops_out[k].tag = (a[i - 1] == b[j - 1]) ? OP_EQ : OP_SUB;
            ops_out[k].pos = i - 1;
            ops_out[k].ch = b[j - 1];
            k++; i--; j--;
        } else {
            ops_out[k].tag = OP_INS;
            ops_out[k].pos = i;
            ops_out[k].ch = b[j - 1];
            k++; j--;
        }
    }
    free(D);
    /* flip back into forward order (anchor position ascending) */
    for (int32_t x = 0, y = k - 1; x < y; x++, y--) {
        elem_op_t t = ops_out[x]; ops_out[x] = ops_out[y]; ops_out[y] = t;
    }
    return k;
}

/* K3 -- indel-aware column profile with a canonical tie-breaker.
** Match/substitute columns: majority vote, ties broken by the fixed order
** A<C<G<T. Insertion columns: accepted only with a strict majority; a tie
** means DROP the column. */

/* Simple hash entry for ins_votes: key (pos, text) -> count. The number of
** distinct entries is bounded by the number of core members, so a linear
** table is fine (a StarcodeV core is typically tens to low hundreds of
** members). */
typedef struct { int32_t pos; char *text; int32_t count; } ins_entry_t;

static ins_entry_t *ins_find_or_add(ins_entry_t *tab, int32_t *ntab, int32_t cap,
                                     int32_t pos, const char *text, int32_t tlen) {
    for (int32_t t = 0; t < *ntab; t++) {
        if (tab[t].pos == pos && (int32_t)strlen(tab[t].text) == tlen &&
            memcmp(tab[t].text, text, (size_t)tlen) == 0) {
            return &tab[t];
        }
    }
    if (*ntab >= cap) return NULL; /* should never happen, cap is chosen conservatively */
    ins_entry_t *e = &tab[*ntab];
    e->pos = pos;
    e->text = (char *)malloc((size_t)tlen + 1);
    memcpy(e->text, text, (size_t)tlen); e->text[tlen] = 0;
    e->count = 0;
    (*ntab)++;
    return e;
}

static void consensus_from_anchor(const char *anchor, int32_t la,
                                   const int32_t *members, int32_t nmem,
                                   const char *const *seqs,
                                   char **out_str, int32_t *out_tieflag) {
    int32_t (*sub_votes)[4] = (int32_t (*)[4])calloc((size_t)la, sizeof(int32_t[4]));
    int32_t *del_votes = (int32_t *)calloc((size_t)la, sizeof(int32_t));
    int32_t m = nmem;

    /* ins_votes table: at most (la+maxlen) entries per member, sized generously */
    int32_t ins_cap = 0;
    for (int32_t k = 0; k < nmem; k++) ins_cap += (int32_t)strlen(seqs[members[k]]) + 2;
    ins_entry_t *ins_tab = (ins_entry_t *)malloc((size_t)(ins_cap > 0 ? ins_cap : 1) * sizeof(ins_entry_t));
    int32_t n_ins_tab = 0;

    for (int32_t k = 0; k < nmem; k++) {
        const char *s = seqs[members[k]];
        int32_t ls = (int32_t)strlen(s);
        int32_t cap_ops = la + ls + 2;
        elem_op_t *ops = (elem_op_t *)malloc((size_t)cap_ops * sizeof(elem_op_t));
        int32_t nops = align_ops(anchor, la, s, ls, ops);

        int32_t run_pos = -1, run_len = 0, run_cap = 64;
        char *run_buf = (char *)malloc((size_t)run_cap);
        for (int32_t oi = 0; oi < nops; oi++) {
            elem_op_t *op = &ops[oi];
            if (op->tag == OP_INS) {
                if (run_pos != op->pos && run_len > 0) {
                    ins_entry_t *e = ins_find_or_add(ins_tab, &n_ins_tab, ins_cap, run_pos, run_buf, run_len);
                    if (e) e->count++;
                    run_len = 0;
                }
                run_pos = op->pos;
                if (run_len + 1 > run_cap) { run_cap *= 2; run_buf = (char *)realloc(run_buf, (size_t)run_cap); }
                run_buf[run_len++] = op->ch;
            } else {
                if (run_len > 0) {
                    ins_entry_t *e = ins_find_or_add(ins_tab, &n_ins_tab, ins_cap, run_pos, run_buf, run_len);
                    if (e) e->count++;
                    run_len = 0;
                }
                if (op->tag == OP_EQ || op->tag == OP_SUB) {
                    int bi = base_idx(op->ch);
                    if (bi >= 0) sub_votes[op->pos][bi]++;
                } else { /* OP_DEL */
                    del_votes[op->pos]++;
                }
            }
        }
        if (run_len > 0) {
            ins_entry_t *e = ins_find_or_add(ins_tab, &n_ins_tab, ins_cap, run_pos, run_buf, run_len);
            if (e) e->count++;
        }
        free(run_buf);
        free(ops);
    }

    /* Pick the best insertion per position: strict majority (2*c > m); ties
    ** between different texts at the same position are broken canonically
    ** (largest count, then shortest text, then lexicographically smallest)
    ** -- i.e. the maximum of the tuple (count, -length, text). */
    int32_t *best_ins_for_pos = (int32_t *)malloc((size_t)(la + 1) * sizeof(int32_t));
    for (int32_t p = 0; p <= la; p++) best_ins_for_pos[p] = -1;
    for (int32_t t = 0; t < n_ins_tab; t++) {
        ins_entry_t *e = &ins_tab[t];
        if (2 * e->count <= m) continue;
        int32_t p = e->pos;
        int32_t cur = best_ins_for_pos[p];
        if (cur < 0) { best_ins_for_pos[p] = t; continue; }
        ins_entry_t *ce = &ins_tab[cur];
        int32_t clen = (int32_t)strlen(ce->text), tlen = (int32_t)strlen(e->text);
        /* compare (count, -len, text) -- the maximum wins */
        int take_new = 0;
        if (e->count != ce->count) take_new = e->count > ce->count;
        else if (tlen != clen) take_new = tlen < clen;
        else take_new = strcmp(e->text, ce->text) < 0;
        if (take_new) best_ins_for_pos[p] = t;
    }

    int32_t tie_col = 0;
    size_t outcap = (size_t)la * 2 + 16;
    char *out = (char *)malloc(outcap);
    size_t olen = 0;
    #define OUT_PUSH_STR(str, slen_) do { \
        if (olen + (size_t)(slen_) + 1 > outcap) { outcap = (olen + (size_t)(slen_) + 16) * 2; out = (char*)realloc(out, outcap); } \
        memcpy(out + olen, (str), (size_t)(slen_)); olen += (size_t)(slen_); } while (0)
    #define OUT_PUSH_CH(c) do { \
        if (olen + 2 > outcap) { outcap *= 2; out = (char*)realloc(out, outcap); } \
        out[olen++] = (c); } while (0)

    for (int32_t p = 0; p < la; p++) {
        if (best_ins_for_pos[p] >= 0) {
            ins_entry_t *e = &ins_tab[best_ins_for_pos[p]];
            OUT_PUSH_STR(e->text, (int32_t)strlen(e->text));
        }
        if (2 * del_votes[p] > m) continue; /* a majority favors deletion -> skip this column */
        int32_t *v = sub_votes[p];
        int32_t vsum = v[0] + v[1] + v[2] + v[3];
        if (vsum == 0) { OUT_PUSH_CH(anchor[p]); continue; }
        int32_t mx = v[0]; int32_t argmax = 0; int32_t ntie = 0;
        for (int32_t b = 0; b < 4; b++) if (v[b] == mx) ntie++;
        for (int32_t b = 1; b < 4; b++) if (v[b] > mx) { mx = v[b]; argmax = b; }
        ntie = 0;
        for (int32_t b = 0; b < 4; b++) if (v[b] == mx) ntie++;
        if (ntie > 1) tie_col++;
        OUT_PUSH_CH(BASES[argmax]);
    }
    if (best_ins_for_pos[la] >= 0) {
        ins_entry_t *e = &ins_tab[best_ins_for_pos[la]];
        OUT_PUSH_STR(e->text, (int32_t)strlen(e->text));
    }
    OUT_PUSH_CH('\0'); olen--; /* null-terminate without counting it towards the length */

    for (int32_t t = 0; t < n_ins_tab; t++) free(ins_tab[t].text);
    free(ins_tab); free(best_ins_for_pos);
    free(sub_votes); free(del_votes);

    *out_str = out;
    *out_tieflag = tie_col;
}

void sv_s2_consensus(const sv_core_list_t *cores, const char *const *seqs,
                      sv_s2_result_t *out) {
    int32_t nc = cores->count;
    char **cons = (char **)malloc((size_t)nc * sizeof(char *));
    int32_t *sizes = (int32_t *)malloc((size_t)nc * sizeof(int32_t));
    int32_t *anchors = (int32_t *)malloc((size_t)nc * sizeof(int32_t));
    int32_t tie_med = 0, tie_col = 0;

    for (int32_t c = 0; c < nc; c++) {
        int32_t nmem = cores->sizes[c];
        const int32_t *mem = cores->members[c];
        int32_t ntie_medoid;
        int32_t anc_i = medoid_pick(mem, nmem, seqs, &ntie_medoid);
        if (ntie_medoid > 1) tie_med++;
        int32_t la = (int32_t)strlen(seqs[anc_i]);
        char *consensus_str; int32_t ct;
        consensus_from_anchor(seqs[anc_i], la, mem, nmem, seqs, &consensus_str, &ct);
        if (ct > 0) tie_col++;
        cons[c] = consensus_str;
        sizes[c] = nmem;
        anchors[c] = anc_i;
    }
    out->cons = cons; out->sizes = sizes; out->anchors = anchors;
    out->count = nc; out->tie_medoid = tie_med; out->tie_column = tie_col;
}

void sv_s2_result_free(sv_s2_result_t *r) {
    for (int32_t c = 0; c < r->count; c++) free(r->cons[c]);
    free(r->cons); free(r->sizes); free(r->anchors);
    r->cons = NULL; r->sizes = NULL; r->anchors = NULL; r->count = 0;
}

/* ================================================================== STAGE 4: candidate search
** K1 -- combined 64-bit decision key (d << 32) | consensus_idx. `min` over
** this key is associative-commutative: reducing it in any order gives the
** same result (on a GPU this is a single atomicMin on an unsigned long long). */
void sv_s4_absorb(const char *const *orph_seqs, const int32_t *orph_lens, int32_t norph,
                   const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                   int32_t tau_abs, int32_t sig_abs, int32_t n_iter, int32_t anchor_occ,
                   int32_t *bd_out, int32_t *bi_out, uint8_t *absorbed_out,
                   sv_work_t *work_out) {
    int32_t ntot = norph + ncons;
    if (anchor_occ < 1) anchor_occ = 1;
    const char **allseq = (const char **)malloc((size_t)ntot * sizeof(char *));
    int32_t *alllens = (int32_t *)malloc((size_t)ntot * sizeof(int32_t));
    uint8_t *is_cons = (uint8_t *)calloc((size_t)ntot, 1);
    for (int32_t i = 0; i < norph; i++) { allseq[i] = orph_seqs[i]; alllens[i] = orph_lens[i]; }
    for (int32_t i = 0; i < ncons; i++) {
        allseq[norph + i] = cons[i]; alllens[norph + i] = cons_lens[i];
        is_cons[norph + i] = 1;
    }

    pg_encoded_t enc;
    pg_encode(allseq, alllens, ntot, &enc);
    pg_signatures_t sig;
    pg_signatures(&enc, &sig);

    /* best key per orphan; sentinel = (tau_abs+1) << 32 | 0xFFFFFFFF */
    int64_t *best = (int64_t *)malloc((size_t)norph * sizeof(int64_t));
    for (int32_t i = 0; i < norph; i++)
        best[i] = ((int64_t)(tau_abs + 1) << 32) | (int64_t)0xFFFFFFFFLL;

    sv_work_t work = {0, 0, 0, 0};
    u64set_t seen;
    u64set_init(&seen, 4096);
    uint64_t *keys = (uint64_t *)malloc((size_t)ntot * anchor_occ * sizeof(uint64_t));
    int32_t *owner = anchor_occ > 1 ? (int32_t *)malloc((size_t)ntot * anchor_occ * sizeof(int32_t)) : NULL;

    for (int32_t it = 0; it < n_iter; it++) {
        pg_pairlist_t pl;
        if (anchor_occ == 1) {
            pg_anchor_bucket(&enc, it, keys);
            pg_candidate_pairs(keys, ntot, PG_BUCKET_CAP, &pl);
        } else {
            int64_t nk = pg_anchor_keys_multi(&enc, it, anchor_occ, keys, owner);
            pg_candidate_pairs_owned(keys, owner, nk, PG_BUCKET_CAP, &pl);
        }
        for (int64_t k = 0; k < pl.count; k++) {
            int64_t x = pl.a[k], y = pl.b[k];
            uint8_t cx = is_cons[x], cy = is_cons[y];
            if (cx == cy) continue;  /* only cross-type pairs matter here */
            int64_t oi = cx ? y : x;             /* orphan side (global index) */
            int64_t ci_global = cx ? x : y;       /* consensus side (global index) */
            int64_t ci = ci_global - norph;       /* local consensus index */
            work.cand++;
            int dh = pg_sig_hamming(&sig, oi, ci_global);
            if (dh > sig_abs) continue;
            work.sig_pass++;
            uint64_t key = ((uint64_t)oi << 32) | (uint64_t)(uint32_t)ci;
            if (!u64set_add(&seen, key)) continue;
            work.edit++;
            int32_t d = sv_levenshtein(orph_seqs[oi], orph_lens[oi], cons[ci], cons_lens[ci], tau_abs);
            if (d <= tau_abs) {
                int64_t cand_key = ((int64_t)d << 32) | ci;
                if (cand_key < best[oi]) best[oi] = cand_key;
            }
        }
        pg_pairlist_free(&pl);
    }
    free(keys);
    free(owner);
    u64set_free(&seen);

    for (int32_t i = 0; i < norph; i++) {
        int32_t bd = (int32_t)(best[i] >> 32);
        int32_t bi = (int32_t)(best[i] & 0xFFFFFFFFLL);
        bd_out[i] = bd;
        uint8_t ab = bd <= tau_abs;
        absorbed_out[i] = ab;
        bi_out[i] = ab ? bi : -1;
    }
    if (work_out) *work_out = work;

    free(best);
    pg_free_signatures(&sig);
    pg_free_encoded(&enc);
    free(allseq); free(alllens); free(is_cons);
}

/* ================================================================== STAGE 3: filtered merge (--approx) */
void sv_s3_merge_cons(const char *const *cons, const int32_t *cons_lens, int32_t ncons,
                       int32_t tau_cons, int32_t sig_abs, int32_t n_iter,
                       int32_t *lab_out, sv_work_t *work_out) {
    pg_encoded_t enc;
    pg_encode(cons, cons_lens, ncons, &enc);
    pg_signatures_t sig;
    pg_signatures(&enc, &sig);

    sv_uf_t uf;
    sv_uf_init(&uf, ncons);

    u64set_t seen;
    u64set_init(&seen, 4096);
    sv_work_t work = {0, 0, 0, 0};
    uint64_t *keys = (uint64_t *)malloc((size_t)ncons * sizeof(uint64_t));

    for (int32_t it = 0; it < n_iter; it++) {
        pg_anchor_bucket(&enc, it, keys);
        pg_pairlist_t pl;
        pg_candidate_pairs(keys, ncons, PG_BUCKET_CAP, &pl);
        work.cand += pl.count;
        for (int64_t k = 0; k < pl.count; k++) {
            int64_t i = pl.a[k], j = pl.b[k];
            int dh = pg_sig_hamming(&sig, i, j);
            if (dh > sig_abs) continue;
            work.sig_pass++;
            uint64_t key = ((uint64_t)i << 32) | (uint64_t)(uint32_t)j;
            if (!u64set_add(&seen, key)) continue;
            work.edit++;
            int32_t d = sv_levenshtein(cons[i], cons_lens[i], cons[j], cons_lens[j], tau_cons);
            if (d <= tau_cons) {
                if (sv_uf_union(&uf, (int32_t)i, (int32_t)j)) work.uni++;
            }
        }
        pg_pairlist_free(&pl);
    }
    free(keys);
    u64set_free(&seen);
    sv_uf_labels(&uf, lab_out);
    if (work_out) *work_out = work;
    sv_uf_free(&uf);
    pg_free_signatures(&sig);
    pg_free_encoded(&enc);
}

/* ================================================================== final assembly */
void sv_assemble(const sv_core_list_t *cores,
                  const int32_t *orph_idx, int32_t norph,
                  const int32_t *bi4, const uint8_t *absorbed,
                  const int32_t *lab_s3, int32_t n_reads,
                  int32_t *out) {
    for (int32_t i = 0; i < n_reads; i++) out[i] = -1;
    int32_t nc = cores->count;
    int32_t maxlab = -1;
    for (int32_t c = 0; c < nc; c++) {
        int32_t lab = lab_s3[c];
        if (lab > maxlab) maxlab = lab;
        const int32_t *mem = cores->members[c];
        for (int32_t k = 0; k < cores->sizes[c]; k++) out[mem[k]] = lab;
    }
    int32_t nxt = maxlab + 1;
    for (int32_t k = 0; k < norph; k++) {
        int32_t g = orph_idx[k];
        if (absorbed[k]) {
            out[g] = lab_s3[bi4[k]];
        } else {
            out[g] = nxt++;
        }
    }
}
