/*
** pairgen.h -- StarcodeV's candidate-pair generator (anchor hash + signature).
**
** Takes the role of src/trie.c and src/trie.h in Starcode
** (gui11aume/starcode v1.4); CONVENTIONS_STARCODEV.md §2 explains why. A trie
** answers "which strings are within radius tau of X?" with a node-by-node
** tree walk. This module answers it with anchor hashing plus a q-gram
** signature pre-filter, and every unit of work (hash one string, compare one
** pair) is independent of the others, so it can also run in parallel.
**
** Contents: encode, signatures, anchor keys, candidate pairs. No runtime RNG
** (the anchor table comes from fixed LCG constants, computed once), no
** floating point, no dependency on input order.
*/
#ifndef _PAIRGEN_HEADER
#define _PAIRGEN_HEADER

#include <stdint.h>
#include <stddef.h>

/* ---------------------------------------------------------------- constants */
#define PG_W_ANCHOR   4        /* anchor length */
#define PG_L_HASH     6        /* number of bases hashed after the anchor */
#define PG_Q_SIG      3        /* q-gram signature length */
#define PG_SIG_BLOCK  22       /* position-block width for signature blocking */
#define PG_SIG_BITS   384      /* total signature bits */
#define PG_SIG_BYTES  (PG_SIG_BITS / 8)
#define PG_BUCKET_CAP 4096     /* buckets larger than this are skipped (safety fuse) */

/* ---------------------------------------------------------------- encoding */
typedef struct {
    int32_t  n;        /* number of strings */
    int32_t  maxlen;   /* length of the longest string */
    int32_t  glen;     /* maxlen - PG_W_ANCHOR + 1 (width of the G column) */
    uint8_t *T;         /* n x maxlen, base codes 0-3 (A/C/G/T), or 255 past the string's end */
    int32_t *L;         /* n, length of each string */
    int32_t *G;         /* n x glen, rolling anchor w-gram code, -1 where invalid */
} pg_encoded_t;

typedef struct {
    int32_t  n;
    uint8_t *S;         /* n x PG_SIG_BYTES, binary q-gram signature */
} pg_signatures_t;

/* Dynamic candidate-pair list (a[i] < b[i] always, normalized). */
typedef struct {
    int64_t *a;
    int64_t *b;
    int64_t  count;
    int64_t  cap;
} pg_pairlist_t;

/* Anchor table: a permutation of 4^PG_W_ANCHOR values from a fixed-constant
** LCG (Numerical Recipes) -- a CONSTANT, not a runtime RNG. Built once by
** pg_init_anchors(); the values are identical on every machine and every
** run. */
void            pg_init_anchors(void);
const int32_t  *pg_anchor_table(int32_t *len_out);

/* seqs: array of char* (need not be null-terminated at a common length), lens: length of each seq. */
void pg_encode(const char *const *seqs, const int32_t *lens, int32_t n,
               pg_encoded_t *out);
void pg_free_encoded(pg_encoded_t *enc);

void pg_signatures(const pg_encoded_t *enc, pg_signatures_t *out);
void pg_free_signatures(pg_signatures_t *sig);

/* Hamming popcount between two signature rows. */
int pg_sig_hamming(const pg_signatures_t *sig, int64_t a, int64_t b);

/* Deterministic anchor hash for iteration `it`. keys_out must be allocated
** by the caller, sized enc->n. A string with no anchor gets a unique
** sentinel. */
void pg_anchor_bucket(const pg_encoded_t *enc, int32_t it, uint64_t *keys_out);

/* All pairs within a bucket, not just the adjacent ones after sorting: the
** all-pairs version finds many more true pairs per iteration. */
void pg_candidate_pairs(const uint64_t *keys, int32_t n, int32_t cap,
                         pg_pairlist_t *out);

/* One key per occurrence of the anchor, up to max_occ occurrences per
** string. An indel right after one occurrence changes that key, but a
** later occurrence past the indel can still match, so the pair is not lost.
** keys_out/owner_out need room for n * max_occ entries. Strings without the
** anchor emit nothing. Returns the number of keys written. With max_occ = 1
** the candidate pairs are the same as pg_anchor_bucket + pg_candidate_pairs. */
int64_t pg_anchor_keys_multi(const pg_encoded_t *enc, int32_t it, int32_t max_occ,
                             uint64_t *keys_out, int32_t *owner_out);

/* Pairs of different owners that share a key. The same pair can come out
** more than once (two shared occurrences); callers already de-duplicate. */
void pg_candidate_pairs_owned(const uint64_t *keys, const int32_t *owner, int64_t count,
                              int32_t cap, pg_pairlist_t *out);

void pg_pairlist_init(pg_pairlist_t *pl);
void pg_pairlist_push(pg_pairlist_t *pl, int64_t a, int64_t b);
void pg_pairlist_free(pg_pairlist_t *pl);

#endif /* _PAIRGEN_HEADER */
