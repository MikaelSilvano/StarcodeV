/*
** pairgen.h -- StarcodeV's candidate-pair generator (anchor-hash + signature).
**
** Replaces the role of src/trie.c and src/trie.h from the original
** Starcode (gui11aume/starcode v1.4). See CONVENTIONS_STARCODEV.md §2 for
** the rationale behind the replacement: a trie answers "which strings are
** within radius tau of X?" via a sequential, node-by-node tree walk; this
** module answers the same question with deterministic anchor hashing plus
** a q-gram signature pre-filter, whose unit of work (hash one string,
** compare one pair) is independent across strings/pairs -- a direct port
** target for a GPU kernel (p1_pairgen.cu).
**
** An exact mirror of starcodev.py (encode/signatures/anchor_bucket/
** candidate_pairs). The restrictions from §4.2 of the StarcodeV design doc
** apply here too: no runtime RNG (the anchor table comes from a fixed LCG
** constant, computed once), no floating-point reduction, no dependency on
** input order.
*/
#ifndef _PAIRGEN_HEADER
#define _PAIRGEN_HEADER

#include <stdint.h>
#include <stddef.h>

/* ---------------------------------------------------------------- constants */
#define PG_W_ANCHOR   4        /* anchor length (StarcodeV design doc §7) */
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

/* All pairs within a bucket (not just adjacent ones) -- see the calibration
** note in starcodev.py: 14x more effective per iteration than only
** considering adjacent pairs after sorting. */
void pg_candidate_pairs(const uint64_t *keys, int32_t n, int32_t cap,
                         pg_pairlist_t *out);

void pg_pairlist_init(pg_pairlist_t *pl);
void pg_pairlist_push(pg_pairlist_t *pl, int64_t a, int64_t b);
void pg_pairlist_free(pg_pairlist_t *pl);

#endif /* _PAIRGEN_HEADER */
