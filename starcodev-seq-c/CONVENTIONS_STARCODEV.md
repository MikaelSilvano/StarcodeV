# StarcodeV conventions: following the `gui11aume/starcode` structure

This document maps the sequential C implementation of StarcodeV onto the
file layout and code style of the original **gui11aume/starcode v1.4** repo
(`src/main-starcode.c`, `starcode.c/h`, `trie.c/h`, `view.c/h`). The goal:
anyone already familiar with the original repo should be able to find the
equivalent piece in StarcodeV right away, and anywhere this implementation has
**no** direct equivalent, that's called out explicitly instead of left as
a gap the reader has to guess at.

## 1. File correspondence table

| Original file (`gui11aume/starcode`) | StarcodeV file | Why |
|---|---|---|
| `src/main-starcode.c` | `main-starcodev.c` | CLI driver: `getopt_long` with a `long_options[]` table, the `USAGE` string, `say_usage()`/`say_version()`, `SIGSEGV_handler()` — the pattern is kept identical. StarcodeV's own options (`--tau-core`, `--tau-abs`, `--tau-cons`, `--auto-tau`, `--anchor-occ`, `--audit-sample`, `--exact-max`, `--approx`, ...) are added without changing the meaning of the Starcode options that still apply (`-i -o -t -d -c -q -v -h`). |
| `src/starcode.h` | `starcodev.h` | Public declarations for stages 1-4, union-find, `canonical_partition`, and the constants (`TAU_CORE`, `TAU_ABS`, `TAU_CONS`, etc.) — the equivalent of `starcode.h` declaring `int starcode(...)` and `output_t`/`cluster_t`. |
| `src/starcode.c` | `starcodev.c` | Core logic: `s1_cores`, `s2_consensus`, `s4_absorb` (candidate search), `s3_merge_cons` (filtered merge for `--approx`), `assemble`. The equivalent of `starcode.c`, which orchestrates the trie + clustering calls. |
| `src/trie.c` + `src/trie.h` | **`pairgen.c` + `pairgen.h`** (NOT "trie.c/h") | See §2 — this is a design replacement, not a rename. |
| `src/view.c` + `src/view.h` | `output_writer.c` + `output_writer.h` | Output writer: default mode, `--print-clusters`, `--seq-id`, `--tidy`. The line format is kept identical to original Starcode, so any script that already parses Starcode's output can read StarcodeV's output unchanged. |
| *(not present in the original repo)* | `autocal.c` + `autocal.h` | Separability-gap rule (`cal_ceiling_below_peak`), the filtered histogram used by `--approx --auto-tau`, and the S1 core-coverage diagnostic. Starcode has no counterpart: its τ is either given by the user or comes from a length heuristic. |
| *(not present in the original repo)* | `exactcons.c` + `exactcons.h` | Exact work on the consensus set: bit-parallel edit distance, the exact distance histogram H, the exact S3 merge, the S4 nearest-consensus check and the optional S1 audit. |
| `Makefile` | `Makefile` | A `make seq` target: plain gcc, no CUDA. |

## 2. Why `trie.c/h` was replaced rather than ported

The trie in original Starcode answers one question: **"which reads are
within Levenshtein radius τ of read X?"** — answered with a tree walk
that prunes a branch as soon as the running distance exceeds the
threshold (branch-and-bound over an `A/C/G/T` + `PAD` trie). That
structure is fundamentally sequential, node by node: each step of the
walk depends on the result of the previous step on the same branch, so
there's no direct mapping to independent parallel units of work without
a substantial restructuring.

StarcodeV replaces it with
**anchor-hash + signature candidate-pair generation**: every read is
hashed at the same set of anchor positions (the `ANCHORS` table, a fixed
LCG constant — not a runtime RNG), reads that land in the same hash
bucket get grouped together, and **every pair within a bucket** becomes
a candidate, further filtered by q-gram signature Hamming distance
before exact Levenshtein verification. The unit of work here — "hash one
read at one iteration", "compare one candidate pair" — is independent
across reads/pairs, which also means the work can be spread over
many threads.

Since there's no 1:1 data-structure correspondence (tree vs. hash table),
this file was given the new name **`pairgen.c/h`** (pair generator)
rather than forcing it to stay `trie.c/h` — the old name would mislead a
reader looking for a tree structure that simply isn't here. `pairgen.c`
takes over the trie's role as the **source of candidate pairs**, but
final verification is still exact edit distance (Levenshtein with
`score_cutoff`), same as original Starcode — precision is guaranteed by
construction; only recall (whether every true edge gets found) is
probabilistic, and that's measured via the iteration sweep
(`SV_N_ITER_S1`) and checked with the `--audit-sample` option (see VALIDATION.md).

## 3. Determinism rules

1. **NO RNG** — the anchor-position schedule comes from a constant
   table (`pg_init_anchors()`, a fixed-constant LCG from Numerical
   Recipes), computed once at startup, identical on every machine and
   every run. There is no `rand()`/`srand()` anywhere on a decision
   path.
2. **NO floating-point reduction** — every decision quantity (edit
   distance, Hamming signature, consensus vote count) is `int`/`int64_t`.
   The only `double` on a decision path is the margin (`--margin`, default
   `CAL_MARGIN = 0.25`), used purely as a one-shot multiplier right before
   it's rounded down to an `int` — never as an accumulator.
3. **NO order-dependent accumulator that affects a decision** — every
   decision (medoid, absorption argmin, consensus vote) uses a combined
   key `(value << 32) | index` reduced with `min`, which is
   associative and commutative, so the visiting order never matters (on a GPU it is a single `atomicMin`).
4. **NO dependency on input order** — the K5 union-find rule always
   picks the smallest-indexed root; the K5 canonical relabeling sorts
   clusters by their lexicographically-smallest member. This is tested
   directly in §7 (determinism under input permutation).

## 4. Canonical tie-breakers K1-K5

| Code | Context | Rule |
|---|---|---|
| K1 | Absorption argmin (S4) | Combined key `(distance << 32) \| consensus_index`; `min` over that key. |
| K2 | Medoid (S2) | Among reads tied for minimum total distance, pick the lexicographically smallest read; if the reads are identical, the smallest global index. |
| K3 | Consensus-column voting (S2) | A substitution-vote tie is broken by the fixed order `A < C < G < T` (smallest `argmax` index). An insertion-vote tie is **dropped** (a conservative policy — never resolved by a random pick). |
| K4 | Medoid subsampling for large cores (S2) | Take the `M_MAX=64` lexicographically-smallest members, not the first `M_MAX` in input order. |
| K5 | Union-find & canonical partition | Union-find root = smallest index; final cluster id = lexicographic rank of the cluster's smallest member. |

## 5. Output format

`output_writer.c` reproduces original Starcode's `view.c` in exactly the
same four modes, so the output lines follow the format Starcode users
already know:

- **default**: one line per cluster, `centroid<TAB>size` (in ascending
  canonical label order).
- **`--print-clusters`**: `centroid<TAB>size<TAB>member1,member2,...`.
- **`--seq-id`**: the member column holds 1-based input line numbers
  instead of read strings.
- **`--tidy`**: one line per **input** read (not per cluster):
  `read<TAB>centroid`.

## 6. Tie-breaks in the consensus alignment

`align_ops()` (the backtrace inside `consensus_from_anchor`, S2) has to choose
between edit paths of equal cost, for example a homopolymer deletion
`XAAAAY` -> `XAAY`, where the deletion can sit at any position of the run. The
rule is fixed: delete before diagonal (match/replace) before insert, so the
deletion is always placed at the first position of the run and all members of
a core vote on the same column. The choice is arbitrary, but because it is
fixed the consensus depends on nothing but the reads themselves.

## 7. Determinism under input permutation

Each data set in `validation/determinism.csv` (L40_p05, L60_p10, D3, D10, D11)
was run again on the original input and on 3 shuffled copies of the input
lines. The re-run is byte-identical and the partition of the shuffled inputs
is identical, as a set of read-string groups (line indices shift under
permutation, so they are not compared): 15 of 15 combinations match. This is
rule §3.4 at work: the smallest-rooted union-find and the lexicographic
canonical relabelling make the result independent of input order, whatever
the order in which anchor buckets are processed in S1 and S4.

Results against Starcode are in VALIDATION.md.
