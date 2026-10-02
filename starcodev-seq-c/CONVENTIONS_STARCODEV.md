# Sequential C Conventions: StarcodeV follows the `gui11aume/starcode` structure

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
| `src/main-starcode.c` | `main-starcodev.c` | CLI driver: `getopt_long` with a `long_options[]` table, the `USAGE` string, `say_usage()`/`say_version()`, `SIGSEGV_handler()` — the pattern is kept identical. New options (`--tau-core`, `--tau-abs`, `--tau-cons`, `--auto-tau`) are added without changing the meaning of the legacy options that still apply (`-i -o -t -d -c -q -v -h`). |
| `src/starcode.h` | `starcodev.h` | Public declarations for stages 1-4, union-find, `canonical_partition`, and the constants (`TAU_CORE`, `TAU_ABS`, `TAU_CONS`, etc.) — the equivalent of `starcode.h` declaring `int starcode(...)` and `output_t`/`cluster_t`. |
| `src/starcode.c` | `starcodev.c` | Core logic: `s1_cores`, `s2_consensus`, `s3_absorb`, `s4_merge_cons`, `assemble`. The equivalent of `starcode.c`, which orchestrates the trie + clustering calls. |
| `src/trie.c` + `src/trie.h` | **`pairgen.c` + `pairgen.h`** (NOT "trie.c/h") | See §2 — this is a design replacement, not a rename. |
| `src/view.c` + `src/view.h` | `output_writer.c` + `output_writer.h` | Output writer: default mode, `--print-clusters`, `--seq-id`, `--tidy`. The line format is kept identical to original Starcode, so any script that already parses Starcode's output can read StarcodeV's output unchanged. |
| *(not present in the original repo)* | `autocal.c` + `autocal.h` | A new module: separability-gap detection and `auto_tau()`. Nothing like it exists in Starcode v1.4, because Starcode doesn't do automatic threshold calibration (τ is either given by the user or comes from a different heuristic via `-d auto`). |
| *(not present in the original repo)* | `exactcons.c` + `exactcons.h` | Added in v2. Exact work on the consensus set: bit-parallel edit distance, the exact distance histogram, S4 without the candidate filter, the S3 nearest-consensus check and the optional S1 audit. See CHANGES_v2.md. |
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

The StarcodeV design (see the design document §4.6/§7) replaces it with
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
(`N_ITER_S1`, etc. — see VALIDATION_v2.md and the `S1 audit` option).

## 3. §4.2 restrictions binding on this implementation

1. **NO RNG** — the anchor-position schedule comes from a constant
   table (`pg_init_anchors()`, a fixed-constant LCG from Numerical
   Recipes), computed once at startup, identical on every machine and
   every run. There is no `rand()`/`srand()` anywhere on a decision
   path.
2. **NO floating-point reduction** — every decision quantity (edit
   distance, Hamming signature, consensus vote count) is `int`/`int64_t`.
   The only `double` anywhere in this code is `CAL_MARGIN = 0.25` (default of `--margin`, fixed in v1.1; see AUTOTAU_REVISION.md) in
   `cal_auto_tau()`, used purely as a one-shot multiplier right before
   it's rounded down to an `int` — never as an accumulator.
3. **NO order-dependent accumulator that affects a decision** — every
   decision (medoid, absorption argmin, consensus vote) uses a combined
   key `(value << 32) | index` reduced with `min`, which is
   associative-commutative — the direct mirror of a GPU `atomicMin`.
4. **NO dependency on input order** — the K5 union-find rule always
   picks the smallest-indexed root; the K5 canonical relabeling sorts
   clusters by their lexicographically-smallest member. This is tested
   directly in §7 (determinism under input permutation).

## 4. Canonical tie-breakers K1-K5

| Code | Context | Rule |
|---|---|---|
| K1 | Absorption argmin (S3) | Combined key `(distance << 32) \| consensus_index`; `min` over that key. |
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

## 6. Validation status: StarcodeV against Starcode

The check that matters is whether StarcodeV keeps everything Starcode got
right and only adds to it. Two partitions are compared with an overlap graph
(one node per cluster of each partition, an edge when two clusters share a
read). Every connected component is *identical*, a *merge* (several Starcode
clusters inside one StarcodeV cluster), a *split* (one Starcode cluster divided
among several StarcodeV clusters) or *crossing* (reads moved between groups).
Merge and identical are the intended outcomes; a split or a crossing would mean
StarcodeV broke something Starcode had kept together.

Microsoft clustered-nanopore data (269,709 reads, 9,984 oligos), Starcode v1.4
at -d 8 against StarcodeV v2 default:

| | Clusters | Purity | Whole oligos | ARI | Recovery γ = 0.5 | Strict γ = 0.9 |
|---|---|---|---|---|---|---|
| Starcode -c | 74,435 | 100% | 129 | 0.769 | 89.73% | 5.93% |
| StarcodeV v2 | 11,163 | 100% | 9,491 | 0.999 | 98.08% | 97.63% |

Overlap graph, Starcode -c against StarcodeV: 1,503 identical, 9,660 merge,
0 split, 0 crossing. The same holds for Starcode -s and message passing. No
Starcode cluster is divided and no read crosses into another group, so every
oligo Starcode recovers is recovered by StarcodeV too.

On the synthetic sets (seed 7 and seed 2026, 15 data sets in all) StarcodeV v2
shows 0 mixed clusters and 0 crossing components against Starcode -c. It does
leave 9 split components, all in two regimes with very shallow or very noisy
reads (D3 and D7), where the S1 filters miss some true d <= 8 edges; the split
direction is the intended failure mode. VALIDATION_v2.md has the full tables.

### 6.1 Tie-breaks in the consensus alignment

`align_ops()` (the backtrace inside `consensus_from_anchor`, S2) has to choose
between edit paths of equal cost, for example a homopolymer deletion
`XAAAAY` -> `XAAY`, where the deletion can sit at any position of the run. The
rule is fixed: delete before diagonal (match/replace) before insert, so the
deletion is always placed at the first position of the run and all members of
a core vote on the same column. The choice is arbitrary, but because it is
fixed the consensus does not depend on anything but the reads themselves.

## 7. Input-permutation determinism test

Each data set in `validation/determinism.csv` (L40_p05, L60_p10, D3, D10, D11)
was run again on the original input and on 3 shuffled copies of the input
lines. The re-run is byte-identical and the partition of the shuffled inputs
is identical, as a set of read-string groups (line indices shift under
permutation, so they are not compared): 15 of 15 combinations match. This
confirms restriction §3.4: the smallest-rooted union-find and the
lexicographic canonical relabelling make the final result independent of input
order, whatever the order in which anchor buckets are processed in S1 and S3.
