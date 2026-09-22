# Sequential C Port Conventions: StarcodeV follows the `gui11aume/starcode` structure

This document maps the sequential C port of the Python reference
(`starcodev.py` + `autocal.py`, already validated across 11 synthetic
regimes in `CORRECTNESS_REPORT_GENERALIZATION.md`) onto the file layout
and code style of the original **gui11aume/starcode v1.4** repo
(`src/main-starcode.c`, `starcode.c/h`, `trie.c/h`, `view.c/h`). The goal:
anyone already familiar with the original repo should be able to find the
equivalent piece in this port right away, and anywhere this port has
**no** direct equivalent, that's called out explicitly instead of left as
a gap the reader has to guess at.

## 1. File correspondence table

| Original file (`gui11aume/starcode`) | This port's file | Why |
|---|---|---|
| `src/main-starcode.c` | `main-starcodev.c` | CLI driver: `getopt_long` with a `long_options[]` table, the `USAGE` string, `say_usage()`/`say_version()`, `SIGSEGV_handler()` — the pattern is kept identical. New options (`--tau-core`, `--tau-abs`, `--tau-cons`, `--auto-tau`) are added without changing the meaning of the legacy options that still apply (`-i -o -t -d -c -q -v -h`). |
| `src/starcode.h` | `starcodev.h` | Public declarations for stages 1-4, union-find, `canonical_partition`, and the constants (`TAU_CORE`, `TAU_ABS`, `TAU_CONS`, etc.) — the equivalent of `starcode.h` declaring `int starcode(...)` and `output_t`/`cluster_t`. |
| `src/starcode.c` | `starcodev.c` | Core logic: `s1_cores`, `s2_consensus`, `s3_absorb`, `s4_merge_cons`, `assemble`. The equivalent of `starcode.c`, which orchestrates the trie + clustering calls. |
| `src/trie.c` + `src/trie.h` | **`pairgen.c` + `pairgen.h`** (NOT "trie.c/h") | See §2 — this is a design replacement, not a rename. |
| `src/view.c` + `src/view.h` | `output_writer.c` + `output_writer.h` | Output writer: default mode, `--print-clusters`, `--seq-id`, `--tidy`. The line format is kept identical to original Starcode so it's directly compatible with the Python evaluation harness (`eval_generalisasi.py`, `compare.py`) without touching the parser on the Python side. |
| *(not present in the original repo)* | `autocal.c` + `autocal.h` | A new module: separability-gap detection and `auto_tau()`. Nothing like it exists in Starcode v1.4, because Starcode doesn't do automatic threshold calibration (τ is either given by the user or comes from a different heuristic via `-d auto`). |
| `Makefile` | `Makefile` | A `make seq` target (plain gcc, no CUDA) added alongside this project's existing `make gpu`/`make shim` targets from the CUDA port. |

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
across reads/pairs, which is exactly what makes it a direct port target
for a GPU kernel (`p1_pairgen.cu` in this project's CUDA implementation).

Since there's no 1:1 data-structure correspondence (tree vs. hash table),
this file was given the new name **`pairgen.c/h`** (pair generator)
rather than forcing it to stay `trie.c/h` — the old name would mislead a
reader looking for a tree structure that simply isn't here. `pairgen.c`
takes over the trie's role as the **source of candidate pairs**, but
final verification is still exact edit distance (Levenshtein with
`score_cutoff`), same as original Starcode — precision is guaranteed by
construction; only recall (whether every true edge gets found) is
probabilistic, and that's measured via the iteration sweep
(`N_ITER_S1`, etc. — see the calibration comments in `starcodev.py`).

## 3. §4.2 restrictions binding on the C port (identical to the Python version)

1. **NO RNG** — the anchor-position schedule comes from a constant
   table (`pg_init_anchors()`, a fixed-constant LCG from Numerical
   Recipes), computed once at startup, identical on every machine and
   every run. There is no `rand()`/`srand()` anywhere on a decision
   path.
2. **NO floating-point reduction** — every decision quantity (edit
   distance, Hamming signature, consensus vote count) is `int`/`int64_t`.
   The only `double` anywhere in this port is `CAL_MARGIN = 0.25` in
   `cal_auto_tau()`, used purely as a one-shot multiplier right before
   it's rounded down to an `int` — never as an accumulator.
3. **NO order-dependent accumulator that affects a decision** — every
   decision (medoid, absorption argmin, consensus vote) uses a combined
   key `(value << 32) | index` reduced with `min`, which is
   associative-commutative — the direct mirror of a GPU `atomicMin`.
4. **NO dependency on input order** — the K5 union-find rule always
   picks the smallest-indexed root; the K5 canonical relabeling sorts
   clusters by their lexicographically-smallest member. This is tested
   directly in §7 (the C determinism test under input permutation).

## 4. Canonical tie-breakers K1-K5 (kept exactly as in the Python reference)

| Code | Context | Rule |
|---|---|---|
| K1 | Absorption argmin (S3) | Combined key `(distance << 32) \| consensus_index`; `min` over that key. |
| K2 | Medoid (S2) | Among reads tied for minimum total distance, pick the lexicographically smallest read; if the reads are identical, the smallest global index. |
| K3 | Consensus-column voting (S2) | A substitution-vote tie is broken by the fixed order `A < C < G < T` (smallest `argmax` index). An insertion-vote tie is **dropped** (a conservative policy — never resolved by a random pick). |
| K4 | Medoid subsampling for large cores (S2) | Take the `M_MAX=64` lexicographically-smallest members, not the first `M_MAX` in input order. |
| K5 | Union-find & canonical partition | Union-find root = smallest index; final cluster id = lexicographic rank of the cluster's smallest member. |

## 5. Output format — harness compatibility

`output_writer.c` reproduces original Starcode's `view.c` in exactly the
same four modes, because `eval_generalisasi.py`/`compare.py` parse the
output lines against this format:

- **default**: one line per cluster, `centroid<TAB>size` (in ascending
  canonical label order).
- **`--print-clusters`**: `centroid<TAB>size<TAB>member1,member2,...`.
- **`--seq-id`**: the member column holds 1-based input line numbers
  instead of read strings.
- **`--tidy`**: one line per **input** read (not per cluster):
  `read<TAB>centroid`.

## 6. Validation status

This port was checked for **correctness** (the C canonical partition
compared as a set of groups against the Python canonical partition) on
regimes D1-baseline, D7-depth-3, D9-length-150, D10-length-40 — chosen
to cover the baseline case, low copy depth, a core-coverage failure case
(D9), and a fixed-τ failure case that `cal_auto_tau()` resolves (D10).
Full results in `correctness_c_vs_python.csv`:

| Regime | Reads | Py clusters | C clusters | Identical groups | Reads affected | % different |
|---|---|---|---|---|---|---|
| D1-baseline | 3274 | 500 | 500 | 500/500 | 0 | 0.000% |
| D7-depth-3 | 1500 | 579 | 579 | 579/579 | 0 | 0.000% |
| D9-length-150 | 4000 | 3540 | 3540 | 3540/3540 | 0 | 0.000% |
| D10-length-40 | 3910 | 470 | 471 | 469/470 | 32 | 0.818% |

D1, D7, D9: the C canonical partition is **bit-for-bit identical** to
Python (0 differing groups). D10 has **1 differing group out of
470/471** (32 of 3910 reads, 0.818%), traced to one residual alignment
tie-break — see §6.1.

### 6.1 Where the residual D10 difference comes from, and why it isn't a logic bug

`align_ops()` (the Needleman-Wunsch backtrace inside
`consensus_from_anchor`, S2) hits path ambiguity whenever several edit
paths tie on distance (e.g. a homopolymer deletion `XAAAAY`→`XAAY` can
be explained by two different deletion positions at the same distance).
The Python reference calls
`rapidfuzz.distance.Levenshtein.opcodes()`, computed via a bit-parallel
algorithm (Myers/Hyyrö) — **not** a fixed-priority DP-table backtrace.
Empirical investigation (see the session history) found that a
backtrace priority order of **delete before diagonal (match/replace)
before insert** in `align_ops()`'s backward pass matches `rapidfuzz` on
most of the ambiguous cases, and brought the number of Stage-2 cores
whose consensus differs from Python on D10 down from 19/502 to 15/502
(≈3%) — that fix is applied directly in `align_ops()` (see the comment
in the function body). Since `rapidfuzz` doesn't guarantee that any
single fixed-priority DP backtrace can reproduce its tie-break in EVERY
ambiguous column, the remaining ~3% of cores with a one-base-different
consensus is a **documented limitation**, not a hidden bug.

On D10, one consensus-base difference (core anchor
`CGGTCTCGCACATTAACTCACGTCTACACTCGAAAATTG` vs the C version's
`CGGTCTCGCACATTAACTTCACGTCTACACTCGAAAATTG`, position 16) shifts the
inter-consensus Levenshtein distance right around the Stage 4
`tau_cons=20` threshold, breaking a union-find chain differently from
Python for 5 cores that should merge into 1 component (Python) but
split into 2 components (C). The impact is fully contained to 32 reads
(0.818% of D10), doesn't spread to the other regimes (D1/D7/D9: 0
differences), and doesn't change this report's generalization
conclusions.

## 7. Input-permutation determinism test

Every regime (D1, D7, D9, D10) was tested with **3 random permutations**
of the input lines (`random.shuffle`, seeds 1/2/3) run through the same
C binary. The cluster partition was compared as a set of **read-string**
groups (not line indices, since those shift under permutation) against
the partition on the original ordering. Results
(`determinism_c_permutation.csv`): **all 12 combinations (4 regimes × 3
seeds) match 100%** — 0 differing groups in every comparison. This
confirms restriction §3.4: the smallest-rooted K5 union-find and the
lexicographic canonical relabeling make the final result independent of
input order, regardless of the anchor-hash processing order in Stage
1/3/4.
