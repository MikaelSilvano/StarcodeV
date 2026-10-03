# StarcodeV (sequential C)

StarcodeV clusters noisy DNA-storage reads by the oligo they came from. It
keeps Starcode's exact distance-8 read graph as its first stage and then works
on one consensus sequence per core, which carries far less noise than a single
read, so that larger thresholds can be applied safely.

## Build and run

```sh
make seq                    # builds ./starcodev_seq (gcc, C99, no other dependencies)
make test                   # 5-line smoke test
./starcodev_seq -i reads.txt -o clusters.txt --print-clusters
```

Input: one read per line, optionally `read<TAB>count`, as for Starcode.
Output: the Starcode formats (default, `--print-clusters`, `--seq-id`,
`--tidy`), see CONVENTIONS_STARCODEV.md §5.

## Pipeline

| Stage | What it does | Where |
|---|---|---|
| S1 exact cores | connected components of the read graph at d <= tau_core = 8; candidate pairs from an anchor hash and a q-gram signature filter, every edge checked with the exact distance | `sv_s1_cores`, `pairgen.c` |
| S2 consensus | per core: medoid as the frame, column vote over the aligned members (bases, deletions, insertions; strict majority) | `sv_s2_consensus` |
| H exact distances | exact distance between every pair of consensus sequences (bit-parallel Myers); the gap below the histogram peak gives the ceiling C | `ec_all_pairs`, `cal_ceiling_below_peak` |
| S3 consensus merge | union-find over all consensus pairs with d <= tau_cons = 20, on the exact distances | `ec_s3_merge` |
| S4 orphan absorption | a read left alone by S1 joins the cluster of its nearest consensus within tau_abs = 16; filtered search first, then an exact check unless 2a < c_post proves the pick | `sv_s4_absorb`, `ec_s4_resolve` |

Thresholds at or above the measured ceiling C would join different oligos, so
tau_abs and tau_cons are lowered to floor(0.75 (C - 1)) when they are not below
C. With `--auto-tau` both are set to that value.

Every decision that can join two oligos uses an exact distance. A case that
cannot be settled (no consensus within tau_abs, or two components equally
close) leaves the read on its own: a split costs one extra cluster, a wrong
merge corrupts a cluster.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--tau-core` (`-d`) | 8 | S1 threshold |
| `--tau-cons` | 20 | S3 threshold |
| `--tau-abs` | 16 | S4 threshold |
| `--auto-tau` | off | set tau_cons = tau_abs = floor((C - 1)(1 - margin)) |
| `--margin` | 0.25 | margin used by `--auto-tau` |
| `--auto-tau-hist FILE` | | write the exact consensus histogram (TSV: d, pairs) |
| `--sig-abs` | 96 | signature filter of the S3/S4 candidate searches |
| `--anchor-occ` | 4 | anchor occurrences per read in the S1/S4 candidate search |
| `--audit-sample N` | 0 | check N reads exactly against S1 and report missed edges (no labels needed) |
| `--exact-max` | 50000 | above this many cores, S3/S4 use the filtered searches only |
| `--approx` | off | filtered S3/S4 only, one anchor occurrence per read; faster, but without the exact checks |
| `-q` | | quiet; otherwise the stage log goes to stderr |

`-c` is accepted and changes nothing (StarcodeV always uses connected
components); `-s` and `-r` are rejected because StarcodeV has no sphere or
message-passing mode.

## Files

| File | Content |
|---|---|
| `main-starcodev.c` | command-line driver and pipeline |
| `starcodev.c/.h` | union-find, canonical partition, S1, S2, filtered S3/S4 searches, assembly |
| `pairgen.c/.h` | encoding, anchor keys, signatures, candidate pairs |
| `exactcons.c/.h` | exact distance, H, exact S3, S4 check, S1 audit |
| `autocal.c/.h` | gap rule, filtered histogram for `--approx --auto-tau`, core coverage |
| `output_writer.c/.h` | input reading, de-duplication, Starcode output formats |
| `CONVENTIONS_STARCODEV.md` | mapping to the Starcode source tree, determinism rules, tie-breakers |
| `VALIDATION.md`, `validation/` | results on the Microsoft data and 15 synthetic sets |
