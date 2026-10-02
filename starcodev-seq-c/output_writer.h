/*
** output_writer.h -- StarcodeV's output writer, format-compatible with
** original Starcode (gui11aume/starcode v1.4, -c / connected-components
** mode) so scripts that already parse Starcode's output can read this
** binary's output unchanged.
**
** Mirrors the role of src/view.c/.h from original Starcode -- see
** CONVENTIONS_STARCODEV.md §1 and §5. Four modes, each checked directly
** against the starcode-v1.4 binary's actual output (see the manual test
** notes at the bottom of output_writer.c):
**   default            centroid<TAB>size
**   --print-clusters   centroid<TAB>size<TAB>centroid,member2,...
**   --seq-id           adds a 1-based id column (or, without
**                       --print-clusters: centroid<TAB>size<TAB>centroid<TAB>id1,id2,...)
**   --tidy             one line per INPUT read (not per cluster): read<TAB>centroid
**                       (cannot be combined with --print-clusters/--seq-id --
**                       matches the original binary's own validation)
**
** "size" here means the TOTAL COUNT (sum of counts) of a cluster's
** members, following original Starcode's canonical->count semantics --
** NOT the number of distinct reads.
*/
#ifndef _OUTPUT_WRITER_HEADER
#define _OUTPUT_WRITER_HEADER

#include <stdint.h>
#include <stdio.h>

/* Raw input rows, in their original file order (needed for --seq-id and
** --tidy, both of which must report the original input line numbers/order,
** not our internal lexicographic ordering). */
typedef struct {
    char    **seqs;      /* n_lines, the raw read string per input line */
    int32_t  *counts;    /* n_lines, count column (defaults to 1 when absent) */
    int32_t   n_lines;
} sv_out_rawinput_t;

/* Distinct reads, sorted lexicographically, with a list of 1-based raw line ids mapped to
** each distinct read. */
typedef struct {
    char    **seqs;         /* n, allocated (owned by this struct) */
    int32_t  *lens;         /* n */
    int64_t  *counts;       /* n, total count (sum over all identical raw lines) */
    int32_t **ids;          /* n, each element: a list of 1-based ids, ids[i][0..nids[i]-1] */
    int32_t  *nids;         /* n */
    int32_t   n;
} sv_out_uniqinput_t;

int  sv_out_read_input(FILE *in, sv_out_rawinput_t *raw_out);
void sv_out_free_raw(sv_out_rawinput_t *raw);

/* Build the deduplicated representation from raw input (dedup + count
** aggregation + id collection). */
void sv_out_build_unique(const sv_out_rawinput_t *raw, sv_out_uniqinput_t *out);
void sv_out_free_unique(sv_out_uniqinput_t *u);

typedef struct {
    int show_clusters;   /* --print-clusters */
    int show_seqid;      /* --seq-id */
    int tidy;            /* --tidy */
} sv_out_opts_t;

/* labels: sized u->n, canonical label (from sv_canonical_partition) per
** distinct read u->seqs[i]. raw is only needed for --tidy mode (to recover
** the original input line order). Writes to `out`. */
void sv_out_write(FILE *out, const sv_out_uniqinput_t *u, const int32_t *labels,
                   const sv_out_rawinput_t *raw, const sv_out_opts_t *opts);

#endif /* _OUTPUT_WRITER_HEADER */
