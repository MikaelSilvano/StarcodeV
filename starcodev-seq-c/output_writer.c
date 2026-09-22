/*
** output_writer.c -- see output_writer.h. The output format was checked
** by hand against the real starcode-v1.4 binary (bioconda build), example
** (tau=2, 5 input lines):
**
**   $ starcode -i t.txt -c -d 2 -t 1
**   TTTTGGGGCCCCAAAA\t6
**   AAAACCCCGGGGTTTT\t5
**   GATTACAGATTACAGA\t4
**
**   $ starcode -i t.txt -c -d 2 --print-clusters -t 1
**   TTTTGGGGCCCCAAAA\t6\tTTTTGGGGCCCCAAAA,TTTTGGGGCCCCAAAT
**   AAAACCCCGGGGTTTT\t5\tAAAACCCCGGGGTTTT,AAAACCCCGGGGTTTA
**   GATTACAGATTACAGA\t4\tGATTACAGATTACAGA
**
**   $ starcode -i t.txt -c -d 2 --print-clusters --seq-id -t 1
**   TTTTGGGGCCCCAAAA\t6\tTTTTGGGGCCCCAAAA,TTTTGGGGCCCCAAAT\t3,4
**   ...
**
**   $ starcode -i t.txt -d 2 --tidy -t 1   (default mode: message passing)
**   AAAACCCCGGGGTTTT\tAAAACCCCGGGGTTTT
**   AAAACCCCGGGGTTTA\tAAAACCCCGGGGTTTA
**   TTTTGGGGCCCCAAAA\tTTTTGGGGCCCCAAAA
**   TTTTGGGGCCCCAAAT\tTTTTGGGGCCCCAAAA
**   GATTACAGATTACAGA\tGATTACAGATTACAGA
**
** StarcodeV only has a connected-components-equivalent algorithm
** (union-find), so we reproduce the lines above exactly for
** default/print-clusters/seq-id mode. StarcodeV's --tidy is driven by the
** same S1-S4 assembled labels (not message-passing like that last example
** above) -- the line format (read<TAB>centroid) stays identical, only the
** underlying clustering values come from a different algorithm. This is
** documented in CONVENTIONS_STARCODEV.md §5.
*/
#include "output_writer.h"
#include "starcodev.h"

#include <stdlib.h>
#include <string.h>

#define SV_OUT_LINE_CAP 1048576  /* 1 MiB per line -- enough for any reasonably long read */

int sv_out_read_input(FILE *in, sv_out_rawinput_t *raw_out) {
    char *buf = (char *)malloc(SV_OUT_LINE_CAP);
    size_t cap = 4096;
    char **seqs = (char **)malloc(cap * sizeof(char *));
    int32_t *counts = (int32_t *)malloc(cap * sizeof(int32_t));
    int32_t n = 0;

    while (fgets(buf, SV_OUT_LINE_CAP, in) != NULL) {
        size_t len = strlen(buf);
        while (len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r')) buf[--len] = '\0';
        if (len == 0) continue;  /* skip blank lines (matches original Starcode's behavior) */

        char *tab = strchr(buf, '\t');
        int32_t count = 1;
        size_t seqlen = len;
        if (tab != NULL) {
            *tab = '\0';
            seqlen = (size_t)(tab - buf);
            count = atoi(tab + 1);
            if (count <= 0) count = 1;
        }
        if (n >= (int32_t)cap) {
            cap *= 2;
            seqs = (char **)realloc(seqs, cap * sizeof(char *));
            counts = (int32_t *)realloc(counts, cap * sizeof(int32_t));
        }
        seqs[n] = (char *)malloc(seqlen + 1);
        memcpy(seqs[n], buf, seqlen + 1);
        counts[n] = count;
        n++;
    }
    free(buf);
    raw_out->seqs = seqs;
    raw_out->counts = counts;
    raw_out->n_lines = n;
    return n;
}

void sv_out_free_raw(sv_out_rawinput_t *raw) {
    for (int32_t i = 0; i < raw->n_lines; i++) free(raw->seqs[i]);
    free(raw->seqs); free(raw->counts);
    raw->seqs = NULL; raw->counts = NULL; raw->n_lines = 0;
}

typedef struct { const char *seq; int32_t line_idx; } sv_out_kv_t;
static int sv_out_kv_cmp(const void *pa, const void *pb) {
    const sv_out_kv_t *a = (const sv_out_kv_t *)pa, *b = (const sv_out_kv_t *)pb;
    int c = strcmp(a->seq, b->seq);
    if (c != 0) return c;
    return a->line_idx < b->line_idx ? -1 : (a->line_idx > b->line_idx ? 1 : 0);
}

void sv_out_build_unique(const sv_out_rawinput_t *raw, sv_out_uniqinput_t *out) {
    int32_t nl = raw->n_lines;
    sv_out_kv_t *kv = (sv_out_kv_t *)malloc((size_t)nl * sizeof(sv_out_kv_t));
    for (int32_t i = 0; i < nl; i++) { kv[i].seq = raw->seqs[i]; kv[i].line_idx = i; }
    qsort(kv, (size_t)nl, sizeof(sv_out_kv_t), sv_out_kv_cmp);

    /* count the number of distinct groups */
    int32_t nu = 0;
    for (int32_t i = 0; i < nl; i++) if (i == 0 || strcmp(kv[i].seq, kv[i - 1].seq) != 0) nu++;

    char **useqs = (char **)malloc((size_t)nu * sizeof(char *));
    int32_t *ulens = (int32_t *)malloc((size_t)nu * sizeof(int32_t));
    int64_t *ucounts = (int64_t *)calloc((size_t)nu, sizeof(int64_t));
    int32_t **uids = (int32_t **)malloc((size_t)nu * sizeof(int32_t *));
    int32_t *unids = (int32_t *)calloc((size_t)nu, sizeof(int32_t));

    int32_t g = -1;
    for (int32_t i = 0; i < nl; i++) {
        if (i == 0 || strcmp(kv[i].seq, kv[i - 1].seq) != 0) {
            g++;
            useqs[g] = strdup(kv[i].seq);
            ulens[g] = (int32_t)strlen(kv[i].seq);
        }
        unids[g]++;
    }
    for (int32_t k = 0; k < nu; k++) uids[k] = (int32_t *)malloc((size_t)unids[k] * sizeof(int32_t));
    int32_t *fillpos = (int32_t *)calloc((size_t)nu, sizeof(int32_t));
    g = -1;
    for (int32_t i = 0; i < nl; i++) {
        if (i == 0 || strcmp(kv[i].seq, kv[i - 1].seq) != 0) g++;
        uids[g][fillpos[g]++] = kv[i].line_idx + 1;  /* 1-based, ascending original file order */
        ucounts[g] += raw->counts[kv[i].line_idx];
    }
    free(fillpos);
    free(kv);

    out->seqs = useqs; out->lens = ulens; out->counts = ucounts;
    out->ids = uids; out->nids = unids; out->n = nu;
}

void sv_out_free_unique(sv_out_uniqinput_t *u) {
    for (int32_t i = 0; i < u->n; i++) { free(u->seqs[i]); free(u->ids[i]); }
    free(u->seqs); free(u->lens); free(u->counts); free(u->ids); free(u->nids);
    u->seqs = NULL; u->lens = NULL; u->counts = NULL; u->ids = NULL; u->nids = NULL; u->n = 0;
}

/* ------------------------------------------------------------------ */
typedef struct { int32_t label; int32_t *members; int32_t nmem; char *rep; } sv_out_cluster_t;

static int sv_out_cluster_cmp(const void *pa, const void *pb) {
    const sv_out_cluster_t *a = (const sv_out_cluster_t *)pa, *b = (const sv_out_cluster_t *)pb;
    return strcmp(a->rep, b->rep);
}
static int sv_out_member_cmp_ctx(const int32_t *a, const int32_t *b, char *const *seqs) {
    return strcmp(seqs[*a], seqs[*b]);
}
static char *const *g_sv_out_seqs;
static int sv_out_member_cmp(const void *pa, const void *pb) {
    return sv_out_member_cmp_ctx((const int32_t *)pa, (const int32_t *)pb, g_sv_out_seqs);
}

void sv_out_write(FILE *out, const sv_out_uniqinput_t *u, const int32_t *labels,
                   const sv_out_rawinput_t *raw, const sv_out_opts_t *opts) {
    int32_t n = u->n;

    if (opts->tidy) {
        /* one line per UNIQUE read -- repeated for each raw input line that
        ** maps to it, in the original input line order (matches original
        ** Starcode's --tidy behavior: one line per input, not per distinct
        ** read). */
        /* rep[label] = canonical representative string (lexicographically smallest). */
        int32_t maxlab = -1;
        for (int32_t i = 0; i < n; i++) if (labels[i] > maxlab) maxlab = labels[i];
        char **rep = (char **)calloc((size_t)(maxlab + 1), sizeof(char *));
        for (int32_t i = 0; i < n; i++) {
            int32_t l = labels[i];
            if (rep[l] == NULL || strcmp(u->seqs[i], rep[l]) < 0) rep[l] = u->seqs[i];
        }
        /* map raw line id -> distinct-read index */
        int32_t *line_to_uniq = (int32_t *)malloc((size_t)raw->n_lines * sizeof(int32_t));
        for (int32_t i = 0; i < n; i++)
            for (int32_t k = 0; k < u->nids[i]; k++)
                line_to_uniq[u->ids[i][k] - 1] = i;
        for (int32_t li = 0; li < raw->n_lines; li++) {
            int32_t ui = line_to_uniq[li];
            fprintf(out, "%s\t%s\n", raw->seqs[li], rep[labels[ui]]);
        }
        free(line_to_uniq); free(rep);
        return;
    }

    /* cluster mode (default / --print-clusters / --seq-id) */
    int32_t maxlab = -1;
    for (int32_t i = 0; i < n; i++) if (labels[i] > maxlab) maxlab = labels[i];
    int32_t nlab = maxlab + 1;
    int32_t *cnt = (int32_t *)calloc((size_t)nlab, sizeof(int32_t));
    for (int32_t i = 0; i < n; i++) cnt[labels[i]]++;

    sv_out_cluster_t *cl = (sv_out_cluster_t *)malloc((size_t)nlab * sizeof(sv_out_cluster_t));
    for (int32_t l = 0; l < nlab; l++) {
        cl[l].label = l;
        cl[l].members = (int32_t *)malloc((size_t)cnt[l] * sizeof(int32_t));
        cl[l].nmem = 0;
        cl[l].rep = NULL;
    }
    for (int32_t i = 0; i < n; i++) {
        sv_out_cluster_t *c = &cl[labels[i]];
        c->members[c->nmem++] = i;
    }
    g_sv_out_seqs = u->seqs;
    for (int32_t l = 0; l < nlab; l++) {
        qsort(cl[l].members, (size_t)cl[l].nmem, sizeof(int32_t), sv_out_member_cmp);
        cl[l].rep = u->seqs[cl[l].members[0]];  /* lexicographically smallest (K5) */
    }
    qsort(cl, (size_t)nlab, sizeof(sv_out_cluster_t), sv_out_cluster_cmp);

    for (int32_t l = 0; l < nlab; l++) {
        sv_out_cluster_t *c = &cl[l];
        int64_t total = 0;
        for (int32_t k = 0; k < c->nmem; k++) total += u->counts[c->members[k]];
        fprintf(out, "%s\t%lld", c->rep, (long long)total);

        if (opts->show_clusters) {
            fprintf(out, "\t%s", u->seqs[c->members[0]]);
            for (int32_t k = 1; k < c->nmem; k++) fprintf(out, ",%s", u->seqs[c->members[k]]);
        } else if (opts->show_seqid) {
            fprintf(out, "\t%s", c->rep);
        }

        if (opts->show_seqid) {
            /* Gather every member's ids across the whole cluster, sort them
            ** ASCENDING -- mirrors original Starcode's sort_and_print_ids()
            ** (the id stack is sorted before printing, regardless of the
            ** member traversal order). */
            int32_t total_ids = 0;
            for (int32_t k = 0; k < c->nmem; k++) total_ids += u->nids[c->members[k]];
            int32_t *allids = (int32_t *)malloc((size_t)total_ids * sizeof(int32_t));
            int32_t p = 0;
            for (int32_t k = 0; k < c->nmem; k++) {
                int32_t ui = c->members[k];
                for (int32_t j = 0; j < u->nids[ui]; j++) allids[p++] = u->ids[ui][j];
            }
            /* insertion sort -- the number of ids per cluster is small in practice */
            for (int32_t i = 1; i < total_ids; i++) {
                int32_t key = allids[i], j = i - 1;
                while (j >= 0 && allids[j] > key) { allids[j + 1] = allids[j]; j--; }
                allids[j + 1] = key;
            }
            fputc('\t', out);
            for (int32_t i = 0; i < total_ids; i++) fprintf(out, "%s%d", i ? "," : "", allids[i]);
            free(allids);
        }
        fputc('\n', out);
    }

    for (int32_t l = 0; l < nlab; l++) free(cl[l].members);
    free(cl); free(cnt);
}
