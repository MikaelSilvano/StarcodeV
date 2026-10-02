/*
** main-starcodev.c -- StarcodeV's CLI driver, using the same getopt_long
** style as src/main-starcode.c from original Starcode (gui11aume/starcode
** v1.4). See CONVENTIONS_STARCODEV.md §1: the legacy options (-i -o -t -d
** -c -q -v -h, --print-clusters, --seq-id, --tidy) keep their original
** meaning and names; the new options (--tau-core, --tau-abs, --tau-cons,
** --auto-tau, --sig-abs, --margin, and in v2 --anchor-occ, --exact-max,
** --audit-sample, --legacy) are added on top without clobbering them.
**
** v2 pipeline order: S1 -> S2 -> H (exact histogram) -> S4 -> S3. S4 used to
** run last; it doesn't read the orphans, so moving it up changes nothing by
** itself, but it lets S3 check its picks against the S4 components.
**
** StarcodeV does not implement original Starcode's message-passing (-r)
** or sphere (-s) clustering -- its only underlying algorithm is
** connected-components via union-find (S1+S4), so -c/--connected-comp is
** accepted for argument compatibility but doesn't change behavior (it's
** always connected-components). -s/--sphere and -r/--cluster-ratio are
** rejected with an explicit error message rather than silently ignored.
*/
#include <getopt.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "starcodev.h"
#include "pairgen.h"
#include "autocal.h"
#include "exactcons.h"
#include "output_writer.h"

#define ERRM "starcodev error:"
#define STARCODEV_MAX_TAU 32

char *USAGE =
"\n"
"Usage:  starcodev [options]\n"
"\n"
"  general options:\n"
"    -t --threads: number of concurrent threads (unused, this build is sequential; accepted for compatibility)\n"
"    -q --quiet: quiet output (default verbose)\n"
"    -v --version: display version and exit\n"
"    -h --help: display this help and exit\n"
"\n"
"  cluster options: (algorithm: StarcodeV's 4-stage connected components)\n"
"    -c --connected-comp: accepted for compatibility (always on)\n"
"       --tau-core: Levenshtein threshold for Stage 1 core formation (default 8)\n"
"       --tau-abs: Stage 3 orphan-absorption threshold (default 16, or automatic with --auto-tau)\n"
"       --tau-cons: Stage 4 consensus-merge threshold (default 20, or automatic with --auto-tau)\n"
"       --auto-tau: pick tau-abs & tau-cons automatically from the separability gap (see autocal.h)\n"
"       --sig-abs: Hamming signature pre-filter threshold for Stage 3/4 (default 96)\n"
"       --margin: safety margin for --auto-tau, tau = floor((ceiling-1)*(1-margin)) (default 0.25)\n"
"       --auto-tau-hist: file to write the consensus distance histogram to (TSV: d, pairs)\n"
"\n"
"  v2 options\n"
"       --anchor-occ: anchor occurrences used per read in S1/S3 candidate search (default 4, 1 = v1)\n"
"       --exact-max: largest number of cores handled with exact S3/S4 (default 50000)\n"
"       --audit-sample: reads checked exactly against S1 (default 0 = off)\n"
"       --legacy: run the v1.1 pipeline unchanged (filtered S3/S4, single anchor occurrence)\n"
"\n"
"  input/output options (single file, default)\n"
"    -i --input: input file (default stdin)\n"
"    -o --output: output file (default stdout)\n"
"\n"
"  output format options\n"
"       --print-clusters: outputs cluster compositions\n"
"       --seq-id: print sequence id numbers (1-based)\n"
"       --tidy: print each sequence and its centroid\n"
"\n";

#define VERSION "starcodev-v2.0 (sequential C)"

void say_usage(void) { fprintf(stderr, "%s\n", USAGE); }
void say_version(void) { fprintf(stderr, VERSION "\n"); }

void SIGSEGV_handler(int sig) {
    fprintf(stderr, "Error: signal %d (starcodev)\n", sig);
    exit(1);
}

int main(int argc, char **argv) {
    signal(SIGSEGV, SIGSEGV_handler);

    static int cl_flag = 0;   /* --print-clusters */
    static int id_flag = 0;   /* --seq-id */
    static int tidy_flag = 0; /* --tidy */
    static int vb_flag = 1;   /* verbose, turned off with -q */
    static int cp_flag = 0;   /* -c, accepted, doesn't change behavior */
    static int auto_flag = 0; /* --auto-tau */

    int tau_core = SV_TAU_CORE;
    int tau_abs = SV_TAU_ABS;
    int tau_cons = SV_TAU_CONS;
    int sig_abs = SV_SIG_ABS;
    double margin = CAL_MARGIN;
    char *hist_path = NULL;
    static int legacy_flag = 0;  /* --legacy */
    int anchor_occ = 4;
    int exact_max = 50000;
    int audit_sample = 0;
    int threads = -1;  /* accepted, ignored (this build is sequential) */
    (void)threads;

    char *const UNSET = "unset";
    char *input = UNSET;
    char *output = UNSET;

    if (argc == 1 && isatty(0)) {
        say_usage();
        return EXIT_SUCCESS;
    }

    int c;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"print-clusters", no_argument,       &cl_flag,   1},
            {"seq-id",         no_argument,       &id_flag,   1},
            {"tidy",           no_argument,       &tidy_flag, 1},
            {"quiet",          no_argument,       &vb_flag,   0},
            {"connected-comp", no_argument,       &cp_flag,   1},
            {"auto-tau",       no_argument,       &auto_flag, 1},
            {"version",        no_argument,             0, 'v'},
            {"help",           no_argument,             0, 'h'},
            {"input",          required_argument,       0, 'i'},
            {"output",         required_argument,       0, 'o'},
            {"threads",        required_argument,       0, 't'},
            {"tau-core",       required_argument,       0, 1001},
            {"tau-abs",        required_argument,       0, 1002},
            {"tau-cons",       required_argument,       0, 1003},
            {"sig-abs",        required_argument,       0, 1004},
            {"margin",         required_argument,       0, 1005},
            {"auto-tau-hist",  required_argument,       0, 1006},
            {"legacy",         no_argument,       &legacy_flag, 1},
            {"anchor-occ",     required_argument,       0, 1007},
            {"exact-max",      required_argument,       0, 1008},
            {"audit-sample",   required_argument,       0, 1009},
            {"dist",           required_argument,       0, 'd'},  /* alias for --tau-core, for compatibility */
            {"sphere",         no_argument,             0, 's'},
            {"cluster-ratio",  required_argument,       0, 'r'},
            {0, 0, 0, 0}
        };

        c = getopt_long(argc, argv, "hi:o:qct:v", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
        case 0:
            break;  /* the flag was already set directly by getopt_long */
        case 'v':
            say_version();
            return EXIT_SUCCESS;
        case 'h':
            say_version();
            say_usage();
            return EXIT_SUCCESS;
        case 'i':
            if (input == UNSET) input = optarg;
            else { fprintf(stderr, "%s --input set more than once\n", ERRM); say_usage(); return EXIT_FAILURE; }
            break;
        case 'o':
            if (output == UNSET) output = optarg;
            else { fprintf(stderr, "%s --output set more than once\n", ERRM); say_usage(); return EXIT_FAILURE; }
            break;
        case 'q':
            vb_flag = 0;
            break;
        case 'c':
            cp_flag = 1;
            break;
        case 't':
            threads = atoi(optarg);
            break;
        case 'd':
        case 1001:
            tau_core = atoi(optarg);
            if (tau_core < 0 || tau_core > STARCODEV_MAX_TAU) {
                fprintf(stderr, "%s --tau-core/--dist must be in [0,%d]\n", ERRM, STARCODEV_MAX_TAU);
                return EXIT_FAILURE;
            }
            break;
        case 1002:
            tau_abs = atoi(optarg);
            break;
        case 1003:
            tau_cons = atoi(optarg);
            break;
        case 1004:
            sig_abs = atoi(optarg);
            break;
        case 1005:
            margin = atof(optarg);
            if (margin < 0.0 || margin >= 1.0) {
                fprintf(stderr, "%s --margin must be in [0,1)\n", ERRM);
                return EXIT_FAILURE;
            }
            break;
        case 1006:
            hist_path = optarg;
            break;
        case 1007:
            anchor_occ = atoi(optarg);
            if (anchor_occ < 1 || anchor_occ > 64) {
                fprintf(stderr, "%s --anchor-occ must be in [1,64]\n", ERRM);
                return EXIT_FAILURE;
            }
            break;
        case 1008:
            exact_max = atoi(optarg);
            break;
        case 1009:
            audit_sample = atoi(optarg);
            break;
        case 's':
            fprintf(stderr, "%s --sphere is not supported by StarcodeV "
                    "(connected-components only)\n", ERRM);
            return EXIT_FAILURE;
        case 'r':
            fprintf(stderr, "%s --cluster-ratio is not supported by StarcodeV "
                    "(there is no message-passing algorithm)\n", ERRM);
            return EXIT_FAILURE;
        default:
            say_usage();
            return EXIT_FAILURE;
        }
    }

    if (optind < argc) {
        if (optind == argc - 1 && input == UNSET) {
            input = argv[optind];
        } else {
            fprintf(stderr, "%s too many options\n", ERRM);
            say_usage();
            return EXIT_FAILURE;
        }
    }

    if (tidy_flag && (cl_flag || id_flag)) {
        fprintf(stderr, "%s --tidy flag is not compatible with options "
                "--print-clusters and --seq-id\n", ERRM);
        say_usage();
        return EXIT_FAILURE;
    }

    FILE *fin = stdin;
    if (input != UNSET) {
        fin = fopen(input, "r");
        if (!fin) { fprintf(stderr, "%s cannot open input file %s\n", ERRM, input); return EXIT_FAILURE; }
    }
    FILE *fout = stdout;
    if (output != UNSET) {
        fout = fopen(output, "w");
        if (!fout) { fprintf(stderr, "%s cannot open output file %s\n", ERRM, output); return EXIT_FAILURE; }
    }

    if (vb_flag) fprintf(stderr, "running %s (sequential)\n", VERSION);

    /* ---------------------------------------------------------- read input */
    sv_out_rawinput_t raw;
    sv_out_read_input(fin, &raw);
    if (fin != stdin) fclose(fin);
    if (raw.n_lines == 0) {
        fprintf(stderr, "%s no input sequences\n", ERRM);
        return EXIT_FAILURE;
    }

    sv_out_uniqinput_t uniq;
    sv_out_build_unique(&raw, &uniq);
    int32_t n = uniq.n;
    if (vb_flag) fprintf(stderr, "%d input lines, %d distinct reads\n", raw.n_lines, n);

    const char **seqs_c = (const char **)malloc((size_t)n * sizeof(char *));
    for (int32_t i = 0; i < n; i++) seqs_c[i] = uniq.seqs[i];

    /* ---------------------------------------------------------- pipeline */
    if (legacy_flag) anchor_occ = 1;
    int32_t *lab1 = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    sv_work_t w1; int32_t iters1;
    sv_s1_cores((const char *const *)seqs_c, uniq.lens, n, tau_core,
                SV_N_ITER_S1, SV_THETA_HIGH, anchor_occ, lab1, &w1, &iters1);

    sv_core_list_t cores;
    int32_t *orph_idx; int32_t norph;
    sv_split_cores_orphans(lab1, n, &cores, &orph_idx, &norph);

    if (vb_flag) {
        cal_coverage_t cv = cal_core_coverage(lab1, n);
        fprintf(stderr, "S1: %d cores, %d orphans, %.2f%% of distinct reads in cores (anchor occurrences: %d)\n",
                cv.core_count, cv.orphan_count, 100.0 * cv.frac_reads_in_cores, anchor_occ);
    }
    if (audit_sample > 0) {
        ec_audit_t au;
        ec_s1_audit((const char *const *)seqs_c, uniq.lens, n, lab1, tau_core, audit_sample, &au);
        if (vb_flag)
            fprintf(stderr, "S1 audit: %d sampled reads, %lld exact neighbours within %d, "
                    "%lld of them in another S1 component (%d sampled reads affected)\n",
                    au.sampled, (long long)au.neighbours, tau_core, (long long)au.cross, au.reads_with_cross);
    }

    sv_s2_result_t s2;
    sv_s2_consensus(&cores, (const char *const *)seqs_c, &s2);
    int32_t K = s2.count;

    const char **orph_seqs = (const char **)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(char *));
    int32_t *orph_lens = (int32_t *)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(int32_t));
    for (int32_t k = 0; k < norph; k++) {
        orph_seqs[k] = seqs_c[orph_idx[k]];
        orph_lens[k] = uniq.lens[orph_idx[k]];
    }
    int32_t *cons_lens = (int32_t *)malloc((size_t)(K > 0 ? K : 1) * sizeof(int32_t));
    for (int32_t k = 0; k < K; k++) cons_lens[k] = (int32_t)strlen(s2.cons[k]);

    int exact = !legacy_flag && K >= 2 && K <= exact_max;
    if (!legacy_flag && K > exact_max && vb_flag)
        fprintf(stderr, "warning: %d cores > --exact-max %d, falling back to filtered S3/S4 "
                "(mode=approx, no purity guarantee)\n", K, exact_max);

    int32_t *bd = (int32_t *)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(int32_t));
    int32_t *bi = (int32_t *)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(int32_t));
    uint8_t *absorbed = (uint8_t *)calloc((size_t)(norph > 0 ? norph : 1), 1);
    int32_t *lab4 = (int32_t *)malloc((size_t)(K > 0 ? K : 1) * sizeof(int32_t));
    for (int32_t k = 0; k < norph; k++) bi[k] = -1;

    if (exact) {
        /* ---- H: exact distance histogram over all consensus pairs */
        ec_peq_t *peq = (ec_peq_t *)malloc((size_t)K * sizeof(ec_peq_t));
        for (int32_t k = 0; k < K; k++) ec_peq_build(s2.cons[k], cons_lens[k], &peq[k]);
        int32_t dmax = cal_estimate_dmax(cons_lens, K, 0.5);
        int64_t *hist = (int64_t *)calloc((size_t)dmax + 2, sizeof(int64_t));
        ec_edges_t edges = {0};

        /* Edges are only needed up to 2*tau (S3 neighbour search) -- storing
        ** everything up to dmax can run into millions of pairs on long reads.
        ** With --auto-tau tau isn't known yet, so that case takes two passes. */
        int32_t keep;
        if (auto_flag) {
            ec_all_pairs(peq, (const char *const *)s2.cons, cons_lens, K, dmax, -1, hist, NULL);
        } else {
            int32_t tmax = tau_abs > tau_cons ? tau_abs : tau_cons;
            keep = 2 * tmax + 1 < dmax ? 2 * tmax + 1 : dmax;
            ec_all_pairs(peq, (const char *const *)s2.cons, cons_lens, K, dmax, keep, hist, &edges);
        }

        int32_t peak = -1, gap_start = 0, gap_len = 0, ceiling = -1; int kind = 0;
        if (!cal_ceiling_below_peak(hist, dmax, CAL_MIN_GAP, &peak, &gap_start, &gap_len, &ceiling, &kind))
            ceiling = dmax + 1;   /* no pair within dmax at all */
        int32_t tau_safe = (int32_t)((double)(ceiling - 1) * (1.0 - margin));
        if (tau_safe > ceiling - 1) tau_safe = ceiling - 1;
        if (tau_safe < 2) tau_safe = 2;

        if (vb_flag) {
            int64_t np = 0; for (int32_t d = 0; d <= dmax; d++) np += hist[d];
            if (kind == 1)
                fprintf(stderr, "H: %d consensus, %lld pairs within %d, peak d=%d, gap d=%d..%d, ceiling C=%d\n",
                        K, (long long)np, dmax, peak, gap_start, gap_start + gap_len - 1, ceiling);
            else
                fprintf(stderr, "H: %d consensus, %lld pairs within %d, no gap below the peak, ceiling C=%d\n",
                        K, (long long)np, dmax, ceiling);
        }
        if (hist_path) {
            FILE *fh = fopen(hist_path, "w");
            if (fh) { for (int32_t d = 0; d <= dmax; d++) fprintf(fh, "%d\t%lld\n", d, (long long)hist[d]); fclose(fh); }
        }

        if (auto_flag) {
            tau_abs = tau_cons = tau_safe;
            int32_t want = 2 * tau_safe + 1;
            keep = want < dmax ? want : dmax;
            ec_all_pairs(peq, (const char *const *)s2.cons, cons_lens, K, dmax, keep, hist, &edges);
            if (vb_flag) fprintf(stderr, "auto-tau: tau_abs=tau_cons=%d (margin %.2f)\n", tau_safe, margin);
        } else {
            /* a fixed threshold at or above the ceiling would join different oligos */
            if (tau_cons >= ceiling) {
                if (vb_flag) fprintf(stderr, "warning: tau_cons=%d >= ceiling %d, lowered to %d\n", tau_cons, ceiling, tau_safe);
                tau_cons = tau_safe;
            }
            if (tau_abs >= ceiling) {
                if (vb_flag) fprintf(stderr, "warning: tau_abs=%d >= ceiling %d, lowered to %d\n", tau_abs, ceiling, tau_safe);
                tau_abs = tau_safe;
            }
        }

        /* ---- S4 before S3: S3 needs the components and the gap between them */
        int32_t c_post;
        int64_t unions = ec_merge(&edges, K, tau_cons, keep, lab4, &c_post);
        if (vb_flag)
            fprintf(stderr, "S4: %d consensus -> %lld components at tau_cons=%d, closest pair across components d%s%d\n",
                    K, (long long)(K - unions), tau_cons, c_post > keep ? ">=" : "=", c_post);

        /* ---- S3: filtered candidates first, then the exact check */
        if (norph > 0) {
            int32_t *bd0 = (int32_t *)malloc((size_t)norph * sizeof(int32_t));
            int32_t *bi0 = (int32_t *)malloc((size_t)norph * sizeof(int32_t));
            uint8_t *ab0 = (uint8_t *)calloc((size_t)norph, 1);
            sv_work_t w3 = {0,0,0,0};
            sv_s3_absorb(orph_seqs, orph_lens, norph, (const char *const *)s2.cons, cons_lens, K,
                         tau_abs, sig_abs, SV_N_ITER_S3, anchor_occ, bd0, bi0, ab0, &w3);
            ec_s3_stats_t st;
            ec_s3_resolve(orph_seqs, orph_lens, norph, peq, (const char *const *)s2.cons, cons_lens, K,
                          lab4, &edges, keep, c_post, tau_abs, bd0, bi0, bi, absorbed, &st);
            for (int32_t k = 0; k < norph; k++) bd[k] = bd0[k];
            if (vb_flag) {
                int64_t na = 0; for (int32_t k = 0; k < norph; k++) na += absorbed[k];
                fprintf(stderr, "S3: %lld of %d orphans absorbed at tau_abs=%d; %lld proven by 2d < %d, "
                        "%lld checked against neighbours, %lld full scans; %lld moved to another component, "
                        "%lld found only by the exact check, %lld left alone on a tie\n",
                        (long long)na, norph, tau_abs, (long long)st.fast, c_post, (long long)st.local,
                        (long long)st.full, (long long)st.moved, (long long)st.added, (long long)st.tie_rejected);
            }
            free(bd0); free(bi0); free(ab0);
        }
        for (int32_t k = 0; k < K; k++) ec_peq_free(&peq[k]);
        free(peq); free(hist); ec_edges_free(&edges);
    } else {
        /* ---- v1.1 path (also used past --exact-max) */
        if (auto_flag && K > 0) {
            int32_t tau_out = -1, ceiling_out = -1;
            cal_diagnostics_t dg;
            cal_auto_tau((const char *const *)s2.cons, cons_lens, K, margin, SV_N_ITER_S3, sig_abs,
                         &tau_out, &ceiling_out, &dg);
            if (hist_path && dg.hist) {
                FILE *fh = fopen(hist_path, "w");
                if (fh) { for (int32_t d = 0; d <= dg.dmax; d++) fprintf(fh, "%d\t%lld\n", d, (long long)dg.hist[d]); fclose(fh); }
            }
            if (tau_out > 0) {
                tau_abs = tau_out;
                tau_cons = tau_out;
                if (vb_flag) {
                    if (dg.kind == 1)
                        fprintf(stderr, "auto-tau: %lld candidate pairs, peak d=%d, gap d=%d..%d, ceiling=%d -> tau_abs=tau_cons=%d (margin %.2f)\n",
                                (long long)dg.n_pairs, dg.peak, dg.gap_start, dg.gap_start + dg.gap_len - 1, ceiling_out, tau_out, margin);
                    else
                        fprintf(stderr, "auto-tau: %lld candidate pairs, peak d=%d, no gap below the peak (no same-oligo mode); "
                                "ceiling = smallest distance %d -> tau_abs=tau_cons=%d (margin %.2f)\n",
                                (long long)dg.n_pairs, dg.peak, ceiling_out, tau_out, margin);
                }
            } else if (vb_flag) {
                fprintf(stderr, "auto-tau: empty histogram, falling back to default tau_abs=%d tau_cons=%d\n",
                        tau_abs, tau_cons);
            }
            free(dg.hist);
        }
        sv_work_t w3 = {0,0,0,0};
        if (norph > 0 && K > 0)
            sv_s3_absorb(orph_seqs, orph_lens, norph, (const char *const *)s2.cons, cons_lens, K,
                         tau_abs, sig_abs, SV_N_ITER_S3, anchor_occ, bd, bi, absorbed, &w3);
        sv_work_t w4 = {0,0,0,0};
        if (K > 0)
            sv_s4_merge_cons((const char *const *)s2.cons, cons_lens, K, tau_cons, sig_abs, SV_N_ITER_S4, lab4, &w4);
    }

    int32_t *final_lab = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    sv_assemble(&cores, orph_idx, norph, bi, absorbed, lab4, n, final_lab);

    int32_t *canon = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    sv_canonical_partition(final_lab, seqs_c, n, canon);

    if (vb_flag) {
        int32_t maxc = -1;
        for (int32_t i = 0; i < n; i++) if (canon[i] > maxc) maxc = canon[i];
        fprintf(stderr, "done: %d clusters from %d distinct reads\n", maxc + 1, n);
    }

    sv_out_opts_t opts = { .show_clusters = cl_flag, .show_seqid = id_flag, .tidy = tidy_flag };
    sv_out_write(fout, &uniq, canon, &raw, &opts);

    if (fout != stdout) fclose(fout);

    free(seqs_c); free(lab1); free(orph_idx);
    free(orph_seqs); free(orph_lens); free(cons_lens);
    free(bd); free(bi); free(absorbed); free(lab4); free(final_lab); free(canon);
    sv_core_list_free(&cores);
    sv_s2_result_free(&s2);
    sv_out_free_unique(&uniq);
    sv_out_free_raw(&raw);

    return EXIT_SUCCESS;
}
