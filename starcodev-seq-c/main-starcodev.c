/*
** main-starcodev.c -- StarcodeV's CLI driver, using the same getopt_long
** style as src/main-starcode.c from original Starcode (gui11aume/starcode
** v1.4). See CONVENTIONS_STARCODEV.md §1: the legacy options (-i -o -t -d
** -c -q -v -h, --print-clusters, --seq-id, --tidy) keep their original
** meaning and names; the new options (--tau-core, --tau-abs, --tau-cons,
** --auto-tau, --sig-abs, --margin) are added on top without clobbering
** them.
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
"       --margin: fraction of margin below the non-sibling ceiling for --auto-tau (default 0.75)\n"
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

#define VERSION "starcodev-v1.0 (sequential C port of the validated starcodev.py reference)"

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
    double margin = 0.75;
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
    int32_t *lab1 = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    sv_work_t w1; int32_t iters1;
    sv_s1_cores((const char *const *)seqs_c, uniq.lens, n, tau_core,
                SV_N_ITER_S1, SV_THETA_HIGH, lab1, &w1, &iters1);

    sv_core_list_t cores;
    int32_t *orph_idx; int32_t norph;
    sv_split_cores_orphans(lab1, n, &cores, &orph_idx, &norph);

    sv_s2_result_t s2;
    sv_s2_consensus(&cores, (const char *const *)seqs_c, &s2);

    const char **orph_seqs = (const char **)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(char *));
    int32_t *orph_lens = (int32_t *)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(int32_t));
    for (int32_t k = 0; k < norph; k++) {
        orph_seqs[k] = seqs_c[orph_idx[k]];
        orph_lens[k] = uniq.lens[orph_idx[k]];
    }

    int32_t *cons_lens = (int32_t *)malloc((size_t)(s2.count > 0 ? s2.count : 1) * sizeof(int32_t));
    for (int32_t k = 0; k < s2.count; k++) cons_lens[k] = (int32_t)strlen(s2.cons[k]);

    if (auto_flag && s2.count > 0) {
        int32_t tau_out = -1, ceiling_out = -1;
        cal_auto_tau((const char *const *)s2.cons, cons_lens, s2.count, margin, SV_N_ITER_S3,
                     &tau_out, &ceiling_out, NULL);
        if (tau_out > 0) {
            tau_abs = tau_out;
            tau_cons = tau_out;
            if (vb_flag) fprintf(stderr, "auto-tau: tau_abs=tau_cons=%d (ceiling=%d)\n", tau_out, ceiling_out);
        } else if (vb_flag) {
            fprintf(stderr, "auto-tau: no gap found, falling back to default tau_abs=%d tau_cons=%d\n",
                    tau_abs, tau_cons);
        }
    }

    int32_t *bd = (int32_t *)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(int32_t));
    int32_t *bi = (int32_t *)malloc((size_t)(norph > 0 ? norph : 1) * sizeof(int32_t));
    uint8_t *absorbed = (uint8_t *)calloc((size_t)(norph > 0 ? norph : 1), 1);
    sv_work_t w3 = {0,0,0,0};
    if (norph > 0 && s2.count > 0) {
        sv_s3_absorb(orph_seqs, orph_lens, norph, (const char *const *)s2.cons, cons_lens, s2.count,
                     tau_abs, sig_abs, SV_N_ITER_S3, bd, bi, absorbed, &w3);
    } else {
        for (int32_t k = 0; k < norph; k++) bi[k] = -1;
    }

    int32_t *lab4 = (int32_t *)malloc((size_t)(s2.count > 0 ? s2.count : 1) * sizeof(int32_t));
    sv_work_t w4 = {0,0,0,0};
    if (s2.count > 0) {
        sv_s4_merge_cons((const char *const *)s2.cons, cons_lens, s2.count, tau_cons, sig_abs, SV_N_ITER_S4, lab4, &w4);
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
