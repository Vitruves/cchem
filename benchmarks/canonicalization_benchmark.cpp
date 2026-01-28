/**
 * C/C++ Benchmark for SMILES canonicalization
 *
 * Compares native performance of:
 * - cchem (this library)
 * - RDKit (C++ API)
 * - OpenBabel (C++ API)
 *
 * Build:
 *   mkdir build && cd build
 *   cmake .. -DBUILD_BENCHMARKS=ON
 *   make canonicalization_benchmark
 *
 * Usage:
 *   ./canonicalization_benchmark -f molecules.csv -s smiles -n 8
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <getopt.h>
#include <sys/time.h>

// cchem
extern "C" {
#include "cchem/canonicalizer/parser.h"
#include "cchem/canonicalizer/canon.h"
}

// Conditional includes for RDKit and OpenBabel
#ifdef HAVE_RDKIT
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#endif

#ifdef HAVE_OPENBABEL
// OpenBabel uses std::binary_function which was removed in C++17
// Provide a compatibility shim for C++17/20
#if __cplusplus >= 201703L
namespace std {
    template<typename Arg1, typename Arg2, typename Result>
    struct binary_function {
        typedef Arg1 first_argument_type;
        typedef Arg2 second_argument_type;
        typedef Result result_type;
    };
}
#endif
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/op.h>
#endif

// ============================================================================
// Utilities
// ============================================================================

static double get_time_ms(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000.0 + tv.tv_usec / 1000.0;
}

typedef struct {
    char **smiles;
    size_t count;
    size_t capacity;
} smiles_list_t;

static smiles_list_t *smiles_list_create(void) {
    smiles_list_t *list = (smiles_list_t *)malloc(sizeof(smiles_list_t));
    list->count = 0;
    list->capacity = 1024;
    list->smiles = (char **)malloc(list->capacity * sizeof(char *));
    return list;
}

static void smiles_list_add(smiles_list_t *list, const char *smi) {
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        list->smiles = (char **)realloc(list->smiles, list->capacity * sizeof(char *));
    }
    list->smiles[list->count++] = strdup(smi);
}

static void smiles_list_free(smiles_list_t *list) {
    for (size_t i = 0; i < list->count; i++) {
        free(list->smiles[i]);
    }
    free(list->smiles);
    free(list);
}

static smiles_list_t *load_smiles_csv(const char *filepath, const char *smiles_col, size_t limit) {
    FILE *f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open file %s\n", filepath);
        return NULL;
    }

    smiles_list_t *list = smiles_list_create();
    char line[65536];
    int smiles_idx = -1;
    int line_num = 0;

    while (fgets(line, sizeof(line), f)) {
        // Remove newline
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[--len] = '\0';
        }

        if (line_num == 0) {
            // Parse header to find smiles column
            char *token = strtok(line, ",");
            int idx = 0;
            while (token) {
                if (strcmp(token, smiles_col) == 0) {
                    smiles_idx = idx;
                    break;
                }
                token = strtok(NULL, ",");
                idx++;
            }
            if (smiles_idx < 0) {
                fprintf(stderr, "Error: Column '%s' not found in CSV\n", smiles_col);
                fclose(f);
                smiles_list_free(list);
                return NULL;
            }
        } else {
            // Parse data row
            char *token = strtok(line, ",");
            int idx = 0;
            while (token) {
                if (idx == smiles_idx) {
                    if (strlen(token) > 0) {
                        smiles_list_add(list, token);
                    }
                    break;
                }
                token = strtok(NULL, ",");
                idx++;
            }
        }

        line_num++;
        if (limit > 0 && list->count >= limit) break;
    }

    fclose(f);
    return list;
}

// ============================================================================
// Benchmark result
// ============================================================================

typedef struct {
    const char *library;
    double total_time_ms;
    size_t molecules_processed;
    size_t molecules_failed;
    double throughput;      // mol/s
    double avg_time_us;     // microseconds per molecule
    int threads;
} benchmark_result_t;

static void print_result(const benchmark_result_t *r) {
    printf("  %-20s %d threads: %10.0f mol/s  (%.2f us/mol, %zu failed)\n",
           r->library, r->threads, r->throughput, r->avg_time_us, r->molecules_failed);
}

// ============================================================================
// cchem benchmark
// ============================================================================

typedef struct {
    char **smiles;
    size_t start;
    size_t end;
    size_t success;
    size_t failed;
} cchem_thread_args_t;

static void *cchem_worker(void *arg) {
    cchem_thread_args_t *args = (cchem_thread_args_t *)arg;
    args->success = 0;
    args->failed = 0;

    char error_buf[256];
    for (size_t i = args->start; i < args->end; i++) {
        molecule_t *mol = smiles_to_molecule(args->smiles[i], error_buf, sizeof(error_buf));
        if (mol) {
            char *canonical = molecule_to_canonical_smiles(mol, NULL);
            if (canonical && strlen(canonical) > 0) {
                args->success++;
                free(canonical);
            } else {
                args->failed++;
            }
            molecule_free(mol);
        } else {
            args->failed++;
        }
    }

    return NULL;
}

static benchmark_result_t benchmark_cchem(char **smiles, size_t count, int threads) {
    benchmark_result_t result = {0};
    result.library = "cchem";
    result.threads = threads;

    double start = get_time_ms();

    if (threads == 1) {
        // Single-threaded
        size_t success = 0, failed = 0;
        char error_buf[256];
        for (size_t i = 0; i < count; i++) {
            molecule_t *mol = smiles_to_molecule(smiles[i], error_buf, sizeof(error_buf));
            if (mol) {
                char *canonical = molecule_to_canonical_smiles(mol, NULL);
                if (canonical && strlen(canonical) > 0) {
                    success++;
                    free(canonical);
                } else {
                    failed++;
                }
                molecule_free(mol);
            } else {
                failed++;
            }
        }
        result.molecules_processed = success;
        result.molecules_failed = failed;
    } else {
        // Multi-threaded
        pthread_t *thread_ids = (pthread_t *)malloc(threads * sizeof(pthread_t));
        cchem_thread_args_t *thread_args = (cchem_thread_args_t *)malloc(threads * sizeof(cchem_thread_args_t));

        size_t per_thread = count / threads;
        size_t remainder = count % threads;

        size_t offset = 0;
        for (int t = 0; t < threads; t++) {
            thread_args[t].smiles = smiles;
            thread_args[t].start = offset;
            thread_args[t].end = offset + per_thread + (t < (int)remainder ? 1 : 0);
            offset = thread_args[t].end;
            pthread_create(&thread_ids[t], NULL, cchem_worker, &thread_args[t]);
        }

        size_t total_success = 0, total_failed = 0;
        for (int t = 0; t < threads; t++) {
            pthread_join(thread_ids[t], NULL);
            total_success += thread_args[t].success;
            total_failed += thread_args[t].failed;
        }

        result.molecules_processed = total_success;
        result.molecules_failed = total_failed;

        free(thread_ids);
        free(thread_args);
    }

    double elapsed = get_time_ms() - start;
    result.total_time_ms = elapsed;
    result.throughput = (result.molecules_processed / elapsed) * 1000.0;
    result.avg_time_us = (elapsed / result.molecules_processed) * 1000.0;

    return result;
}

// ============================================================================
// RDKit benchmark
// ============================================================================

#ifdef HAVE_RDKIT

typedef struct {
    char **smiles;
    size_t start;
    size_t end;
    size_t success;
    size_t failed;
} rdkit_thread_args_t;

static void *rdkit_worker(void *arg) {
    rdkit_thread_args_t *args = (rdkit_thread_args_t *)arg;
    args->success = 0;
    args->failed = 0;

    for (size_t i = args->start; i < args->end; i++) {
        try {
            RDKit::ROMol *mol = RDKit::SmilesToMol(args->smiles[i]);
            if (mol) {
                std::string canonical = RDKit::MolToSmiles(*mol);
                if (!canonical.empty()) {
                    args->success++;
                } else {
                    args->failed++;
                }
                delete mol;
            } else {
                args->failed++;
            }
        } catch (...) {
            args->failed++;
        }
    }

    return NULL;
}

static benchmark_result_t benchmark_rdkit(char **smiles, size_t count, int threads) {
    benchmark_result_t result = {0};
    result.library = "RDKit";
    result.threads = threads;

    double start = get_time_ms();

    if (threads == 1) {
        size_t success = 0, failed = 0;
        for (size_t i = 0; i < count; i++) {
            try {
                RDKit::ROMol *mol = RDKit::SmilesToMol(smiles[i]);
                if (mol) {
                    std::string canonical = RDKit::MolToSmiles(*mol);
                    if (!canonical.empty()) {
                        success++;
                    } else {
                        failed++;
                    }
                    delete mol;
                } else {
                    failed++;
                }
            } catch (...) {
                failed++;
            }
        }
        result.molecules_processed = success;
        result.molecules_failed = failed;
    } else {
        pthread_t *thread_ids = (pthread_t *)malloc(threads * sizeof(pthread_t));
        rdkit_thread_args_t *thread_args = (rdkit_thread_args_t *)malloc(threads * sizeof(rdkit_thread_args_t));

        size_t per_thread = count / threads;
        size_t remainder = count % threads;

        size_t offset = 0;
        for (int t = 0; t < threads; t++) {
            thread_args[t].smiles = smiles;
            thread_args[t].start = offset;
            thread_args[t].end = offset + per_thread + (t < (int)remainder ? 1 : 0);
            offset = thread_args[t].end;
            pthread_create(&thread_ids[t], NULL, rdkit_worker, &thread_args[t]);
        }

        size_t total_success = 0, total_failed = 0;
        for (int t = 0; t < threads; t++) {
            pthread_join(thread_ids[t], NULL);
            total_success += thread_args[t].success;
            total_failed += thread_args[t].failed;
        }

        result.molecules_processed = total_success;
        result.molecules_failed = total_failed;

        free(thread_ids);
        free(thread_args);
    }

    double elapsed = get_time_ms() - start;
    result.total_time_ms = elapsed;
    result.throughput = (result.molecules_processed / elapsed) * 1000.0;
    result.avg_time_us = (elapsed / result.molecules_processed) * 1000.0;

    return result;
}

#endif // HAVE_RDKIT

// ============================================================================
// OpenBabel benchmark
// ============================================================================

#ifdef HAVE_OPENBABEL

typedef struct {
    char **smiles;
    size_t start;
    size_t end;
    size_t success;
    size_t failed;
} openbabel_thread_args_t;

static void *openbabel_worker(void *arg) {
    openbabel_thread_args_t *args = (openbabel_thread_args_t *)arg;
    args->success = 0;
    args->failed = 0;

    OpenBabel::OBConversion conv;
    conv.SetInFormat("smi");
    conv.SetOutFormat("can");

    for (size_t i = args->start; i < args->end; i++) {
        try {
            OpenBabel::OBMol mol;
            if (conv.ReadString(&mol, args->smiles[i])) {
                std::string canonical = conv.WriteString(&mol, true);
                if (!canonical.empty()) {
                    args->success++;
                } else {
                    args->failed++;
                }
            } else {
                args->failed++;
            }
        } catch (...) {
            args->failed++;
        }
    }

    return NULL;
}

static benchmark_result_t benchmark_openbabel(char **smiles, size_t count, int threads) {
    benchmark_result_t result = {0};
    result.library = "OpenBabel";

    // OpenBabel is not thread-safe, always run single-threaded
    if (threads > 1) {
        printf("  Warning: OpenBabel runs single-threaded (not thread-safe)\n");
    }
    result.threads = 1;

    double start = get_time_ms();

    {
        size_t success = 0, failed = 0;
        OpenBabel::OBConversion conv;
        conv.SetInFormat("smi");
        conv.SetOutFormat("can");

        for (size_t i = 0; i < count; i++) {
            try {
                OpenBabel::OBMol mol;
                if (conv.ReadString(&mol, smiles[i])) {
                    std::string canonical = conv.WriteString(&mol, true);
                    if (!canonical.empty()) {
                        success++;
                    } else {
                        failed++;
                    }
                } else {
                    failed++;
                }
            } catch (...) {
                failed++;
            }
        }
        result.molecules_processed = success;
        result.molecules_failed = failed;
    }

    double elapsed = get_time_ms() - start;
    result.total_time_ms = elapsed;
    result.throughput = (result.molecules_processed / elapsed) * 1000.0;
    result.avg_time_us = (elapsed / result.molecules_processed) * 1000.0;

    return result;
}

#endif // HAVE_OPENBABEL

// ============================================================================
// Main
// ============================================================================

static void print_usage(const char *prog) {
    printf("Usage: %s [options]\n", prog);
    printf("\nOptions:\n");
    printf("  -f, --file <path>       Input CSV file with SMILES (required)\n");
    printf("  -s, --smiles-col <name> Column name containing SMILES (default: smiles)\n");
    printf("  -n, --threads <num>     Number of threads (default: 1)\n");
    printf("  -i, --iterations <num>  Number of iterations (default: 3)\n");
    printf("  -l, --limit <num>       Limit number of molecules\n");
    printf("  -h, --help              Show this help\n");
    printf("\nAvailable libraries:\n");
    printf("  - cchem: Yes\n");
#ifdef HAVE_RDKIT
    printf("  - RDKit: Yes\n");
#else
    printf("  - RDKit: No (rebuild with -DHAVE_RDKIT=ON)\n");
#endif
#ifdef HAVE_OPENBABEL
    printf("  - OpenBabel: Yes\n");
#else
    printf("  - OpenBabel: No (rebuild with -DHAVE_OPENBABEL=ON)\n");
#endif
}

int main(int argc, char **argv) {
    const char *filepath = NULL;
    const char *smiles_col = "smiles";
    int threads = 1;
    int iterations = 3;
    size_t limit = 0;

    static struct option long_options[] = {
        {"file",       required_argument, 0, 'f'},
        {"smiles-col", required_argument, 0, 's'},
        {"threads",    required_argument, 0, 'n'},
        {"iterations", required_argument, 0, 'i'},
        {"limit",      required_argument, 0, 'l'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "f:s:n:i:l:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'f': filepath = optarg; break;
            case 's': smiles_col = optarg; break;
            case 'n': threads = atoi(optarg); break;
            case 'i': iterations = atoi(optarg); break;
            case 'l': limit = (size_t)atol(optarg); break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (!filepath) {
        fprintf(stderr, "Error: Input file required (-f)\n");
        print_usage(argv[0]);
        return 1;
    }

    // Load SMILES
    printf("Loading SMILES from %s...\n", filepath);
    smiles_list_t *smiles = load_smiles_csv(filepath, smiles_col, limit);
    if (!smiles) {
        return 1;
    }
    printf("Loaded %zu SMILES\n\n", smiles->count);

    // Print header
    printf("=============================================================================\n");
    printf("C/C++ CANONICALIZATION BENCHMARK\n");
    printf("=============================================================================\n");
    printf("Libraries: cchem");
#ifdef HAVE_RDKIT
    printf(", RDKit");
#endif
#ifdef HAVE_OPENBABEL
    printf(", OpenBabel");
#endif
    printf("\n");
    printf("Molecules: %zu\n", smiles->count);
    printf("Iterations: %d\n", iterations);
    printf("Threads: 1 and %d\n", threads);
    printf("=============================================================================\n\n");

    // Aggregate results
    double cchem_1t_throughput = 0, cchem_mt_throughput = 0;
#ifdef HAVE_RDKIT
    double rdkit_1t_throughput = 0, rdkit_mt_throughput = 0;
#endif
#ifdef HAVE_OPENBABEL
    double openbabel_1t_throughput = 0, openbabel_mt_throughput = 0;
#endif

    for (int iter = 0; iter < iterations; iter++) {
        printf("--- Iteration %d/%d ---\n", iter + 1, iterations);

        // Single-threaded
        printf("\nSingle-threaded:\n");

        benchmark_result_t r = benchmark_cchem(smiles->smiles, smiles->count, 1);
        print_result(&r);
        cchem_1t_throughput += r.throughput;

#ifdef HAVE_RDKIT
        r = benchmark_rdkit(smiles->smiles, smiles->count, 1);
        print_result(&r);
        rdkit_1t_throughput += r.throughput;
#endif

#ifdef HAVE_OPENBABEL
        r = benchmark_openbabel(smiles->smiles, smiles->count, 1);
        print_result(&r);
        openbabel_1t_throughput += r.throughput;
#endif

        // Multi-threaded
        if (threads > 1) {
            printf("\nMulti-threaded (n=%d):\n", threads);

            r = benchmark_cchem(smiles->smiles, smiles->count, threads);
            print_result(&r);
            cchem_mt_throughput += r.throughput;

#ifdef HAVE_RDKIT
            r = benchmark_rdkit(smiles->smiles, smiles->count, threads);
            print_result(&r);
            rdkit_mt_throughput += r.throughput;
#endif

#ifdef HAVE_OPENBABEL
            r = benchmark_openbabel(smiles->smiles, smiles->count, threads);
            print_result(&r);
            openbabel_mt_throughput += r.throughput;
#endif
        }

        printf("\n");
    }

    // Print summary
    printf("=============================================================================\n");
    printf("SUMMARY (averaged over %d iterations)\n", iterations);
    printf("=============================================================================\n\n");

    printf("%-20s %15s %15s %10s\n", "Library", "1-thread", threads > 1 ? "N-thread" : "", "Speedup");
    printf("-----------------------------------------------------------------------------\n");

    cchem_1t_throughput /= iterations;
    cchem_mt_throughput /= iterations;
    printf("%-20s %12.0f/s", "cchem", cchem_1t_throughput);
    if (threads > 1) {
        printf(" %12.0f/s %9.2fx", cchem_mt_throughput, cchem_mt_throughput / cchem_1t_throughput);
    }
    printf("\n");

#ifdef HAVE_RDKIT
    rdkit_1t_throughput /= iterations;
    rdkit_mt_throughput /= iterations;
    printf("%-20s %12.0f/s", "RDKit", rdkit_1t_throughput);
    if (threads > 1) {
        printf(" %12.0f/s %9.2fx", rdkit_mt_throughput, rdkit_mt_throughput / rdkit_1t_throughput);
    }
    printf("\n");
#endif

#ifdef HAVE_OPENBABEL
    openbabel_1t_throughput /= iterations;
    openbabel_mt_throughput /= iterations;
    printf("%-20s %12.0f/s", "OpenBabel", openbabel_1t_throughput);
    if (threads > 1) {
        printf(" %12.0f/s %9.2fx", openbabel_mt_throughput, openbabel_mt_throughput / openbabel_1t_throughput);
    }
    printf("\n");
#endif

    printf("-----------------------------------------------------------------------------\n\n");

    // Print comparison vs RDKit
#ifdef HAVE_RDKIT
    printf("cchem vs RDKit:\n");
    printf("  Single-threaded: %.2fx faster\n", cchem_1t_throughput / rdkit_1t_throughput);
    if (threads > 1) {
        printf("  Multi-threaded:  %.2fx faster\n", cchem_mt_throughput / rdkit_mt_throughput);
    }
    printf("\n");
#endif

    smiles_list_free(smiles);
    return 0;
}
