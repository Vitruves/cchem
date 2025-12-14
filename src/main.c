/**
 * @file main.c
 * @brief cchem CLI - SMILES canonicalizer command line interface
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <limits.h>
#include "cchem/cchem.h"
#include "cchem/descriptors.h"
#include "cchem/csv.h"
#include "cchem/progress.h"
#include "cchem/threading.h"
#include "cchem/memory.h"
#include "cchem/depictor/depictor.h"
#include "cchem/splitter.h"
#include <stdatomic.h>

#define VERSION "1.0.0"

static void print_version(void) {
    printf("cchem version %s\n", VERSION);
    printf("SMILES canonicalization tool\n");
}

static void print_usage(const char* prog_name) {
    printf("Usage: %s <command> [options]\n\n", prog_name);
    printf("Commands:\n");
    printf("  canonicalize   Canonicalize SMILES strings\n");
    printf("  compute        Compute molecular descriptors\n");
    printf("  depict         Generate 2D/3D molecular structure images\n");
    printf("  split          Split dataset into train/test/validation sets\n");
    printf("  validate       Validate SMILES syntax\n");
    printf("  version        Print version information\n");
    printf("  help           Print this help message\n");
    printf("\n");
    printf("For command-specific help:\n");
    printf("  %s <command> --help\n", prog_name);
}

static void print_canonicalize_usage(const char* prog_name) {
    printf("Usage: %s canonicalize [options]\n\n", prog_name);
    printf("Canonicalize SMILES strings (aromatized by default)\n\n");
    printf("Options:\n");
    printf("  -S, --smiles <SMILES>   Single SMILES string to canonicalize\n");
    printf("  -f, --file <file.csv>   Input CSV file\n");
    printf("  -s, --col <name>        Input SMILES column name (required with -f)\n");
    printf("  -c, --outcol <name>     Output column name for canonical SMILES (default: canonical_smiles)\n");
    printf("  -o, --output <file.csv> Output CSV file (required with -f)\n");
    printf("  -n, --ncpu <N>          Number of CPU cores for parallel processing (default: auto)\n");
    printf("  --sanitize <opts>       Apply sanitization before canonicalization\n");
    printf("                          Values: \"complete\" or comma-separated list of:\n");
    printf("                            unsalt          - Remove salts, keep largest fragment\n");
    printf("                            aromatize       - Perceive and apply aromaticity (default)\n");
    printf("                            kekulize        - Convert aromatic to Kekule form\n");
    printf("                            neutralize      - Neutralize charges\n");
    printf("                            normalize       - Normalize functional groups (nitro, etc.)\n");
    printf("                            remove-stereo   - Remove stereochemistry\n");
    printf("                            remove-isotopes - Remove isotope labels\n");
    printf("                            remove-h        - Remove explicit hydrogens\n");
    printf("                            validate        - Validate structure\n");
    printf("  --list-tautomers        List all tautomeric forms instead of canonicalizing\n");
    printf("  -v, --verbose           Verbose output\n");
    printf("  -h, --help              Print this help message\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s canonicalize -S \"c1ccccc1\"\n", prog_name);
    printf("  %s canonicalize -S \"[Na+].CC(=O)[O-]\" --sanitize unsalt\n", prog_name);
    printf("  %s canonicalize -S \"C1=CC=CC=C1\" --sanitize kekulize\n", prog_name);
    printf("  %s canonicalize -S \"CC(=O)NC\" --list-tautomers\n", prog_name);
    printf("  %s canonicalize -f input.csv -s smiles -c canonical -o output.csv -n 4\n", prog_name);
}

static void print_validate_usage(const char* prog_name) {
    printf("Usage: %s validate [options]\n\n", prog_name);
    printf("Validate SMILES syntax\n\n");
    printf("Options:\n");
    printf("  -S, --smiles <SMILES>   SMILES string to validate\n");
    printf("  -v, --verbose           Verbose output\n");
    printf("  -h, --help              Print this help message\n");
}

static void print_compute_usage(const char* prog_name) {
    printf("Usage: %s compute [options]\n\n", prog_name);
    printf("Compute molecular descriptors from SMILES\n\n");
    printf("Options:\n");
    printf("  -S, --smiles <SMILES>       Single SMILES string\n");
    printf("  -f, --file <file.csv>       Input CSV file\n");
    printf("  -s, --col <name>            Input SMILES column name (required with -f)\n");
    printf("  -d, --descriptors <list>    Comma-separated descriptor names or \"all\" (default: all)\n");
    printf("  -o, --output <file.csv>     Output CSV file (required with -f)\n");
    printf("  -n, --ncpu <N>              Number of CPU cores (default: auto)\n");
    printf("  --no-canonicalization       Skip canonicalization, parse input SMILES directly (faster)\n");
    printf("  -l, --list                  List available descriptors\n");
    printf("  -v, --verbose               Verbose output\n");
    printf("  -h, --help                  Print this help message\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s compute -S \"CCO\" -d CarbonCount,HydrogenCount\n", prog_name);
    printf("  %s compute -S \"c1ccccc1\" -d all\n", prog_name);
    printf("  %s compute -f input.csv -s smiles -d CarbonCount,RingCount -o output.csv\n", prog_name);
    printf("  %s compute -f input.csv -s smiles -d all -o out.csv --no-canonicalization\n", prog_name);
    printf("  %s compute --list\n", prog_name);
}

/* Canonicalize command */
static int cmd_canonicalize(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"smiles",        required_argument, 0, 'S'},
        {"file",          required_argument, 0, 'f'},
        {"col",           required_argument, 0, 's'},
        {"outcol",        required_argument, 0, 'c'},
        {"output",        required_argument, 0, 'o'},
        {"ncpu",          required_argument, 0, 'n'},
        {"sanitize",      required_argument, 0, 1000},
        {"list-tautomers", no_argument,      0, 1001},
        {"verbose",       no_argument,       0, 'v'},
        {"help",          no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* smiles = NULL;
    char* input_file = NULL;
    char* input_col = NULL;
    char* output_col = "canonical_smiles";
    char* output_file = NULL;
    int ncpu = 0;  /* 0 = auto */
    bool verbose = false;
    char* sanitize_str = NULL;
    bool list_tautomers = false;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "S:f:s:c:o:n:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'S':
                smiles = optarg;
                break;
            case 'f':
                input_file = optarg;
                break;
            case 's':
                input_col = optarg;
                break;
            case 'c':
                output_col = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'n':
                ncpu = atoi(optarg);
                break;
            case 1000:
                sanitize_str = optarg;
                break;
            case 1001:
                list_tautomers = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_canonicalize_usage(argv[0]);
                return 0;
            default:
                print_canonicalize_usage(argv[0]);
                return 1;
        }
    }

    /* Single SMILES mode */
    if (smiles) {
        char error_buf[256];

        /* Parse sanitization flags if provided */
        sanitize_flags_t sanitize_flags = SANITIZE_NONE;
        if (sanitize_str) {
            if (sanitize_parse_flags(sanitize_str, &sanitize_flags) != CCHEM_OK) {
                fprintf(stderr, "Error: Invalid sanitize option: %s\n", sanitize_str);
                fprintf(stderr, "Valid options: complete, unsalt, aromatize, kekulize, neutralize,\n");
                fprintf(stderr, "               normalize, remove-stereo, remove-isotopes, remove-h, validate\n");
                return 1;
            }
        }

        /* Handle list-tautomers mode */
        if (list_tautomers) {
            molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
            if (!mol) {
                fprintf(stderr, "Error: %s\n", error_buf);
                return 1;
            }

            /* Apply sanitization if requested */
            if (sanitize_flags != SANITIZE_NONE) {
                sanitize_options_t san_opts = SANITIZE_OPTIONS_DEFAULT;
                san_opts.flags = sanitize_flags;
                if (molecule_sanitize(mol, &san_opts, error_buf, sizeof(error_buf)) != CCHEM_OK) {
                    fprintf(stderr, "Error: %s\n", error_buf);
                    molecule_free(mol);
                    return 1;
                }
            }

            tautomer_result_t result;
            tautomer_options_t taut_opts = TAUTOMER_OPTIONS_DEFAULT;

            if (tautomer_enumerate(mol, &taut_opts, &result) != CCHEM_OK) {
                fprintf(stderr, "Error: Tautomer enumeration failed\n");
                molecule_free(mol);
                return 1;
            }

            if (verbose) {
                printf("Found %d tautomer(s):\n", result.num_tautomers);
            }

            for (int i = 0; i < result.num_tautomers; i++) {
                if (verbose) {
                    printf("%d. %s%s\n", i + 1, result.smiles[i],
                           i == result.canonical_idx ? " (canonical)" : "");
                } else {
                    printf("%s\n", result.smiles[i]);
                }
            }

            tautomer_result_free(&result);
            molecule_free(mol);
            return 0;
        }

        /* Standard canonicalization with optional sanitization */
        if (sanitize_flags != SANITIZE_NONE) {
            /* Parse, sanitize, then canonicalize */
            molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
            if (!mol) {
                fprintf(stderr, "Error: %s\n", error_buf);
                return 1;
            }

            /* Apply sanitization */
            sanitize_options_t san_opts = SANITIZE_OPTIONS_DEFAULT;
            san_opts.flags = sanitize_flags;

            if (molecule_sanitize(mol, &san_opts, error_buf, sizeof(error_buf)) != CCHEM_OK) {
                fprintf(stderr, "Error: %s\n", error_buf);
                molecule_free(mol);
                return 1;
            }

            /* Generate canonical SMILES */
            char* canonical = molecule_to_canonical_smiles(mol, NULL);
            molecule_free(mol);

            if (canonical) {
                printf("%s\n", canonical);
                free(canonical);
                return 0;
            } else {
                fprintf(stderr, "Error: Failed to generate canonical SMILES\n");
                return 1;
            }
        }

        /* Standard canonicalization (no sanitization) */
        char* canonical = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));

        if (canonical) {
            printf("%s\n", canonical);
            free(canonical);
            return 0;
        } else {
            fprintf(stderr, "Error: %s\n", error_buf);
            return 1;
        }
    }

    /* CSV batch mode */
    if (input_file) {
        if (!input_col) {
            fprintf(stderr, "Error: Input column name required (-s/--col)\n");
            return 1;
        }
        if (!output_file) {
            fprintf(stderr, "Error: Output file required (-o/--output)\n");
            return 1;
        }

        if (verbose) {
            printf("Processing: %s\n", input_file);
            printf("Input column: %s\n", input_col);
            printf("Output column: %s\n", output_col);
            printf("Output file: %s\n", output_file);
            printf("CPU cores: %s\n", ncpu > 0 ? "specified" : "auto");
            printf("\n");
        }

        csv_batch_context_t* ctx = csv_batch_context_create(
            input_file, output_file, input_col, output_col, ncpu);

        if (!ctx) {
            fprintf(stderr, "Error: Failed to create batch context\n");
            return 1;
        }

        cchem_status_t status = csv_batch_canonicalize(ctx);

        if (status == CCHEM_OK) {
            size_t total, success, errors;
            csv_batch_get_stats(ctx, &total, &success, &errors);

            if (verbose) {
                printf("\n");
            }
            printf("Processed: %zu/%zu (%.1f%% success)\n",
                   success, total, total > 0 ? (double)success / total * 100.0 : 0.0);

            if (errors > 0) {
                printf("Errors: %zu\n", errors);
            }

            printf("Output written to: %s\n", output_file);

            csv_batch_context_free(ctx);
            return 0;
        } else {
            fprintf(stderr, "Error: Batch processing failed\n");
            csv_batch_context_free(ctx);
            return 1;
        }
    }

    /* No input specified */
    fprintf(stderr, "Error: No input specified. Use -S for single SMILES or -f for CSV file.\n");
    print_canonicalize_usage(argv[0]);
    return 1;
}

/* Validate command */
static int cmd_validate(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"smiles",  required_argument, 0, 'S'},
        {"verbose", no_argument,       0, 'v'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* smiles = NULL;
    bool verbose = false;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "S:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'S':
                smiles = optarg;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_validate_usage(argv[0]);
                return 0;
            default:
                print_validate_usage(argv[0]);
                return 1;
        }
    }

    if (!smiles) {
        fprintf(stderr, "Error: SMILES string required (-S/--smiles)\n");
        return 1;
    }

    char error_buf[256];
    cchem_status_t status = cchem_validate_smiles(smiles, error_buf, sizeof(error_buf));

    if (status == CCHEM_OK) {
        if (verbose) {
            printf("Valid SMILES: %s\n", smiles);
        } else {
            printf("valid\n");
        }
        return 0;
    } else {
        if (verbose) {
            printf("Invalid SMILES: %s\n", smiles);
            printf("Error: %s\n", error_buf);
        } else {
            printf("invalid: %s\n", error_buf);
        }
        return 1;
    }
}


/* SMILES column data loaded via mmap (more efficient than csv_document_read) */
typedef struct {
    char* buffer;           /* Contiguous buffer holding all SMILES strings */
    char** smiles;          /* Array of pointers into buffer */
    char** header_fields;   /* Header field names (for output) */
    int num_header_fields;
    int num_rows;
    int smiles_col_idx;
} smiles_data_t;

/* Free SMILES data */
static void smiles_data_free(smiles_data_t* data) {
    if (data) {
        if (data->buffer) free(data->buffer);
        if (data->smiles) free(data->smiles);
        if (data->header_fields) {
            for (int i = 0; i < data->num_header_fields; i++) {
                if (data->header_fields[i]) free(data->header_fields[i]);
            }
            free(data->header_fields);
        }
        free(data);
    }
}

/* Load SMILES column from CSV using mmap (much faster than csv_document_read)
 * Note: Currently unused but available for future optimization */
__attribute__((unused))
static smiles_data_t* smiles_data_load(const char* filename, const char* smiles_col_name) {
    csv_mmap_reader_t* reader = csv_mmap_reader_create(filename);
    if (!reader) return NULL;

    smiles_data_t* data = (smiles_data_t*)calloc(1, sizeof(smiles_data_t));
    if (!data) {
        csv_mmap_reader_free(reader);
        return NULL;
    }

    /* First pass: read header and count rows */
    size_t line_len;
    const char* line = csv_mmap_reader_next_line(reader, &line_len);
    if (!line) {
        csv_mmap_reader_free(reader);
        smiles_data_free(data);
        return NULL;
    }

    /* Parse header */
    const char* header_fields[256];
    int num_header = csv_mmap_parse_line(line, line_len, ',', header_fields, 256);
    if (num_header <= 0) {
        csv_mmap_reader_free(reader);
        smiles_data_free(data);
        return NULL;
    }

    /* Copy header fields (need to persist) */
    data->num_header_fields = num_header;
    data->header_fields = (char**)calloc(num_header, sizeof(char*));
    for (int i = 0; i < num_header; i++) {
        data->header_fields[i] = strdup(header_fields[i]);
    }

    /* Find SMILES column */
    data->smiles_col_idx = -1;
    for (int i = 0; i < num_header; i++) {
        if (strcmp(header_fields[i], smiles_col_name) == 0) {
            data->smiles_col_idx = i;
            break;
        }
    }
    if (data->smiles_col_idx < 0) {
        csv_mmap_reader_free(reader);
        smiles_data_free(data);
        return NULL;
    }

    /* Count rows and total SMILES length */
    size_t total_smiles_len = 0;
    int row_count = 0;

    while ((line = csv_mmap_reader_next_line(reader, &line_len))) {
        const char* fields[256];
        int n = csv_mmap_parse_line(line, line_len, ',', fields, 256);
        if (n > data->smiles_col_idx && fields[data->smiles_col_idx]) {
            total_smiles_len += strlen(fields[data->smiles_col_idx]) + 1;
        } else {
            total_smiles_len += 1;  /* Empty string */
        }
        row_count++;
    }

    data->num_rows = row_count;

    /* Allocate contiguous buffer for all SMILES strings */
    data->buffer = (char*)malloc(total_smiles_len);
    data->smiles = (char**)malloc(row_count * sizeof(char*));
    if (!data->buffer || !data->smiles) {
        csv_mmap_reader_free(reader);
        smiles_data_free(data);
        return NULL;
    }

    /* Second pass: read header again to reset position */
    csv_mmap_reader_free(reader);
    reader = csv_mmap_reader_create(filename);
    csv_mmap_reader_next_line(reader, &line_len);  /* Skip header */

    /* Extract SMILES strings into contiguous buffer */
    char* buf_ptr = data->buffer;
    int row_idx = 0;
    while ((line = csv_mmap_reader_next_line(reader, &line_len)) && row_idx < row_count) {
        const char* fields[256];
        int n = csv_mmap_parse_line(line, line_len, ',', fields, 256);

        const char* smiles = (n > data->smiles_col_idx) ? fields[data->smiles_col_idx] : "";
        size_t slen = strlen(smiles);

        data->smiles[row_idx] = buf_ptr;
        memcpy(buf_ptr, smiles, slen);
        buf_ptr[slen] = '\0';
        buf_ptr += slen + 1;
        row_idx++;
    }

    csv_mmap_reader_free(reader);
    return data;
}

/* Fixed-width field for pre-allocated buffer (24 bytes covers most descriptor values) */
#define DESC_VALUE_WIDTH 24

/* Fast integer to string - avoids snprintf overhead for integers */
static inline void fast_i64toa(int64_t val, char* buf, int buf_size) {
    if (buf_size < 2) return;

    char tmp[24];
    int i = 0;
    bool neg = val < 0;
    if (neg) val = -val;

    do {
        tmp[i++] = '0' + (val % 10);
        val /= 10;
    } while (val > 0 && i < 22);

    int j = 0;
    if (neg && j < buf_size - 1) buf[j++] = '-';
    while (i > 0 && j < buf_size - 1) {
        buf[j++] = tmp[--i];
    }
    buf[j] = '\0';
}

/* Fast double to string with 6 significant figures (like %.6g) */
static inline void fast_dtoa(double val, char* buf, int buf_size) {
    if (buf_size < 2) { buf[0] = '\0'; return; }

    /* Handle special cases - replace NaN/Inf with 0 for clean CSV output */
    if (val != val) { /* NaN */
        buf[0] = '0'; buf[1] = '\0';
        return;
    }
    /* Check for infinity (works with -ffast-math) */
    if (val > 1e300 || val < -1e300) {
        buf[0] = '0'; buf[1] = '\0';
        return;
    }
    if (val == 0.0) {
        buf[0] = '0'; buf[1] = '\0';
        return;
    }

    int pos = 0;
    if (val < 0) {
        buf[pos++] = '-';
        val = -val;
    }

    /* For very small or very large numbers, fall back to snprintf */
    if (val < 1e-4 || val >= 1e7) {
        int n = snprintf(buf, buf_size, "%.6g", (pos > 0) ? -val : val);
        /* Ensure the result is valid (snprintf might produce nan/inf in edge cases) */
        if (n > 0 && n < buf_size) {
            char c = buf[pos > 0 ? 1 : 0];  /* First char after optional sign */
            if (c == 'n' || c == 'i' || c == 'N' || c == 'I') {
                /* Replace nan/inf with 0 */
                buf[0] = '0'; buf[1] = '\0';
            }
        }
        return;
    }

    /* Fixed-point formatting for normal range */
    int64_t int_part = (int64_t)val;
    double frac_part = val - int_part;

    /* Write integer part */
    char tmp[24];
    int i = 0;
    int64_t ip = int_part;
    if (ip == 0) {
        tmp[i++] = '0';
    } else {
        while (ip > 0 && i < 20) {
            tmp[i++] = '0' + (ip % 10);
            ip /= 10;
        }
    }
    while (i > 0 && pos < buf_size - 1) {
        buf[pos++] = tmp[--i];
    }

    /* Write fractional part (up to 6 digits) */
    if (frac_part > 0.0000005 && pos < buf_size - 2) {
        buf[pos++] = '.';
        int frac_digits = 0;
        while (frac_part > 0.0000005 && frac_digits < 6 && pos < buf_size - 1) {
            frac_part *= 10;
            int digit = (int)frac_part;
            buf[pos++] = '0' + digit;
            frac_part -= digit;
            frac_digits++;
        }
        /* Remove trailing zeros */
        while (pos > 1 && buf[pos-1] == '0') pos--;
        if (buf[pos-1] == '.') pos--;
    }
    buf[pos] = '\0';
}


/* ============================================================================
 * Pipeline Streaming Data Structures
 * ============================================================================ */

/* Input work item - from reader to workers */
typedef struct {
    int row_idx;
    char* smiles;                 /* Owned, freed after processing */
    char** original_fields;       /* Owned array of owned strings */
    int num_original_fields;
} work_item_t;

/* Output work item - from workers to writer */
typedef struct {
    int row_idx;                  /* For potential ordering */
    char** output_fields;         /* Ready-to-write field array (owned) */
    int num_fields;
    bool success;
} result_item_t;

/* Pipeline context */
typedef struct {
    /* Input */
    csv_mmap_reader_t* reader;
    int smiles_col_idx;
    int num_input_fields;
    char** header_fields;         /* Original header field names */

    /* Output */
    FILE* output_file;
    char** desc_names;
    int num_descs;

    /* Queues */
    bounded_queue_t* input_queue;
    bounded_queue_t* output_queue;

    /* Threads */
    pthread_t reader_thread;
    pthread_t writer_thread;
    pthread_t* worker_threads;
    int num_workers;

    /* Progress */
    atomic_size_t rows_read;
    atomic_size_t rows_processed;
    atomic_size_t rows_written;
    size_t total_rows;

    /* Control */
    bool skip_canon;
    char delimiter;
} pipeline_ctx_t;

/* Free a work item */
static void work_item_free(work_item_t* item) {
    if (!item) return;
    free(item->smiles);
    if (item->original_fields) {
        for (int i = 0; i < item->num_original_fields; i++) {
            free(item->original_fields[i]);
        }
        free(item->original_fields);
    }
    free(item);
}

/* Free a result item */
static void result_item_free(result_item_t* result) {
    if (!result) return;
    if (result->output_fields) {
        for (int i = 0; i < result->num_fields; i++) {
            free(result->output_fields[i]);
        }
        free(result->output_fields);
    }
    free(result);
}

/* Reader thread - reads CSV rows and pushes to input queue */
static void* pipeline_reader_thread(void* arg) {
    pipeline_ctx_t* ctx = (pipeline_ctx_t*)arg;

    size_t line_len;
    const char* line;
    int row_idx = 0;
    const char* fields[256];

    while ((line = csv_mmap_reader_next_line(ctx->reader, &line_len))) {
        int nfields = csv_mmap_parse_line(line, line_len, ctx->delimiter, fields, 256);

        /* Create work item */
        work_item_t* item = (work_item_t*)malloc(sizeof(work_item_t));
        if (!item) continue;

        item->row_idx = row_idx++;
        item->num_original_fields = nfields;

        /* Copy SMILES string */
        if (ctx->smiles_col_idx < nfields && fields[ctx->smiles_col_idx]) {
            item->smiles = strdup(fields[ctx->smiles_col_idx]);
        } else {
            item->smiles = strdup("");
        }

        /* Copy all original fields */
        item->original_fields = (char**)calloc(nfields, sizeof(char*));
        if (item->original_fields) {
            for (int i = 0; i < nfields; i++) {
                item->original_fields[i] = fields[i] ? strdup(fields[i]) : strdup("");
            }
        }

        /* Push to input queue (blocks if full) */
        if (!queue_push(ctx->input_queue, item)) {
            work_item_free(item);
            break;
        }
        atomic_fetch_add(&ctx->rows_read, 1);
    }

    /* Signal EOF */
    queue_close(ctx->input_queue);
    return NULL;
}

/* Worker thread - processes SMILES and computes descriptors */
static void* pipeline_worker_thread(void* arg) {
    pipeline_ctx_t* ctx = (pipeline_ctx_t*)arg;

    /* Thread-local molecule for reuse */
    molecule_t* tl_mol = molecule_create_with_capacity(256, 256);

    /* Cache descriptor definitions */
    const descriptor_def_t* cached_defs[MAX_DESCRIPTORS];
    int cached_num_defs = descriptor_get_all(cached_defs, MAX_DESCRIPTORS);
    int cached_total = descriptor_count();
    bool compute_all = (ctx->num_descs >= cached_total);

    char error_buf[256];
    work_item_t* item;

    while ((item = (work_item_t*)queue_pop(ctx->input_queue)) != NULL) {
        /* Create result with all fields */
        int out_num_fields = item->num_original_fields + ctx->num_descs;
        result_item_t* result = (result_item_t*)malloc(sizeof(result_item_t));
        if (!result) {
            work_item_free(item);
            continue;
        }

        result->row_idx = item->row_idx;
        result->num_fields = out_num_fields;
        result->output_fields = (char**)calloc(out_num_fields, sizeof(char*));
        result->success = false;

        if (!result->output_fields) {
            result_item_free(result);
            work_item_free(item);
            continue;
        }

        /* Copy original fields */
        for (int i = 0; i < item->num_original_fields; i++) {
            result->output_fields[i] = item->original_fields[i] ?
                                       strdup(item->original_fields[i]) : strdup("");
        }

        /* Compute descriptors */
        molecule_t* mol = NULL;
        cchem_status_t status;

        if (item->smiles && item->smiles[0] != '\0') {
            if (ctx->skip_canon) {
                status = smiles_to_molecule_reuse(tl_mol, item->smiles, error_buf, sizeof(error_buf));
                mol = (status == CCHEM_OK) ? tl_mol : NULL;
            } else {
                char* canonical = smiles_canonicalize(item->smiles, NULL, error_buf, sizeof(error_buf));
                if (canonical) {
                    status = smiles_to_molecule_reuse(tl_mol, canonical, error_buf, sizeof(error_buf));
                    mol = (status == CCHEM_OK) ? tl_mol : NULL;
                    free(canonical);
                }
            }
        }

        if (mol) {
            result->success = true;

            if (compute_all) {
                /* Batch computation */
                descriptor_value_t all_values[MAX_DESCRIPTORS];
                int num_computed = descriptors_compute_all(mol, NULL, all_values, MAX_DESCRIPTORS);

                int n = (ctx->num_descs < num_computed) ? ctx->num_descs : num_computed;
                if (n > cached_num_defs) n = cached_num_defs;

                for (int i = 0; i < n; i++) {
                    char val_buf[DESC_VALUE_WIDTH];
                    if (cached_defs[i] == NULL) {
                        strcpy(val_buf, "0");
                    } else if (cached_defs[i]->value_type == DESC_VALUE_INT) {
                        fast_i64toa(all_values[i].i, val_buf, DESC_VALUE_WIDTH);
                    } else {
                        fast_dtoa(all_values[i].d, val_buf, DESC_VALUE_WIDTH);
                    }
                    result->output_fields[item->num_original_fields + i] = strdup(val_buf);
                }
                /* Fill remaining with "0" */
                for (int i = n; i < ctx->num_descs; i++) {
                    result->output_fields[item->num_original_fields + i] = strdup("0");
                }
            } else {
                /* Individual computation for subset */
                for (int i = 0; i < ctx->num_descs; i++) {
                    char val_buf[DESC_VALUE_WIDTH];
                    strcpy(val_buf, "0");

                    const descriptor_def_t* def = descriptor_get(ctx->desc_names[i]);
                    if (def) {
                        descriptor_value_t value;
                        memset(&value, 0, sizeof(value));
                        if (def->compute(mol, &value) == CCHEM_OK) {
                            if (def->value_type == DESC_VALUE_INT) {
                                fast_i64toa(value.i, val_buf, DESC_VALUE_WIDTH);
                            } else {
                                fast_dtoa(value.d, val_buf, DESC_VALUE_WIDTH);
                            }
                        }
                    }
                    result->output_fields[item->num_original_fields + i] = strdup(val_buf);
                }
            }
        } else {
            /* Failed to parse - fill with "0" */
            for (int i = 0; i < ctx->num_descs; i++) {
                result->output_fields[item->num_original_fields + i] = strdup("0");
            }
        }

        atomic_fetch_add(&ctx->rows_processed, 1);

        /* Push to output queue */
        if (!queue_push(ctx->output_queue, result)) {
            result_item_free(result);
        }

        work_item_free(item);
    }

    molecule_free(tl_mol);
    return NULL;
}

/* Writer thread - writes results to CSV */
static void* pipeline_writer_thread(void* arg) {
    pipeline_ctx_t* ctx = (pipeline_ctx_t*)arg;

    result_item_t* result;
    char line_buffer[65536];  /* 64KB line buffer */

    while ((result = (result_item_t*)queue_pop(ctx->output_queue)) != NULL) {
        /* Format and write row */
        int pos = 0;
        for (int i = 0; i < result->num_fields; i++) {
            if (i > 0) {
                line_buffer[pos++] = ctx->delimiter;
            }
            const char* field = result->output_fields[i] ? result->output_fields[i] : "";

            /* Check if field needs quoting */
            bool needs_quote = false;
            for (const char* p = field; *p; p++) {
                if (*p == ctx->delimiter || *p == '"' || *p == '\n' || *p == '\r') {
                    needs_quote = true;
                    break;
                }
            }

            if (needs_quote) {
                line_buffer[pos++] = '"';
                for (const char* p = field; *p && pos < 65530; p++) {
                    if (*p == '"') {
                        line_buffer[pos++] = '"';
                    }
                    line_buffer[pos++] = *p;
                }
                line_buffer[pos++] = '"';
            } else {
                size_t flen = strlen(field);
                if (pos + flen < 65530) {
                    memcpy(line_buffer + pos, field, flen);
                    pos += flen;
                }
            }
        }
        line_buffer[pos++] = '\n';
        line_buffer[pos] = '\0';

        fwrite(line_buffer, 1, pos, ctx->output_file);
        atomic_fetch_add(&ctx->rows_written, 1);

        result_item_free(result);
    }

    return NULL;
}

/* Pipeline compute - streaming with constant memory */
static int cmd_compute_pipeline(const char* input_file, const char* output_file,
                                const char* input_col, char** desc_names,
                                int num_descs, int ncpu, bool skip_canon, bool verbose) {
    pipeline_ctx_t ctx = {0};
    ctx.delimiter = ',';
    ctx.skip_canon = skip_canon;
    ctx.desc_names = desc_names;
    ctx.num_descs = num_descs;

    /* Count total rows for progress */
    ctx.total_rows = csv_count_rows(input_file);
    if (ctx.total_rows <= 1) {
        fprintf(stderr, "Error: Input file is empty or has only header\n");
        return 1;
    }
    ctx.total_rows--;  /* Exclude header */

    /* Create mmap reader */
    ctx.reader = csv_mmap_reader_create(input_file);
    if (!ctx.reader) {
        fprintf(stderr, "Error: Failed to open input file: %s\n", input_file);
        return 1;
    }

    /* Read and parse header */
    size_t header_len;
    const char* header_line = csv_mmap_reader_next_line(ctx.reader, &header_len);
    if (!header_line) {
        fprintf(stderr, "Error: Failed to read header\n");
        csv_mmap_reader_free(ctx.reader);
        return 1;
    }

    const char* header_fields[256];
    ctx.num_input_fields = csv_mmap_parse_line(header_line, header_len, ctx.delimiter, header_fields, 256);

    /* Find SMILES column */
    ctx.smiles_col_idx = -1;
    for (int i = 0; i < ctx.num_input_fields; i++) {
        if (strcmp(header_fields[i], input_col) == 0) {
            ctx.smiles_col_idx = i;
            break;
        }
    }
    if (ctx.smiles_col_idx < 0) {
        fprintf(stderr, "Error: Column '%s' not found in input file\n", input_col);
        csv_mmap_reader_free(ctx.reader);
        return 1;
    }

    /* Copy header field names */
    ctx.header_fields = (char**)calloc(ctx.num_input_fields, sizeof(char*));
    for (int i = 0; i < ctx.num_input_fields; i++) {
        ctx.header_fields[i] = strdup(header_fields[i]);
    }

    /* Open output file */
    ctx.output_file = fopen(output_file, "w");
    if (!ctx.output_file) {
        fprintf(stderr, "Error: Failed to create output file: %s\n", output_file);
        for (int i = 0; i < ctx.num_input_fields; i++) free(ctx.header_fields[i]);
        free(ctx.header_fields);
        csv_mmap_reader_free(ctx.reader);
        return 1;
    }
    setvbuf(ctx.output_file, NULL, _IOFBF, 262144);  /* 256KB buffer */

    /* Write output header */
    int out_num_fields = ctx.num_input_fields + num_descs;
    const char** out_header = (const char**)calloc(out_num_fields, sizeof(char*));
    for (int i = 0; i < ctx.num_input_fields; i++) {
        out_header[i] = ctx.header_fields[i];
    }
    for (int i = 0; i < num_descs; i++) {
        out_header[ctx.num_input_fields + i] = desc_names[i];
    }
    csv_write_header_line(ctx.output_file, out_header, out_num_fields, ctx.delimiter);
    free(out_header);

    /* Determine thread count */
    ctx.num_workers = (ncpu > 0) ? ncpu : parallel_get_num_cores();

    if (verbose) {
        printf("Pipeline mode: %zu molecules, %d workers\n", ctx.total_rows, ctx.num_workers);
    }

    /* Create bounded queues */
    int queue_size = ctx.num_workers * 2;
    ctx.input_queue = queue_create(queue_size);
    ctx.output_queue = queue_create(queue_size);
    if (!ctx.input_queue || !ctx.output_queue) {
        fprintf(stderr, "Error: Failed to create queues\n");
        if (ctx.input_queue) queue_free(ctx.input_queue);
        if (ctx.output_queue) queue_free(ctx.output_queue);
        fclose(ctx.output_file);
        for (int i = 0; i < ctx.num_input_fields; i++) free(ctx.header_fields[i]);
        free(ctx.header_fields);
        csv_mmap_reader_free(ctx.reader);
        return 1;
    }

    /* Initialize atomics */
    atomic_store(&ctx.rows_read, 0);
    atomic_store(&ctx.rows_processed, 0);
    atomic_store(&ctx.rows_written, 0);

    printf("Computing %d descriptors for %zu molecules...\n", num_descs, ctx.total_rows);

    /* Create progress bar */
    progress_config_t prog_config = PROGRESS_CONFIG_DEFAULT;
    prog_config.prefix = "Processing";
    progress_t* progress = progress_create(ctx.total_rows, &prog_config);

    /* Start reader thread */
    pthread_create(&ctx.reader_thread, NULL, pipeline_reader_thread, &ctx);

    /* Start worker threads */
    ctx.worker_threads = (pthread_t*)malloc(ctx.num_workers * sizeof(pthread_t));
    for (int i = 0; i < ctx.num_workers; i++) {
        pthread_create(&ctx.worker_threads[i], NULL, pipeline_worker_thread, &ctx);
    }

    /* Start writer thread */
    pthread_create(&ctx.writer_thread, NULL, pipeline_writer_thread, &ctx);

    /* Progress monitoring (main thread) */
    while (atomic_load(&ctx.rows_written) < ctx.total_rows) {
        size_t written = atomic_load(&ctx.rows_written);
        progress_update(progress, written);
        usleep(50000);  /* 50ms */
    }

    /* Wait for reader to finish */
    pthread_join(ctx.reader_thread, NULL);

    /* Wait for all workers to finish */
    for (int i = 0; i < ctx.num_workers; i++) {
        pthread_join(ctx.worker_threads[i], NULL);
    }

    /* Signal output queue EOF and wait for writer */
    queue_close(ctx.output_queue);
    pthread_join(ctx.writer_thread, NULL);

    progress_finish(progress);
    progress_free(progress);

    /* Cleanup */
    queue_free(ctx.input_queue);
    queue_free(ctx.output_queue);
    free(ctx.worker_threads);
    fclose(ctx.output_file);
    for (int i = 0; i < ctx.num_input_fields; i++) free(ctx.header_fields[i]);
    free(ctx.header_fields);
    csv_mmap_reader_free(ctx.reader);

    size_t rows_written = atomic_load(&ctx.rows_written);
    fprintf(stderr, "\033[KProcessed: %zu/%zu (%.1f%% success)\n",
           rows_written, ctx.total_rows,
           ctx.total_rows > 0 ? (double)rows_written / ctx.total_rows * 100.0 : 0.0);
    fprintf(stderr, "Output written to: %s\n", output_file);

    return 0;
}

/* Compute command */
static int cmd_compute(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"smiles",              required_argument, 0, 'S'},
        {"file",                required_argument, 0, 'f'},
        {"col",                 required_argument, 0, 's'},
        {"descriptors",         required_argument, 0, 'd'},
        {"output",              required_argument, 0, 'o'},
        {"ncpu",                required_argument, 0, 'n'},
        {"no-canonicalization", no_argument,       0, 1000},
        {"list",                no_argument,       0, 'l'},
        {"verbose",             no_argument,       0, 'v'},
        {"help",                no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* smiles = NULL;
    char* input_file = NULL;
    char* input_col = NULL;
    char* desc_list = "all";
    char* output_file = NULL;
    int ncpu = 0;  /* 0 = auto */
    bool list_descriptors = false;
    bool verbose = false;
    bool skip_canonicalization = false;
    (void)argc;

    int opt;
    int option_index = 0;

    /* Initialize descriptor registry */
    descriptors_init();

    while ((opt = getopt_long(argc, argv, "S:f:s:d:o:n:lvh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'S':
                smiles = optarg;
                break;
            case 'f':
                input_file = optarg;
                break;
            case 's':
                input_col = optarg;
                break;
            case 'd':
                desc_list = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'n':
                ncpu = atoi(optarg);
                if (ncpu < 0) ncpu = 0;
                break;
            case 1000:
                skip_canonicalization = true;
                break;
            case 'l':
                list_descriptors = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_compute_usage(argv[0]);
                return 0;
            default:
                print_compute_usage(argv[0]);
                return 1;
        }
    }

    /* List descriptors mode */
    if (list_descriptors) {
        descriptor_list_all();
        return 0;
    }

    /* Parse descriptor names */
    char* desc_names[MAX_DESCRIPTORS];
    int num_descs = descriptor_parse_names(desc_list, desc_names, MAX_DESCRIPTORS);
    if (num_descs <= 0) {
        fprintf(stderr, "Error: Invalid descriptor specification: %s\n", desc_list);
        fprintf(stderr, "Use --list to see available descriptors\n");
        return 1;
    }

    /* Single SMILES mode */
    if (smiles) {
        char error_buf[256];
        molecule_t* mol = NULL;
        char* canonical = NULL;

        if (skip_canonicalization) {
            /* Parse directly without canonicalization */
            mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
            if (!mol) {
                fprintf(stderr, "Error: %s\n", error_buf);
                descriptor_free_names(desc_names, num_descs);
                return 1;
            }
        } else {
            /* Canonicalize first */
            canonical = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));
            if (!canonical) {
                fprintf(stderr, "Error: %s\n", error_buf);
                descriptor_free_names(desc_names, num_descs);
                return 1;
            }

            /* Parse into molecule */
            mol = smiles_to_molecule(canonical, error_buf, sizeof(error_buf));
            if (!mol) {
                fprintf(stderr, "Error: %s\n", error_buf);
                free(canonical);
                descriptor_free_names(desc_names, num_descs);
                return 1;
            }
        }

        /* Compute and print descriptors */
        if (verbose) {
            printf("SMILES: %s\n", smiles);
            if (canonical) {
                printf("Canonical: %s\n", canonical);
            } else {
                printf("(canonicalization skipped)\n");
            }
            printf("\n");
        }

        for (int i = 0; i < num_descs; i++) {
            const descriptor_def_t* def = descriptor_get(desc_names[i]);
            if (!def) continue;

            descriptor_value_t value;
            memset(&value, 0, sizeof(value));  /* Zero-initialize entire union */
            cchem_status_t status = def->compute(mol, &value);

            if (status == CCHEM_OK) {
                char val_buf[64];
                descriptor_format_value(def, &value, val_buf, sizeof(val_buf));
                printf("%s: %s\n", def->name, val_buf);
            } else {
                printf("%s: ERROR\n", def->name);
            }
        }

        molecule_free(mol);
        free(canonical);
        descriptor_free_names(desc_names, num_descs);
        return 0;
    }

    /* CSV batch mode - use pipeline streaming for constant memory usage */
    if (input_file) {
        if (!input_col) {
            fprintf(stderr, "Error: Input column name required (-s/--col)\n");
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }
        if (!output_file) {
            fprintf(stderr, "Error: Output file required (-o/--output)\n");
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        if (verbose) {
            printf("Processing: %s\n", input_file);
            printf("Input column: %s\n", input_col);
            printf("Descriptors: %s\n", desc_list);
            printf("Output file: %s\n", output_file);
            printf("CPU cores: %s\n", ncpu > 0 ? "specified" : "auto");
            if (skip_canonicalization) {
                printf("Canonicalization: skipped\n");
            }
            printf("\n");
        }

        /* Use pipeline streaming for constant memory regardless of input size */
        int result = cmd_compute_pipeline(input_file, output_file, input_col,
                                          desc_names, num_descs, ncpu,
                                          skip_canonicalization, verbose);
        descriptor_free_names(desc_names, num_descs);
        return result;
    }

    /* No input specified */
    fprintf(stderr, "Error: No input specified. Use -S for single SMILES or -f for CSV file.\n");
    print_compute_usage(argv[0]);
    descriptor_free_names(desc_names, num_descs);
    return 1;
}

static void print_depict_usage(const char* prog_name) {
    printf("Usage: %s depict [options]\n\n", prog_name);
    printf("Generate publication-quality molecular structure images from SMILES\n\n");
    printf("Required:\n");
    printf("  -S, --smiles <SMILES>   SMILES string to depict\n");
    printf("  -o, --output <file>     Output file (.png, .jpg, .svg)\n");
    printf("\n");
    printf("Render Style:\n");
    printf("  -s, --style <style>     wireframe (default), sticks, balls-sticks, spacefill, surface\n");
    printf("  -m, --mode <mode>       2d (default) or 3d (MMFF94 optimized geometry)\n");
    printf("\n");
    printf("Dimensions:\n");
    printf("  -W, --width <pixels>    Image width (default: 800)\n");
    printf("  -H, --height <pixels>   Image height (default: 800)\n");
    printf("  --margin <px>           Margin (default: 50)\n");
    printf("  --scale <factor>        Resolution multiplier (e.g., 2.0 for 2x)\n");
    printf("\n");
    printf("Bond Rendering:\n");
    printf("  --bond-length <px>      Bond length (default: 35)\n");
    printf("  --bond-width <px>       Bond line width (default: 1.8)\n");
    printf("  --heteroatom-gap <0-1>  Gap at heteroatoms (default: 1.0, 0 to disable)\n");
    printf("  --line-cap <style>      round, butt (default), square\n");
    printf("\n");
    printf("Atom Labels:\n");
    printf("  --show-carbons          Show 'C' labels on carbons\n");
    printf("  --show-hydrogens        Show hydrogen labels\n");
    printf("  --terminal-carbons      Show CH3 on terminal carbons\n");
    printf("  --font-size <scale>     Font size scale (default: 3.5)\n");
    printf("\n");
    printf("3D/Surface Options:\n");
    printf("  --atom-filling          Draw atoms as filled CPK spheres\n");
    printf("  --proportional-atoms    Scale atoms by VDW radius (default: on)\n");
    printf("  --no-proportional       Disable proportional sizing\n");
    printf("  --surface-color <mode>  uniform, atom, polarity (for -s surface)\n");
    printf("  --max-iter <N>          3D optimization iterations (default: 500)\n");
    printf("\n");
    printf("Other:\n");
    printf("  --toggle-aromaticity    Draw aromatic circles\n");
    printf("  --quality <1-100>       JPEG quality (default: 95)\n");
    printf("  -v, --verbose           Verbose output\n");
    printf("  -h, --help              Print this help\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s depict -S \"c1ccccc1\" -o benzene.png\n", prog_name);
    printf("  %s depict -S \"CCO\" -o ethanol.svg\n", prog_name);
    printf("  %s depict -S \"CCO\" -o ethanol.png -m 3d -s balls-sticks\n", prog_name);
    printf("  %s depict -S \"CCO\" -o ethanol.png --scale 2\n", prog_name);
}

/* Depict command */
static int cmd_depict(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"smiles",         required_argument, 0, 'S'},
        {"output",         required_argument, 0, 'o'},
        {"mode",           required_argument, 0, 'm'},
        {"style",          required_argument, 0, 's'},
        {"width",          required_argument, 0, 'W'},
        {"height",         required_argument, 0, 'H'},
        {"bond-length",    required_argument, 0, 1000},
        {"bond-width",     required_argument, 0, 1001},
        {"margin",         required_argument, 0, 1002},
        {"show-carbons",   no_argument,       0, 1003},
        {"show-hydrogens", no_argument,       0, 1004},
        {"toggle-aromaticity", no_argument,   0, 1005},
        {"atom-filling",   no_argument,       0, 1006},
        {"font-size",      required_argument, 0, 1007},
        {"quality",        required_argument, 0, 1008},
        {"max-iter",       required_argument, 0, 1009},
        {"proportional-atoms", no_argument,   0, 1010},
        {"no-proportional", no_argument,      0, 1011},
        {"surface-color",  required_argument, 0, 1012},
        {"modern",         no_argument,       0, 1013},
        {"heteroatom-gap", required_argument, 0, 1014},
        {"scale",          required_argument, 0, 1015},
        {"line-cap",       required_argument, 0, 1016},
        {"terminal-carbons", no_argument,     0, 1017},
        {"verbose",        no_argument,       0, 'v'},
        {"help",           no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* smiles = NULL;
    char* output_file = NULL;
    bool verbose = false;

    depictor_options_t options = DEPICTOR_OPTIONS_DEFAULT;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "S:o:m:s:W:H:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'S':
                smiles = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'm':
                if (strcmp(optarg, "2d") == 0 || strcmp(optarg, "2D") == 0) {
                    options.mode = DEPICT_MODE_2D;
                } else if (strcmp(optarg, "3d") == 0 || strcmp(optarg, "3D") == 0) {
                    options.mode = DEPICT_MODE_3D;
                } else {
                    fprintf(stderr, "Error: Invalid mode '%s'. Use '2d' or '3d'.\n", optarg);
                    return 1;
                }
                break;
            case 's':  /* --style */
                if (strcmp(optarg, "wireframe") == 0) {
                    options.render_style = RENDER_STYLE_WIREFRAME;
                } else if (strcmp(optarg, "sticks") == 0) {
                    options.render_style = RENDER_STYLE_STICKS;
                } else if (strcmp(optarg, "balls-sticks") == 0 || strcmp(optarg, "balls") == 0) {
                    options.render_style = RENDER_STYLE_BALLS_AND_STICKS;
                } else if (strcmp(optarg, "spacefill") == 0 || strcmp(optarg, "cpk") == 0) {
                    options.render_style = RENDER_STYLE_SPACEFILL;
                } else if (strcmp(optarg, "surface") == 0) {
                    options.render_style = RENDER_STYLE_SURFACE;
                } else {
                    fprintf(stderr, "Error: Invalid style '%s'. Use wireframe, sticks, balls-sticks, spacefill, or surface.\n", optarg);
                    return 1;
                }
                break;
            case 'W':
                options.width = atoi(optarg);
                if (options.width < 50 || options.width > 4096) {
                    fprintf(stderr, "Error: Width must be between 50 and 4096.\n");
                    return 1;
                }
                break;
            case 'H':
                options.height = atoi(optarg);
                if (options.height < 50 || options.height > 4096) {
                    fprintf(stderr, "Error: Height must be between 50 and 4096.\n");
                    return 1;
                }
                break;
            case 1000:  /* --bond-length */
                options.bond_length = atof(optarg);
                break;
            case 1001:  /* --bond-width */
                options.bond_width = atof(optarg);
                break;
            case 1002:  /* --margin */
                options.margin = atoi(optarg);
                break;
            case 1003:  /* --show-carbons */
                options.show_carbons = true;
                break;
            case 1004:  /* --show-hydrogens */
                options.show_hydrogens = true;
                break;
            case 1005:  /* --toggle-aromaticity */
                options.draw_aromatic_circles = true;
                break;
            case 1006:  /* --atom-filling */
                options.atom_filling = true;
                break;
            case 1007:  /* --font-size */
                options.font_size = atof(optarg);
                if (options.font_size < 0.5) options.font_size = 0.5;
                if (options.font_size > 10.0) options.font_size = 10.0;
                break;
            case 1008:  /* --quality */
                options.jpeg_quality = atoi(optarg);
                if (options.jpeg_quality < 1) options.jpeg_quality = 1;
                if (options.jpeg_quality > 100) options.jpeg_quality = 100;
                break;
            case 1009:  /* --max-iter */
                options.max_iterations = atoi(optarg);
                break;
            case 1010:  /* --proportional-atoms */
                options.proportional_atoms = true;
                break;
            case 1011:  /* --no-proportional */
                options.proportional_atoms = false;
                break;
            case 1012:  /* --surface-color */
                if (strcmp(optarg, "uniform") == 0 || strcmp(optarg, "blue") == 0) {
                    options.surface_color = SURFACE_COLOR_UNIFORM;
                } else if (strcmp(optarg, "atom") == 0 || strcmp(optarg, "element") == 0) {
                    options.surface_color = SURFACE_COLOR_ATOM;
                } else if (strcmp(optarg, "polarity") == 0 || strcmp(optarg, "charge") == 0) {
                    options.surface_color = SURFACE_COLOR_POLARITY;
                } else {
                    fprintf(stderr, "Error: Invalid surface color '%s'. Use uniform, atom, or polarity.\n", optarg);
                    return 1;
                }
                break;
            case 1013:  /* --modern - now default, kept for compatibility */
                /* Already default, no changes needed */
                break;
            case 1014:  /* --heteroatom-gap */
                options.heteroatom_gap = atof(optarg);
                if (options.heteroatom_gap < 0.0) options.heteroatom_gap = 0.0;
                if (options.heteroatom_gap > 1.0) options.heteroatom_gap = 1.0;
                break;
            case 1015:  /* --scale */
                options.scale_factor = atof(optarg);
                if (options.scale_factor < 0.1) options.scale_factor = 0.1;
                if (options.scale_factor > 10.0) options.scale_factor = 10.0;
                break;
            case 1016:  /* --line-cap */
                if (strcmp(optarg, "round") == 0) {
                    options.line_cap = LINE_CAP_ROUND;
                } else if (strcmp(optarg, "butt") == 0) {
                    options.line_cap = LINE_CAP_BUTT;
                } else if (strcmp(optarg, "square") == 0) {
                    options.line_cap = LINE_CAP_SQUARE;
                } else {
                    fprintf(stderr, "Error: Invalid line cap '%s'. Use round, butt, or square.\n", optarg);
                    return 1;
                }
                break;
            case 1017:  /* --terminal-carbons */
                options.terminal_carbon_labels = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_depict_usage(argv[0]);
                return 0;
            default:
                print_depict_usage(argv[0]);
                return 1;
        }
    }

    if (!smiles) {
        fprintf(stderr, "Error: SMILES string required (-S/--smiles)\n");
        print_depict_usage(argv[0]);
        return 1;
    }

    if (!output_file) {
        fprintf(stderr, "Error: Output file required (-o/--output)\n");
        print_depict_usage(argv[0]);
        return 1;
    }

    /* Auto-detect output format from file extension */
    const char* ext = strrchr(output_file, '.');
    if (ext) {
        if (strcasecmp(ext, ".svg") == 0) {
            options.format = IMG_FORMAT_SVG;
        } else if (strcasecmp(ext, ".png") == 0) {
            options.format = IMG_FORMAT_PNG;
        } else if (strcasecmp(ext, ".jpg") == 0 || strcasecmp(ext, ".jpeg") == 0) {
            options.format = IMG_FORMAT_JPEG;
        }
    }

    /* Get render style name */
    const char* style_name = "sticks";
    switch (options.render_style) {
        case RENDER_STYLE_WIREFRAME: style_name = "wireframe"; break;
        case RENDER_STYLE_STICKS: style_name = "sticks"; break;
        case RENDER_STYLE_BALLS_AND_STICKS: style_name = "balls-and-sticks"; break;
        case RENDER_STYLE_SPACEFILL: style_name = "spacefill"; break;
        case RENDER_STYLE_SURFACE: style_name = "surface"; break;
    }

    char error_buf[256];
    depict_info_t info = {0};
    cchem_status_t status;

    if (verbose) {
        status = depict_smiles_verbose(smiles, output_file, &options,
                                       &info, error_buf, sizeof(error_buf));
    } else {
        status = depict_smiles(smiles, output_file, &options,
                               error_buf, sizeof(error_buf));
    }

    if (status == CCHEM_OK) {
        if (verbose) {
            printf("=== Depiction Summary ===\n");
            printf("Input SMILES:     %s\n", smiles);
            printf("Canonical SMILES: %s\n", info.canonical_smiles);
            printf("Output file:      %s\n", output_file);
            printf("Mode:             %s\n", options.mode == DEPICT_MODE_2D ? "2D" : "3D (MMFF94)");
            printf("Render style:     %s\n", style_name);
            printf("Dimensions:       %dx%d\n", options.width, options.height);
            printf("Molecule:         %d atoms, %d bonds, %d rings\n",
                   info.num_atoms, info.num_bonds, info.num_rings);
            if (options.mode == DEPICT_MODE_3D) {
                printf("MMFF94 Energy:    %.2f -> %.2f kcal/mol (delta: %.2f)\n",
                       info.energy_initial, info.energy_final,
                       info.energy_initial - info.energy_final);
            }
            printf("=========================\n");
        }
        return 0;
    } else {
        fprintf(stderr, "Error: %s\n", error_buf);
        return 1;
    }
}

static void print_split_usage(const char* prog_name) {
    printf("Usage: %s split [options]\n\n", prog_name);
    printf("Split dataset into train/test/validation sets\n\n");
    printf("Options:\n");
    printf("  -f, --file <file.csv>         Input CSV file (required)\n");
    printf("  -o, --output <f1,f2,...>      Comma-separated output file names (required)\n");
    printf("  -s, --col <name>              SMILES column name (required)\n");
    printf("  --split-ratios <ratios>       Split ratios, e.g., \"80,10,10\" or \"0.8,0.1,0.1\"\n");
    printf("                                (default: 80,20 for 2 outputs, 80,10,10 for 3)\n");
    printf("  --splitting-method <method>   Splitting method: 'random' or 'scaffold'\n");
    printf("                                (default: random)\n");
    printf("  --stratified                  For scaffold split: distribute each scaffold\n");
    printf("                                proportionally across all outputs (default: off)\n");
    printf("  -n, --ncpu <N>                Number of CPU cores (default: all)\n");
    printf("  --seed <N>                    Random seed for reproducibility\n");
    printf("  -v, --verbose                 Verbose output\n");
    printf("  -h, --help                    Print this help message\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,val.csv,test.csv --split-ratios 80,10,10\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv --splitting-method scaffold\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv --splitting-method scaffold --stratified\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv --split-ratios 8,2 --seed 42\n", prog_name);
}

/* Split command */
static int cmd_split(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"file",             required_argument, 0, 'f'},
        {"output",           required_argument, 0, 'o'},
        {"col",              required_argument, 0, 's'},
        {"split-ratios",     required_argument, 0, 1000},
        {"splitting-method", required_argument, 0, 1001},
        {"stratified",       no_argument,       0, 1003},
        {"ncpu",             required_argument, 0, 'n'},
        {"seed",             required_argument, 0, 1002},
        {"verbose",          no_argument,       0, 'v'},
        {"help",             no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* input_file = NULL;
    char* output_str = NULL;
    char* smiles_col = NULL;
    char* ratio_str = NULL;
    char* method_str = "random";
    int ncpu = 0;
    unsigned int seed = 0;
    bool verbose = false;
    bool stratified = false;
    (void)argc;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "f:o:s:n:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'f':
                input_file = optarg;
                break;
            case 'o':
                output_str = optarg;
                break;
            case 's':
                smiles_col = optarg;
                break;
            case 1000:  /* --split-ratios */
                ratio_str = optarg;
                break;
            case 1001:  /* --splitting-method */
                method_str = optarg;
                break;
            case 'n':
                ncpu = atoi(optarg);
                if (ncpu < 0) ncpu = 0;
                break;
            case 1002:  /* --seed */
                seed = (unsigned int)atoi(optarg);
                break;
            case 1003:  /* --stratified */
                stratified = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_split_usage(argv[0]);
                return 0;
            default:
                print_split_usage(argv[0]);
                return 1;
        }
    }

    /* Validate required arguments */
    if (!input_file) {
        fprintf(stderr, "Error: Input file required (-f/--file)\n");
        return 1;
    }
    if (!output_str) {
        fprintf(stderr, "Error: Output files required (-o/--output)\n");
        return 1;
    }
    if (!smiles_col) {
        fprintf(stderr, "Error: SMILES column name required (-s/--col)\n");
        return 1;
    }

    /* Parse output files */
    char* output_files_buf[SPLITTER_MAX_SPLITS];
    int num_outputs = 0;

    char* output_copy = strdup(output_str);
    if (!output_copy) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }

    char* token = strtok(output_copy, ",");
    while (token && num_outputs < SPLITTER_MAX_SPLITS) {
        output_files_buf[num_outputs++] = strdup(token);
        token = strtok(NULL, ",");
    }
    free(output_copy);

    if (num_outputs < 2) {
        fprintf(stderr, "Error: At least 2 output files required\n");
        for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
        return 1;
    }

    /* Parse splitting method */
    split_method_t method = SPLIT_METHOD_RANDOM;
    if (strcmp(method_str, "scaffold") == 0) {
        method = SPLIT_METHOD_SCAFFOLD;
    } else if (strcmp(method_str, "random") != 0) {
        fprintf(stderr, "Error: Unknown splitting method '%s'. Use 'random' or 'scaffold'.\n", method_str);
        for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
        return 1;
    }

    /* Set up split options */
    split_options_t options = SPLIT_OPTIONS_DEFAULT;
    options.method = method;
    options.num_splits = num_outputs;
    options.num_threads = ncpu;
    options.seed = seed;
    options.verbose = verbose;
    options.stratified = stratified;

    /* Parse ratios */
    if (ratio_str) {
        int num_ratios;
        if (splitter_parse_ratios(ratio_str, options.ratios, &num_ratios) != CCHEM_OK) {
            fprintf(stderr, "Error: Invalid split ratios '%s'\n", ratio_str);
            for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
            return 1;
        }
        if (num_ratios != num_outputs) {
            fprintf(stderr, "Error: Number of ratios (%d) must match number of output files (%d)\n",
                    num_ratios, num_outputs);
            for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
            return 1;
        }
    } else {
        /* Default ratios based on number of outputs */
        if (num_outputs == 2) {
            options.ratios[0] = 0.8;
            options.ratios[1] = 0.2;
        } else if (num_outputs == 3) {
            options.ratios[0] = 0.8;
            options.ratios[1] = 0.1;
            options.ratios[2] = 0.1;
        } else {
            /* Equal splits for more than 3 */
            for (int i = 0; i < num_outputs; i++) {
                options.ratios[i] = 1.0 / num_outputs;
            }
        }
    }

    /* Run the split */
    char error_buf[512];
    const char* output_files[SPLITTER_MAX_SPLITS];
    for (int i = 0; i < num_outputs; i++) {
        output_files[i] = output_files_buf[i];
    }

    cchem_status_t status = splitter_split_csv(
        input_file, output_files, num_outputs, smiles_col,
        &options, error_buf, sizeof(error_buf));

    /* Cleanup */
    for (int i = 0; i < num_outputs; i++) {
        free(output_files_buf[i]);
    }

    if (status != CCHEM_OK) {
        fprintf(stderr, "Error: %s\n", error_buf);
        return 1;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char* command = argv[1];

    /* Shift arguments for subcommand parsing */
    argc--;
    argv++;

    if (strcmp(command, "canonicalize") == 0) {
        return cmd_canonicalize(argc, argv);
    } else if (strcmp(command, "compute") == 0) {
        return cmd_compute(argc, argv);
    } else if (strcmp(command, "depict") == 0) {
        return cmd_depict(argc, argv);
    } else if (strcmp(command, "split") == 0) {
        return cmd_split(argc, argv);
    } else if (strcmp(command, "validate") == 0) {
        return cmd_validate(argc, argv);
    } else if (strcmp(command, "version") == 0 || strcmp(command, "--version") == 0) {
        print_version();
        return 0;
    } else if (strcmp(command, "help") == 0 || strcmp(command, "--help") == 0 || strcmp(command, "-h") == 0) {
        print_usage(argv[0]);
        return 0;
    } else {
        fprintf(stderr, "Unknown command: %s\n", command);
        print_usage(argv[0]);
        return 1;
    }
}

/* Library initialization */
const char* cchem_version(void) {
    return CCHEM_VERSION_STRING;
}

cchem_status_t cchem_init(void) {
    /* Nothing to initialize currently */
    return CCHEM_OK;
}

void cchem_cleanup(void) {
    /* Nothing to cleanup currently */
}

char* cchem_canonicalize(const char* smiles, char* error_buf, size_t error_buf_size) {
    return smiles_canonicalize(smiles, NULL, error_buf, error_buf_size);
}

bool cchem_smiles_equal(const char* smiles1, const char* smiles2) {
    return smiles_are_equivalent(smiles1, smiles2);
}

cchem_status_t cchem_validate_smiles(const char* smiles, char* error_buf, size_t error_buf_size) {
    return smiles_validate(smiles, error_buf, error_buf_size);
}

cchem_status_t cchem_process_csv(const char* input_file, const char* output_file,
                                 const char* smiles_column, const char* output_column,
                                 int num_threads) {
    csv_batch_context_t* ctx = csv_batch_context_create(
        input_file, output_file, smiles_column, output_column, num_threads);

    if (!ctx) return CCHEM_ERROR_MEMORY;

    cchem_status_t status = csv_batch_canonicalize(ctx);
    csv_batch_context_free(ctx);

    return status;
}
