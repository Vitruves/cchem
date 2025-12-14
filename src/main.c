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
#include "cchem/parallel.h"
#include "cchem/depictor/depictor.h"
#include "cchem/arena.h"
#include "cchem/splitter.h"

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

/* Descriptor batch processing types */
typedef struct {
    int row_idx;
    const char* smiles;  /* Pointer to source data (not owned) */
    char** desc_names;
    int num_descs;
    bool skip_canonicalization;
} desc_task_arg_t;

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

typedef struct {
    int row_idx;
    char* buffer;         /* Single pre-allocated buffer: num_values * DESC_VALUE_WIDTH */
    int num_values;
    bool success;
} desc_task_result_t;

/* Task argument for parallel CSV row formatting */
typedef struct {
    int row_idx;
    csv_bulk_writer_t* bulk;       /* Shared bulk writer (thread-safe per-row writes) */
    const char** fields;           /* Array of field pointers for this row */
    int num_fields;
} csv_format_task_t;

/* Get pointer to i-th value in result buffer */
static inline char* result_value_ptr(desc_task_result_t* r, int i) {
    return r->buffer + (i * DESC_VALUE_WIDTH);
}

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

/* Worker function for descriptor computation */
static void* desc_batch_worker(void* arg) {
    desc_task_arg_t* task = (desc_task_arg_t*)arg;

    /* Thread-local molecule pool to avoid malloc/free per molecule */
    static __thread molecule_t* tl_mol = NULL;

    /* Use thread-local cache for descriptor defs to avoid repeated lookups */
    static __thread const descriptor_def_t* cached_defs[MAX_DESCRIPTORS];
    static __thread int cached_num_defs = 0;
    static __thread int cached_total = 0;
    static __thread bool tl_initialized = false;

    if (!tl_initialized) {
        memset(cached_defs, 0, sizeof(cached_defs));
        cached_total = descriptor_count();
        cached_num_defs = descriptor_get_all(cached_defs, MAX_DESCRIPTORS);
        tl_initialized = true;
    }

    /* Initialize thread-local molecule on first use */
    if (!tl_mol) {
        tl_mol = molecule_create_with_capacity(256, 256);  /* Pre-size for typical molecules */
    }

    desc_task_result_t* result = (desc_task_result_t*)malloc(sizeof(desc_task_result_t));
    if (!result) return NULL;

    result->row_idx = task->row_idx;
    result->num_values = task->num_descs;

    /* Single allocation for all values */
    result->buffer = (char*)malloc(task->num_descs * DESC_VALUE_WIDTH);
    if (!result->buffer) {
        free(result);
        return NULL;
    }
    memset(result->buffer, 0, task->num_descs * DESC_VALUE_WIDTH);

    if (!task->smiles || task->smiles[0] == '\0') {
        result->success = false;
        return result;
    }

    char error_buf[256];
    molecule_t* mol = NULL;
    cchem_status_t status;

    if (task->skip_canonicalization) {
        /* Parse directly without canonicalization - reuse thread-local molecule */
        status = smiles_to_molecule_reuse(tl_mol, task->smiles, error_buf, sizeof(error_buf));
        mol = (status == CCHEM_OK) ? tl_mol : NULL;
    } else {
        /* Canonicalize SMILES first */
        char* canonical = smiles_canonicalize(task->smiles, NULL, error_buf, sizeof(error_buf));
        if (!canonical) {
            result->success = false;
            return result;
        }

        /* Parse molecule - reuse thread-local molecule */
        status = smiles_to_molecule_reuse(tl_mol, canonical, error_buf, sizeof(error_buf));
        mol = (status == CCHEM_OK) ? tl_mol : NULL;
        free(canonical);
    }

    if (!mol) {
        result->success = false;
        return result;
    }

    /* Compute descriptors - use batch for all, individual for subset */
    result->success = true;

    /* Check if computing all descriptors (use optimized batch path) */
    bool compute_all = (task->num_descs >= cached_total);

    if (compute_all) {
        /* Use cached defs (already loaded in thread-local storage) */
        /* Batch computation - zero-initialize to avoid valgrind warnings */
        descriptor_value_t all_values[MAX_DESCRIPTORS];
        memset(all_values, 0, sizeof(all_values));
        int num_computed = descriptors_compute_all(mol, NULL, all_values, MAX_DESCRIPTORS);

        /* Format with correct types - fast formatting for both ints and doubles */
        int n = (task->num_descs < num_computed) ? task->num_descs : num_computed;
        if (n > cached_num_defs) n = cached_num_defs;
        for (int i = 0; i < n; i++) {
            char* dest = result_value_ptr(result, i);
            if (cached_defs[i]->value_type == DESC_VALUE_INT) {
                fast_i64toa(all_values[i].i, dest, DESC_VALUE_WIDTH);
            } else {
                fast_dtoa(all_values[i].d, dest, DESC_VALUE_WIDTH);
            }
        }
    } else {
        /* Individual computation for subset of descriptors */
        for (int i = 0; i < task->num_descs; i++) {
            const descriptor_def_t* def = descriptor_get(task->desc_names[i]);
            char* dest = result_value_ptr(result, i);
            if (def) {
                descriptor_value_t value = {0};  /* Zero-initialize */
                if (def->compute(mol, &value) == CCHEM_OK) {
                    if (def->value_type == DESC_VALUE_INT) {
                        fast_i64toa(value.i, dest, DESC_VALUE_WIDTH);
                    } else {
                        fast_dtoa(value.d, dest, DESC_VALUE_WIDTH);
                    }
                }
            }
        }
    }

    /* Note: mol points to tl_mol, do NOT free - it's reused */
    return result;
}

/* Worker function for parallel CSV row formatting */
static void* csv_format_worker(void* arg) {
    csv_format_task_t* task = (csv_format_task_t*)arg;

    csv_bulk_format_row(task->bulk, task->row_idx, task->fields, task->num_fields);

    /* Return task pointer as result (for cleanup tracking) */
    return task;
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

            descriptor_value_t value = {0};  /* Zero-initialize */
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

    /* CSV batch mode */
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

        /* Determine number of threads */
        int actual_ncpu = (ncpu > 0) ? ncpu : parallel_get_num_cores();

        if (verbose) {
            printf("Processing: %s\n", input_file);
            printf("Input column: %s\n", input_col);
            printf("Descriptors: %s\n", desc_list);
            printf("Output file: %s\n", output_file);
            printf("CPU cores: %d\n", actual_ncpu);
            if (skip_canonicalization) {
                printf("Canonicalization: skipped\n");
            }
            printf("\n");
        }

        /* Read entire input CSV */
        csv_document_t* in_doc = csv_document_create();
        if (!in_doc || csv_document_read(in_doc, input_file, true) != CSV_OK) {
            fprintf(stderr, "Error: Failed to read input file: %s\n", input_file);
            if (in_doc) csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Find SMILES column */
        int smiles_col_idx = csv_find_column(&in_doc->header, input_col);
        if (smiles_col_idx < 0) {
            fprintf(stderr, "Error: Column '%s' not found in input file\n", input_col);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        size_t total_rows = (size_t)in_doc->num_rows;

        /* Create thread pool */
        thread_pool_t* pool = thread_pool_create(actual_ncpu);
        if (!pool) {
            fprintf(stderr, "Error: Failed to create thread pool\n");
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Pre-allocate task capacity */
        if (thread_pool_ensure_capacity(pool, (int)total_rows + 1) < 0) {
            fprintf(stderr, "Error: Failed to allocate task capacity\n");
            thread_pool_free(pool);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Show processing info */
        printf("Computing %d descriptors for %zu molecules...\n", num_descs, total_rows);

        /* Create progress bar */
        progress_config_t prog_config = PROGRESS_CONFIG_DEFAULT;
        prog_config.prefix = "Computing";
        progress_t* progress = progress_create(total_rows, &prog_config);

        /* Allocate task arguments */
        desc_task_arg_t* task_args = (desc_task_arg_t*)calloc(total_rows, sizeof(desc_task_arg_t));
        if (!task_args) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            progress_free(progress);
            thread_pool_free(pool);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Submit all tasks - avoid strdup, document stays alive during processing */
        for (int row_idx = 0; row_idx < in_doc->num_rows; row_idx++) {
            csv_row_t* row = &in_doc->rows[row_idx];
            const char* row_smiles = (smiles_col_idx < row->num_fields) ?
                                     row->fields[smiles_col_idx] : NULL;

            task_args[row_idx].row_idx = row_idx;
            task_args[row_idx].smiles = row_smiles ? row_smiles : "";  /* No strdup - doc stays alive */
            task_args[row_idx].desc_names = desc_names;
            task_args[row_idx].num_descs = num_descs;
            task_args[row_idx].skip_canonicalization = skip_canonicalization;

            thread_pool_submit(pool, desc_batch_worker, &task_args[row_idx]);
        }

        /* Wait for completion with progress updates */
        while (thread_pool_num_completed(pool) < (int)total_rows) {
            int completed = thread_pool_num_completed(pool);
            progress_update(progress, completed);
            usleep(50000);  /* 50ms */
        }

        progress_finish(progress);
        progress_free(progress);

        /* ================================================================
         * Phase 2: Parallel CSV formatting + bulk writing
         * ================================================================ */

        /* Prepare output header: original columns + descriptor columns */
        int out_num_fields = in_doc->header.num_fields + num_descs;
        const char** header_fields = (const char**)calloc(out_num_fields, sizeof(char*));

        for (int i = 0; i < in_doc->header.num_fields; i++) {
            header_fields[i] = in_doc->header.fields[i];
        }
        for (int i = 0; i < num_descs; i++) {
            header_fields[in_doc->header.num_fields + i] = desc_names[i];
        }

        /* Estimate average line size for bulk writer pre-allocation */
        size_t avg_line_size = (size_t)(out_num_fields * 12);  /* ~12 chars per field average */

        /* Create bulk writer for parallel row formatting */
        csv_bulk_writer_t* bulk = csv_bulk_writer_create((int)total_rows, avg_line_size);
        if (!bulk) {
            fprintf(stderr, "Error: Failed to create bulk writer\n");
            free(header_fields);
            free(task_args);
            thread_pool_free(pool);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Allocate field buffers for each row (for parallel access) */
        const char*** all_row_fields = (const char***)calloc(total_rows, sizeof(const char**));
        csv_format_task_t* format_tasks = (csv_format_task_t*)calloc(total_rows, sizeof(csv_format_task_t));
        if (!all_row_fields || !format_tasks) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            if (all_row_fields) free(all_row_fields);
            if (format_tasks) free(format_tasks);
            csv_bulk_writer_free(bulk);
            free(header_fields);
            free(task_args);
            thread_pool_free(pool);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Save descriptor results before reusing pool (need them for both value access and cleanup) */
        desc_task_result_t** desc_results = (desc_task_result_t**)calloc(total_rows, sizeof(desc_task_result_t*));
        if (!desc_results) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            free(all_row_fields);
            free(format_tasks);
            csv_bulk_writer_free(bulk);
            free(header_fields);
            free(task_args);
            thread_pool_free(pool);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }
        for (size_t i = 0; i < total_rows; i++) {
            desc_results[i] = (desc_task_result_t*)thread_pool_get_result(pool, (int)i);
        }

        /* Count successes and prepare row field arrays */
        size_t success_count = 0;
        for (int row_idx = 0; row_idx < in_doc->num_rows; row_idx++) {
            csv_row_t* row = &in_doc->rows[row_idx];
            const char** row_fields = (const char**)calloc(out_num_fields, sizeof(const char*));
            if (!row_fields) {
                /* Cleanup on failure */
                for (int j = 0; j < row_idx; j++) {
                    free((void*)all_row_fields[j]);
                }
                for (size_t j = 0; j < total_rows; j++) {
                    if (desc_results[j]) {
                        free(desc_results[j]->buffer);
                        free(desc_results[j]);
                    }
                }
                free(desc_results);
                free(all_row_fields);
                free(format_tasks);
                csv_bulk_writer_free(bulk);
                free(header_fields);
                free(task_args);
                thread_pool_free(pool);
                csv_document_free(in_doc);
                descriptor_free_names(desc_names, num_descs);
                return 1;
            }
            all_row_fields[row_idx] = row_fields;

            /* Copy original fields */
            for (int i = 0; i < in_doc->header.num_fields; i++) {
                row_fields[i] = (i < row->num_fields && row->fields[i]) ?
                                row->fields[i] : "";
            }

            /* Get computed result and fill descriptor values */
            desc_task_result_t* result = desc_results[row_idx];
            if (result && result->buffer) {
                for (int i = 0; i < num_descs; i++) {
                    const char* raw_value = (i < result->num_values) ? result_value_ptr(result, i) : "";
                    row_fields[in_doc->header.num_fields + i] = csv_sanitize_numeric(raw_value);
                }
                if (result->success) {
                    success_count++;
                }
            } else {
                for (int i = 0; i < num_descs; i++) {
                    row_fields[in_doc->header.num_fields + i] = "0";
                }
            }
        }

        /* Clear thread pool completed tasks for reuse */
        thread_pool_clear_completed(pool);

        /* Ensure capacity for format tasks */
        if (thread_pool_ensure_capacity(pool, (int)total_rows + 1) < 0) {
            fprintf(stderr, "Error: Failed to allocate format task capacity\n");
            for (size_t i = 0; i < total_rows; i++) free((void*)all_row_fields[i]);
            for (size_t i = 0; i < total_rows; i++) {
                if (desc_results[i]) { free(desc_results[i]->buffer); free(desc_results[i]); }
            }
            free(desc_results);
            free(all_row_fields);
            free(format_tasks);
            csv_bulk_writer_free(bulk);
            free(header_fields);
            free(task_args);
            thread_pool_free(pool);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Show CSV writing progress */
        printf("Writing CSV output (%d rows)...\n", (int)total_rows);

        /* Submit parallel CSV formatting tasks */
        for (int row_idx = 0; row_idx < (int)total_rows; row_idx++) {
            format_tasks[row_idx].row_idx = row_idx;
            format_tasks[row_idx].bulk = bulk;
            format_tasks[row_idx].fields = all_row_fields[row_idx];
            format_tasks[row_idx].num_fields = out_num_fields;
            thread_pool_submit(pool, csv_format_worker, &format_tasks[row_idx]);
        }

        /* Wait for all formatting to complete */
        thread_pool_wait_all(pool);

        /* Open output file and write header + bulk data */
        FILE* out_file = fopen(output_file, "w");
        if (!out_file) {
            fprintf(stderr, "Error: Failed to create output file: %s\n", output_file);
            for (size_t i = 0; i < total_rows; i++) free((void*)all_row_fields[i]);
            for (size_t i = 0; i < total_rows; i++) {
                if (desc_results[i]) { free(desc_results[i]->buffer); free(desc_results[i]); }
            }
            free(desc_results);
            free(all_row_fields);
            free(format_tasks);
            csv_bulk_writer_free(bulk);
            free(header_fields);
            free(task_args);
            thread_pool_free(pool);
            csv_document_free(in_doc);
            descriptor_free_names(desc_names, num_descs);
            return 1;
        }

        /* Use large buffer for output file */
        setvbuf(out_file, NULL, _IOFBF, 262144);  /* 256KB buffer */

        /* Write header */
        csv_write_header_line(out_file, header_fields, out_num_fields, ',');

        /* Bulk write all formatted rows */
        csv_status_t write_status = csv_bulk_write_all(bulk, out_file);
        fclose(out_file);

        if (write_status != CSV_OK) {
            fprintf(stderr, "Error: Failed to write output file\n");
        }

        /* Cleanup descriptor results (saved before pool reuse) */
        for (size_t i = 0; i < total_rows; i++) {
            if (desc_results[i]) {
                free(desc_results[i]->buffer);
                free(desc_results[i]);
            }
        }
        free(desc_results);

        /* Cleanup row field arrays */
        for (size_t i = 0; i < total_rows; i++) {
            free((void*)all_row_fields[i]);
        }
        free(all_row_fields);
        free(format_tasks);
        csv_bulk_writer_free(bulk);
        free(header_fields);
        free(task_args);
        thread_pool_free(pool);
        csv_document_free(in_doc);
        descriptor_free_names(desc_names, num_descs);

        printf("Processed: %zu/%zu (%.1f%% success)\n",
               success_count, total_rows,
               total_rows > 0 ? (double)success_count / total_rows * 100.0 : 0.0);
        printf("Output written to: %s\n", output_file);

        return 0;
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
