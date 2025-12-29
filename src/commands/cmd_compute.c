/**
 * @file cmd_compute.c
 * @brief Compute command implementation with pipeline streaming
 */

#include "cchem/utils/commands.h"

void print_compute_usage(const char* prog_name) {
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

int cmd_compute(int argc, char* argv[]) {
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
