/**
 * @file cmd_canonicalize.c
 * @brief Canonicalize command implementation
 */

#include "cchem/utils/commands.h"

/* ============================================================================
 * Parquet Batch Canonicalization
 * ============================================================================ */

#ifdef HAVE_PARQUET

/* Work item for parquet processing */
typedef struct {
    int row_idx;
    char* smiles;
    char** original_fields;
    int num_fields;
} pq_work_item_t;

/* Result item for parquet processing */
typedef struct {
    int row_idx;
    char** output_fields;
    int num_fields;
    bool success;
} pq_result_item_t;

/* Parquet canonicalization pipeline context */
typedef struct {
    parquet_stream_reader_t* reader;
    int smiles_col_idx;
    int num_input_fields;

    parquet_writer_t* writer;
    char** output_col_names;
    int num_output_cols;

    bounded_queue_t* input_queue;
    bounded_queue_t* output_queue;

    pthread_t reader_thread;
    pthread_t writer_thread;
    pthread_t* worker_threads;
    int num_workers;

    atomic_size_t rows_read;
    atomic_size_t rows_processed;
    atomic_size_t rows_written;
    size_t total_rows;

    const char* output_col_name;
    int output_col_idx;
} pq_pipeline_ctx_t;

static void pq_work_item_free(pq_work_item_t* item) {
    if (!item) return;
    free(item->smiles);
    if (item->original_fields) {
        for (int i = 0; i < item->num_fields; i++) {
            free(item->original_fields[i]);
        }
        free(item->original_fields);
    }
    free(item);
}

static void pq_result_item_free(pq_result_item_t* result) {
    if (!result) return;
    if (result->output_fields) {
        for (int i = 0; i < result->num_fields; i++) {
            free(result->output_fields[i]);
        }
        free(result->output_fields);
    }
    free(result);
}

/* Reader thread for parquet */
static void* pq_canon_reader_thread(void* arg) {
    pq_pipeline_ctx_t* ctx = (pq_pipeline_ctx_t*)arg;

    const char* fields[256];
    int row_idx = 0;

    int nfields;
    while ((nfields = parquet_stream_reader_next(ctx->reader, fields, 256)) > 0) {
        pq_work_item_t* item = (pq_work_item_t*)malloc(sizeof(pq_work_item_t));
        if (!item) continue;

        item->row_idx = row_idx++;
        item->num_fields = nfields;

        /* Copy SMILES */
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

        if (!queue_push(ctx->input_queue, item)) {
            pq_work_item_free(item);
            break;
        }
        atomic_fetch_add(&ctx->rows_read, 1);
    }

    queue_close(ctx->input_queue);
    return NULL;
}

/* Worker thread for parquet canonicalization */
static void* pq_canon_worker_thread(void* arg) {
    pq_pipeline_ctx_t* ctx = (pq_pipeline_ctx_t*)arg;

    char error_buf[256];
    pq_work_item_t* item;

    while ((item = (pq_work_item_t*)queue_pop(ctx->input_queue)) != NULL) {
        pq_result_item_t* result = (pq_result_item_t*)malloc(sizeof(pq_result_item_t));
        if (!result) {
            pq_work_item_free(item);
            continue;
        }

        result->row_idx = item->row_idx;
        result->num_fields = ctx->num_output_cols;
        result->output_fields = (char**)calloc(ctx->num_output_cols, sizeof(char*));
        result->success = false;

        if (!result->output_fields) {
            pq_result_item_free(result);
            pq_work_item_free(item);
            continue;
        }

        /* Copy original fields */
        for (int i = 0; i < item->num_fields && i < ctx->num_output_cols - 1; i++) {
            result->output_fields[i] = item->original_fields[i] ?
                                       strdup(item->original_fields[i]) : strdup("");
        }

        /* Canonicalize SMILES */
        char* canonical = NULL;
        if (item->smiles && item->smiles[0] != '\0') {
            canonical = smiles_canonicalize(item->smiles, NULL, error_buf, sizeof(error_buf));
        }

        /* Set output column (last column is the canonical SMILES) */
        result->output_fields[ctx->output_col_idx] = canonical ? canonical : strdup("");
        result->success = (canonical != NULL);

        atomic_fetch_add(&ctx->rows_processed, 1);

        if (!queue_push(ctx->output_queue, result)) {
            pq_result_item_free(result);
        }

        pq_work_item_free(item);
    }

    return NULL;
}

/* Writer thread for parquet */
static void* pq_canon_writer_thread(void* arg) {
    pq_pipeline_ctx_t* ctx = (pq_pipeline_ctx_t*)arg;

    pq_result_item_t* result;

    while ((result = (pq_result_item_t*)queue_pop(ctx->output_queue)) != NULL) {
        parquet_writer_write_row(ctx->writer,
                                  (const char**)result->output_fields,
                                  result->num_fields);
        atomic_fetch_add(&ctx->rows_written, 1);
        pq_result_item_free(result);
    }

    return NULL;
}

static int cmd_canonicalize_parquet(const char* input_file, const char* output_file,
                                     const char* input_col, const char* output_col,
                                     int ncpu, bool verbose) {
    pq_pipeline_ctx_t ctx = {0};

    /* Count rows */
    ctx.total_rows = parquet_count_rows(input_file);
    if (ctx.total_rows == 0) {
        fprintf(stderr, "Error: Input file is empty or invalid\n");
        return 1;
    }

    /* Create parquet reader */
    ctx.reader = parquet_stream_reader_create(input_file);
    if (!ctx.reader) {
        fprintf(stderr, "Error: Failed to open parquet file: %s\n", input_file);
        return 1;
    }

    ctx.num_input_fields = parquet_stream_reader_num_columns(ctx.reader);

    /* Find SMILES column */
    ctx.smiles_col_idx = parquet_stream_reader_find_column(ctx.reader, input_col);
    if (ctx.smiles_col_idx < 0) {
        fprintf(stderr, "Error: Column '%s' not found in parquet file\n", input_col);
        parquet_stream_reader_free(ctx.reader);
        return 1;
    }

    /* Build output column names (original columns + canonical column) */
    ctx.num_output_cols = ctx.num_input_fields + 1;
    ctx.output_col_names = (char**)calloc(ctx.num_output_cols, sizeof(char*));
    if (!ctx.output_col_names) {
        parquet_stream_reader_free(ctx.reader);
        return 1;
    }

    for (int i = 0; i < ctx.num_input_fields; i++) {
        ctx.output_col_names[i] = strdup(parquet_stream_reader_column_name(ctx.reader, i));
    }
    ctx.output_col_names[ctx.num_input_fields] = strdup(output_col);
    ctx.output_col_idx = ctx.num_input_fields;
    ctx.output_col_name = output_col;

    /* Determine output format (parquet or CSV based on extension) */
    bool output_parquet = is_parquet_file(output_file);

    if (output_parquet) {
        /* All columns are strings for canonicalize (pass NULL for types = all STRING) */
        ctx.writer = parquet_writer_create(output_file,
                                            (const char**)ctx.output_col_names,
                                            NULL,
                                            ctx.num_output_cols);
        if (!ctx.writer) {
            fprintf(stderr, "Error: Failed to create output parquet file\n");
            for (int i = 0; i < ctx.num_output_cols; i++) free(ctx.output_col_names[i]);
            free(ctx.output_col_names);
            parquet_stream_reader_free(ctx.reader);
            return 1;
        }
    } else {
        /* Fall back to CSV output */
        fprintf(stderr, "Error: Output must be a parquet file when input is parquet\n");
        for (int i = 0; i < ctx.num_output_cols; i++) free(ctx.output_col_names[i]);
        free(ctx.output_col_names);
        parquet_stream_reader_free(ctx.reader);
        return 1;
    }

    /* Set up thread count */
    ctx.num_workers = (ncpu > 0) ? ncpu : parallel_get_num_cores();

    if (verbose) {
        printf("Parquet mode: %zu molecules, %d workers\n", ctx.total_rows, ctx.num_workers);
    }

    /* Create queues */
    int queue_size = ctx.num_workers * 2;
    ctx.input_queue = queue_create(queue_size);
    ctx.output_queue = queue_create(queue_size);
    if (!ctx.input_queue || !ctx.output_queue) {
        fprintf(stderr, "Error: Failed to create queues\n");
        if (ctx.input_queue) queue_free(ctx.input_queue);
        if (ctx.output_queue) queue_free(ctx.output_queue);
        parquet_writer_free(ctx.writer);
        for (int i = 0; i < ctx.num_output_cols; i++) free(ctx.output_col_names[i]);
        free(ctx.output_col_names);
        parquet_stream_reader_free(ctx.reader);
        return 1;
    }

    /* Initialize atomics */
    atomic_store(&ctx.rows_read, 0);
    atomic_store(&ctx.rows_processed, 0);
    atomic_store(&ctx.rows_written, 0);

    printf("Canonicalizing %zu molecules...\n", ctx.total_rows);

    /* Create progress bar */
    progress_config_t prog_config = PROGRESS_CONFIG_DEFAULT;
    prog_config.prefix = "Processing";
    progress_t* progress = progress_create(ctx.total_rows, &prog_config);

    /* Start threads */
    pthread_create(&ctx.reader_thread, NULL, pq_canon_reader_thread, &ctx);

    ctx.worker_threads = (pthread_t*)malloc(ctx.num_workers * sizeof(pthread_t));
    for (int i = 0; i < ctx.num_workers; i++) {
        pthread_create(&ctx.worker_threads[i], NULL, pq_canon_worker_thread, &ctx);
    }

    pthread_create(&ctx.writer_thread, NULL, pq_canon_writer_thread, &ctx);

    /* Monitor progress */
    while (atomic_load(&ctx.rows_written) < ctx.total_rows) {
        size_t written = atomic_load(&ctx.rows_written);
        progress_update(progress, written);
        usleep(50000);
    }

    /* Wait for completion */
    pthread_join(ctx.reader_thread, NULL);
    for (int i = 0; i < ctx.num_workers; i++) {
        pthread_join(ctx.worker_threads[i], NULL);
    }
    queue_close(ctx.output_queue);
    pthread_join(ctx.writer_thread, NULL);

    progress_finish(progress);
    progress_free(progress);

    /* Cleanup */
    queue_free(ctx.input_queue);
    queue_free(ctx.output_queue);
    free(ctx.worker_threads);
    parquet_writer_free(ctx.writer);
    for (int i = 0; i < ctx.num_output_cols; i++) free(ctx.output_col_names[i]);
    free(ctx.output_col_names);
    parquet_stream_reader_free(ctx.reader);

    size_t rows_written = atomic_load(&ctx.rows_written);
    fprintf(stderr, "\033[KProcessed: %zu/%zu (%.1f%% success)\n",
           rows_written, ctx.total_rows,
           ctx.total_rows > 0 ? (double)rows_written / ctx.total_rows * 100.0 : 0.0);
    fprintf(stderr, "Output written to: %s\n", output_file);

    return 0;
}

/* ============================================================================
 * CSV to Parquet Canonicalization
 * ============================================================================ */

/* CSV to Parquet pipeline context */
typedef struct {
    csv_mmap_reader_t* reader;
    int smiles_col_idx;
    int num_input_fields;
    char** header_fields;

    parquet_writer_t* writer;
    char** output_col_names;
    int num_output_cols;

    bounded_queue_t* input_queue;
    bounded_queue_t* output_queue;

    pthread_t reader_thread;
    pthread_t writer_thread;
    pthread_t* worker_threads;
    int num_workers;

    atomic_size_t rows_read;
    atomic_size_t rows_processed;
    atomic_size_t rows_written;
    size_t total_rows;

    const char* output_col_name;
    int output_col_idx;
    char delimiter;
} csv_to_pq_canon_ctx_t;

/* CSV reader thread for csv-to-parquet canonicalization */
static void* csv_to_pq_canon_reader_thread(void* arg) {
    csv_to_pq_canon_ctx_t* ctx = (csv_to_pq_canon_ctx_t*)arg;

    size_t line_len;
    const char* line;
    int row_idx = 0;
    const char* fields[256];

    while ((line = csv_mmap_reader_next_line(ctx->reader, &line_len))) {
        int nfields = csv_mmap_parse_line(line, line_len, ctx->delimiter, fields, 256);

        pq_work_item_t* item = (pq_work_item_t*)malloc(sizeof(pq_work_item_t));
        if (!item) continue;

        item->row_idx = row_idx++;
        item->num_fields = nfields;

        /* Copy SMILES */
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

        if (!queue_push(ctx->input_queue, item)) {
            pq_work_item_free(item);
            break;
        }
        atomic_fetch_add(&ctx->rows_read, 1);
    }

    queue_close(ctx->input_queue);
    return NULL;
}

/* Worker thread for csv-to-parquet canonicalization */
static void* csv_to_pq_canon_worker_thread(void* arg) {
    csv_to_pq_canon_ctx_t* ctx = (csv_to_pq_canon_ctx_t*)arg;

    char error_buf[256];
    pq_work_item_t* item;

    while ((item = (pq_work_item_t*)queue_pop(ctx->input_queue)) != NULL) {
        pq_result_item_t* result = (pq_result_item_t*)malloc(sizeof(pq_result_item_t));
        if (!result) {
            pq_work_item_free(item);
            continue;
        }

        result->row_idx = item->row_idx;
        result->num_fields = ctx->num_output_cols;
        result->output_fields = (char**)calloc(ctx->num_output_cols, sizeof(char*));
        result->success = false;

        if (!result->output_fields) {
            pq_result_item_free(result);
            pq_work_item_free(item);
            continue;
        }

        /* Copy original fields */
        for (int i = 0; i < item->num_fields && i < ctx->num_output_cols - 1; i++) {
            result->output_fields[i] = item->original_fields[i] ?
                                       strdup(item->original_fields[i]) : strdup("");
        }

        /* Canonicalize SMILES */
        char* canonical = NULL;
        if (item->smiles && item->smiles[0] != '\0') {
            canonical = smiles_canonicalize(item->smiles, NULL, error_buf, sizeof(error_buf));
        }

        /* Set output column (last column is the canonical SMILES) */
        result->output_fields[ctx->output_col_idx] = canonical ? canonical : strdup("");
        result->success = (canonical != NULL);

        atomic_fetch_add(&ctx->rows_processed, 1);

        if (!queue_push(ctx->output_queue, result)) {
            pq_result_item_free(result);
        }

        pq_work_item_free(item);
    }

    return NULL;
}

/* Writer thread for csv-to-parquet canonicalization */
static void* csv_to_pq_canon_writer_thread(void* arg) {
    csv_to_pq_canon_ctx_t* ctx = (csv_to_pq_canon_ctx_t*)arg;

    pq_result_item_t* result;

    while ((result = (pq_result_item_t*)queue_pop(ctx->output_queue)) != NULL) {
        parquet_writer_write_row(ctx->writer,
                                  (const char**)result->output_fields,
                                  result->num_fields);
        atomic_fetch_add(&ctx->rows_written, 1);
        pq_result_item_free(result);
    }

    return NULL;
}

static int cmd_canonicalize_csv_to_parquet(const char* input_file, const char* output_file,
                                            const char* input_col, const char* output_col,
                                            int ncpu, bool verbose) {
    csv_to_pq_canon_ctx_t ctx = {0};
    ctx.delimiter = ',';

    /* Count rows */
    ctx.total_rows = csv_count_rows(input_file);
    if (ctx.total_rows <= 1) {
        fprintf(stderr, "Error: Input file is empty or has only header\n");
        return 1;
    }
    ctx.total_rows--;  /* Exclude header */

    /* Create CSV reader */
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

    /* Build output column names (original columns + canonical column) */
    ctx.num_output_cols = ctx.num_input_fields + 1;
    ctx.output_col_names = (char**)calloc(ctx.num_output_cols, sizeof(char*));
    if (!ctx.output_col_names) {
        for (int i = 0; i < ctx.num_input_fields; i++) free(ctx.header_fields[i]);
        free(ctx.header_fields);
        csv_mmap_reader_free(ctx.reader);
        return 1;
    }

    for (int i = 0; i < ctx.num_input_fields; i++) {
        ctx.output_col_names[i] = strdup(ctx.header_fields[i]);
    }
    ctx.output_col_names[ctx.num_input_fields] = strdup(output_col);
    ctx.output_col_idx = ctx.num_input_fields;
    ctx.output_col_name = output_col;

    /* Create parquet writer (all columns are strings, pass NULL for types) */
    ctx.writer = parquet_writer_create(output_file,
                                        (const char**)ctx.output_col_names,
                                        NULL,
                                        ctx.num_output_cols);
    if (!ctx.writer) {
        fprintf(stderr, "Error: Failed to create output parquet file: %s\n", output_file);
        for (int i = 0; i < ctx.num_output_cols; i++) free(ctx.output_col_names[i]);
        free(ctx.output_col_names);
        for (int i = 0; i < ctx.num_input_fields; i++) free(ctx.header_fields[i]);
        free(ctx.header_fields);
        csv_mmap_reader_free(ctx.reader);
        return 1;
    }

    /* Set up thread count */
    ctx.num_workers = (ncpu > 0) ? ncpu : parallel_get_num_cores();

    if (verbose) {
        printf("CSV to Parquet mode: %zu molecules, %d workers\n", ctx.total_rows, ctx.num_workers);
    }

    /* Create queues */
    int queue_size = ctx.num_workers * 2;
    ctx.input_queue = queue_create(queue_size);
    ctx.output_queue = queue_create(queue_size);
    if (!ctx.input_queue || !ctx.output_queue) {
        fprintf(stderr, "Error: Failed to create queues\n");
        if (ctx.input_queue) queue_free(ctx.input_queue);
        if (ctx.output_queue) queue_free(ctx.output_queue);
        parquet_writer_free(ctx.writer);
        for (int i = 0; i < ctx.num_output_cols; i++) free(ctx.output_col_names[i]);
        free(ctx.output_col_names);
        for (int i = 0; i < ctx.num_input_fields; i++) free(ctx.header_fields[i]);
        free(ctx.header_fields);
        csv_mmap_reader_free(ctx.reader);
        return 1;
    }

    /* Initialize atomics */
    atomic_store(&ctx.rows_read, 0);
    atomic_store(&ctx.rows_processed, 0);
    atomic_store(&ctx.rows_written, 0);

    printf("Canonicalizing %zu molecules (CSV -> Parquet)...\n", ctx.total_rows);

    /* Create progress bar */
    progress_config_t prog_config = PROGRESS_CONFIG_DEFAULT;
    prog_config.prefix = "Processing";
    progress_t* progress = progress_create(ctx.total_rows, &prog_config);

    /* Start threads */
    pthread_create(&ctx.reader_thread, NULL, csv_to_pq_canon_reader_thread, &ctx);

    ctx.worker_threads = (pthread_t*)malloc(ctx.num_workers * sizeof(pthread_t));
    for (int i = 0; i < ctx.num_workers; i++) {
        pthread_create(&ctx.worker_threads[i], NULL, csv_to_pq_canon_worker_thread, &ctx);
    }

    pthread_create(&ctx.writer_thread, NULL, csv_to_pq_canon_writer_thread, &ctx);

    /* Monitor progress */
    while (atomic_load(&ctx.rows_written) < ctx.total_rows) {
        size_t written = atomic_load(&ctx.rows_written);
        progress_update(progress, written);
        usleep(50000);
    }

    /* Wait for completion */
    pthread_join(ctx.reader_thread, NULL);
    for (int i = 0; i < ctx.num_workers; i++) {
        pthread_join(ctx.worker_threads[i], NULL);
    }
    queue_close(ctx.output_queue);
    pthread_join(ctx.writer_thread, NULL);

    progress_finish(progress);
    progress_free(progress);

    /* Cleanup */
    queue_free(ctx.input_queue);
    queue_free(ctx.output_queue);
    free(ctx.worker_threads);
    parquet_writer_free(ctx.writer);
    for (int i = 0; i < ctx.num_output_cols; i++) free(ctx.output_col_names[i]);
    free(ctx.output_col_names);
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

#endif /* HAVE_PARQUET */

/* ============================================================================
 * Canonicalize Command
 * ============================================================================ */

void print_canonicalize_usage(const char* prog_name) {
    printf("Usage: %s canonicalize [options]\n\n", prog_name);
    printf("Canonicalize SMILES strings (aromatized by default)\n\n");
    printf("Options:\n");
    printf("  -S, --smiles <SMILES>   Single SMILES string to canonicalize\n");
    printf("  -f, --file <file>       Input file (CSV or Parquet)\n");
    printf("  -s, --col <name>        Input SMILES column name (required with -f)\n");
    printf("  -c, --outcol <name>     Output column name for canonical SMILES (default: canonical_smiles)\n");
    printf("  -o, --output <file>     Output file (CSV or Parquet, required with -f)\n");
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
    printf("Supported file formats:\n");
    printf("  .csv                    Comma-separated values\n");
    printf("  .parquet, .pq           Apache Parquet (requires matching output format)\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s canonicalize -S \"c1ccccc1\"\n", prog_name);
    printf("  %s canonicalize -S \"[Na+].CC(=O)[O-]\" --sanitize unsalt\n", prog_name);
    printf("  %s canonicalize -S \"C1=CC=CC=C1\" --sanitize kekulize\n", prog_name);
    printf("  %s canonicalize -S \"CC(=O)NC\" --list-tautomers\n", prog_name);
    printf("  %s canonicalize -f input.csv -s smiles -c canonical -o output.csv -n 4\n", prog_name);
    printf("  %s canonicalize -f input.parquet -s smiles -o output.parquet\n", prog_name);
}

int cmd_canonicalize(int argc, char* argv[]) {
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

    /* Batch mode (CSV or Parquet) */
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

#ifdef HAVE_PARQUET
        bool input_is_parquet = is_parquet_file(input_file);
        bool output_is_parquet = is_parquet_file(output_file);

        if (input_is_parquet && output_is_parquet) {
            /* Parquet -> Parquet */
            return cmd_canonicalize_parquet(input_file, output_file, input_col,
                                             output_col, ncpu, verbose);
        } else if (!input_is_parquet && output_is_parquet) {
            /* CSV -> Parquet */
            return cmd_canonicalize_csv_to_parquet(input_file, output_file, input_col,
                                                    output_col, ncpu, verbose);
        } else if (input_is_parquet && !output_is_parquet) {
            fprintf(stderr, "Error: Cannot write CSV output from Parquet input. Use .parquet extension for output.\n");
            return 1;
        }
#endif

        /* CSV batch processing */
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
