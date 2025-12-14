/**
 * @file parallel.c
 * @brief Parallel processing utilities using pthreads
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <pthread.h>
#include "cchem/parallel.h"
#include "cchem/csv.h"
#include "cchem/progress.h"
#include "cchem/canonicalizer/canon.h"
#include "cchem/canonicalizer/parser.h"

int parallel_get_num_cores(void) {
#ifdef _SC_NPROCESSORS_ONLN
    long cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (cores > 0) return (int)cores;
#endif
    return 4;  /* Default fallback */
}

/* Thread pool worker */
static void* thread_pool_worker(void* arg) {
    thread_pool_t* pool = (thread_pool_t*)arg;

    while (1) {
        pthread_mutex_lock(&pool->mutex);

        /* Wait for task or shutdown */
        while (pool->next_task >= pool->num_tasks && !pool->shutdown) {
            pthread_cond_wait(&pool->task_available, &pool->mutex);
        }

        if (pool->shutdown && pool->next_task >= pool->num_tasks) {
            pthread_mutex_unlock(&pool->mutex);
            break;
        }

        /* Get next task and copy data locally to avoid race with realloc */
        int task_id = pool->next_task++;
        parallel_task_fn fn = pool->tasks[task_id].fn;
        void* task_arg = pool->tasks[task_id].arg;
        pthread_mutex_unlock(&pool->mutex);

        if (task_id < pool->num_tasks) {
            /* Execute task with locally copied data */
            void* result = fn(task_arg);

            pthread_mutex_lock(&pool->mutex);
            pool->tasks[task_id].result = result;
            pool->tasks[task_id].completed = true;
            pool->completed_tasks++;
            pthread_cond_signal(&pool->task_completed);
            pthread_mutex_unlock(&pool->mutex);
        }
    }

    return NULL;
}

/* Pre-allocate task capacity to avoid reallocation during execution */
int thread_pool_ensure_capacity(thread_pool_t* pool, int capacity) {
    if (!pool || capacity <= 0) return -1;

    pthread_mutex_lock(&pool->mutex);

    if (capacity > pool->tasks_capacity) {
        void* new_tasks = realloc(pool->tasks, capacity * sizeof(*pool->tasks));
        if (!new_tasks) {
            pthread_mutex_unlock(&pool->mutex);
            return -1;
        }
        /* Zero the new portion */
        memset((char*)new_tasks + pool->tasks_capacity * sizeof(*pool->tasks),
               0, (capacity - pool->tasks_capacity) * sizeof(*pool->tasks));
        pool->tasks = new_tasks;
        pool->tasks_capacity = capacity;
    }

    pthread_mutex_unlock(&pool->mutex);
    return 0;
}

thread_pool_t* thread_pool_create(int num_threads) {
    if (num_threads <= 0) {
        num_threads = parallel_get_num_cores();
    }

    thread_pool_t* pool = (thread_pool_t*)calloc(1, sizeof(thread_pool_t));
    if (!pool) return NULL;

    pool->num_threads = num_threads;
    pool->threads = (pthread_t*)calloc(num_threads, sizeof(pthread_t));
    pool->thread_active = (bool*)calloc(num_threads, sizeof(bool));

    pool->tasks_capacity = 1024;
    pool->tasks = calloc(pool->tasks_capacity, sizeof(*pool->tasks));

    if (!pool->threads || !pool->thread_active || !pool->tasks) {
        thread_pool_free(pool);
        return NULL;
    }

    pthread_mutex_init(&pool->mutex, NULL);
    pthread_cond_init(&pool->task_available, NULL);
    pthread_cond_init(&pool->task_completed, NULL);

    pool->shutdown = false;
    pool->running = true;
    pool->num_tasks = 0;
    pool->next_task = 0;
    pool->completed_tasks = 0;

    /* Create worker threads */
    for (int i = 0; i < num_threads; i++) {
        if (pthread_create(&pool->threads[i], NULL, thread_pool_worker, pool) == 0) {
            pool->thread_active[i] = true;
        }
    }

    return pool;
}

void thread_pool_free(thread_pool_t* pool) {
    if (!pool) return;

    /* Signal shutdown */
    pthread_mutex_lock(&pool->mutex);
    pool->shutdown = true;
    pthread_cond_broadcast(&pool->task_available);
    pthread_mutex_unlock(&pool->mutex);

    /* Wait for threads */
    for (int i = 0; i < pool->num_threads; i++) {
        if (pool->thread_active[i]) {
            pthread_join(pool->threads[i], NULL);
        }
    }

    pthread_mutex_destroy(&pool->mutex);
    pthread_cond_destroy(&pool->task_available);
    pthread_cond_destroy(&pool->task_completed);

    if (pool->threads) free(pool->threads);
    if (pool->thread_active) free(pool->thread_active);
    if (pool->tasks) free(pool->tasks);

    free(pool);
}

int thread_pool_submit(thread_pool_t* pool, parallel_task_fn fn, void* arg) {
    if (!pool || !fn) return -1;

    pthread_mutex_lock(&pool->mutex);

    /* Ensure capacity */
    if (pool->num_tasks >= pool->tasks_capacity) {
        int new_cap = pool->tasks_capacity * 2;
        void* new_tasks = realloc(pool->tasks, new_cap * sizeof(*pool->tasks));
        if (!new_tasks) {
            pthread_mutex_unlock(&pool->mutex);
            return -1;
        }
        pool->tasks = new_tasks;
        pool->tasks_capacity = new_cap;
    }

    int task_id = pool->num_tasks++;
    pool->tasks[task_id].fn = fn;
    pool->tasks[task_id].arg = arg;
    pool->tasks[task_id].result = NULL;
    pool->tasks[task_id].completed = false;
    pool->tasks[task_id].has_error = false;

    pthread_cond_signal(&pool->task_available);
    pthread_mutex_unlock(&pool->mutex);

    return task_id;
}

void thread_pool_wait_all(thread_pool_t* pool) {
    if (!pool) return;

    pthread_mutex_lock(&pool->mutex);
    while (pool->completed_tasks < pool->num_tasks) {
        pthread_cond_wait(&pool->task_completed, &pool->mutex);
    }
    pthread_mutex_unlock(&pool->mutex);
}

void* thread_pool_get_result(thread_pool_t* pool, int task_id) {
    if (!pool || task_id < 0 || task_id >= pool->num_tasks) return NULL;
    return pool->tasks[task_id].result;
}

bool thread_pool_task_completed(thread_pool_t* pool, int task_id) {
    if (!pool || task_id < 0 || task_id >= pool->num_tasks) return false;
    return pool->tasks[task_id].completed;
}

int thread_pool_num_completed(thread_pool_t* pool) {
    if (!pool) return 0;

    pthread_mutex_lock(&pool->mutex);
    int completed = pool->completed_tasks;
    pthread_mutex_unlock(&pool->mutex);

    return completed;
}

void thread_pool_clear_completed(thread_pool_t* pool) {
    if (!pool) return;

    pthread_mutex_lock(&pool->mutex);

    /* Wait for all current tasks to complete first */
    while (pool->completed_tasks < pool->num_tasks) {
        pthread_cond_wait(&pool->task_completed, &pool->mutex);
    }

    /* Reset task queue counters for reuse (keeps allocated memory) */
    pool->num_tasks = 0;
    pool->next_task = 0;
    pool->completed_tasks = 0;

    /* Zero out the task array to clear old results */
    if (pool->tasks && pool->tasks_capacity > 0) {
        memset(pool->tasks, 0, pool->tasks_capacity * sizeof(*pool->tasks));
    }

    pthread_mutex_unlock(&pool->mutex);
}

/* Batch processing types */
typedef struct {
    int row_idx;
    char* smiles;
    csv_batch_context_t* ctx;
} batch_task_arg_t;

typedef struct {
    char* canonical_smiles;
    char* error_msg;
    bool success;
    int row_idx;
} batch_task_result_t;

/* Worker function for canonicalization */
void* csv_batch_worker(void* arg) {
    batch_task_arg_t* task = (batch_task_arg_t*)arg;
    batch_task_result_t* result = (batch_task_result_t*)calloc(1, sizeof(batch_task_result_t));
    if (!result) return NULL;

    result->row_idx = task->row_idx;

    if (!task->smiles || task->smiles[0] == '\0') {
        result->success = false;
        result->error_msg = strdup("Empty SMILES");
        result->canonical_smiles = strdup("");
    } else {
        char error_buf[256];
        char* canonical = smiles_canonicalize(task->smiles, NULL, error_buf, sizeof(error_buf));

        if (canonical) {
            result->success = true;
            result->canonical_smiles = canonical;
            result->error_msg = NULL;
        } else {
            result->success = false;
            result->canonical_smiles = strdup("");
            result->error_msg = strdup(error_buf);
        }
    }

    return result;
}

csv_batch_context_t* csv_batch_context_create(const char* input_file,
                                               const char* output_file,
                                               const char* smiles_column,
                                               const char* output_column,
                                               int num_threads) {
    if (!input_file || !output_file || !smiles_column || !output_column) {
        return NULL;
    }

    csv_batch_context_t* ctx = (csv_batch_context_t*)calloc(1, sizeof(csv_batch_context_t));
    if (!ctx) return NULL;

    ctx->input_file = input_file;
    ctx->output_file = output_file;
    ctx->smiles_column = smiles_column;
    ctx->output_column = output_column;
    ctx->num_threads = (num_threads > 0) ? num_threads : parallel_get_num_cores();
    ctx->has_header = true;

    pthread_mutex_init(&ctx->progress_mutex, NULL);

    /* Count rows */
    ctx->total_rows = csv_count_rows(input_file);
    if (ctx->total_rows > 0) ctx->total_rows--;  /* Subtract header */

    ctx->processed_rows = 0;
    ctx->error_rows = 0;

    return ctx;
}

void csv_batch_context_free(csv_batch_context_t* ctx) {
    if (!ctx) return;

    pthread_mutex_destroy(&ctx->progress_mutex);

    if (ctx->results) {
        for (size_t i = 0; i < ctx->total_rows; i++) {
            if (ctx->results[i]) free(ctx->results[i]);
        }
        free(ctx->results);
    }

    if (ctx->errors) {
        for (size_t i = 0; i < ctx->total_rows; i++) {
            if (ctx->errors[i]) free(ctx->errors[i]);
        }
        free(ctx->errors);
    }

    if (ctx->success) free(ctx->success);

    free(ctx);
}

cchem_status_t csv_batch_canonicalize(csv_batch_context_t* ctx) {
    if (!ctx) return CCHEM_ERROR_INVALID_INPUT;

    /* Allocate results arrays */
    ctx->results = (char**)calloc(ctx->total_rows, sizeof(char*));
    ctx->success = (bool*)calloc(ctx->total_rows, sizeof(bool));
    ctx->errors = (char**)calloc(ctx->total_rows, sizeof(char*));

    if (!ctx->results || !ctx->success || !ctx->errors) {
        return CCHEM_ERROR_MEMORY;
    }

    /* Open input file */
    csv_reader_t* reader = csv_reader_create(ctx->input_file);
    if (!reader) {
        return CCHEM_ERROR_FILE_IO;
    }

    /* Read header */
    csv_row_t header;
    memset(&header, 0, sizeof(header));

    if (csv_reader_next(reader, &header) != CSV_OK) {
        csv_reader_free(reader);
        return CCHEM_ERROR_FILE_IO;
    }

    int smiles_col = csv_find_column(&header, ctx->smiles_column);
    if (smiles_col < 0) {
        fprintf(stderr, "Error: Column '%s' not found in input file\n", ctx->smiles_column);
        csv_reader_free(reader);
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Create thread pool */
    thread_pool_t* pool = thread_pool_create(ctx->num_threads);
    if (!pool) {
        csv_reader_free(reader);
        return CCHEM_ERROR_THREAD;
    }

    /* Pre-allocate capacity to avoid reallocation race conditions */
    if (thread_pool_ensure_capacity(pool, ctx->total_rows + 1) < 0) {
        thread_pool_free(pool);
        csv_reader_free(reader);
        return CCHEM_ERROR_MEMORY;
    }

    /* Create progress bar */
    progress_config_t prog_config = PROGRESS_CONFIG_DEFAULT;
    prog_config.prefix = "Canonicalizing";
    progress_t* progress = progress_create(ctx->total_rows, &prog_config);

    /* Read all rows and submit tasks */
    batch_task_arg_t* task_args = (batch_task_arg_t*)calloc(ctx->total_rows, sizeof(batch_task_arg_t));
    if (!task_args) {
        progress_free(progress);
        thread_pool_free(pool);
        csv_reader_free(reader);
        return CCHEM_ERROR_MEMORY;
    }

    /* Read all SMILES and store in document for output */
    csv_document_t* doc = csv_document_create();
    if (!doc) {
        free(task_args);
        progress_free(progress);
        thread_pool_free(pool);
        csv_reader_free(reader);
        return CCHEM_ERROR_MEMORY;
    }

    /* Copy header and add output column */
    doc->has_header = true;
    doc->header.capacity = header.num_fields + 2;
    doc->header.fields = (char**)calloc(doc->header.capacity, sizeof(char*));
    for (int i = 0; i < header.num_fields; i++) {
        doc->header.fields[doc->header.num_fields++] = strdup(header.fields[i]);
    }
    doc->header.fields[doc->header.num_fields++] = strdup(ctx->output_column);

    csv_row_t row;
    memset(&row, 0, sizeof(row));
    int row_idx = 0;

    while (csv_reader_next(reader, &row) == CSV_OK && row_idx < (int)ctx->total_rows) {
        /* Store row in document */
        if (doc->num_rows >= doc->rows_capacity) {
            int new_cap = doc->rows_capacity * 2;
            csv_row_t* new_rows = (csv_row_t*)realloc(doc->rows, new_cap * sizeof(csv_row_t));
            if (new_rows) {
                doc->rows = new_rows;
                doc->rows_capacity = new_cap;
            }
        }

        doc->rows[doc->num_rows].capacity = row.num_fields + 2;
        doc->rows[doc->num_rows].num_fields = 0;
        doc->rows[doc->num_rows].fields = (char**)calloc(doc->rows[doc->num_rows].capacity, sizeof(char*));

        for (int i = 0; i < row.num_fields; i++) {
            doc->rows[doc->num_rows].fields[doc->rows[doc->num_rows].num_fields++] = strdup(row.fields[i]);
        }
        doc->num_rows++;

        /* Submit task */
        const char* smiles = csv_row_get_field(&row, smiles_col);
        task_args[row_idx].row_idx = row_idx;
        task_args[row_idx].smiles = smiles ? strdup(smiles) : strdup("");
        task_args[row_idx].ctx = ctx;

        thread_pool_submit(pool, csv_batch_worker, &task_args[row_idx]);

        row_idx++;
    }

    csv_reader_free(reader);

    /* Wait for completion with progress updates */
    while (thread_pool_num_completed(pool) < (int)ctx->total_rows) {
        int completed = thread_pool_num_completed(pool);
        progress_update(progress, completed);
        usleep(50000);  /* 50ms */
    }

    /* Collect results */
    for (int i = 0; i < (int)ctx->total_rows; i++) {
        batch_task_result_t* result = (batch_task_result_t*)thread_pool_get_result(pool, i);
        if (result) {
            int idx = result->row_idx;
            ctx->results[idx] = result->canonical_smiles;
            ctx->success[idx] = result->success;
            ctx->errors[idx] = result->error_msg;

            if (result->success) {
                ctx->processed_rows++;
            } else {
                ctx->error_rows++;
            }

            /* Add to document row */
            if (idx < doc->num_rows) {
                doc->rows[idx].fields[doc->rows[idx].num_fields++] =
                    strdup(result->canonical_smiles ? result->canonical_smiles : "");
            }

            free(result);
        } else {
            /* Result allocation failed - record as error */
            ctx->results[i] = strdup("");
            ctx->success[i] = false;
            ctx->errors[i] = strdup("Memory allocation failed");
            ctx->error_rows++;

            if (i < doc->num_rows) {
                doc->rows[i].fields[doc->rows[i].num_fields++] = strdup("");
            }
        }

        if (task_args[i].smiles) free(task_args[i].smiles);
    }

    progress_finish(progress);
    progress_free(progress);
    thread_pool_free(pool);
    free(task_args);

    /* Write output */
    csv_status_t csv_status = csv_document_write(doc, ctx->output_file);
    csv_document_free(doc);

    /* Free header */
    for (int i = 0; i < header.num_fields; i++) {
        if (header.fields[i]) free(header.fields[i]);
    }
    if (header.fields) free(header.fields);

    if (csv_status != CSV_OK) {
        return CCHEM_ERROR_FILE_IO;
    }

    return CCHEM_OK;
}

void csv_batch_get_stats(const csv_batch_context_t* ctx,
                         size_t* total, size_t* success, size_t* errors) {
    if (!ctx) return;

    if (total) *total = ctx->total_rows;
    if (success) *success = ctx->processed_rows;
    if (errors) *errors = ctx->error_rows;
}
