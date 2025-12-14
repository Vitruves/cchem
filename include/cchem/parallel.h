/**
 * @file parallel.h
 * @brief Parallel processing utilities using pthreads
 */

#ifndef CCHEM_PARALLEL_H
#define CCHEM_PARALLEL_H

#include <stddef.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <pthread.h>
#include "canonicalizer/types.h"

/* Task function type */
typedef void* (*parallel_task_fn)(void* arg);

/* Thread pool */
typedef struct {
    pthread_t* threads;
    int num_threads;
    bool* thread_active;

    /* Task queue */
    struct {
        parallel_task_fn fn;
        void* arg;
        void* result;
        bool completed;
        bool has_error;
        char error_msg[256];
    }* tasks;
    int num_tasks;
    int tasks_capacity;
    int next_task;
    atomic_int completed_tasks;  /* Atomic for lock-free progress reads */

    /* Synchronization */
    pthread_mutex_t mutex;
    pthread_cond_t task_available;
    pthread_cond_t task_completed;

    /* State */
    bool shutdown;
    bool running;
} thread_pool_t;

/* Create thread pool */
thread_pool_t* thread_pool_create(int num_threads);

/* Free thread pool */
void thread_pool_free(thread_pool_t* pool);

/* Pre-allocate task capacity to avoid reallocation during execution */
int thread_pool_ensure_capacity(thread_pool_t* pool, int capacity);

/* Submit task to pool */
int thread_pool_submit(thread_pool_t* pool, parallel_task_fn fn, void* arg);

/* Wait for all tasks to complete */
void thread_pool_wait_all(thread_pool_t* pool);

/* Get task result */
void* thread_pool_get_result(thread_pool_t* pool, int task_id);

/* Check if task completed */
bool thread_pool_task_completed(thread_pool_t* pool, int task_id);

/* Get number of completed tasks */
int thread_pool_num_completed(thread_pool_t* pool);

/* Clear completed tasks */
void thread_pool_clear_completed(thread_pool_t* pool);

/* CSV batch processing context */
typedef struct {
    const char* input_file;
    const char* output_file;
    const char* smiles_column;
    const char* output_column;
    int num_threads;
    bool has_header;

    /* Progress tracking */
    size_t total_rows;
    size_t processed_rows;
    size_t error_rows;
    pthread_mutex_t progress_mutex;

    /* Results */
    char** results;          /* Canonical SMILES for each row */
    bool* success;           /* Success flag for each row */
    char** errors;           /* Error messages for failures */
} csv_batch_context_t;

/* Create batch processing context */
csv_batch_context_t* csv_batch_context_create(const char* input_file,
                                               const char* output_file,
                                               const char* smiles_column,
                                               const char* output_column,
                                               int num_threads);

/* Free batch processing context */
void csv_batch_context_free(csv_batch_context_t* ctx);

/* Run batch canonicalization */
cchem_status_t csv_batch_canonicalize(csv_batch_context_t* ctx);

/* Get processing statistics */
void csv_batch_get_stats(const csv_batch_context_t* ctx,
                         size_t* total, size_t* success, size_t* errors);

/* Worker function for batch processing */
void* csv_batch_worker(void* arg);

/* Get number of available CPU cores */
int parallel_get_num_cores(void);

#endif /* CCHEM_PARALLEL_H */
