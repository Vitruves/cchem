/**
 * @file threading.h
 * @brief Threading utilities for cchem
 *
 * This module provides:
 * - Thread pool for parallel task execution
 * - Bounded queue for producer-consumer patterns
 * - CSV batch processing utilities
 */

#ifndef CCHEM_THREADING_H
#define CCHEM_THREADING_H

#include <stddef.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <pthread.h>
#include "canonicalizer/types.h"

/* ============================================================================
 * Bounded Queue (Thread-Safe Ring Buffer)
 *
 * Used for producer-consumer patterns in pipeline streaming.
 * Supports blocking push/pop with close semantics for EOF signaling.
 * ============================================================================ */

/**
 * Bounded thread-safe queue (ring buffer)
 */
typedef struct {
    void** items;
    int capacity;
    int head;
    int tail;
    int count;
    pthread_mutex_t mutex;
    pthread_cond_t not_full;
    pthread_cond_t not_empty;
    bool closed;
} bounded_queue_t;

/**
 * Create a bounded queue with specified capacity
 * @param capacity Maximum number of items in queue
 * @return New queue or NULL on failure
 */
bounded_queue_t* queue_create(int capacity);

/**
 * Free queue and all resources
 * Does NOT free items still in queue - caller must drain first
 */
void queue_free(bounded_queue_t* q);

/**
 * Push item to queue (blocks if full)
 * @param q Queue
 * @param item Item to push (must not be NULL)
 * @return true if pushed, false if queue is closed
 */
bool queue_push(bounded_queue_t* q, void* item);

/**
 * Pop item from queue (blocks if empty)
 * @param q Queue
 * @return Item or NULL if queue is closed and empty
 */
void* queue_pop(bounded_queue_t* q);

/**
 * Try to pop without blocking
 * @param q Queue
 * @param out Output pointer for item
 * @return true if item popped, false if empty
 */
bool queue_try_pop(bounded_queue_t* q, void** out);

/**
 * Close queue - signals EOF to all consumers
 * Blocked pushers will return false
 * Blocked poppers will return NULL after queue drains
 */
void queue_close(bounded_queue_t* q);

/**
 * Check if queue is closed
 */
bool queue_is_closed(bounded_queue_t* q);

/**
 * Get current item count (approximate, may change)
 */
int queue_count(bounded_queue_t* q);

/* ============================================================================
 * Thread Pool
 *
 * Work-stealing thread pool for parallel task execution.
 * Supports dynamic task submission and result retrieval.
 * ============================================================================ */

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

/* Get number of available CPU cores */
int parallel_get_num_cores(void);

/* ============================================================================
 * CSV Batch Processing
 *
 * High-level utilities for parallel CSV processing.
 * ============================================================================ */

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

#endif /* CCHEM_THREADING_H */
