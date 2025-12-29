/**
 * @file parquet.h
 * @brief Parquet file reading and writing utilities
 *
 * This module provides parquet file I/O using the carquet library.
 * API is designed to mirror csv.h for easy integration.
 */

#ifndef CCHEM_PARQUET_H
#define CCHEM_PARQUET_H

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

/* Status codes */
typedef enum {
    PARQUET_OK = 0,
    PARQUET_ERROR_IO = -1,
    PARQUET_ERROR_MEMORY = -2,
    PARQUET_ERROR_PARSE = -3,
    PARQUET_ERROR_INVALID = -4,
    PARQUET_ERROR_COLUMN_NOT_FOUND = -5,
    PARQUET_ERROR_TYPE_MISMATCH = -6,
    PARQUET_ERROR_END_OF_DATA = -7
} parquet_status_t;

/* Forward declarations for carquet types */
struct carquet_reader;
struct carquet_writer;
struct carquet_batch_reader;
struct carquet_row_batch;
struct carquet_schema;

/* ============================================================================
 * Parquet Row (for row-by-row access)
 * ============================================================================ */

/**
 * Parquet row - holds string representations of field values
 */
typedef struct {
    char** fields;          /* Array of string field values */
    int num_fields;
    int capacity;
} parquet_row_t;

/* Create parquet row */
parquet_row_t* parquet_row_create(void);

/* Free parquet row */
void parquet_row_free(parquet_row_t* row);

/* Clear row (reuse memory) */
void parquet_row_clear(parquet_row_t* row);

/* Add field to row */
parquet_status_t parquet_row_add_field(parquet_row_t* row, const char* field);

/* Get field by index */
const char* parquet_row_get_field(const parquet_row_t* row, int idx);

/* ============================================================================
 * Parquet Reader
 * ============================================================================ */

/**
 * Parquet reader state
 */
typedef struct parquet_reader {
    struct carquet_reader* reader;
    struct carquet_batch_reader* batch_reader;
    struct carquet_row_batch* current_batch;
    int64_t current_row_in_batch;
    int64_t total_rows;
    int64_t rows_read;
    int32_t num_columns;
    char** column_names;        /* Cached column names */
    int* string_column_indices; /* Indices of string columns */
    int num_string_columns;
    bool at_eof;
    char error_msg[256];
} parquet_reader_t;

/**
 * Create parquet reader
 * @param filename Path to parquet file
 * @return Reader or NULL on failure
 */
parquet_reader_t* parquet_reader_create(const char* filename);

/**
 * Free parquet reader
 */
void parquet_reader_free(parquet_reader_t* reader);

/**
 * Read next row into provided row structure
 * @param reader Parquet reader
 * @param row Output row (must be pre-allocated)
 * @return PARQUET_OK on success, PARQUET_ERROR_END_OF_DATA at EOF
 */
parquet_status_t parquet_reader_next(parquet_reader_t* reader, parquet_row_t* row);

/**
 * Check if at end of file
 */
bool parquet_reader_at_end(const parquet_reader_t* reader);

/**
 * Get total row count
 */
int64_t parquet_reader_num_rows(const parquet_reader_t* reader);

/**
 * Get number of columns
 */
int32_t parquet_reader_num_columns(const parquet_reader_t* reader);

/**
 * Get column name by index
 */
const char* parquet_reader_column_name(const parquet_reader_t* reader, int idx);

/**
 * Find column index by name
 * @return Column index or -1 if not found
 */
int parquet_reader_find_column(const parquet_reader_t* reader, const char* name);

/**
 * Get error message
 */
const char* parquet_reader_get_error(const parquet_reader_t* reader);

/* ============================================================================
 * Parquet Writer
 * ============================================================================ */

/**
 * Column data types for parquet writer
 */
typedef enum {
    PARQUET_COL_STRING = 0,    /* String/text data (BYTE_ARRAY) */
    PARQUET_COL_DOUBLE = 1,    /* Floating point (DOUBLE) */
    PARQUET_COL_INT64 = 2      /* Integer (INT64) */
} parquet_col_type_t;

/**
 * Parquet writer state
 */
typedef struct parquet_writer {
    struct carquet_writer* writer;
    struct carquet_schema* schema;
    char** column_names;
    parquet_col_type_t* column_types;  /* Type of each column */
    int32_t num_columns;
    int64_t rows_written;

    /* Row buffer for batch writing */
    char*** row_buffer;         /* Array of rows, each row is array of strings */
    int buffer_size;
    int buffer_capacity;

    char error_msg[256];
} parquet_writer_t;

/**
 * Create parquet writer with typed columns
 * @param filename Output file path
 * @param column_names Array of column names
 * @param column_types Array of column types (NULL = all strings)
 * @param num_columns Number of columns
 * @return Writer or NULL on failure
 */
parquet_writer_t* parquet_writer_create(const char* filename,
                                         const char** column_names,
                                         const parquet_col_type_t* column_types,
                                         int num_columns);

/**
 * Free parquet writer (flushes and closes file)
 */
void parquet_writer_free(parquet_writer_t* writer);

/**
 * Write a single row
 * @param writer Parquet writer
 * @param fields Array of field values (as strings)
 * @param num_fields Number of fields (must match num_columns)
 * @return PARQUET_OK on success
 */
parquet_status_t parquet_writer_write_row(parquet_writer_t* writer,
                                           const char** fields,
                                           int num_fields);

/**
 * Flush buffered data to file
 */
parquet_status_t parquet_writer_flush(parquet_writer_t* writer);

/**
 * Get error message
 */
const char* parquet_writer_get_error(const parquet_writer_t* writer);

/* ============================================================================
 * Streaming Parquet Reader (for pipeline processing)
 * ============================================================================ */

/**
 * Streaming parquet reader - optimized for sequential row access
 * Similar to csv_mmap_reader but for parquet files
 */
typedef struct parquet_stream_reader {
    struct carquet_reader* reader;
    struct carquet_batch_reader* batch_reader;
    struct carquet_row_batch* current_batch;

    int64_t current_row_in_batch;
    int64_t batch_num_rows;
    int64_t total_rows;
    int64_t rows_read;

    int32_t num_columns;
    char** column_names;

    /* Thread-local string buffers for value conversion */
    char** field_buffers;
    size_t* field_buffer_sizes;

    bool at_eof;
    char error_msg[256];
} parquet_stream_reader_t;

/**
 * Create streaming parquet reader
 */
parquet_stream_reader_t* parquet_stream_reader_create(const char* filename);

/**
 * Free streaming parquet reader
 */
void parquet_stream_reader_free(parquet_stream_reader_t* reader);

/**
 * Get next row as array of string pointers
 * @param reader Stream reader
 * @param fields Pre-allocated array of char* (size = num_columns)
 * @param max_fields Maximum number of fields to read
 * @return Number of fields read, or -1 at EOF/error
 *
 * NOTE: Field pointers are valid until next call to this function.
 *       Caller must copy if needed.
 */
int parquet_stream_reader_next(parquet_stream_reader_t* reader,
                                const char** fields,
                                int max_fields);

/**
 * Check if at end of data
 */
bool parquet_stream_reader_at_end(const parquet_stream_reader_t* reader);

/**
 * Get total row count
 */
int64_t parquet_stream_reader_num_rows(const parquet_stream_reader_t* reader);

/**
 * Get number of columns
 */
int32_t parquet_stream_reader_num_columns(const parquet_stream_reader_t* reader);

/**
 * Get column name by index
 */
const char* parquet_stream_reader_column_name(const parquet_stream_reader_t* reader, int idx);

/**
 * Find column index by name
 */
int parquet_stream_reader_find_column(const parquet_stream_reader_t* reader, const char* name);

/**
 * Get error message
 */
const char* parquet_stream_reader_get_error(const parquet_stream_reader_t* reader);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Check if filename has parquet extension (.parquet or .pq)
 */
bool is_parquet_file(const char* filename);

/**
 * Count rows in parquet file (fast, reads only metadata)
 */
int64_t parquet_count_rows(const char* filename);

#endif /* CCHEM_PARQUET_H */
