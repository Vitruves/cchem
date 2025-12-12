/**
 * @file csv.h
 * @brief CSV file reading and writing
 */

#ifndef CCHEM_CSV_H
#define CCHEM_CSV_H

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>

/* Status codes */
typedef enum {
    CSV_OK = 0,
    CSV_ERROR_IO = -1,
    CSV_ERROR_MEMORY = -2,
    CSV_ERROR_PARSE = -3,
    CSV_ERROR_INVALID = -4
} csv_status_t;

/* CSV row */
typedef struct {
    char** fields;
    int num_fields;
    int capacity;
} csv_row_t;

/* CSV document */
typedef struct {
    csv_row_t* rows;
    int num_rows;
    int rows_capacity;
    csv_row_t header;
    bool has_header;
    char delimiter;
    char quote_char;
} csv_document_t;

/* CSV reader state */
typedef struct {
    FILE* file;
    char* line_buffer;
    size_t line_buffer_size;
    char delimiter;
    char quote_char;
    int current_line;
    bool at_eof;
    char error_msg[256];
} csv_reader_t;

/* CSV writer state */
typedef struct {
    FILE* file;
    char delimiter;
    char quote_char;
    bool needs_quotes;
    char error_msg[256];
} csv_writer_t;

/* Create CSV row */
csv_row_t* csv_row_create(void);

/* Free CSV row */
void csv_row_free(csv_row_t* row);

/* Add field to row */
csv_status_t csv_row_add_field(csv_row_t* row, const char* field);

/* Get field by index */
const char* csv_row_get_field(const csv_row_t* row, int idx);

/* Get field by column name (requires header) */
const char* csv_row_get_field_by_name(const csv_row_t* row,
                                       const csv_row_t* header,
                                       const char* name);

/* Find column index by name */
int csv_find_column(const csv_row_t* header, const char* name);

/* Create CSV document */
csv_document_t* csv_document_create(void);

/* Free CSV document */
void csv_document_free(csv_document_t* doc);

/* Read entire CSV file into document */
csv_status_t csv_document_read(csv_document_t* doc, const char* filename,
                               bool first_row_header);

/* Write document to file */
csv_status_t csv_document_write(const csv_document_t* doc, const char* filename);

/* Create CSV reader */
csv_reader_t* csv_reader_create(const char* filename);

/* Free CSV reader */
void csv_reader_free(csv_reader_t* reader);

/* Read next row */
csv_status_t csv_reader_next(csv_reader_t* reader, csv_row_t* row);

/* Check if at end */
bool csv_reader_at_end(const csv_reader_t* reader);

/* Get current line number */
int csv_reader_get_line(const csv_reader_t* reader);

/* Get error message */
const char* csv_reader_get_error(const csv_reader_t* reader);

/* Create CSV writer */
csv_writer_t* csv_writer_create(const char* filename);

/* Free CSV writer */
void csv_writer_free(csv_writer_t* writer);

/* Write row to file */
csv_status_t csv_writer_write_row(csv_writer_t* writer, const csv_row_t* row);

/* Write fields to file */
csv_status_t csv_writer_write_fields(csv_writer_t* writer,
                                     const char** fields, int num_fields);

/* Get error message */
const char* csv_writer_get_error(const csv_writer_t* writer);

/* Utility: count rows in file */
int csv_count_rows(const char* filename);

/* Utility: sanitize numeric field - returns "0" if empty or non-numeric
 * This is used to ensure CSV output has no empty cells for numeric columns.
 * Returns the original field if valid, or "0" if empty/non-numeric.
 * Note: Returns static "0" string for invalid values - do not free.
 */
const char* csv_sanitize_numeric(const char* field);

/* ============================================================================
 * Memory-Mapped CSV Reader (High Performance)
 *
 * Uses mmap() to map entire file into memory, eliminating syscall overhead.
 * Ideal for batch processing millions of SMILES strings.
 * ============================================================================ */

/* Memory-mapped CSV reader state */
typedef struct csv_mmap_reader {
    char* data;              /* mmap'd file data */
    size_t size;             /* File size */
    size_t pos;              /* Current read position */
    char delimiter;
    char quote_char;
    int current_line;
    bool at_eof;
    int fd;                  /* File descriptor (for munmap) */
    char error_msg[256];
} csv_mmap_reader_t;

/* Create memory-mapped CSV reader
 * Returns NULL on failure (check errno)
 * Much faster than csv_reader_create for large files
 */
csv_mmap_reader_t* csv_mmap_reader_create(const char* filename);

/* Free memory-mapped CSV reader (unmaps file) */
void csv_mmap_reader_free(csv_mmap_reader_t* reader);

/* Get next line as zero-copy pointer (does NOT allocate)
 * Returns pointer to start of line, sets *len to line length
 * Returns NULL at EOF
 * IMPORTANT: The returned pointer is NOT null-terminated!
 *            Line includes any trailing \r\n
 */
const char* csv_mmap_reader_next_line(csv_mmap_reader_t* reader, size_t* len);

/* Parse a line into fields (zero-copy where possible)
 * fields: pre-allocated array of char* pointers
 * max_fields: capacity of fields array
 * line: pointer from csv_mmap_reader_next_line
 * len: length from csv_mmap_reader_next_line
 * Returns number of fields parsed, or -1 on error
 *
 * NOTE: For simple CSV (no quotes), fields point directly into mmap'd data.
 *       For quoted fields, uses thread-local buffer.
 */
int csv_mmap_parse_line(const char* line, size_t len, char delimiter,
                        const char** fields, int max_fields);

/* Check if at end */
bool csv_mmap_reader_at_end(const csv_mmap_reader_t* reader);

/* Get current line number */
int csv_mmap_reader_get_line(const csv_mmap_reader_t* reader);

/* Get error message */
const char* csv_mmap_reader_get_error(const csv_mmap_reader_t* reader);

#endif /* CCHEM_CSV_H */
