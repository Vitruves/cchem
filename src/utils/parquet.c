/**
 * @file parquet.c
 * @brief Parquet file reading and writing implementation
 */

#include "cchem/utils/parquet.h"
#include <carquet/carquet.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Parquet Row Implementation
 * ============================================================================ */

parquet_row_t* parquet_row_create(void) {
    parquet_row_t* row = (parquet_row_t*)calloc(1, sizeof(parquet_row_t));
    if (!row) return NULL;

    row->capacity = 16;
    row->fields = (char**)calloc(row->capacity, sizeof(char*));
    if (!row->fields) {
        free(row);
        return NULL;
    }
    row->num_fields = 0;
    return row;
}

void parquet_row_free(parquet_row_t* row) {
    if (!row) return;
    if (row->fields) {
        for (int i = 0; i < row->num_fields; i++) {
            free(row->fields[i]);
        }
        free(row->fields);
    }
    free(row);
}

void parquet_row_clear(parquet_row_t* row) {
    if (!row) return;
    for (int i = 0; i < row->num_fields; i++) {
        free(row->fields[i]);
        row->fields[i] = NULL;
    }
    row->num_fields = 0;
}

parquet_status_t parquet_row_add_field(parquet_row_t* row, const char* field) {
    if (!row) return PARQUET_ERROR_INVALID;

    /* Grow capacity if needed */
    if (row->num_fields >= row->capacity) {
        int new_capacity = row->capacity * 2;
        char** new_fields = (char**)realloc(row->fields, new_capacity * sizeof(char*));
        if (!new_fields) return PARQUET_ERROR_MEMORY;
        row->fields = new_fields;
        row->capacity = new_capacity;
    }

    row->fields[row->num_fields] = field ? strdup(field) : strdup("");
    if (!row->fields[row->num_fields]) return PARQUET_ERROR_MEMORY;
    row->num_fields++;
    return PARQUET_OK;
}

const char* parquet_row_get_field(const parquet_row_t* row, int idx) {
    if (!row || idx < 0 || idx >= row->num_fields) return NULL;
    return row->fields[idx];
}

/* ============================================================================
 * Parquet Reader Implementation
 * ============================================================================ */

parquet_reader_t* parquet_reader_create(const char* filename) {
    if (!filename) return NULL;

    parquet_reader_t* reader = (parquet_reader_t*)calloc(1, sizeof(parquet_reader_t));
    if (!reader) return NULL;

    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Open parquet file */
    reader->reader = carquet_reader_open(filename, NULL, &err);
    if (!reader->reader) {
        snprintf(reader->error_msg, sizeof(reader->error_msg),
                 "Failed to open parquet file: %s", err.message);
        free(reader);
        return NULL;
    }

    /* Get metadata */
    reader->total_rows = carquet_reader_num_rows(reader->reader);
    reader->num_columns = carquet_reader_num_columns(reader->reader);

    /* Cache column names */
    const carquet_schema_t* schema = carquet_reader_schema(reader->reader);
    reader->column_names = (char**)calloc(reader->num_columns, sizeof(char*));
    if (!reader->column_names) {
        snprintf(reader->error_msg, sizeof(reader->error_msg), "Memory allocation failed");
        carquet_reader_close(reader->reader);
        free(reader);
        return NULL;
    }

    for (int i = 0; i < reader->num_columns; i++) {
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, i + 1);
        if (node) {
            reader->column_names[i] = strdup(carquet_schema_node_name(node));
        } else {
            char buf[32];
            snprintf(buf, sizeof(buf), "column_%d", i);
            reader->column_names[i] = strdup(buf);
        }
    }

    /* Create batch reader */
    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = 10000;

    reader->batch_reader = carquet_batch_reader_create(reader->reader, &config, &err);
    if (!reader->batch_reader) {
        snprintf(reader->error_msg, sizeof(reader->error_msg),
                 "Failed to create batch reader: %s", err.message);
        for (int i = 0; i < reader->num_columns; i++) free(reader->column_names[i]);
        free(reader->column_names);
        carquet_reader_close(reader->reader);
        free(reader);
        return NULL;
    }

    reader->current_batch = NULL;
    reader->current_row_in_batch = 0;
    reader->rows_read = 0;
    reader->at_eof = false;

    return reader;
}

void parquet_reader_free(parquet_reader_t* reader) {
    if (!reader) return;

    if (reader->current_batch) {
        carquet_row_batch_free(reader->current_batch);
    }
    if (reader->batch_reader) {
        carquet_batch_reader_free(reader->batch_reader);
    }
    if (reader->reader) {
        carquet_reader_close(reader->reader);
    }
    if (reader->column_names) {
        for (int i = 0; i < reader->num_columns; i++) {
            free(reader->column_names[i]);
        }
        free(reader->column_names);
    }
    free(reader->string_column_indices);
    free(reader);
}

/* Helper to convert parquet value to string */
static char* parquet_value_to_string(const void* data, int64_t row_idx,
                                      carquet_physical_type_t type,
                                      const uint8_t* null_bitmap) {
    /* Check for null */
    if (null_bitmap) {
        bool is_null = !(null_bitmap[row_idx / 8] & (1 << (row_idx % 8)));
        if (is_null) {
            return strdup("");
        }
    }

    char buf[64];

    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN: {
            const uint8_t* vals = (const uint8_t*)data;
            snprintf(buf, sizeof(buf), "%s", vals[row_idx] ? "true" : "false");
            break;
        }
        case CARQUET_PHYSICAL_INT32: {
            const int32_t* vals = (const int32_t*)data;
            snprintf(buf, sizeof(buf), "%d", vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_INT64: {
            const int64_t* vals = (const int64_t*)data;
            snprintf(buf, sizeof(buf), "%lld", (long long)vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_FLOAT: {
            const float* vals = (const float*)data;
            snprintf(buf, sizeof(buf), "%.6g", vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_DOUBLE: {
            const double* vals = (const double*)data;
            snprintf(buf, sizeof(buf), "%.6g", vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_BYTE_ARRAY: {
            const carquet_byte_array_t* vals = (const carquet_byte_array_t*)data;
            /* Return string directly (assume UTF-8) */
            char* str = (char*)malloc(vals[row_idx].length + 1);
            if (str) {
                memcpy(str, vals[row_idx].data, vals[row_idx].length);
                str[vals[row_idx].length] = '\0';
            }
            return str;
        }
        default:
            return strdup("");
    }

    return strdup(buf);
}

parquet_status_t parquet_reader_next(parquet_reader_t* reader, parquet_row_t* row) {
    if (!reader || !row) return PARQUET_ERROR_INVALID;
    if (reader->at_eof) return PARQUET_ERROR_END_OF_DATA;

    parquet_row_clear(row);

    /* Load next batch if needed */
    if (!reader->current_batch ||
        reader->current_row_in_batch >= carquet_row_batch_num_rows(reader->current_batch)) {

        if (reader->current_batch) {
            carquet_row_batch_free(reader->current_batch);
            reader->current_batch = NULL;
        }

        carquet_status_t status = carquet_batch_reader_next(reader->batch_reader,
                                                             &reader->current_batch);
        if (status != CARQUET_OK || !reader->current_batch) {
            reader->at_eof = true;
            return PARQUET_ERROR_END_OF_DATA;
        }
        reader->current_row_in_batch = 0;
    }

    /* Get schema for type info */
    const carquet_schema_t* schema = carquet_reader_schema(reader->reader);

    /* Extract values for current row */
    for (int col = 0; col < reader->num_columns; col++) {
        const void* data;
        const uint8_t* null_bitmap;
        int64_t num_values;

        carquet_status_t status = carquet_row_batch_column(reader->current_batch, col,
                                                            &data, &null_bitmap, &num_values);
        if (status != CARQUET_OK) {
            parquet_row_add_field(row, "");
            continue;
        }

        /* Get column type */
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, col + 1);
        carquet_physical_type_t ptype = node ? carquet_schema_node_physical_type(node)
                                              : CARQUET_PHYSICAL_BYTE_ARRAY;

        char* value = parquet_value_to_string(data, reader->current_row_in_batch,
                                               ptype, null_bitmap);
        parquet_row_add_field(row, value ? value : "");
        free(value);
    }

    reader->current_row_in_batch++;
    reader->rows_read++;

    return PARQUET_OK;
}

bool parquet_reader_at_end(const parquet_reader_t* reader) {
    return reader ? reader->at_eof : true;
}

int64_t parquet_reader_num_rows(const parquet_reader_t* reader) {
    return reader ? reader->total_rows : 0;
}

int32_t parquet_reader_num_columns(const parquet_reader_t* reader) {
    return reader ? reader->num_columns : 0;
}

const char* parquet_reader_column_name(const parquet_reader_t* reader, int idx) {
    if (!reader || idx < 0 || idx >= reader->num_columns) return NULL;
    return reader->column_names[idx];
}

int parquet_reader_find_column(const parquet_reader_t* reader, const char* name) {
    if (!reader || !name) return -1;
    for (int i = 0; i < reader->num_columns; i++) {
        if (reader->column_names[i] && strcmp(reader->column_names[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

const char* parquet_reader_get_error(const parquet_reader_t* reader) {
    return reader ? reader->error_msg : "Invalid reader";
}

/* ============================================================================
 * Parquet Writer Implementation
 * ============================================================================ */

#define PARQUET_WRITE_BUFFER_SIZE 10000

parquet_writer_t* parquet_writer_create(const char* filename,
                                         const char** column_names,
                                         const parquet_col_type_t* column_types,
                                         int num_columns) {
    if (!filename || !column_names || num_columns <= 0) return NULL;

    parquet_writer_t* writer = (parquet_writer_t*)calloc(1, sizeof(parquet_writer_t));
    if (!writer) return NULL;

    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Create schema */
    writer->schema = carquet_schema_create(&err);
    if (!writer->schema) {
        snprintf(writer->error_msg, sizeof(writer->error_msg),
                 "Failed to create schema: %s", err.message);
        free(writer);
        return NULL;
    }

    writer->num_columns = num_columns;
    writer->column_names = (char**)calloc(num_columns, sizeof(char*));
    writer->column_types = (parquet_col_type_t*)calloc(num_columns, sizeof(parquet_col_type_t));
    if (!writer->column_names || !writer->column_types) {
        free(writer->column_names);
        free(writer->column_types);
        carquet_schema_free(writer->schema);
        free(writer);
        return NULL;
    }

    /* Add columns to schema with appropriate types */
    for (int i = 0; i < num_columns; i++) {
        writer->column_names[i] = strdup(column_names[i]);
        writer->column_types[i] = column_types ? column_types[i] : PARQUET_COL_STRING;

        carquet_physical_type_t ptype;
        switch (writer->column_types[i]) {
            case PARQUET_COL_DOUBLE:
                ptype = CARQUET_PHYSICAL_DOUBLE;
                break;
            case PARQUET_COL_INT64:
                ptype = CARQUET_PHYSICAL_INT64;
                break;
            case PARQUET_COL_STRING:
            default:
                ptype = CARQUET_PHYSICAL_BYTE_ARRAY;
                break;
        }

        carquet_status_t status = carquet_schema_add_column(
            writer->schema,
            column_names[i],
            ptype,
            NULL,
            CARQUET_REPETITION_REQUIRED,
            0
        );

        if (status != CARQUET_OK) {
            snprintf(writer->error_msg, sizeof(writer->error_msg),
                     "Failed to add column '%s' to schema", column_names[i]);
            for (int j = 0; j <= i; j++) free(writer->column_names[j]);
            free(writer->column_names);
            free(writer->column_types);
            carquet_schema_free(writer->schema);
            free(writer);
            return NULL;
        }
    }

    /* Create writer with ZSTD compression (better compression ratio) */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.compression_level = 3;  /* ZSTD level 1-22, default 3 for good balance */

    writer->writer = carquet_writer_create(filename, writer->schema, &opts, &err);
    if (!writer->writer) {
        snprintf(writer->error_msg, sizeof(writer->error_msg),
                 "Failed to create parquet writer: %s", err.message);
        for (int i = 0; i < num_columns; i++) free(writer->column_names[i]);
        free(writer->column_names);
        free(writer->column_types);
        carquet_schema_free(writer->schema);
        free(writer);
        return NULL;
    }

    /* Initialize row buffer */
    writer->buffer_capacity = PARQUET_WRITE_BUFFER_SIZE;
    writer->buffer_size = 0;
    writer->row_buffer = (char***)calloc(writer->buffer_capacity, sizeof(char**));
    if (!writer->row_buffer) {
        carquet_writer_abort(writer->writer);
        for (int i = 0; i < num_columns; i++) free(writer->column_names[i]);
        free(writer->column_names);
        free(writer->column_types);
        carquet_schema_free(writer->schema);
        free(writer);
        return NULL;
    }

    writer->rows_written = 0;

    return writer;
}

static void free_row_buffer(parquet_writer_t* writer) {
    if (!writer || !writer->row_buffer) return;

    for (int i = 0; i < writer->buffer_size; i++) {
        if (writer->row_buffer[i]) {
            for (int j = 0; j < writer->num_columns; j++) {
                free(writer->row_buffer[i][j]);
            }
            free(writer->row_buffer[i]);
        }
    }
    writer->buffer_size = 0;
}

parquet_status_t parquet_writer_flush(parquet_writer_t* writer) {
    if (!writer || writer->buffer_size == 0) return PARQUET_OK;

    /* Write buffered data column by column */
    for (int col = 0; col < writer->num_columns; col++) {
        carquet_status_t status;

        switch (writer->column_types[col]) {
            case PARQUET_COL_DOUBLE: {
                /* Convert strings to doubles and write */
                double* values = (double*)malloc(writer->buffer_size * sizeof(double));
                if (!values) return PARQUET_ERROR_MEMORY;

                for (int row = 0; row < writer->buffer_size; row++) {
                    const char* val = writer->row_buffer[row][col];
                    if (val && val[0] != '\0') {
                        values[row] = strtod(val, NULL);
                    } else {
                        values[row] = 0.0;  /* Default for empty values */
                    }
                }

                status = carquet_writer_write_batch(
                    writer->writer, col, values, writer->buffer_size, NULL, NULL);
                free(values);
                break;
            }

            case PARQUET_COL_INT64: {
                /* Convert strings to int64 and write */
                int64_t* values = (int64_t*)malloc(writer->buffer_size * sizeof(int64_t));
                if (!values) return PARQUET_ERROR_MEMORY;

                for (int row = 0; row < writer->buffer_size; row++) {
                    const char* val = writer->row_buffer[row][col];
                    if (val && val[0] != '\0') {
                        values[row] = strtoll(val, NULL, 10);
                    } else {
                        values[row] = 0;  /* Default for empty values */
                    }
                }

                status = carquet_writer_write_batch(
                    writer->writer, col, values, writer->buffer_size, NULL, NULL);
                free(values);
                break;
            }

            case PARQUET_COL_STRING:
            default: {
                /* Write as byte arrays (strings) */
                carquet_byte_array_t* values = (carquet_byte_array_t*)malloc(
                    writer->buffer_size * sizeof(carquet_byte_array_t));
                if (!values) return PARQUET_ERROR_MEMORY;

                for (int row = 0; row < writer->buffer_size; row++) {
                    const char* val = writer->row_buffer[row][col];
                    if (val && val[0] != '\0') {
                        values[row].data = (uint8_t*)val;
                        values[row].length = (int32_t)strlen(val);
                    } else {
                        values[row].data = (uint8_t*)"";
                        values[row].length = 0;
                    }
                }

                status = carquet_writer_write_batch(
                    writer->writer, col, values, writer->buffer_size, NULL, NULL);
                free(values);
                break;
            }
        }

        if (status != CARQUET_OK) {
            snprintf(writer->error_msg, sizeof(writer->error_msg),
                     "Failed to write column %d", col);
            return PARQUET_ERROR_IO;
        }
    }

    writer->rows_written += writer->buffer_size;
    free_row_buffer(writer);

    return PARQUET_OK;
}

void parquet_writer_free(parquet_writer_t* writer) {
    if (!writer) return;

    /* Flush remaining data */
    parquet_writer_flush(writer);

    /* Close writer (writes footer) */
    if (writer->writer) {
        (void)carquet_writer_close(writer->writer);
    }

    /* Cleanup */
    free_row_buffer(writer);
    free(writer->row_buffer);

    if (writer->column_names) {
        for (int i = 0; i < writer->num_columns; i++) {
            free(writer->column_names[i]);
        }
        free(writer->column_names);
    }

    free(writer->column_types);

    if (writer->schema) {
        carquet_schema_free(writer->schema);
    }

    free(writer);
}

parquet_status_t parquet_writer_write_row(parquet_writer_t* writer,
                                           const char** fields,
                                           int num_fields) {
    if (!writer || !fields) return PARQUET_ERROR_INVALID;
    if (num_fields != writer->num_columns) return PARQUET_ERROR_INVALID;

    /* Flush if buffer is full */
    if (writer->buffer_size >= writer->buffer_capacity) {
        parquet_status_t status = parquet_writer_flush(writer);
        if (status != PARQUET_OK) return status;
    }

    /* Add row to buffer */
    char** row = (char**)calloc(writer->num_columns, sizeof(char*));
    if (!row) return PARQUET_ERROR_MEMORY;

    for (int i = 0; i < writer->num_columns; i++) {
        row[i] = fields[i] ? strdup(fields[i]) : strdup("");
        if (!row[i]) {
            for (int j = 0; j < i; j++) free(row[j]);
            free(row);
            return PARQUET_ERROR_MEMORY;
        }
    }

    writer->row_buffer[writer->buffer_size++] = row;

    return PARQUET_OK;
}

const char* parquet_writer_get_error(const parquet_writer_t* writer) {
    return writer ? writer->error_msg : "Invalid writer";
}

/* ============================================================================
 * Streaming Parquet Reader Implementation
 * ============================================================================ */

parquet_stream_reader_t* parquet_stream_reader_create(const char* filename) {
    if (!filename) return NULL;

    parquet_stream_reader_t* reader = (parquet_stream_reader_t*)calloc(1, sizeof(parquet_stream_reader_t));
    if (!reader) return NULL;

    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Open parquet file */
    reader->reader = carquet_reader_open(filename, NULL, &err);
    if (!reader->reader) {
        snprintf(reader->error_msg, sizeof(reader->error_msg),
                 "Failed to open parquet file: %s", err.message);
        free(reader);
        return NULL;
    }

    /* Get metadata */
    reader->total_rows = carquet_reader_num_rows(reader->reader);
    reader->num_columns = carquet_reader_num_columns(reader->reader);

    /* Cache column names */
    const carquet_schema_t* schema = carquet_reader_schema(reader->reader);
    reader->column_names = (char**)calloc(reader->num_columns, sizeof(char*));
    if (!reader->column_names) {
        carquet_reader_close(reader->reader);
        free(reader);
        return NULL;
    }

    for (int i = 0; i < reader->num_columns; i++) {
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, i + 1);
        if (node) {
            reader->column_names[i] = strdup(carquet_schema_node_name(node));
        } else {
            char buf[32];
            snprintf(buf, sizeof(buf), "column_%d", i);
            reader->column_names[i] = strdup(buf);
        }
    }

    /* Allocate field buffers for value conversion */
    reader->field_buffers = (char**)calloc(reader->num_columns, sizeof(char*));
    reader->field_buffer_sizes = (size_t*)calloc(reader->num_columns, sizeof(size_t));
    if (!reader->field_buffers || !reader->field_buffer_sizes) {
        for (int i = 0; i < reader->num_columns; i++) free(reader->column_names[i]);
        free(reader->column_names);
        free(reader->field_buffers);
        free(reader->field_buffer_sizes);
        carquet_reader_close(reader->reader);
        free(reader);
        return NULL;
    }

    /* Pre-allocate field buffers */
    for (int i = 0; i < reader->num_columns; i++) {
        reader->field_buffer_sizes[i] = 256;
        reader->field_buffers[i] = (char*)malloc(reader->field_buffer_sizes[i]);
        if (!reader->field_buffers[i]) {
            for (int j = 0; j < i; j++) free(reader->field_buffers[j]);
            for (int j = 0; j < reader->num_columns; j++) free(reader->column_names[j]);
            free(reader->column_names);
            free(reader->field_buffers);
            free(reader->field_buffer_sizes);
            carquet_reader_close(reader->reader);
            free(reader);
            return NULL;
        }
    }

    /* Create batch reader */
    carquet_batch_reader_config_t config;
    carquet_batch_reader_config_init(&config);
    config.batch_size = 10000;

    reader->batch_reader = carquet_batch_reader_create(reader->reader, &config, &err);
    if (!reader->batch_reader) {
        snprintf(reader->error_msg, sizeof(reader->error_msg),
                 "Failed to create batch reader: %s", err.message);
        for (int i = 0; i < reader->num_columns; i++) {
            free(reader->column_names[i]);
            free(reader->field_buffers[i]);
        }
        free(reader->column_names);
        free(reader->field_buffers);
        free(reader->field_buffer_sizes);
        carquet_reader_close(reader->reader);
        free(reader);
        return NULL;
    }

    reader->current_batch = NULL;
    reader->current_row_in_batch = 0;
    reader->batch_num_rows = 0;
    reader->rows_read = 0;
    reader->at_eof = false;

    return reader;
}

void parquet_stream_reader_free(parquet_stream_reader_t* reader) {
    if (!reader) return;

    if (reader->current_batch) {
        carquet_row_batch_free(reader->current_batch);
    }
    if (reader->batch_reader) {
        carquet_batch_reader_free(reader->batch_reader);
    }
    if (reader->reader) {
        carquet_reader_close(reader->reader);
    }
    if (reader->column_names) {
        for (int i = 0; i < reader->num_columns; i++) {
            free(reader->column_names[i]);
        }
        free(reader->column_names);
    }
    if (reader->field_buffers) {
        for (int i = 0; i < reader->num_columns; i++) {
            free(reader->field_buffers[i]);
        }
        free(reader->field_buffers);
    }
    free(reader->field_buffer_sizes);
    free(reader);
}

/* Helper to convert value to string in buffer */
static const char* value_to_string_buf(const void* data, int64_t row_idx,
                                        carquet_physical_type_t type,
                                        const uint8_t* null_bitmap,
                                        char** buf, size_t* buf_size) {
    /* Check for null */
    if (null_bitmap) {
        bool is_null = !(null_bitmap[row_idx / 8] & (1 << (row_idx % 8)));
        if (is_null) {
            (*buf)[0] = '\0';
            return *buf;
        }
    }

    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN: {
            const uint8_t* vals = (const uint8_t*)data;
            snprintf(*buf, *buf_size, "%s", vals[row_idx] ? "true" : "false");
            break;
        }
        case CARQUET_PHYSICAL_INT32: {
            const int32_t* vals = (const int32_t*)data;
            snprintf(*buf, *buf_size, "%d", vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_INT64: {
            const int64_t* vals = (const int64_t*)data;
            snprintf(*buf, *buf_size, "%lld", (long long)vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_FLOAT: {
            const float* vals = (const float*)data;
            snprintf(*buf, *buf_size, "%.6g", vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_DOUBLE: {
            const double* vals = (const double*)data;
            snprintf(*buf, *buf_size, "%.6g", vals[row_idx]);
            break;
        }
        case CARQUET_PHYSICAL_BYTE_ARRAY: {
            const carquet_byte_array_t* vals = (const carquet_byte_array_t*)data;
            size_t needed = vals[row_idx].length + 1;
            if (needed > *buf_size) {
                char* new_buf = (char*)realloc(*buf, needed);
                if (new_buf) {
                    *buf = new_buf;
                    *buf_size = needed;
                }
            }
            if (vals[row_idx].length > 0 && vals[row_idx].data) {
                memcpy(*buf, vals[row_idx].data, vals[row_idx].length);
                (*buf)[vals[row_idx].length] = '\0';
            } else {
                (*buf)[0] = '\0';
            }
            break;
        }
        default:
            (*buf)[0] = '\0';
            break;
    }

    return *buf;
}

int parquet_stream_reader_next(parquet_stream_reader_t* reader,
                                const char** fields,
                                int max_fields) {
    if (!reader || !fields) return -1;
    if (reader->at_eof) return -1;

    /* Load next batch if needed */
    if (!reader->current_batch ||
        reader->current_row_in_batch >= reader->batch_num_rows) {

        if (reader->current_batch) {
            carquet_row_batch_free(reader->current_batch);
            reader->current_batch = NULL;
        }

        carquet_status_t status = carquet_batch_reader_next(reader->batch_reader,
                                                             &reader->current_batch);
        if (status != CARQUET_OK || !reader->current_batch) {
            reader->at_eof = true;
            return -1;
        }
        reader->current_row_in_batch = 0;
        reader->batch_num_rows = carquet_row_batch_num_rows(reader->current_batch);
    }

    /* Get schema for type info */
    const carquet_schema_t* schema = carquet_reader_schema(reader->reader);

    /* Extract values for current row */
    int num_fields = (max_fields < reader->num_columns) ? max_fields : reader->num_columns;

    for (int col = 0; col < num_fields; col++) {
        const void* data;
        const uint8_t* null_bitmap;
        int64_t num_values;

        carquet_status_t status = carquet_row_batch_column(reader->current_batch, col,
                                                            &data, &null_bitmap, &num_values);
        if (status != CARQUET_OK) {
            reader->field_buffers[col][0] = '\0';
            fields[col] = reader->field_buffers[col];
            continue;
        }

        /* Get column type */
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, col + 1);
        carquet_physical_type_t ptype = node ? carquet_schema_node_physical_type(node)
                                              : CARQUET_PHYSICAL_BYTE_ARRAY;

        fields[col] = value_to_string_buf(data, reader->current_row_in_batch,
                                           ptype, null_bitmap,
                                           &reader->field_buffers[col],
                                           &reader->field_buffer_sizes[col]);
    }

    reader->current_row_in_batch++;
    reader->rows_read++;

    return num_fields;
}

bool parquet_stream_reader_at_end(const parquet_stream_reader_t* reader) {
    return reader ? reader->at_eof : true;
}

int64_t parquet_stream_reader_num_rows(const parquet_stream_reader_t* reader) {
    return reader ? reader->total_rows : 0;
}

int32_t parquet_stream_reader_num_columns(const parquet_stream_reader_t* reader) {
    return reader ? reader->num_columns : 0;
}

const char* parquet_stream_reader_column_name(const parquet_stream_reader_t* reader, int idx) {
    if (!reader || idx < 0 || idx >= reader->num_columns) return NULL;
    return reader->column_names[idx];
}

int parquet_stream_reader_find_column(const parquet_stream_reader_t* reader, const char* name) {
    if (!reader || !name) return -1;
    for (int i = 0; i < reader->num_columns; i++) {
        if (reader->column_names[i] && strcmp(reader->column_names[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

const char* parquet_stream_reader_get_error(const parquet_stream_reader_t* reader) {
    return reader ? reader->error_msg : "Invalid reader";
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

bool is_parquet_file(const char* filename) {
    if (!filename) return false;

    size_t len = strlen(filename);
    if (len < 4) return false;

    /* Check for .parquet extension */
    if (len >= 8 && strcasecmp(filename + len - 8, ".parquet") == 0) {
        return true;
    }

    /* Check for .pq extension */
    if (len >= 3 && strcasecmp(filename + len - 3, ".pq") == 0) {
        return true;
    }

    return false;
}

int64_t parquet_count_rows(const char* filename) {
    if (!filename) return 0;

    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Open reader just to get row count */
    carquet_reader_t* reader = carquet_reader_open(filename, NULL, &err);
    if (!reader) {
        return 0;
    }

    int64_t num_rows = carquet_reader_num_rows(reader);
    carquet_reader_close(reader);

    return num_rows;
}
