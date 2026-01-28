/**
 * @file csv.c
 * @brief CSV file reading and writing implementation
 */

#include "cchem/compat.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
/* Windows file operations */
#define open _open
#define close _close
#define fstat _fstat
#define stat _stat
#define O_RDONLY _O_RDONLY
/* Windows memory mapping */
static void* win_mmap(size_t length, int fd) {
    HANDLE hFile = (HANDLE)_get_osfhandle(fd);
    HANDLE hMapping = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!hMapping) return NULL;
    void* ptr = MapViewOfFile(hMapping, FILE_MAP_READ, 0, 0, length);
    CloseHandle(hMapping);
    return ptr;
}
static void win_munmap(void* addr, size_t length) {
    (void)length;
    UnmapViewOfFile(addr);
}
#define mmap(addr, len, prot, flags, fd, off) win_mmap(len, fd)
#define munmap(addr, len) win_munmap(addr, len)
#define MAP_FAILED NULL
#define PROT_READ 0
#define MAP_PRIVATE 0
#else
#include <strings.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#endif

#include "cchem/utils/csv.h"

#define INITIAL_ROW_CAPACITY 16
#define INITIAL_LINE_BUFFER_SIZE 4096

csv_row_t* csv_row_create(void) {
    csv_row_t* row = (csv_row_t*)calloc(1, sizeof(csv_row_t));
    if (!row) return NULL;

    row->capacity = INITIAL_ROW_CAPACITY;
    row->fields = (char**)calloc(row->capacity, sizeof(char*));
    if (!row->fields) {
        free(row);
        return NULL;
    }

    row->num_fields = 0;
    return row;
}

void csv_row_free(csv_row_t* row) {
    if (!row) return;

    for (int i = 0; i < row->num_fields; i++) {
        if (row->fields[i]) free(row->fields[i]);
    }

    if (row->fields) free(row->fields);
    free(row);
}

static void csv_row_clear(csv_row_t* row) {
    if (!row) return;

    for (int i = 0; i < row->num_fields; i++) {
        if (row->fields[i]) {
            free(row->fields[i]);
            row->fields[i] = NULL;
        }
    }
    row->num_fields = 0;
}

csv_status_t csv_row_add_field(csv_row_t* row, const char* field) {
    if (!row) return CSV_ERROR_INVALID;

    if (row->num_fields >= row->capacity) {
        int new_cap = row->capacity * 2;
        char** new_fields = (char**)realloc(row->fields, new_cap * sizeof(char*));
        if (!new_fields) return CSV_ERROR_MEMORY;
        row->fields = new_fields;
        row->capacity = new_cap;
    }

    row->fields[row->num_fields] = field ? strdup(field) : strdup("");
    if (!row->fields[row->num_fields]) return CSV_ERROR_MEMORY;

    row->num_fields++;
    return CSV_OK;
}

const char* csv_row_get_field(const csv_row_t* row, int idx) {
    if (!row || idx < 0 || idx >= row->num_fields) return NULL;
    return row->fields[idx];
}

const char* csv_row_get_field_by_name(const csv_row_t* row,
                                       const csv_row_t* header,
                                       const char* name) {
    int idx = csv_find_column(header, name);
    if (idx < 0) return NULL;
    return csv_row_get_field(row, idx);
}

int csv_find_column(const csv_row_t* header, const char* name) {
    if (!header || !name) return -1;

    for (int i = 0; i < header->num_fields; i++) {
        if (header->fields[i] && strcmp(header->fields[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

csv_document_t* csv_document_create(void) {
    csv_document_t* doc = (csv_document_t*)calloc(1, sizeof(csv_document_t));
    if (!doc) return NULL;

    doc->rows_capacity = 64;
    doc->rows = (csv_row_t*)calloc(doc->rows_capacity, sizeof(csv_row_t));
    if (!doc->rows) {
        free(doc);
        return NULL;
    }

    doc->delimiter = ',';
    doc->quote_char = '"';
    doc->has_header = false;
    doc->num_rows = 0;

    memset(&doc->header, 0, sizeof(csv_row_t));

    return doc;
}

void csv_document_free(csv_document_t* doc) {
    if (!doc) return;

    /* Free header */
    for (int i = 0; i < doc->header.num_fields; i++) {
        if (doc->header.fields[i]) free(doc->header.fields[i]);
    }
    if (doc->header.fields) free(doc->header.fields);

    /* Free rows */
    for (int r = 0; r < doc->num_rows; r++) {
        for (int i = 0; i < doc->rows[r].num_fields; i++) {
            if (doc->rows[r].fields[i]) free(doc->rows[r].fields[i]);
        }
        if (doc->rows[r].fields) free(doc->rows[r].fields);
    }
    if (doc->rows) free(doc->rows);

    free(doc);
}

csv_reader_t* csv_reader_create(const char* filename) {
    if (!filename) return NULL;

    csv_reader_t* reader = (csv_reader_t*)calloc(1, sizeof(csv_reader_t));
    if (!reader) return NULL;

    reader->file = fopen(filename, "r");
    if (!reader->file) {
        snprintf(reader->error_msg, sizeof(reader->error_msg),
                "Failed to open file: %s", filename);
        free(reader);
        return NULL;
    }

    /* Use large buffer to reduce syscalls */
    setvbuf(reader->file, NULL, _IOFBF, 262144);  /* 256KB buffer */

    reader->line_buffer_size = INITIAL_LINE_BUFFER_SIZE;
    reader->line_buffer = (char*)calloc(reader->line_buffer_size, sizeof(char));
    if (!reader->line_buffer) {
        fclose(reader->file);
        free(reader);
        return NULL;
    }

    reader->delimiter = ',';
    reader->quote_char = '"';
    reader->current_line = 0;
    reader->at_eof = false;
    reader->error_msg[0] = '\0';

    return reader;
}

void csv_reader_free(csv_reader_t* reader) {
    if (!reader) return;

    if (reader->file) fclose(reader->file);
    if (reader->line_buffer) free(reader->line_buffer);

    free(reader);
}

/* Parse a single line into fields */
static csv_status_t parse_line(csv_reader_t* reader, csv_row_t* row) {
    csv_row_clear(row);

    if (!row->fields) {
        row->capacity = INITIAL_ROW_CAPACITY;
        row->fields = (char**)calloc(row->capacity, sizeof(char*));
        if (!row->fields) return CSV_ERROR_MEMORY;
    }

    const char* line = reader->line_buffer;
    size_t len = strlen(line);

    /* Remove trailing newline/carriage return */
    while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
        len--;
    }

    if (len == 0) {
        return CSV_OK;
    }

    /* Parse fields */
    size_t pos = 0;
    char* field_buf = (char*)malloc(len + 1);
    if (!field_buf) return CSV_ERROR_MEMORY;

    while (pos <= len) {
        size_t field_len = 0;

        if (pos < len && line[pos] == reader->quote_char) {
            /* Quoted field */
            pos++;  /* Skip opening quote */

            while (pos < len) {
                if (line[pos] == reader->quote_char) {
                    if (pos + 1 < len && line[pos + 1] == reader->quote_char) {
                        /* Escaped quote */
                        field_buf[field_len++] = reader->quote_char;
                        pos += 2;
                    } else {
                        /* End of quoted field */
                        pos++;  /* Skip closing quote */
                        break;
                    }
                } else {
                    field_buf[field_len++] = line[pos++];
                }
            }

            /* Skip to delimiter or end */
            while (pos < len && line[pos] != reader->delimiter) {
                pos++;
            }
        } else {
            /* Unquoted field */
            while (pos < len && line[pos] != reader->delimiter) {
                field_buf[field_len++] = line[pos++];
            }
        }

        field_buf[field_len] = '\0';

        /* Add field */
        csv_status_t status = csv_row_add_field(row, field_buf);
        if (status != CSV_OK) {
            free(field_buf);
            return status;
        }

        /* Skip delimiter */
        if (pos < len && line[pos] == reader->delimiter) {
            pos++;
            /* Handle trailing delimiter */
            if (pos == len) {
                csv_row_add_field(row, "");
            }
        } else {
            break;
        }
    }

    free(field_buf);
    return CSV_OK;
}

csv_status_t csv_reader_next(csv_reader_t* reader, csv_row_t* row) {
    if (!reader || !row) return CSV_ERROR_INVALID;

    if (reader->at_eof || !reader->file) {
        return CSV_ERROR_IO;
    }

    /* Read line */
    if (fgets(reader->line_buffer, (int)reader->line_buffer_size, reader->file) == NULL) {
        reader->at_eof = true;
        return CSV_ERROR_IO;
    }

    reader->current_line++;

    /* Handle very long lines */
    while (strchr(reader->line_buffer, '\n') == NULL && !feof(reader->file)) {
        size_t current_len = strlen(reader->line_buffer);
        size_t new_size = reader->line_buffer_size * 2;
        char* new_buf = (char*)realloc(reader->line_buffer, new_size);
        if (!new_buf) return CSV_ERROR_MEMORY;

        reader->line_buffer = new_buf;
        reader->line_buffer_size = new_size;

        if (fgets(reader->line_buffer + current_len,
                  (int)(new_size - current_len), reader->file) == NULL) {
            break;
        }
    }

    return parse_line(reader, row);
}

bool csv_reader_at_end(const csv_reader_t* reader) {
    if (!reader) return true;
    return reader->at_eof;
}

int csv_reader_get_line(const csv_reader_t* reader) {
    if (!reader) return 0;
    return reader->current_line;
}

const char* csv_reader_get_error(const csv_reader_t* reader) {
    if (!reader) return "NULL reader";
    return reader->error_msg;
}

csv_writer_t* csv_writer_create(const char* filename) {
    if (!filename) return NULL;

    csv_writer_t* writer = (csv_writer_t*)calloc(1, sizeof(csv_writer_t));
    if (!writer) return NULL;

    writer->file = fopen(filename, "w");
    if (!writer->file) {
        snprintf(writer->error_msg, sizeof(writer->error_msg),
                "Failed to create file: %s", filename);
        free(writer);
        return NULL;
    }

    /* Use large buffer to reduce syscalls */
    setvbuf(writer->file, NULL, _IOFBF, 262144);  /* 256KB buffer */

    writer->delimiter = ',';
    writer->quote_char = '"';
    writer->needs_quotes = true;
    writer->error_msg[0] = '\0';

    return writer;
}

void csv_writer_free(csv_writer_t* writer) {
    if (!writer) return;

    if (writer->file) fclose(writer->file);
    free(writer);
}

const char* csv_sanitize_numeric(const char* field) {
    /* Handle NULL or empty strings */
    if (!field || field[0] == '\0') {
        return "0";
    }

    /* Skip leading whitespace */
    const char* p = field;
    while (*p == ' ' || *p == '\t') p++;

    /* Check for empty after whitespace */
    if (*p == '\0') {
        return "0";
    }

    /* Check for valid numeric format:
     * [+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?
     * Also accepts: nan, inf, -inf
     */

    /* Handle special values - fast inline checks instead of strcasecmp */
    /* Check for nan/NaN/-nan (3-4 chars) */
    if ((p[0] == 'n' || p[0] == 'N') && (p[1] == 'a' || p[1] == 'A') &&
        (p[2] == 'n' || p[2] == 'N') && p[3] == '\0') {
        return "0";
    }
    if (p[0] == '-' && (p[1] == 'n' || p[1] == 'N') && (p[2] == 'a' || p[2] == 'A') &&
        (p[3] == 'n' || p[3] == 'N') && p[4] == '\0') {
        return "0";
    }
    /* Check for inf/Inf/+inf/-inf (3-4 chars) */
    const char* ip = p;
    if (*ip == '+' || *ip == '-') ip++;
    if ((ip[0] == 'i' || ip[0] == 'I') && (ip[1] == 'n' || ip[1] == 'N') &&
        (ip[2] == 'f' || ip[2] == 'F')) {
        if (ip[3] == '\0') return "0";  /* inf */
        /* Check for infinity */
        if ((ip[3] == 'i' || ip[3] == 'I') && (ip[4] == 'n' || ip[4] == 'N') &&
            (ip[5] == 'i' || ip[5] == 'I') && (ip[6] == 't' || ip[6] == 'T') &&
            (ip[7] == 'y' || ip[7] == 'Y') && ip[8] == '\0') {
            return "0";  /* infinity */
        }
    }

    /* Check for valid number format */
    bool has_digit = false;

    /* Optional sign */
    if (*p == '+' || *p == '-') p++;

    /* Integer part */
    while (*p >= '0' && *p <= '9') {
        has_digit = true;
        p++;
    }

    /* Optional decimal point and fractional part */
    if (*p == '.') {
        p++;
        while (*p >= '0' && *p <= '9') {
            has_digit = true;
            p++;
        }
    }

    /* Optional exponent */
    if (*p == 'e' || *p == 'E') {
        p++;
        /* Optional sign in exponent */
        if (*p == '+' || *p == '-') p++;
        /* Must have at least one digit in exponent */
        bool exp_digit = false;
        while (*p >= '0' && *p <= '9') {
            exp_digit = true;
            p++;
        }
        if (!exp_digit) {
            return "0";  /* Invalid: 'e' without digits */
        }
    }

    /* Skip trailing whitespace */
    while (*p == ' ' || *p == '\t') p++;

    /* Must be at end of string and have at least one digit */
    if (*p != '\0' || !has_digit) {
        return "0";
    }

    /* Valid numeric value - return original */
    return field;
}

/* Fast buffered row writing - avoid per-character fputc overhead */
csv_status_t csv_writer_write_row(csv_writer_t* writer, const csv_row_t* row) {
    if (!writer || !row || !writer->file) return CSV_ERROR_INVALID;

    /* Use thread-local static buffer for speed (avoid malloc per row) */
    static __thread char line_buf[65536];
    char* buf_end = line_buf + sizeof(line_buf) - 2;
    char* p = line_buf;

    for (int i = 0; i < row->num_fields; i++) {
        if (i > 0 && p < buf_end) {
            *p++ = writer->delimiter;
        }

        const char* field = row->fields[i] ? row->fields[i] : "";

        /* Check if quoting needed while copying (single pass) */
        bool needs_quote = false;
        for (const char* f = field; *f; f++) {
            if (*f == writer->delimiter || *f == '"' || *f == '\n' || *f == '\r') {
                needs_quote = true;
                break;
            }
        }

        if (needs_quote) {
            if (p < buf_end) *p++ = writer->quote_char;
            for (const char* f = field; *f && p < buf_end; f++) {
                if (*f == writer->quote_char && p < buf_end - 1) {
                    *p++ = writer->quote_char;
                }
                *p++ = *f;
            }
            if (p < buf_end) *p++ = writer->quote_char;
        } else {
            /* Fast copy for unquoted fields */
            while (*field && p < buf_end) {
                *p++ = *field++;
            }
        }
    }

    *p++ = '\n';

    /* Single write call for entire line */
    fwrite(line_buf, 1, p - line_buf, writer->file);
    return CSV_OK;
}

csv_status_t csv_writer_write_fields(csv_writer_t* writer,
                                     const char** fields, int num_fields) {
    if (!writer || !fields || !writer->file) return CSV_ERROR_INVALID;

    /* Direct write - avoid csv_row_t wrapper overhead */
    static __thread char line_buf[65536];
    char* buf_end = line_buf + sizeof(line_buf) - 2;
    char* p = line_buf;

    for (int i = 0; i < num_fields; i++) {
        if (i > 0 && p < buf_end) {
            *p++ = writer->delimiter;
        }

        const char* field = fields[i] ? fields[i] : "";

        /* Fast path: check first char for common numeric cases (no quoting needed) */
        char first = *field;
        bool maybe_numeric = (first == '-' || first == '+' || (first >= '0' && first <= '9') || first == '.');

        if (maybe_numeric) {
            /* Fast copy - numeric fields never need quoting */
            while (*field && p < buf_end) {
                *p++ = *field++;
            }
        } else {
            /* Check if quoting needed */
            bool needs_quote = false;
            for (const char* f = field; *f; f++) {
                if (*f == writer->delimiter || *f == '"' || *f == '\n' || *f == '\r') {
                    needs_quote = true;
                    break;
                }
            }

            if (needs_quote) {
                if (p < buf_end) *p++ = writer->quote_char;
                for (const char* f = field; *f && p < buf_end; f++) {
                    if (*f == writer->quote_char && p < buf_end - 1) {
                        *p++ = writer->quote_char;
                    }
                    *p++ = *f;
                }
                if (p < buf_end) *p++ = writer->quote_char;
            } else {
                while (*field && p < buf_end) {
                    *p++ = *field++;
                }
            }
        }
    }

    *p++ = '\n';
    fwrite(line_buf, 1, p - line_buf, writer->file);
    return CSV_OK;
}

const char* csv_writer_get_error(const csv_writer_t* writer) {
    if (!writer) return "NULL writer";
    return writer->error_msg;
}

int csv_count_rows(const char* filename) {
    if (!filename) return 0;

    FILE* f = fopen(filename, "r");
    if (!f) return 0;

    int count = 0;
    char buf[4096];

    while (fgets(buf, sizeof(buf), f)) {
        /* Skip empty lines */
        char* p = buf;
        while (*p && isspace(*p)) p++;
        if (*p) count++;
    }

    fclose(f);
    return count;
}

csv_status_t csv_document_read(csv_document_t* doc, const char* filename,
                               bool first_row_header) {
    if (!doc || !filename) return CSV_ERROR_INVALID;

    csv_reader_t* reader = csv_reader_create(filename);
    if (!reader) return CSV_ERROR_IO;

    csv_row_t row;
    memset(&row, 0, sizeof(row));

    bool first = true;
    csv_status_t status;

    while ((status = csv_reader_next(reader, &row)) == CSV_OK) {
        if (first && first_row_header) {
            /* Copy to header */
            doc->header.capacity = row.capacity;
            doc->header.num_fields = row.num_fields;
            doc->header.fields = row.fields;
            row.fields = NULL;
            row.num_fields = 0;
            row.capacity = 0;
            doc->has_header = true;
            first = false;
            continue;
        }
        first = false;

        /* Ensure capacity */
        if (doc->num_rows >= doc->rows_capacity) {
            int new_cap = doc->rows_capacity * 2;
            csv_row_t* new_rows = (csv_row_t*)realloc(doc->rows, new_cap * sizeof(csv_row_t));
            if (!new_rows) {
                csv_reader_free(reader);
                return CSV_ERROR_MEMORY;
            }
            doc->rows = new_rows;
            doc->rows_capacity = new_cap;
        }

        /* Move row data */
        doc->rows[doc->num_rows].capacity = row.capacity;
        doc->rows[doc->num_rows].num_fields = row.num_fields;
        doc->rows[doc->num_rows].fields = row.fields;
        doc->num_rows++;

        row.fields = NULL;
        row.num_fields = 0;
        row.capacity = 0;
    }

    csv_reader_free(reader);
    return CSV_OK;
}

csv_status_t csv_document_write(const csv_document_t* doc, const char* filename) {
    if (!doc || !filename) return CSV_ERROR_INVALID;

    csv_writer_t* writer = csv_writer_create(filename);
    if (!writer) return CSV_ERROR_IO;

    csv_status_t status;

    /* Write header if present */
    if (doc->has_header && doc->header.num_fields > 0) {
        status = csv_writer_write_row(writer, &doc->header);
        if (status != CSV_OK) {
            csv_writer_free(writer);
            return status;
        }
    }

    /* Write rows */
    for (int i = 0; i < doc->num_rows; i++) {
        status = csv_writer_write_row(writer, &doc->rows[i]);
        if (status != CSV_OK) {
            csv_writer_free(writer);
            return status;
        }
    }

    csv_writer_free(writer);
    return CSV_OK;
}

/* ============================================================================
 * Memory-Mapped CSV Reader Implementation
 * ============================================================================ */

csv_mmap_reader_t* csv_mmap_reader_create(const char* filename) {
    if (!filename) return NULL;

    csv_mmap_reader_t* reader = (csv_mmap_reader_t*)calloc(1, sizeof(csv_mmap_reader_t));
    if (!reader) return NULL;

    reader->fd = open(filename, O_RDONLY);
    if (reader->fd < 0) {
        snprintf(reader->error_msg, sizeof(reader->error_msg),
                "Failed to open file: %s", filename);
        free(reader);
        return NULL;
    }

    /* Get file size */
    struct stat st;
    if (fstat(reader->fd, &st) < 0) {
        snprintf(reader->error_msg, sizeof(reader->error_msg),
                "Failed to stat file: %s", filename);
        close(reader->fd);
        free(reader);
        return NULL;
    }

    reader->size = st.st_size;

    if (reader->size == 0) {
        /* Empty file - valid but nothing to read */
        reader->data = NULL;
        reader->at_eof = true;
    } else {
        /* Memory map the file */
        reader->data = mmap(NULL, reader->size, PROT_READ, MAP_PRIVATE, reader->fd, 0);
        if (reader->data == MAP_FAILED) {
            snprintf(reader->error_msg, sizeof(reader->error_msg),
                    "Failed to mmap file: %s", filename);
            close(reader->fd);
            free(reader);
            return NULL;
        }

        /* Advise kernel we'll read sequentially */
        madvise(reader->data, reader->size, MADV_SEQUENTIAL);
    }

    reader->pos = 0;
    reader->delimiter = ',';
    reader->quote_char = '"';
    reader->current_line = 0;
    reader->at_eof = (reader->size == 0);
    reader->error_msg[0] = '\0';

    return reader;
}

void csv_mmap_reader_free(csv_mmap_reader_t* reader) {
    if (!reader) return;

    if (reader->data && reader->data != MAP_FAILED) {
        munmap(reader->data, reader->size);
    }
    if (reader->fd >= 0) {
        close(reader->fd);
    }
    free(reader);
}

const char* csv_mmap_reader_next_line(csv_mmap_reader_t* reader, size_t* len) {
    if (!reader || reader->at_eof || reader->pos >= reader->size) {
        reader->at_eof = true;
        if (len) *len = 0;
        return NULL;
    }

    const char* start = reader->data + reader->pos;
    const char* end = reader->data + reader->size;
    const char* p = start;

    /* Find end of line */
    while (p < end && *p != '\n') {
        p++;
    }

    /* Calculate line length (excluding \n but may include \r) */
    size_t line_len = p - start;

    /* Move position past newline */
    if (p < end) {
        reader->pos = (p - reader->data) + 1;  /* Skip \n */
    } else {
        reader->pos = reader->size;
        reader->at_eof = true;
    }

    reader->current_line++;

    /* Strip trailing \r if present */
    if (line_len > 0 && start[line_len - 1] == '\r') {
        line_len--;
    }

    if (len) *len = line_len;
    return start;
}

/* Thread-local buffer for quoted field parsing */
static __thread char tl_field_buf[65536];
static __thread size_t tl_field_pos = 0;

int csv_mmap_parse_line(const char* line, size_t len, char delimiter,
                        const char** fields, int max_fields) {
    if (!line || !fields || max_fields <= 0) return -1;

    if (len == 0) {
        return 0;  /* Empty line = no fields */
    }

    int num_fields = 0;
    size_t pos = 0;
    tl_field_pos = 0;  /* Reset thread-local buffer */

    while (pos <= len && num_fields < max_fields) {
        if (pos < len && line[pos] == '"') {
            /* Quoted field - need to copy to handle escaped quotes */
            pos++;  /* Skip opening quote */
            size_t field_start = tl_field_pos;

            while (pos < len) {
                if (line[pos] == '"') {
                    if (pos + 1 < len && line[pos + 1] == '"') {
                        /* Escaped quote */
                        tl_field_buf[tl_field_pos++] = '"';
                        pos += 2;
                    } else {
                        /* End of quoted field */
                        pos++;
                        break;
                    }
                } else {
                    tl_field_buf[tl_field_pos++] = line[pos++];
                }
            }

            /* Null-terminate the field */
            tl_field_buf[tl_field_pos++] = '\0';
            fields[num_fields++] = tl_field_buf + field_start;

            /* Skip to delimiter or end */
            while (pos < len && line[pos] != delimiter) pos++;

        } else {
            /* Unquoted field - zero-copy pointer directly into mmap'd data */
            const char* field_start = line + pos;
            size_t field_len = 0;

            while (pos < len && line[pos] != delimiter) {
                pos++;
                field_len++;
            }

            /* For unquoted fields, we need to null-terminate, so copy to buffer */
            if (tl_field_pos + field_len + 1 < sizeof(tl_field_buf)) {
                memcpy(tl_field_buf + tl_field_pos, field_start, field_len);
                tl_field_buf[tl_field_pos + field_len] = '\0';
                fields[num_fields++] = tl_field_buf + tl_field_pos;
                tl_field_pos += field_len + 1;
            }
        }

        /* Skip delimiter */
        if (pos < len && line[pos] == delimiter) {
            pos++;
            /* Handle trailing delimiter = empty field */
            if (pos == len && num_fields < max_fields) {
                tl_field_buf[tl_field_pos] = '\0';
                fields[num_fields++] = tl_field_buf + tl_field_pos;
                tl_field_pos++;
            }
        } else {
            break;
        }
    }

    return num_fields;
}

bool csv_mmap_reader_at_end(const csv_mmap_reader_t* reader) {
    if (!reader) return true;
    return reader->at_eof;
}

int csv_mmap_reader_get_line(const csv_mmap_reader_t* reader) {
    if (!reader) return 0;
    return reader->current_line;
}

const char* csv_mmap_reader_get_error(const csv_mmap_reader_t* reader) {
    if (!reader) return "NULL reader";
    return reader->error_msg;
}

/* ============================================================================
 * Bulk CSV Writer Implementation
 * ============================================================================ */

FILE* csv_writer_get_file(csv_writer_t* writer) {
    return writer ? writer->file : NULL;
}

csv_bulk_writer_t* csv_bulk_writer_create(int num_rows, size_t avg_line_size) {
    if (num_rows <= 0) return NULL;

    csv_bulk_writer_t* bulk = (csv_bulk_writer_t*)calloc(1, sizeof(csv_bulk_writer_t));
    if (!bulk) return NULL;

    bulk->lines = (csv_line_buffer_t*)calloc(num_rows, sizeof(csv_line_buffer_t));
    if (!bulk->lines) {
        free(bulk);
        return NULL;
    }

    bulk->num_lines = num_rows;
    bulk->delimiter = ',';
    bulk->quote_char = '"';

    /* Pre-allocate each line buffer - use calloc to zero-initialize */
    size_t initial_cap = (avg_line_size > 0) ? avg_line_size : 1024;
    for (int i = 0; i < num_rows; i++) {
        bulk->lines[i].data = (char*)calloc(initial_cap, 1);
        if (!bulk->lines[i].data) {
            /* Cleanup on failure */
            for (int j = 0; j < i; j++) {
                free(bulk->lines[j].data);
            }
            free(bulk->lines);
            free(bulk);
            return NULL;
        }
        bulk->lines[i].capacity = initial_cap;
        bulk->lines[i].len = 0;
    }

    return bulk;
}

void csv_bulk_writer_free(csv_bulk_writer_t* bulk) {
    if (!bulk) return;

    if (bulk->lines) {
        for (int i = 0; i < bulk->num_lines; i++) {
            if (bulk->lines[i].data) {
                free(bulk->lines[i].data);
            }
        }
        free(bulk->lines);
    }
    free(bulk);
}

/* Thread-safe row formatting into pre-allocated buffer */
csv_status_t csv_bulk_format_row(csv_bulk_writer_t* bulk, int row_idx,
                                  const char** fields, int num_fields) {
    if (!bulk || row_idx < 0 || row_idx >= bulk->num_lines || !fields) {
        return CSV_ERROR_INVALID;
    }

    csv_line_buffer_t* line = &bulk->lines[row_idx];

    /* Calculate required size */
    size_t needed = 0;
    for (int i = 0; i < num_fields; i++) {
        const char* field = fields[i] ? fields[i] : "";
        size_t flen = strlen(field);
        needed += flen + 3;  /* field + possible quotes + delimiter */

        /* Check if quoting needed */
        for (const char* p = field; *p; p++) {
            if (*p == bulk->quote_char) needed++;  /* Escaped quotes */
        }
    }
    needed += 2;  /* newline + null terminator */

    /* Ensure capacity */
    if (needed > line->capacity) {
        size_t new_cap = needed * 2;
        char* new_data = (char*)realloc(line->data, new_cap);
        if (!new_data) return CSV_ERROR_MEMORY;
        /* Zero-initialize the new space */
        memset(new_data + line->capacity, 0, new_cap - line->capacity);
        line->data = new_data;
        line->capacity = new_cap;
    }

    /* Format row */
    char* p = line->data;
    char* end = line->data + line->capacity - 2;

    for (int i = 0; i < num_fields; i++) {
        if (i > 0 && p < end) {
            *p++ = bulk->delimiter;
        }

        const char* field = fields[i] ? fields[i] : "";

        /* Fast path: check if numeric (no quoting needed) */
        char first = *field;
        bool maybe_numeric = (first == '-' || first == '+' ||
                              (first >= '0' && first <= '9') ||
                              first == '.' || first == '\0');

        if (maybe_numeric) {
            /* Fast copy for numeric/empty fields */
            while (*field && p < end) {
                *p++ = *field++;
            }
        } else {
            /* Check if quoting needed */
            bool needs_quote = false;
            for (const char* f = field; *f; f++) {
                if (*f == bulk->delimiter || *f == '"' || *f == '\n' || *f == '\r') {
                    needs_quote = true;
                    break;
                }
            }

            if (needs_quote) {
                if (p < end) *p++ = bulk->quote_char;
                for (const char* f = field; *f && p < end; f++) {
                    if (*f == bulk->quote_char && p < end - 1) {
                        *p++ = bulk->quote_char;
                    }
                    *p++ = *f;
                }
                if (p < end) *p++ = bulk->quote_char;
            } else {
                while (*field && p < end) {
                    *p++ = *field++;
                }
            }
        }
    }

    *p++ = '\n';
    line->len = p - line->data;

    return CSV_OK;
}

/* Bulk write all formatted lines to file */
csv_status_t csv_bulk_write_all(csv_bulk_writer_t* bulk, FILE* file) {
    if (!bulk || !file) return CSV_ERROR_INVALID;

    /* Calculate total size for potential single write */
    size_t total_size = 0;
    for (int i = 0; i < bulk->num_lines; i++) {
        total_size += bulk->lines[i].len;
    }

    /* For large datasets, use chunked writing to avoid memory issues */
    const size_t CHUNK_SIZE = 64 * 1024 * 1024;  /* 64MB chunks */

    if (total_size < CHUNK_SIZE) {
        /* Small enough - allocate single buffer and write once */
        char* buffer = (char*)malloc(total_size);
        if (buffer) {
            char* p = buffer;
            for (int i = 0; i < bulk->num_lines; i++) {
                memcpy(p, bulk->lines[i].data, bulk->lines[i].len);
                p += bulk->lines[i].len;
            }
            size_t written = fwrite(buffer, 1, total_size, file);
            free(buffer);
            if (written != total_size) return CSV_ERROR_IO;
            return CSV_OK;
        }
        /* Fall through to chunked write if malloc fails */
    }

    /* Chunked writing for large files or if buffer allocation failed */
    char* chunk = (char*)malloc(CHUNK_SIZE);
    if (!chunk) {
        /* Last resort: write line by line */
        for (int i = 0; i < bulk->num_lines; i++) {
            if (fwrite(bulk->lines[i].data, 1, bulk->lines[i].len, file) != bulk->lines[i].len) {
                return CSV_ERROR_IO;
            }
        }
        return CSV_OK;
    }

    size_t chunk_pos = 0;
    for (int i = 0; i < bulk->num_lines; i++) {
        size_t line_len = bulk->lines[i].len;

        /* If line doesn't fit in remaining chunk, flush */
        if (chunk_pos + line_len > CHUNK_SIZE) {
            if (chunk_pos > 0) {
                if (fwrite(chunk, 1, chunk_pos, file) != chunk_pos) {
                    free(chunk);
                    return CSV_ERROR_IO;
                }
                chunk_pos = 0;
            }

            /* If single line > chunk size, write directly */
            if (line_len > CHUNK_SIZE) {
                if (fwrite(bulk->lines[i].data, 1, line_len, file) != line_len) {
                    free(chunk);
                    return CSV_ERROR_IO;
                }
                continue;
            }
        }

        memcpy(chunk + chunk_pos, bulk->lines[i].data, line_len);
        chunk_pos += line_len;
    }

    /* Flush remaining chunk */
    if (chunk_pos > 0) {
        if (fwrite(chunk, 1, chunk_pos, file) != chunk_pos) {
            free(chunk);
            return CSV_ERROR_IO;
        }
    }

    free(chunk);
    return CSV_OK;
}

/* Bulk write with ninja-style progress */
csv_status_t csv_bulk_write_all_progress(csv_bulk_writer_t* bulk, FILE* file) {
    if (!bulk || !file) return CSV_ERROR_INVALID;

    int total = bulk->num_lines;
    int update_interval = total > 1000 ? total / 100 : 10;  /* ~100 updates max */
    if (update_interval < 1) update_interval = 1;

    for (int i = 0; i < total; i++) {
        if (fwrite(bulk->lines[i].data, 1, bulk->lines[i].len, file) != bulk->lines[i].len) {
            fprintf(stderr, "\r\033[K");
            return CSV_ERROR_IO;
        }

        /* Unified progress: [current/total] Writing */
        if (i % update_interval == 0 || i == total - 1) {
            fprintf(stderr, "\r[%d/%d] Writing\033[K", i + 1, total);
            fflush(stderr);
        }
    }

    fprintf(stderr, "\r");  /* Return to start of line */
    fflush(stderr);
    return CSV_OK;
}

/* Write header line directly to file */
csv_status_t csv_write_header_line(FILE* file, const char** fields,
                                    int num_fields, char delimiter) {
    if (!file || !fields || num_fields <= 0) return CSV_ERROR_INVALID;

    char line_buf[65536];
    char* p = line_buf;
    char* end = line_buf + sizeof(line_buf) - 2;

    for (int i = 0; i < num_fields; i++) {
        if (i > 0 && p < end) {
            *p++ = delimiter;
        }

        const char* field = fields[i] ? fields[i] : "";

        /* Check if quoting needed */
        bool needs_quote = false;
        for (const char* f = field; *f; f++) {
            if (*f == delimiter || *f == '"' || *f == '\n' || *f == '\r') {
                needs_quote = true;
                break;
            }
        }

        if (needs_quote) {
            if (p < end) *p++ = '"';
            for (const char* f = field; *f && p < end; f++) {
                if (*f == '"' && p < end - 1) {
                    *p++ = '"';
                }
                *p++ = *f;
            }
            if (p < end) *p++ = '"';
        } else {
            while (*field && p < end) {
                *p++ = *field++;
            }
        }
    }

    *p++ = '\n';
    size_t len = p - line_buf;

    if (fwrite(line_buf, 1, len, file) != len) {
        return CSV_ERROR_IO;
    }

    return CSV_OK;
}
