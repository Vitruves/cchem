/**
 * @file commands.h
 * @brief CLI command declarations and shared utilities
 */

#ifndef CCHEM_COMMANDS_H
#define CCHEM_COMMANDS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "cchem/compat.h"

#ifdef _WIN32
#include <io.h>
#include <getopt.h>  /* from vcpkg */
#include <windows.h>
#define isatty _isatty
#define fileno _fileno
#define usleep(us) Sleep((us) / 1000)
#include <pthread.h>  /* from vcpkg pthreads */
#else
#include <getopt.h>
#include <unistd.h>
#include <pthread.h>
#endif

#include "cchem/cchem.h"
#include "cchem/utils/descriptors.h"
#include "cchem/utils/csv.h"
#include "cchem/utils/progress.h"
#include "cchem/utils/threading.h"
#include "cchem/utils/memory.h"
#ifdef HAVE_CAIRO
#include "cchem/depictor/depictor.h"
#endif
#include "cchem/splitter/splitter.h"

#ifdef HAVE_PARQUET
#include "cchem/utils/parquet.h"
#endif

#define VERSION "1.1.1"

/* ============================================================================
 * Command Entry Points
 * ============================================================================ */

/**
 * @brief Execute the canonicalize command
 * @param argc Argument count (after subcommand)
 * @param argv Argument vector (argv[0] is "canonicalize")
 * @return Exit code (0 on success)
 */
int cmd_canonicalize(int argc, char* argv[]);

/**
 * @brief Execute the validate command
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit code (0 on success)
 */
int cmd_validate(int argc, char* argv[]);

/**
 * @brief Execute the compute command
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit code (0 on success)
 */
int cmd_compute(int argc, char* argv[]);

#ifdef HAVE_CAIRO
/**
 * @brief Execute the depict command
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit code (0 on success)
 */
int cmd_depict(int argc, char* argv[]);
#endif

/**
 * @brief Execute the split command
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit code (0 on success)
 */
int cmd_split(int argc, char* argv[]);

/* ============================================================================
 * Usage Printers
 * ============================================================================ */

void print_canonicalize_usage(const char* prog_name);
void print_validate_usage(const char* prog_name);
void print_compute_usage(const char* prog_name);
#ifdef HAVE_CAIRO
void print_depict_usage(const char* prog_name);
#endif
void print_split_usage(const char* prog_name);

/* ============================================================================
 * Shared Utilities for Compute Pipeline
 * ============================================================================ */

/* Fixed-width field for pre-allocated buffer (24 bytes covers most descriptor values) */
#define DESC_VALUE_WIDTH 24

/**
 * @brief Fast integer to string conversion (avoids snprintf overhead)
 * @param val Integer value to convert
 * @param buf Output buffer
 * @param buf_size Buffer size
 */
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

/**
 * @brief Fast double to string with 6 significant figures (like %.6g)
 * @param val Double value to convert
 * @param buf Output buffer
 * @param buf_size Buffer size
 */
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

#endif /* CCHEM_COMMANDS_H */
