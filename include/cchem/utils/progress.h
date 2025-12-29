/**
 * @file progress.h
 * @brief Progress bar with rate estimation
 */

#ifndef CCHEM_PROGRESS_H
#define CCHEM_PROGRESS_H

#include <stddef.h>
#include <stdbool.h>
#include <time.h>

/* Progress bar configuration */
typedef struct {
    int bar_width;           /* Width of progress bar in characters */
    bool show_percent;       /* Show percentage */
    bool show_rate;          /* Show items per second */
    bool show_eta;           /* Show estimated time remaining */
    bool show_count;         /* Show current/total count */
    const char* prefix;      /* Prefix string */
    const char* suffix;      /* Suffix string */
    char fill_char;          /* Fill character */
    char empty_char;         /* Empty character */
    int update_interval_ms;  /* Minimum time between updates */
} progress_config_t;

/* Default progress configuration */
extern const progress_config_t PROGRESS_CONFIG_DEFAULT;

/* Progress bar state */
typedef struct {
    size_t total;            /* Total items */
    size_t current;          /* Current item */
    time_t start_time;       /* Start time */
    double start_time_ms;    /* Start time in milliseconds for precise timing */
    double last_update_ms;   /* Last update time in milliseconds */

    /* Rate calculation */
    double rate;             /* Items per second */
    double* rate_samples;    /* Rate samples for smoothing */
    int num_samples;
    int sample_idx;
    size_t last_count;       /* Count at last rate sample */
    double last_sample_time_ms; /* Last rate sample time in milliseconds */

    /* Configuration */
    progress_config_t config;

    /* State */
    bool finished;
    bool shown;
} progress_t;

/* Create progress bar */
progress_t* progress_create(size_t total, const progress_config_t* config);

/* Free progress bar */
void progress_free(progress_t* progress);

/* Initialize progress bar with defaults */
void progress_init(progress_t* progress, size_t total);

/* Update progress */
void progress_update(progress_t* progress, size_t current);

/* Increment progress by 1 */
void progress_increment(progress_t* progress);

/* Increment progress by n */
void progress_increment_by(progress_t* progress, size_t n);

/* Finish progress (move to 100% and newline) */
void progress_finish(progress_t* progress);

/* Reset progress bar */
void progress_reset(progress_t* progress, size_t new_total);

/* Get current rate (items/sec) */
double progress_get_rate(const progress_t* progress);

/* Get elapsed time in seconds */
double progress_get_elapsed(const progress_t* progress);

/* Get estimated time remaining in seconds */
double progress_get_eta(const progress_t* progress);

/* Get percentage complete */
double progress_get_percent(const progress_t* progress);

/* Force redraw */
void progress_redraw(progress_t* progress);

/* Set prefix string */
void progress_set_prefix(progress_t* progress, const char* prefix);

/* Set suffix string */
void progress_set_suffix(progress_t* progress, const char* suffix);

/* Format time duration as string */
void progress_format_time(double seconds, char* buf, size_t buf_size);

/* Format rate as string */
void progress_format_rate(double rate, char* buf, size_t buf_size);

#endif /* CCHEM_PROGRESS_H */
