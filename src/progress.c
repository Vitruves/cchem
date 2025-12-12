/**
 * @file progress.c
 * @brief Progress bar with rate estimation implementation
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "cchem/progress.h"

#define RATE_SAMPLES 10
#define MIN_UPDATE_INTERVAL_MS 100

/* ANSI color codes */
#define COLOR_RESET   "\033[0m"
#define COLOR_GREEN   "\033[32m"
#define COLOR_CYAN    "\033[36m"
#define COLOR_DIM     "\033[2m"
#define COLOR_BOLD    "\033[1m"

/* Unicode block characters for smooth progress bar */
static const char* BLOCK_FULL = "█";
static const char* BLOCK_EMPTY = "░";
static const char* PARTIAL_BLOCKS[] = {
    " ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"
};

const progress_config_t PROGRESS_CONFIG_DEFAULT = {
    .bar_width = 30,
    .show_percent = true,
    .show_rate = true,
    .show_eta = true,
    .show_count = true,
    .prefix = "",
    .suffix = "",
    .fill_char = '=',
    .empty_char = '-',
    .update_interval_ms = MIN_UPDATE_INTERVAL_MS
};

static double get_time_ms(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000.0 + tv.tv_usec / 1000.0;
}

progress_t* progress_create(size_t total, const progress_config_t* config) {
    progress_t* progress = (progress_t*)calloc(1, sizeof(progress_t));
    if (!progress) return NULL;

    double now_ms = get_time_ms();
    progress->total = total;
    progress->current = 0;
    progress->start_time = time(NULL);
    progress->start_time_ms = now_ms;
    progress->last_update_ms = now_ms;
    progress->rate = 0.0;
    progress->last_count = 0;
    progress->last_sample_time_ms = now_ms;
    progress->finished = false;
    progress->shown = false;

    progress->rate_samples = (double*)calloc(RATE_SAMPLES, sizeof(double));
    if (!progress->rate_samples) {
        free(progress);
        return NULL;
    }
    progress->num_samples = 0;
    progress->sample_idx = 0;

    if (config) {
        progress->config = *config;
    } else {
        progress->config = PROGRESS_CONFIG_DEFAULT;
    }

    return progress;
}

void progress_free(progress_t* progress) {
    if (!progress) return;
    if (progress->rate_samples) free(progress->rate_samples);
    free(progress);
}

void progress_init(progress_t* progress, size_t total) {
    if (!progress) return;

    double now_ms = get_time_ms();
    progress->total = total;
    progress->current = 0;
    progress->start_time = time(NULL);
    progress->start_time_ms = now_ms;
    progress->last_update_ms = now_ms;
    progress->rate = 0.0;
    progress->last_count = 0;
    progress->last_sample_time_ms = now_ms;
    progress->finished = false;
    progress->shown = false;

    if (progress->rate_samples) {
        memset(progress->rate_samples, 0, RATE_SAMPLES * sizeof(double));
    }
    progress->num_samples = 0;
    progress->sample_idx = 0;
}

void progress_format_time(double seconds, char* buf, size_t buf_size) {
    if (seconds < 0) seconds = 0;

    int total_secs = (int)seconds;
    int hours = total_secs / 3600;
    int mins = (total_secs % 3600) / 60;
    int secs = total_secs % 60;

    if (hours > 0) {
        snprintf(buf, buf_size, "%d:%02d:%02d", hours, mins, secs);
    } else {
        snprintf(buf, buf_size, "%d:%02d", mins, secs);
    }
}

void progress_format_rate(double rate, char* buf, size_t buf_size) {
    if (rate < 1.0) {
        snprintf(buf, buf_size, "%.2f it/s", rate);
    } else if (rate < 1000.0) {
        snprintf(buf, buf_size, "%.1f it/s", rate);
    } else if (rate < 1000000.0) {
        snprintf(buf, buf_size, "%.1fk it/s", rate / 1000.0);
    } else {
        snprintf(buf, buf_size, "%.1fM it/s", rate / 1000000.0);
    }
}

static void update_rate(progress_t* progress) {
    double now_ms = get_time_ms();
    double elapsed_since_sample_ms = now_ms - progress->last_sample_time_ms;
    double elapsed_total_ms = now_ms - progress->start_time_ms;

    /* Use shorter interval initially (100ms), then 500ms for stability */
    double sample_interval = (progress->num_samples < 3) ? 100.0 : 500.0;

    if (elapsed_since_sample_ms >= sample_interval && progress->current > progress->last_count) {
        size_t count_diff = progress->current - progress->last_count;
        double sample_rate = (count_diff / elapsed_since_sample_ms) * 1000.0;  /* items/sec */

        /* Add to rolling average */
        progress->rate_samples[progress->sample_idx] = sample_rate;
        progress->sample_idx = (progress->sample_idx + 1) % RATE_SAMPLES;
        if (progress->num_samples < RATE_SAMPLES) progress->num_samples++;

        /* Calculate average */
        double sum = 0;
        for (int i = 0; i < progress->num_samples; i++) {
            sum += progress->rate_samples[i];
        }
        progress->rate = sum / progress->num_samples;

        progress->last_count = progress->current;
        progress->last_sample_time_ms = now_ms;
    }

    /* Fallback: calculate overall rate if no samples yet but enough time elapsed */
    if (progress->rate == 0.0 && elapsed_total_ms >= 50.0 && progress->current > 0) {
        progress->rate = (progress->current / elapsed_total_ms) * 1000.0;
    }
}

static void draw_progress(progress_t* progress) {
    if (!progress) return;

    double percent = progress_get_percent(progress);
    double fill_exact = percent / 100.0 * progress->config.bar_width;
    int filled = (int)fill_exact;
    int partial_idx = (int)((fill_exact - filled) * 8);  /* 8 partial block levels */
    if (filled > progress->config.bar_width) filled = progress->config.bar_width;
    if (partial_idx < 0) partial_idx = 0;
    if (partial_idx > 8) partial_idx = 8;

    /* Clear entire line and move to start */
    fprintf(stderr, "\r\033[K");

    /* Prefix */
    if (progress->config.prefix && progress->config.prefix[0]) {
        fprintf(stderr, "%s ", progress->config.prefix);
    }

    /* Progress bar with Unicode blocks and colors */
    fprintf(stderr, COLOR_DIM "│" COLOR_RESET);

    /* Filled portion (green) */
    fprintf(stderr, COLOR_GREEN);
    for (int i = 0; i < filled; i++) {
        fprintf(stderr, "%s", BLOCK_FULL);
    }

    /* Partial block for smooth transition */
    if (filled < progress->config.bar_width && partial_idx > 0) {
        fprintf(stderr, "%s", PARTIAL_BLOCKS[partial_idx]);
        filled++;  /* Account for partial block position */
    }
    fprintf(stderr, COLOR_RESET);

    /* Empty portion (dim) */
    fprintf(stderr, COLOR_DIM);
    for (int i = filled; i < progress->config.bar_width; i++) {
        fprintf(stderr, "%s", BLOCK_EMPTY);
    }
    fprintf(stderr, COLOR_RESET);

    fprintf(stderr, COLOR_DIM "│" COLOR_RESET);

    /* Percentage */
    if (progress->config.show_percent) {
        fprintf(stderr, " " COLOR_BOLD "%5.1f%%" COLOR_RESET, percent);
    }

    /* Count */
    if (progress->config.show_count) {
        fprintf(stderr, " %zu/%zu", progress->current, progress->total);
    }

    /* Rate */
    if (progress->config.show_rate && progress->rate > 0) {
        char rate_buf[32];
        progress_format_rate(progress->rate, rate_buf, sizeof(rate_buf));
        fprintf(stderr, " " COLOR_CYAN "[%s]" COLOR_RESET, rate_buf);
    }

    /* ETA */
    if (progress->config.show_eta && progress->rate > 0 && !progress->finished) {
        double eta = progress_get_eta(progress);
        if (eta > 0 && eta < 86400 * 365) {  /* Less than a year */
            char eta_buf[32];
            progress_format_time(eta, eta_buf, sizeof(eta_buf));
            fprintf(stderr, " ETA: %s", eta_buf);
        }
    }

    /* Elapsed time for finished */
    if (progress->finished) {
        double elapsed = progress_get_elapsed(progress);
        char elapsed_buf[32];
        progress_format_time(elapsed, elapsed_buf, sizeof(elapsed_buf));
        fprintf(stderr, " " COLOR_GREEN "[%s]" COLOR_RESET, elapsed_buf);
    }

    /* Suffix */
    if (progress->config.suffix && progress->config.suffix[0]) {
        fprintf(stderr, " %s", progress->config.suffix);
    }

    fflush(stderr);
    progress->shown = true;
}

void progress_update(progress_t* progress, size_t current) {
    if (!progress || progress->finished) return;

    progress->current = current;
    if (progress->current > progress->total) {
        progress->current = progress->total;
    }

    /* Rate limiting */
    double now_ms = get_time_ms();
    double elapsed_ms = now_ms - progress->last_update_ms;

    if (elapsed_ms >= progress->config.update_interval_ms || progress->current == progress->total) {
        update_rate(progress);
        draw_progress(progress);
        progress->last_update_ms = now_ms;
    }
}

void progress_increment(progress_t* progress) {
    if (!progress) return;
    progress_update(progress, progress->current + 1);
}

void progress_increment_by(progress_t* progress, size_t n) {
    if (!progress) return;
    progress_update(progress, progress->current + n);
}

void progress_finish(progress_t* progress) {
    if (!progress || progress->finished) return;

    progress->current = progress->total;
    progress->finished = true;

    /* Calculate final rate if not already computed */
    double elapsed_ms = get_time_ms() - progress->start_time_ms;
    if (elapsed_ms > 0 && progress->total > 0) {
        double final_rate = (progress->total / elapsed_ms) * 1000.0;
        /* Use final rate if higher confidence (full run) or no rate computed */
        if (progress->rate == 0.0 || progress->num_samples < 3) {
            progress->rate = final_rate;
        }
    }

    draw_progress(progress);
    fprintf(stderr, "\n");
    fflush(stderr);
}

void progress_reset(progress_t* progress, size_t new_total) {
    if (!progress) return;
    progress_init(progress, new_total);
}

double progress_get_rate(const progress_t* progress) {
    if (!progress) return 0.0;
    return progress->rate;
}

double progress_get_elapsed(const progress_t* progress) {
    if (!progress) return 0.0;

    double now_ms = get_time_ms();
    return (now_ms - progress->start_time_ms) / 1000.0;  /* Convert to seconds */
}

double progress_get_eta(const progress_t* progress) {
    if (!progress || progress->rate <= 0) return 0.0;

    size_t remaining = progress->total - progress->current;
    return (double)remaining / progress->rate;
}

double progress_get_percent(const progress_t* progress) {
    if (!progress || progress->total == 0) return 0.0;
    return (double)progress->current / (double)progress->total * 100.0;
}

void progress_redraw(progress_t* progress) {
    if (progress) {
        draw_progress(progress);
    }
}

void progress_set_prefix(progress_t* progress, const char* prefix) {
    if (progress) {
        progress->config.prefix = prefix ? prefix : "";
    }
}

void progress_set_suffix(progress_t* progress, const char* suffix) {
    if (progress) {
        progress->config.suffix = suffix ? suffix : "";
    }
}
