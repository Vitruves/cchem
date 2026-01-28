/**
 * @file progress.c
 * @brief Modern gradient progress bar with rate estimation
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
struct timeval {
    long tv_sec;
    long tv_usec;
};
static int gettimeofday(struct timeval* tv, void* tz) {
    (void)tz;
    FILETIME ft;
    GetSystemTimeAsFileTime(&ft);
    unsigned long long t = ((unsigned long long)ft.dwHighDateTime << 32) | ft.dwLowDateTime;
    t -= 116444736000000000ULL;  /* Convert to Unix epoch */
    t /= 10;  /* Convert to microseconds */
    tv->tv_sec = (long)(t / 1000000ULL);
    tv->tv_usec = (long)(t % 1000000ULL);
    return 0;
}
#else
#include <sys/time.h>
#endif

#include "cchem/utils/progress.h"

#define RATE_SAMPLES 10
#define MIN_UPDATE_INTERVAL_MS 50



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
        snprintf(buf, buf_size, "%02d:%02d", mins, secs);
    }
}

void progress_format_rate(double rate, char* buf, size_t buf_size) {
    if (rate < 1.0) {
        snprintf(buf, buf_size, "%5.2f/s", rate);
    } else if (rate < 1000.0) {
        snprintf(buf, buf_size, "%5.1f/s", rate);
    } else if (rate < 1000000.0) {
        snprintf(buf, buf_size, "%5.1fK/s", rate / 1000.0);
    } else {
        snprintf(buf, buf_size, "%5.1fM/s", rate / 1000000.0);
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
        double sample_rate = (count_diff / elapsed_since_sample_ms) * 1000.0;

        progress->rate_samples[progress->sample_idx] = sample_rate;
        progress->sample_idx = (progress->sample_idx + 1) % RATE_SAMPLES;
        if (progress->num_samples < RATE_SAMPLES) progress->num_samples++;

        double sum = 0;
        for (int i = 0; i < progress->num_samples; i++) {
            sum += progress->rate_samples[i];
        }
        progress->rate = sum / progress->num_samples;

        progress->last_count = progress->current;
        progress->last_sample_time_ms = now_ms;
    }

    /* Fallback: calculate overall rate if no samples yet */
    if (progress->rate == 0.0 && elapsed_total_ms >= 50.0 && progress->current > 0) {
        progress->rate = (progress->current / elapsed_total_ms) * 1000.0;
    }
}

static void draw_progress(progress_t* progress) {
    if (!progress) return;

    double pct = progress_get_percent(progress);
    double elapsed_s = progress_get_elapsed(progress);
    double eta_s = progress_get_eta(progress);

    char elapsed_str[16], eta_str[16], rate_str[16];
    progress_format_time(elapsed_s, elapsed_str, sizeof(elapsed_str));
    progress_format_time(eta_s, eta_str, sizeof(eta_str));
    progress_format_rate(progress->rate, rate_str, sizeof(rate_str));

    /* Build output in buffer to avoid flicker */
    char output[512];
    int pos = 0;

    /* Unified ninja-style: [current/total] prefix percent rate ETA [elapsed] */
    pos += snprintf(output + pos, sizeof(output) - pos,
                    "\r[%zu/%zu]", progress->current, progress->total);

    /* Prefix */
    if (progress->config.prefix && progress->config.prefix[0]) {
        pos += snprintf(output + pos, sizeof(output) - pos,
                        " %s", progress->config.prefix);
    }

    /* Percentage */
    if (progress->config.show_percent) {
        pos += snprintf(output + pos, sizeof(output) - pos, " %5.1f%%", pct);
    }

    /* Rate */
    if (progress->config.show_rate && progress->rate > 0) {
        pos += snprintf(output + pos, sizeof(output) - pos, " %s", rate_str);
    }

    /* ETA */
    if (progress->config.show_eta && progress->rate > 0 && !progress->finished) {
        if (eta_s > 0 && eta_s < 86400 * 365) {
            pos += snprintf(output + pos, sizeof(output) - pos, " ETA %s", eta_str);
        }
    }

    /* Elapsed time */
    pos += snprintf(output + pos, sizeof(output) - pos, " [%s]", elapsed_str);

    /* Clear to end of line */
    pos += snprintf(output + pos, sizeof(output) - pos, "\033[K");

    /* Single write to stderr */
    fputs(output, stderr);
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

    /* Calculate final rate */
    double elapsed_ms = get_time_ms() - progress->start_time_ms;
    if (elapsed_ms > 0 && progress->total > 0) {
        double final_rate = (progress->total / elapsed_ms) * 1000.0;
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
    return (now_ms - progress->start_time_ms) / 1000.0;
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
