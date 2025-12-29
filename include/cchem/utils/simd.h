/**
 * @file simd.h
 * @brief Portable SIMD utilities for ARM NEON and x86 SSE/AVX
 *
 * Provides common vector operations that compile to efficient SIMD on both platforms.
 * Falls back to scalar loops when SIMD is unavailable.
 */

#ifndef CCHEM_SIMD_H
#define CCHEM_SIMD_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

/* Detect SIMD capability */
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #define CCHEM_SIMD_NEON 1
    #include <arm_neon.h>
#elif defined(__SSE2__)
    #define CCHEM_SIMD_SSE2 1
    #include <emmintrin.h>
    #if defined(__SSE4_1__)
        #define CCHEM_SIMD_SSE41 1
        #include <smmintrin.h>
    #endif
    #if defined(__AVX__)
        #define CCHEM_SIMD_AVX 1
        #include <immintrin.h>
        #if defined(__FMA__)
            #define CCHEM_SIMD_FMA 1
        #endif
    #endif
#endif

/* ============================================================================
 * Horizontal sum of double array (most common operation)
 * ============================================================================ */

static inline double simd_sum_double(const double* arr, int n) {
#if defined(CCHEM_SIMD_AVX)
    /* AVX: process 4 doubles at a time */
    __m256d sum_vec = _mm256_setzero_pd();
    int i = 0;
    for (; i <= n - 4; i += 4) {
        __m256d v = _mm256_loadu_pd(&arr[i]);
        sum_vec = _mm256_add_pd(sum_vec, v);
    }
    /* Horizontal sum of 4 doubles */
    __m128d lo = _mm256_castpd256_pd128(sum_vec);
    __m128d hi = _mm256_extractf128_pd(sum_vec, 1);
    __m128d sum128 = _mm_add_pd(lo, hi);
    sum128 = _mm_hadd_pd(sum128, sum128);
    double sum = _mm_cvtsd_f64(sum128);
    /* Handle remainder */
    for (; i < n; i++) sum += arr[i];
    return sum;

#elif defined(CCHEM_SIMD_SSE2)
    /* SSE2: process 2 doubles at a time */
    __m128d sum_vec = _mm_setzero_pd();
    int i = 0;
    for (; i <= n - 2; i += 2) {
        __m128d v = _mm_loadu_pd(&arr[i]);
        sum_vec = _mm_add_pd(sum_vec, v);
    }
    /* Horizontal sum */
    __m128d shuf = _mm_shuffle_pd(sum_vec, sum_vec, 1);
    sum_vec = _mm_add_pd(sum_vec, shuf);
    double sum = _mm_cvtsd_f64(sum_vec);
    for (; i < n; i++) sum += arr[i];
    return sum;

#elif defined(CCHEM_SIMD_NEON)
    /* NEON: process 2 doubles at a time */
    float64x2_t sum_vec = vdupq_n_f64(0.0);
    int i = 0;
    for (; i <= n - 2; i += 2) {
        float64x2_t v = vld1q_f64(&arr[i]);
        sum_vec = vaddq_f64(sum_vec, v);
    }
    double sum = vgetq_lane_f64(sum_vec, 0) + vgetq_lane_f64(sum_vec, 1);
    for (; i < n; i++) sum += arr[i];
    return sum;

#else
    /* Scalar fallback */
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += arr[i];
    return sum;
#endif
}

/* ============================================================================
 * Sum of absolute values
 * ============================================================================ */

static inline double simd_sum_abs_double(const double* arr, int n) {
#if defined(CCHEM_SIMD_AVX)
    __m256d sum_vec = _mm256_setzero_pd();
    __m256d sign_mask = _mm256_set1_pd(-0.0);
    int i = 0;
    for (; i <= n - 4; i += 4) {
        __m256d v = _mm256_loadu_pd(&arr[i]);
        v = _mm256_andnot_pd(sign_mask, v);  /* abs */
        sum_vec = _mm256_add_pd(sum_vec, v);
    }
    __m128d lo = _mm256_castpd256_pd128(sum_vec);
    __m128d hi = _mm256_extractf128_pd(sum_vec, 1);
    __m128d sum128 = _mm_add_pd(lo, hi);
    sum128 = _mm_hadd_pd(sum128, sum128);
    double sum = _mm_cvtsd_f64(sum128);
    for (; i < n; i++) sum += fabs(arr[i]);
    return sum;

#elif defined(CCHEM_SIMD_SSE2)
    __m128d sum_vec = _mm_setzero_pd();
    __m128d sign_mask = _mm_set1_pd(-0.0);
    int i = 0;
    for (; i <= n - 2; i += 2) {
        __m128d v = _mm_loadu_pd(&arr[i]);
        v = _mm_andnot_pd(sign_mask, v);
        sum_vec = _mm_add_pd(sum_vec, v);
    }
    __m128d shuf = _mm_shuffle_pd(sum_vec, sum_vec, 1);
    sum_vec = _mm_add_pd(sum_vec, shuf);
    double sum = _mm_cvtsd_f64(sum_vec);
    for (; i < n; i++) sum += fabs(arr[i]);
    return sum;

#elif defined(CCHEM_SIMD_NEON)
    float64x2_t sum_vec = vdupq_n_f64(0.0);
    int i = 0;
    for (; i <= n - 2; i += 2) {
        float64x2_t v = vld1q_f64(&arr[i]);
        v = vabsq_f64(v);
        sum_vec = vaddq_f64(sum_vec, v);
    }
    double sum = vgetq_lane_f64(sum_vec, 0) + vgetq_lane_f64(sum_vec, 1);
    for (; i < n; i++) sum += fabs(arr[i]);
    return sum;

#else
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += fabs(arr[i]);
    return sum;
#endif
}

/* ============================================================================
 * Sum of squares
 * ============================================================================ */

static inline double simd_sum_sq_double(const double* arr, int n) {
#if defined(CCHEM_SIMD_AVX)
    __m256d sum_vec = _mm256_setzero_pd();
    int i = 0;
    for (; i <= n - 4; i += 4) {
        __m256d v = _mm256_loadu_pd(&arr[i]);
#if defined(CCHEM_SIMD_FMA)
        sum_vec = _mm256_fmadd_pd(v, v, sum_vec);  /* sum += v*v (FMA) */
#else
        __m256d sq = _mm256_mul_pd(v, v);          /* sq = v*v */
        sum_vec = _mm256_add_pd(sum_vec, sq);      /* sum += sq */
#endif
    }
    __m128d lo = _mm256_castpd256_pd128(sum_vec);
    __m128d hi = _mm256_extractf128_pd(sum_vec, 1);
    __m128d sum128 = _mm_add_pd(lo, hi);
    sum128 = _mm_hadd_pd(sum128, sum128);
    double sum = _mm_cvtsd_f64(sum128);
    for (; i < n; i++) sum += arr[i] * arr[i];
    return sum;

#elif defined(CCHEM_SIMD_SSE2)
    __m128d sum_vec = _mm_setzero_pd();
    int i = 0;
    for (; i <= n - 2; i += 2) {
        __m128d v = _mm_loadu_pd(&arr[i]);
        __m128d sq = _mm_mul_pd(v, v);
        sum_vec = _mm_add_pd(sum_vec, sq);
    }
    __m128d shuf = _mm_shuffle_pd(sum_vec, sum_vec, 1);
    sum_vec = _mm_add_pd(sum_vec, shuf);
    double sum = _mm_cvtsd_f64(sum_vec);
    for (; i < n; i++) sum += arr[i] * arr[i];
    return sum;

#elif defined(CCHEM_SIMD_NEON)
    float64x2_t sum_vec = vdupq_n_f64(0.0);
    int i = 0;
    for (; i <= n - 2; i += 2) {
        float64x2_t v = vld1q_f64(&arr[i]);
        sum_vec = vfmaq_f64(sum_vec, v, v);
    }
    double sum = vgetq_lane_f64(sum_vec, 0) + vgetq_lane_f64(sum_vec, 1);
    for (; i < n; i++) sum += arr[i] * arr[i];
    return sum;

#else
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += arr[i] * arr[i];
    return sum;
#endif
}

/* ============================================================================
 * Max value
 * ============================================================================ */

static inline double simd_max_double(const double* arr, int n) {
    if (n <= 0) return 0.0;

#if defined(CCHEM_SIMD_AVX)
    __m256d max_vec = _mm256_set1_pd(arr[0]);
    int i = 1;
    for (; i <= n - 4; i += 4) {
        __m256d v = _mm256_loadu_pd(&arr[i]);
        max_vec = _mm256_max_pd(max_vec, v);
    }
    /* Reduce */
    __m128d lo = _mm256_castpd256_pd128(max_vec);
    __m128d hi = _mm256_extractf128_pd(max_vec, 1);
    __m128d max128 = _mm_max_pd(lo, hi);
    __m128d shuf = _mm_shuffle_pd(max128, max128, 1);
    max128 = _mm_max_pd(max128, shuf);
    double max_val = _mm_cvtsd_f64(max128);
    for (; i < n; i++) if (arr[i] > max_val) max_val = arr[i];
    return max_val;

#elif defined(CCHEM_SIMD_SSE2)
    __m128d max_vec = _mm_set1_pd(arr[0]);
    int i = 1;
    for (; i <= n - 2; i += 2) {
        __m128d v = _mm_loadu_pd(&arr[i]);
        max_vec = _mm_max_pd(max_vec, v);
    }
    __m128d shuf = _mm_shuffle_pd(max_vec, max_vec, 1);
    max_vec = _mm_max_pd(max_vec, shuf);
    double max_val = _mm_cvtsd_f64(max_vec);
    for (; i < n; i++) if (arr[i] > max_val) max_val = arr[i];
    return max_val;

#elif defined(CCHEM_SIMD_NEON)
    float64x2_t max_vec = vdupq_n_f64(arr[0]);
    int i = 1;
    for (; i <= n - 2; i += 2) {
        float64x2_t v = vld1q_f64(&arr[i]);
        max_vec = vmaxq_f64(max_vec, v);
    }
    double max_val = vmaxvq_f64(max_vec);
    for (; i < n; i++) if (arr[i] > max_val) max_val = arr[i];
    return max_val;

#else
    double max_val = arr[0];
    for (int i = 1; i < n; i++) if (arr[i] > max_val) max_val = arr[i];
    return max_val;
#endif
}

/* ============================================================================
 * Min value
 * ============================================================================ */

static inline double simd_min_double(const double* arr, int n) {
    if (n <= 0) return 0.0;

#if defined(CCHEM_SIMD_AVX)
    __m256d min_vec = _mm256_set1_pd(arr[0]);
    int i = 1;
    for (; i <= n - 4; i += 4) {
        __m256d v = _mm256_loadu_pd(&arr[i]);
        min_vec = _mm256_min_pd(min_vec, v);
    }
    __m128d lo = _mm256_castpd256_pd128(min_vec);
    __m128d hi = _mm256_extractf128_pd(min_vec, 1);
    __m128d min128 = _mm_min_pd(lo, hi);
    __m128d shuf = _mm_shuffle_pd(min128, min128, 1);
    min128 = _mm_min_pd(min128, shuf);
    double min_val = _mm_cvtsd_f64(min128);
    for (; i < n; i++) if (arr[i] < min_val) min_val = arr[i];
    return min_val;

#elif defined(CCHEM_SIMD_SSE2)
    __m128d min_vec = _mm_set1_pd(arr[0]);
    int i = 1;
    for (; i <= n - 2; i += 2) {
        __m128d v = _mm_loadu_pd(&arr[i]);
        min_vec = _mm_min_pd(min_vec, v);
    }
    __m128d shuf = _mm_shuffle_pd(min_vec, min_vec, 1);
    min_vec = _mm_min_pd(min_vec, shuf);
    double min_val = _mm_cvtsd_f64(min_vec);
    for (; i < n; i++) if (arr[i] < min_val) min_val = arr[i];
    return min_val;

#elif defined(CCHEM_SIMD_NEON)
    float64x2_t min_vec = vdupq_n_f64(arr[0]);
    int i = 1;
    for (; i <= n - 2; i += 2) {
        float64x2_t v = vld1q_f64(&arr[i]);
        min_vec = vminq_f64(min_vec, v);
    }
    double min_val = vminvq_f64(min_vec);
    for (; i < n; i++) if (arr[i] < min_val) min_val = arr[i];
    return min_val;

#else
    double min_val = arr[0];
    for (int i = 1; i < n; i++) if (arr[i] < min_val) min_val = arr[i];
    return min_val;
#endif
}

/* ============================================================================
 * Dot product
 * ============================================================================ */

static inline double simd_dot_double(const double* a, const double* b, int n) {
#if defined(CCHEM_SIMD_AVX)
    __m256d sum_vec = _mm256_setzero_pd();
    int i = 0;
    for (; i <= n - 4; i += 4) {
        __m256d va = _mm256_loadu_pd(&a[i]);
        __m256d vb = _mm256_loadu_pd(&b[i]);
        sum_vec = _mm256_fmadd_pd(va, vb, sum_vec);
    }
    __m128d lo = _mm256_castpd256_pd128(sum_vec);
    __m128d hi = _mm256_extractf128_pd(sum_vec, 1);
    __m128d sum128 = _mm_add_pd(lo, hi);
    sum128 = _mm_hadd_pd(sum128, sum128);
    double sum = _mm_cvtsd_f64(sum128);
    for (; i < n; i++) sum += a[i] * b[i];
    return sum;

#elif defined(CCHEM_SIMD_SSE2)
    __m128d sum_vec = _mm_setzero_pd();
    int i = 0;
    for (; i <= n - 2; i += 2) {
        __m128d va = _mm_loadu_pd(&a[i]);
        __m128d vb = _mm_loadu_pd(&b[i]);
        __m128d prod = _mm_mul_pd(va, vb);
        sum_vec = _mm_add_pd(sum_vec, prod);
    }
    __m128d shuf = _mm_shuffle_pd(sum_vec, sum_vec, 1);
    sum_vec = _mm_add_pd(sum_vec, shuf);
    double sum = _mm_cvtsd_f64(sum_vec);
    for (; i < n; i++) sum += a[i] * b[i];
    return sum;

#elif defined(CCHEM_SIMD_NEON)
    float64x2_t sum_vec = vdupq_n_f64(0.0);
    int i = 0;
    for (; i <= n - 2; i += 2) {
        float64x2_t va = vld1q_f64(&a[i]);
        float64x2_t vb = vld1q_f64(&b[i]);
        sum_vec = vfmaq_f64(sum_vec, va, vb);
    }
    double sum = vgetq_lane_f64(sum_vec, 0) + vgetq_lane_f64(sum_vec, 1);
    for (; i < n; i++) sum += a[i] * b[i];
    return sum;

#else
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
#endif
}

/* ============================================================================
 * Population count for fingerprint operations
 * ============================================================================ */

static inline int simd_popcount_bytes(const uint8_t* data, int n) {
    int count = 0;

#if defined(__POPCNT__) || defined(__ARM_FEATURE_CRC32)
    /* Use hardware popcount */
    int i = 0;
    #if defined(__LP64__) || defined(_WIN64)
    for (; i <= n - 8; i += 8) {
        uint64_t v;
        memcpy(&v, &data[i], 8);
        count += __builtin_popcountll(v);
    }
    #endif
    for (; i <= n - 4; i += 4) {
        uint32_t v;
        memcpy(&v, &data[i], 4);
        count += __builtin_popcount(v);
    }
    for (; i < n; i++) {
        count += __builtin_popcount(data[i]);
    }
#else
    /* Lookup table fallback */
    static const uint8_t popcount_table[256] = {
        0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
    };
    for (int i = 0; i < n; i++) {
        count += popcount_table[data[i]];
    }
#endif
    return count;
}

#endif /* CCHEM_SIMD_H */
