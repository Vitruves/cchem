/**
 * @file moments.c
 * @brief Property Moment Descriptors - Statistical moments of atomic properties
 *
 * Higher-order statistical moments for 7 atomic properties:
 * - Mass (atomic weight)
 * - EN (electronegativity)
 * - Pol (polarizability)
 * - IP (ionization potential)
 * - EA (electron affinity)
 * - VdW (van der Waals volume)
 * - Val (valence electrons)
 *
 * For each property, compute:
 * - Skewness (3rd standardized moment)
 * - Kurtosis (4th standardized moment - 3, excess kurtosis)
 * - Median (middle value)
 * - IQR (interquartile range, Q3 - Q1)
 * - Entropy (Shannon entropy of binned distribution)
 * - Gini (Gini coefficient of inequality)
 *
 * Total: 7 properties × 6 statistics = 42 descriptors
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"

/* Number of statistics per property */
#define NUM_STATS 6

/* Number of properties */
#define NUM_PROPERTIES 7

/* Total moments descriptors */
#define NUM_MOMENTS_DESCRIPTORS (NUM_PROPERTIES * NUM_STATS)

/* Maximum atoms for stack allocation */
#define MAX_ATOMS_STACK 256

/* Number of bins for entropy calculation */
#define NUM_ENTROPY_BINS 10

/* ============================================================================
 * Atomic Property Tables (from standard references)
 * ============================================================================ */

/* Atomic Mass (Da) */
static const double PROP_MASS[] = {
    [ELEM_H]  = 1.008,   [ELEM_C]  = 12.011,  [ELEM_N]  = 14.007,
    [ELEM_O]  = 15.999,  [ELEM_F]  = 18.998,  [ELEM_P]  = 30.974,
    [ELEM_S]  = 32.065,  [ELEM_Cl] = 35.453,  [ELEM_Br] = 79.904,
    [ELEM_I]  = 126.904, [ELEM_Si] = 28.086,  [ELEM_B]  = 10.811,
    [ELEM_Se] = 78.96,   [ELEM_As] = 74.922,  [ELEM_Li] = 6.941,
    [ELEM_Na] = 22.990,  [ELEM_K]  = 39.098,  [ELEM_Mg] = 24.305,
    [ELEM_Ca] = 40.078,  [ELEM_Zn] = 65.380,  [ELEM_Fe] = 55.845,
    [ELEM_Cu] = 63.546,
};

/* Pauling Electronegativity */
static const double PROP_EN[] = {
    [ELEM_H]  = 2.20,  [ELEM_C]  = 2.55,  [ELEM_N]  = 3.04,
    [ELEM_O]  = 3.44,  [ELEM_F]  = 3.98,  [ELEM_P]  = 2.19,
    [ELEM_S]  = 2.58,  [ELEM_Cl] = 3.16,  [ELEM_Br] = 2.96,
    [ELEM_I]  = 2.66,  [ELEM_Si] = 1.90,  [ELEM_B]  = 2.04,
    [ELEM_Se] = 2.55,  [ELEM_As] = 2.18,  [ELEM_Li] = 0.98,
    [ELEM_Na] = 0.93,  [ELEM_K]  = 0.82,  [ELEM_Mg] = 1.31,
    [ELEM_Ca] = 1.00,  [ELEM_Zn] = 1.65,  [ELEM_Fe] = 1.83,
    [ELEM_Cu] = 1.90,
};

/* Atomic Polarizability (Å³) */
static const double PROP_POL[] = {
    [ELEM_H]  = 0.667, [ELEM_C]  = 1.76,  [ELEM_N]  = 1.10,
    [ELEM_O]  = 0.802, [ELEM_F]  = 0.557, [ELEM_P]  = 3.63,
    [ELEM_S]  = 2.90,  [ELEM_Cl] = 2.18,  [ELEM_Br] = 3.05,
    [ELEM_I]  = 4.70,  [ELEM_Si] = 5.38,  [ELEM_B]  = 3.03,
    [ELEM_Se] = 3.77,  [ELEM_As] = 4.31,  [ELEM_Li] = 24.3,
    [ELEM_Na] = 24.1,  [ELEM_K]  = 43.4,  [ELEM_Mg] = 10.6,
    [ELEM_Ca] = 22.8,  [ELEM_Zn] = 5.75,  [ELEM_Fe] = 8.4,
    [ELEM_Cu] = 6.2,
};

/* First Ionization Potential (eV) */
static const double PROP_IP[] = {
    [ELEM_H]  = 13.598, [ELEM_C]  = 11.260, [ELEM_N]  = 14.534,
    [ELEM_O]  = 13.618, [ELEM_F]  = 17.422, [ELEM_P]  = 10.486,
    [ELEM_S]  = 10.360, [ELEM_Cl] = 12.967, [ELEM_Br] = 11.814,
    [ELEM_I]  = 10.451, [ELEM_Si] = 8.151,  [ELEM_B]  = 8.298,
    [ELEM_Se] = 9.752,  [ELEM_As] = 9.815,  [ELEM_Li] = 5.392,
    [ELEM_Na] = 5.139,  [ELEM_K]  = 4.341,  [ELEM_Mg] = 7.646,
    [ELEM_Ca] = 6.113,  [ELEM_Zn] = 9.394,  [ELEM_Fe] = 7.902,
    [ELEM_Cu] = 7.726,
};

/* Electron Affinity (eV) */
static const double PROP_EA[] = {
    [ELEM_H]  = 0.754,  [ELEM_C]  = 1.262,  [ELEM_N]  = -0.07,
    [ELEM_O]  = 1.461,  [ELEM_F]  = 3.401,  [ELEM_P]  = 0.746,
    [ELEM_S]  = 2.077,  [ELEM_Cl] = 3.612,  [ELEM_Br] = 3.364,
    [ELEM_I]  = 3.059,  [ELEM_Si] = 1.389,  [ELEM_B]  = 0.277,
    [ELEM_Se] = 2.021,  [ELEM_As] = 0.804,  [ELEM_Li] = 0.618,
    [ELEM_Na] = 0.548,  [ELEM_K]  = 0.501,  [ELEM_Mg] = 0.0,
    [ELEM_Ca] = 0.024,  [ELEM_Zn] = 0.0,    [ELEM_Fe] = 0.151,
    [ELEM_Cu] = 1.235,
};

/* Van der Waals Volume (Å³) */
static const double PROP_VDW[] = {
    [ELEM_H]  = 7.24,   [ELEM_C]  = 20.58,  [ELEM_N]  = 15.60,
    [ELEM_O]  = 14.71,  [ELEM_F]  = 13.31,  [ELEM_P]  = 24.43,
    [ELEM_S]  = 24.43,  [ELEM_Cl] = 22.45,  [ELEM_Br] = 26.52,
    [ELEM_I]  = 32.52,  [ELEM_Si] = 38.79,  [ELEM_B]  = 17.87,
    [ELEM_Se] = 28.73,  [ELEM_As] = 26.52,  [ELEM_Li] = 22.28,
    [ELEM_Na] = 38.79,  [ELEM_K]  = 76.73,  [ELEM_Mg] = 25.08,
    [ELEM_Ca] = 53.81,  [ELEM_Zn] = 21.19,  [ELEM_Fe] = 25.72,
    [ELEM_Cu] = 25.72,
};

/* Valence Electrons */
static const double PROP_VAL[] = {
    [ELEM_H]  = 1,  [ELEM_C]  = 4,  [ELEM_N]  = 5,
    [ELEM_O]  = 6,  [ELEM_F]  = 7,  [ELEM_P]  = 5,
    [ELEM_S]  = 6,  [ELEM_Cl] = 7,  [ELEM_Br] = 7,
    [ELEM_I]  = 7,  [ELEM_Si] = 4,  [ELEM_B]  = 3,
    [ELEM_Se] = 6,  [ELEM_As] = 5,  [ELEM_Li] = 1,
    [ELEM_Na] = 1,  [ELEM_K]  = 1,  [ELEM_Mg] = 2,
    [ELEM_Ca] = 2,  [ELEM_Zn] = 2,  [ELEM_Fe] = 2,
    [ELEM_Cu] = 1,
};

/* Default values for each property */
static const double PROP_DEFAULTS[] = {
    12.011,  /* Mass - carbon */
    2.55,    /* EN - carbon */
    1.76,    /* Pol - carbon */
    11.26,   /* IP - carbon */
    1.26,    /* EA - carbon */
    20.58,   /* VdW - carbon */
    4.0,     /* Val - carbon */
};

/* Property tables array */
static const double* PROP_TABLES[] = {
    PROP_MASS, PROP_EN, PROP_POL, PROP_IP, PROP_EA, PROP_VDW, PROP_VAL
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_property(element_t elem, int prop_idx) {
    if (elem <= 0 || elem >= 128) return PROP_DEFAULTS[prop_idx];
    double v = PROP_TABLES[prop_idx][elem];
    return (v != 0.0 || elem == ELEM_C) ? v : PROP_DEFAULTS[prop_idx];
}

/* Comparison function for qsort */
static int compare_doubles(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* ============================================================================
 * Statistical Functions
 * ============================================================================ */

/* Compute skewness (3rd standardized moment) */
static double compute_skewness(const double* values, int n, double mean, double std) {
    if (n < 3 || std < 1e-10) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double z = (values[i] - mean) / std;
        sum += z * z * z;
    }
    return sum / n;
}

/* Compute excess kurtosis (4th standardized moment - 3) */
static double compute_kurtosis(const double* values, int n, double mean, double std) {
    if (n < 4 || std < 1e-10) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double z = (values[i] - mean) / std;
        double z2 = z * z;
        sum += z2 * z2;
    }
    return (sum / n) - 3.0;  /* Excess kurtosis */
}

/* Compute median from sorted array */
static double compute_median(const double* sorted, int n) {
    if (n == 0) return 0.0;
    if (n % 2 == 0) {
        return (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
    } else {
        return sorted[n/2];
    }
}

/* Compute quartile (0.25 for Q1, 0.75 for Q3) from sorted array */
static double compute_quartile(const double* sorted, int n, double q) {
    if (n == 0) return 0.0;
    if (n == 1) return sorted[0];

    double idx = q * (n - 1);
    int lo = (int)floor(idx);
    int hi = (int)ceil(idx);

    if (lo == hi || hi >= n) return sorted[lo];

    double frac = idx - lo;
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

/* Compute IQR (interquartile range) */
static double compute_iqr(const double* sorted, int n) {
    return compute_quartile(sorted, n, 0.75) - compute_quartile(sorted, n, 0.25);
}

/* Compute Shannon entropy using binned distribution */
static double compute_entropy(const double* values, int n, double min_val, double max_val) {
    if (n == 0 || max_val - min_val < 1e-10) return 0.0;

    int bins[NUM_ENTROPY_BINS] = {0};
    double bin_width = (max_val - min_val) / NUM_ENTROPY_BINS;

    for (int i = 0; i < n; i++) {
        int bin = (int)((values[i] - min_val) / bin_width);
        if (bin >= NUM_ENTROPY_BINS) bin = NUM_ENTROPY_BINS - 1;
        if (bin < 0) bin = 0;
        bins[bin]++;
    }

    double entropy = 0.0;
    for (int i = 0; i < NUM_ENTROPY_BINS; i++) {
        if (bins[i] > 0) {
            double p = (double)bins[i] / n;
            entropy -= p * log(p);
        }
    }

    /* Normalize to [0, 1] */
    double max_entropy = log((double)NUM_ENTROPY_BINS);
    return max_entropy > 0 ? entropy / max_entropy : 0.0;
}

/* Compute Gini coefficient */
static double compute_gini(const double* sorted, int n) {
    if (n <= 1) return 0.0;

    double sum = 0.0;
    double total = 0.0;

    for (int i = 0; i < n; i++) {
        /* Shift values to be non-negative for Gini */
        double v = sorted[i] - sorted[0] + 1.0;
        sum += (2.0 * (i + 1) - n - 1) * v;
        total += v;
    }

    if (total < 1e-10) return 0.0;
    return sum / (n * total);
}

/* ============================================================================
 * Main Computation
 * ============================================================================ */

int descriptors_compute_moments_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_MOMENTS_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    /* Count heavy atoms */
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    }

    if (n_heavy == 0) return NUM_MOMENTS_DESCRIPTORS;

    /* Allocate arrays - use size_t to avoid overflow warnings */
    size_t n = (size_t)n_heavy;
    size_t arr_size = n * sizeof(double);
    bool heap_alloc = (n_heavy > MAX_ATOMS_STACK);
    double* prop_values;
    double* sorted;

    if (heap_alloc) {
        prop_values = (double*)malloc(arr_size);
        sorted = (double*)malloc(arr_size);
        if (!prop_values || !sorted) {
            free(prop_values);
            free(sorted);
            return -1;
        }
    } else {
        prop_values = (double*)alloca(arr_size);
        sorted = (double*)alloca(arr_size);
    }

    int out_idx = 0;

    /* Process each property */
    for (int p = 0; p < NUM_PROPERTIES; p++) {
        /* Collect property values for heavy atoms */
        int idx = 0;
        double sum = 0.0;
        double min_val = 1e30, max_val = -1e30;

        for (int i = 0; i < mol->num_atoms; i++) {
            if (mol->atoms[i].element != ELEM_H) {
                double v = get_property(mol->atoms[i].element, p);
                prop_values[idx++] = v;
                sum += v;
                if (v < min_val) min_val = v;
                if (v > max_val) max_val = v;
            }
        }

        /* Compute mean and std */
        double mean = sum / (double)n;
        double var_sum = 0.0;
        for (size_t i = 0; i < n; i++) {
            double diff = prop_values[i] - mean;
            var_sum += diff * diff;
        }
        double std = sqrt(var_sum / (double)n);

        /* Create sorted copy for quantile-based stats */
        for (size_t i = 0; i < n; i++) sorted[i] = prop_values[i];
        qsort(sorted, n, sizeof(double), compare_doubles);

        /* Compute statistics */
        values[out_idx++].d = compute_skewness(prop_values, n_heavy, mean, std);
        values[out_idx++].d = compute_kurtosis(prop_values, n_heavy, mean, std);
        values[out_idx++].d = compute_median(sorted, n_heavy);
        values[out_idx++].d = compute_iqr(sorted, n_heavy);
        values[out_idx++].d = compute_entropy(prop_values, n_heavy, min_val, max_val);
        values[out_idx++].d = compute_gini(sorted, n_heavy);
    }

    if (heap_alloc) {
        free(prop_values);
        free(sorted);
    }

    return NUM_MOMENTS_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* moments_cached_mol = NULL;
static _Thread_local uint64_t moments_cached_gen = 0;
static _Thread_local descriptor_value_t moments_cached_values[NUM_MOMENTS_DESCRIPTORS];

static inline void ensure_moments_computed(const molecule_t* mol) {
    uint64_t current_gen = descriptor_cache_generation();
    if (moments_cached_mol != mol || moments_cached_gen != current_gen) {
        descriptors_compute_moments_all(mol, moments_cached_values);
        moments_cached_mol = mol;
        moments_cached_gen = current_gen;
    }
}

#define DEFINE_MOMENTS_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_moments_computed(mol); \
    value->d = moments_cached_values[idx].d; \
    return CCHEM_OK; \
}

/* Mass moments: idx 0-5 */
DEFINE_MOMENTS_FUNC(mass_skew, 0)
DEFINE_MOMENTS_FUNC(mass_kurt, 1)
DEFINE_MOMENTS_FUNC(mass_median, 2)
DEFINE_MOMENTS_FUNC(mass_iqr, 3)
DEFINE_MOMENTS_FUNC(mass_entropy, 4)
DEFINE_MOMENTS_FUNC(mass_gini, 5)

/* EN moments: idx 6-11 */
DEFINE_MOMENTS_FUNC(en_skew, 6)
DEFINE_MOMENTS_FUNC(en_kurt, 7)
DEFINE_MOMENTS_FUNC(en_median, 8)
DEFINE_MOMENTS_FUNC(en_iqr, 9)
DEFINE_MOMENTS_FUNC(en_entropy, 10)
DEFINE_MOMENTS_FUNC(en_gini, 11)

/* Polarizability moments: idx 12-17 */
DEFINE_MOMENTS_FUNC(pol_skew, 12)
DEFINE_MOMENTS_FUNC(pol_kurt, 13)
DEFINE_MOMENTS_FUNC(pol_median, 14)
DEFINE_MOMENTS_FUNC(pol_iqr, 15)
DEFINE_MOMENTS_FUNC(pol_entropy, 16)
DEFINE_MOMENTS_FUNC(pol_gini, 17)

/* IP moments: idx 18-23 */
DEFINE_MOMENTS_FUNC(ip_skew, 18)
DEFINE_MOMENTS_FUNC(ip_kurt, 19)
DEFINE_MOMENTS_FUNC(ip_median, 20)
DEFINE_MOMENTS_FUNC(ip_iqr, 21)
DEFINE_MOMENTS_FUNC(ip_entropy, 22)
DEFINE_MOMENTS_FUNC(ip_gini, 23)

/* EA moments: idx 24-29 */
DEFINE_MOMENTS_FUNC(ea_skew, 24)
DEFINE_MOMENTS_FUNC(ea_kurt, 25)
DEFINE_MOMENTS_FUNC(ea_median, 26)
DEFINE_MOMENTS_FUNC(ea_iqr, 27)
DEFINE_MOMENTS_FUNC(ea_entropy, 28)
DEFINE_MOMENTS_FUNC(ea_gini, 29)

/* VdW moments: idx 30-35 */
DEFINE_MOMENTS_FUNC(vdw_skew, 30)
DEFINE_MOMENTS_FUNC(vdw_kurt, 31)
DEFINE_MOMENTS_FUNC(vdw_median, 32)
DEFINE_MOMENTS_FUNC(vdw_iqr, 33)
DEFINE_MOMENTS_FUNC(vdw_entropy, 34)
DEFINE_MOMENTS_FUNC(vdw_gini, 35)

/* Valence moments: idx 36-41 */
DEFINE_MOMENTS_FUNC(val_skew, 36)
DEFINE_MOMENTS_FUNC(val_kurt, 37)
DEFINE_MOMENTS_FUNC(val_median, 38)
DEFINE_MOMENTS_FUNC(val_iqr, 39)
DEFINE_MOMENTS_FUNC(val_entropy, 40)
DEFINE_MOMENTS_FUNC(val_gini, 41)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_MOMENT(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_moments(void) {
    /* Mass moments */
    REGISTER_MOMENT("Mass_Skew", "Skewness of atomic mass distribution", desc_mass_skew);
    REGISTER_MOMENT("Mass_Kurt", "Kurtosis of atomic mass distribution", desc_mass_kurt);
    REGISTER_MOMENT("Mass_Median", "Median atomic mass", desc_mass_median);
    REGISTER_MOMENT("Mass_IQR", "IQR of atomic mass", desc_mass_iqr);
    REGISTER_MOMENT("Mass_Entropy", "Entropy of atomic mass distribution", desc_mass_entropy);
    REGISTER_MOMENT("Mass_Gini", "Gini coefficient of atomic mass", desc_mass_gini);

    /* EN moments */
    REGISTER_MOMENT("EN_Skew", "Skewness of electronegativity distribution", desc_en_skew);
    REGISTER_MOMENT("EN_Kurt", "Kurtosis of electronegativity distribution", desc_en_kurt);
    REGISTER_MOMENT("EN_Median", "Median electronegativity", desc_en_median);
    REGISTER_MOMENT("EN_IQR", "IQR of electronegativity", desc_en_iqr);
    REGISTER_MOMENT("EN_Entropy", "Entropy of electronegativity distribution", desc_en_entropy);
    REGISTER_MOMENT("EN_Gini", "Gini coefficient of electronegativity", desc_en_gini);

    /* Polarizability moments */
    REGISTER_MOMENT("Pol_Skew", "Skewness of polarizability distribution", desc_pol_skew);
    REGISTER_MOMENT("Pol_Kurt", "Kurtosis of polarizability distribution", desc_pol_kurt);
    REGISTER_MOMENT("Pol_Median", "Median polarizability", desc_pol_median);
    REGISTER_MOMENT("Pol_IQR", "IQR of polarizability", desc_pol_iqr);
    REGISTER_MOMENT("Pol_Entropy", "Entropy of polarizability distribution", desc_pol_entropy);
    REGISTER_MOMENT("Pol_Gini", "Gini coefficient of polarizability", desc_pol_gini);

    /* IP moments */
    REGISTER_MOMENT("IP_Skew", "Skewness of ionization potential distribution", desc_ip_skew);
    REGISTER_MOMENT("IP_Kurt", "Kurtosis of ionization potential distribution", desc_ip_kurt);
    REGISTER_MOMENT("IP_Median", "Median ionization potential", desc_ip_median);
    REGISTER_MOMENT("IP_IQR", "IQR of ionization potential", desc_ip_iqr);
    REGISTER_MOMENT("IP_Entropy", "Entropy of ionization potential distribution", desc_ip_entropy);
    REGISTER_MOMENT("IP_Gini", "Gini coefficient of ionization potential", desc_ip_gini);

    /* EA moments */
    REGISTER_MOMENT("EA_Skew", "Skewness of electron affinity distribution", desc_ea_skew);
    REGISTER_MOMENT("EA_Kurt", "Kurtosis of electron affinity distribution", desc_ea_kurt);
    REGISTER_MOMENT("EA_Median", "Median electron affinity", desc_ea_median);
    REGISTER_MOMENT("EA_IQR", "IQR of electron affinity", desc_ea_iqr);
    REGISTER_MOMENT("EA_Entropy", "Entropy of electron affinity distribution", desc_ea_entropy);
    REGISTER_MOMENT("EA_Gini", "Gini coefficient of electron affinity", desc_ea_gini);

    /* VdW moments */
    REGISTER_MOMENT("VdW_Skew", "Skewness of VdW volume distribution", desc_vdw_skew);
    REGISTER_MOMENT("VdW_Kurt", "Kurtosis of VdW volume distribution", desc_vdw_kurt);
    REGISTER_MOMENT("VdW_Median", "Median VdW volume", desc_vdw_median);
    REGISTER_MOMENT("VdW_IQR", "IQR of VdW volume", desc_vdw_iqr);
    REGISTER_MOMENT("VdW_Entropy", "Entropy of VdW volume distribution", desc_vdw_entropy);
    REGISTER_MOMENT("VdW_Gini", "Gini coefficient of VdW volume", desc_vdw_gini);

    /* Valence moments */
    REGISTER_MOMENT("Val_Skew", "Skewness of valence electron distribution", desc_val_skew);
    REGISTER_MOMENT("Val_Kurt", "Kurtosis of valence electron distribution", desc_val_kurt);
    REGISTER_MOMENT("Val_Median", "Median valence electrons", desc_val_median);
    REGISTER_MOMENT("Val_IQR", "IQR of valence electrons", desc_val_iqr);
    REGISTER_MOMENT("Val_Entropy", "Entropy of valence electron distribution", desc_val_entropy);
    REGISTER_MOMENT("Val_Gini", "Gini coefficient of valence electrons", desc_val_gini);
}
