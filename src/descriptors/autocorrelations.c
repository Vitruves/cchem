/**
 * @file autocorrelations.c
 * @brief Broto-Moreau Autocorrelation Descriptors
 *
 * Implements 2D autocorrelation descriptors based on topological distance.
 * ATS_d(p) = Σ p_i * p_j for all pairs (i,j) where topological distance d(i,j) = d
 *
 * Properties computed:
 * - ATSm: Atomic mass weighted
 * - ATSv: Van der Waals volume weighted
 * - ATSe: Electronegativity weighted (Sanderson)
 * - ATSp: Polarizability weighted
 * - ATSi: Ionization potential weighted
 * - ATSc: Charge weighted (Gasteiger partial charges)
 *
 * For each property, distances 0-8 are computed (lag 0 = sum of squared properties).
 * Total: 6 properties × 9 lags = 54 descriptors
 *
 * Optimized for speed using BFS-based distance matrix computation.
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Physical Constants for Atomic Properties
 * ============================================================================ */

/* Maximum topological distance (lag) to compute */
#define MAX_LAG 8

/* Number of properties */
#define NUM_PROPERTIES 6

/* Number of lags (0-8) */
#define NUM_LAGS 9

/* Autocorrelation types:
 * ATS:  Broto-Moreau - 6 props × 9 lags = 54
 * AATS: Average ATS  - 6 props × 8 lags = 48 (no lag 0)
 * ATSC: Centered     - 6 props × 8 lags = 48 (no lag 0)
 * MATS: Moran        - 6 props × 8 lags = 48 (no lag 0)
 * GATS: Geary        - 6 props × 8 lags = 48 (no lag 0)
 * Total: 54 + 48*4 = 246 descriptors
 */
#define NUM_ATS_DESCRIPTORS  54
#define NUM_AATS_DESCRIPTORS 48
#define NUM_ATSC_DESCRIPTORS 48
#define NUM_MATS_DESCRIPTORS 48
#define NUM_GATS_DESCRIPTORS 48
#define NUM_AUTOCORR_DESCRIPTORS 246

/* Maximum atoms for stack-allocated distance matrix */
#define MAX_ATOMS_STACK 256

/* Atomic masses (Da) */
static const double ATOMIC_MASS[] = {
    [ELEM_H]  = 1.008,
    [ELEM_C]  = 12.011,
    [ELEM_N]  = 14.007,
    [ELEM_O]  = 15.999,
    [ELEM_F]  = 18.998,
    [ELEM_P]  = 30.974,
    [ELEM_S]  = 32.065,
    [ELEM_Cl] = 35.453,
    [ELEM_Br] = 79.904,
    [ELEM_I]  = 126.904,
    [ELEM_Si] = 28.086,
    [ELEM_B]  = 10.811,
    [ELEM_Se] = 78.96,
    [ELEM_As] = 74.922,
};

/* Van der Waals volumes (Å³) - Bondi radii converted to volume */
static const double VDW_VOLUME[] = {
    [ELEM_H]  = 7.24,
    [ELEM_C]  = 20.58,
    [ELEM_N]  = 15.60,
    [ELEM_O]  = 14.71,
    [ELEM_F]  = 13.31,
    [ELEM_P]  = 24.43,
    [ELEM_S]  = 24.43,
    [ELEM_Cl] = 22.45,
    [ELEM_Br] = 26.52,
    [ELEM_I]  = 32.52,
    [ELEM_Si] = 38.79,
    [ELEM_B]  = 29.64,
    [ELEM_Se] = 28.73,
    [ELEM_As] = 26.52,
};

/* Sanderson electronegativity */
static const double ELECTRONEGATIVITY[] = {
    [ELEM_H]  = 2.59,
    [ELEM_C]  = 2.75,
    [ELEM_N]  = 3.19,
    [ELEM_O]  = 3.65,
    [ELEM_F]  = 4.00,
    [ELEM_P]  = 2.52,
    [ELEM_S]  = 2.96,
    [ELEM_Cl] = 3.48,
    [ELEM_Br] = 3.22,
    [ELEM_I]  = 2.78,
    [ELEM_Si] = 2.14,
    [ELEM_B]  = 2.28,
    [ELEM_Se] = 2.76,
    [ELEM_As] = 2.26,
};

/* Atomic polarizability (Å³) */
static const double POLARIZABILITY[] = {
    [ELEM_H]  = 0.387,
    [ELEM_C]  = 1.76,
    [ELEM_N]  = 1.10,
    [ELEM_O]  = 0.802,
    [ELEM_F]  = 0.557,
    [ELEM_P]  = 3.63,
    [ELEM_S]  = 2.90,
    [ELEM_Cl] = 2.18,
    [ELEM_Br] = 3.05,
    [ELEM_I]  = 4.7,
    [ELEM_Si] = 5.38,
    [ELEM_B]  = 3.03,
    [ELEM_Se] = 3.77,
    [ELEM_As] = 4.31,
};

/* First ionization potential (eV) */
static const double IONIZATION_POTENTIAL[] = {
    [ELEM_H]  = 13.598,
    [ELEM_C]  = 11.260,
    [ELEM_N]  = 14.534,
    [ELEM_O]  = 13.618,
    [ELEM_F]  = 17.422,
    [ELEM_P]  = 10.486,
    [ELEM_S]  = 10.360,
    [ELEM_Cl] = 12.967,
    [ELEM_Br] = 11.814,
    [ELEM_I]  = 10.451,
    [ELEM_Si] = 8.151,
    [ELEM_B]  = 8.298,
    [ELEM_Se] = 9.752,
    [ELEM_As] = 9.815,
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_atomic_mass(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ATOMIC_MASS)/sizeof(ATOMIC_MASS[0])))
        return 12.011;
    double m = ATOMIC_MASS[elem];
    return (m > 0.0) ? m : 12.011;
}

static inline double get_vdw_volume(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(VDW_VOLUME)/sizeof(VDW_VOLUME[0])))
        return 20.58;
    double v = VDW_VOLUME[elem];
    return (v > 0.0) ? v : 20.58;
}

static inline double get_electronegativity(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ELECTRONEGATIVITY)/sizeof(ELECTRONEGATIVITY[0])))
        return 2.75;
    double e = ELECTRONEGATIVITY[elem];
    return (e > 0.0) ? e : 2.75;
}

static inline double get_polarizability(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(POLARIZABILITY)/sizeof(POLARIZABILITY[0])))
        return 1.76;
    double p = POLARIZABILITY[elem];
    return (p > 0.0) ? p : 1.76;
}

static inline double get_ionization_potential(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(IONIZATION_POTENTIAL)/sizeof(IONIZATION_POTENTIAL[0])))
        return 11.26;
    double ip = IONIZATION_POTENTIAL[elem];
    return (ip > 0.0) ? ip : 11.26;
}

/* ============================================================================
 * Distance Matrix Computation (BFS-based, O(n²) for sparse molecular graphs)
 * ============================================================================ */

/**
 * Compute topological distance matrix using BFS from each atom.
 * Only considers heavy atoms (non-H).
 *
 * @param mol Molecule
 * @param dist Output distance matrix (n_heavy × n_heavy), -1 for disconnected
 * @param heavy_idx Output mapping from heavy atom index to original atom index
 * @param n_heavy Output number of heavy atoms
 * @param max_atoms Maximum atoms to process
 */
static void compute_distance_matrix(const molecule_t* mol,
                                    int* dist,
                                    int* heavy_idx,
                                    int* n_heavy,
                                    int max_atoms) {
    /* First pass: identify heavy atoms */
    int nh = 0;
    for (int i = 0; i < mol->num_atoms && nh < max_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_idx[nh++] = i;
        }
    }
    *n_heavy = nh;

    if (nh == 0) return;

    /* Create reverse mapping: original index -> heavy index */
    int* reverse_map = (int*)alloca(mol->num_atoms * sizeof(int));
    memset(reverse_map, -1, mol->num_atoms * sizeof(int));
    for (int i = 0; i < nh; i++) {
        reverse_map[heavy_idx[i]] = i;
    }

    /* Initialize distance matrix to -1 (disconnected) */
    for (int i = 0; i < nh * nh; i++) {
        dist[i] = -1;
    }

    /* BFS queue (stack allocated for small molecules) */
    int* queue = (int*)alloca(nh * sizeof(int));

    /* BFS from each heavy atom */
    for (int src = 0; src < nh; src++) {
        int src_orig = heavy_idx[src];
        dist[src * nh + src] = 0;  /* Distance to self is 0 */

        int head = 0, tail = 0;
        queue[tail++] = src_orig;

        /* Visited array for this BFS */
        bool* visited = (bool*)alloca(mol->num_atoms * sizeof(bool));
        memset(visited, 0, mol->num_atoms * sizeof(bool));
        visited[src_orig] = true;

        while (head < tail) {
            int curr = queue[head++];
            int curr_heavy = reverse_map[curr];
            int curr_dist = (curr_heavy >= 0) ? dist[src * nh + curr_heavy] : -1;

            const atom_t* atom = &mol->atoms[curr];
            for (int j = 0; j < atom->num_neighbors; j++) {
                int neighbor = atom->neighbors[j];
                if (visited[neighbor]) continue;
                visited[neighbor] = true;

                int neighbor_heavy = reverse_map[neighbor];

                if (mol->atoms[neighbor].element == ELEM_H) {
                    /* Skip H but continue BFS through it (for bridging) */
                    /* Actually, in molecular graphs H is terminal, so skip */
                    continue;
                }

                /* Set distance for heavy neighbor */
                if (neighbor_heavy >= 0 && curr_heavy >= 0) {
                    dist[src * nh + neighbor_heavy] = curr_dist + 1;
                }
                queue[tail++] = neighbor;
            }
        }
    }
}

/* ============================================================================
 * Gasteiger Partial Charges (simplified PEOE)
 * ============================================================================ */

/* Gasteiger parameters: a, b, c for EN = a + b*q + c*q² */
typedef struct {
    double a, b, c;
} gasteiger_params_t;

static const gasteiger_params_t GASTEIGER_PARAMS[] = {
    [ELEM_H]  = {7.17, 6.24, -0.56},
    [ELEM_C]  = {7.98, 9.18, 1.88},
    [ELEM_N]  = {11.54, 10.82, 1.36},
    [ELEM_O]  = {14.18, 12.92, 1.39},
    [ELEM_F]  = {14.66, 13.85, 2.31},
    [ELEM_P]  = {8.90, 8.24, 1.62},
    [ELEM_S]  = {10.14, 9.13, 1.38},
    [ELEM_Cl] = {11.00, 9.69, 1.35},
    [ELEM_Br] = {10.08, 8.47, 1.16},
    [ELEM_I]  = {9.90, 7.96, 0.96},
};

static void compute_gasteiger_charges(const molecule_t* mol, double* charges) {
    const int n = mol->num_atoms;
    const int max_iter = 6;

    /* Initialize charges to formal charge */
    for (int i = 0; i < n; i++) {
        charges[i] = mol->atoms[i].charge;
    }

    /* Iterative charge equilibration */
    double* delta = (double*)alloca(n * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        memset(delta, 0, n * sizeof(double));

        /* Damping factor decreases with iteration */
        double damping = pow(0.5, iter + 1);

        /* For each bond, transfer charge based on EN difference */
        for (int b = 0; b < mol->num_bonds; b++) {
            const bond_t* bond = &mol->bonds[b];
            int i = bond->atom1;
            int j = bond->atom2;

            element_t ei = mol->atoms[i].element;
            element_t ej = mol->atoms[j].element;

            /* Get Gasteiger params */
            gasteiger_params_t pi = {7.98, 9.18, 1.88};  /* Default to C */
            gasteiger_params_t pj = {7.98, 9.18, 1.88};

            if (ei > 0 && ei < (int)(sizeof(GASTEIGER_PARAMS)/sizeof(GASTEIGER_PARAMS[0]))) {
                if (GASTEIGER_PARAMS[ei].a > 0) pi = GASTEIGER_PARAMS[ei];
            }
            if (ej > 0 && ej < (int)(sizeof(GASTEIGER_PARAMS)/sizeof(GASTEIGER_PARAMS[0]))) {
                if (GASTEIGER_PARAMS[ej].a > 0) pj = GASTEIGER_PARAMS[ej];
            }

            /* Compute EN at current charge (clamped for stability) */
            double qi = fmax(-1.0, fmin(1.0, charges[i]));
            double qj = fmax(-1.0, fmin(1.0, charges[j]));
            double eni = pi.a + pi.b * qi + pi.c * qi * qi;
            double enj = pj.a + pj.b * qj + pj.c * qj * qj;

            /* Charge transfer proportional to EN difference */
            double dq = (enj - eni) / (pi.b + pj.b) * damping;

            /* Clamp delta to prevent instability */
            dq = fmax(-0.1, fmin(0.1, dq));

            delta[i] += dq;
            delta[j] -= dq;
        }

        /* Apply deltas with clamping */
        for (int i = 0; i < n; i++) {
            charges[i] += delta[i];
            /* Clamp to reasonable range */
            charges[i] = fmax(-2.0, fmin(2.0, charges[i]));
        }
    }
}

/* ============================================================================
 * Autocorrelation Computation
 * ============================================================================ */

/**
 * Compute Broto-Moreau autocorrelation for a given property at distance d.
 * ATS_d(p) = Σ p_i * p_j for all pairs (i,j) where d(i,j) = d
 */
static double compute_ats(const double* props, const int* dist, int n, int d) {
    double sum = 0.0;

    if (d == 0) {
        /* Lag 0: sum of squared properties */
        for (int i = 0; i < n; i++) {
            sum += props[i] * props[i];
        }
    } else {
        /* Lag > 0: sum over pairs at distance d */
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (dist[i * n + j] == d) {
                    sum += props[i] * props[j];
                }
            }
        }
        sum *= 2.0;  /* Count both (i,j) and (j,i) */
    }

    return sum;
}

/**
 * Count pairs at distance d
 */
static int count_pairs_at_distance(const int* dist, int n, int d) {
    int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i * n + j] == d) count++;
        }
    }
    return count;
}

/**
 * Compute mean of properties
 */
static double compute_mean(const double* props, int n) {
    if (n == 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += props[i];
    return sum / n;
}

/**
 * Compute variance of properties
 */
static double compute_variance(const double* props, int n, double mean) {
    if (n <= 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double d = props[i] - mean;
        sum += d * d;
    }
    return sum / n;
}

/**
 * AATS: Average ATS - normalized by number of pairs
 * AATS_d = ATS_d / n_pairs_d
 */
static double compute_aats(const double* props, const int* dist, int n, int d) {
    if (d == 0) return 0.0;  /* AATS not defined for lag 0 */

    int n_pairs = count_pairs_at_distance(dist, n, d);
    if (n_pairs == 0) return 0.0;

    double ats = compute_ats(props, dist, n, d) / 2.0;  /* Undo the *2 */
    return ats / n_pairs;
}

/**
 * ATSC: Centered ATS - uses centered properties
 * ATSC_d = Σ (p_i - mean) * (p_j - mean) for d(i,j) = d
 */
static double compute_atsc(const double* props, const int* dist, int n, int d, double mean) {
    if (d == 0) return 0.0;  /* ATSC not defined for lag 0 */

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i * n + j] == d) {
                sum += (props[i] - mean) * (props[j] - mean);
            }
        }
    }
    return sum * 2.0;
}

/**
 * MATS: Moran autocorrelation coefficient
 * MATS_d = [n * Σ(p_i - mean)(p_j - mean)] / [Σ(p_i - mean)^2 * n_pairs_d]
 */
static double compute_mats(const double* props, const int* dist, int n, int d, double mean, double var) {
    if (d == 0 || var < 1e-10) return 0.0;

    int n_pairs = count_pairs_at_distance(dist, n, d);
    if (n_pairs == 0) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i * n + j] == d) {
                sum += (props[i] - mean) * (props[j] - mean);
            }
        }
    }

    /* Moran's I formula */
    return (n * sum) / (var * n * n_pairs);
}

/**
 * GATS: Geary autocorrelation coefficient
 * GATS_d = [(n-1) * Σ(p_i - p_j)^2] / [2 * n_pairs_d * Σ(p_i - mean)^2]
 */
static double compute_gats(const double* props, const int* dist, int n, int d, double var) {
    if (d == 0 || n <= 1 || var < 1e-10) return 0.0;

    int n_pairs = count_pairs_at_distance(dist, n, d);
    if (n_pairs == 0) return 0.0;

    double sum_sq_diff = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i * n + j] == d) {
                double diff = props[i] - props[j];
                sum_sq_diff += diff * diff;
            }
        }
    }

    /* Geary's C formula */
    return ((n - 1) * sum_sq_diff) / (2.0 * n_pairs * var * n);
}

/* ============================================================================
 * Batch Computation (thread-safe)
 * ============================================================================ */

int descriptors_compute_autocorr_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    /* Initialize all to 0 */
    for (int i = 0; i < NUM_AUTOCORR_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    /* Count heavy atoms */
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    }

    if (n_heavy == 0) return NUM_AUTOCORR_DESCRIPTORS;

    /* Use stack allocation for small molecules, heap for large */
    int* dist_matrix;
    int* heavy_idx;
    double* props_mass;
    double* props_vol;
    double* props_en;
    double* props_pol;
    double* props_ip;
    double* props_charge;
    double* gasteiger;

    bool heap_alloc = (n_heavy > MAX_ATOMS_STACK);

    if (heap_alloc) {
        dist_matrix = (int*)malloc(n_heavy * n_heavy * sizeof(int));
        heavy_idx = (int*)malloc(n_heavy * sizeof(int));
        props_mass = (double*)malloc(n_heavy * sizeof(double));
        props_vol = (double*)malloc(n_heavy * sizeof(double));
        props_en = (double*)malloc(n_heavy * sizeof(double));
        props_pol = (double*)malloc(n_heavy * sizeof(double));
        props_ip = (double*)malloc(n_heavy * sizeof(double));
        props_charge = (double*)malloc(n_heavy * sizeof(double));
        gasteiger = (double*)malloc(mol->num_atoms * sizeof(double));

        if (!dist_matrix || !heavy_idx || !props_mass || !props_vol ||
            !props_en || !props_pol || !props_ip || !props_charge || !gasteiger) {
            free(dist_matrix); free(heavy_idx); free(props_mass); free(props_vol);
            free(props_en); free(props_pol); free(props_ip); free(props_charge);
            free(gasteiger);
            return -1;
        }
    } else {
        dist_matrix = (int*)alloca(n_heavy * n_heavy * sizeof(int));
        heavy_idx = (int*)alloca(n_heavy * sizeof(int));
        props_mass = (double*)alloca(n_heavy * sizeof(double));
        props_vol = (double*)alloca(n_heavy * sizeof(double));
        props_en = (double*)alloca(n_heavy * sizeof(double));
        props_pol = (double*)alloca(n_heavy * sizeof(double));
        props_ip = (double*)alloca(n_heavy * sizeof(double));
        props_charge = (double*)alloca(n_heavy * sizeof(double));
        gasteiger = (double*)alloca(mol->num_atoms * sizeof(double));
    }

    /* Compute distance matrix */
    int actual_heavy;
    compute_distance_matrix(mol, dist_matrix, heavy_idx, &actual_heavy, n_heavy);

    if (actual_heavy == 0) {
        if (heap_alloc) {
            free(dist_matrix); free(heavy_idx); free(props_mass); free(props_vol);
            free(props_en); free(props_pol); free(props_ip); free(props_charge);
            free(gasteiger);
        }
        return NUM_AUTOCORR_DESCRIPTORS;
    }

    /* Compute Gasteiger charges */
    compute_gasteiger_charges(mol, gasteiger);

    /* Extract properties for heavy atoms */
    for (int i = 0; i < actual_heavy; i++) {
        int orig_idx = heavy_idx[i];
        element_t elem = mol->atoms[orig_idx].element;

        props_mass[i] = get_atomic_mass(elem);
        props_vol[i] = get_vdw_volume(elem);
        props_en[i] = get_electronegativity(elem);
        props_pol[i] = get_polarizability(elem);
        props_ip[i] = get_ionization_potential(elem);
        props_charge[i] = gasteiger[orig_idx];
    }

    /* Compute means and variances for centered/normalized autocorrelations */
    double mean_mass = compute_mean(props_mass, actual_heavy);
    double mean_vol = compute_mean(props_vol, actual_heavy);
    double mean_en = compute_mean(props_en, actual_heavy);
    double mean_pol = compute_mean(props_pol, actual_heavy);
    double mean_ip = compute_mean(props_ip, actual_heavy);
    double mean_charge = compute_mean(props_charge, actual_heavy);

    double var_mass = compute_variance(props_mass, actual_heavy, mean_mass);
    double var_vol = compute_variance(props_vol, actual_heavy, mean_vol);
    double var_en = compute_variance(props_en, actual_heavy, mean_en);
    double var_pol = compute_variance(props_pol, actual_heavy, mean_pol);
    double var_ip = compute_variance(props_ip, actual_heavy, mean_ip);
    double var_charge = compute_variance(props_charge, actual_heavy, mean_charge);

    /* Store property arrays and stats for iteration */
    const double* props[6] = {props_mass, props_vol, props_en, props_pol, props_ip, props_charge};
    double means[6] = {mean_mass, mean_vol, mean_en, mean_pol, mean_ip, mean_charge};
    double vars[6] = {var_mass, var_vol, var_en, var_pol, var_ip, var_charge};

    int idx = 0;

    /* ATS: Broto-Moreau autocorrelations (6 properties × 9 lags) = 54 */
    for (int p = 0; p < NUM_PROPERTIES; p++) {
        for (int d = 0; d <= MAX_LAG; d++) {
            values[idx++].d = compute_ats(props[p], dist_matrix, actual_heavy, d);
        }
    }

    /* AATS: Average ATS (6 properties × 8 lags, no lag 0) = 48 */
    for (int p = 0; p < NUM_PROPERTIES; p++) {
        for (int d = 1; d <= MAX_LAG; d++) {
            values[idx++].d = compute_aats(props[p], dist_matrix, actual_heavy, d);
        }
    }

    /* ATSC: Centered autocorrelations (6 properties × 8 lags, no lag 0) = 48 */
    for (int p = 0; p < NUM_PROPERTIES; p++) {
        for (int d = 1; d <= MAX_LAG; d++) {
            values[idx++].d = compute_atsc(props[p], dist_matrix, actual_heavy, d, means[p]);
        }
    }

    /* MATS: Moran autocorrelation (6 properties × 8 lags, no lag 0) = 48 */
    for (int p = 0; p < NUM_PROPERTIES; p++) {
        for (int d = 1; d <= MAX_LAG; d++) {
            values[idx++].d = compute_mats(props[p], dist_matrix, actual_heavy, d, means[p], vars[p]);
        }
    }

    /* GATS: Geary autocorrelation (6 properties × 8 lags, no lag 0) = 48 */
    for (int p = 0; p < NUM_PROPERTIES; p++) {
        for (int d = 1; d <= MAX_LAG; d++) {
            values[idx++].d = compute_gats(props[p], dist_matrix, actual_heavy, d, vars[p]);
        }
    }

    if (heap_alloc) {
        free(dist_matrix); free(heavy_idx); free(props_mass); free(props_vol);
        free(props_en); free(props_pol); free(props_ip); free(props_charge);
        free(gasteiger);
    }

    return NUM_AUTOCORR_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Thread-Local Caching
 * ============================================================================ */

/* Offset for each autocorrelation type in the batch output */
#define ATS_OFFSET  0
#define AATS_OFFSET NUM_ATS_DESCRIPTORS
#define ATSC_OFFSET (AATS_OFFSET + NUM_AATS_DESCRIPTORS)
#define MATS_OFFSET (ATSC_OFFSET + NUM_ATSC_DESCRIPTORS)
#define GATS_OFFSET (MATS_OFFSET + NUM_MATS_DESCRIPTORS)

/* Thread-local cache to avoid recomputing for same molecule */
static _Thread_local const molecule_t* cached_mol = NULL;
static _Thread_local descriptor_value_t cached_values[NUM_AUTOCORR_DESCRIPTORS];

static inline void ensure_autocorr_computed(const molecule_t* mol) {
    if (cached_mol != mol) {
        descriptors_compute_autocorr_all(mol, cached_values);
        cached_mol = mol;
    }
}

/* Macro for generating individual ATS descriptor functions (has lag 0) */
#define DEFINE_ATS_FUNC(name, prop_idx, lag) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_autocorr_computed(mol); \
    value->d = cached_values[ATS_OFFSET + (prop_idx) * 9 + (lag)].d; \
    return CCHEM_OK; \
}

/* Macro for generating AATS descriptor functions (lags 1-8 only) */
#define DEFINE_AATS_FUNC(name, prop_idx, lag) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_autocorr_computed(mol); \
    value->d = cached_values[AATS_OFFSET + (prop_idx) * 8 + ((lag) - 1)].d; \
    return CCHEM_OK; \
}

/* Macro for generating ATSC descriptor functions (lags 1-8 only) */
#define DEFINE_ATSC_FUNC(name, prop_idx, lag) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_autocorr_computed(mol); \
    value->d = cached_values[ATSC_OFFSET + (prop_idx) * 8 + ((lag) - 1)].d; \
    return CCHEM_OK; \
}

/* Macro for generating MATS descriptor functions (lags 1-8 only) */
#define DEFINE_MATS_FUNC(name, prop_idx, lag) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_autocorr_computed(mol); \
    value->d = cached_values[MATS_OFFSET + (prop_idx) * 8 + ((lag) - 1)].d; \
    return CCHEM_OK; \
}

/* Macro for generating GATS descriptor functions (lags 1-8 only) */
#define DEFINE_GATS_FUNC(name, prop_idx, lag) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_autocorr_computed(mol); \
    value->d = cached_values[GATS_OFFSET + (prop_idx) * 8 + ((lag) - 1)].d; \
    return CCHEM_OK; \
}

/* ATSm: Mass-weighted (property index 0) */
DEFINE_ATS_FUNC(ats_m0, 0, 0)
DEFINE_ATS_FUNC(ats_m1, 0, 1)
DEFINE_ATS_FUNC(ats_m2, 0, 2)
DEFINE_ATS_FUNC(ats_m3, 0, 3)
DEFINE_ATS_FUNC(ats_m4, 0, 4)
DEFINE_ATS_FUNC(ats_m5, 0, 5)
DEFINE_ATS_FUNC(ats_m6, 0, 6)
DEFINE_ATS_FUNC(ats_m7, 0, 7)
DEFINE_ATS_FUNC(ats_m8, 0, 8)

/* ATSv: Volume-weighted (property index 1) */
DEFINE_ATS_FUNC(ats_v0, 1, 0)
DEFINE_ATS_FUNC(ats_v1, 1, 1)
DEFINE_ATS_FUNC(ats_v2, 1, 2)
DEFINE_ATS_FUNC(ats_v3, 1, 3)
DEFINE_ATS_FUNC(ats_v4, 1, 4)
DEFINE_ATS_FUNC(ats_v5, 1, 5)
DEFINE_ATS_FUNC(ats_v6, 1, 6)
DEFINE_ATS_FUNC(ats_v7, 1, 7)
DEFINE_ATS_FUNC(ats_v8, 1, 8)

/* ATSe: Electronegativity-weighted (property index 2) */
DEFINE_ATS_FUNC(ats_e0, 2, 0)
DEFINE_ATS_FUNC(ats_e1, 2, 1)
DEFINE_ATS_FUNC(ats_e2, 2, 2)
DEFINE_ATS_FUNC(ats_e3, 2, 3)
DEFINE_ATS_FUNC(ats_e4, 2, 4)
DEFINE_ATS_FUNC(ats_e5, 2, 5)
DEFINE_ATS_FUNC(ats_e6, 2, 6)
DEFINE_ATS_FUNC(ats_e7, 2, 7)
DEFINE_ATS_FUNC(ats_e8, 2, 8)

/* ATSp: Polarizability-weighted (property index 3) */
DEFINE_ATS_FUNC(ats_p0, 3, 0)
DEFINE_ATS_FUNC(ats_p1, 3, 1)
DEFINE_ATS_FUNC(ats_p2, 3, 2)
DEFINE_ATS_FUNC(ats_p3, 3, 3)
DEFINE_ATS_FUNC(ats_p4, 3, 4)
DEFINE_ATS_FUNC(ats_p5, 3, 5)
DEFINE_ATS_FUNC(ats_p6, 3, 6)
DEFINE_ATS_FUNC(ats_p7, 3, 7)
DEFINE_ATS_FUNC(ats_p8, 3, 8)

/* ATSi: Ionization potential-weighted (property index 4) */
DEFINE_ATS_FUNC(ats_i0, 4, 0)
DEFINE_ATS_FUNC(ats_i1, 4, 1)
DEFINE_ATS_FUNC(ats_i2, 4, 2)
DEFINE_ATS_FUNC(ats_i3, 4, 3)
DEFINE_ATS_FUNC(ats_i4, 4, 4)
DEFINE_ATS_FUNC(ats_i5, 4, 5)
DEFINE_ATS_FUNC(ats_i6, 4, 6)
DEFINE_ATS_FUNC(ats_i7, 4, 7)
DEFINE_ATS_FUNC(ats_i8, 4, 8)

/* ATSc: Charge-weighted (property index 5) */
DEFINE_ATS_FUNC(ats_c0, 5, 0)
DEFINE_ATS_FUNC(ats_c1, 5, 1)
DEFINE_ATS_FUNC(ats_c2, 5, 2)
DEFINE_ATS_FUNC(ats_c3, 5, 3)
DEFINE_ATS_FUNC(ats_c4, 5, 4)
DEFINE_ATS_FUNC(ats_c5, 5, 5)
DEFINE_ATS_FUNC(ats_c6, 5, 6)
DEFINE_ATS_FUNC(ats_c7, 5, 7)
DEFINE_ATS_FUNC(ats_c8, 5, 8)

/* ============================================================================
 * AATS: Average Autocorrelations (lags 1-8 only)
 * ============================================================================ */

/* AATSm: Mass-weighted (property index 0) */
DEFINE_AATS_FUNC(aats_m1, 0, 1)
DEFINE_AATS_FUNC(aats_m2, 0, 2)
DEFINE_AATS_FUNC(aats_m3, 0, 3)
DEFINE_AATS_FUNC(aats_m4, 0, 4)
DEFINE_AATS_FUNC(aats_m5, 0, 5)
DEFINE_AATS_FUNC(aats_m6, 0, 6)
DEFINE_AATS_FUNC(aats_m7, 0, 7)
DEFINE_AATS_FUNC(aats_m8, 0, 8)

/* AATSv: Volume-weighted (property index 1) */
DEFINE_AATS_FUNC(aats_v1, 1, 1)
DEFINE_AATS_FUNC(aats_v2, 1, 2)
DEFINE_AATS_FUNC(aats_v3, 1, 3)
DEFINE_AATS_FUNC(aats_v4, 1, 4)
DEFINE_AATS_FUNC(aats_v5, 1, 5)
DEFINE_AATS_FUNC(aats_v6, 1, 6)
DEFINE_AATS_FUNC(aats_v7, 1, 7)
DEFINE_AATS_FUNC(aats_v8, 1, 8)

/* AATSe: Electronegativity-weighted (property index 2) */
DEFINE_AATS_FUNC(aats_e1, 2, 1)
DEFINE_AATS_FUNC(aats_e2, 2, 2)
DEFINE_AATS_FUNC(aats_e3, 2, 3)
DEFINE_AATS_FUNC(aats_e4, 2, 4)
DEFINE_AATS_FUNC(aats_e5, 2, 5)
DEFINE_AATS_FUNC(aats_e6, 2, 6)
DEFINE_AATS_FUNC(aats_e7, 2, 7)
DEFINE_AATS_FUNC(aats_e8, 2, 8)

/* AATSp: Polarizability-weighted (property index 3) */
DEFINE_AATS_FUNC(aats_p1, 3, 1)
DEFINE_AATS_FUNC(aats_p2, 3, 2)
DEFINE_AATS_FUNC(aats_p3, 3, 3)
DEFINE_AATS_FUNC(aats_p4, 3, 4)
DEFINE_AATS_FUNC(aats_p5, 3, 5)
DEFINE_AATS_FUNC(aats_p6, 3, 6)
DEFINE_AATS_FUNC(aats_p7, 3, 7)
DEFINE_AATS_FUNC(aats_p8, 3, 8)

/* AATSi: Ionization potential-weighted (property index 4) */
DEFINE_AATS_FUNC(aats_i1, 4, 1)
DEFINE_AATS_FUNC(aats_i2, 4, 2)
DEFINE_AATS_FUNC(aats_i3, 4, 3)
DEFINE_AATS_FUNC(aats_i4, 4, 4)
DEFINE_AATS_FUNC(aats_i5, 4, 5)
DEFINE_AATS_FUNC(aats_i6, 4, 6)
DEFINE_AATS_FUNC(aats_i7, 4, 7)
DEFINE_AATS_FUNC(aats_i8, 4, 8)

/* AATSc: Charge-weighted (property index 5) */
DEFINE_AATS_FUNC(aats_c1, 5, 1)
DEFINE_AATS_FUNC(aats_c2, 5, 2)
DEFINE_AATS_FUNC(aats_c3, 5, 3)
DEFINE_AATS_FUNC(aats_c4, 5, 4)
DEFINE_AATS_FUNC(aats_c5, 5, 5)
DEFINE_AATS_FUNC(aats_c6, 5, 6)
DEFINE_AATS_FUNC(aats_c7, 5, 7)
DEFINE_AATS_FUNC(aats_c8, 5, 8)

/* ============================================================================
 * ATSC: Centered Autocorrelations (lags 1-8 only)
 * ============================================================================ */

/* ATSCm: Mass-weighted (property index 0) */
DEFINE_ATSC_FUNC(atsc_m1, 0, 1)
DEFINE_ATSC_FUNC(atsc_m2, 0, 2)
DEFINE_ATSC_FUNC(atsc_m3, 0, 3)
DEFINE_ATSC_FUNC(atsc_m4, 0, 4)
DEFINE_ATSC_FUNC(atsc_m5, 0, 5)
DEFINE_ATSC_FUNC(atsc_m6, 0, 6)
DEFINE_ATSC_FUNC(atsc_m7, 0, 7)
DEFINE_ATSC_FUNC(atsc_m8, 0, 8)

/* ATSCv: Volume-weighted (property index 1) */
DEFINE_ATSC_FUNC(atsc_v1, 1, 1)
DEFINE_ATSC_FUNC(atsc_v2, 1, 2)
DEFINE_ATSC_FUNC(atsc_v3, 1, 3)
DEFINE_ATSC_FUNC(atsc_v4, 1, 4)
DEFINE_ATSC_FUNC(atsc_v5, 1, 5)
DEFINE_ATSC_FUNC(atsc_v6, 1, 6)
DEFINE_ATSC_FUNC(atsc_v7, 1, 7)
DEFINE_ATSC_FUNC(atsc_v8, 1, 8)

/* ATSCe: Electronegativity-weighted (property index 2) */
DEFINE_ATSC_FUNC(atsc_e1, 2, 1)
DEFINE_ATSC_FUNC(atsc_e2, 2, 2)
DEFINE_ATSC_FUNC(atsc_e3, 2, 3)
DEFINE_ATSC_FUNC(atsc_e4, 2, 4)
DEFINE_ATSC_FUNC(atsc_e5, 2, 5)
DEFINE_ATSC_FUNC(atsc_e6, 2, 6)
DEFINE_ATSC_FUNC(atsc_e7, 2, 7)
DEFINE_ATSC_FUNC(atsc_e8, 2, 8)

/* ATSCp: Polarizability-weighted (property index 3) */
DEFINE_ATSC_FUNC(atsc_p1, 3, 1)
DEFINE_ATSC_FUNC(atsc_p2, 3, 2)
DEFINE_ATSC_FUNC(atsc_p3, 3, 3)
DEFINE_ATSC_FUNC(atsc_p4, 3, 4)
DEFINE_ATSC_FUNC(atsc_p5, 3, 5)
DEFINE_ATSC_FUNC(atsc_p6, 3, 6)
DEFINE_ATSC_FUNC(atsc_p7, 3, 7)
DEFINE_ATSC_FUNC(atsc_p8, 3, 8)

/* ATSCi: Ionization potential-weighted (property index 4) */
DEFINE_ATSC_FUNC(atsc_i1, 4, 1)
DEFINE_ATSC_FUNC(atsc_i2, 4, 2)
DEFINE_ATSC_FUNC(atsc_i3, 4, 3)
DEFINE_ATSC_FUNC(atsc_i4, 4, 4)
DEFINE_ATSC_FUNC(atsc_i5, 4, 5)
DEFINE_ATSC_FUNC(atsc_i6, 4, 6)
DEFINE_ATSC_FUNC(atsc_i7, 4, 7)
DEFINE_ATSC_FUNC(atsc_i8, 4, 8)

/* ATSCc: Charge-weighted (property index 5) */
DEFINE_ATSC_FUNC(atsc_c1, 5, 1)
DEFINE_ATSC_FUNC(atsc_c2, 5, 2)
DEFINE_ATSC_FUNC(atsc_c3, 5, 3)
DEFINE_ATSC_FUNC(atsc_c4, 5, 4)
DEFINE_ATSC_FUNC(atsc_c5, 5, 5)
DEFINE_ATSC_FUNC(atsc_c6, 5, 6)
DEFINE_ATSC_FUNC(atsc_c7, 5, 7)
DEFINE_ATSC_FUNC(atsc_c8, 5, 8)

/* ============================================================================
 * MATS: Moran Autocorrelations (lags 1-8 only)
 * ============================================================================ */

/* MATSm: Mass-weighted (property index 0) */
DEFINE_MATS_FUNC(mats_m1, 0, 1)
DEFINE_MATS_FUNC(mats_m2, 0, 2)
DEFINE_MATS_FUNC(mats_m3, 0, 3)
DEFINE_MATS_FUNC(mats_m4, 0, 4)
DEFINE_MATS_FUNC(mats_m5, 0, 5)
DEFINE_MATS_FUNC(mats_m6, 0, 6)
DEFINE_MATS_FUNC(mats_m7, 0, 7)
DEFINE_MATS_FUNC(mats_m8, 0, 8)

/* MATSv: Volume-weighted (property index 1) */
DEFINE_MATS_FUNC(mats_v1, 1, 1)
DEFINE_MATS_FUNC(mats_v2, 1, 2)
DEFINE_MATS_FUNC(mats_v3, 1, 3)
DEFINE_MATS_FUNC(mats_v4, 1, 4)
DEFINE_MATS_FUNC(mats_v5, 1, 5)
DEFINE_MATS_FUNC(mats_v6, 1, 6)
DEFINE_MATS_FUNC(mats_v7, 1, 7)
DEFINE_MATS_FUNC(mats_v8, 1, 8)

/* MATSe: Electronegativity-weighted (property index 2) */
DEFINE_MATS_FUNC(mats_e1, 2, 1)
DEFINE_MATS_FUNC(mats_e2, 2, 2)
DEFINE_MATS_FUNC(mats_e3, 2, 3)
DEFINE_MATS_FUNC(mats_e4, 2, 4)
DEFINE_MATS_FUNC(mats_e5, 2, 5)
DEFINE_MATS_FUNC(mats_e6, 2, 6)
DEFINE_MATS_FUNC(mats_e7, 2, 7)
DEFINE_MATS_FUNC(mats_e8, 2, 8)

/* MATSp: Polarizability-weighted (property index 3) */
DEFINE_MATS_FUNC(mats_p1, 3, 1)
DEFINE_MATS_FUNC(mats_p2, 3, 2)
DEFINE_MATS_FUNC(mats_p3, 3, 3)
DEFINE_MATS_FUNC(mats_p4, 3, 4)
DEFINE_MATS_FUNC(mats_p5, 3, 5)
DEFINE_MATS_FUNC(mats_p6, 3, 6)
DEFINE_MATS_FUNC(mats_p7, 3, 7)
DEFINE_MATS_FUNC(mats_p8, 3, 8)

/* MATSi: Ionization potential-weighted (property index 4) */
DEFINE_MATS_FUNC(mats_i1, 4, 1)
DEFINE_MATS_FUNC(mats_i2, 4, 2)
DEFINE_MATS_FUNC(mats_i3, 4, 3)
DEFINE_MATS_FUNC(mats_i4, 4, 4)
DEFINE_MATS_FUNC(mats_i5, 4, 5)
DEFINE_MATS_FUNC(mats_i6, 4, 6)
DEFINE_MATS_FUNC(mats_i7, 4, 7)
DEFINE_MATS_FUNC(mats_i8, 4, 8)

/* MATSc: Charge-weighted (property index 5) */
DEFINE_MATS_FUNC(mats_c1, 5, 1)
DEFINE_MATS_FUNC(mats_c2, 5, 2)
DEFINE_MATS_FUNC(mats_c3, 5, 3)
DEFINE_MATS_FUNC(mats_c4, 5, 4)
DEFINE_MATS_FUNC(mats_c5, 5, 5)
DEFINE_MATS_FUNC(mats_c6, 5, 6)
DEFINE_MATS_FUNC(mats_c7, 5, 7)
DEFINE_MATS_FUNC(mats_c8, 5, 8)

/* ============================================================================
 * GATS: Geary Autocorrelations (lags 1-8 only)
 * ============================================================================ */

/* GATSm: Mass-weighted (property index 0) */
DEFINE_GATS_FUNC(gats_m1, 0, 1)
DEFINE_GATS_FUNC(gats_m2, 0, 2)
DEFINE_GATS_FUNC(gats_m3, 0, 3)
DEFINE_GATS_FUNC(gats_m4, 0, 4)
DEFINE_GATS_FUNC(gats_m5, 0, 5)
DEFINE_GATS_FUNC(gats_m6, 0, 6)
DEFINE_GATS_FUNC(gats_m7, 0, 7)
DEFINE_GATS_FUNC(gats_m8, 0, 8)

/* GATSv: Volume-weighted (property index 1) */
DEFINE_GATS_FUNC(gats_v1, 1, 1)
DEFINE_GATS_FUNC(gats_v2, 1, 2)
DEFINE_GATS_FUNC(gats_v3, 1, 3)
DEFINE_GATS_FUNC(gats_v4, 1, 4)
DEFINE_GATS_FUNC(gats_v5, 1, 5)
DEFINE_GATS_FUNC(gats_v6, 1, 6)
DEFINE_GATS_FUNC(gats_v7, 1, 7)
DEFINE_GATS_FUNC(gats_v8, 1, 8)

/* GATSe: Electronegativity-weighted (property index 2) */
DEFINE_GATS_FUNC(gats_e1, 2, 1)
DEFINE_GATS_FUNC(gats_e2, 2, 2)
DEFINE_GATS_FUNC(gats_e3, 2, 3)
DEFINE_GATS_FUNC(gats_e4, 2, 4)
DEFINE_GATS_FUNC(gats_e5, 2, 5)
DEFINE_GATS_FUNC(gats_e6, 2, 6)
DEFINE_GATS_FUNC(gats_e7, 2, 7)
DEFINE_GATS_FUNC(gats_e8, 2, 8)

/* GATSp: Polarizability-weighted (property index 3) */
DEFINE_GATS_FUNC(gats_p1, 3, 1)
DEFINE_GATS_FUNC(gats_p2, 3, 2)
DEFINE_GATS_FUNC(gats_p3, 3, 3)
DEFINE_GATS_FUNC(gats_p4, 3, 4)
DEFINE_GATS_FUNC(gats_p5, 3, 5)
DEFINE_GATS_FUNC(gats_p6, 3, 6)
DEFINE_GATS_FUNC(gats_p7, 3, 7)
DEFINE_GATS_FUNC(gats_p8, 3, 8)

/* GATSi: Ionization potential-weighted (property index 4) */
DEFINE_GATS_FUNC(gats_i1, 4, 1)
DEFINE_GATS_FUNC(gats_i2, 4, 2)
DEFINE_GATS_FUNC(gats_i3, 4, 3)
DEFINE_GATS_FUNC(gats_i4, 4, 4)
DEFINE_GATS_FUNC(gats_i5, 4, 5)
DEFINE_GATS_FUNC(gats_i6, 4, 6)
DEFINE_GATS_FUNC(gats_i7, 4, 7)
DEFINE_GATS_FUNC(gats_i8, 4, 8)

/* GATSc: Charge-weighted (property index 5) */
DEFINE_GATS_FUNC(gats_c1, 5, 1)
DEFINE_GATS_FUNC(gats_c2, 5, 2)
DEFINE_GATS_FUNC(gats_c3, 5, 3)
DEFINE_GATS_FUNC(gats_c4, 5, 4)
DEFINE_GATS_FUNC(gats_c5, 5, 5)
DEFINE_GATS_FUNC(gats_c6, 5, 6)
DEFINE_GATS_FUNC(gats_c7, 5, 7)
DEFINE_GATS_FUNC(gats_c8, 5, 8)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_AUTOCORR_DESC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_autocorrelations(void) {
    /* ATSm: Mass-weighted autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSm0", "Mass autocorr lag 0", desc_ats_m0);
    REGISTER_AUTOCORR_DESC("ATSm1", "Mass autocorr lag 1", desc_ats_m1);
    REGISTER_AUTOCORR_DESC("ATSm2", "Mass autocorr lag 2", desc_ats_m2);
    REGISTER_AUTOCORR_DESC("ATSm3", "Mass autocorr lag 3", desc_ats_m3);
    REGISTER_AUTOCORR_DESC("ATSm4", "Mass autocorr lag 4", desc_ats_m4);
    REGISTER_AUTOCORR_DESC("ATSm5", "Mass autocorr lag 5", desc_ats_m5);
    REGISTER_AUTOCORR_DESC("ATSm6", "Mass autocorr lag 6", desc_ats_m6);
    REGISTER_AUTOCORR_DESC("ATSm7", "Mass autocorr lag 7", desc_ats_m7);
    REGISTER_AUTOCORR_DESC("ATSm8", "Mass autocorr lag 8", desc_ats_m8);

    /* ATSv: Volume-weighted autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSv0", "Volume autocorr lag 0", desc_ats_v0);
    REGISTER_AUTOCORR_DESC("ATSv1", "Volume autocorr lag 1", desc_ats_v1);
    REGISTER_AUTOCORR_DESC("ATSv2", "Volume autocorr lag 2", desc_ats_v2);
    REGISTER_AUTOCORR_DESC("ATSv3", "Volume autocorr lag 3", desc_ats_v3);
    REGISTER_AUTOCORR_DESC("ATSv4", "Volume autocorr lag 4", desc_ats_v4);
    REGISTER_AUTOCORR_DESC("ATSv5", "Volume autocorr lag 5", desc_ats_v5);
    REGISTER_AUTOCORR_DESC("ATSv6", "Volume autocorr lag 6", desc_ats_v6);
    REGISTER_AUTOCORR_DESC("ATSv7", "Volume autocorr lag 7", desc_ats_v7);
    REGISTER_AUTOCORR_DESC("ATSv8", "Volume autocorr lag 8", desc_ats_v8);

    /* ATSe: Electronegativity-weighted autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSe0", "Electronegativity autocorr lag 0", desc_ats_e0);
    REGISTER_AUTOCORR_DESC("ATSe1", "Electronegativity autocorr lag 1", desc_ats_e1);
    REGISTER_AUTOCORR_DESC("ATSe2", "Electronegativity autocorr lag 2", desc_ats_e2);
    REGISTER_AUTOCORR_DESC("ATSe3", "Electronegativity autocorr lag 3", desc_ats_e3);
    REGISTER_AUTOCORR_DESC("ATSe4", "Electronegativity autocorr lag 4", desc_ats_e4);
    REGISTER_AUTOCORR_DESC("ATSe5", "Electronegativity autocorr lag 5", desc_ats_e5);
    REGISTER_AUTOCORR_DESC("ATSe6", "Electronegativity autocorr lag 6", desc_ats_e6);
    REGISTER_AUTOCORR_DESC("ATSe7", "Electronegativity autocorr lag 7", desc_ats_e7);
    REGISTER_AUTOCORR_DESC("ATSe8", "Electronegativity autocorr lag 8", desc_ats_e8);

    /* ATSp: Polarizability-weighted autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSp0", "Polarizability autocorr lag 0", desc_ats_p0);
    REGISTER_AUTOCORR_DESC("ATSp1", "Polarizability autocorr lag 1", desc_ats_p1);
    REGISTER_AUTOCORR_DESC("ATSp2", "Polarizability autocorr lag 2", desc_ats_p2);
    REGISTER_AUTOCORR_DESC("ATSp3", "Polarizability autocorr lag 3", desc_ats_p3);
    REGISTER_AUTOCORR_DESC("ATSp4", "Polarizability autocorr lag 4", desc_ats_p4);
    REGISTER_AUTOCORR_DESC("ATSp5", "Polarizability autocorr lag 5", desc_ats_p5);
    REGISTER_AUTOCORR_DESC("ATSp6", "Polarizability autocorr lag 6", desc_ats_p6);
    REGISTER_AUTOCORR_DESC("ATSp7", "Polarizability autocorr lag 7", desc_ats_p7);
    REGISTER_AUTOCORR_DESC("ATSp8", "Polarizability autocorr lag 8", desc_ats_p8);

    /* ATSi: Ionization potential-weighted autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSi0", "Ionization potential autocorr lag 0", desc_ats_i0);
    REGISTER_AUTOCORR_DESC("ATSi1", "Ionization potential autocorr lag 1", desc_ats_i1);
    REGISTER_AUTOCORR_DESC("ATSi2", "Ionization potential autocorr lag 2", desc_ats_i2);
    REGISTER_AUTOCORR_DESC("ATSi3", "Ionization potential autocorr lag 3", desc_ats_i3);
    REGISTER_AUTOCORR_DESC("ATSi4", "Ionization potential autocorr lag 4", desc_ats_i4);
    REGISTER_AUTOCORR_DESC("ATSi5", "Ionization potential autocorr lag 5", desc_ats_i5);
    REGISTER_AUTOCORR_DESC("ATSi6", "Ionization potential autocorr lag 6", desc_ats_i6);
    REGISTER_AUTOCORR_DESC("ATSi7", "Ionization potential autocorr lag 7", desc_ats_i7);
    REGISTER_AUTOCORR_DESC("ATSi8", "Ionization potential autocorr lag 8", desc_ats_i8);

    /* ATSc: Charge-weighted autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSc0", "Charge autocorr lag 0", desc_ats_c0);
    REGISTER_AUTOCORR_DESC("ATSc1", "Charge autocorr lag 1", desc_ats_c1);
    REGISTER_AUTOCORR_DESC("ATSc2", "Charge autocorr lag 2", desc_ats_c2);
    REGISTER_AUTOCORR_DESC("ATSc3", "Charge autocorr lag 3", desc_ats_c3);
    REGISTER_AUTOCORR_DESC("ATSc4", "Charge autocorr lag 4", desc_ats_c4);
    REGISTER_AUTOCORR_DESC("ATSc5", "Charge autocorr lag 5", desc_ats_c5);
    REGISTER_AUTOCORR_DESC("ATSc6", "Charge autocorr lag 6", desc_ats_c6);
    REGISTER_AUTOCORR_DESC("ATSc7", "Charge autocorr lag 7", desc_ats_c7);
    REGISTER_AUTOCORR_DESC("ATSc8", "Charge autocorr lag 8", desc_ats_c8);

    /* ========== AATS: Average Autocorrelations ========== */

    /* AATSm: Mass-weighted average autocorrelations */
    REGISTER_AUTOCORR_DESC("AATSm1", "Avg mass autocorr lag 1", desc_aats_m1);
    REGISTER_AUTOCORR_DESC("AATSm2", "Avg mass autocorr lag 2", desc_aats_m2);
    REGISTER_AUTOCORR_DESC("AATSm3", "Avg mass autocorr lag 3", desc_aats_m3);
    REGISTER_AUTOCORR_DESC("AATSm4", "Avg mass autocorr lag 4", desc_aats_m4);
    REGISTER_AUTOCORR_DESC("AATSm5", "Avg mass autocorr lag 5", desc_aats_m5);
    REGISTER_AUTOCORR_DESC("AATSm6", "Avg mass autocorr lag 6", desc_aats_m6);
    REGISTER_AUTOCORR_DESC("AATSm7", "Avg mass autocorr lag 7", desc_aats_m7);
    REGISTER_AUTOCORR_DESC("AATSm8", "Avg mass autocorr lag 8", desc_aats_m8);

    /* AATSv: Volume-weighted average autocorrelations */
    REGISTER_AUTOCORR_DESC("AATSv1", "Avg volume autocorr lag 1", desc_aats_v1);
    REGISTER_AUTOCORR_DESC("AATSv2", "Avg volume autocorr lag 2", desc_aats_v2);
    REGISTER_AUTOCORR_DESC("AATSv3", "Avg volume autocorr lag 3", desc_aats_v3);
    REGISTER_AUTOCORR_DESC("AATSv4", "Avg volume autocorr lag 4", desc_aats_v4);
    REGISTER_AUTOCORR_DESC("AATSv5", "Avg volume autocorr lag 5", desc_aats_v5);
    REGISTER_AUTOCORR_DESC("AATSv6", "Avg volume autocorr lag 6", desc_aats_v6);
    REGISTER_AUTOCORR_DESC("AATSv7", "Avg volume autocorr lag 7", desc_aats_v7);
    REGISTER_AUTOCORR_DESC("AATSv8", "Avg volume autocorr lag 8", desc_aats_v8);

    /* AATSe: Electronegativity-weighted average autocorrelations */
    REGISTER_AUTOCORR_DESC("AATSe1", "Avg EN autocorr lag 1", desc_aats_e1);
    REGISTER_AUTOCORR_DESC("AATSe2", "Avg EN autocorr lag 2", desc_aats_e2);
    REGISTER_AUTOCORR_DESC("AATSe3", "Avg EN autocorr lag 3", desc_aats_e3);
    REGISTER_AUTOCORR_DESC("AATSe4", "Avg EN autocorr lag 4", desc_aats_e4);
    REGISTER_AUTOCORR_DESC("AATSe5", "Avg EN autocorr lag 5", desc_aats_e5);
    REGISTER_AUTOCORR_DESC("AATSe6", "Avg EN autocorr lag 6", desc_aats_e6);
    REGISTER_AUTOCORR_DESC("AATSe7", "Avg EN autocorr lag 7", desc_aats_e7);
    REGISTER_AUTOCORR_DESC("AATSe8", "Avg EN autocorr lag 8", desc_aats_e8);

    /* AATSp: Polarizability-weighted average autocorrelations */
    REGISTER_AUTOCORR_DESC("AATSp1", "Avg polariz autocorr lag 1", desc_aats_p1);
    REGISTER_AUTOCORR_DESC("AATSp2", "Avg polariz autocorr lag 2", desc_aats_p2);
    REGISTER_AUTOCORR_DESC("AATSp3", "Avg polariz autocorr lag 3", desc_aats_p3);
    REGISTER_AUTOCORR_DESC("AATSp4", "Avg polariz autocorr lag 4", desc_aats_p4);
    REGISTER_AUTOCORR_DESC("AATSp5", "Avg polariz autocorr lag 5", desc_aats_p5);
    REGISTER_AUTOCORR_DESC("AATSp6", "Avg polariz autocorr lag 6", desc_aats_p6);
    REGISTER_AUTOCORR_DESC("AATSp7", "Avg polariz autocorr lag 7", desc_aats_p7);
    REGISTER_AUTOCORR_DESC("AATSp8", "Avg polariz autocorr lag 8", desc_aats_p8);

    /* AATSi: Ionization potential-weighted average autocorrelations */
    REGISTER_AUTOCORR_DESC("AATSi1", "Avg IP autocorr lag 1", desc_aats_i1);
    REGISTER_AUTOCORR_DESC("AATSi2", "Avg IP autocorr lag 2", desc_aats_i2);
    REGISTER_AUTOCORR_DESC("AATSi3", "Avg IP autocorr lag 3", desc_aats_i3);
    REGISTER_AUTOCORR_DESC("AATSi4", "Avg IP autocorr lag 4", desc_aats_i4);
    REGISTER_AUTOCORR_DESC("AATSi5", "Avg IP autocorr lag 5", desc_aats_i5);
    REGISTER_AUTOCORR_DESC("AATSi6", "Avg IP autocorr lag 6", desc_aats_i6);
    REGISTER_AUTOCORR_DESC("AATSi7", "Avg IP autocorr lag 7", desc_aats_i7);
    REGISTER_AUTOCORR_DESC("AATSi8", "Avg IP autocorr lag 8", desc_aats_i8);

    /* AATSc: Charge-weighted average autocorrelations */
    REGISTER_AUTOCORR_DESC("AATSc1", "Avg charge autocorr lag 1", desc_aats_c1);
    REGISTER_AUTOCORR_DESC("AATSc2", "Avg charge autocorr lag 2", desc_aats_c2);
    REGISTER_AUTOCORR_DESC("AATSc3", "Avg charge autocorr lag 3", desc_aats_c3);
    REGISTER_AUTOCORR_DESC("AATSc4", "Avg charge autocorr lag 4", desc_aats_c4);
    REGISTER_AUTOCORR_DESC("AATSc5", "Avg charge autocorr lag 5", desc_aats_c5);
    REGISTER_AUTOCORR_DESC("AATSc6", "Avg charge autocorr lag 6", desc_aats_c6);
    REGISTER_AUTOCORR_DESC("AATSc7", "Avg charge autocorr lag 7", desc_aats_c7);
    REGISTER_AUTOCORR_DESC("AATSc8", "Avg charge autocorr lag 8", desc_aats_c8);

    /* ========== ATSC: Centered Autocorrelations ========== */

    /* ATSCm: Mass-weighted centered autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSCm1", "Centered mass autocorr lag 1", desc_atsc_m1);
    REGISTER_AUTOCORR_DESC("ATSCm2", "Centered mass autocorr lag 2", desc_atsc_m2);
    REGISTER_AUTOCORR_DESC("ATSCm3", "Centered mass autocorr lag 3", desc_atsc_m3);
    REGISTER_AUTOCORR_DESC("ATSCm4", "Centered mass autocorr lag 4", desc_atsc_m4);
    REGISTER_AUTOCORR_DESC("ATSCm5", "Centered mass autocorr lag 5", desc_atsc_m5);
    REGISTER_AUTOCORR_DESC("ATSCm6", "Centered mass autocorr lag 6", desc_atsc_m6);
    REGISTER_AUTOCORR_DESC("ATSCm7", "Centered mass autocorr lag 7", desc_atsc_m7);
    REGISTER_AUTOCORR_DESC("ATSCm8", "Centered mass autocorr lag 8", desc_atsc_m8);

    /* ATSCv: Volume-weighted centered autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSCv1", "Centered volume autocorr lag 1", desc_atsc_v1);
    REGISTER_AUTOCORR_DESC("ATSCv2", "Centered volume autocorr lag 2", desc_atsc_v2);
    REGISTER_AUTOCORR_DESC("ATSCv3", "Centered volume autocorr lag 3", desc_atsc_v3);
    REGISTER_AUTOCORR_DESC("ATSCv4", "Centered volume autocorr lag 4", desc_atsc_v4);
    REGISTER_AUTOCORR_DESC("ATSCv5", "Centered volume autocorr lag 5", desc_atsc_v5);
    REGISTER_AUTOCORR_DESC("ATSCv6", "Centered volume autocorr lag 6", desc_atsc_v6);
    REGISTER_AUTOCORR_DESC("ATSCv7", "Centered volume autocorr lag 7", desc_atsc_v7);
    REGISTER_AUTOCORR_DESC("ATSCv8", "Centered volume autocorr lag 8", desc_atsc_v8);

    /* ATSCe: Electronegativity-weighted centered autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSCe1", "Centered EN autocorr lag 1", desc_atsc_e1);
    REGISTER_AUTOCORR_DESC("ATSCe2", "Centered EN autocorr lag 2", desc_atsc_e2);
    REGISTER_AUTOCORR_DESC("ATSCe3", "Centered EN autocorr lag 3", desc_atsc_e3);
    REGISTER_AUTOCORR_DESC("ATSCe4", "Centered EN autocorr lag 4", desc_atsc_e4);
    REGISTER_AUTOCORR_DESC("ATSCe5", "Centered EN autocorr lag 5", desc_atsc_e5);
    REGISTER_AUTOCORR_DESC("ATSCe6", "Centered EN autocorr lag 6", desc_atsc_e6);
    REGISTER_AUTOCORR_DESC("ATSCe7", "Centered EN autocorr lag 7", desc_atsc_e7);
    REGISTER_AUTOCORR_DESC("ATSCe8", "Centered EN autocorr lag 8", desc_atsc_e8);

    /* ATSCp: Polarizability-weighted centered autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSCp1", "Centered polariz autocorr lag 1", desc_atsc_p1);
    REGISTER_AUTOCORR_DESC("ATSCp2", "Centered polariz autocorr lag 2", desc_atsc_p2);
    REGISTER_AUTOCORR_DESC("ATSCp3", "Centered polariz autocorr lag 3", desc_atsc_p3);
    REGISTER_AUTOCORR_DESC("ATSCp4", "Centered polariz autocorr lag 4", desc_atsc_p4);
    REGISTER_AUTOCORR_DESC("ATSCp5", "Centered polariz autocorr lag 5", desc_atsc_p5);
    REGISTER_AUTOCORR_DESC("ATSCp6", "Centered polariz autocorr lag 6", desc_atsc_p6);
    REGISTER_AUTOCORR_DESC("ATSCp7", "Centered polariz autocorr lag 7", desc_atsc_p7);
    REGISTER_AUTOCORR_DESC("ATSCp8", "Centered polariz autocorr lag 8", desc_atsc_p8);

    /* ATSCi: Ionization potential-weighted centered autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSCi1", "Centered IP autocorr lag 1", desc_atsc_i1);
    REGISTER_AUTOCORR_DESC("ATSCi2", "Centered IP autocorr lag 2", desc_atsc_i2);
    REGISTER_AUTOCORR_DESC("ATSCi3", "Centered IP autocorr lag 3", desc_atsc_i3);
    REGISTER_AUTOCORR_DESC("ATSCi4", "Centered IP autocorr lag 4", desc_atsc_i4);
    REGISTER_AUTOCORR_DESC("ATSCi5", "Centered IP autocorr lag 5", desc_atsc_i5);
    REGISTER_AUTOCORR_DESC("ATSCi6", "Centered IP autocorr lag 6", desc_atsc_i6);
    REGISTER_AUTOCORR_DESC("ATSCi7", "Centered IP autocorr lag 7", desc_atsc_i7);
    REGISTER_AUTOCORR_DESC("ATSCi8", "Centered IP autocorr lag 8", desc_atsc_i8);

    /* ATSCc: Charge-weighted centered autocorrelations */
    REGISTER_AUTOCORR_DESC("ATSCc1", "Centered charge autocorr lag 1", desc_atsc_c1);
    REGISTER_AUTOCORR_DESC("ATSCc2", "Centered charge autocorr lag 2", desc_atsc_c2);
    REGISTER_AUTOCORR_DESC("ATSCc3", "Centered charge autocorr lag 3", desc_atsc_c3);
    REGISTER_AUTOCORR_DESC("ATSCc4", "Centered charge autocorr lag 4", desc_atsc_c4);
    REGISTER_AUTOCORR_DESC("ATSCc5", "Centered charge autocorr lag 5", desc_atsc_c5);
    REGISTER_AUTOCORR_DESC("ATSCc6", "Centered charge autocorr lag 6", desc_atsc_c6);
    REGISTER_AUTOCORR_DESC("ATSCc7", "Centered charge autocorr lag 7", desc_atsc_c7);
    REGISTER_AUTOCORR_DESC("ATSCc8", "Centered charge autocorr lag 8", desc_atsc_c8);

    /* ========== MATS: Moran Autocorrelations ========== */

    /* MATSm: Mass-weighted Moran autocorrelations */
    REGISTER_AUTOCORR_DESC("MATSm1", "Moran mass autocorr lag 1", desc_mats_m1);
    REGISTER_AUTOCORR_DESC("MATSm2", "Moran mass autocorr lag 2", desc_mats_m2);
    REGISTER_AUTOCORR_DESC("MATSm3", "Moran mass autocorr lag 3", desc_mats_m3);
    REGISTER_AUTOCORR_DESC("MATSm4", "Moran mass autocorr lag 4", desc_mats_m4);
    REGISTER_AUTOCORR_DESC("MATSm5", "Moran mass autocorr lag 5", desc_mats_m5);
    REGISTER_AUTOCORR_DESC("MATSm6", "Moran mass autocorr lag 6", desc_mats_m6);
    REGISTER_AUTOCORR_DESC("MATSm7", "Moran mass autocorr lag 7", desc_mats_m7);
    REGISTER_AUTOCORR_DESC("MATSm8", "Moran mass autocorr lag 8", desc_mats_m8);

    /* MATSv: Volume-weighted Moran autocorrelations */
    REGISTER_AUTOCORR_DESC("MATSv1", "Moran volume autocorr lag 1", desc_mats_v1);
    REGISTER_AUTOCORR_DESC("MATSv2", "Moran volume autocorr lag 2", desc_mats_v2);
    REGISTER_AUTOCORR_DESC("MATSv3", "Moran volume autocorr lag 3", desc_mats_v3);
    REGISTER_AUTOCORR_DESC("MATSv4", "Moran volume autocorr lag 4", desc_mats_v4);
    REGISTER_AUTOCORR_DESC("MATSv5", "Moran volume autocorr lag 5", desc_mats_v5);
    REGISTER_AUTOCORR_DESC("MATSv6", "Moran volume autocorr lag 6", desc_mats_v6);
    REGISTER_AUTOCORR_DESC("MATSv7", "Moran volume autocorr lag 7", desc_mats_v7);
    REGISTER_AUTOCORR_DESC("MATSv8", "Moran volume autocorr lag 8", desc_mats_v8);

    /* MATSe: Electronegativity-weighted Moran autocorrelations */
    REGISTER_AUTOCORR_DESC("MATSe1", "Moran EN autocorr lag 1", desc_mats_e1);
    REGISTER_AUTOCORR_DESC("MATSe2", "Moran EN autocorr lag 2", desc_mats_e2);
    REGISTER_AUTOCORR_DESC("MATSe3", "Moran EN autocorr lag 3", desc_mats_e3);
    REGISTER_AUTOCORR_DESC("MATSe4", "Moran EN autocorr lag 4", desc_mats_e4);
    REGISTER_AUTOCORR_DESC("MATSe5", "Moran EN autocorr lag 5", desc_mats_e5);
    REGISTER_AUTOCORR_DESC("MATSe6", "Moran EN autocorr lag 6", desc_mats_e6);
    REGISTER_AUTOCORR_DESC("MATSe7", "Moran EN autocorr lag 7", desc_mats_e7);
    REGISTER_AUTOCORR_DESC("MATSe8", "Moran EN autocorr lag 8", desc_mats_e8);

    /* MATSp: Polarizability-weighted Moran autocorrelations */
    REGISTER_AUTOCORR_DESC("MATSp1", "Moran polariz autocorr lag 1", desc_mats_p1);
    REGISTER_AUTOCORR_DESC("MATSp2", "Moran polariz autocorr lag 2", desc_mats_p2);
    REGISTER_AUTOCORR_DESC("MATSp3", "Moran polariz autocorr lag 3", desc_mats_p3);
    REGISTER_AUTOCORR_DESC("MATSp4", "Moran polariz autocorr lag 4", desc_mats_p4);
    REGISTER_AUTOCORR_DESC("MATSp5", "Moran polariz autocorr lag 5", desc_mats_p5);
    REGISTER_AUTOCORR_DESC("MATSp6", "Moran polariz autocorr lag 6", desc_mats_p6);
    REGISTER_AUTOCORR_DESC("MATSp7", "Moran polariz autocorr lag 7", desc_mats_p7);
    REGISTER_AUTOCORR_DESC("MATSp8", "Moran polariz autocorr lag 8", desc_mats_p8);

    /* MATSi: Ionization potential-weighted Moran autocorrelations */
    REGISTER_AUTOCORR_DESC("MATSi1", "Moran IP autocorr lag 1", desc_mats_i1);
    REGISTER_AUTOCORR_DESC("MATSi2", "Moran IP autocorr lag 2", desc_mats_i2);
    REGISTER_AUTOCORR_DESC("MATSi3", "Moran IP autocorr lag 3", desc_mats_i3);
    REGISTER_AUTOCORR_DESC("MATSi4", "Moran IP autocorr lag 4", desc_mats_i4);
    REGISTER_AUTOCORR_DESC("MATSi5", "Moran IP autocorr lag 5", desc_mats_i5);
    REGISTER_AUTOCORR_DESC("MATSi6", "Moran IP autocorr lag 6", desc_mats_i6);
    REGISTER_AUTOCORR_DESC("MATSi7", "Moran IP autocorr lag 7", desc_mats_i7);
    REGISTER_AUTOCORR_DESC("MATSi8", "Moran IP autocorr lag 8", desc_mats_i8);

    /* MATSc: Charge-weighted Moran autocorrelations */
    REGISTER_AUTOCORR_DESC("MATSc1", "Moran charge autocorr lag 1", desc_mats_c1);
    REGISTER_AUTOCORR_DESC("MATSc2", "Moran charge autocorr lag 2", desc_mats_c2);
    REGISTER_AUTOCORR_DESC("MATSc3", "Moran charge autocorr lag 3", desc_mats_c3);
    REGISTER_AUTOCORR_DESC("MATSc4", "Moran charge autocorr lag 4", desc_mats_c4);
    REGISTER_AUTOCORR_DESC("MATSc5", "Moran charge autocorr lag 5", desc_mats_c5);
    REGISTER_AUTOCORR_DESC("MATSc6", "Moran charge autocorr lag 6", desc_mats_c6);
    REGISTER_AUTOCORR_DESC("MATSc7", "Moran charge autocorr lag 7", desc_mats_c7);
    REGISTER_AUTOCORR_DESC("MATSc8", "Moran charge autocorr lag 8", desc_mats_c8);

    /* ========== GATS: Geary Autocorrelations ========== */

    /* GATSm: Mass-weighted Geary autocorrelations */
    REGISTER_AUTOCORR_DESC("GATSm1", "Geary mass autocorr lag 1", desc_gats_m1);
    REGISTER_AUTOCORR_DESC("GATSm2", "Geary mass autocorr lag 2", desc_gats_m2);
    REGISTER_AUTOCORR_DESC("GATSm3", "Geary mass autocorr lag 3", desc_gats_m3);
    REGISTER_AUTOCORR_DESC("GATSm4", "Geary mass autocorr lag 4", desc_gats_m4);
    REGISTER_AUTOCORR_DESC("GATSm5", "Geary mass autocorr lag 5", desc_gats_m5);
    REGISTER_AUTOCORR_DESC("GATSm6", "Geary mass autocorr lag 6", desc_gats_m6);
    REGISTER_AUTOCORR_DESC("GATSm7", "Geary mass autocorr lag 7", desc_gats_m7);
    REGISTER_AUTOCORR_DESC("GATSm8", "Geary mass autocorr lag 8", desc_gats_m8);

    /* GATSv: Volume-weighted Geary autocorrelations */
    REGISTER_AUTOCORR_DESC("GATSv1", "Geary volume autocorr lag 1", desc_gats_v1);
    REGISTER_AUTOCORR_DESC("GATSv2", "Geary volume autocorr lag 2", desc_gats_v2);
    REGISTER_AUTOCORR_DESC("GATSv3", "Geary volume autocorr lag 3", desc_gats_v3);
    REGISTER_AUTOCORR_DESC("GATSv4", "Geary volume autocorr lag 4", desc_gats_v4);
    REGISTER_AUTOCORR_DESC("GATSv5", "Geary volume autocorr lag 5", desc_gats_v5);
    REGISTER_AUTOCORR_DESC("GATSv6", "Geary volume autocorr lag 6", desc_gats_v6);
    REGISTER_AUTOCORR_DESC("GATSv7", "Geary volume autocorr lag 7", desc_gats_v7);
    REGISTER_AUTOCORR_DESC("GATSv8", "Geary volume autocorr lag 8", desc_gats_v8);

    /* GATSe: Electronegativity-weighted Geary autocorrelations */
    REGISTER_AUTOCORR_DESC("GATSe1", "Geary EN autocorr lag 1", desc_gats_e1);
    REGISTER_AUTOCORR_DESC("GATSe2", "Geary EN autocorr lag 2", desc_gats_e2);
    REGISTER_AUTOCORR_DESC("GATSe3", "Geary EN autocorr lag 3", desc_gats_e3);
    REGISTER_AUTOCORR_DESC("GATSe4", "Geary EN autocorr lag 4", desc_gats_e4);
    REGISTER_AUTOCORR_DESC("GATSe5", "Geary EN autocorr lag 5", desc_gats_e5);
    REGISTER_AUTOCORR_DESC("GATSe6", "Geary EN autocorr lag 6", desc_gats_e6);
    REGISTER_AUTOCORR_DESC("GATSe7", "Geary EN autocorr lag 7", desc_gats_e7);
    REGISTER_AUTOCORR_DESC("GATSe8", "Geary EN autocorr lag 8", desc_gats_e8);

    /* GATSp: Polarizability-weighted Geary autocorrelations */
    REGISTER_AUTOCORR_DESC("GATSp1", "Geary polariz autocorr lag 1", desc_gats_p1);
    REGISTER_AUTOCORR_DESC("GATSp2", "Geary polariz autocorr lag 2", desc_gats_p2);
    REGISTER_AUTOCORR_DESC("GATSp3", "Geary polariz autocorr lag 3", desc_gats_p3);
    REGISTER_AUTOCORR_DESC("GATSp4", "Geary polariz autocorr lag 4", desc_gats_p4);
    REGISTER_AUTOCORR_DESC("GATSp5", "Geary polariz autocorr lag 5", desc_gats_p5);
    REGISTER_AUTOCORR_DESC("GATSp6", "Geary polariz autocorr lag 6", desc_gats_p6);
    REGISTER_AUTOCORR_DESC("GATSp7", "Geary polariz autocorr lag 7", desc_gats_p7);
    REGISTER_AUTOCORR_DESC("GATSp8", "Geary polariz autocorr lag 8", desc_gats_p8);

    /* GATSi: Ionization potential-weighted Geary autocorrelations */
    REGISTER_AUTOCORR_DESC("GATSi1", "Geary IP autocorr lag 1", desc_gats_i1);
    REGISTER_AUTOCORR_DESC("GATSi2", "Geary IP autocorr lag 2", desc_gats_i2);
    REGISTER_AUTOCORR_DESC("GATSi3", "Geary IP autocorr lag 3", desc_gats_i3);
    REGISTER_AUTOCORR_DESC("GATSi4", "Geary IP autocorr lag 4", desc_gats_i4);
    REGISTER_AUTOCORR_DESC("GATSi5", "Geary IP autocorr lag 5", desc_gats_i5);
    REGISTER_AUTOCORR_DESC("GATSi6", "Geary IP autocorr lag 6", desc_gats_i6);
    REGISTER_AUTOCORR_DESC("GATSi7", "Geary IP autocorr lag 7", desc_gats_i7);
    REGISTER_AUTOCORR_DESC("GATSi8", "Geary IP autocorr lag 8", desc_gats_i8);

    /* GATSc: Charge-weighted Geary autocorrelations */
    REGISTER_AUTOCORR_DESC("GATSc1", "Geary charge autocorr lag 1", desc_gats_c1);
    REGISTER_AUTOCORR_DESC("GATSc2", "Geary charge autocorr lag 2", desc_gats_c2);
    REGISTER_AUTOCORR_DESC("GATSc3", "Geary charge autocorr lag 3", desc_gats_c3);
    REGISTER_AUTOCORR_DESC("GATSc4", "Geary charge autocorr lag 4", desc_gats_c4);
    REGISTER_AUTOCORR_DESC("GATSc5", "Geary charge autocorr lag 5", desc_gats_c5);
    REGISTER_AUTOCORR_DESC("GATSc6", "Geary charge autocorr lag 6", desc_gats_c6);
    REGISTER_AUTOCORR_DESC("GATSc7", "Geary charge autocorr lag 7", desc_gats_c7);
    REGISTER_AUTOCORR_DESC("GATSc8", "Geary charge autocorr lag 8", desc_gats_c8);
}
