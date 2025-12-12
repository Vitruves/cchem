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

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/descriptors.h"
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

/* Total descriptors: 6 properties × 9 lags (0-8) */
#define NUM_AUTOCORR_DESCRIPTORS 54

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

    /* Compute autocorrelations for each property and lag */
    int idx = 0;

    /* ATSm: Mass-weighted (lags 0-8) */
    for (int d = 0; d <= MAX_LAG; d++) {
        values[idx++].d = compute_ats(props_mass, dist_matrix, actual_heavy, d);
    }

    /* ATSv: Volume-weighted (lags 0-8) */
    for (int d = 0; d <= MAX_LAG; d++) {
        values[idx++].d = compute_ats(props_vol, dist_matrix, actual_heavy, d);
    }

    /* ATSe: Electronegativity-weighted (lags 0-8) */
    for (int d = 0; d <= MAX_LAG; d++) {
        values[idx++].d = compute_ats(props_en, dist_matrix, actual_heavy, d);
    }

    /* ATSp: Polarizability-weighted (lags 0-8) */
    for (int d = 0; d <= MAX_LAG; d++) {
        values[idx++].d = compute_ats(props_pol, dist_matrix, actual_heavy, d);
    }

    /* ATSi: Ionization potential-weighted (lags 0-8) */
    for (int d = 0; d <= MAX_LAG; d++) {
        values[idx++].d = compute_ats(props_ip, dist_matrix, actual_heavy, d);
    }

    /* ATSc: Charge-weighted (lags 0-8) */
    for (int d = 0; d <= MAX_LAG; d++) {
        values[idx++].d = compute_ats(props_charge, dist_matrix, actual_heavy, d);
    }

    if (heap_alloc) {
        free(dist_matrix); free(heavy_idx); free(props_mass); free(props_vol);
        free(props_en); free(props_pol); free(props_ip); free(props_charge);
        free(gasteiger);
    }

    return NUM_AUTOCORR_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions
 * ============================================================================ */

/* Macro for generating individual ATS descriptor functions */
#define DEFINE_ATS_FUNC(name, prop_idx, lag) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    descriptor_value_t all_values[NUM_AUTOCORR_DESCRIPTORS]; \
    descriptors_compute_autocorr_all(mol, all_values); \
    value->d = all_values[(prop_idx) * 9 + (lag)].d; \
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
}
