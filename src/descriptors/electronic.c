/**
 * @file electronic.c
 * @brief Electronic molecular descriptors
 *
 * Ultra-fast electronic descriptors based on Gasteiger-Marsili partial charges,
 * Kier-Hall E-States, and electronegativity-weighted topological indices.
 *
 * Group A: Gasteiger-Marsili Partial Charges (PEOE)
 * Group B: Electronegativity Weighted Topology
 * Group C: Electrotopological States (E-States)
 * Group D: Polarizability and Refractivity
 * Group E: Electronic Connectivity & Dipoles
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"
#include "cchem/simd.h"

/* ============================================================================
 * Physical Constants
 * ============================================================================ */

/* Gasteiger-Marsili parameters: a, b, c for chi = a + b*q + c*q^2 */
typedef struct {
    double a;  /* Electronegativity at q=0 */
    double b;  /* Linear term */
    double c;  /* Quadratic term */
} gm_params_t;

static const gm_params_t GM_PARAMS[] = {
    [ELEM_H]  = {7.17,  6.24, -0.56},
    [ELEM_C]  = {7.98,  9.18,  1.88},
    [ELEM_N]  = {11.54, 10.82, 1.36},
    [ELEM_O]  = {14.18, 12.92, 1.39},
    [ELEM_F]  = {14.66, 13.85, 2.31},
    [ELEM_P]  = {8.90,  8.24,  0.96},
    [ELEM_S]  = {10.14, 9.13,  1.38},
    [ELEM_Cl] = {11.00, 9.69,  1.35},
    [ELEM_Br] = {10.08, 8.47,  1.16},
    [ELEM_I]  = {9.90,  7.96,  0.96},
    [ELEM_Si] = {7.30,  6.56,  0.70},
    [ELEM_B]  = {6.00,  5.50,  0.50},
    [ELEM_Se] = {10.00, 8.50,  1.20},
    [ELEM_As] = {9.00,  7.80,  1.00},
};

/* Pauling electronegativities */
static const double ELECTRONEGATIVITY[] = {
    [ELEM_H]  = 2.20,
    [ELEM_C]  = 2.55,
    [ELEM_N]  = 3.04,
    [ELEM_O]  = 3.44,
    [ELEM_F]  = 3.98,
    [ELEM_P]  = 2.19,
    [ELEM_S]  = 2.58,
    [ELEM_Cl] = 3.16,
    [ELEM_Br] = 2.96,
    [ELEM_I]  = 2.66,
    [ELEM_Si] = 1.90,
    [ELEM_B]  = 2.04,
    [ELEM_Se] = 2.55,
    [ELEM_As] = 2.18,
};

/* Atomic polarizabilities in Å³ */
static const double POLARIZABILITY[] = {
    [ELEM_H]  = 0.667,
    [ELEM_C]  = 1.76,
    [ELEM_N]  = 1.10,
    [ELEM_O]  = 0.802,
    [ELEM_F]  = 0.557,
    [ELEM_P]  = 3.63,
    [ELEM_S]  = 2.90,
    [ELEM_Cl] = 2.18,
    [ELEM_Br] = 3.05,
    [ELEM_I]  = 4.70,
    [ELEM_Si] = 5.38,
    [ELEM_B]  = 3.03,
    [ELEM_Se] = 3.77,
    [ELEM_As] = 4.31,
};

/* Molar refractivity in Å³ */
static const double MOLAR_REFRACTIVITY[] = {
    [ELEM_H]  = 1.057,
    [ELEM_C]  = 2.503,
    [ELEM_N]  = 2.262,
    [ELEM_O]  = 1.607,
    [ELEM_F]  = 0.997,
    [ELEM_P]  = 6.920,
    [ELEM_S]  = 7.365,
    [ELEM_Cl] = 5.853,
    [ELEM_Br] = 8.927,
    [ELEM_I]  = 13.900,
    [ELEM_Si] = 6.840,
    [ELEM_B]  = 3.710,
    [ELEM_Se] = 7.990,
    [ELEM_As] = 6.280,
};

/* Intrinsic state values (Kier-Hall) for E-State calculation */
static const double INTRINSIC_STATE[] = {
    [ELEM_H]  = 0.0,
    [ELEM_C]  = 2.0,
    [ELEM_N]  = 3.0,
    [ELEM_O]  = 4.0,
    [ELEM_F]  = 5.0,
    [ELEM_P]  = 2.0,
    [ELEM_S]  = 3.0,
    [ELEM_Cl] = 4.0,
    [ELEM_Br] = 4.0,
    [ELEM_I]  = 4.0,
    [ELEM_Si] = 2.0,
    [ELEM_B]  = 1.5,
    [ELEM_Se] = 3.0,
    [ELEM_As] = 2.5,
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline gm_params_t get_gm_params(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(GM_PARAMS)/sizeof(GM_PARAMS[0]))) {
        return GM_PARAMS[ELEM_C];
    }
    gm_params_t p = GM_PARAMS[elem];
    if (p.a == 0.0 && p.b == 0.0 && p.c == 0.0) {
        return GM_PARAMS[ELEM_C];
    }
    return p;
}

static inline double get_electronegativity(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ELECTRONEGATIVITY)/sizeof(ELECTRONEGATIVITY[0])))
        return 2.55;
    double chi = ELECTRONEGATIVITY[elem];
    return (chi > 0.0) ? chi : 2.55;
}

static inline double get_polarizability(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(POLARIZABILITY)/sizeof(POLARIZABILITY[0])))
        return 1.76;
    double alpha = POLARIZABILITY[elem];
    return (alpha > 0.0) ? alpha : 1.76;
}

static inline double get_molar_refractivity(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(MOLAR_REFRACTIVITY)/sizeof(MOLAR_REFRACTIVITY[0])))
        return 2.503;
    double mr = MOLAR_REFRACTIVITY[elem];
    return (mr > 0.0) ? mr : 2.503;
}

static inline double get_intrinsic_state(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(INTRINSIC_STATE)/sizeof(INTRINSIC_STATE[0])))
        return 2.0;
    return INTRINSIC_STATE[elem];
}

/* ============================================================================
 * Gasteiger-Marsili PEOE Charge Calculation
 * ============================================================================ */

#define MAX_GM_ATOMS 1024
#define GM_ITERATIONS 6
#define GM_DAMPING 0.5

static void compute_gasteiger_charges(const molecule_t* mol, double* charges) {
    if (mol->num_atoms > MAX_GM_ATOMS) {
        /* Fall back to simple electronegativity-based charges */
        double sum_chi = 0.0;
        for (int i = 0; i < mol->num_atoms; i++) {
            sum_chi += get_electronegativity(mol->atoms[i].element);
        }
        double mean_chi = sum_chi / mol->num_atoms;
        for (int i = 0; i < mol->num_atoms; i++) {
            charges[i] = (get_electronegativity(mol->atoms[i].element) - mean_chi) / 10.0;
        }
        return;
    }

    gm_params_t params[MAX_GM_ATOMS];

    /* Initialize */
    for (int i = 0; i < mol->num_atoms; i++) {
        charges[i] = 0.0;
        params[i] = get_gm_params(mol->atoms[i].element);
    }

    /* PEOE iterations */
    double damping = 1.0;
    for (int iter = 0; iter < GM_ITERATIONS; iter++) {
        damping *= GM_DAMPING;

        /* For each bond, transfer charge based on electronegativity difference */
        for (int b = 0; b < mol->num_bonds; b++) {
            int i = mol->bonds[b].atom1;
            int j = mol->bonds[b].atom2;

            double chi_i = params[i].a + params[i].b * charges[i] + params[i].c * charges[i] * charges[i];
            double chi_j = params[j].a + params[j].b * charges[j] + params[j].c * charges[j] * charges[j];

            double delta_chi = chi_j - chi_i;
            double transfer;

            if (delta_chi > 0) {
                /* Charge flows from i to j */
                transfer = delta_chi / (params[i].a + params[i].b + params[i].c);
            } else {
                /* Charge flows from j to i */
                transfer = delta_chi / (params[j].a + params[j].b + params[j].c);
            }

            transfer *= damping;
            charges[i] += transfer;
            charges[j] -= transfer;
        }
    }
}

/* ============================================================================
 * E-State Calculation
 * ============================================================================ */

static void compute_estates(const molecule_t* mol, double* estates) {
    /* Initialize with intrinsic state values */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Intrinsic state: I = (2/n)^2 * (delta_v + 1) / delta */
        /* Simplified: use element-based intrinsic state */
        double I = get_intrinsic_state(atom->element);

        /* Modify by valence state */
        int degree = atom->num_neighbors;
        if (degree > 0) {
            int h_count = atom->implicit_h_count;
            for (int j = 0; j < atom->num_neighbors; j++) {
                if (mol->atoms[atom->neighbors[j]].element == ELEM_H) h_count++;
            }
            double delta = (double)(degree + h_count);
            double delta_v = get_electronegativity(atom->element) - 2.0;
            I = ((delta_v + 1.0) / delta) * I;
        }

        estates[i] = I;
    }

    /* Perturbation from neighbors: E_i = I_i + sum((I_i - I_j) / r_ij^2) */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        double I_i = estates[i];
        double perturbation = 0.0;

        for (int j = 0; j < atom->num_neighbors; j++) {
            int neighbor = atom->neighbors[j];
            double I_j = get_intrinsic_state(mol->atoms[neighbor].element);
            /* r_ij = 1 for direct neighbors (topological distance) */
            perturbation += (I_i - I_j);
        }

        estates[i] = I_i + perturbation;
    }
}

/* ============================================================================
 * Batch Computation (thread-safe: stack-allocated arrays)
 * ============================================================================ */

#define NUM_ELECTRONIC_DESCRIPTORS 30

/* Batch compute ALL electronic descriptors - computes charges/estates once */
int descriptors_compute_electronic_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    int n = mol->num_atoms;
    int idx = 0;

    /* Compute Gasteiger charges once */
    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    /* Compute E-States once */
    double estates[MAX_GM_ATOMS];
    compute_estates(mol, estates);

    /* Pre-compute charge statistics */
    double sum_pos = 0.0, sum_neg = 0.0;
    double max_pos = 0.0, max_neg = 0.0;
    double sum_q = 0.0, sum_q_sq = 0.0;
    double sum_abs_q = 0.0;

    for (int i = 0; i < n; i++) {
        double q = charges[i];
        if (q > 0) { sum_pos += q; if (q > max_pos) max_pos = q; }
        else { sum_neg += q; if (q < max_neg) max_neg = q; }
        sum_q += q;
        sum_q_sq += q * q;
        sum_abs_q += fabs(q);
    }
    double mean_q = sum_q / n;
    double var_q = (sum_q_sq / n) - (mean_q * mean_q);
    double total_abs = sum_pos - sum_neg;

    /* Pre-compute E-State statistics using SIMD */
    double sum_e = simd_sum_double(estates, n);
    double sum_e_sq = simd_sum_sq_double(estates, n);
    double max_e = simd_max_double(estates, n);
    double min_e = simd_min_double(estates, n);
    double mean_e = sum_e / n;
    double var_e = (sum_e_sq / n) - (mean_e * mean_e);

    /* Pre-compute electronegativity/polarizability sums */
    double sum_chi = 0.0, sum_pol = 0.0, sum_mr = 0.0, sum_heavy_pol = 0.0;
    double sum_chi_path = 0.0, sum_chi_asym = 0.0;
    double sum_chi_branch = 0.0, sum_chi_ring = 0.0;
    double sum_bond_chi = 0.0;
    double sum_elec_chi = 0.0, sum_elec_topo = 0.0, sum_polar_surface = 0.0;
    double sum_chi_pos = 0.0;

    for (int i = 0; i < n; i++) {
        const atom_t* atom = &mol->atoms[i];
        element_t elem = atom->element;
        double chi = get_electronegativity(elem);
        double pol = get_polarizability(elem);
        double mr = get_molar_refractivity(elem);

        sum_chi += chi;
        sum_pol += pol;
        sum_mr += mr;
        sum_chi_pos += chi * i;

        if (elem != ELEM_H) {
            sum_heavy_pol += pol;

            /* Branching */
            int degree = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                if (mol->atoms[atom->neighbors[j]].element != ELEM_H) degree++;
            }
            if (degree >= 3) sum_chi_branch += chi * (degree - 2);
        }

        /* Ring index */
        if (atom->ring_count > 0) sum_chi_ring += chi * atom->ring_count;

        /* Electronic topo sum */
        if (atom->num_neighbors > 0) sum_elec_topo += chi / atom->num_neighbors;

        /* Polar surface proxy */
        if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S || elem == ELEM_P) {
            sum_polar_surface += chi * pol;
        }
    }

    /* Bond-based computations */
    for (int b = 0; b < mol->num_bonds; b++) {
        const atom_t* a1 = &mol->atoms[mol->bonds[b].atom1];
        const atom_t* a2 = &mol->atoms[mol->bonds[b].atom2];
        double chi1 = get_electronegativity(a1->element);
        double chi2 = get_electronegativity(a2->element);

        sum_chi_path += sqrt(chi1 * chi2);
        sum_chi_asym += fabs(chi1 - chi2);
        sum_bond_chi += (chi1 + chi2) / 2.0;

        int deg1 = a1->num_neighbors + a1->implicit_h_count;
        int deg2 = a2->num_neighbors + a2->implicit_h_count;
        if (deg1 > 0 && deg2 > 0) {
            sum_elec_chi += 1.0 / sqrt(deg1 * deg2) * sqrt(chi1 * chi2) / 2.55;
        }
    }

    /* Dipole proxy */
    double center_chi = (sum_chi > 0) ? sum_chi_pos / sum_chi : 0.0;
    double dipole_var = 0.0;
    for (int i = 0; i < n; i++) {
        double chi = get_electronegativity(mol->atoms[i].element);
        double diff = i - center_chi;
        dipole_var += chi * diff * diff;
    }
    dipole_var /= n;

    /* Group A: Gasteiger-Marsili Partial Charges (8) */
    values[idx++].d = sum_pos;                                              /* TotalPositiveCharge */
    values[idx++].d = sum_neg;                                              /* TotalNegativeCharge */
    values[idx++].d = max_pos;                                              /* MaxPositiveCharge */
    values[idx++].d = max_neg;                                              /* MaxNegativeCharge */
    values[idx++].d = max_pos - max_neg;                                    /* ChargeRangePEOE */
    values[idx++].d = sum_abs_q / n;                                        /* MeanAbsCharge */
    values[idx++].d = var_q;                                                /* ChargeVariance */
    values[idx++].d = (total_abs > 0.001) ? sum_pos / total_abs : 0.5;     /* RelativePositiveCharge */

    /* Group B: Electronegativity Weighted Topology (6) */
    values[idx++].d = sum_chi;                                              /* SumChi */
    values[idx++].d = sum_chi_path;                                         /* ChiPathSum */
    values[idx++].d = sum_chi_asym;                                         /* ChiAsymmetry */
    values[idx++].d = sum_chi_branch;                                       /* ChiBranching */
    values[idx++].d = sum_chi_ring;                                         /* ChiRingIndex */
    values[idx++].d = (mol->num_bonds > 0) ? sum_bond_chi / mol->num_bonds : 0.0; /* MeanBondChi */

    /* Group C: Electrotopological States (6) */
    values[idx++].d = sum_e;                                                /* SumEState */
    values[idx++].d = max_e;                                                /* MaxEState */
    values[idx++].d = min_e;                                                /* MinEState */
    values[idx++].d = max_e - min_e;                                        /* EStateRange */
    values[idx++].d = mean_e;                                               /* MeanEState */
    values[idx++].d = var_e;                                                /* EStateVariance */

    /* Group D: Polarizability and Refractivity (5) */
    values[idx++].d = sum_pol;                                              /* TotalPolarizability */
    values[idx++].d = sum_mr;                                               /* TotalMolarRefractivity */
    values[idx++].d = sum_pol / n;                                          /* MeanPolarizability */
    values[idx++].d = (sum_mr > 0.001) ? sum_pol / sum_mr : 0.0;           /* PolMRRatio */
    values[idx++].d = sum_heavy_pol;                                        /* HeavyPolarizability */

    /* Group E: Electronic Connectivity & Dipoles (5) */
    values[idx++].d = sum_elec_chi;                                         /* ElectronicChi */
    values[idx++].d = dipole_var;                                           /* DipoleProxy */
    values[idx++].d = sum_elec_topo;                                        /* ElectronicTopoSum */
    values[idx++].d = sum_polar_surface;                                    /* PolarSurfaceProxy */
    values[idx++].d = max_pos * fabs(max_neg);                             /* ChargeSeparation */

    return NUM_ELECTRONIC_DESCRIPTORS;
}

/* ============================================================================
 * Group A: Gasteiger-Marsili Partial Charges (8 descriptors)
 * ============================================================================ */

/* 1. Total Positive Charge */
static cchem_status_t desc_total_positive_charge(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (charges[i] > 0) sum += charges[i];
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 2. Total Negative Charge */
static cchem_status_t desc_total_negative_charge(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (charges[i] < 0) sum += charges[i];
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 3. Maximum Positive Charge */
static cchem_status_t desc_max_positive_charge(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double max_pos = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (charges[i] > max_pos) max_pos = charges[i];
    }
    value->d = max_pos;
    return CCHEM_OK;
}

/* 4. Maximum Negative Charge */
static cchem_status_t desc_max_negative_charge(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double max_neg = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (charges[i] < max_neg) max_neg = charges[i];
    }
    value->d = max_neg;
    return CCHEM_OK;
}

/* 5. Charge Range */
static cchem_status_t desc_charge_range(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double max_q = charges[0], min_q = charges[0];
    for (int i = 1; i < mol->num_atoms; i++) {
        if (charges[i] > max_q) max_q = charges[i];
        if (charges[i] < min_q) min_q = charges[i];
    }
    value->d = max_q - min_q;
    return CCHEM_OK;
}

/* 6. Mean Absolute Charge */
static cchem_status_t desc_mean_abs_charge(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += fabs(charges[i]);
    }
    value->d = sum / mol->num_atoms;
    return CCHEM_OK;
}

/* 7. Charge Variance */
static cchem_status_t desc_charge_variance(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double sum = 0.0, sum_sq = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += charges[i];
        sum_sq += charges[i] * charges[i];
    }
    double mean = sum / mol->num_atoms;
    value->d = (sum_sq / mol->num_atoms) - (mean * mean);
    return CCHEM_OK;
}

/* 8. Relative Positive Charge */
static cchem_status_t desc_relative_positive_charge(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double pos = 0.0, neg = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (charges[i] > 0) pos += charges[i];
        else neg -= charges[i];
    }
    double total = pos + neg;
    value->d = (total > 0.001) ? pos / total : 0.5;
    return CCHEM_OK;
}

/* ============================================================================
 * Group B: Electronegativity Weighted Topology (6 descriptors)
 * ============================================================================ */

/* 9. Sum of Electronegativities */
static cchem_status_t desc_sum_chi(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += get_electronegativity(mol->atoms[i].element);
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 10. Chi-Weighted Path Sum (first-order connectivity) */
static cchem_status_t desc_chi_path_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int b = 0; b < mol->num_bonds; b++) {
        double chi_i = get_electronegativity(mol->atoms[mol->bonds[b].atom1].element);
        double chi_j = get_electronegativity(mol->atoms[mol->bonds[b].atom2].element);
        sum += sqrt(chi_i * chi_j);
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 11. Electronegativity Asymmetry: sum(|chi_i - chi_j|) for bonds */
static cchem_status_t desc_chi_asymmetry(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int b = 0; b < mol->num_bonds; b++) {
        double chi_i = get_electronegativity(mol->atoms[mol->bonds[b].atom1].element);
        double chi_j = get_electronegativity(mol->atoms[mol->bonds[b].atom2].element);
        sum += fabs(chi_i - chi_j);
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 12. Chi-Weighted Branching Index */
static cchem_status_t desc_chi_branching(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        int degree = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) degree++;
        }

        if (degree >= 3) {
            sum += get_electronegativity(atom->element) * (degree - 2);
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 13. Chi-Weighted Ring Index */
static cchem_status_t desc_chi_ring_index(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count > 0) {
            sum += get_electronegativity(atom->element) * atom->ring_count;
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 14. Mean Bond Electronegativity */
static cchem_status_t desc_mean_bond_chi(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    if (mol->num_bonds == 0) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    double sum = 0.0;
    for (int b = 0; b < mol->num_bonds; b++) {
        double chi_i = get_electronegativity(mol->atoms[mol->bonds[b].atom1].element);
        double chi_j = get_electronegativity(mol->atoms[mol->bonds[b].atom2].element);
        sum += (chi_i + chi_j) / 2.0;
    }
    value->d = sum / mol->num_bonds;
    return CCHEM_OK;
}

/* ============================================================================
 * Group C: Electrotopological States (6 descriptors)
 * ============================================================================ */

/* 15. Sum of E-States */
static cchem_status_t desc_sum_estate(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double estates[MAX_GM_ATOMS];
    compute_estates(mol, estates);

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += estates[i];
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 16. Maximum E-State */
static cchem_status_t desc_max_estate(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double estates[MAX_GM_ATOMS];
    compute_estates(mol, estates);

    double max_e = estates[0];
    for (int i = 1; i < mol->num_atoms; i++) {
        if (estates[i] > max_e) max_e = estates[i];
    }
    value->d = max_e;
    return CCHEM_OK;
}

/* 17. Minimum E-State */
static cchem_status_t desc_min_estate(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double estates[MAX_GM_ATOMS];
    compute_estates(mol, estates);

    double min_e = estates[0];
    for (int i = 1; i < mol->num_atoms; i++) {
        if (estates[i] < min_e) min_e = estates[i];
    }
    value->d = min_e;
    return CCHEM_OK;
}

/* 18. E-State Range */
static cchem_status_t desc_estate_range(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double estates[MAX_GM_ATOMS];
    compute_estates(mol, estates);

    double max_e = estates[0], min_e = estates[0];
    for (int i = 1; i < mol->num_atoms; i++) {
        if (estates[i] > max_e) max_e = estates[i];
        if (estates[i] < min_e) min_e = estates[i];
    }
    value->d = max_e - min_e;
    return CCHEM_OK;
}

/* 19. Mean E-State */
static cchem_status_t desc_mean_estate(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double estates[MAX_GM_ATOMS];
    compute_estates(mol, estates);

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += estates[i];
    }
    value->d = sum / mol->num_atoms;
    return CCHEM_OK;
}

/* 20. E-State Variance */
static cchem_status_t desc_estate_variance(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double estates[MAX_GM_ATOMS];
    compute_estates(mol, estates);

    double sum = 0.0, sum_sq = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += estates[i];
        sum_sq += estates[i] * estates[i];
    }
    double mean = sum / mol->num_atoms;
    value->d = (sum_sq / mol->num_atoms) - (mean * mean);
    return CCHEM_OK;
}

/* ============================================================================
 * Group D: Polarizability and Refractivity (5 descriptors)
 * ============================================================================ */

/* 21. Total Polarizability */
static cchem_status_t desc_total_polarizability(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += get_polarizability(mol->atoms[i].element);
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 22. Total Molar Refractivity */
static cchem_status_t desc_total_mr(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += get_molar_refractivity(mol->atoms[i].element);
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 23. Mean Polarizability */
static cchem_status_t desc_mean_polarizability(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum += get_polarizability(mol->atoms[i].element);
    }
    value->d = sum / mol->num_atoms;
    return CCHEM_OK;
}

/* 24. Polarizability/MR Ratio */
static cchem_status_t desc_pol_mr_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum_pol = 0.0, sum_mr = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        sum_pol += get_polarizability(mol->atoms[i].element);
        sum_mr += get_molar_refractivity(mol->atoms[i].element);
    }
    value->d = (sum_mr > 0.001) ? sum_pol / sum_mr : 0.0;
    return CCHEM_OK;
}

/* 25. Heavy Atom Polarizability */
static cchem_status_t desc_heavy_polarizability(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            sum += get_polarizability(mol->atoms[i].element);
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* ============================================================================
 * Group E: Electronic Connectivity & Dipoles (5 descriptors)
 * ============================================================================ */

/* 26. Electronic Connectivity Index (first-order chi) */
static cchem_status_t desc_electronic_chi(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    /* First-order connectivity index with electronegativity weighting */
    double sum = 0.0;
    for (int b = 0; b < mol->num_bonds; b++) {
        const atom_t* a1 = &mol->atoms[mol->bonds[b].atom1];
        const atom_t* a2 = &mol->atoms[mol->bonds[b].atom2];

        int deg1 = a1->num_neighbors + a1->implicit_h_count;
        int deg2 = a2->num_neighbors + a2->implicit_h_count;

        if (deg1 > 0 && deg2 > 0) {
            double chi1 = get_electronegativity(a1->element);
            double chi2 = get_electronegativity(a2->element);
            sum += 1.0 / sqrt(deg1 * deg2) * sqrt(chi1 * chi2) / 2.55;
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 27. Dipole Moment Proxy: sum(chi_i * position_i) variance */
static cchem_status_t desc_dipole_proxy(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    /* Use atom index as 1D position proxy */
    double sum_chi_pos = 0.0;
    double sum_chi = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        double chi = get_electronegativity(mol->atoms[i].element);
        sum_chi_pos += chi * i;
        sum_chi += chi;
    }

    /* Center of electronegativity */
    double center = (sum_chi > 0) ? sum_chi_pos / sum_chi : 0.0;

    /* Variance as dipole proxy */
    double variance = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double chi = get_electronegativity(mol->atoms[i].element);
        double diff = i - center;
        variance += chi * diff * diff;
    }

    value->d = variance / mol->num_atoms;
    return CCHEM_OK;
}

/* 28. Electronic Topological State Sum */
static cchem_status_t desc_electronic_topo_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        double chi = get_electronegativity(atom->element);
        int degree = atom->num_neighbors;
        if (degree > 0) {
            sum += chi / degree;
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 29. Polar Surface Proxy: chi-weighted polar atom sum */
static cchem_status_t desc_polar_surface_proxy(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S || elem == ELEM_P) {
            sum += get_electronegativity(elem) * get_polarizability(elem);
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* 30. Charge Separation Index: product of max_pos and |max_neg| */
static cchem_status_t desc_charge_separation(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double charges[MAX_GM_ATOMS];
    compute_gasteiger_charges(mol, charges);

    double max_pos = 0.0, max_neg = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (charges[i] > max_pos) max_pos = charges[i];
        if (charges[i] < max_neg) max_neg = charges[i];
    }
    value->d = max_pos * fabs(max_neg);
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_ELECTRONIC_DESC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_ELECTRONIC; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_electronic(void) {
    /* Group A: Gasteiger-Marsili Partial Charges */
    REGISTER_ELECTRONIC_DESC("TotalPositiveCharge", "Sum of positive PEOE charges", desc_total_positive_charge);
    REGISTER_ELECTRONIC_DESC("TotalNegativeCharge", "Sum of negative PEOE charges", desc_total_negative_charge);
    REGISTER_ELECTRONIC_DESC("MaxPositiveCharge", "Maximum positive PEOE charge", desc_max_positive_charge);
    REGISTER_ELECTRONIC_DESC("MaxNegativeCharge", "Maximum negative PEOE charge", desc_max_negative_charge);
    REGISTER_ELECTRONIC_DESC("ChargeRangePEOE", "Max - min PEOE charge", desc_charge_range);
    REGISTER_ELECTRONIC_DESC("MeanAbsCharge", "Mean absolute PEOE charge", desc_mean_abs_charge);
    REGISTER_ELECTRONIC_DESC("ChargeVariance", "Variance of PEOE charges", desc_charge_variance);
    REGISTER_ELECTRONIC_DESC("RelativePositiveCharge", "Positive charge fraction", desc_relative_positive_charge);

    /* Group B: Electronegativity Weighted Topology */
    REGISTER_ELECTRONIC_DESC("SumChi", "Sum of Pauling electronegativities", desc_sum_chi);
    REGISTER_ELECTRONIC_DESC("ChiPathSum", "Chi-weighted first-order path sum", desc_chi_path_sum);
    REGISTER_ELECTRONIC_DESC("ChiAsymmetry", "Bond electronegativity asymmetry", desc_chi_asymmetry);
    REGISTER_ELECTRONIC_DESC("ChiBranching", "Chi-weighted branching index", desc_chi_branching);
    REGISTER_ELECTRONIC_DESC("ChiRingIndex", "Chi-weighted ring membership", desc_chi_ring_index);
    REGISTER_ELECTRONIC_DESC("MeanBondChi", "Mean bond electronegativity", desc_mean_bond_chi);

    /* Group C: Electrotopological States */
    REGISTER_ELECTRONIC_DESC("SumEState", "Sum of E-State indices", desc_sum_estate);
    REGISTER_ELECTRONIC_DESC("MaxEState", "Maximum E-State index", desc_max_estate);
    REGISTER_ELECTRONIC_DESC("MinEState", "Minimum E-State index", desc_min_estate);
    REGISTER_ELECTRONIC_DESC("EStateRange", "E-State range", desc_estate_range);
    REGISTER_ELECTRONIC_DESC("MeanEState", "Mean E-State index", desc_mean_estate);
    REGISTER_ELECTRONIC_DESC("EStateVariance", "E-State variance", desc_estate_variance);

    /* Group D: Polarizability and Refractivity */
    REGISTER_ELECTRONIC_DESC("TotalPolarizability", "Sum of atomic polarizabilities", desc_total_polarizability);
    REGISTER_ELECTRONIC_DESC("TotalMolarRefractivity", "Sum of molar refractivity", desc_total_mr);
    REGISTER_ELECTRONIC_DESC("MeanPolarizability", "Mean atomic polarizability", desc_mean_polarizability);
    REGISTER_ELECTRONIC_DESC("PolMRRatio", "Polarizability/MR ratio", desc_pol_mr_ratio);
    REGISTER_ELECTRONIC_DESC("HeavyPolarizability", "Heavy atom polarizability", desc_heavy_polarizability);

    /* Group E: Electronic Connectivity & Dipoles */
    REGISTER_ELECTRONIC_DESC("ElectronicChi", "Electronic connectivity index", desc_electronic_chi);
    REGISTER_ELECTRONIC_DESC("DipoleProxy", "1D dipole moment proxy", desc_dipole_proxy);
    REGISTER_ELECTRONIC_DESC("ElectronicTopoSum", "Electronic topological sum", desc_electronic_topo_sum);
    REGISTER_ELECTRONIC_DESC("PolarSurfaceProxy", "Chi-weighted polar surface", desc_polar_surface_proxy);
    REGISTER_ELECTRONIC_DESC("ChargeSeparation", "Charge separation index", desc_charge_separation);
}
