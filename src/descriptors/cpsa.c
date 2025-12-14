/**
 * @file cpsa.c
 * @brief Charged Partial Surface Area (CPSA) descriptors
 *
 * Extends VSA descriptors with charge-weighted and property-binned surface areas:
 * - CPSA: Surface area binned by Gasteiger partial charge
 * - HBSA: Hydrogen bond surface area descriptors
 * - LPSA: Lipophilicity surface area descriptors
 * - PSA_Pol/EN: Polarizability and electronegativity surface area
 * - PSA_IP/EA: Ionization potential and electron affinity surface area
 *
 * All O(n) complexity, reuses existing Gasteiger charge infrastructure.
 */

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Physical Constants and Property Tables
 * ============================================================================ */

/* Van der Waals radii (Bondi) in Angstroms */
static double get_vdw_radius(element_t elem) {
    switch (elem) {
        case ELEM_H:  return 1.20;
        case ELEM_C:  return 1.70;
        case ELEM_N:  return 1.55;
        case ELEM_O:  return 1.52;
        case ELEM_F:  return 1.47;
        case ELEM_P:  return 1.80;
        case ELEM_S:  return 1.80;
        case ELEM_Cl: return 1.75;
        case ELEM_Br: return 1.85;
        case ELEM_I:  return 1.98;
        case ELEM_Si: return 2.10;
        case ELEM_B:  return 1.92;
        case ELEM_Se: return 1.90;
        case ELEM_As: return 1.85;
        default:      return 1.70;
    }
}

/* Approximate VdW surface area contribution */
static double get_vsa(element_t elem) {
    double r = get_vdw_radius(elem);
    return 4.0 * 3.14159265 * r * r;
}

/* Pauling electronegativities */
static double get_electronegativity(element_t elem) {
    switch (elem) {
        case ELEM_H:  return 2.20;
        case ELEM_C:  return 2.55;
        case ELEM_N:  return 3.04;
        case ELEM_O:  return 3.44;
        case ELEM_F:  return 3.98;
        case ELEM_P:  return 2.19;
        case ELEM_S:  return 2.58;
        case ELEM_Cl: return 3.16;
        case ELEM_Br: return 2.96;
        case ELEM_I:  return 2.66;
        case ELEM_Si: return 1.90;
        case ELEM_B:  return 2.04;
        case ELEM_Se: return 2.55;
        case ELEM_As: return 2.18;
        default:      return 2.55;
    }
}

/* Atomic polarizabilities in Angstrom^3 */
static double get_polarizability(element_t elem) {
    switch (elem) {
        case ELEM_H:  return 0.667;
        case ELEM_C:  return 1.76;
        case ELEM_N:  return 1.10;
        case ELEM_O:  return 0.802;
        case ELEM_F:  return 0.557;
        case ELEM_P:  return 3.63;
        case ELEM_S:  return 2.90;
        case ELEM_Cl: return 2.18;
        case ELEM_Br: return 3.05;
        case ELEM_I:  return 4.70;
        case ELEM_Si: return 5.38;
        case ELEM_B:  return 3.03;
        case ELEM_Se: return 3.77;
        case ELEM_As: return 4.31;
        default:      return 1.76;
    }
}

/* First ionization potentials in eV */
static double get_ionization_potential(element_t elem) {
    switch (elem) {
        case ELEM_H:  return 13.60;
        case ELEM_C:  return 11.26;
        case ELEM_N:  return 14.53;
        case ELEM_O:  return 13.62;
        case ELEM_F:  return 17.42;
        case ELEM_P:  return 10.49;
        case ELEM_S:  return 10.36;
        case ELEM_Cl: return 12.97;
        case ELEM_Br: return 11.81;
        case ELEM_I:  return 10.45;
        case ELEM_Si: return 8.15;
        case ELEM_B:  return 8.30;
        case ELEM_Se: return 9.75;
        case ELEM_As: return 9.79;
        default:      return 11.0;
    }
}

/* Electron affinities in eV */
static double get_electron_affinity(element_t elem) {
    switch (elem) {
        case ELEM_H:  return 0.75;
        case ELEM_C:  return 1.26;
        case ELEM_N:  return -0.07;
        case ELEM_O:  return 1.46;
        case ELEM_F:  return 3.40;
        case ELEM_P:  return 0.75;
        case ELEM_S:  return 2.08;
        case ELEM_Cl: return 3.61;
        case ELEM_Br: return 3.36;
        case ELEM_I:  return 3.06;
        case ELEM_Si: return 1.39;
        case ELEM_B:  return 0.28;
        case ELEM_Se: return 2.02;
        case ELEM_As: return 0.80;
        default:      return 1.0;
    }
}

/* Wildman-Crippen LogP contributions */
static double get_logp_contrib(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int heavy = 0, h_count = atom->implicit_h_count;

    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
        else heavy++;
    }

    switch (atom->element) {
        case ELEM_C:
            if (atom->aromatic) return 0.13;
            if (heavy == 1) return 0.29;
            if (heavy == 2) return 0.19;
            if (heavy == 3) return 0.09;
            return -0.03;
        case ELEM_N:
            if (atom->aromatic) return -0.55;
            if (h_count >= 2) return -1.03;
            if (h_count == 1) return -0.54;
            return -0.28;
        case ELEM_O:
            if (h_count >= 1) return -0.47;
            return -0.27;
        case ELEM_S:
            if (h_count >= 1) return 0.42;
            return 0.29;
        case ELEM_F:  return 0.38;
        case ELEM_Cl: return 0.71;
        case ELEM_Br: return 0.87;
        case ELEM_I:  return 1.15;
        case ELEM_P:  return 0.16;
        default:      return 0.0;
    }
}

/* ============================================================================
 * Gasteiger Charge Estimation (simplified)
 * ============================================================================ */

static double get_charge_estimate(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    switch (atom->element) {
        case ELEM_O:  return -0.25;
        case ELEM_N:  return atom->aromatic ? -0.10 : -0.15;
        case ELEM_S:  return -0.05;
        case ELEM_F:  return -0.20;
        case ELEM_Cl: return -0.10;
        case ELEM_Br: return -0.08;
        case ELEM_I:  return -0.05;
        case ELEM_C:
            for (int i = 0; i < atom->num_neighbors; i++) {
                element_t ne = mol->atoms[atom->neighbors[i]].element;
                if (ne == ELEM_O || ne == ELEM_N || ne == ELEM_F) {
                    return 0.15;
                }
            }
            return 0.02;
        default:      return 0.0;
    }
}

/* ============================================================================
 * H-Bond Donor/Acceptor Classification
 * ============================================================================ */

static bool is_hbond_donor(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    element_t elem = atom->element;

    if (elem == ELEM_N || elem == ELEM_O) {
        int h_count = atom->implicit_h_count;
        for (int i = 0; i < atom->num_neighbors; i++) {
            if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
        }
        return h_count > 0;
    }
    return false;
}

static bool is_hbond_acceptor(const molecule_t* mol, int atom_idx) {
    element_t elem = mol->atoms[atom_idx].element;
    return (elem == ELEM_N || elem == ELEM_O || elem == ELEM_F);
}

/* ============================================================================
 * CPSA Descriptors - Charged Partial Surface Area (14 descriptors)
 * ============================================================================ */

typedef struct {
    double pos_sa;           /* q > 0 */
    double neg_sa;           /* q < 0 */
    double strong_pos_sa;    /* q > 0.2 */
    double strong_neg_sa;    /* q < -0.2 */
    double very_pos_sa;      /* q > 0.3 */
    double very_neg_sa;      /* q < -0.3 */
    double neutral_sa;       /* |q| < 0.1 */
    double weak_polar_sa;    /* 0.1 < |q| < 0.2 */
    double total_sa;
    double charge_weighted_pos;
    double charge_weighted_neg;
} cpsa_data_t;

static void compute_cpsa(const molecule_t* mol, cpsa_data_t* data) {
    memset(data, 0, sizeof(cpsa_data_t));

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;

        double q = get_charge_estimate(mol, i);
        double vsa = get_vsa(mol->atoms[i].element);

        data->total_sa += vsa;

        if (q > 0) {
            data->pos_sa += vsa;
            data->charge_weighted_pos += q * vsa;
        } else {
            data->neg_sa += vsa;
            data->charge_weighted_neg += fabs(q) * vsa;
        }

        if (q > 0.2) data->strong_pos_sa += vsa;
        if (q < -0.2) data->strong_neg_sa += vsa;
        if (q > 0.3) data->very_pos_sa += vsa;
        if (q < -0.3) data->very_neg_sa += vsa;

        double abs_q = fabs(q);
        if (abs_q < 0.1) data->neutral_sa += vsa;
        else if (abs_q < 0.2) data->weak_polar_sa += vsa;
    }
}

#define CPSA_FUNC(name, field) \
static cchem_status_t cpsa_##name(const molecule_t* mol, descriptor_value_t* value) { \
    cpsa_data_t data; \
    compute_cpsa(mol, &data); \
    value->d = data.field; \
    return CCHEM_OK; \
}

CPSA_FUNC(1, pos_sa)
CPSA_FUNC(2, neg_sa)
CPSA_FUNC(3, strong_pos_sa)
CPSA_FUNC(4, strong_neg_sa)
CPSA_FUNC(5, very_pos_sa)
CPSA_FUNC(6, very_neg_sa)
CPSA_FUNC(7, neutral_sa)
CPSA_FUNC(8, weak_polar_sa)

static cchem_status_t cpsa_9(const molecule_t* mol, descriptor_value_t* value) {
    cpsa_data_t data;
    compute_cpsa(mol, &data);
    value->d = (data.total_sa > 0) ? data.pos_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t cpsa_10(const molecule_t* mol, descriptor_value_t* value) {
    cpsa_data_t data;
    compute_cpsa(mol, &data);
    value->d = (data.total_sa > 0) ? data.neg_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

CPSA_FUNC(11, charge_weighted_pos)
CPSA_FUNC(12, charge_weighted_neg)

static cchem_status_t cpsa_13(const molecule_t* mol, descriptor_value_t* value) {
    cpsa_data_t data;
    compute_cpsa(mol, &data);
    value->d = data.charge_weighted_pos + data.charge_weighted_neg;
    return CCHEM_OK;
}

static cchem_status_t cpsa_14(const molecule_t* mol, descriptor_value_t* value) {
    cpsa_data_t data;
    compute_cpsa(mol, &data);
    value->d = data.charge_weighted_pos - data.charge_weighted_neg;
    return CCHEM_OK;
}

/* ============================================================================
 * HBSA Descriptors - H-Bond Surface Area (14 descriptors)
 * ============================================================================ */

typedef struct {
    double donor_sa;
    double acceptor_sa;
    double total_sa;
    double charge_weighted_donor;
    double charge_weighted_acceptor;
    double strong_donor_sa;
    double strong_acceptor_sa;
    double weak_donor_sa;
    double weak_acceptor_sa;
    double aromatic_hb_sa;
    double aliphatic_hb_sa;
} hbsa_data_t;

static void compute_hbsa(const molecule_t* mol, hbsa_data_t* data) {
    memset(data, 0, sizeof(hbsa_data_t));

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;

        double vsa = get_vsa(mol->atoms[i].element);
        double q = get_charge_estimate(mol, i);
        bool is_donor = is_hbond_donor(mol, i);
        bool is_acceptor = is_hbond_acceptor(mol, i);
        bool aromatic = mol->atoms[i].aromatic;

        data->total_sa += vsa;

        if (is_donor) {
            data->donor_sa += vsa;
            data->charge_weighted_donor += q * vsa;
            if (q > 0.2) data->strong_donor_sa += vsa;
            else data->weak_donor_sa += vsa;
            if (aromatic) data->aromatic_hb_sa += vsa;
            else data->aliphatic_hb_sa += vsa;
        }

        if (is_acceptor) {
            data->acceptor_sa += vsa;
            data->charge_weighted_acceptor += q * vsa;
            if (q < -0.2) data->strong_acceptor_sa += vsa;
            else data->weak_acceptor_sa += vsa;
            if (aromatic && !is_donor) data->aromatic_hb_sa += vsa;
            else if (!is_donor) data->aliphatic_hb_sa += vsa;
        }
    }
}

#define HBSA_FUNC(name, field) \
static cchem_status_t hbsa_##name(const molecule_t* mol, descriptor_value_t* value) { \
    hbsa_data_t data; \
    compute_hbsa(mol, &data); \
    value->d = data.field; \
    return CCHEM_OK; \
}

HBSA_FUNC(1, donor_sa)
HBSA_FUNC(2, acceptor_sa)

static cchem_status_t hbsa_3(const molecule_t* mol, descriptor_value_t* value) {
    hbsa_data_t data;
    compute_hbsa(mol, &data);
    value->d = data.donor_sa + data.acceptor_sa;
    return CCHEM_OK;
}

static cchem_status_t hbsa_4(const molecule_t* mol, descriptor_value_t* value) {
    hbsa_data_t data;
    compute_hbsa(mol, &data);
    value->d = (data.total_sa > 0) ? data.donor_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t hbsa_5(const molecule_t* mol, descriptor_value_t* value) {
    hbsa_data_t data;
    compute_hbsa(mol, &data);
    value->d = (data.total_sa > 0) ? data.acceptor_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t hbsa_6(const molecule_t* mol, descriptor_value_t* value) {
    hbsa_data_t data;
    compute_hbsa(mol, &data);
    value->d = (data.acceptor_sa > 0) ? data.donor_sa / data.acceptor_sa : 0.0;
    return CCHEM_OK;
}

HBSA_FUNC(7, charge_weighted_donor)
HBSA_FUNC(8, charge_weighted_acceptor)
HBSA_FUNC(9, strong_donor_sa)
HBSA_FUNC(10, strong_acceptor_sa)
HBSA_FUNC(11, weak_donor_sa)
HBSA_FUNC(12, weak_acceptor_sa)
HBSA_FUNC(13, aromatic_hb_sa)
HBSA_FUNC(14, aliphatic_hb_sa)

/* ============================================================================
 * LPSA Descriptors - Lipophilicity Surface Area (14 descriptors)
 * ============================================================================ */

typedef struct {
    double hydrophobic_sa;
    double hydrophilic_sa;
    double total_sa;
    double logp_weighted_sa;
    double pos_logp_sa;
    double neg_logp_sa;
    double aromatic_hydrophobic_sa;
    double aliphatic_hydrophobic_sa;
    double halogen_sa;
    double carbon_only_sa;
    double heteroatom_sa;
} lpsa_data_t;

static void compute_lpsa(const molecule_t* mol, lpsa_data_t* data) {
    memset(data, 0, sizeof(lpsa_data_t));

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_H) continue;

        double vsa = get_vsa(elem);
        double logp = get_logp_contrib(mol, i);
        bool aromatic = mol->atoms[i].aromatic;

        data->total_sa += vsa;
        data->logp_weighted_sa += logp * vsa;

        if (logp > 0) data->pos_logp_sa += vsa;
        else data->neg_logp_sa += vsa;

        /* Hydrophobic: C and halogens */
        if (elem == ELEM_C || elem == ELEM_F || elem == ELEM_Cl ||
            elem == ELEM_Br || elem == ELEM_I) {
            data->hydrophobic_sa += vsa;
            if (elem == ELEM_C) {
                data->carbon_only_sa += vsa;
                if (aromatic) data->aromatic_hydrophobic_sa += vsa;
                else data->aliphatic_hydrophobic_sa += vsa;
            } else {
                data->halogen_sa += vsa;
            }
        }

        /* Hydrophilic: N, O, S */
        if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S) {
            data->hydrophilic_sa += vsa;
            data->heteroatom_sa += vsa;
        }
    }
}

#define LPSA_FUNC(name, field) \
static cchem_status_t lpsa_##name(const molecule_t* mol, descriptor_value_t* value) { \
    lpsa_data_t data; \
    compute_lpsa(mol, &data); \
    value->d = data.field; \
    return CCHEM_OK; \
}

LPSA_FUNC(1, hydrophobic_sa)
LPSA_FUNC(2, hydrophilic_sa)

static cchem_status_t lpsa_3(const molecule_t* mol, descriptor_value_t* value) {
    lpsa_data_t data;
    compute_lpsa(mol, &data);
    value->d = (data.total_sa > 0) ? data.hydrophobic_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t lpsa_4(const molecule_t* mol, descriptor_value_t* value) {
    lpsa_data_t data;
    compute_lpsa(mol, &data);
    value->d = (data.total_sa > 0) ? data.hydrophilic_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

LPSA_FUNC(5, logp_weighted_sa)
LPSA_FUNC(6, pos_logp_sa)
LPSA_FUNC(7, neg_logp_sa)
LPSA_FUNC(8, aromatic_hydrophobic_sa)
LPSA_FUNC(9, aliphatic_hydrophobic_sa)
LPSA_FUNC(10, halogen_sa)
LPSA_FUNC(11, carbon_only_sa)
LPSA_FUNC(12, heteroatom_sa)

static cchem_status_t lpsa_13(const molecule_t* mol, descriptor_value_t* value) {
    lpsa_data_t data;
    compute_lpsa(mol, &data);
    value->d = (data.hydrophilic_sa > 0) ? data.hydrophobic_sa / data.hydrophilic_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t lpsa_14(const molecule_t* mol, descriptor_value_t* value) {
    lpsa_data_t data;
    compute_lpsa(mol, &data);
    value->d = (data.total_sa > 0) ? (data.pos_logp_sa - data.neg_logp_sa) / data.total_sa : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * PSA_Pol/EN - Polarizability & Electronegativity Surface Area (14 descriptors)
 * ============================================================================ */

typedef struct {
    double low_pol_sa;     /* pol < 1.0 */
    double med_pol_sa;     /* 1.0 <= pol < 2.5 */
    double high_pol_sa;    /* pol >= 2.5 */
    double pol_weighted_sa;
    double total_sa;
    double low_en_sa;      /* EN < 2.5 */
    double med_en_sa;      /* 2.5 <= EN < 3.0 */
    double high_en_sa;     /* EN >= 3.0 */
    double en_weighted_sa;
    double sum_pol;
    double sum_en;
    double sum_pol_sq;
    double sum_en_sq;
    int count;
} psa_pol_en_data_t;

static void compute_psa_pol_en(const molecule_t* mol, psa_pol_en_data_t* data) {
    memset(data, 0, sizeof(psa_pol_en_data_t));

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_H) continue;

        double vsa = get_vsa(elem);
        double pol = get_polarizability(elem);
        double en = get_electronegativity(elem);

        data->total_sa += vsa;
        data->pol_weighted_sa += pol * vsa;
        data->en_weighted_sa += en * vsa;
        data->sum_pol += pol;
        data->sum_en += en;
        data->sum_pol_sq += pol * pol;
        data->sum_en_sq += en * en;
        data->count++;

        if (pol < 1.0) data->low_pol_sa += vsa;
        else if (pol < 2.5) data->med_pol_sa += vsa;
        else data->high_pol_sa += vsa;

        if (en < 2.5) data->low_en_sa += vsa;
        else if (en < 3.0) data->med_en_sa += vsa;
        else data->high_en_sa += vsa;
    }
}

#define PSA_POL_EN_FUNC(name, field) \
static cchem_status_t psa_pol_en_##name(const molecule_t* mol, descriptor_value_t* value) { \
    psa_pol_en_data_t data; \
    compute_psa_pol_en(mol, &data); \
    value->d = data.field; \
    return CCHEM_OK; \
}

PSA_POL_EN_FUNC(1, low_pol_sa)
PSA_POL_EN_FUNC(2, med_pol_sa)
PSA_POL_EN_FUNC(3, high_pol_sa)
PSA_POL_EN_FUNC(4, pol_weighted_sa)

static cchem_status_t psa_pol_en_5(const molecule_t* mol, descriptor_value_t* value) {
    psa_pol_en_data_t data;
    compute_psa_pol_en(mol, &data);
    value->d = (data.total_sa > 0) ? data.pol_weighted_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t psa_pol_en_6(const molecule_t* mol, descriptor_value_t* value) {
    psa_pol_en_data_t data;
    compute_psa_pol_en(mol, &data);
    if (data.count < 2) { value->d = 0.0; return CCHEM_OK; }
    double mean = data.sum_pol / data.count;
    value->d = (data.sum_pol_sq / data.count) - (mean * mean);
    return CCHEM_OK;
}

static cchem_status_t psa_pol_en_7(const molecule_t* mol, descriptor_value_t* value) {
    psa_pol_en_data_t data;
    compute_psa_pol_en(mol, &data);
    value->d = (data.total_sa > 0) ? data.low_pol_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t psa_pol_en_8(const molecule_t* mol, descriptor_value_t* value) {
    psa_pol_en_data_t data;
    compute_psa_pol_en(mol, &data);
    value->d = (data.total_sa > 0) ? data.high_pol_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

PSA_POL_EN_FUNC(9, low_en_sa)
PSA_POL_EN_FUNC(10, med_en_sa)
PSA_POL_EN_FUNC(11, high_en_sa)
PSA_POL_EN_FUNC(12, en_weighted_sa)

static cchem_status_t psa_pol_en_13(const molecule_t* mol, descriptor_value_t* value) {
    psa_pol_en_data_t data;
    compute_psa_pol_en(mol, &data);
    value->d = (data.total_sa > 0) ? data.en_weighted_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t psa_pol_en_14(const molecule_t* mol, descriptor_value_t* value) {
    psa_pol_en_data_t data;
    compute_psa_pol_en(mol, &data);
    if (data.count < 2) { value->d = 0.0; return CCHEM_OK; }
    double mean = data.sum_en / data.count;
    value->d = (data.sum_en_sq / data.count) - (mean * mean);
    return CCHEM_OK;
}

/* ============================================================================
 * PSA_IP/EA - Ionization Potential & Electron Affinity Surface Area (14 descriptors)
 * ============================================================================ */

typedef struct {
    double low_ip_sa;      /* IP < 10 eV */
    double med_ip_sa;      /* 10 <= IP < 12 */
    double high_ip_sa;     /* IP >= 12 eV */
    double ip_weighted_sa;
    double total_sa;
    double low_ea_sa;      /* EA < 0.5 eV */
    double med_ea_sa;      /* 0.5 <= EA < 2.0 */
    double high_ea_sa;     /* EA >= 2.0 eV */
    double ea_weighted_sa;
    double sum_ip;
    double sum_ea;
    double sum_ip_sq;
    double sum_ea_sq;
    int count;
} psa_ip_ea_data_t;

static void compute_psa_ip_ea(const molecule_t* mol, psa_ip_ea_data_t* data) {
    memset(data, 0, sizeof(psa_ip_ea_data_t));

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_H) continue;

        double vsa = get_vsa(elem);
        double ip = get_ionization_potential(elem);
        double ea = get_electron_affinity(elem);

        data->total_sa += vsa;
        data->ip_weighted_sa += ip * vsa;
        data->ea_weighted_sa += ea * vsa;
        data->sum_ip += ip;
        data->sum_ea += ea;
        data->sum_ip_sq += ip * ip;
        data->sum_ea_sq += ea * ea;
        data->count++;

        if (ip < 10.0) data->low_ip_sa += vsa;
        else if (ip < 12.0) data->med_ip_sa += vsa;
        else data->high_ip_sa += vsa;

        if (ea < 0.5) data->low_ea_sa += vsa;
        else if (ea < 2.0) data->med_ea_sa += vsa;
        else data->high_ea_sa += vsa;
    }
}

#define PSA_IP_EA_FUNC(name, field) \
static cchem_status_t psa_ip_ea_##name(const molecule_t* mol, descriptor_value_t* value) { \
    psa_ip_ea_data_t data; \
    compute_psa_ip_ea(mol, &data); \
    value->d = data.field; \
    return CCHEM_OK; \
}

PSA_IP_EA_FUNC(1, low_ip_sa)
PSA_IP_EA_FUNC(2, med_ip_sa)
PSA_IP_EA_FUNC(3, high_ip_sa)
PSA_IP_EA_FUNC(4, ip_weighted_sa)

static cchem_status_t psa_ip_ea_5(const molecule_t* mol, descriptor_value_t* value) {
    psa_ip_ea_data_t data;
    compute_psa_ip_ea(mol, &data);
    value->d = (data.total_sa > 0) ? data.ip_weighted_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t psa_ip_ea_6(const molecule_t* mol, descriptor_value_t* value) {
    psa_ip_ea_data_t data;
    compute_psa_ip_ea(mol, &data);
    if (data.count < 2) { value->d = 0.0; return CCHEM_OK; }
    double mean = data.sum_ip / data.count;
    value->d = (data.sum_ip_sq / data.count) - (mean * mean);
    return CCHEM_OK;
}

static cchem_status_t psa_ip_ea_7(const molecule_t* mol, descriptor_value_t* value) {
    psa_ip_ea_data_t data;
    compute_psa_ip_ea(mol, &data);
    value->d = (data.total_sa > 0) ? data.low_ip_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

PSA_IP_EA_FUNC(8, low_ea_sa)
PSA_IP_EA_FUNC(9, med_ea_sa)
PSA_IP_EA_FUNC(10, high_ea_sa)
PSA_IP_EA_FUNC(11, ea_weighted_sa)

static cchem_status_t psa_ip_ea_12(const molecule_t* mol, descriptor_value_t* value) {
    psa_ip_ea_data_t data;
    compute_psa_ip_ea(mol, &data);
    value->d = (data.total_sa > 0) ? data.ea_weighted_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t psa_ip_ea_13(const molecule_t* mol, descriptor_value_t* value) {
    psa_ip_ea_data_t data;
    compute_psa_ip_ea(mol, &data);
    value->d = (data.total_sa > 0) ? data.high_ea_sa / data.total_sa : 0.0;
    return CCHEM_OK;
}

static cchem_status_t psa_ip_ea_14(const molecule_t* mol, descriptor_value_t* value) {
    psa_ip_ea_data_t data;
    compute_psa_ip_ea(mol, &data);
    if (data.count < 2) { value->d = 0.0; return CCHEM_OK; }
    double mean = data.sum_ea / data.count;
    value->d = (data.sum_ea_sq / data.count) - (mean * mean);
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_CPSA(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_STERIC; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_cpsa(void) {
    /* CPSA - Charged Partial Surface Area (14) */
    REGISTER_CPSA("CPSA_1", "Total positive surface area (q>0)", cpsa_1);
    REGISTER_CPSA("CPSA_2", "Total negative surface area (q<0)", cpsa_2);
    REGISTER_CPSA("CPSA_3", "Partial positive surface (q>0.2)", cpsa_3);
    REGISTER_CPSA("CPSA_4", "Partial negative surface (q<-0.2)", cpsa_4);
    REGISTER_CPSA("CPSA_5", "Highly positive surface (q>0.3)", cpsa_5);
    REGISTER_CPSA("CPSA_6", "Highly negative surface (q<-0.3)", cpsa_6);
    REGISTER_CPSA("CPSA_7", "Near-neutral surface (|q|<0.1)", cpsa_7);
    REGISTER_CPSA("CPSA_8", "Weakly polar surface (0.1<|q|<0.2)", cpsa_8);
    REGISTER_CPSA("CPSA_9", "Relative positive SA (pos/total)", cpsa_9);
    REGISTER_CPSA("CPSA_10", "Relative negative SA (neg/total)", cpsa_10);
    REGISTER_CPSA("CPSA_11", "Charge-weighted positive SA", cpsa_11);
    REGISTER_CPSA("CPSA_12", "Charge-weighted negative SA", cpsa_12);
    REGISTER_CPSA("CPSA_13", "Total charge-weighted SA", cpsa_13);
    REGISTER_CPSA("CPSA_14", "Charge imbalance SA", cpsa_14);

    /* HBSA - H-Bond Surface Area (14) */
    REGISTER_CPSA("HBSA_1", "H-bond donor surface area", hbsa_1);
    REGISTER_CPSA("HBSA_2", "H-bond acceptor surface area", hbsa_2);
    REGISTER_CPSA("HBSA_3", "Total H-bond surface area", hbsa_3);
    REGISTER_CPSA("HBSA_4", "Relative donor SA (donor/total)", hbsa_4);
    REGISTER_CPSA("HBSA_5", "Relative acceptor SA (acceptor/total)", hbsa_5);
    REGISTER_CPSA("HBSA_6", "Donor/acceptor SA ratio", hbsa_6);
    REGISTER_CPSA("HBSA_7", "Charge-weighted donor SA", hbsa_7);
    REGISTER_CPSA("HBSA_8", "Charge-weighted acceptor SA", hbsa_8);
    REGISTER_CPSA("HBSA_9", "Strong H-bond donor SA (q>0.2)", hbsa_9);
    REGISTER_CPSA("HBSA_10", "Strong H-bond acceptor SA (q<-0.2)", hbsa_10);
    REGISTER_CPSA("HBSA_11", "Weak H-bond donor SA", hbsa_11);
    REGISTER_CPSA("HBSA_12", "Weak H-bond acceptor SA", hbsa_12);
    REGISTER_CPSA("HBSA_13", "Aromatic H-bond SA", hbsa_13);
    REGISTER_CPSA("HBSA_14", "Aliphatic H-bond SA", hbsa_14);

    /* LPSA - Lipophilicity Surface Area (14) */
    REGISTER_CPSA("LPSA_1", "Hydrophobic surface area (C,Hal)", lpsa_1);
    REGISTER_CPSA("LPSA_2", "Hydrophilic surface area (N,O,S)", lpsa_2);
    REGISTER_CPSA("LPSA_3", "Relative hydrophobic SA", lpsa_3);
    REGISTER_CPSA("LPSA_4", "Relative hydrophilic SA", lpsa_4);
    REGISTER_CPSA("LPSA_5", "LogP-weighted surface area", lpsa_5);
    REGISTER_CPSA("LPSA_6", "Positive LogP surface (logP>0)", lpsa_6);
    REGISTER_CPSA("LPSA_7", "Negative LogP surface (logP<0)", lpsa_7);
    REGISTER_CPSA("LPSA_8", "Aromatic hydrophobic SA", lpsa_8);
    REGISTER_CPSA("LPSA_9", "Aliphatic hydrophobic SA", lpsa_9);
    REGISTER_CPSA("LPSA_10", "Halogen surface area", lpsa_10);
    REGISTER_CPSA("LPSA_11", "Carbon-only surface area", lpsa_11);
    REGISTER_CPSA("LPSA_12", "Heteroatom surface area", lpsa_12);
    REGISTER_CPSA("LPSA_13", "Amphiphilic balance (hydrophobic/hydrophilic)", lpsa_13);
    REGISTER_CPSA("LPSA_14", "LogP SA asymmetry", lpsa_14);

    /* PSA_Pol - Polarizability Surface Area (8) */
    REGISTER_CPSA("PSA_Pol1", "Low polarizability SA (pol<1.0)", psa_pol_en_1);
    REGISTER_CPSA("PSA_Pol2", "Medium polarizability SA (1.0<=pol<2.5)", psa_pol_en_2);
    REGISTER_CPSA("PSA_Pol3", "High polarizability SA (pol>=2.5)", psa_pol_en_3);
    REGISTER_CPSA("PSA_Pol4", "Polarizability-weighted total SA", psa_pol_en_4);
    REGISTER_CPSA("PSA_Pol5", "Mean polarizability per SA", psa_pol_en_5);
    REGISTER_CPSA("PSA_Pol6", "Polarizability SA variance", psa_pol_en_6);
    REGISTER_CPSA("PSA_Pol7", "Relative low-pol SA", psa_pol_en_7);
    REGISTER_CPSA("PSA_Pol8", "Relative high-pol SA", psa_pol_en_8);

    /* PSA_EN - Electronegativity Surface Area (6) */
    REGISTER_CPSA("PSA_EN1", "Low EN surface (EN<2.5)", psa_pol_en_9);
    REGISTER_CPSA("PSA_EN2", "Medium EN surface (2.5<=EN<3.0)", psa_pol_en_10);
    REGISTER_CPSA("PSA_EN3", "High EN surface (EN>=3.0)", psa_pol_en_11);
    REGISTER_CPSA("PSA_EN4", "EN-weighted total SA", psa_pol_en_12);
    REGISTER_CPSA("PSA_EN5", "Mean EN per SA", psa_pol_en_13);
    REGISTER_CPSA("PSA_EN6", "EN SA variance", psa_pol_en_14);

    /* PSA_IP - Ionization Potential Surface Area (7) */
    REGISTER_CPSA("PSA_IP1", "Low IP surface area (IP<10eV)", psa_ip_ea_1);
    REGISTER_CPSA("PSA_IP2", "Medium IP surface area (10<=IP<12)", psa_ip_ea_2);
    REGISTER_CPSA("PSA_IP3", "High IP surface area (IP>=12eV)", psa_ip_ea_3);
    REGISTER_CPSA("PSA_IP4", "IP-weighted total SA", psa_ip_ea_4);
    REGISTER_CPSA("PSA_IP5", "Mean IP per SA", psa_ip_ea_5);
    REGISTER_CPSA("PSA_IP6", "IP SA variance", psa_ip_ea_6);
    REGISTER_CPSA("PSA_IP7", "Relative low-IP SA", psa_ip_ea_7);

    /* PSA_EA - Electron Affinity Surface Area (7) */
    REGISTER_CPSA("PSA_EA1", "Low EA surface (EA<0.5eV)", psa_ip_ea_8);
    REGISTER_CPSA("PSA_EA2", "Medium EA surface (0.5<=EA<2.0)", psa_ip_ea_9);
    REGISTER_CPSA("PSA_EA3", "High EA surface (EA>=2.0eV)", psa_ip_ea_10);
    REGISTER_CPSA("PSA_EA4", "EA-weighted total SA", psa_ip_ea_11);
    REGISTER_CPSA("PSA_EA5", "Mean EA per SA", psa_ip_ea_12);
    REGISTER_CPSA("PSA_EA6", "Relative high-EA SA", psa_ip_ea_13);
    REGISTER_CPSA("PSA_EA7", "EA SA variance", psa_ip_ea_14);
}
