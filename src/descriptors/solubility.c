/**
 * @file solubility.c
 * @brief CLogS - Calculated aqueous solubility descriptor
 *
 * Global linear model for aqueous solubility (LogS) prediction using
 * 24 molecular features including WCLogP, TPSA, atom counts, and
 * derived ratios.
 *
 * Performance (calibrated on ~20,000 compounds):
 *   - 5-Fold CV R² = 0.59
 *   - Training R² = 0.59
 *   - Pearson R = 0.77
 *   - MAE = 1.08 log units
 *
 * Note: R² > 0.8 requires nonlinear modeling (XGBoost achieves CV R² = 0.88)
 *
 * Reference: Extended ESOL-style model with Ridge regularization
 */

#include <string.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Model Coefficients (Ridge regression, alpha=1.0, trained on 19794 compounds)
 * ============================================================================ */

static const double SOL_INTERCEPT = -0.932933;

/* Primary feature coefficients */
static const double COEF_WCLOGP = -0.152894;
static const double COEF_HEAVYATOMCOUNT = -0.219451;
static const double COEF_ROTATABLEBONDCOUNT = -0.062109;
static const double COEF_AROMATICATOMCOUNT = -0.093821;
static const double COEF_CLASSICALTPSA = -0.002790;
static const double COEF_HBONDDONORCOUNT = 0.277804;
static const double COEF_HBONDACCEPTORCOUNT = 0.134443;
static const double COEF_SP3CARBONCOUNT = -0.200198;
static const double COEF_CARBONCOUNT = 0.422869;
static const double COEF_AROMATICRINGCOUNT = -0.398125;
static const double COEF_ALIPHATICRINGCOUNT = -0.273489;
static const double COEF_HALOGENCOUNT = 0.136953;
static const double COEF_NITROGENCOUNT = 0.042268;
static const double COEF_OXYGENCOUNT = 0.092175;
static const double COEF_SULFURCOUNT = 0.785413;

/* Derived feature coefficients */
static const double COEF_AROMATICRATIO = 4.301442;
static const double COEF_SP3RATIO = 3.092855;
static const double COEF_LOGP_SQ = -0.000763;
static const double COEF_TPSA_PER_HEAVY = -0.110240;
static const double COEF_HBD_RATIO = -0.914875;
static const double COEF_HBA_RATIO = 6.082462;
static const double COEF_LOGP_X_AROM = -0.814372;
static const double COEF_LOGP_X_SIZE = 0.002490;
static const double COEF_LOG_SIZE = -1.783780;

/* ============================================================================
 * WCLogP Atom Type Definitions (from logp.c)
 * ============================================================================ */

typedef enum {
    WC_C1, WC_C2, WC_C3, WC_C4, WC_C5, WC_C6, WC_C7, WC_C8, WC_C9, WC_C10,
    WC_C11, WC_C12, WC_C13, WC_C14, WC_C15, WC_C16, WC_C17, WC_C18, WC_C19, WC_C20,
    WC_C21, WC_C22, WC_C23, WC_C24, WC_C25, WC_C26, WC_C27,
    WC_H1, WC_H2, WC_H3, WC_H4,
    WC_N1, WC_N2, WC_N3, WC_N4, WC_N5, WC_N6, WC_N7, WC_N8, WC_N9, WC_N10,
    WC_N11, WC_N12, WC_N13, WC_N14,
    WC_O1, WC_O2, WC_O3, WC_O4, WC_O5, WC_O6, WC_O7, WC_O8, WC_O9, WC_O10,
    WC_O11, WC_O12,
    WC_S1, WC_S2, WC_S3,
    WC_P,
    WC_F, WC_Cl, WC_Br, WC_I,
    WC_OTHER,
    WC_NUM_TYPES
} wc_atom_type_t;

/* Wildman-Crippen LogP contributions */
static const double WC_LOGP_CONTRIB[WC_NUM_TYPES] = {
    0.0000, 1.0723, 1.1346, 1.0277, 0.8608, 0.9147, 0.9050, 0.8313, 0.7960, 0.5764,
    1.0112, 1.0344, 1.3594, 0.8700, 1.5099, 0.6656, 0.6256, 0.6611, 0.1747, 0.6085,
    0.4483, -0.4099, 0.2052, 0.8438, 0.7803, 0.8701, 0.0743,
    -0.3886, -0.4920, -1.8781, -0.5653,
    -0.1304, -0.0790, 0.2395, 0.1061, 0.4470, 0.3982, -0.3336, 0.2702, 0.8613, 1.4574,
    -0.9497, 0.0000, 0.0000, 0.0146,
    -0.4399, -0.1770, -0.1094, 0.0090, -0.6701, -0.4514, -0.4514, 0.1132, 0.0000, 0.0000,
    0.0000, 0.0000,
    1.5838, 0.5807, 0.2585,
    0.6808,
    -0.1986, 0.1285, 0.1537, -0.2238,
    0.0000
};

static const double WC_INTERCEPT = 1.3924;
static const double WC_COEF_RING_COUNT = -1.0047;
static const double WC_COEF_AROM_RING = -0.2689;
static const double WC_COEF_H_DONOR = 1.2144;
static const double WC_COEF_H_ACCEPTOR = 0.1247;
static const double WC_COEF_ROT_BONDS = 0.0539;

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static int count_heteroatom_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        element_t e = mol->atoms[atom->neighbors[i]].element;
        if (e != ELEM_C && e != ELEM_H) count++;
    }
    return count;
}

static int count_heavy_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) count++;
    }
    return count;
}

static int get_total_h(const molecule_t* mol, const atom_t* atom) {
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h++;
    }
    return h;
}

static bool has_double_bond_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE &&
            mol->atoms[atom->neighbors[i]].element == elem) {
            return true;
        }
    }
    return false;
}

static bool has_triple_bond(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_TRIPLE) return true;
    }
    return false;
}

static int count_double_bonds(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE) count++;
    }
    return count;
}

static bool is_ring_fusion(const molecule_t* mol, const atom_t* atom) {
    if (atom->ring_count < 2) return false;
    int arom_neighbors = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].aromatic) arom_neighbors++;
    }
    return (arom_neighbors >= 3);
}

static bool is_aromatic_substituted(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        const atom_t* nb = &mol->atoms[atom->neighbors[i]];
        if (!nb->aromatic && nb->element != ELEM_H && nb->element != ELEM_C) return true;
        if (nb->element == ELEM_C && !nb->aromatic) return true;
    }
    return false;
}

static bool is_rotatable_bond(const molecule_t* mol, const bond_t* bond) {
    if (bond->type != BOND_SINGLE) return false;
    if (bond->in_ring) return false;
    const atom_t* a1 = &mol->atoms[bond->atom1];
    const atom_t* a2 = &mol->atoms[bond->atom2];
    if (a1->element == ELEM_H || a2->element == ELEM_H) return false;
    int heavy1 = count_heavy_neighbors(mol, a1);
    int heavy2 = count_heavy_neighbors(mol, a2);
    return (heavy1 >= 2 && heavy2 >= 2);
}

/* ============================================================================
 * WCLogP Atom Type Assignment
 * ============================================================================ */

static wc_atom_type_t assign_carbon_type(const molecule_t* mol, const atom_t* atom) {
    int h_count = get_total_h(mol, atom);
    int heavy_neighbors = count_heavy_neighbors(mol, atom);
    int hetero_neighbors = count_heteroatom_neighbors(mol, atom);
    int double_bonds = count_double_bonds(mol, atom);
    bool has_triple = has_triple_bond(mol, atom);

    if (atom->aromatic) {
        if (is_ring_fusion(mol, atom)) return WC_C26;
        if (is_aromatic_substituted(mol, atom)) return WC_C25;
        return WC_C24;
    }
    if (has_double_bond_to(mol, atom, ELEM_O) || has_double_bond_to(mol, atom, ELEM_S)) {
        return WC_C27;
    }
    if (has_triple) return (h_count > 0) ? WC_C22 : WC_C23;
    if (double_bonds > 0) {
        if (hetero_neighbors == 0) {
            if (h_count == 2) return WC_C16;
            if (h_count == 1) return WC_C17;
            return WC_C18;
        } else {
            if (h_count == 1) return WC_C19;
            if (heavy_neighbors - hetero_neighbors > 0) return WC_C20;
            return WC_C21;
        }
    }
    if (hetero_neighbors == 0) {
        if (h_count >= 3) return WC_C2;
        if (h_count == 2) return WC_C3;
        if (h_count == 1) return WC_C4;
        return WC_C5;
    } else {
        int r_count = heavy_neighbors - hetero_neighbors;
        if (h_count == 3) return WC_C6;
        if (h_count == 2 && r_count >= 1) return WC_C7;
        if (h_count == 1 && r_count >= 2) return WC_C8;
        if (h_count == 0 && r_count >= 3) return WC_C9;
        if (h_count == 2) return WC_C10;
        if (h_count == 1 && r_count >= 1) return WC_C11;
        if (h_count == 0 && r_count >= 2) return WC_C12;
        if (h_count == 1) return WC_C13;
        if (h_count == 0 && r_count >= 1) return WC_C14;
        return WC_C15;
    }
}

static wc_atom_type_t assign_hydrogen_type(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        const atom_t* neighbor = &mol->atoms[atom->neighbors[i]];
        if (neighbor->element == ELEM_C) {
            if (neighbor->aromatic) return WC_H4;
            int hetero = count_heteroatom_neighbors(mol, neighbor);
            return (hetero > 0) ? WC_H2 : WC_H1;
        }
        return WC_H3;
    }
    return WC_H1;
}

static wc_atom_type_t assign_nitrogen_type(const molecule_t* mol, const atom_t* atom) {
    int h_count = get_total_h(mol, atom);
    int heavy_neighbors = count_heavy_neighbors(mol, atom);
    bool has_dbl = count_double_bonds(mol, atom) > 0;
    bool has_trpl = has_triple_bond(mol, atom);

    if (has_trpl) return WC_N11;
    if (atom->charge == 1 && has_double_bond_to(mol, atom, ELEM_O)) return WC_N10;
    if (atom->charge == 1 && heavy_neighbors == 4) return WC_N9;
    if (atom->aromatic) return (h_count > 0) ? WC_N5 : WC_N7;

    for (int i = 0; i < atom->num_neighbors; i++) {
        const atom_t* nb = &mol->atoms[atom->neighbors[i]];
        if (nb->element == ELEM_C && has_double_bond_to(mol, nb, ELEM_O)) return WC_N8;
    }

    bool bonded_arom = false;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].aromatic) { bonded_arom = true; break; }
    }
    if (bonded_arom) {
        if (h_count == 2) return WC_N4;
        if (h_count == 1) return WC_N5;
        return WC_N6;
    }
    if (h_count == 2) return WC_N1;
    if (h_count == 1) return WC_N2;
    if (h_count == 0 && !has_dbl) return WC_N3;
    return WC_N14;
}

static wc_atom_type_t assign_oxygen_type(const molecule_t* mol, const atom_t* atom) {
    int h_count = get_total_h(mol, atom);
    int heavy_neighbors = count_heavy_neighbors(mol, atom);
    bool has_dbl = count_double_bonds(mol, atom) > 0;

    if (atom->aromatic) return WC_O4;
    if (has_dbl && heavy_neighbors == 1) {
        for (int i = 0; i < atom->num_neighbors; i++) {
            const atom_t* nb = &mol->atoms[atom->neighbors[i]];
            if (nb->element == ELEM_C) {
                for (int j = 0; j < nb->num_neighbors; j++) {
                    const atom_t* nb2 = &mol->atoms[nb->neighbors[j]];
                    if (nb2 != atom && nb2->element == ELEM_O && get_total_h(mol, nb2) > 0)
                        return WC_O6;
                }
            }
        }
        return WC_O5;
    }
    if (h_count > 0) {
        for (int i = 0; i < atom->num_neighbors; i++) {
            const atom_t* nb = &mol->atoms[atom->neighbors[i]];
            if (nb->element == ELEM_C) {
                if (nb->aromatic) return WC_O2;
                if (has_double_bond_to(mol, nb, ELEM_O)) return WC_O7;
            }
        }
        return WC_O1;
    }
    if (heavy_neighbors == 2) {
        for (int i = 0; i < atom->num_neighbors; i++) {
            const atom_t* nb = &mol->atoms[atom->neighbors[i]];
            if (nb->element == ELEM_C && has_double_bond_to(mol, nb, ELEM_O)) return WC_O8;
        }
        return WC_O3;
    }
    return WC_O12;
}

static wc_atom_type_t assign_sulfur_type(const molecule_t* mol, const atom_t* atom) {
    int h_count = get_total_h(mol, atom);
    int double_bonds = count_double_bonds(mol, atom);
    if (double_bonds > 0) return WC_S3;
    if (h_count > 0) return WC_S1;
    return WC_S2;
}

static wc_atom_type_t assign_atom_type(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    switch (atom->element) {
        case ELEM_C:  return assign_carbon_type(mol, atom);
        case ELEM_H:  return assign_hydrogen_type(mol, atom);
        case ELEM_N:  return assign_nitrogen_type(mol, atom);
        case ELEM_O:  return assign_oxygen_type(mol, atom);
        case ELEM_S:  return assign_sulfur_type(mol, atom);
        case ELEM_P:  return WC_P;
        case ELEM_F:  return WC_F;
        case ELEM_Cl: return WC_Cl;
        case ELEM_Br: return WC_Br;
        case ELEM_I:  return WC_I;
        default:      return WC_OTHER;
    }
}

/* ============================================================================
 * TPSA Calculation (Ertl method)
 * ============================================================================ */

static double get_tpsa_contribution(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    element_t elem = atom->element;

    if (elem != ELEM_N && elem != ELEM_O && elem != ELEM_S && elem != ELEM_P) {
        return 0.0;
    }

    int h_count = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
    }
    int double_bonds = count_double_bonds(mol, atom);

    if (elem == ELEM_N) {
        if (atom->aromatic) return (h_count == 1) ? 15.79 : 12.89;
        if (atom->charge > 0) {
            if (h_count == 3) return 27.64;
            else if (h_count == 2) return 25.59;
            else if (h_count == 1) return 23.79;
            else return 16.61;
        }
        if (double_bonds > 0) return (h_count == 1) ? 23.85 : 12.36;
        if (h_count == 2) return 26.03;
        else if (h_count == 1) return 12.03;
        else return 3.24;
    }
    if (elem == ELEM_O) {
        if (atom->aromatic) return 13.14;
        if (atom->charge < 0) return 23.06;
        if (double_bonds > 0) return 17.07;
        if (h_count >= 1) return 20.23;
        return 9.23;
    }
    if (elem == ELEM_S) {
        if (atom->aromatic) return 28.24;
        if (double_bonds == 2) return 36.28;
        if (double_bonds == 1) return 32.09;
        if (h_count >= 1) return 38.80;
        return 25.30;
    }
    if (elem == ELEM_P) return 13.59;
    return 0.0;
}

/* ============================================================================
 * Feature Extraction and CLogS Calculation
 * ============================================================================ */

typedef struct {
    double wclogp;
    int heavy_atoms;
    int rotatable_bonds;
    int aromatic_atoms;
    double tpsa;
    int h_donors;
    int h_acceptors;
    int sp3_carbons;
    int carbons;
    int aromatic_rings;
    int aliphatic_rings;
    int halogens;
    int nitrogens;
    int oxygens;
    int sulfurs;
} mol_features_t;

static void extract_features(const molecule_t* mol, mol_features_t* f) {
    memset(f, 0, sizeof(mol_features_t));

    f->wclogp = WC_INTERCEPT;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        wc_atom_type_t type = assign_atom_type(mol, i);

        f->wclogp += WC_LOGP_CONTRIB[type];
        f->tpsa += get_tpsa_contribution(mol, i);

        if (atom->element != ELEM_H) {
            f->heavy_atoms++;
            if (atom->aromatic) f->aromatic_atoms++;

            /* H-bond donors: N or O with attached H */
            if ((atom->element == ELEM_O || atom->element == ELEM_N) &&
                get_total_h(mol, atom) > 0) {
                f->h_donors++;
            }

            /* H-bond acceptors: N or O */
            if (atom->element == ELEM_N || atom->element == ELEM_O) {
                f->h_acceptors++;
            }

            /* Element counts */
            switch (atom->element) {
                case ELEM_C:
                    f->carbons++;
                    if (!atom->aromatic && count_double_bonds(mol, atom) == 0 &&
                        !has_triple_bond(mol, atom)) {
                        f->sp3_carbons++;
                    }
                    break;
                case ELEM_N: f->nitrogens++; break;
                case ELEM_O: f->oxygens++; break;
                case ELEM_S: f->sulfurs++; break;
                case ELEM_F:
                case ELEM_Cl:
                case ELEM_Br:
                case ELEM_I:
                    f->halogens++;
                    break;
                default: break;
            }
        }
    }

    /* Add implicit H contributions */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_H && atom->implicit_h_count > 0) {
            wc_atom_type_t h_type = WC_H1;
            if (atom->element != ELEM_C) {
                h_type = WC_H3;
            } else if (atom->aromatic) {
                h_type = WC_H4;
            } else if (count_heteroatom_neighbors(mol, atom) > 0) {
                h_type = WC_H2;
            }
            f->wclogp += atom->implicit_h_count * WC_LOGP_CONTRIB[h_type];
        }
    }

    /* Rotatable bonds */
    for (int i = 0; i < mol->num_bonds; i++) {
        if (is_rotatable_bond(mol, &mol->bonds[i])) f->rotatable_bonds++;
    }

    /* Ring counts */
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].aromatic) f->aromatic_rings++;
        else f->aliphatic_rings++;
    }

    /* WCLogP corrections */
    f->wclogp += mol->num_rings * WC_COEF_RING_COUNT;
    f->wclogp += (f->aromatic_atoms / 5) * WC_COEF_AROM_RING;
    f->wclogp += f->h_donors * WC_COEF_H_DONOR;
    f->wclogp += f->h_acceptors * WC_COEF_H_ACCEPTOR;
    f->wclogp += f->rotatable_bonds * WC_COEF_ROT_BONDS;
}

static double calculate_clogs(const molecule_t* mol) {
    mol_features_t f;
    extract_features(mol, &f);

    /* Derived features */
    double aromatic_ratio = (f.heavy_atoms > 0) ?
        (double)f.aromatic_atoms / f.heavy_atoms : 0.0;
    double sp3_ratio = (f.carbons > 0) ?
        (double)f.sp3_carbons / f.carbons : 0.0;
    double logp_sq = f.wclogp * f.wclogp;
    double tpsa_per_heavy = (f.heavy_atoms > 0) ?
        f.tpsa / f.heavy_atoms : 0.0;
    double hbd_ratio = (f.heavy_atoms > 0) ?
        (double)f.h_donors / f.heavy_atoms : 0.0;
    double hba_ratio = (f.heavy_atoms > 0) ?
        (double)f.h_acceptors / f.heavy_atoms : 0.0;
    double logp_x_arom = f.wclogp * aromatic_ratio;
    double logp_x_size = f.wclogp * f.heavy_atoms;
    double log_size = log(1.0 + f.heavy_atoms);

    /* Calculate CLogS using linear model */
    double logs = SOL_INTERCEPT;

    /* Primary features */
    logs += COEF_WCLOGP * f.wclogp;
    logs += COEF_HEAVYATOMCOUNT * f.heavy_atoms;
    logs += COEF_ROTATABLEBONDCOUNT * f.rotatable_bonds;
    logs += COEF_AROMATICATOMCOUNT * f.aromatic_atoms;
    logs += COEF_CLASSICALTPSA * f.tpsa;
    logs += COEF_HBONDDONORCOUNT * f.h_donors;
    logs += COEF_HBONDACCEPTORCOUNT * f.h_acceptors;
    logs += COEF_SP3CARBONCOUNT * f.sp3_carbons;
    logs += COEF_CARBONCOUNT * f.carbons;
    logs += COEF_AROMATICRINGCOUNT * f.aromatic_rings;
    logs += COEF_ALIPHATICRINGCOUNT * f.aliphatic_rings;
    logs += COEF_HALOGENCOUNT * f.halogens;
    logs += COEF_NITROGENCOUNT * f.nitrogens;
    logs += COEF_OXYGENCOUNT * f.oxygens;
    logs += COEF_SULFURCOUNT * f.sulfurs;

    /* Derived features */
    logs += COEF_AROMATICRATIO * aromatic_ratio;
    logs += COEF_SP3RATIO * sp3_ratio;
    logs += COEF_LOGP_SQ * logp_sq;
    logs += COEF_TPSA_PER_HEAVY * tpsa_per_heavy;
    logs += COEF_HBD_RATIO * hbd_ratio;
    logs += COEF_HBA_RATIO * hba_ratio;
    logs += COEF_LOGP_X_AROM * logp_x_arom;
    logs += COEF_LOGP_X_SIZE * logp_x_size;
    logs += COEF_LOG_SIZE * log_size;

    return logs;
}

/* ============================================================================
 * Descriptor Implementation
 * ============================================================================ */

static cchem_status_t desc_clogs(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->d = calculate_clogs(mol);
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

void descriptors_register_solubility(void) {
    descriptor_def_t def = {0};
    strncpy(def.name, "CLogS", MAX_DESCRIPTOR_NAME - 1);
    strncpy(def.description, "Calculated aqueous solubility (log mol/L)",
            sizeof(def.description) - 1);
    def.category = DESC_CATEGORY_PROPERTIES;
    def.value_type = DESC_VALUE_DOUBLE;
    def.compute = desc_clogs;
    descriptor_register(&def);
}
