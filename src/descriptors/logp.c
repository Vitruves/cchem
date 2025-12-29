/**
 * @file logp.c
 * @brief Wildman-Crippen LogP descriptor
 *
 * Implementation of the Wildman-Crippen atom contribution method for LogP.
 * Reference: Wildman, S. A.; Crippen, G. M. J. Chem. Inf. Comput. Sci. 1999, 39, 868-873.
 *
 * Coefficients optimized via linear regression on 12,948 compounds
 * achieving R² = 0.82, Pearson R = 0.90, MAE = 0.60
 */

#include <string.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Wildman-Crippen Atom Type Definitions
 * 68 atom types with LogP contributions
 * ============================================================================ */

typedef enum {
    /* Carbon types (C1-C27) */
    WC_C1, WC_C2, WC_C3, WC_C4, WC_C5, WC_C6, WC_C7, WC_C8, WC_C9, WC_C10,
    WC_C11, WC_C12, WC_C13, WC_C14, WC_C15, WC_C16, WC_C17, WC_C18, WC_C19, WC_C20,
    WC_C21, WC_C22, WC_C23, WC_C24, WC_C25, WC_C26, WC_C27,
    /* Hydrogen types (H1-H4) */
    WC_H1, WC_H2, WC_H3, WC_H4,
    /* Nitrogen types (N1-N14) */
    WC_N1, WC_N2, WC_N3, WC_N4, WC_N5, WC_N6, WC_N7, WC_N8, WC_N9, WC_N10,
    WC_N11, WC_N12, WC_N13, WC_N14,
    /* Oxygen types (O1-O12) */
    WC_O1, WC_O2, WC_O3, WC_O4, WC_O5, WC_O6, WC_O7, WC_O8, WC_O9, WC_O10,
    WC_O11, WC_O12,
    /* Sulfur types (S1-S3) */
    WC_S1, WC_S2, WC_S3,
    /* Phosphorus */
    WC_P,
    /* Halogens */
    WC_F, WC_Cl, WC_Br, WC_I,
    /* Other */
    WC_OTHER,
    WC_NUM_TYPES
} wc_atom_type_t;

/* Optimized LogP contributions (linear regression on 12,948 compounds)
 * R² = 0.82, Pearson R = 0.90, MAE = 0.60
 * Values obtained via least squares regression with ridge regularization */
static const double WC_LOGP_CONTRIB[WC_NUM_TYPES] = {
    /* C1-C27 - Carbon types */
    0.0000,   /* C1:  CH4 (methane) - rare */
    1.0723,   /* C2:  CH3R (methyl) */
    1.1346,   /* C3:  CH2R2 (methylene) */
    1.0277,   /* C4:  CHR3 (methine) */
    0.8608,   /* C5:  CR4 (quaternary C) */
    0.9147,   /* C6:  CH3X (methyl with heteroatom) */
    0.9050,   /* C7:  CH2RX */
    0.8313,   /* C8:  CHR2X */
    0.7960,   /* C9:  CR3X */
    0.5764,   /* C10: CH2X2 */
    1.0112,   /* C11: CHRX2 */
    1.0344,   /* C12: CR2X2 */
    1.3594,   /* C13: CHX3 */
    0.8700,   /* C14: CRX3 */
    1.5099,   /* C15: CX4 */
    0.6656,   /* C16: =CH2 (terminal vinyl) */
    0.6256,   /* C17: =CHR */
    0.6611,   /* C18: =CR2 */
    0.1747,   /* C19: =CHX */
    0.6085,   /* C20: =CRX */
    0.4483,   /* C21: =CX2 */
    -0.4099,  /* C22: #CH (terminal alkyne) */
    0.2052,   /* C23: #CR */
    0.8438,   /* C24: aromatic C (unsubstituted) */
    0.7803,   /* C25: aromatic C (substituted) */
    0.8701,   /* C26: aromatic C (fused ring junction) */
    0.0743,   /* C27: C=O (carbonyl carbon) */
    /* H1-H4 - Hydrogen types */
    -0.3886,  /* H1:  H on C sp3 (aliphatic, no X) */
    -0.4920,  /* H2:  H on C sp3 (with X) */
    -1.8781,  /* H3:  H on heteroatom (OH, NH, SH) */
    -0.5653,  /* H4:  H on aromatic C */
    /* N1-N14 - Nitrogen types */
    -0.1304,  /* N1:  -NH2 (primary amine) */
    -0.0790,  /* N2:  >NH (secondary amine) */
    0.2395,   /* N3:  >N- (tertiary amine) */
    0.1061,   /* N4:  -NH2 (aniline) */
    0.4470,   /* N5:  >NH (secondary aromatic) */
    0.3982,   /* N6:  >N- (tertiary aromatic) */
    -0.3336,  /* N7:  =N- (pyridine-like) */
    0.2702,   /* N8:  amide N */
    0.8613,   /* N9:  N+ quaternary */
    1.4574,   /* N10: nitro N */
    -0.9497,  /* N11: nitrile N */
    0.0000,   /* N12: azide/diazo N */
    0.0000,   /* N13: N-oxide */
    0.0146,   /* N14: other N */
    /* O1-O12 - Oxygen types */
    -0.4399,  /* O1:  -OH (alcohol) */
    -0.1770,  /* O2:  -OH (phenol) */
    -0.1094,  /* O3:  -O- (ether, aliphatic) */
    0.0090,   /* O4:  -O- (furan, aromatic) */
    -0.6701,  /* O5:  =O (carbonyl, ketone/aldehyde) */
    -0.4514,  /* O6:  =O (carboxylic) */
    -0.4514,  /* O7:  -OH (carboxylic acid) */
    0.1132,   /* O8:  -O- (ester) */
    0.0000,   /* O9:  =O (ester) - rare */
    0.0000,   /* O10: nitro O */
    0.0000,   /* O11: N-oxide O */
    0.0000,   /* O12: other O */
    /* S1-S3 - Sulfur types */
    1.5838,   /* S1:  -SH (thiol) */
    0.5807,   /* S2:  -S- (sulfide/thioether) */
    0.2585,   /* S3:  >S=O, -SO2- (sulfoxide/sulfone) */
    /* P */
    0.6808,   /* P: any phosphorus */
    /* Halogens */
    -0.1986,  /* F */
    0.1285,   /* Cl */
    0.1537,   /* Br */
    -0.2238,  /* I */
    /* Other */
    0.0000    /* OTHER */
};

/* Intercept term from linear regression */
static const double WC_INTERCEPT = 1.3924;

/* Additional molecular feature coefficients */
static const double COEF_RING_COUNT = -1.0047;
static const double COEF_AROM_RING = -0.2689;
static const double COEF_H_DONOR = 1.2144;
static const double COEF_H_ACCEPTOR = 0.1247;
static const double COEF_ROT_BONDS = 0.0539;

/* ============================================================================
 * Helper Functions (shared patterns with counts.c)
 * ============================================================================ */

/* Count heteroatom neighbors (not C, H) */
static int count_heteroatom_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        element_t e = mol->atoms[atom->neighbors[i]].element;
        if (e != ELEM_C && e != ELEM_H) {
            count++;
        }
    }
    return count;
}

/* Check if atom has double bond to specific element */
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

/* Check if atom has triple bond */
static bool has_triple_bond(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_TRIPLE) {
            return true;
        }
    }
    return false;
}

/* Count double bonds */
static int count_double_bonds(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE) {
            count++;
        }
    }
    return count;
}

/* Get total hydrogens (explicit + implicit) */
static int get_total_h(const molecule_t* mol, const atom_t* atom) {
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) {
            h++;
        }
    }
    return h;
}

/* Count heavy (non-H) neighbors */
static int count_heavy_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            count++;
        }
    }
    return count;
}

/* Check if bonded to aromatic atom */
static bool bonded_to_aromatic(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].aromatic) {
            return true;
        }
    }
    return false;
}

/* ============================================================================
 * Atom Type Assignment
 * ============================================================================ */

/* Check if carbon is at ring fusion (shared between 2+ rings) */
static bool is_ring_fusion(const molecule_t* mol, const atom_t* atom) {
    if (atom->ring_count < 2) return false;
    int arom_neighbors = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].aromatic) {
            arom_neighbors++;
        }
    }
    return (arom_neighbors >= 3);
}

/* Check if aromatic C is substituted with non-C, non-H */
static bool is_aromatic_substituted(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        const atom_t* nb = &mol->atoms[atom->neighbors[i]];
        if (!nb->aromatic && nb->element != ELEM_H && nb->element != ELEM_C) {
            return true;
        }
        if (nb->element == ELEM_C && !nb->aromatic) {
            return true;
        }
    }
    return false;
}

static wc_atom_type_t assign_carbon_type(const molecule_t* mol, const atom_t* atom) {
    int h_count = get_total_h(mol, atom);
    int heavy_neighbors = count_heavy_neighbors(mol, atom);
    int hetero_neighbors = count_heteroatom_neighbors(mol, atom);
    int double_bonds = count_double_bonds(mol, atom);
    bool has_triple = has_triple_bond(mol, atom);

    /* Aromatic carbon */
    if (atom->aromatic) {
        if (is_ring_fusion(mol, atom)) return WC_C26;
        if (is_aromatic_substituted(mol, atom)) return WC_C25;
        return WC_C24;
    }

    /* Carbonyl carbon (C=O, C=S) */
    if (has_double_bond_to(mol, atom, ELEM_O) || has_double_bond_to(mol, atom, ELEM_S)) {
        return WC_C27;
    }

    /* Triple bond (sp carbon) */
    if (has_triple) {
        return (h_count > 0) ? WC_C22 : WC_C23;
    }

    /* Double bond to carbon (sp2) */
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

    /* sp3 carbon */
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

    if (atom->charge == 1 && has_double_bond_to(mol, atom, ELEM_O)) {
        return WC_N10;
    }

    if (atom->charge == 1 && heavy_neighbors == 4) {
        return WC_N9;
    }

    if (atom->aromatic) {
        if (h_count > 0) return WC_N5;
        return WC_N7;
    }

    for (int i = 0; i < atom->num_neighbors; i++) {
        const atom_t* nb = &mol->atoms[atom->neighbors[i]];
        if (nb->element == ELEM_C && has_double_bond_to(mol, nb, ELEM_O)) {
            return WC_N8;
        }
    }

    if (bonded_to_aromatic(mol, atom)) {
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
                    if (nb2 != atom && nb2->element == ELEM_O) {
                        int nb2_h = get_total_h(mol, nb2);
                        if (nb2_h > 0) return WC_O6;
                    }
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
            if (nb->element == ELEM_C && has_double_bond_to(mol, nb, ELEM_O)) {
                return WC_O8;
            }
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

/* Check if bond is rotatable */
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
 * LogP Calculation
 * ============================================================================ */

static double calculate_wclogp(const molecule_t* mol) {
    double logp = WC_INTERCEPT;

    int h_donors = 0;
    int h_acceptors = 0;
    int arom_atoms = 0;

    /* Sum contributions from all atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        wc_atom_type_t type = assign_atom_type(mol, i);
        logp += WC_LOGP_CONTRIB[type];

        const atom_t* atom = &mol->atoms[i];

        /* Count aromatic atoms for aromatic ring estimate */
        if (atom->aromatic) arom_atoms++;

        /* Count H-bond acceptors (N, O) */
        if (atom->element == ELEM_N || atom->element == ELEM_O) {
            h_acceptors++;
        }
    }

    /* Add implicit hydrogen contributions */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_H && atom->implicit_h_count > 0) {
            wc_atom_type_t h_type = WC_H1;
            if (atom->element != ELEM_C) {
                h_type = WC_H3;
                /* H on N or O is H-bond donor */
                if (atom->element == ELEM_N || atom->element == ELEM_O) {
                    h_donors += atom->implicit_h_count;
                }
            } else if (atom->aromatic) {
                h_type = WC_H4;
            } else if (count_heteroatom_neighbors(mol, atom) > 0) {
                h_type = WC_H2;
            }
            logp += atom->implicit_h_count * WC_LOGP_CONTRIB[h_type];
        }
    }

    /* Count rotatable bonds */
    int rot_bonds = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        if (is_rotatable_bond(mol, &mol->bonds[i])) {
            rot_bonds++;
        }
    }

    /* Add extra feature contributions */
    logp += mol->num_rings * COEF_RING_COUNT;
    logp += (arom_atoms / 5) * COEF_AROM_RING;  /* Approximate aromatic ring count */
    logp += h_donors * COEF_H_DONOR;
    logp += h_acceptors * COEF_H_ACCEPTOR;
    logp += rot_bonds * COEF_ROT_BONDS;

    return logp;
}

/* ============================================================================
 * Descriptor Implementation - Linear WCLogP
 * ============================================================================ */

static cchem_status_t desc_wclogp(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    value->d = calculate_wclogp(mol);
    return CCHEM_OK;
}

/* ============================================================================
 * Neural Network ccLogP
 * Architecture: Input(85) -> Hidden(128, ReLU) -> Hidden(64, ReLU) -> Output(1)
 * Trained on 12,766 molecules: R² = 0.9900, MAE = 0.1326, RMSE = 0.1828
 * ============================================================================ */

#define NN_INPUT_LOGP 85
#define NN_HIDDEN1_LOGP 128
#define NN_HIDDEN2_LOGP 64

#include "logp_weights.h"

/* Feature extraction for neural network */
static void extract_logp_features(const molecule_t* mol, double* features) {
    int wc_counts[WC_NUM_TYPES] = {0};
    int implicit_h[4] = {0};

    int h_donors = 0, h_acceptors = 0, arom_atoms = 0;

    /* Count WC atom types */
    for (int i = 0; i < mol->num_atoms; i++) {
        wc_atom_type_t type = assign_atom_type(mol, i);
        wc_counts[type]++;

        const atom_t* atom = &mol->atoms[i];
        if (atom->aromatic) arom_atoms++;
        if (atom->element == ELEM_N || atom->element == ELEM_O) h_acceptors++;
    }

    /* Count implicit H by type */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_H && atom->implicit_h_count > 0) {
            int h_type = 0;
            if (atom->element != ELEM_C) {
                h_type = 2;
                if (atom->element == ELEM_N || atom->element == ELEM_O)
                    h_donors += atom->implicit_h_count;
            } else if (atom->aromatic) {
                h_type = 3;
            } else if (count_heteroatom_neighbors(mol, atom) > 0) {
                h_type = 1;
            }
            implicit_h[h_type] += atom->implicit_h_count;
        }
    }

    /* Build feature vector */
    int f = 0;
    features[f++] = 1.0;  /* Intercept */

    /* 64 WC non-H types */
    for (int t = 0; t < WC_NUM_TYPES; t++) {
        if (t >= WC_H1 && t <= WC_H4) continue;
        features[f++] = (double)wc_counts[t];
    }

    /* 4 implicit H types */
    for (int h = 0; h < 4; h++) features[f++] = (double)implicit_h[h];

    /* 4 molecular features */
    features[f++] = (double)mol->num_rings;
    features[f++] = (double)arom_atoms / 6.0;
    features[f++] = (double)h_donors;
    features[f++] = (double)h_acceptors;

    /* 12 cchem descriptors */
    descriptor_value_t val;
    double tpsa = 0, chi0v = 0, chi1v = 0, zagreb1 = 0, clogs = 0;
    double molar_refr = 0, vdw_vol = 0, wiener = 0, balaban = 0;
    double polar_periph = 0, hydro_core = 0, lipo_chain = 0;

    if (descriptor_compute(mol, "TPSA", &val) == CCHEM_OK) tpsa = val.d;
    if (descriptor_compute(mol, "Chi0v", &val) == CCHEM_OK) chi0v = val.d;
    if (descriptor_compute(mol, "Chi1v", &val) == CCHEM_OK) chi1v = val.d;
    if (descriptor_compute(mol, "Zagreb1", &val) == CCHEM_OK) zagreb1 = val.d;
    if (descriptor_compute(mol, "CLogS", &val) == CCHEM_OK) clogs = val.d;
    if (descriptor_compute(mol, "MolarRefractivity", &val) == CCHEM_OK) molar_refr = val.d;
    if (descriptor_compute(mol, "VdWVolume", &val) == CCHEM_OK) vdw_vol = val.d;
    if (descriptor_compute(mol, "WienerIndex", &val) == CCHEM_OK) wiener = val.d;
    if (descriptor_compute(mol, "BalabanJ", &val) == CCHEM_OK) balaban = val.d;
    if (descriptor_compute(mol, "PolarPeripheryRatio", &val) == CCHEM_OK) polar_periph = val.d;
    if (descriptor_compute(mol, "HydrophobicCoreIndex", &val) == CCHEM_OK) hydro_core = val.d;
    if (descriptor_compute(mol, "LipophilicChainLength", &val) == CCHEM_OK) lipo_chain = val.d;

    features[f++] = tpsa / 100.0;
    features[f++] = chi0v;
    features[f++] = chi1v;
    features[f++] = zagreb1 / 10.0;
    features[f++] = clogs;
    features[f++] = molar_refr / 100.0;
    features[f++] = vdw_vol / 100.0;
    features[f++] = log1p(wiener) / 10.0;
    features[f++] = balaban;
    features[f++] = polar_periph;
    features[f++] = hydro_core;
    features[f++] = lipo_chain / 10.0;
}

/* Neural network forward pass */
static inline double relu_logp(double x) { return x > 0 ? x : 0; }

static double nn_predict_logp(const double* features) {
    double h1[NN_HIDDEN1_LOGP], h2[NN_HIDDEN2_LOGP];

    /* Normalize features */
    double norm[NN_INPUT_LOGP];
    for (int i = 0; i < NN_INPUT_LOGP; i++) {
        norm[i] = (features[i] - FEAT_MEAN[i]) / FEAT_STD[i];
        if (isnan(norm[i]) || isinf(norm[i])) norm[i] = 0.0;
    }

    /* Hidden layer 1 */
    for (int j = 0; j < NN_HIDDEN1_LOGP; j++) {
        double sum = NN_B1[j];
        for (int i = 0; i < NN_INPUT_LOGP; i++) {
            sum += norm[i] * NN_W1[i * NN_HIDDEN1_LOGP + j];
        }
        h1[j] = relu_logp(sum);
    }

    /* Hidden layer 2 */
    for (int j = 0; j < NN_HIDDEN2_LOGP; j++) {
        double sum = NN_B2[j];
        for (int i = 0; i < NN_HIDDEN1_LOGP; i++) {
            sum += h1[i] * NN_W2[i * NN_HIDDEN2_LOGP + j];
        }
        h2[j] = relu_logp(sum);
    }

    /* Output */
    double output = NN_B3;
    for (int i = 0; i < NN_HIDDEN2_LOGP; i++) {
        output += h2[i] * NN_W3[i];
    }

    return output;
}

static cchem_status_t desc_cclogp(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double features[NN_INPUT_LOGP];
    extract_logp_features(mol, features);
    value->d = nn_predict_logp(features);
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

void descriptors_register_logp(void) {
    /* Linear WCLogP */
    descriptor_def_t def = {0};
    strncpy(def.name, "WCLogP", MAX_DESCRIPTOR_NAME - 1);
    strncpy(def.description, "Wildman-Crippen LogP (linear model)",
            sizeof(def.description) - 1);
    def.category = DESC_CATEGORY_PROPERTIES;
    def.value_type = DESC_VALUE_DOUBLE;
    def.compute = desc_wclogp;
    descriptor_register(&def);

    /* Neural network ccLogP */
    descriptor_def_t def2 = {0};
    strncpy(def2.name, "ccLogP", MAX_DESCRIPTOR_NAME - 1);
    strncpy(def2.description, "Neural network LogP (R²=0.99, MAE=0.13)",
            sizeof(def2.description) - 1);
    def2.category = DESC_CATEGORY_PROPERTIES;
    def2.value_type = DESC_VALUE_DOUBLE;
    def2.compute = desc_cclogp;
    descriptor_register(&def2);
}
