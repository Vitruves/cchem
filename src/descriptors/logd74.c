/**
 * @file logd74.c
 * @brief ccLogD74 - Distribution coefficient at pH 7.4
 *
 * Neural network model for LogD prediction at physiological pH.
 * Architecture: Input(108) -> Hidden(128, ReLU) -> Hidden(64, ReLU) -> Output(1)
 * Trained on 26,036 molecules from ChEMBL, verified performance:
 *   R² = 0.9427, Pearson R = 0.9719, MAE = 0.2261, RMSE = 0.3194
 *
 * Features (108 total):
 * - 64 Wildman-Crippen atom type counts
 * - 4 implicit hydrogen type counts
 * - 7 ionizable group counts (linear, squared, log-transformed)
 * - 4 molecular topology features
 * - 12 cchem descriptors (TPSA, Chi, Zagreb, CLogS, etc.)
 * - 2 physics-based ionization corrections (Henderson-Hasselbalch)
 */

#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Neural Network Architecture Constants
 * ============================================================================ */

#define NN_INPUT 108
#define NN_HIDDEN1 128
#define NN_HIDDEN2 64

/* ============================================================================
 * Wildman-Crippen Atom Types (same as logp.c)
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
    WC_P, WC_F, WC_Cl, WC_Br, WC_I, WC_OTHER,
    WC_NUM_TYPES
} wc_atom_type_t;

/* WC LogP contributions for calculating base LogP */
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
    0.6808, -0.1986, 0.1285, 0.1537, -0.2238, 0.0000
};

static const double WC_INTERCEPT = 1.3924;
static const double COEF_RING_COUNT = -1.0047;
static const double COEF_AROM_RING = -0.2689;
static const double COEF_H_DONOR = 1.2144;
static const double COEF_H_ACCEPTOR = 0.1247;
static const double COEF_ROT_BONDS = 0.0539;

/* ============================================================================
 * Helper Functions for Atom Typing
 * ============================================================================ */

static int count_heteroatom_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        element_t e = mol->atoms[atom->neighbors[i]].element;
        if (e != ELEM_C && e != ELEM_H) count++;
    }
    return count;
}

static bool has_double_bond_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE &&
            mol->atoms[atom->neighbors[i]].element == elem) return true;
    }
    return false;
}

static bool has_triple_bond(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->bonds[atom->neighbor_bonds[i]].type == BOND_TRIPLE) return true;
    }
    return false;
}

static int count_double_bonds(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->bonds[atom->neighbor_bonds[i]].type == BOND_DOUBLE) count++;
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

static int count_heavy_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) count++;
    }
    return count;
}

static bool bonded_to_aromatic(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].aromatic) return true;
    }
    return false;
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

/* ============================================================================
 * Atom Type Assignment Functions
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
    if (has_double_bond_to(mol, atom, ELEM_O) || has_double_bond_to(mol, atom, ELEM_S)) return WC_C27;
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
            return (count_heteroatom_neighbors(mol, neighbor) > 0) ? WC_H2 : WC_H1;
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
                    if (nb2 != atom && nb2->element == ELEM_O && get_total_h(mol, nb2) > 0) return WC_O6;
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
    if (count_double_bonds(mol, atom) > 0) return WC_S3;
    if (get_total_h(mol, atom) > 0) return WC_S1;
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

static bool is_rotatable_bond(const molecule_t* mol, const bond_t* bond) {
    if (bond->type != BOND_SINGLE || bond->in_ring) return false;
    const atom_t* a1 = &mol->atoms[bond->atom1];
    const atom_t* a2 = &mol->atoms[bond->atom2];
    if (a1->element == ELEM_H || a2->element == ELEM_H) return false;
    return (count_heavy_neighbors(mol, a1) >= 2 && count_heavy_neighbors(mol, a2) >= 2);
}

/* ============================================================================
 * Ionizable Group Detection
 * ============================================================================ */

typedef struct {
    int n_cooh, n_phenol, n_amine_pri, n_amine_sec, n_amine_ter;
    int n_pyridine, n_imidazole;
} ionizable_t;

static bool is_guanidine_n(const molecule_t* mol, const atom_t* n_atom) {
    for (int i = 0; i < n_atom->num_neighbors; i++) {
        const atom_t* c = &mol->atoms[n_atom->neighbors[i]];
        if (c->element != ELEM_C) continue;
        bool has_dbl_n = false;
        int n_count = 0;
        for (int j = 0; j < c->num_neighbors; j++) {
            const atom_t* nb = &mol->atoms[c->neighbors[j]];
            if (nb->element == ELEM_N) {
                n_count++;
                if (mol->bonds[c->neighbor_bonds[j]].type == BOND_DOUBLE) has_dbl_n = true;
            }
        }
        if (has_dbl_n && n_count >= 3) return true;
    }
    return false;
}

static bool is_imidazole_n(const molecule_t* mol, const atom_t* n_atom) {
    if (!n_atom->aromatic || n_atom->ring_count == 0) return false;
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (ring->size != 5) continue;
        bool in_ring = false;
        int n_in_ring = 0;
        for (int i = 0; i < ring->size; i++) {
            if (&mol->atoms[ring->atoms[i]] == n_atom) in_ring = true;
            if (mol->atoms[ring->atoms[i]].element == ELEM_N) n_in_ring++;
        }
        if (in_ring && n_in_ring == 2) return true;
    }
    return false;
}

static void detect_ionizable(const molecule_t* mol, ionizable_t* g) {
    memset(g, 0, sizeof(*g));

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        if (atom->element == ELEM_N) {
            int h_count = get_total_h(mol, atom);
            int heavy_nb = count_heavy_neighbors(mol, atom);
            bool has_dbl = count_double_bonds(mol, atom) > 0;

            if (atom->charge == 1 && heavy_nb == 4) continue;
            if (is_guanidine_n(mol, atom)) continue;
            if (is_imidazole_n(mol, atom)) { g->n_imidazole++; continue; }
            if (atom->aromatic) { if (h_count == 0 && !has_dbl) g->n_pyridine++; continue; }

            bool is_amide = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* nb = &mol->atoms[atom->neighbors[j]];
                if (nb->element == ELEM_C && has_double_bond_to(mol, nb, ELEM_O)) { is_amide = true; break; }
            }
            if (is_amide) continue;
            if (bonded_to_aromatic(mol, atom)) continue;

            if (h_count == 2 && !has_dbl) g->n_amine_pri++;
            else if (h_count == 1 && !has_dbl) g->n_amine_sec++;
            else if (h_count == 0 && !has_dbl && heavy_nb == 3) g->n_amine_ter++;
        }

        if (atom->element == ELEM_O && get_total_h(mol, atom) > 0) {
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* c = &mol->atoms[atom->neighbors[j]];
                if (c->element == ELEM_C && has_double_bond_to(mol, c, ELEM_O)) { g->n_cooh++; goto next; }
                if (c->element == ELEM_C && c->aromatic) { g->n_phenol++; goto next; }
            }
        }
        next:;
    }
}

/* ============================================================================
 * WCLogP Calculation
 * ============================================================================ */

static double calculate_wclogp(const molecule_t* mol, int* h_donors_out, int* h_acceptors_out,
                               int* arom_atoms_out, int* rot_bonds_out) {
    double logp = WC_INTERCEPT;
    int h_donors = 0, h_acceptors = 0, arom_atoms = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        wc_atom_type_t type = assign_atom_type(mol, i);
        logp += WC_LOGP_CONTRIB[type];
        const atom_t* atom = &mol->atoms[i];
        if (atom->aromatic) arom_atoms++;
        if (atom->element == ELEM_N || atom->element == ELEM_O) h_acceptors++;
    }

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_H && atom->implicit_h_count > 0) {
            wc_atom_type_t h_type = WC_H1;
            if (atom->element != ELEM_C) {
                h_type = WC_H3;
                if (atom->element == ELEM_N || atom->element == ELEM_O)
                    h_donors += atom->implicit_h_count;
            } else if (atom->aromatic) {
                h_type = WC_H4;
            } else if (count_heteroatom_neighbors(mol, atom) > 0) {
                h_type = WC_H2;
            }
            logp += atom->implicit_h_count * WC_LOGP_CONTRIB[h_type];
        }
    }

    int rot_bonds = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        if (is_rotatable_bond(mol, &mol->bonds[i])) rot_bonds++;
    }

    logp += mol->num_rings * COEF_RING_COUNT;
    logp += (arom_atoms / 5) * COEF_AROM_RING;
    logp += h_donors * COEF_H_DONOR;
    logp += h_acceptors * COEF_H_ACCEPTOR;
    logp += rot_bonds * COEF_ROT_BONDS;

    *h_donors_out = h_donors;
    *h_acceptors_out = h_acceptors;
    *arom_atoms_out = arom_atoms;
    *rot_bonds_out = rot_bonds;

    return logp;
}

/* ============================================================================
 * Pre-trained Neural Network Weights
 * ============================================================================ */

#include "logd74_weights.h"

/* ============================================================================
 * Feature Extraction
 * ============================================================================ */

static void extract_features(const molecule_t* mol, double* features) {
    int wc_counts[WC_NUM_TYPES] = {0};
    int implicit_h[4] = {0};

    int h_donors = 0, h_acceptors = 0, arom_atoms = 0, rot_bonds = 0;
    (void)calculate_wclogp(mol, &h_donors, &h_acceptors, &arom_atoms, &rot_bonds);

    /* Count WC atom types */
    for (int i = 0; i < mol->num_atoms; i++) {
        wc_atom_type_t type = assign_atom_type(mol, i);
        wc_counts[type]++;
    }

    /* Count implicit H by type */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_H && atom->implicit_h_count > 0) {
            int h_type = 0;
            if (atom->element != ELEM_C) h_type = 2;
            else if (atom->aromatic) h_type = 3;
            else if (count_heteroatom_neighbors(mol, atom) > 0) h_type = 1;
            implicit_h[h_type] += atom->implicit_h_count;
        }
    }

    ionizable_t ion;
    detect_ionizable(mol, &ion);

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

    /* 7 ionization features */
    double ionf[7] = {
        (double)ion.n_cooh, (double)ion.n_amine_pri, (double)ion.n_amine_sec,
        (double)ion.n_amine_ter, (double)ion.n_phenol, (double)ion.n_pyridine,
        (double)ion.n_imidazole
    };
    for (int i = 0; i < 7; i++) features[f++] = ionf[i];

    /* 4 molecular features */
    features[f++] = (double)mol->num_rings;
    features[f++] = (double)arom_atoms / 6.0;
    features[f++] = (double)h_donors;
    features[f++] = (double)h_acceptors;

    /* 7 ionization^2 */
    for (int i = 0; i < 7; i++) features[f++] = ionf[i] * ionf[i];

    /* 7 log-transformed ionization */
    for (int i = 0; i < 7; i++) features[f++] = log1p(ionf[i]);

    /* 12 cchem descriptors - call real descriptors */
    descriptor_value_t val;
    double tpsa = 0, chi0v = 0, chi1v = 0, zagreb1 = 0, clogs = 0;
    double acidic_count = 0, basic_count = 0;
    double mean_acidic_pka = 0, mean_basic_pka = 0;
    double polar_periphery = 0, hydrophobic_core = 0, lipophilic_chain = 0;

    if (descriptor_compute(mol, "TPSA", &val) == CCHEM_OK) tpsa = val.d;
    if (descriptor_compute(mol, "Chi0v", &val) == CCHEM_OK) chi0v = val.d;
    if (descriptor_compute(mol, "Chi1v", &val) == CCHEM_OK) chi1v = val.d;
    if (descriptor_compute(mol, "Zagreb1", &val) == CCHEM_OK) zagreb1 = val.d;
    if (descriptor_compute(mol, "CLogS", &val) == CCHEM_OK) clogs = val.d;
    if (descriptor_compute(mol, "AcidicGroupCount", &val) == CCHEM_OK) acidic_count = (double)val.i;
    if (descriptor_compute(mol, "BasicGroupCount", &val) == CCHEM_OK) basic_count = (double)val.i;
    if (descriptor_compute(mol, "MeanAcidicPka", &val) == CCHEM_OK) mean_acidic_pka = val.d;
    if (descriptor_compute(mol, "MeanBasicPka", &val) == CCHEM_OK) mean_basic_pka = val.d;
    if (descriptor_compute(mol, "PolarPeripheryRatio", &val) == CCHEM_OK) polar_periphery = val.d;
    if (descriptor_compute(mol, "HydrophobicCoreIndex", &val) == CCHEM_OK) hydrophobic_core = val.d;
    if (descriptor_compute(mol, "LipophilicChainLength", &val) == CCHEM_OK) lipophilic_chain = val.d;

    features[f++] = tpsa / 100.0;
    features[f++] = chi0v;
    features[f++] = chi1v;
    features[f++] = zagreb1 / 10.0;
    features[f++] = clogs;
    features[f++] = acidic_count;
    features[f++] = basic_count;
    features[f++] = (mean_acidic_pka - 7.4) / 5.0;
    features[f++] = (mean_basic_pka - 7.4) / 5.0;
    features[f++] = polar_periphery;
    features[f++] = hydrophobic_core;
    features[f++] = lipophilic_chain / 10.0;

    /* 2 physics-based ionization corrections (using estimated pKa) */
    double acidic_corr = 0, basic_corr = 0;
    double pH = 7.4;

    /* Acidic groups: LogD shift = log10(fraction neutral) */
    if (acidic_count > 0 && mean_acidic_pka > 0 && mean_acidic_pka < 14) {
        double frac_ionized = 1.0 / (1.0 + pow(10.0, mean_acidic_pka - pH));
        acidic_corr = log10(1.0 - frac_ionized + 1e-6);
    }

    /* Basic groups: LogD shift = log10(fraction neutral) */
    if (basic_count > 0 && mean_basic_pka > 0 && mean_basic_pka < 14) {
        double frac_ionized = 1.0 / (1.0 + pow(10.0, pH - mean_basic_pka));
        basic_corr = log10(1.0 - frac_ionized + 1e-6);
    }

    features[f++] = acidic_corr;
    features[f++] = basic_corr;
}

/* ============================================================================
 * Neural Network Forward Pass
 * ============================================================================ */

static inline double relu(double x) { return x > 0 ? x : 0; }

static double nn_predict(const double* features) {
    double h1[NN_HIDDEN1] = {0}, h2[NN_HIDDEN2] = {0};

    /* Normalize features */
    double norm[NN_INPUT] = {0};
    for (int i = 0; i < NN_INPUT; i++) {
        norm[i] = (features[i] - FEAT_MEAN[i]) / FEAT_STD[i];
        if (isnan(norm[i]) || isinf(norm[i])) norm[i] = 0.0;
    }

    /* Hidden layer 1 */
    for (int j = 0; j < NN_HIDDEN1; j++) {
        double sum = NN_B1[j];
        for (int i = 0; i < NN_INPUT; i++) {
            sum += norm[i] * NN_W1[i * NN_HIDDEN1 + j];
        }
        h1[j] = relu(sum);
    }

    /* Hidden layer 2 */
    for (int j = 0; j < NN_HIDDEN2; j++) {
        double sum = NN_B2[j];
        for (int i = 0; i < NN_HIDDEN1; i++) {
            sum += h1[i] * NN_W2[i * NN_HIDDEN2 + j];
        }
        h2[j] = relu(sum);
    }

    /* Output */
    double output = NN_B3;
    for (int i = 0; i < NN_HIDDEN2; i++) {
        output += h2[i] * NN_W3[i];
    }

    return output;
}

/* ============================================================================
 * Descriptor Compute Function
 * ============================================================================ */

static cchem_status_t compute_cclogd74(const molecule_t* mol, descriptor_value_t* value) {
    double features[NN_INPUT] = {0};
    extract_features(mol, features);
    value->d = nn_predict(features);
    return CCHEM_OK;
}

/* ============================================================================
 * Descriptor Registration
 * ============================================================================ */

void descriptors_register_logd74(void) {
    descriptor_def_t def = {
        .name = "ccLogD74",
        .description = "Distribution coefficient at pH 7.4 (neural network model)",
        .category = DESC_CATEGORY_PROPERTIES,
        .value_type = DESC_VALUE_DOUBLE,
        .compute = compute_cclogd74,
        .user_data = NULL,
        .registered = false
    };
    descriptor_register(&def);
}
