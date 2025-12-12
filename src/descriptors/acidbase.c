/**
 * @file acidbase.c
 * @brief Acid-base functional group descriptors with pKa-based classification
 *
 * This module provides descriptors for counting acidic and basic functional groups
 * classified by their approximate pKa ranges:
 *
 * ACIDIC GROUPS:
 * - SuperAcidicCount:       pKa < -1   (middle: -3)   - Sulfonic acids, superacids
 * - StronglyAcidicCount:    pKa -1 to 1 (middle: 0)   - Phosphoric, polyhaloacetic acids
 * - ModeratelyAcidicCount:  pKa 1 to 5  (middle: 3)   - Carboxylic acids, tetrazoles
 * - WeaklyAcidicCount:      pKa 5 to 14 (middle: 9.5) - Phenols, thiols, sulfonamides
 *
 * BASIC GROUPS:
 * - SuperBasicCount:        conj. pKa > 15 (middle: 17)  - Alkoxides, amide anions
 * - StronglyBasicCount:     conj. pKa 10-14 (middle: 12) - Aliphatic amines, guanidines
 * - ModeratelyBasicCount:   conj. pKa 5-10  (middle: 7.5) - Pyridines, imidazoles
 * - WeaklyBasicCount:       conj. pKa 0-5   (middle: 2.5) - Anilines, amides, ethers
 *
 * Aggregate descriptors:
 * - TotalAcidicFunctions:   Sum of all acidic groups
 * - TotalBasicFunctions:    Sum of all basic groups
 * - MeanAcidicPotential:    Weighted mean pKa of acidic groups
 * - MeanBasicPotential:     Weighted mean pKa of basic groups
 */

#include <string.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* pKa midpoints for each category */
#define PKA_SUPER_ACIDIC_MID    -3.0
#define PKA_STRONG_ACIDIC_MID    0.0
#define PKA_MODERATE_ACIDIC_MID  3.0
#define PKA_WEAK_ACIDIC_MID      9.5

#define PKA_SUPER_BASIC_MID     17.0
#define PKA_STRONG_BASIC_MID    12.0
#define PKA_MODERATE_BASIC_MID   7.5
#define PKA_WEAK_BASIC_MID       2.5

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

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

/* Get total H count (implicit + explicit) */
static int get_total_h(const molecule_t* mol, const atom_t* atom) {
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) {
            h++;
        }
    }
    return h;
}

/* Check if atom has double bond to specific element */
static bool has_double_bond_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE) {
            if (mol->atoms[atom->neighbors[i]].element == elem) {
                return true;
            }
        }
    }
    return false;
}

/* Check if atom has triple bond to specific element */
static bool has_triple_bond_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_TRIPLE) {
            if (mol->atoms[atom->neighbors[i]].element == elem) {
                return true;
            }
        }
    }
    return false;
}

/* Count double bonds to specific element */
static int count_double_bonds_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE) {
            if (mol->atoms[atom->neighbors[i]].element == elem) {
                count++;
            }
        }
    }
    return count;
}

/* Get neighbor of specific element */
static const atom_t* get_neighbor_of_element(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == elem) {
            return &mol->atoms[atom->neighbors[i]];
        }
    }
    return NULL;
}

/* Count neighbors of specific element */
static int count_neighbors_of_element(const molecule_t* mol, const atom_t* atom, element_t elem) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == elem) {
            count++;
        }
    }
    return count;
}

/* Count halogens bonded to atom */
static int count_bonded_halogens(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        element_t e = mol->atoms[atom->neighbors[i]].element;
        if (e == ELEM_F || e == ELEM_Cl || e == ELEM_Br || e == ELEM_I) {
            count++;
        }
    }
    return count;
}

/* Count fluorines bonded to atom (strongest EWG) */
static int count_bonded_fluorines(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_F) {
            count++;
        }
    }
    return count;
}

/* Check if atom is in aromatic ring of given size */
static bool is_in_aromatic_ring_of_size(const molecule_t* mol, int atom_idx, int ring_size) {
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (ring->aromatic && ring->size == ring_size) {
            for (int j = 0; j < ring->size; j++) {
                if (ring->atoms[j] == atom_idx) {
                    return true;
                }
            }
        }
    }
    return false;
}

/* Count nitrogens in the same aromatic ring as given atom */
static int count_ring_nitrogens(const molecule_t* mol, int atom_idx) {
    int max_n = 0;
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (!ring->aromatic) continue;

        bool atom_in_ring = false;
        int n_count = 0;

        for (int j = 0; j < ring->size; j++) {
            if (ring->atoms[j] == atom_idx) atom_in_ring = true;
            if (mol->atoms[ring->atoms[j]].element == ELEM_N) n_count++;
        }

        if (atom_in_ring && n_count > max_n) {
            max_n = n_count;
        }
    }
    return max_n;
}

/* Check if atom is in a ring containing NH nitrogen */
static bool is_in_ring_with_nh(const molecule_t* mol, int atom_idx) {
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (!ring->aromatic) continue;

        bool atom_in_ring = false;
        bool has_nh = false;

        for (int j = 0; j < ring->size; j++) {
            if (ring->atoms[j] == atom_idx) atom_in_ring = true;
            const atom_t* ra = &mol->atoms[ring->atoms[j]];
            if (ra->element == ELEM_N && get_total_h(mol, ra) >= 1) {
                has_nh = true;
            }
        }

        if (atom_in_ring && has_nh) return true;
    }
    return false;
}

/* Check if phenol (hydroxyl on aromatic ring) */
static bool is_phenol(const molecule_t* mol, const atom_t* o_atom) {
    if (o_atom->element != ELEM_O) return false;
    if (o_atom->aromatic) return false;

    int h = get_total_h(mol, o_atom);
    if (h < 1) return false;

    for (int i = 0; i < o_atom->num_neighbors; i++) {
        const atom_t* neighbor = &mol->atoms[o_atom->neighbors[i]];
        if (neighbor->element == ELEM_C && neighbor->aromatic) {
            return true;
        }
    }
    return false;
}

/* Check if thiophenol (thiol on aromatic ring) */
static bool is_thiophenol(const molecule_t* mol, const atom_t* s_atom) {
    if (s_atom->element != ELEM_S) return false;
    if (s_atom->aromatic) return false;

    int h = get_total_h(mol, s_atom);
    if (h < 1) return false;

    for (int i = 0; i < s_atom->num_neighbors; i++) {
        const atom_t* neighbor = &mol->atoms[s_atom->neighbors[i]];
        if (neighbor->element == ELEM_C && neighbor->aromatic) {
            return true;
        }
    }
    return false;
}

/* Count electron-withdrawing groups on aromatic ring */
static int count_ewg_on_ring(const molecule_t* mol, const atom_t* ring_atom) {
    int ewg_count = 0;

    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (!ring->aromatic) continue;

        bool atom_in_ring = false;
        for (int j = 0; j < ring->size; j++) {
            if (&mol->atoms[ring->atoms[j]] == ring_atom) {
                atom_in_ring = true;
                break;
            }
        }
        if (!atom_in_ring) continue;

        for (int j = 0; j < ring->size; j++) {
            const atom_t* ra = &mol->atoms[ring->atoms[j]];
            for (int k = 0; k < ra->num_neighbors; k++) {
                const atom_t* sub = &mol->atoms[ra->neighbors[k]];

                /* Skip ring atoms */
                bool sub_in_ring = false;
                for (int m = 0; m < ring->size; m++) {
                    if (ring->atoms[m] == ra->neighbors[k]) {
                        sub_in_ring = true;
                        break;
                    }
                }
                if (sub_in_ring) continue;

                /* Nitro group: N+ with O */
                if (sub->element == ELEM_N && sub->charge > 0) {
                    ewg_count += 2;  /* Strong EWG */
                }
                /* Carbonyl */
                else if (sub->element == ELEM_C && has_double_bond_to(mol, sub, ELEM_O)) {
                    ewg_count++;
                }
                /* Cyano */
                else if (sub->element == ELEM_C && has_triple_bond_to(mol, sub, ELEM_N)) {
                    ewg_count += 2;  /* Strong EWG */
                }
                /* Sulfonyl */
                else if (sub->element == ELEM_S && count_double_bonds_to(mol, sub, ELEM_O) >= 2) {
                    ewg_count += 2;
                }
                /* Halogens */
                else if (sub->element == ELEM_F) {
                    ewg_count++;  /* F is strong */
                }
                else if (sub->element == ELEM_Cl || sub->element == ELEM_Br || sub->element == ELEM_I) {
                    /* Halogens are weak EWG by induction */
                }
                /* Trifluoromethyl */
                else if (sub->element == ELEM_C) {
                    int f_count = count_bonded_fluorines(mol, sub);
                    if (f_count >= 3) ewg_count += 2;
                }
            }
        }
    }
    return ewg_count;
}

/* Count electron-donating groups on aromatic ring */
static int count_edg_on_ring(const molecule_t* mol, const atom_t* ring_atom) {
    int edg_count = 0;

    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (!ring->aromatic) continue;

        bool atom_in_ring = false;
        for (int j = 0; j < ring->size; j++) {
            if (&mol->atoms[ring->atoms[j]] == ring_atom) {
                atom_in_ring = true;
                break;
            }
        }
        if (!atom_in_ring) continue;

        for (int j = 0; j < ring->size; j++) {
            const atom_t* ra = &mol->atoms[ring->atoms[j]];
            for (int k = 0; k < ra->num_neighbors; k++) {
                const atom_t* sub = &mol->atoms[ra->neighbors[k]];

                /* Skip ring atoms */
                bool sub_in_ring = false;
                for (int m = 0; m < ring->size; m++) {
                    if (ring->atoms[m] == ra->neighbors[k]) {
                        sub_in_ring = true;
                        break;
                    }
                }
                if (sub_in_ring) continue;

                /* Amino group: NH2, NHR, NR2 */
                if (sub->element == ELEM_N && !sub->aromatic && sub->charge == 0) {
                    if (!has_double_bond_to(mol, sub, ELEM_C) &&
                        !has_double_bond_to(mol, sub, ELEM_O)) {
                        edg_count++;
                    }
                }
                /* Alkoxy/hydroxy */
                else if (sub->element == ELEM_O && !sub->aromatic) {
                    if (!has_double_bond_to(mol, sub, ELEM_C)) {
                        edg_count++;
                    }
                }
                /* Alkyl groups */
                else if (sub->element == ELEM_C && !sub->aromatic) {
                    /* Check it's not a carbonyl carbon */
                    if (!has_double_bond_to(mol, sub, ELEM_O) &&
                        !has_triple_bond_to(mol, sub, ELEM_N)) {
                        /* Simple alkyl is weak EDG */
                    }
                }
            }
        }
    }
    return edg_count;
}

/* Check if nitrogen is part of a guanidine N=C(N)N */
static bool is_guanidine_nitrogen(const molecule_t* mol, const atom_t* n_atom) {
    for (int i = 0; i < n_atom->num_neighbors; i++) {
        const atom_t* neighbor = &mol->atoms[n_atom->neighbors[i]];
        if (neighbor->element == ELEM_C) {
            /* Check if C has =N and another N */
            bool has_double_n = has_double_bond_to(mol, neighbor, ELEM_N);
            int n_count = count_neighbors_of_element(mol, neighbor, ELEM_N);
            if (has_double_n && n_count >= 3) {
                return true;
            }
        }
    }
    return false;
}

/* Check if nitrogen is part of an amidine R-C(=N)-N */
static bool is_amidine_nitrogen(const molecule_t* mol, const atom_t* n_atom) {
    for (int i = 0; i < n_atom->num_neighbors; i++) {
        const atom_t* neighbor = &mol->atoms[n_atom->neighbors[i]];
        if (neighbor->element == ELEM_C) {
            bool has_double_n = has_double_bond_to(mol, neighbor, ELEM_N);
            int n_count = count_neighbors_of_element(mol, neighbor, ELEM_N);
            if (has_double_n && n_count >= 2 && n_count < 3) {
                return true;
            }
        }
    }
    return false;
}

/* ============================================================================
 * Acidic Group Detection Functions
 * ============================================================================ */

/**
 * Super Acidic Groups (pKa < -1):
 * - Sulfonic acids: RS(=O)(=O)OH (~-2 to -3)
 * - Perfluorosulfonic acids (~-6 to -14)
 * - Fluorosulfonic acid, chlorosulfonic acid
 */
static int count_super_acidic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Sulfonic acid: S with 2 =O and one -OH */
        if (atom->element == ELEM_S) {
            int double_o = count_double_bonds_to(mol, atom, ELEM_O);

            if (double_o >= 2) {
                /* Check for -OH on sulfur */
                for (int j = 0; j < atom->num_neighbors; j++) {
                    const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                    if (neighbor->element == ELEM_O) {
                        int bond_idx = atom->neighbor_bonds[j];
                        bond_type_t bt = mol->bonds[bond_idx].type;
                        if (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) {
                            int o_h = get_total_h(mol, neighbor);
                            if (o_h > 0) {
                                count++;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    return count;
}

/**
 * Strongly Acidic Groups (pKa -1 to 1):
 * - Phosphoric acids: P(=O)(OH)x (~2.1 first pKa, but di/tri-phosphates stronger)
 * - Phosphonic acids: RP(=O)(OH)2 (~1.5-2.5)
 * - Trifluoroacetic acid: CF3-COOH (~0.5)
 * - Polyhaloacetic acids: CX3-COOH, CX2H-COOH (~0.7-1.3)
 * - Sulfinic acids: RS(=O)OH (~2)
 * - Oxalic acid (first pKa ~1.25)
 * - Picric acid (2,4,6-trinitrophenol ~0.4)
 */
static int count_strongly_acidic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Phosphoric/Phosphonic acids */
        if (atom->element == ELEM_P) {
            bool has_double_o = has_double_bond_to(mol, atom, ELEM_O);

            int oh_count = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                if (neighbor->element == ELEM_O) {
                    int bond_idx = atom->neighbor_bonds[j];
                    bond_type_t bt = mol->bonds[bond_idx].type;
                    if (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) {
                        int o_h = get_total_h(mol, neighbor);
                        if (o_h > 0) {
                            oh_count++;
                        }
                    }
                }
            }

            if (has_double_o && oh_count > 0) {
                count += oh_count;
            }
        }

        /* Sulfinic acid: RS(=O)OH (one =O, one -OH) */
        if (atom->element == ELEM_S) {
            int double_o = count_double_bonds_to(mol, atom, ELEM_O);

            if (double_o == 1) {
                for (int j = 0; j < atom->num_neighbors; j++) {
                    const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                    if (neighbor->element == ELEM_O) {
                        int bond_idx = atom->neighbor_bonds[j];
                        bond_type_t bt = mol->bonds[bond_idx].type;
                        if (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) {
                            int o_h = get_total_h(mol, neighbor);
                            if (o_h > 0) {
                                count++;
                                break;
                            }
                        }
                    }
                }
            }
        }

        /* Carboxylic acid alpha to multiple halogens or strong EWG */
        if (atom->element == ELEM_C && has_double_bond_to(mol, atom, ELEM_O)) {
            bool has_oh = false;
            const atom_t* alpha_c = NULL;

            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                int bond_idx = atom->neighbor_bonds[j];
                bond_type_t bt = mol->bonds[bond_idx].type;

                if (neighbor->element == ELEM_O &&
                    (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN)) {
                    int o_h = get_total_h(mol, neighbor);
                    int o_heavy = count_heavy_neighbors(mol, neighbor);
                    if (o_h > 0 && o_heavy == 1) {
                        has_oh = true;
                    }
                }
                if (neighbor->element == ELEM_C &&
                    (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN)) {
                    alpha_c = neighbor;
                }
            }

            if (has_oh && alpha_c) {
                int halogens = count_bonded_halogens(mol, alpha_c);
                int fluorines = count_bonded_fluorines(mol, alpha_c);

                /* Trifluoro or 2+ halogens -> strongly acidic */
                if (fluorines >= 3 || halogens >= 2) {
                    count++;
                }
            }

            /* Alpha-oxo carboxylic (glyoxylic/pyruvic type) */
            if (has_oh && alpha_c && has_double_bond_to(mol, alpha_c, ELEM_O)) {
                count++;  /* Alpha-keto acid is more acidic */
            }
        }

        /* Highly activated phenols (trinitrophenol/picric acid type) */
        if (atom->element == ELEM_O && is_phenol(mol, atom)) {
            const atom_t* ar_c = NULL;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                if (neighbor->element == ELEM_C && neighbor->aromatic) {
                    ar_c = neighbor;
                    break;
                }
            }

            if (ar_c) {
                int ewg = count_ewg_on_ring(mol, ar_c);
                /* 3+ strong EWG (like trinitro) makes phenol super acidic */
                if (ewg >= 6) {  /* 3 nitro groups = 6 EWG units */
                    count++;
                }
            }
        }
    }

    return count;
}

/**
 * Moderately Acidic Groups (pKa 1 to 5):
 * - Carboxylic acids: R-COOH (~4-5)
 * - Haloacetic acids (mono): XCH2-COOH (~2.9)
 * - Dicarboxylic acids: malonic (~2.85), succinic (~4.2)
 * - Tetrazoles (~4.5-4.9, bioisostere of COOH)
 * - Strongly activated phenols (2+ EWG, pKa 4-7)
 * - Squaric acid derivatives (~1-3)
 * - Hydroxamic acids in some cases (~8-9, but placed here when activated)
 */
static int count_moderately_acidic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Carboxylic acid: C with =O and -OH */
        if (atom->element == ELEM_C && !atom->aromatic) {
            bool has_carbonyl_o = has_double_bond_to(mol, atom, ELEM_O);
            if (!has_carbonyl_o) continue;

            bool has_oh = false;
            bool is_ester = false;
            const atom_t* alpha_c = NULL;
            bool alpha_is_carbonyl = false;

            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                int bond_idx = atom->neighbor_bonds[j];
                bond_type_t bt = mol->bonds[bond_idx].type;

                if (neighbor->element == ELEM_O && mol->bonds[bond_idx].type != BOND_DOUBLE) {
                    int o_h = get_total_h(mol, neighbor);
                    int o_heavy = count_heavy_neighbors(mol, neighbor);

                    if (o_h > 0 && o_heavy == 1) {
                        has_oh = true;
                    } else if (o_h == 0 && o_heavy == 2) {
                        is_ester = true;
                    }
                }
                if (neighbor->element == ELEM_C &&
                    (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN)) {
                    alpha_c = neighbor;
                    if (has_double_bond_to(mol, neighbor, ELEM_O)) {
                        alpha_is_carbonyl = true;
                    }
                }
            }

            if (has_oh && !is_ester) {
                int halogens = 0;
                int fluorines = 0;
                if (alpha_c) {
                    halogens = count_bonded_halogens(mol, alpha_c);
                    fluorines = count_bonded_fluorines(mol, alpha_c);
                }

                /* Skip strongly acidic (handled above) */
                if (fluorines >= 3 || halogens >= 2) continue;
                if (alpha_is_carbonyl) continue;  /* Alpha-keto handled above */

                count++;
            }
        }

        /* Tetrazole: 5-membered aromatic ring with 4 N */
        if (atom->element == ELEM_N && atom->aromatic) {
            if (is_in_aromatic_ring_of_size(mol, i, 5)) {
                int ring_n = count_ring_nitrogens(mol, i);
                if (ring_n >= 4) {
                    /* Only count once per tetrazole (count NH) */
                    int h = get_total_h(mol, atom);
                    if (h >= 1) {
                        count++;
                    }
                }
            }
        }

        /* Activated phenols (2+ EWG but not super strong) */
        if (atom->element == ELEM_O && is_phenol(mol, atom)) {
            const atom_t* ar_c = NULL;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                if (neighbor->element == ELEM_C && neighbor->aromatic) {
                    ar_c = neighbor;
                    break;
                }
            }

            if (ar_c) {
                int ewg = count_ewg_on_ring(mol, ar_c);
                /* 2-5 EWG units -> moderately acidic phenol */
                if (ewg >= 2 && ewg < 6) {
                    count++;
                }
            }
        }

        /* 1,3-Dicarbonyl systems (barbituric acid type) */
        if (atom->element == ELEM_N && !atom->aromatic) {
            int h = get_total_h(mol, atom);
            if (h < 1) continue;

            int carbonyl_neighbors = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                    carbonyl_neighbors++;
                }
            }

            /* Uracil/barbiturate type: NH between two C=O */
            if (carbonyl_neighbors >= 2) {
                /* Check if in ring (cyclic imide) */
                if (atom->ring_count > 0) {
                    count++;  /* Cyclic imides are more acidic (~4-5) */
                }
            }
        }
    }

    return count;
}

/**
 * Weakly Acidic Groups (pKa 5 to 14):
 * - Simple phenols (~10)
 * - Mildly activated phenols (~7-9)
 * - Thiols: R-SH (~10-11)
 * - Thiophenols: Ar-SH (~6.5)
 * - Sulfonamides: R-SO2-NH2 (~10)
 * - Simple imides: -CO-NH-CO- (~8-10 for linear)
 * - Beta-diketones: -CO-CH2-CO- (~9)
 * - Beta-ketoesters, malonates (~10-13)
 * - Indoles, pyrroles (NH) (~15-17, weak but counted)
 * - Hydroxamic acids: R-CO-NH-OH (~8-9)
 * - Oximes: R2C=N-OH (~10-12)
 * - Triazoles (1,2,3 and 1,2,4) NH (~9-10)
 * - Benzimidazoles NH (~12-13)
 * - Cyanoacetates, malononitriles activated CH (~9-11)
 */
static int count_weakly_acidic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Simple phenols (0-1 EWG) - EDG make phenols less acidic */
        if (atom->element == ELEM_O && !atom->aromatic) {
            if (is_phenol(mol, atom)) {
                const atom_t* ar_c = NULL;
                for (int j = 0; j < atom->num_neighbors; j++) {
                    const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                    if (neighbor->element == ELEM_C && neighbor->aromatic) {
                        ar_c = neighbor;
                        break;
                    }
                }

                if (ar_c) {
                    int ewg = count_ewg_on_ring(mol, ar_c);
                    int edg = count_edg_on_ring(mol, ar_c);

                    /* Skip activated phenols (handled in moderate/strong) */
                    if (ewg >= 2) continue;

                    /* Phenols with strong EDG are less acidic (pKa ~10-11) */
                    /* Still count them as weakly acidic but with higher pKa */
                    if (edg >= 2) {
                        count++;  /* EDG-substituted phenol (pKa ~10-11) */
                    } else {
                        count++;  /* Simple phenol (pKa ~10) */
                    }
                }
            }
        }

        /* Thiols: R-SH (aliphatic ~10-11) and thiophenols (Ar-SH ~6.5) */
        if (atom->element == ELEM_S && !atom->aromatic) {
            int h = get_total_h(mol, atom);
            int heavy = count_heavy_neighbors(mol, atom);

            if (h >= 1 && heavy == 1) {
                if (is_thiophenol(mol, atom)) {
                    count++;  /* Thiophenol (pKa ~6.5) */
                } else {
                    count++;  /* Aliphatic thiol (pKa ~10-11) */
                }
            }
        }

        /* Sulfonamide: SO2-NH */
        if (atom->element == ELEM_S) {
            int double_o = count_double_bonds_to(mol, atom, ELEM_O);

            if (double_o >= 2) {
                for (int j = 0; j < atom->num_neighbors; j++) {
                    const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                    if (neighbor->element == ELEM_N) {
                        int n_h = get_total_h(mol, neighbor);
                        if (n_h >= 1) {
                            count++;
                            break;
                        }
                    }
                }
            }
        }

        /* Linear imide: -CO-NH-CO- (not cyclic, those are in moderate) */
        if (atom->element == ELEM_N && !atom->aromatic && atom->ring_count == 0) {
            int h = get_total_h(mol, atom);
            if (h < 1) continue;

            int carbonyl_count = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                    carbonyl_count++;
                }
            }

            if (carbonyl_count >= 2) {
                count++;
            }
        }

        /* Pyrrole/indole NH (aromatic N with H) */
        if (atom->element == ELEM_N && atom->aromatic) {
            int h = get_total_h(mol, atom);
            if (h >= 1) {
                /* Check it's not a tetrazole (handled in moderate) */
                int ring_n = count_ring_nitrogens(mol, i);
                if (ring_n < 4) {
                    count++;
                }
            }
        }

        /* Beta-diketone: C(=O)-CH-C(=O) or activated methylene */
        if (atom->element == ELEM_C && !atom->aromatic) {
            int h = get_total_h(mol, atom);
            if (h < 1) continue;

            int ewg_count = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];

                /* Carbonyl neighbor */
                if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                    ewg_count++;
                }
                /* Cyano neighbor */
                if (neighbor->element == ELEM_C && has_triple_bond_to(mol, neighbor, ELEM_N)) {
                    ewg_count++;
                }
                /* Ester carbonyl (malonate type) */
                if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                    /* Check if also has OR (ester) */
                    for (int k = 0; k < neighbor->num_neighbors; k++) {
                        const atom_t* nn = &mol->atoms[neighbor->neighbors[k]];
                        if (nn->element == ELEM_O) {
                            int bond_idx = neighbor->neighbor_bonds[k];
                            if (mol->bonds[bond_idx].type != BOND_DOUBLE) {
                                int o_h = get_total_h(mol, nn);
                                if (o_h == 0) {
                                    /* It's an ester, already counted carbonyl */
                                    break;
                                }
                            }
                        }
                    }
                }
                /* Nitro neighbor */
                if (neighbor->element == ELEM_N && neighbor->charge > 0) {
                    ewg_count++;
                }
                /* Sulfonyl neighbor */
                if (neighbor->element == ELEM_S && count_double_bonds_to(mol, neighbor, ELEM_O) >= 2) {
                    ewg_count++;
                }
            }

            /* 2+ EWG makes methylene acidic */
            if (ewg_count >= 2) {
                count++;
            }
        }

        /* Hydroxamic acid: R-CO-NH-OH or R-CO-N(R)-OH */
        if (atom->element == ELEM_N && !atom->aromatic) {
            bool has_carbonyl = false;
            bool has_oh = false;

            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];

                if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                    has_carbonyl = true;
                }
                if (neighbor->element == ELEM_O) {
                    int o_h = get_total_h(mol, neighbor);
                    if (o_h >= 1) {
                        has_oh = true;
                    }
                }
            }

            if (has_carbonyl && has_oh) {
                count++;  /* Hydroxamic acid */
            }
        }

        /* Oxime: C=N-OH */
        if (atom->element == ELEM_N && !atom->aromatic) {
            bool has_double_c = has_double_bond_to(mol, atom, ELEM_C);

            if (has_double_c) {
                for (int j = 0; j < atom->num_neighbors; j++) {
                    const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                    if (neighbor->element == ELEM_O) {
                        int o_h = get_total_h(mol, neighbor);
                        if (o_h >= 1) {
                            count++;  /* Oxime */
                            break;
                        }
                    }
                }
            }
        }

        /* Boronic acid: R-B(OH)2 (pKa ~8-9) */
        if (atom->element == ELEM_B) {
            int oh_count = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                if (neighbor->element == ELEM_O) {
                    int o_h = get_total_h(mol, neighbor);
                    if (o_h >= 1) {
                        oh_count++;
                    }
                }
            }
            if (oh_count >= 2) {
                count++;  /* Boronic acid */
            }
        }
    }

    return count;
}

/* ============================================================================
 * Basic Group Detection Functions
 * ============================================================================ */

/**
 * Super Basic Groups (conjugate acid pKa > 15):
 * - Alkoxides (charged): R-O(-) (~16-18)
 * - Amide anions (charged): R2N(-) (~35)
 * - Carbanions (charged)
 */
static int count_super_basic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Alkoxide: O with negative charge */
        if (atom->element == ELEM_O && atom->charge < 0) {
            int heavy = count_heavy_neighbors(mol, atom);
            if (heavy == 1) {
                const atom_t* neighbor = get_neighbor_of_element(mol, atom, ELEM_C);
                /* Exclude carboxylate (resonance stabilized) */
                if (neighbor && !has_double_bond_to(mol, neighbor, ELEM_O)) {
                    count++;
                }
            }
        }

        /* Amide anion: N with negative charge (not nitrile anion) */
        if (atom->element == ELEM_N && atom->charge < 0) {
            bool is_nitrile = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                int bond_idx = atom->neighbor_bonds[j];
                if (mol->bonds[bond_idx].type == BOND_TRIPLE) {
                    is_nitrile = true;
                    break;
                }
            }
            if (!is_nitrile) {
                count++;
            }
        }

        /* Carbanion */
        if (atom->element == ELEM_C && atom->charge < 0) {
            count++;
        }
    }

    return count;
}

/**
 * Strongly Basic Groups (conjugate acid pKa 10-14):
 * - Primary aliphatic amines: R-NH2 (~10.6)
 * - Secondary aliphatic amines: R2-NH (~11)
 * - Tertiary aliphatic amines: R3-N (~10)
 * - Cyclic aliphatic amines: piperidine (~11.2), pyrrolidine (~11.3)
 * - Guanidines: NC(=N)N (~13.6)
 * - Amidines: RC(=N)N (~12)
 * - DBU/DBN type bicyclic amidines (~12-13)
 * - Imidazoles (=N-) (~7, borderline, but counted here)
 */
static int count_strongly_basic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Non-aromatic nitrogen bases */
        if (atom->element == ELEM_N && !atom->aromatic && atom->charge == 0) {
            int heavy = count_heavy_neighbors(mol, atom);

            /* Identify nitrogen type */
            bool is_amide = false;
            bool is_carbamate = false;
            bool is_urea = false;
            bool is_sulfonamide = false;
            bool is_guanidine_n = is_guanidine_nitrogen(mol, atom);
            bool is_amidine_n = is_amidine_nitrogen(mol, atom);
            bool attached_to_aromatic = false;

            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];

                /* Attached to aromatic (aniline) */
                if (neighbor->aromatic) {
                    attached_to_aromatic = true;
                }

                /* Amide: N-C=O */
                if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                    /* Check if also has OR (carbamate) or NR (urea) */
                    int n_neighbor_n = count_neighbors_of_element(mol, neighbor, ELEM_N);
                    if (n_neighbor_n >= 2) {
                        is_urea = true;
                    } else {
                        is_amide = true;
                    }

                    for (int k = 0; k < neighbor->num_neighbors; k++) {
                        const atom_t* nn = &mol->atoms[neighbor->neighbors[k]];
                        if (nn->element == ELEM_O && nn != neighbor) {
                            int bond_idx = neighbor->neighbor_bonds[k];
                            if (mol->bonds[bond_idx].type != BOND_DOUBLE) {
                                is_carbamate = true;
                            }
                        }
                    }
                }

                /* Sulfonamide: N-SO2 */
                if (neighbor->element == ELEM_S) {
                    int so = count_double_bonds_to(mol, neighbor, ELEM_O);
                    if (so >= 2) is_sulfonamide = true;
                }
            }

            /* Guanidines and amidines are strongly basic */
            if (is_guanidine_n || is_amidine_n) {
                count++;
                continue;
            }

            /* Skip weakly basic nitrogen types */
            if (is_amide || is_carbamate || is_urea || is_sulfonamide) continue;
            if (attached_to_aromatic) continue;

            /* Aliphatic amine (primary, secondary, tertiary) */
            if (heavy >= 1 && heavy <= 3) {
                count++;
            }
        }

        /* Imidazole =N- (not NH) */
        if (atom->element == ELEM_N && atom->aromatic && atom->charge == 0) {
            int h = get_total_h(mol, atom);

            if (h == 0 && is_in_aromatic_ring_of_size(mol, i, 5)) {
                /* Check if in imidazole-type ring (has NH in same ring) */
                if (is_in_ring_with_nh(mol, i)) {
                    count++;  /* Imidazole =N- is basic (~7) */
                }
            }
        }
    }

    return count;
}

/**
 * Moderately Basic Groups (conjugate acid pKa 5-10):
 * - Pyridines (~5.2)
 * - Alkyl-substituted pyridines (~6-7.5)
 * - Aminopyridines (4-amino ~9.2, 2-amino ~6.7)
 * - Quinolines (~4.9), isoquinolines (~5.1)
 * - Morpholine (~8.4)
 * - Thiomorpholine (~8.7)
 * - 1,2,4-Triazoles (=N-) (~2.5-10)
 * - Pyrazoles (=N-) (~2.5)
 * - Oxazoles, thiazoles (~0.8-1.2, borderline)
 */
static int count_moderately_basic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Aromatic nitrogen (pyridine-like) without H */
        if (atom->element == ELEM_N && atom->aromatic && atom->charge == 0) {
            int h = get_total_h(mol, atom);

            if (h == 0) {
                /* Check ring size and type */
                bool is_6_ring = is_in_aromatic_ring_of_size(mol, i, 6);
                bool is_5_ring = is_in_aromatic_ring_of_size(mol, i, 5);

                /* 6-membered aromatic N (pyridine/quinoline type) */
                if (is_6_ring) {
                    count++;
                }
                /* 5-membered but NOT imidazole type (no NH in ring) */
                else if (is_5_ring) {
                    bool has_nh_in_ring = is_in_ring_with_nh(mol, i);
                    if (!has_nh_in_ring) {
                        /* Oxazole, thiazole, pyrazole =N- type */
                        count++;
                    }
                    /* Imidazole handled in strong */
                }
            }
        }

        /* Morpholine/thiomorpholine (saturated heterocycle with N) */
        if (atom->element == ELEM_N && !atom->aromatic && atom->charge == 0) {
            if (atom->ring_count > 0) {
                /* Check if it's a 6-membered saturated ring with O or S */
                for (int r = 0; r < mol->num_rings; r++) {
                    const ring_t* ring = &mol->rings[r];
                    if (ring->aromatic || ring->size != 6) continue;

                    bool has_this_n = false;
                    bool has_o_or_s = false;

                    for (int j = 0; j < ring->size; j++) {
                        if (ring->atoms[j] == i) has_this_n = true;
                        element_t e = mol->atoms[ring->atoms[j]].element;
                        if (e == ELEM_O || e == ELEM_S) has_o_or_s = true;
                    }

                    if (has_this_n && has_o_or_s) {
                        /* Check it's not already counted as strong base */
                        bool attached_to_aromatic = false;
                        bool is_amide = false;

                        for (int j = 0; j < atom->num_neighbors; j++) {
                            const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
                            if (neighbor->aromatic) attached_to_aromatic = true;
                            if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                                is_amide = true;
                            }
                        }

                        if (!attached_to_aromatic && !is_amide) {
                            count++;
                        }
                        break;
                    }
                }
            }
        }
    }

    return count;
}

/**
 * Weakly Basic Groups (conjugate acid pKa 0-5):
 * - Anilines (~4.6)
 * - Substituted anilines (EWG decrease, EDG increase basicity)
 * - N,N-Dialkylanilines (~5.1)
 * - Amides (as bases) (~-0.5)
 * - Ureas (~0.1)
 * - Ethers (R-O-R) (~-3.5)
 * - Esters (~-6)
 * - Ketones/aldehydes (~-7)
 * - Thioethers (~-6)
 * - Nitriles (~-10)
 * - Furan, thiophene, pyrrole (ring atoms)
 */
static int count_weakly_basic(const molecule_t* mol) {
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Aniline: N on aromatic ring */
        if (atom->element == ELEM_N && !atom->aromatic && atom->charge == 0) {
            bool attached_to_aromatic = false;
            bool is_amide = false;
            bool is_guanidine_n = is_guanidine_nitrogen(mol, atom);
            bool is_amidine_n = is_amidine_nitrogen(mol, atom);

            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];

                if (neighbor->element == ELEM_C && neighbor->aromatic) {
                    attached_to_aromatic = true;
                }
                if (neighbor->element == ELEM_C && has_double_bond_to(mol, neighbor, ELEM_O)) {
                    is_amide = true;
                }
            }

            /* Skip guanidines/amidines (strong) */
            if (is_guanidine_n || is_amidine_n) continue;

            if (attached_to_aromatic && !is_amide) {
                count++;  /* Aniline */
            } else if (is_amide) {
                count++;  /* Amide N as weak base */
            }
        }

        /* Ether oxygen (R-O-R) */
        if (atom->element == ELEM_O && !atom->aromatic && atom->charge == 0) {
            int h = get_total_h(mol, atom);
            int heavy = count_heavy_neighbors(mol, atom);

            bool is_carbonyl = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                int bond_idx = atom->neighbor_bonds[j];
                if (mol->bonds[bond_idx].type == BOND_DOUBLE) {
                    is_carbonyl = true;
                    break;
                }
            }

            if (!is_carbonyl && h == 0 && heavy == 2) {
                count++;  /* Ether */
            }
        }

        /* Carbonyl oxygen (ketone/aldehyde/ester) */
        if (atom->element == ELEM_O && !atom->aromatic) {
            for (int j = 0; j < atom->num_neighbors; j++) {
                int bond_idx = atom->neighbor_bonds[j];
                if (mol->bonds[bond_idx].type == BOND_DOUBLE) {
                    const atom_t* c = &mol->atoms[atom->neighbors[j]];
                    if (c->element == ELEM_C) {
                        count++;
                        break;
                    }
                }
            }
        }

        /* Thioether (R-S-R) */
        if (atom->element == ELEM_S && !atom->aromatic && atom->charge == 0) {
            int h = get_total_h(mol, atom);
            int heavy = count_heavy_neighbors(mol, atom);

            /* Check not sulfoxide/sulfone */
            int double_o = count_double_bonds_to(mol, atom, ELEM_O);

            if (h == 0 && heavy == 2 && double_o == 0) {
                count++;  /* Thioether */
            }
        }

        /* Nitrile nitrogen */
        if (atom->element == ELEM_N && !atom->aromatic) {
            for (int j = 0; j < atom->num_neighbors; j++) {
                int bond_idx = atom->neighbor_bonds[j];
                if (mol->bonds[bond_idx].type == BOND_TRIPLE) {
                    const atom_t* c = &mol->atoms[atom->neighbors[j]];
                    if (c->element == ELEM_C) {
                        count++;
                        break;
                    }
                }
            }
        }

        /* Aromatic O (furan) and S (thiophene) */
        if ((atom->element == ELEM_O || atom->element == ELEM_S) && atom->aromatic) {
            count++;  /* Very weak base */
        }
    }

    return count;
}

/* ============================================================================
 * Descriptor Compute Functions
 * ============================================================================ */

static cchem_status_t desc_super_acidic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_super_acidic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_strongly_acidic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_strongly_acidic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_moderately_acidic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_moderately_acidic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_weakly_acidic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_weakly_acidic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_total_acidic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_super_acidic(mol) + count_strongly_acidic(mol) +
               count_moderately_acidic(mol) + count_weakly_acidic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_mean_acidic_potential(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int n_super = count_super_acidic(mol);
    int n_strong = count_strongly_acidic(mol);
    int n_moderate = count_moderately_acidic(mol);
    int n_weak = count_weakly_acidic(mol);

    int total = n_super + n_strong + n_moderate + n_weak;

    if (total == 0) {
        value->d = 0.0;
    } else {
        double sum = n_super * PKA_SUPER_ACIDIC_MID +
                     n_strong * PKA_STRONG_ACIDIC_MID +
                     n_moderate * PKA_MODERATE_ACIDIC_MID +
                     n_weak * PKA_WEAK_ACIDIC_MID;
        value->d = sum / total;
    }

    return CCHEM_OK;
}

static cchem_status_t desc_super_basic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_super_basic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_strongly_basic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_strongly_basic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_moderately_basic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_moderately_basic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_weakly_basic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_weakly_basic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_total_basic_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = count_super_basic(mol) + count_strongly_basic(mol) +
               count_moderately_basic(mol) + count_weakly_basic(mol);
    return CCHEM_OK;
}

static cchem_status_t desc_mean_basic_potential(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int n_super = count_super_basic(mol);
    int n_strong = count_strongly_basic(mol);
    int n_moderate = count_moderately_basic(mol);
    int n_weak = count_weakly_basic(mol);

    int total = n_super + n_strong + n_moderate + n_weak;

    if (total == 0) {
        value->d = 0.0;
    } else {
        double sum = n_super * PKA_SUPER_BASIC_MID +
                     n_strong * PKA_STRONG_BASIC_MID +
                     n_moderate * PKA_MODERATE_BASIC_MID +
                     n_weak * PKA_WEAK_BASIC_MID;
        value->d = sum / total;
    }

    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_ACIDBASE_COUNT(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

#define REGISTER_ACIDBASE_DOUBLE(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_acidbase(void) {
    /* Acidic group counts */
    REGISTER_ACIDBASE_COUNT("SuperAcidicCount", "Super acidic groups (pKa < -1): sulfonic acids", desc_super_acidic_count);
    REGISTER_ACIDBASE_COUNT("StronglyAcidicCount", "Strongly acidic groups (pKa -1 to 1): phosphoric, haloacetic", desc_strongly_acidic_count);
    REGISTER_ACIDBASE_COUNT("ModeratelyAcidicCount", "Moderately acidic groups (pKa 1-5): carboxylic acids, tetrazoles", desc_moderately_acidic_count);
    REGISTER_ACIDBASE_COUNT("WeaklyAcidicCount", "Weakly acidic groups (pKa 5-14): phenols, thiols, sulfonamides", desc_weakly_acidic_count);
    REGISTER_ACIDBASE_COUNT("TotalAcidicFunctions", "Total count of all acidic functional groups", desc_total_acidic_count);
    REGISTER_ACIDBASE_DOUBLE("MeanAcidicPotential", "Mean pKa of acidic groups weighted by count", desc_mean_acidic_potential);

    /* Basic group counts */
    REGISTER_ACIDBASE_COUNT("SuperBasicCount", "Super basic groups (conj. pKa > 15): alkoxides, amide anions", desc_super_basic_count);
    REGISTER_ACIDBASE_COUNT("StronglyBasicCount", "Strongly basic groups (conj. pKa 10-14): aliphatic amines, guanidines", desc_strongly_basic_count);
    REGISTER_ACIDBASE_COUNT("ModeratelyBasicCount", "Moderately basic groups (conj. pKa 5-10): pyridines, morpholine", desc_moderately_basic_count);
    REGISTER_ACIDBASE_COUNT("WeaklyBasicCount", "Weakly basic groups (conj. pKa 0-5): anilines, amides, ethers", desc_weakly_basic_count);
    REGISTER_ACIDBASE_COUNT("TotalBasicFunctions", "Total count of all basic functional groups", desc_total_basic_count);
    REGISTER_ACIDBASE_DOUBLE("MeanBasicPotential", "Mean pKa of basic groups weighted by count", desc_mean_basic_potential);
}
