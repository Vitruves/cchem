/**
 * @file functional.c
 * @brief Comprehensive functional group descriptors for medicinal chemistry
 *
 * This module provides detection and counting of functional groups crucial
 * for QSAR, drug-likeness assessment, and medicinal chemistry ML models.
 *
 * Categories:
 * - Carbonyl variants (aldehydes, ketones, acyl halides, anhydrides)
 * - Nitrogen functional groups (imines, oximes, hydrazines, azides, etc.)
 * - Sulfur functional groups (sulfides, disulfides, sulfoxides, etc.)
 * - Oxygen functional groups (epoxides, peroxides, acetals, etc.)
 * - Mixed heteroatom groups (ureas, carbamates, isocyanates, etc.)
 * - Heterocyclic scaffolds (privileged medicinal chemistry rings)
 * - Aliphatic ring systems
 * - Drug-relevant features
 *
 * All functions are O(n) where n is atoms or bonds.
 */

#include <string.h>
#include <stdbool.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/* Count heavy (non-H) neighbors */
static int heavy_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) count++;
    }
    return count;
}

/* Get total H count (implicit + explicit) */
static int total_h(const molecule_t* mol, const atom_t* atom) {
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h++;
    }
    return h;
}

/* Check if atom has double bond to specific element */
static bool has_double_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE &&
            mol->atoms[atom->neighbors[i]].element == elem) {
            return true;
        }
    }
    return false;
}

/* Check if atom has single bond to specific element */
static bool has_single_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        bond_type_t bt = mol->bonds[bond_idx].type;
        if ((bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) &&
            mol->atoms[atom->neighbors[i]].element == elem) {
            return true;
        }
    }
    return false;
}

/* Check if atom has triple bond to specific element - reserved for future use */
__attribute__((unused))
static bool has_triple_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_TRIPLE &&
            mol->atoms[atom->neighbors[i]].element == elem) {
            return true;
        }
    }
    return false;
}

/* Count neighbors of specific element */
static int count_neighbors_elem(const molecule_t* mol, const atom_t* atom, element_t elem) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == elem) count++;
    }
    return count;
}

/* Count neighbors of specific element via single bond */
static int count_single_neighbors_elem(const molecule_t* mol, const atom_t* atom, element_t elem) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        bond_type_t bt = mol->bonds[bond_idx].type;
        if ((bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) &&
            mol->atoms[atom->neighbors[i]].element == elem) {
            count++;
        }
    }
    return count;
}

/* Check if ring contains only specific elements - reserved for future use */
__attribute__((unused))
static bool ring_contains_only(const molecule_t* mol, const ring_t* ring,
                               element_t e1, element_t e2, element_t e3) {
    for (int i = 0; i < ring->size; i++) {
        element_t e = mol->atoms[ring->atoms[i]].element;
        if (e != e1 && e != e2 && e != e3) return false;
    }
    return true;
}

/* Count atoms of element in ring */
static int ring_count_element(const molecule_t* mol, const ring_t* ring, element_t elem) {
    int count = 0;
    for (int i = 0; i < ring->size; i++) {
        if (mol->atoms[ring->atoms[i]].element == elem) count++;
    }
    return count;
}

/* Check if all atoms in ring are aromatic - reserved for future use */
__attribute__((unused))
static bool ring_all_aromatic(const molecule_t* mol, const ring_t* ring) {
    for (int i = 0; i < ring->size; i++) {
        if (!mol->atoms[ring->atoms[i]].aromatic) return false;
    }
    return true;
}

/* Check if all atoms in ring are non-aromatic */
static bool ring_all_aliphatic(const molecule_t* mol, const ring_t* ring) {
    for (int i = 0; i < ring->size; i++) {
        if (mol->atoms[ring->atoms[i]].aromatic) return false;
    }
    return true;
}

/* ============================================================================
 * Carbonyl Variants
 * ============================================================================ */

/* Aldehyde count: CHO groups (R-CHO where R is C or H, includes formaldehyde) */
static cchem_status_t desc_aldehyde_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!has_double_to(mol, atom, ELEM_O)) continue;
        int h = total_h(mol, atom);
        /* Aldehyde: C(=O)H with at least 1 H on carbonyl C (h >= 1) */
        /* Formaldehyde H2C=O has h=2, acetaldehyde CH3-CHO has h=1 */
        if (h >= 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Ketone count: C(=O)C where C is not in carboxyl, ester, or amide */
static cchem_status_t desc_ketone_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!has_double_to(mol, atom, ELEM_O)) continue;

        /* Check it's a ketone: 2 C neighbors, no O/N single bonds */
        int c_count = count_single_neighbors_elem(mol, atom, ELEM_C);
        int o_count = count_single_neighbors_elem(mol, atom, ELEM_O);
        int n_count = count_single_neighbors_elem(mol, atom, ELEM_N);
        int h = total_h(mol, atom);

        if (c_count == 2 && o_count == 0 && n_count == 0 && h == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Acyl halide count: C(=O)X where X = F, Cl, Br, I */
static cchem_status_t desc_acyl_halide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!has_double_to(mol, atom, ELEM_O)) continue;

        if (has_single_to(mol, atom, ELEM_F) ||
            has_single_to(mol, atom, ELEM_Cl) ||
            has_single_to(mol, atom, ELEM_Br) ||
            has_single_to(mol, atom, ELEM_I)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Acid anhydride count: C(=O)-O-C(=O) */
static cchem_status_t desc_anhydride_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    /* Look for bridging O between two carbonyls */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;
        if (atom->aromatic) continue;

        int carbonyl_c = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_C && has_double_to(mol, nb, ELEM_O)) {
                carbonyl_c++;
            }
        }
        if (carbonyl_c == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Nitrogen Functional Groups
 * ============================================================================ */

/* Imine count: C=N-C */
static cchem_status_t desc_imine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_DOUBLE) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if ((a1->element == ELEM_C && a2->element == ELEM_N && !a2->aromatic) ||
            (a1->element == ELEM_N && a2->element == ELEM_C && !a1->aromatic)) {
            const atom_t* n_atom = (a1->element == ELEM_N) ? a1 : a2;
            /* Imine: N has C neighbor (not just =C) */
            if (has_single_to(mol, n_atom, ELEM_C) || total_h(mol, n_atom) > 0) {
                count++;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Oxime count: C=N-OH */
static cchem_status_t desc_oxime_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->aromatic) continue;
        if (!has_double_to(mol, atom, ELEM_C)) continue;

        /* Check for O-H attached to N */
        for (int j = 0; j < atom->num_neighbors; j++) {
            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_O && total_h(mol, nb) > 0) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Hydrazine count: N-N with at least one H on each */
static cchem_status_t desc_hydrazine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        bond_type_t bt = bond->type;
        if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_N && a2->element == ELEM_N &&
            !a1->aromatic && !a2->aromatic) {
            if (total_h(mol, a1) > 0 && total_h(mol, a2) > 0) {
                count++;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Hydrazone count: C=N-N */
static cchem_status_t desc_hydrazone_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->aromatic) continue;
        if (!has_double_to(mol, atom, ELEM_C)) continue;
        if (has_single_to(mol, atom, ELEM_N)) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Azo count: N=N */
static cchem_status_t desc_azo_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_DOUBLE) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_N && a2->element == ELEM_N) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Azide count: N=N=N or N#N=N */
static cchem_status_t desc_azide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;

        /* Central N of azide: two N neighbors via double/triple bonds */
        int n_neighbors = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if ((bt == BOND_DOUBLE || bt == BOND_TRIPLE) &&
                mol->atoms[atom->neighbors[j]].element == ELEM_N) {
                n_neighbors++;
            }
        }
        if (n_neighbors == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Diazo count: C=N=N (diazo compound) */
static cchem_status_t desc_diazo_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;

        /* Central N: double to C and double to N */
        if (has_double_to(mol, atom, ELEM_C) && has_double_to(mol, atom, ELEM_N)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Guanidine count: NC(=N)N */
static cchem_status_t desc_guanidine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        int n_single = count_single_neighbors_elem(mol, atom, ELEM_N);
        bool n_double = has_double_to(mol, atom, ELEM_N);

        if (n_single == 2 && n_double) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Amidine count: C(=N)-N */
static cchem_status_t desc_amidine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        bool n_double = has_double_to(mol, atom, ELEM_N);
        int n_single = count_single_neighbors_elem(mol, atom, ELEM_N);

        if (n_double && n_single == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Quaternary ammonium count: N+ with 4 C neighbors */
static cchem_status_t desc_quat_ammonium_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->charge <= 0) continue;

        int c_count = count_neighbors_elem(mol, atom, ELEM_C);
        if (c_count == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* N-oxide count: N+ with O- neighbor */
static cchem_status_t desc_n_oxide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->charge <= 0) continue;

        for (int j = 0; j < atom->num_neighbors; j++) {
            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_O && nb->charge < 0) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Sulfur Functional Groups
 * ============================================================================ */

/* Sulfide count: C-S-C (excluding sulfoxides and sulfones) */
static cchem_status_t desc_sulfide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;
        if (atom->aromatic) continue;

        /* Check for double bonds to O (would be sulfoxide/sulfone) */
        bool has_double_o = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            if (mol->bonds[bond_idx].type == BOND_DOUBLE &&
                mol->atoms[atom->neighbors[j]].element == ELEM_O) {
                has_double_o = true;
                break;
            }
        }
        if (has_double_o) continue;  /* Skip sulfoxides/sulfones */

        int c_count = count_single_neighbors_elem(mol, atom, ELEM_C);
        int h = total_h(mol, atom);
        if (c_count == 2 && h == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Disulfide count: S-S bonds */
static cchem_status_t desc_disulfide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        bond_type_t bt = bond->type;
        if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_S && a2->element == ELEM_S) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Sulfoxide count: S=O (but not SO2) */
static cchem_status_t desc_sulfoxide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        int double_o = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            if (mol->bonds[bond_idx].type == BOND_DOUBLE &&
                mol->atoms[atom->neighbors[j]].element == ELEM_O) {
                double_o++;
            }
        }
        if (double_o == 1) count++;  /* Exactly one =O */
    }
    value->i = count;
    return CCHEM_OK;
}

/* Sulfonamide count: SO2-N */
static cchem_status_t desc_sulfonamide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        int double_o = 0;
        bool has_n = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            element_t elem = mol->atoms[atom->neighbors[j]].element;
            if (mol->bonds[bond_idx].type == BOND_DOUBLE && elem == ELEM_O) {
                double_o++;
            }
            if (elem == ELEM_N) has_n = true;
        }
        if (double_o == 2 && has_n) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Sulfonic acid count: SO3H */
static cchem_status_t desc_sulfonic_acid_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        int double_o = 0;
        bool has_oh = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (mol->bonds[bond_idx].type == BOND_DOUBLE && nb->element == ELEM_O) {
                double_o++;
            }
            if (nb->element == ELEM_O && total_h(mol, nb) > 0) {
                has_oh = true;
            }
        }
        if (double_o >= 2 && has_oh) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Sulfonate ester count: SO2-O-C */
static cchem_status_t desc_sulfonate_ester_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        int double_o = 0;
        int o_c_bridge = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt == BOND_DOUBLE && nb->element == ELEM_O) {
                double_o++;
            }
            if ((bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) &&
                nb->element == ELEM_O && total_h(mol, nb) == 0 &&
                has_single_to(mol, nb, ELEM_C)) {
                o_c_bridge++;
            }
        }
        if (double_o >= 2 && o_c_bridge >= 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Oxygen Functional Groups
 * ============================================================================ */

/* Epoxide count: 3-membered ring with O */
static cchem_status_t desc_epoxide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 3) continue;
        if (ring_count_element(mol, ring, ELEM_O) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 2) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Peroxide count: O-O bonds */
static cchem_status_t desc_peroxide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        bond_type_t bt = bond->type;
        if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_O && a2->element == ELEM_O) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Acetal count: C with 2 O-C neighbors (R-CH(OR)2) */
static cchem_status_t desc_acetal_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->aromatic) continue;

        int or_count = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_O && has_single_to(mol, nb, ELEM_C)) {
                or_count++;
            }
        }
        if (or_count == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Hemiacetal count: C with OH and OR (R-CH(OH)(OR)) */
static cchem_status_t desc_hemiacetal_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->aromatic) continue;

        int or_count = 0;
        int oh_count = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_O) {
                if (total_h(mol, nb) > 0) oh_count++;
                else if (has_single_to(mol, nb, ELEM_C)) or_count++;
            }
        }
        if (or_count == 1 && oh_count == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Enol count: C=C-OH */
static cchem_status_t desc_enol_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->aromatic) continue;

        /* Check for C=C */
        bool has_c_double = has_double_to(mol, atom, ELEM_C);
        if (!has_c_double) continue;

        /* Check for OH on this carbon */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_O && total_h(mol, nb) > 0) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Mixed Heteroatom Groups
 * ============================================================================ */

/* Urea count: NC(=O)N */
static cchem_status_t desc_urea_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        bool o_double = has_double_to(mol, atom, ELEM_O);
        int n_single = count_single_neighbors_elem(mol, atom, ELEM_N);

        if (o_double && n_single == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Thiourea count: NC(=S)N */
static cchem_status_t desc_thiourea_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        bool s_double = has_double_to(mol, atom, ELEM_S);
        int n_single = count_single_neighbors_elem(mol, atom, ELEM_N);

        if (s_double && n_single == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Carbamate count: NC(=O)O-C */
static cchem_status_t desc_carbamate_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        bool o_double = has_double_to(mol, atom, ELEM_O);
        int n_single = count_single_neighbors_elem(mol, atom, ELEM_N);

        if (!o_double || n_single == 0) continue;

        /* Check for O-C (ester-like O) */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_O && total_h(mol, nb) == 0 &&
                has_single_to(mol, nb, ELEM_C)) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Isocyanate count: N=C=O */
static cchem_status_t desc_isocyanate_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        if (has_double_to(mol, atom, ELEM_N) && has_double_to(mol, atom, ELEM_O)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Isothiocyanate count: N=C=S */
static cchem_status_t desc_isothiocyanate_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        if (has_double_to(mol, atom, ELEM_N) && has_double_to(mol, atom, ELEM_S)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Phosphate count: P with O neighbors */
static cchem_status_t desc_phosphate_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_P) continue;

        int o_count = count_neighbors_elem(mol, atom, ELEM_O);
        if (o_count >= 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Phosphonate count: P with C and O neighbors */
static cchem_status_t desc_phosphonate_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_P) continue;

        int o_count = count_neighbors_elem(mol, atom, ELEM_O);
        int c_count = count_neighbors_elem(mol, atom, ELEM_C);
        if (o_count >= 2 && c_count >= 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Heterocyclic Scaffolds - Privileged Structures in Medicinal Chemistry
 * ============================================================================ */

/* Benzene ring count (6-membered aromatic carbocycle) */
static cchem_status_t desc_benzene_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 6) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_C) == 6) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Pyridine count (6-membered aromatic with 1 N) */
static cchem_status_t desc_pyridine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 6) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 5) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Pyrimidine count (6-membered aromatic with 2 N at 1,3) */
static cchem_status_t desc_pyrimidine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 6) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 2 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Pyrrole count (5-membered aromatic with 1 N) */
static cchem_status_t desc_pyrrole_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Furan count (5-membered aromatic with 1 O) */
static cchem_status_t desc_furan_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_O) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Thiophene count (5-membered aromatic with 1 S) */
static cchem_status_t desc_thiophene_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_S) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Imidazole count (5-membered aromatic with 2 N) */
static cchem_status_t desc_imidazole_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 2 &&
            ring_count_element(mol, ring, ELEM_C) == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Oxazole count (5-membered aromatic with 1 N, 1 O) */
static cchem_status_t desc_oxazole_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_O) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Thiazole count (5-membered aromatic with 1 N, 1 S) */
static cchem_status_t desc_thiazole_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_S) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Triazole count (5-membered aromatic with 3 N) */
static cchem_status_t desc_triazole_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 3 &&
            ring_count_element(mol, ring, ELEM_C) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Tetrazole count (5-membered aromatic with 4 N) */
static cchem_status_t desc_tetrazole_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring->aromatic) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 4 &&
            ring_count_element(mol, ring, ELEM_C) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Aliphatic Ring Systems
 * ============================================================================ */

/* Cyclopropane count */
static cchem_status_t desc_cyclopropane_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 3) continue;
        if (ring_all_aliphatic(mol, ring) &&
            ring_count_element(mol, ring, ELEM_C) == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Cyclobutane count */
static cchem_status_t desc_cyclobutane_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 4) continue;
        if (ring_all_aliphatic(mol, ring) &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Cyclopentane count */
static cchem_status_t desc_cyclopentane_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (ring_all_aliphatic(mol, ring) &&
            ring_count_element(mol, ring, ELEM_C) == 5) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Cyclohexane count */
static cchem_status_t desc_cyclohexane_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 6) continue;
        if (ring_all_aliphatic(mol, ring) &&
            ring_count_element(mol, ring, ELEM_C) == 6) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Piperidine count (6-membered aliphatic with 1 N) */
static cchem_status_t desc_piperidine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 6) continue;
        if (!ring_all_aliphatic(mol, ring)) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 5) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Piperazine count (6-membered aliphatic with 2 N) */
static cchem_status_t desc_piperazine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 6) continue;
        if (!ring_all_aliphatic(mol, ring)) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 2 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Morpholine count (6-membered aliphatic with 1 N, 1 O) */
static cchem_status_t desc_morpholine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 6) continue;
        if (!ring_all_aliphatic(mol, ring)) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_O) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Pyrrolidine count (5-membered aliphatic with 1 N) */
static cchem_status_t desc_pyrrolidine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring_all_aliphatic(mol, ring)) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Tetrahydrofuran count (5-membered aliphatic with 1 O) */
static cchem_status_t desc_thf_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 5) continue;
        if (!ring_all_aliphatic(mol, ring)) continue;
        if (ring_count_element(mol, ring, ELEM_O) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Aziridine count (3-membered ring with 1 N) */
static cchem_status_t desc_aziridine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->size != 3) continue;
        if (ring_count_element(mol, ring, ELEM_N) == 1 &&
            ring_count_element(mol, ring, ELEM_C) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Drug-Relevant Features
 * ============================================================================ */

/* Lactam count: cyclic amide (C(=O)N in ring) */
static cchem_status_t desc_lactam_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->ring_count == 0) continue;
        if (!has_double_to(mol, atom, ELEM_O)) continue;

        /* Check for N in same ring */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_N && nb->ring_count > 0) {
                /* Check if they share a ring */
                for (int r = 0; r < mol->num_rings; r++) {
                    const ring_t* ring = &mol->rings[r];
                    bool has_c = false, has_n = false;
                    for (int k = 0; k < ring->size; k++) {
                        if (ring->atoms[k] == i) has_c = true;
                        if (ring->atoms[k] == atom->neighbors[j]) has_n = true;
                    }
                    if (has_c && has_n) {
                        count++;
                        goto next_atom;
                    }
                }
            }
        }
        next_atom:;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Lactone count: cyclic ester (C(=O)O in ring) */
static cchem_status_t desc_lactone_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->ring_count == 0) continue;
        if (!has_double_to(mol, atom, ELEM_O)) continue;

        /* Check for O in same ring via single bond */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_O && nb->ring_count > 0 && total_h(mol, nb) == 0) {
                /* Check if they share a ring */
                for (int r = 0; r < mol->num_rings; r++) {
                    const ring_t* ring = &mol->rings[r];
                    bool has_c = false, has_o = false;
                    for (int k = 0; k < ring->size; k++) {
                        if (ring->atoms[k] == i) has_c = true;
                        if (ring->atoms[k] == atom->neighbors[j]) has_o = true;
                    }
                    if (has_c && has_o) {
                        count++;
                        goto next_atom2;
                    }
                }
            }
        }
        next_atom2:;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Beta-lactam count: 4-membered lactam */
static cchem_status_t desc_beta_lactam_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (ring->size != 4) continue;

        /* Look for C(=O)-N pattern in ring */
        for (int j = 0; j < ring->size; j++) {
            const atom_t* atom = &mol->atoms[ring->atoms[j]];
            if (atom->element != ELEM_C) continue;
            if (!has_double_to(mol, atom, ELEM_O)) continue;

            /* Check if N is adjacent in ring */
            int prev = (j + ring->size - 1) % ring->size;
            int next = (j + 1) % ring->size;

            if (mol->atoms[ring->atoms[prev]].element == ELEM_N ||
                mol->atoms[ring->atoms[next]].element == ELEM_N) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Vinyl group count: C=C terminal */
static cchem_status_t desc_vinyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_DOUBLE) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element != ELEM_C || a2->element != ELEM_C) continue;
        if (a1->aromatic || a2->aromatic) continue;

        /* At least one end is =CH2 */
        if (total_h(mol, a1) == 2 || total_h(mol, a2) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Allyl group count: C=C-C where terminal C has substituent */
static cchem_status_t desc_allyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_DOUBLE) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element != ELEM_C || a2->element != ELEM_C) continue;
        if (a1->aromatic || a2->aromatic) continue;

        /* Check for CH2 attached to either end */
        for (int j = 0; j < a1->num_neighbors; j++) {
            int bidx = a1->neighbor_bonds[j];
            if (mol->bonds[bidx].type != BOND_SINGLE) continue;
            const atom_t* nb = &mol->atoms[a1->neighbors[j]];
            if (nb->element == ELEM_C && !nb->aromatic && total_h(mol, nb) == 2) {
                count++;
                break;
            }
        }
        for (int j = 0; j < a2->num_neighbors; j++) {
            int bidx = a2->neighbor_bonds[j];
            if (mol->bonds[bidx].type != BOND_SINGLE) continue;
            const atom_t* nb = &mol->atoms[a2->neighbors[j]];
            if (nb->element == ELEM_C && !nb->aromatic && total_h(mol, nb) == 2) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Benzyl group count: aromatic C attached to CH2 */
static cchem_status_t desc_benzyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!atom->aromatic) continue;

        /* Check for CH2 neighbor */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (nb->element == ELEM_C && !nb->aromatic && total_h(mol, nb) == 2) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Phenyl group count: aromatic 6-ring attached to something */
static cchem_status_t desc_phenyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;

    /* Count aromatic carbons that have non-aromatic neighbors */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!atom->aromatic) continue;

        for (int j = 0; j < atom->num_neighbors; j++) {
            const atom_t* nb = &mol->atoms[atom->neighbors[j]];
            if (!nb->aromatic && nb->element != ELEM_H) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_FUNC_DESC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_COUNTS; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_functional(void) {
    /* Carbonyl variants */
    REGISTER_FUNC_DESC("AldehydeCount", "Aldehyde groups (-CHO)", desc_aldehyde_count);
    REGISTER_FUNC_DESC("KetoneCount", "Ketone groups (C(=O)C)", desc_ketone_count);
    REGISTER_FUNC_DESC("AcylHalideCount", "Acyl halide groups (C(=O)X)", desc_acyl_halide_count);
    REGISTER_FUNC_DESC("AnhydrideCount", "Acid anhydride groups", desc_anhydride_count);

    /* Nitrogen functional groups */
    REGISTER_FUNC_DESC("ImineCount", "Imine groups (C=N)", desc_imine_count);
    REGISTER_FUNC_DESC("OximeCount", "Oxime groups (C=N-OH)", desc_oxime_count);
    REGISTER_FUNC_DESC("HydrazineCount", "Hydrazine groups (N-N with H)", desc_hydrazine_count);
    REGISTER_FUNC_DESC("HydrazoneCount", "Hydrazone groups (C=N-N)", desc_hydrazone_count);
    REGISTER_FUNC_DESC("AzoCount", "Azo groups (N=N)", desc_azo_count);
    REGISTER_FUNC_DESC("AzideCount", "Azide groups (N=N=N)", desc_azide_count);
    REGISTER_FUNC_DESC("DiazoCount", "Diazo groups (C=N=N)", desc_diazo_count);
    REGISTER_FUNC_DESC("GuanidineCount", "Guanidine groups (NC(=N)N)", desc_guanidine_count);
    REGISTER_FUNC_DESC("AmidineCount", "Amidine groups (C(=N)N)", desc_amidine_count);
    REGISTER_FUNC_DESC("QuatAmmoniumCount", "Quaternary ammonium (N+)", desc_quat_ammonium_count);
    REGISTER_FUNC_DESC("NOxideCount", "N-oxide groups", desc_n_oxide_count);

    /* Sulfur functional groups */
    REGISTER_FUNC_DESC("SulfideCount", "Sulfide linkages (C-S-C)", desc_sulfide_count);
    REGISTER_FUNC_DESC("DisulfideCount", "Disulfide bonds (S-S)", desc_disulfide_count);
    REGISTER_FUNC_DESC("SulfoxideCount", "Sulfoxide groups (S=O)", desc_sulfoxide_count);
    REGISTER_FUNC_DESC("SulfonamideCount", "Sulfonamide groups (SO2N)", desc_sulfonamide_count);
    REGISTER_FUNC_DESC("SulfonicAcidCount", "Sulfonic acid groups (SO3H)", desc_sulfonic_acid_count);
    REGISTER_FUNC_DESC("SulfonateEsterCount", "Sulfonate ester groups", desc_sulfonate_ester_count);

    /* Oxygen functional groups */
    REGISTER_FUNC_DESC("EpoxideCount", "Epoxide rings (3-ring with O)", desc_epoxide_count);
    REGISTER_FUNC_DESC("PeroxideCount", "Peroxide bonds (O-O)", desc_peroxide_count);
    REGISTER_FUNC_DESC("AcetalCount", "Acetal groups (C(OR)2)", desc_acetal_count);
    REGISTER_FUNC_DESC("HemiacetalCount", "Hemiacetal groups (C(OH)(OR))", desc_hemiacetal_count);
    REGISTER_FUNC_DESC("EnolCount", "Enol groups (C=C-OH)", desc_enol_count);

    /* Mixed heteroatom groups */
    REGISTER_FUNC_DESC("UreaCount", "Urea groups (NC(=O)N)", desc_urea_count);
    REGISTER_FUNC_DESC("ThioureaCount", "Thiourea groups (NC(=S)N)", desc_thiourea_count);
    REGISTER_FUNC_DESC("CarbamateCount", "Carbamate groups (NC(=O)O)", desc_carbamate_count);
    REGISTER_FUNC_DESC("IsocyanateCount", "Isocyanate groups (N=C=O)", desc_isocyanate_count);
    REGISTER_FUNC_DESC("IsothiocyanateCount", "Isothiocyanate groups (N=C=S)", desc_isothiocyanate_count);
    REGISTER_FUNC_DESC("PhosphateCount", "Phosphate groups (PO4)", desc_phosphate_count);
    REGISTER_FUNC_DESC("PhosphonateCount", "Phosphonate groups (P(O)(O)C)", desc_phosphonate_count);

    /* Aromatic heterocyclic scaffolds */
    REGISTER_FUNC_DESC("BenzeneCount", "Benzene rings (C6 aromatic)", desc_benzene_count);
    REGISTER_FUNC_DESC("PyridineCount", "Pyridine rings (6-ring, 1N)", desc_pyridine_count);
    REGISTER_FUNC_DESC("PyrimidineCount", "Pyrimidine rings (6-ring, 2N)", desc_pyrimidine_count);
    REGISTER_FUNC_DESC("PyrroleCount", "Pyrrole rings (5-ring, 1N)", desc_pyrrole_count);
    REGISTER_FUNC_DESC("FuranCount", "Furan rings (5-ring, 1O)", desc_furan_count);
    REGISTER_FUNC_DESC("ThiopheneCount", "Thiophene rings (5-ring, 1S)", desc_thiophene_count);
    REGISTER_FUNC_DESC("ImidazoleCount", "Imidazole rings (5-ring, 2N)", desc_imidazole_count);
    REGISTER_FUNC_DESC("OxazoleCount", "Oxazole rings (5-ring, N+O)", desc_oxazole_count);
    REGISTER_FUNC_DESC("ThiazoleCount", "Thiazole rings (5-ring, N+S)", desc_thiazole_count);
    REGISTER_FUNC_DESC("TriazoleCount", "Triazole rings (5-ring, 3N)", desc_triazole_count);
    REGISTER_FUNC_DESC("TetrazoleCount", "Tetrazole rings (5-ring, 4N)", desc_tetrazole_count);

    /* Aliphatic ring systems */
    REGISTER_FUNC_DESC("CyclopropaneCount", "Cyclopropane rings", desc_cyclopropane_count);
    REGISTER_FUNC_DESC("CyclobutaneCount", "Cyclobutane rings", desc_cyclobutane_count);
    REGISTER_FUNC_DESC("CyclopentaneCount", "Cyclopentane rings", desc_cyclopentane_count);
    REGISTER_FUNC_DESC("CyclohexaneCount", "Cyclohexane rings", desc_cyclohexane_count);
    REGISTER_FUNC_DESC("PiperidineCount", "Piperidine rings (6-ring, 1N sat)", desc_piperidine_count);
    REGISTER_FUNC_DESC("PiperazineCount", "Piperazine rings (6-ring, 2N sat)", desc_piperazine_count);
    REGISTER_FUNC_DESC("MorpholineCount", "Morpholine rings (6-ring, N+O sat)", desc_morpholine_count);
    REGISTER_FUNC_DESC("PyrrolidineCount", "Pyrrolidine rings (5-ring, 1N sat)", desc_pyrrolidine_count);
    REGISTER_FUNC_DESC("THFCount", "Tetrahydrofuran rings (5-ring, 1O sat)", desc_thf_count);
    REGISTER_FUNC_DESC("AziridineCount", "Aziridine rings (3-ring, 1N)", desc_aziridine_count);

    /* Drug-relevant features */
    REGISTER_FUNC_DESC("LactamCount", "Lactam rings (cyclic amide)", desc_lactam_count);
    REGISTER_FUNC_DESC("LactoneCount", "Lactone rings (cyclic ester)", desc_lactone_count);
    REGISTER_FUNC_DESC("BetaLactamCount", "Beta-lactam rings (4-ring lactam)", desc_beta_lactam_count);
    REGISTER_FUNC_DESC("VinylCount", "Vinyl groups (C=CH2)", desc_vinyl_count);
    REGISTER_FUNC_DESC("AllylCount", "Allyl groups (C=C-CH2)", desc_allyl_count);
    REGISTER_FUNC_DESC("BenzylCount", "Benzyl groups (Ar-CH2)", desc_benzyl_count);
    REGISTER_FUNC_DESC("PhenylCount", "Phenyl attachments (Ar-R)", desc_phenyl_count);
}
