/**
 * @file pharmacophore.c
 * @brief Pharmacophore and molecular complexity descriptors
 *
 * Pharmacophore features critical for drug discovery:
 * - Lipophilic centers
 * - Positive ionizable (basic nitrogen)
 * - Negative ionizable (acidic groups)
 * - Halogen bond donors
 * - Aromatic π-stacking
 *
 * Molecular complexity measures:
 * - Bertz complexity index
 * - Bottcher complexity
 * - Information content
 * - Stereocenter complexity
 * - Drug-likeness scores
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Helper functions
 * ============================================================================ */

/* Count heavy (non-H) neighbors */
static int heavy_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            count++;
        }
    }
    return count;
}

/* Total hydrogens attached (explicit + implicit) */
static int total_h(const molecule_t* mol, const atom_t* atom) {
    int explicit_h = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) {
            explicit_h++;
        }
    }
    return explicit_h + atom->implicit_h_count;
}

/* Check if atom has double bond to element */
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

/* Check if atom has single bond to element */
static bool has_single_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_SINGLE &&
            mol->atoms[atom->neighbors[i]].element == elem) {
            return true;
        }
    }
    return false;
}

/* Check if atom is sp3 carbon */
static bool is_sp3_carbon(const molecule_t* mol, const atom_t* atom) {
    if (atom->element != ELEM_C) return false;
    if (atom->aromatic) return false;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE ||
            mol->bonds[bond_idx].type == BOND_TRIPLE) {
            return false;
        }
    }
    return true;
}

/* Check if ring is aromatic */
static bool ring_is_aromatic(const molecule_t* mol, const ring_t* ring) {
    for (int i = 0; i < ring->size; i++) {
        if (!mol->atoms[ring->atoms[i]].aromatic) {
            return false;
        }
    }
    return true;
}

/* Check if ring contains element - reserved for future use */
__attribute__((unused))
static int ring_count_element(const molecule_t* mol, const ring_t* ring, element_t elem) {
    int count = 0;
    for (int i = 0; i < ring->size; i++) {
        if (mol->atoms[ring->atoms[i]].element == elem) {
            count++;
        }
    }
    return count;
}

/* ============================================================================
 * Pharmacophore Point Descriptors
 * ============================================================================ */

/* Lipophilic center count: sp3 carbons with only C/H neighbors */
static cchem_status_t desc_lipophilic_centers(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (!is_sp3_carbon(mol, atom)) continue;

        /* Check all neighbors are C or H */
        bool lipophilic = true;
        for (int j = 0; j < atom->num_neighbors; j++) {
            element_t elem = mol->atoms[atom->neighbors[j]].element;
            if (elem != ELEM_C && elem != ELEM_H) {
                lipophilic = false;
                break;
            }
        }
        if (lipophilic) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Basic nitrogen count (protonatable at pH 7): amines without EWG, guanidines, amidines */
static cchem_status_t desc_basic_nitrogen_count(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;

        int h = total_h(mol, atom);
        int heavy = heavy_neighbors(mol, atom);

        /* Skip if it has no lone pair available (N+, N-oxide) */
        if (atom->charge > 0) continue;  /* Already charged */

        /* Aromatic N with H (like pyrrole) - less basic */
        if (atom->aromatic && h >= 1) continue;

        /* Check for EWG: carbonyl, sulfonyl, nitro attached */
        bool has_ewg = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
            if (neighbor->element == ELEM_C) {
                /* Check for C=O (amide) */
                if (has_double_to(mol, neighbor, ELEM_O)) {
                    has_ewg = true;
                    break;
                }
            } else if (neighbor->element == ELEM_S) {
                /* Sulfonamide */
                if (has_double_to(mol, neighbor, ELEM_O)) {
                    has_ewg = true;
                    break;
                }
            } else if (neighbor->element == ELEM_N && has_double_to(mol, neighbor, ELEM_O)) {
                /* Nitro */
                has_ewg = true;
                break;
            }
        }
        if (has_ewg) continue;

        /* Primary, secondary, tertiary amines; sp2 N in rings like pyridine */
        if (h >= 1 || (!atom->aromatic && heavy <= 3) ||
            (atom->aromatic && h == 0 && heavy == 2)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Acidic oxygen count: carboxylic, sulfonic, phosphoric acids */
static cchem_status_t desc_acidic_oxygen_count(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;
        if (atom->charge == -1) {
            count++;  /* Already deprotonated */
            continue;
        }

        int h = total_h(mol, atom);
        if (h == 0) continue;  /* No acidic H */

        /* Check if attached to C, S, or P with =O */
        for (int j = 0; j < atom->num_neighbors; j++) {
            const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
            if ((neighbor->element == ELEM_C || neighbor->element == ELEM_S ||
                 neighbor->element == ELEM_P) && has_double_to(mol, neighbor, ELEM_O)) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Halogen bond donor count: C-X where X = Cl, Br, I (perfluorinated preferred) */
static cchem_status_t desc_halogen_bond_donor(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        /* Only Cl, Br, I can be XB donors (F is too electronegative) */
        if (atom->element != ELEM_Cl && atom->element != ELEM_Br &&
            atom->element != ELEM_I) continue;

        /* Must be attached to carbon */
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element == ELEM_C) {
                count++;
                break;
            }
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Aromatic ring count (π-stacking capable) */
static cchem_status_t desc_aromatic_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (ring_is_aromatic(mol, &mol->rings[i])) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Positive charge center count (formal positive charges) */
static cchem_status_t desc_positive_charge_centers(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].charge > 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Negative charge center count (formal negative charges) */
static cchem_status_t desc_negative_charge_centers(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].charge < 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Total pharmacophore points (HBD + HBA + Lipophilic + Charged) */
static cchem_status_t desc_total_pharmacophore_points(const molecule_t* mol, descriptor_value_t* value) {
    int hbd = 0, hba = 0, lipo = 0, charged = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        int h = total_h(mol, atom);

        /* HBD: OH, NH */
        if ((atom->element == ELEM_O || atom->element == ELEM_N) && h > 0) {
            hbd++;
        }

        /* HBA: O, N with lone pair */
        if (atom->element == ELEM_O || atom->element == ELEM_N) {
            hba++;
        }

        /* Lipophilic: C attached only to C/H */
        if (is_sp3_carbon(mol, atom)) {
            bool all_ch = true;
            for (int j = 0; j < atom->num_neighbors; j++) {
                element_t e = mol->atoms[atom->neighbors[j]].element;
                if (e != ELEM_C && e != ELEM_H) {
                    all_ch = false;
                    break;
                }
            }
            if (all_ch) lipo++;
        }

        /* Charged */
        if (atom->charge != 0) charged++;
    }

    value->i = hbd + hba + lipo + charged;
    return CCHEM_OK;
}

/* Pharmacophore density (points per heavy atom) */
static cchem_status_t desc_pharmacophore_density(const molecule_t* mol, descriptor_value_t* value) {
    descriptor_value_t total;
    desc_total_pharmacophore_points(mol, &total);

    int heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) heavy++;
    }

    value->d = heavy > 0 ? (double)total.i / heavy : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Molecular Complexity Descriptors
 * ============================================================================ */

/* Bertz complexity index: CT = log2(1 + Σ(bond orders) + Σ(atom types)) */
static cchem_status_t desc_bertz_complexity(const molecule_t* mol, descriptor_value_t* value) {
    if (mol->num_atoms == 0) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    /* Count distinct atom types */
    int type_counts[128] = {0};  /* element * 4 + hybridization */
    int bond_sum = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        int hyb = 0;  /* 0=sp3, 1=sp2, 2=sp */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            if (mol->bonds[bond_idx].type == BOND_DOUBLE) hyb = 1;
            if (mol->bonds[bond_idx].type == BOND_TRIPLE) hyb = 2;
            if (mol->bonds[bond_idx].type == BOND_AROMATIC) hyb = 1;
        }

        int type_idx = atom->element * 4 + hyb;
        if (type_idx < 128) type_counts[type_idx]++;
    }

    for (int i = 0; i < mol->num_bonds; i++) {
        switch (mol->bonds[i].type) {
            case BOND_SINGLE: bond_sum += 1; break;
            case BOND_DOUBLE: bond_sum += 2; break;
            case BOND_TRIPLE: bond_sum += 3; break;
            case BOND_AROMATIC: bond_sum += 1.5; break;
            default: bond_sum += 1; break;
        }
    }

    /* Information content of atom type distribution */
    int total_types = 0;
    double type_entropy = 0.0;
    for (int i = 0; i < 128; i++) {
        if (type_counts[i] > 0) {
            total_types += type_counts[i];
        }
    }

    if (total_types > 1) {
        for (int i = 0; i < 128; i++) {
            if (type_counts[i] > 0) {
                double p = (double)type_counts[i] / total_types;
                type_entropy -= p * log2(p);
            }
        }
    }

    value->d = log2(1.0 + bond_sum) + type_entropy * total_types;
    return CCHEM_OK;
}

/* Molecular flexibility index: rotatable bonds / total bonds */
static cchem_status_t desc_flexibility_index(const molecule_t* mol, descriptor_value_t* value) {
    if (mol->num_bonds == 0) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    int rotatable = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_SINGLE) continue;
        if (bond->in_ring) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        /* Skip terminal bonds */
        if (heavy_neighbors(mol, a1) <= 1 || heavy_neighbors(mol, a2) <= 1) continue;

        rotatable++;
    }

    value->d = (double)rotatable / mol->num_bonds;
    return CCHEM_OK;
}

/* Sp3 fraction */
static cchem_status_t desc_sp3_fraction(const molecule_t* mol, descriptor_value_t* value) {
    int sp3 = 0, total = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        total++;

        if (is_sp3_carbon(mol, atom)) sp3++;
    }

    value->d = total > 0 ? (double)sp3 / total : 0.0;
    return CCHEM_OK;
}

/* Stereocenter count (tetrahedral + double bond) */
static cchem_status_t desc_stereocenter_count(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;

    /* Tetrahedral stereocenters: sp3 C with 4 different substituents (simplified) */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (!is_sp3_carbon(mol, atom)) continue;
        if (heavy_neighbors(mol, atom) < 3) continue;  /* Need at least 3 heavy + H */

        /* Check for chirality marking */
        if (atom->chirality != CHIRALITY_NONE) {
            count++;
        }
    }

    /* E/Z stereocenters: double bonds with different substituents */
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_DOUBLE) continue;

        /* Check for stereochemistry specification */
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        /* C=C with at least 2 substituents on each */
        if (a1->element == ELEM_C && a2->element == ELEM_C) {
            if (a1->num_neighbors >= 2 && a2->num_neighbors >= 2) {
                /* Could be stereogenic - count it */
                count++;
            }
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Ring complexity: sum of ring sizes / num rings */
static cchem_status_t desc_ring_complexity(const molecule_t* mol, descriptor_value_t* value) {
    if (mol->num_rings == 0) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    int sum = 0;
    int hetero_rings = 0;
    int fused = 0;

    for (int i = 0; i < mol->num_rings; i++) {
        sum += mol->rings[i].size;

        /* Count heterocyclic rings */
        bool has_hetero = false;
        for (int j = 0; j < mol->rings[i].size; j++) {
            element_t e = mol->atoms[mol->rings[i].atoms[j]].element;
            if (e != ELEM_C) {
                has_hetero = true;
                break;
            }
        }
        if (has_hetero) hetero_rings++;

        /* Check for ring fusion (shares atoms with other rings) */
        for (int j = i + 1; j < mol->num_rings; j++) {
            for (int k = 0; k < mol->rings[i].size; k++) {
                for (int l = 0; l < mol->rings[j].size; l++) {
                    if (mol->rings[i].atoms[k] == mol->rings[j].atoms[l]) {
                        fused++;
                        goto next_ring;
                    }
                }
            }
            next_ring:;
        }
    }

    /* Complexity = avg size + hetero bonus + fusion bonus */
    value->d = (double)sum / mol->num_rings + 0.5 * hetero_rings + 0.3 * fused;
    return CCHEM_OK;
}

/* Information content: Shannon entropy of element distribution */
static cchem_status_t desc_element_entropy(const molecule_t* mol, descriptor_value_t* value) {
    int counts[128] = {0};
    int total = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int e = mol->atoms[i].element;
        if (e < 128) {
            counts[e]++;
            total++;
        }
    }

    if (total <= 1) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    double entropy = 0.0;
    for (int i = 0; i < 128; i++) {
        if (counts[i] > 0) {
            double p = (double)counts[i] / total;
            entropy -= p * log2(p);
        }
    }

    value->d = entropy;
    return CCHEM_OK;
}

/* Heteroatom diversity: number of distinct heteroatom types */
static cchem_status_t desc_heteroatom_diversity(const molecule_t* mol, descriptor_value_t* value) {
    bool seen[128] = {false};
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (e == ELEM_C || e == ELEM_H) continue;
        if (e < 128 && !seen[e]) {
            seen[e] = true;
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Drug-likeness Descriptors
 * ============================================================================ */

/* Rule of 5 violations (Lipinski) */
static cchem_status_t desc_ro5_violations(const molecule_t* mol, descriptor_value_t* value) {
    int violations = 0;

    /* Calculate MW */
    double mw = 0.0;
    int hbd = 0, hba = 0, rotatable = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Approximate MW */
        switch (atom->element) {
            case ELEM_H: mw += 1.008; break;
            case ELEM_C: mw += 12.011; break;
            case ELEM_N: mw += 14.007; break;
            case ELEM_O: mw += 15.999; break;
            case ELEM_F: mw += 18.998; break;
            case ELEM_Cl: mw += 35.453; break;
            case ELEM_Br: mw += 79.904; break;
            case ELEM_I: mw += 126.904; break;
            case ELEM_S: mw += 32.065; break;
            case ELEM_P: mw += 30.974; break;
            default: mw += 10.0; break;
        }

        int h = total_h(mol, atom);

        /* HBD: OH, NH */
        if ((atom->element == ELEM_O || atom->element == ELEM_N) && h > 0) {
            hbd++;
        }

        /* HBA: O, N */
        if (atom->element == ELEM_O || atom->element == ELEM_N) {
            hba++;
        }
    }

    /* Rotatable bonds */
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_SINGLE) continue;
        if (bond->in_ring) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (heavy_neighbors(mol, a1) > 1 && heavy_neighbors(mol, a2) > 1) {
            rotatable++;
        }
    }

    /* Check violations */
    if (mw > 500) violations++;
    if (hbd > 5) violations++;
    if (hba > 10) violations++;
    if (rotatable > 10) violations++;  /* Veber extension */

    value->i = violations;
    return CCHEM_OK;
}

/* Lead-likeness violations (MW<450, LogP<4.2, HBD<=4, HBA<=8, Rotatable<=10) */
static cchem_status_t desc_lead_violations(const molecule_t* mol, descriptor_value_t* value) {
    int violations = 0;

    double mw = 0.0;
    int hbd = 0, hba = 0, rotatable = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        switch (atom->element) {
            case ELEM_H: mw += 1.008; break;
            case ELEM_C: mw += 12.011; break;
            case ELEM_N: mw += 14.007; break;
            case ELEM_O: mw += 15.999; break;
            case ELEM_F: mw += 18.998; break;
            case ELEM_Cl: mw += 35.453; break;
            case ELEM_Br: mw += 79.904; break;
            case ELEM_I: mw += 126.904; break;
            case ELEM_S: mw += 32.065; break;
            case ELEM_P: mw += 30.974; break;
            default: mw += 10.0; break;
        }

        int h = total_h(mol, atom);
        if ((atom->element == ELEM_O || atom->element == ELEM_N) && h > 0) hbd++;
        if (atom->element == ELEM_O || atom->element == ELEM_N) hba++;
    }

    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_SINGLE || bond->in_ring) continue;
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];
        if (heavy_neighbors(mol, a1) > 1 && heavy_neighbors(mol, a2) > 1) rotatable++;
    }

    if (mw > 450) violations++;
    if (hbd > 4) violations++;
    if (hba > 8) violations++;
    if (rotatable > 10) violations++;

    value->i = violations;
    return CCHEM_OK;
}

/* CNS MPO score (simplified: MW, PSA, HBD, basic N) */
static cchem_status_t desc_cns_mpo(const molecule_t* mol, descriptor_value_t* value) {
    double score = 0.0;

    double mw = 0.0;
    double psa = 0.0;
    int hbd = 0;
    int basic_n = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* MW */
        switch (atom->element) {
            case ELEM_H: mw += 1.008; break;
            case ELEM_C: mw += 12.011; break;
            case ELEM_N: mw += 14.007; break;
            case ELEM_O: mw += 15.999; break;
            case ELEM_F: mw += 18.998; break;
            case ELEM_Cl: mw += 35.453; break;
            case ELEM_Br: mw += 79.904; break;
            case ELEM_I: mw += 126.904; break;
            case ELEM_S: mw += 32.065; break;
            case ELEM_P: mw += 30.974; break;
            default: mw += 10.0; break;
        }

        int h = total_h(mol, atom);

        /* Simplified PSA */
        if (atom->element == ELEM_N) {
            psa += 26.0;  /* N contribution */
            if (h > 0) hbd++;
        }
        if (atom->element == ELEM_O) {
            psa += 20.0;  /* O contribution */
            if (h > 0) hbd++;
        }

        /* Basic N */
        if (atom->element == ELEM_N && !atom->aromatic && h > 0) {
            bool is_amide = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* nbr = &mol->atoms[atom->neighbors[j]];
                if (nbr->element == ELEM_C && has_double_to(mol, nbr, ELEM_O)) {
                    is_amide = true;
                    break;
                }
            }
            if (!is_amide) basic_n++;
        }
    }

    /* Score components (0-1 each, sum to max ~6) */
    /* MW: optimal < 360 */
    score += mw < 360 ? 1.0 : (mw < 500 ? 0.5 : 0.0);

    /* PSA: optimal 40-90 */
    if (psa >= 40 && psa <= 90) score += 1.0;
    else if (psa < 40 || (psa > 90 && psa <= 120)) score += 0.5;

    /* HBD: optimal <=1 */
    score += hbd <= 1 ? 1.0 : (hbd <= 2 ? 0.5 : 0.0);

    /* Basic N: optimal <=1 */
    score += basic_n <= 1 ? 1.0 : (basic_n <= 2 ? 0.5 : 0.0);

    value->d = score;
    return CCHEM_OK;
}

/* Heavy atom count */
static cchem_status_t desc_heavy_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* Ring count */
static cchem_status_t desc_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    value->i = mol->num_rings;
    return CCHEM_OK;
}

/* Heteroatom count */
static cchem_status_t desc_heteroatom_count(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (e != ELEM_C && e != ELEM_H) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Additional Pharmacophore Descriptors
 * ============================================================================ */

/* Donor-Acceptor balance: (HBD - HBA) / (HBD + HBA + 1) */
static cchem_status_t desc_donor_acceptor_balance(const molecule_t* mol, descriptor_value_t* value) {
    int hbd = 0, hba = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        int h = total_h(mol, atom);

        if ((atom->element == ELEM_O || atom->element == ELEM_N) && h > 0) hbd++;
        if (atom->element == ELEM_O || atom->element == ELEM_N) hba++;
    }

    value->d = (double)(hbd - hba) / (hbd + hba + 1);
    return CCHEM_OK;
}

/* Polar surface area fraction */
static cchem_status_t desc_polar_surface_fraction(const molecule_t* mol, descriptor_value_t* value) {
    double polar = 0.0;
    double total = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        double r = 0.0;
        switch (atom->element) {
            case ELEM_C: r = 1.70; break;
            case ELEM_N: r = 1.55; polar += 26.0; break;
            case ELEM_O: r = 1.52; polar += 20.0; break;
            case ELEM_S: r = 1.80; polar += 25.0; break;
            case ELEM_P: r = 1.80; polar += 25.0; break;
            case ELEM_F: r = 1.47; break;
            case ELEM_Cl: r = 1.75; break;
            case ELEM_Br: r = 1.85; break;
            case ELEM_I: r = 1.98; break;
            default: r = 1.70; break;
        }
        total += 4.0 * 3.14159 * r * r;
    }

    value->d = total > 0 ? polar / total : 0.0;
    return CCHEM_OK;
}

/* Aromatic to aliphatic ring ratio */
static cchem_status_t desc_aromatic_aliphatic_ratio(const molecule_t* mol, descriptor_value_t* value) {
    int aromatic = 0, aliphatic = 0;

    for (int i = 0; i < mol->num_rings; i++) {
        if (ring_is_aromatic(mol, &mol->rings[i])) {
            aromatic++;
        } else {
            aliphatic++;
        }
    }

    value->d = (aliphatic > 0) ? (double)aromatic / aliphatic : (aromatic > 0 ? (double)aromatic : 0.0);
    return CCHEM_OK;
}

/* Functional group diversity: count of distinct functional group types */
static cchem_status_t desc_functional_diversity(const molecule_t* mol, descriptor_value_t* value) {
    int types = 0;

    /* Check for various functional groups */
    bool has_alcohol = false, has_amine = false, has_carboxyl = false;
    bool has_amide = false, has_ester = false, has_ketone = false;
    bool has_aldehyde = false, has_ether = false, has_nitro = false;
    bool has_sulfone = false, has_halide = false;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        int h = total_h(mol, atom);
        int heavy = heavy_neighbors(mol, atom);

        /* Alcohol: C-OH */
        if (atom->element == ELEM_O && h == 1 && heavy == 1) {
            for (int j = 0; j < atom->num_neighbors; j++) {
                if (mol->atoms[atom->neighbors[j]].element == ELEM_C) {
                    const atom_t* c = &mol->atoms[atom->neighbors[j]];
                    if (!has_double_to(mol, c, ELEM_O)) {
                        has_alcohol = true;
                    }
                }
            }
        }

        /* Amine: C-NH2, C-NHR, C-NR2 */
        if (atom->element == ELEM_N && !atom->aromatic) {
            bool is_amide = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                const atom_t* nbr = &mol->atoms[atom->neighbors[j]];
                if (nbr->element == ELEM_C && has_double_to(mol, nbr, ELEM_O)) {
                    is_amide = true;
                    has_amide = true;
                    break;
                }
            }
            if (!is_amide) has_amine = true;
        }

        /* Carboxyl, ester, aldehyde, ketone */
        if (atom->element == ELEM_C && has_double_to(mol, atom, ELEM_O)) {
            if (has_single_to(mol, atom, ELEM_O)) {
                /* Could be carboxyl or ester */
                for (int j = 0; j < atom->num_neighbors; j++) {
                    const atom_t* o = &mol->atoms[atom->neighbors[j]];
                    if (o->element == ELEM_O) {
                        int bond_idx = atom->neighbor_bonds[j];
                        if (mol->bonds[bond_idx].type == BOND_SINGLE) {
                            if (total_h(mol, o) > 0) has_carboxyl = true;
                            else has_ester = true;
                        }
                    }
                }
            } else if (has_single_to(mol, atom, ELEM_N)) {
                has_amide = true;
            } else if (h >= 1) {
                has_aldehyde = true;
            } else if (heavy >= 2) {
                has_ketone = true;
            }
        }

        /* Ether */
        if (atom->element == ELEM_O && h == 0 && heavy == 2) {
            bool is_ester = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                if (has_double_to(mol, &mol->atoms[atom->neighbors[j]], ELEM_O)) {
                    is_ester = true;
                    break;
                }
            }
            if (!is_ester) has_ether = true;
        }

        /* Nitro */
        if (atom->element == ELEM_N && atom->charge == 1) {
            int o_count = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                if (mol->atoms[atom->neighbors[j]].element == ELEM_O) o_count++;
            }
            if (o_count >= 2) has_nitro = true;
        }

        /* Sulfone/sulfoxide */
        if (atom->element == ELEM_S && has_double_to(mol, atom, ELEM_O)) {
            has_sulfone = true;
        }

        /* Halide */
        if (atom->element == ELEM_F || atom->element == ELEM_Cl ||
            atom->element == ELEM_Br || atom->element == ELEM_I) {
            has_halide = true;
        }
    }

    if (has_alcohol) types++;
    if (has_amine) types++;
    if (has_carboxyl) types++;
    if (has_amide) types++;
    if (has_ester) types++;
    if (has_ketone) types++;
    if (has_aldehyde) types++;
    if (has_ether) types++;
    if (has_nitro) types++;
    if (has_sulfone) types++;
    if (has_halide) types++;

    value->i = types;
    return CCHEM_OK;
}

/* Chiral center count (marked stereocenters) */
static cchem_status_t desc_chiral_centers(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].chirality != CHIRALITY_NONE) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* Bridge atom count (atoms in multiple rings) */
static cchem_status_t desc_bridge_atoms(const molecule_t* mol, descriptor_value_t* value) {
    if (mol->num_rings < 2) {
        value->i = 0;
        return CCHEM_OK;
    }

    int ring_membership[1024] = {0};  /* Count rings each atom belongs to */

    for (int i = 0; i < mol->num_rings; i++) {
        for (int j = 0; j < mol->rings[i].size; j++) {
            int atom_idx = mol->rings[i].atoms[j];
            if (atom_idx < 1024) {
                ring_membership[atom_idx]++;
            }
        }
    }

    int count = 0;
    for (int i = 0; i < mol->num_atoms && i < 1024; i++) {
        if (ring_membership[i] >= 2) count++;
    }

    value->i = count;
    return CCHEM_OK;
}

/* Spiro center count (atoms connecting two rings via single bonds only) */
static cchem_status_t desc_spiro_centers(const molecule_t* mol, descriptor_value_t* value) {
    if (mol->num_rings < 2) {
        value->i = 0;
        return CCHEM_OK;
    }

    int count = 0;

    /* Find atoms that are the only shared atom between two rings */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count == 0) continue;

        /* Count how many rings this atom is in */
        int in_rings = 0;
        for (int r = 0; r < mol->num_rings; r++) {
            for (int j = 0; j < mol->rings[r].size; j++) {
                if (mol->rings[r].atoms[j] == i) {
                    in_rings++;
                    break;
                }
            }
        }

        if (in_rings < 2) continue;

        /* Check if this is the only shared atom (spiro) by checking neighbors */
        /* A spiro center's ring bonds don't share any other atoms */
        bool is_spiro = true;
        for (int j = 0; j < atom->num_neighbors; j++) {
            const atom_t* nbr = &mol->atoms[atom->neighbors[j]];
            if (nbr->ring_count > 0) {
                /* Check if this neighbor is also in multiple rings */
                int nbr_rings = 0;
                for (int r = 0; r < mol->num_rings; r++) {
                    for (int k = 0; k < mol->rings[r].size; k++) {
                        if (mol->rings[r].atoms[k] == atom->neighbors[j]) {
                            nbr_rings++;
                            break;
                        }
                    }
                }
                if (nbr_rings >= 2) {
                    is_spiro = false;
                    break;
                }
            }
        }

        if (is_spiro && in_rings >= 2) count++;
    }

    value->i = count;
    return CCHEM_OK;
}

/* Macrocycle presence (ring size >= 12) */
static cchem_status_t desc_macrocycle_count(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size >= 12) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration Macro
 * ============================================================================ */

#define REGISTER_PHARM_DESC_INT(name_str, desc_str, func) \
    do { \
        descriptor_def_t def = {0}; \
        strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
        strncpy(def.description, desc_str, sizeof(def.description) - 1); \
        def.category = DESC_CATEGORY_CUSTOM; \
        def.value_type = DESC_VALUE_INT; \
        def.compute = func; \
        descriptor_register(&def); \
    } while (0)

#define REGISTER_PHARM_DESC_DOUBLE(name_str, desc_str, func) \
    do { \
        descriptor_def_t def = {0}; \
        strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
        strncpy(def.description, desc_str, sizeof(def.description) - 1); \
        def.category = DESC_CATEGORY_CUSTOM; \
        def.value_type = DESC_VALUE_DOUBLE; \
        def.compute = func; \
        descriptor_register(&def); \
    } while (0)

/* ============================================================================
 * Public Registration Function
 * ============================================================================ */

void descriptors_register_pharmacophore(void) {
    /* Pharmacophore points */
    REGISTER_PHARM_DESC_INT("LipophilicCenters", "Lipophilic centers (sp3 C with only C/H neighbors)", desc_lipophilic_centers);
    REGISTER_PHARM_DESC_INT("BasicNitrogenCount", "Basic nitrogens (protonatable at pH 7)", desc_basic_nitrogen_count);
    REGISTER_PHARM_DESC_INT("AcidicOxygenCount", "Acidic oxygens (carboxylic, sulfonic, phosphoric)", desc_acidic_oxygen_count);
    REGISTER_PHARM_DESC_INT("HalogenBondDonors", "Halogen bond donors (C-Cl, C-Br, C-I)", desc_halogen_bond_donor);
    REGISTER_PHARM_DESC_INT("AromaticRings", "Aromatic rings (π-stacking capable)", desc_aromatic_ring_count);
    REGISTER_PHARM_DESC_INT("PositiveChargeCenters", "Formal positive charge centers", desc_positive_charge_centers);
    REGISTER_PHARM_DESC_INT("NegativeChargeCenters", "Formal negative charge centers", desc_negative_charge_centers);
    REGISTER_PHARM_DESC_INT("TotalPharmacophorePoints", "Total pharmacophore points", desc_total_pharmacophore_points);
    REGISTER_PHARM_DESC_DOUBLE("PharmacophoreDensity", "Pharmacophore points per heavy atom", desc_pharmacophore_density);

    /* Molecular complexity */
    REGISTER_PHARM_DESC_DOUBLE("BertzComplexity", "Bertz molecular complexity index", desc_bertz_complexity);
    REGISTER_PHARM_DESC_DOUBLE("FlexibilityIndex", "Molecular flexibility (rotatable/total bonds)", desc_flexibility_index);
    REGISTER_PHARM_DESC_DOUBLE("Sp3Fraction", "Fraction of sp3 carbons", desc_sp3_fraction);
    REGISTER_PHARM_DESC_INT("StereocenterCount", "Potential stereocenters (tetrahedral + E/Z)", desc_stereocenter_count);
    REGISTER_PHARM_DESC_DOUBLE("RingComplexity", "Ring system complexity score", desc_ring_complexity);
    REGISTER_PHARM_DESC_DOUBLE("ElementEntropy", "Shannon entropy of element distribution", desc_element_entropy);
    REGISTER_PHARM_DESC_INT("HeteroatomDiversity", "Number of distinct heteroatom types", desc_heteroatom_diversity);

    /* Drug-likeness */
    REGISTER_PHARM_DESC_INT("Ro5Violations", "Lipinski Rule of 5 violations", desc_ro5_violations);
    REGISTER_PHARM_DESC_INT("LeadViolations", "Lead-likeness violations", desc_lead_violations);
    REGISTER_PHARM_DESC_DOUBLE("CNS_MPO", "CNS multiparameter optimization score", desc_cns_mpo);

    /* Additional features */
    REGISTER_PHARM_DESC_INT("HeavyAtomCount2", "Heavy atom count (non-H)", desc_heavy_atom_count);
    REGISTER_PHARM_DESC_INT("RingCount2", "Total ring count", desc_ring_count);
    REGISTER_PHARM_DESC_INT("HeteroatomCount", "Total heteroatom count", desc_heteroatom_count);
    REGISTER_PHARM_DESC_DOUBLE("DonorAcceptorBalance", "HBD-HBA balance ratio", desc_donor_acceptor_balance);
    REGISTER_PHARM_DESC_DOUBLE("PolarSurfaceFraction", "Polar/total surface area fraction", desc_polar_surface_fraction);
    REGISTER_PHARM_DESC_DOUBLE("AromaticAliphaticRatio", "Aromatic/aliphatic ring ratio", desc_aromatic_aliphatic_ratio);
    REGISTER_PHARM_DESC_INT("FunctionalDiversity", "Distinct functional group types", desc_functional_diversity);
    REGISTER_PHARM_DESC_INT("ChiralCenters", "Marked chiral centers", desc_chiral_centers);
    REGISTER_PHARM_DESC_INT("BridgeAtoms", "Atoms in multiple rings (fused)", desc_bridge_atoms);
    REGISTER_PHARM_DESC_INT("SpiroCenters", "Spiro junction atoms", desc_spiro_centers);
    REGISTER_PHARM_DESC_INT("MacrocycleCount", "Macrocycles (ring size >= 12)", desc_macrocycle_count);
}
