/**
 * @file stereo.c
 * @brief Stereochemistry detection and canonicalization
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/stereo.h"

stereo_info_t* stereo_info_create(void) {
    stereo_info_t* info = (stereo_info_t*)calloc(1, sizeof(stereo_info_t));
    if (!info) return NULL;

    info->centers_capacity = 32;
    info->centers = (stereo_center_t*)calloc(info->centers_capacity, sizeof(stereo_center_t));

    info->double_bonds_capacity = 32;
    info->double_bonds = (stereo_double_bond_t*)calloc(info->double_bonds_capacity,
                                                        sizeof(stereo_double_bond_t));

    if (!info->centers || !info->double_bonds) {
        stereo_info_free(info);
        return NULL;
    }

    return info;
}

void stereo_info_free(stereo_info_t* info) {
    if (!info) return;
    if (info->centers) free(info->centers);
    if (info->double_bonds) free(info->double_bonds);
    free(info);
}

/* Check if atom could be a tetrahedral stereocenter */
bool stereo_is_potential_center(const molecule_t* mol, int atom_idx) {
    if (!mol) return false;

    const atom_t* atom = molecule_get_atom_const(mol, atom_idx);
    if (!atom) return false;

    /* Must have 4 distinct neighbors (including implicit H) */
    int total_neighbors = atom->num_neighbors + atom->implicit_h_count;
    if (total_neighbors != 4) return false;

    /* Check for 4 distinct substituents */
    if (atom->num_neighbors < 3) return false;  /* Need at least 3 explicit neighbors */

    /* Typical stereocenter elements */
    if (atom->element != ELEM_C && atom->element != ELEM_N &&
        atom->element != ELEM_P && atom->element != ELEM_S &&
        atom->element != ELEM_Si) {
        return false;
    }

    /* For now, just check number of neighbors */
    return true;
}

/* Check if bond could be E/Z stereocenter */
bool stereo_is_potential_double_bond(const molecule_t* mol, int bond_idx) {
    if (!mol || bond_idx < 0 || bond_idx >= mol->num_bonds) return false;

    const bond_t* bond = &mol->bonds[bond_idx];

    /* Must be double bond */
    if (bond->type != BOND_DOUBLE) return false;

    /* Can't be in ring (ring double bonds are usually fixed) */
    if (bond->in_ring) return false;

    /* Check both atoms have at least one other substituent */
    const atom_t* atom1 = molecule_get_atom_const(mol, bond->atom1);
    const atom_t* atom2 = molecule_get_atom_const(mol, bond->atom2);

    if (!atom1 || !atom2) return false;

    /* Each atom needs at least one other neighbor */
    if (atom1->num_neighbors < 2 || atom2->num_neighbors < 2) return false;

    return true;
}

cchem_status_t stereo_detect_tetrahedral(molecule_t* mol, stereo_info_t* info) {
    if (!mol || !info) return CCHEM_ERROR_INVALID_INPUT;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (!stereo_is_potential_center(mol, i)) continue;

        atom_t* atom = molecule_get_atom(mol, i);

        /* Create stereocenter entry */
        if (info->num_centers >= info->centers_capacity) {
            int new_cap = info->centers_capacity * 2;
            stereo_center_t* new_centers = (stereo_center_t*)realloc(info->centers,
                                            new_cap * sizeof(stereo_center_t));
            if (!new_centers) return CCHEM_ERROR_MEMORY;
            info->centers = new_centers;
            info->centers_capacity = new_cap;
        }

        stereo_center_t* center = &info->centers[info->num_centers];
        memset(center, 0, sizeof(stereo_center_t));

        center->type = STEREO_CENTER_TETRAHEDRAL;
        center->center_atom = i;
        center->chirality = atom->chirality;
        center->is_specified = (atom->chirality != CHIRALITY_NONE);
        center->can_be_stereo = true;

        /* Store neighbors in order */
        center->num_atoms = 0;
        for (int j = 0; j < atom->num_stereo_neighbors && center->num_atoms < 4; j++) {
            center->atoms[center->num_atoms++] = atom->stereo_neighbors[j];
        }

        /* Add implicit hydrogen as "virtual" neighbor if needed */
        if (center->num_atoms < 4 && atom->implicit_h_count > 0) {
            center->atoms[center->num_atoms++] = -1;  /* -1 indicates implicit H */
        }

        info->num_centers++;
    }

    return CCHEM_OK;
}

cchem_status_t stereo_detect_double_bonds(molecule_t* mol, stereo_info_t* info) {
    if (!mol || !info) return CCHEM_ERROR_INVALID_INPUT;

    for (int i = 0; i < mol->num_bonds; i++) {
        if (!stereo_is_potential_double_bond(mol, i)) continue;

        bond_t* bond = &mol->bonds[i];

        /* Create double bond stereo entry */
        if (info->num_double_bonds >= info->double_bonds_capacity) {
            int new_cap = info->double_bonds_capacity * 2;
            stereo_double_bond_t* new_db = (stereo_double_bond_t*)realloc(info->double_bonds,
                                            new_cap * sizeof(stereo_double_bond_t));
            if (!new_db) return CCHEM_ERROR_MEMORY;
            info->double_bonds = new_db;
            info->double_bonds_capacity = new_cap;
        }

        stereo_double_bond_t* db = &info->double_bonds[info->num_double_bonds];
        memset(db, 0, sizeof(stereo_double_bond_t));

        db->bond_idx = i;
        db->atom1 = bond->atom1;
        db->atom2 = bond->atom2;
        db->stereo = bond->stereo;
        db->is_specified = (bond->stereo != STEREO_NONE || bond->stereo_type != BOND_NONE);

        /* Find reference atoms for E/Z determination */
        atom_t* a1 = molecule_get_atom(mol, bond->atom1);
        atom_t* a2 = molecule_get_atom(mol, bond->atom2);

        /* Pick first non-double-bond neighbor as reference */
        db->ref_atom1 = -1;
        db->ref_atom2 = -1;

        for (int j = 0; j < a1->num_neighbors; j++) {
            if (a1->neighbors[j] != bond->atom2) {
                db->ref_atom1 = a1->neighbors[j];
                break;
            }
        }

        for (int j = 0; j < a2->num_neighbors; j++) {
            if (a2->neighbors[j] != bond->atom1) {
                db->ref_atom2 = a2->neighbors[j];
                break;
            }
        }

        info->num_double_bonds++;
    }

    return CCHEM_OK;
}

cchem_status_t stereo_detect_centers(molecule_t* mol, stereo_info_t* info) {
    if (!mol || !info) return CCHEM_ERROR_INVALID_INPUT;

    info->num_centers = 0;
    info->num_double_bonds = 0;

    cchem_status_t status = stereo_detect_tetrahedral(mol, info);
    if (status != CCHEM_OK) return status;

    status = stereo_detect_double_bonds(mol, info);
    if (status != CCHEM_OK) return status;

    return CCHEM_OK;
}

/* Determine chirality based on neighbor ordering */
chirality_t stereo_calc_chirality(const molecule_t* mol, int center_atom,
                                  const int* neighbors, int num_neighbors) {
    (void)mol;
    (void)center_atom;
    (void)neighbors;
    (void)num_neighbors;

    /* For proper implementation, would need to compute the parity of the
     * permutation from input order to canonical order */
    return CHIRALITY_CW;  /* Placeholder */
}

void stereo_invert_chirality(atom_t* atom) {
    if (!atom) return;

    switch (atom->chirality) {
        case CHIRALITY_CW:
            atom->chirality = CHIRALITY_CCW;
            break;
        case CHIRALITY_CCW:
            atom->chirality = CHIRALITY_CW;
            break;
        case CHIRALITY_TH1:
            atom->chirality = CHIRALITY_TH2;
            break;
        case CHIRALITY_TH2:
            atom->chirality = CHIRALITY_TH1;
            break;
        default:
            break;
    }
}

/* Permutation parity calculation */
static int permutation_parity(const int* perm, int n) {
    int* visited = (int*)calloc(n, sizeof(int));
    if (!visited) return 0;

    int parity = 0;

    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            int cycle_len = 0;
            int j = i;

            while (!visited[j]) {
                visited[j] = 1;
                j = perm[j];
                cycle_len++;
            }

            /* Each cycle of length k contributes k-1 to parity */
            parity += (cycle_len - 1);
        }
    }

    free(visited);
    return parity % 2;
}

chirality_t stereo_permute_chirality(chirality_t orig, const int* old_order,
                                     const int* new_order, int n) {
    if (orig == CHIRALITY_NONE) return CHIRALITY_NONE;

    /* Calculate permutation from old to new order */
    int* perm = (int*)calloc(n, sizeof(int));
    if (!perm) return orig;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (old_order[i] == new_order[j]) {
                perm[i] = j;
                break;
            }
        }
    }

    int parity = permutation_parity(perm, n);
    free(perm);

    /* Odd parity inverts chirality */
    if (parity % 2 == 1) {
        if (orig == CHIRALITY_CW) return CHIRALITY_CCW;
        if (orig == CHIRALITY_CCW) return CHIRALITY_CW;
    }

    return orig;
}

chirality_t stereo_get_canonical_chirality(const molecule_t* mol, int atom_idx,
                                           int from_atom, const int* output_order) {
    const atom_t* atom = molecule_get_atom_const(mol, atom_idx);
    if (!atom || atom->chirality == CHIRALITY_NONE) return CHIRALITY_NONE;

    (void)output_order;  /* Unused currently */

    /*
     * Get original neighbor order from parsing.
     * For bracket atoms like [C@H], the order in SMILES is:
     * 1. From-atom (implicit, the atom we came from in the chain)
     * 2. H (explicitly written in bracket notation)
     * 3. Branch atoms (in parentheses)
     * 4. Main chain continuation
     *
     * stereo_neighbors stores: [from, branch1, branch2, ...]
     * We need to insert H at position 1 (after from-atom) for correct stereo.
     */
    int orig_neighbors[4];
    int num_orig = 0;

    if (atom->implicit_h_count > 0 && atom->num_stereo_neighbors > 0) {
        /* For [C@H] notation: from-atom first, then H, then rest */
        orig_neighbors[num_orig++] = atom->stereo_neighbors[0];  /* from-atom */
        orig_neighbors[num_orig++] = -1;  /* H at position 1 */
        for (int i = 1; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
            orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
        }
    } else {
        /* No implicit H - just use stereo_neighbors as-is */
        for (int i = 0; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
            orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
        }
    }

    /*
     * Build canonical neighbor order based on SMILES output order.
     * In SMILES output:
     * 1. From-atom is first (implicit, we came from there in DFS)
     * 2. For [C@H], H is written next (explicit in bracket)
     * 3. Other neighbors follow, sorted by canonical rank (lower rank = branches first)
     * 4. Highest rank neighbor is main chain (written last, no parentheses)
     */
    int canon_neighbors[4];
    int num_canon = 0;

    /* From-atom first */
    if (from_atom >= 0) {
        canon_neighbors[num_canon++] = from_atom;
    }

    /* H second (for atoms with implicit H) */
    if (atom->implicit_h_count > 0) {
        canon_neighbors[num_canon++] = -1;
    }

    /* Collect other neighbors (excluding from_atom) and sort by canonical rank */
    int other_neighbors[4];
    int num_other = 0;
    for (int i = 0; i < atom->num_neighbors && num_other < 4; i++) {
        int neighbor = atom->neighbors[i];
        if (neighbor != from_atom) {
            other_neighbors[num_other++] = neighbor;
        }
    }

    /* Sort other neighbors by canonical rank (ascending) */
    for (int i = 0; i < num_other - 1; i++) {
        for (int j = i + 1; j < num_other; j++) {
            int rank_i = mol->atoms[other_neighbors[i]].canon_rank;
            int rank_j = mol->atoms[other_neighbors[j]].canon_rank;
            if (rank_i > rank_j) {
                int tmp = other_neighbors[i];
                other_neighbors[i] = other_neighbors[j];
                other_neighbors[j] = tmp;
            }
        }
    }

    /* Add sorted other neighbors to canonical order */
    for (int i = 0; i < num_other && num_canon < 4; i++) {
        canon_neighbors[num_canon++] = other_neighbors[i];
    }

#ifdef STEREO_DEBUG
    fprintf(stderr, "stereo_get_canonical_chirality: atom %d, from_atom=%d, orig_chirality=%d\n",
            atom_idx, from_atom, atom->chirality);
    fprintf(stderr, "  orig_neighbors: ");
    for (int i = 0; i < num_orig; i++) fprintf(stderr, "%d ", orig_neighbors[i]);
    fprintf(stderr, "\n  canon_neighbors: ");
    for (int i = 0; i < num_canon; i++) fprintf(stderr, "%d ", canon_neighbors[i]);
    fprintf(stderr, "\n");
#endif

    return stereo_permute_chirality(atom->chirality, orig_neighbors, canon_neighbors, num_canon);
}

cchem_status_t stereo_canonicalize(molecule_t* mol, stereo_info_t* info,
                                   const int* canon_order) {
    if (!mol || !info) return CCHEM_ERROR_INVALID_INPUT;

    /*
     * Note: We don't pre-calculate canonical chirality here because
     * the correct chirality depends on the DFS traversal order (from_atom),
     * which isn't known until SMILES writing time.
     *
     * The actual chirality calculation happens in smiles_write_atom()
     * via stereo_get_canonical_chirality() with the actual from_atom.
     */
    (void)canon_order;  /* Unused - chirality computed during SMILES writing */

    return CCHEM_OK;
}
