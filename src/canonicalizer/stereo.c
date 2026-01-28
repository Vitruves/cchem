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

#if 0  /* Debug */
    fprintf(stderr, "  permutation: ");
    for (int i = 0; i < n; i++) fprintf(stderr, "%d ", perm[i]);
    fprintf(stderr, "-> parity=%d\n", parity);
#endif

    free(perm);

    chirality_t result = orig;
    /* Odd parity inverts chirality */
    if (parity % 2 == 1) {
        if (orig == CHIRALITY_CW) result = CHIRALITY_CCW;
        if (orig == CHIRALITY_CCW) result = CHIRALITY_CW;
    }

    return result;
}

chirality_t stereo_get_canonical_chirality(const molecule_t* mol, int atom_idx,
                                           int from_atom, const bool* visited) {
    const atom_t* atom = molecule_get_atom_const(mol, atom_idx);
    if (!atom || atom->chirality == CHIRALITY_NONE) return CHIRALITY_NONE;

    /*
     * Get original neighbor order from parsing.
     * For bracket atoms like [C@H], the order in SMILES is:
     * 1. From-atom (implicit, the atom we came from in the chain) - if present
     * 2. H (explicitly written in bracket notation)
     * 3. Branch atoms / ring closures
     * 4. Main chain continuation
     *
     * stereo_neighbors stores: [from, branch1, branch2, ...] if there was a from-atom
     * For atoms at the START of SMILES (no from-atom), stereo_neighbors just has
     * the neighbors in SMILES order.
     */
    int orig_neighbors[4];
    int num_orig = 0;

    /* Check if this atom had a from-atom in the original SMILES */
    bool had_from_atom = false;
    if (atom->num_stereo_neighbors > 0) {
        /* Use from_atom parameter: if from_atom >= 0, there's a from-atom in canon order,
         * so there was likely one in original too. If from_atom < 0, this is start atom. */
        had_from_atom = (from_atom >= 0) || (atom->stereo_neighbors[0] < atom->index);
    }

    if (atom->implicit_h_count > 0 && atom->num_stereo_neighbors > 0) {
        if (had_from_atom) {
            /* For [C@H] notation with from-atom: from-atom first, then H, then rest */
            orig_neighbors[num_orig++] = atom->stereo_neighbors[0];  /* from-atom */
            orig_neighbors[num_orig++] = -1;  /* H at position 1 */
            for (int i = 1; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
                orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
            }
        } else {
            /* For [C@H] at START of SMILES (no from-atom): H first, then neighbors */
            orig_neighbors[num_orig++] = -1;  /* H first in bracket notation */
            for (int i = 0; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
                orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
            }
        }
    } else {
        /* No implicit H - just use stereo_neighbors as-is */
        for (int i = 0; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
            orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
        }
    }

    /*
     * Build canonical neighbor order based on SMILES output order.
     * IMPORTANT: In SMILES output, the order is:
     * 1. From-atom (implicit, we came from there in DFS)
     * 2. H (if implicit, written explicitly in bracket)
     * 3. Ring closures (to VISITED atoms) - these come before branches in SMILES syntax
     * 4. Branches/chain (to UNVISITED atoms), sorted by rank
     *
     * This is critical for correct stereo: ring closures appear in SMILES
     * immediately after the atom symbol, before any branch parentheses.
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

    /* Collect neighbors excluding from_atom, split into ring closures and unvisited */
    int ring_closures[4];
    int num_ring = 0;
    int unvisited_neighbors[4];
    int num_unvisited = 0;

    for (int i = 0; i < atom->num_neighbors; i++) {
        int neighbor = atom->neighbors[i];
        if (neighbor == from_atom) continue;

        if (visited && visited[neighbor]) {
            /* Already visited = ring closure */
            ring_closures[num_ring++] = neighbor;
        } else {
            /* Not visited = will be branch or chain */
            unvisited_neighbors[num_unvisited++] = neighbor;
        }
    }

    /* Sort ring closures by rank (they appear in rank order in the output) */
    for (int i = 0; i < num_ring - 1; i++) {
        for (int j = i + 1; j < num_ring; j++) {
            if (mol->atoms[ring_closures[i]].canon_rank > mol->atoms[ring_closures[j]].canon_rank) {
                int tmp = ring_closures[i];
                ring_closures[i] = ring_closures[j];
                ring_closures[j] = tmp;
            }
        }
    }

    /* Sort unvisited neighbors by rank */
    for (int i = 0; i < num_unvisited - 1; i++) {
        for (int j = i + 1; j < num_unvisited; j++) {
            if (mol->atoms[unvisited_neighbors[i]].canon_rank > mol->atoms[unvisited_neighbors[j]].canon_rank) {
                int tmp = unvisited_neighbors[i];
                unvisited_neighbors[i] = unvisited_neighbors[j];
                unvisited_neighbors[j] = tmp;
            }
        }
    }

    /* Add ring closures first (they appear before branches in SMILES) */
    for (int i = 0; i < num_ring && num_canon < 4; i++) {
        canon_neighbors[num_canon++] = ring_closures[i];
    }

    /* Then add unvisited neighbors (branches and chain) */
    for (int i = 0; i < num_unvisited && num_canon < 4; i++) {
        canon_neighbors[num_canon++] = unvisited_neighbors[i];
    }

#if 0  /* Debug output */
    fprintf(stderr, "stereo_get_canonical_chirality: atom %d, from_atom=%d, orig_chirality=%d\n",
            atom_idx, from_atom, atom->chirality);
    fprintf(stderr, "  had_from_atom=%d, stereo_neighbors[0]=%d, atom->index=%d\n",
            had_from_atom, atom->num_stereo_neighbors > 0 ? atom->stereo_neighbors[0] : -999, atom->index);
    fprintf(stderr, "  stereo_neighbors: ");
    for (int i = 0; i < atom->num_stereo_neighbors; i++) fprintf(stderr, "%d(r%d) ", atom->stereo_neighbors[i], mol->atoms[atom->stereo_neighbors[i]].canon_rank);
    fprintf(stderr, "\n  orig_neighbors[%d]: ", num_orig);
    for (int i = 0; i < num_orig; i++) fprintf(stderr, "%d ", orig_neighbors[i]);
    fprintf(stderr, "\n  canon_neighbors[%d]: ", num_canon);
    for (int i = 0; i < num_canon; i++) fprintf(stderr, "%d ", canon_neighbors[i]);
    fprintf(stderr, "\n");
#endif

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
