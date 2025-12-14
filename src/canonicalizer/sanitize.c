/**
 * @file sanitize.c
 * @brief Molecular structure sanitization and normalization
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "cchem/canonicalizer/sanitize.h"
#include "cchem/canonicalizer/canon.h"
#include "cchem/canonicalizer/parser.h"
#include "cchem/canonicalizer/smiles_writer.h"
#include "cchem/canonicalizer/ring_finder.h"

/* Default sanitization options */
const sanitize_options_t SANITIZE_OPTIONS_DEFAULT = {
    .flags = SANITIZE_COMPLETE,
    .keep_all_organic = false,
    .min_fragment_atoms = 3,
    .remove_water = true,
    .neutralize_acids = true,
    .neutralize_bases = true,
    .preserve_quaternary_n = true,
    .normalize_nitro = true,
    .normalize_sulfoxide = true,
    .normalize_phosphate = true,
    .remove_all_h = false,
    .add_explicit_h = false,
    .check_valence = true,
    .check_charges = false
};

/* Strict sanitization options */
const sanitize_options_t SANITIZE_OPTIONS_STRICT = {
    .flags = SANITIZE_ALL,
    .keep_all_organic = false,
    .min_fragment_atoms = 1,
    .remove_water = true,
    .neutralize_acids = true,
    .neutralize_bases = true,
    .preserve_quaternary_n = false,
    .normalize_nitro = true,
    .normalize_sulfoxide = true,
    .normalize_phosphate = true,
    .remove_all_h = true,
    .add_explicit_h = false,
    .check_valence = true,
    .check_charges = true
};

/* Minimal sanitization options */
const sanitize_options_t SANITIZE_OPTIONS_MINIMAL = {
    .flags = SANITIZE_MINIMAL,
    .keep_all_organic = false,
    .min_fragment_atoms = 3,
    .remove_water = false,
    .neutralize_acids = false,
    .neutralize_bases = false,
    .preserve_quaternary_n = true,
    .normalize_nitro = false,
    .normalize_sulfoxide = false,
    .normalize_phosphate = false,
    .remove_all_h = false,
    .add_explicit_h = false,
    .check_valence = false,
    .check_charges = false
};

/* Default tautomer options */
const tautomer_options_t TAUTOMER_OPTIONS_DEFAULT = {
    .max_tautomers = 100,
    .include_input = true,
    .canonical_only = false,
    .consider_rings = true,
    .keto_enol = true,
    .amide_imidic = true,
    .lactam_lactim = true,
    .nitroso_oxime = true,
    .phosphate = false
};

/* Common inorganic ions (single atoms or small fragments) */
static const element_t SALT_CATIONS[] = {
    ELEM_Na, ELEM_K, ELEM_Li, ELEM_Ca, ELEM_Mg, ELEM_Zn,
    ELEM_Fe, ELEM_Cu, ELEM_Mn, ELEM_Co, ELEM_Ni, ELEM_Al,
    ELEM_Ag, ELEM_Ba, ELEM_Sr, ELEM_Cs, ELEM_Rb
};
#define NUM_SALT_CATIONS (sizeof(SALT_CATIONS) / sizeof(SALT_CATIONS[0]))

static const element_t SALT_ANIONS[] = {
    ELEM_Cl, ELEM_Br, ELEM_I, ELEM_F
};
#define NUM_SALT_ANIONS (sizeof(SALT_ANIONS) / sizeof(SALT_ANIONS[0]))

/* Check if element is a common salt cation */
static bool is_salt_cation(element_t elem) {
    for (size_t i = 0; i < NUM_SALT_CATIONS; i++) {
        if (elem == SALT_CATIONS[i]) return true;
    }
    return false;
}

/* Check if element is a common salt anion */
static bool is_salt_anion(element_t elem) {
    for (size_t i = 0; i < NUM_SALT_ANIONS; i++) {
        if (elem == SALT_ANIONS[i]) return true;
    }
    return false;
}

/* Count heavy atoms in a fragment */
static int count_fragment_heavy_atoms(const molecule_t* mol, int fragment_idx) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->fragment_ids[i] == fragment_idx &&
            mol->atoms[i].element != ELEM_H) {
            count++;
        }
    }
    return count;
}

/* Check if fragment contains carbon */
static bool fragment_contains_carbon(const molecule_t* mol, int fragment_idx) {
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->fragment_ids[i] == fragment_idx &&
            mol->atoms[i].element == ELEM_C) {
            return true;
        }
    }
    return false;
}

/* Check if fragment is water (H2O or HOH) */
static bool fragment_is_water(const molecule_t* mol, int fragment_idx) {
    int oxygen_count = 0;
    int hydrogen_count = 0;
    int other_count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->fragment_ids[i] == fragment_idx) {
            if (mol->atoms[i].element == ELEM_O) {
                oxygen_count++;
            } else if (mol->atoms[i].element == ELEM_H) {
                hydrogen_count++;
            } else {
                other_count++;
            }
            /* Also count implicit hydrogens */
            hydrogen_count += mol->atoms[i].implicit_h_count;
        }
    }

    return (oxygen_count == 1 && hydrogen_count == 2 && other_count == 0);
}

bool molecule_fragment_is_salt(const molecule_t* mol, int fragment_idx) {
    if (!mol || !mol->fragment_ids) return false;
    if (fragment_idx < 0 || fragment_idx >= mol->num_fragments) return false;

    int atom_count = 0;
    int charged_count = 0;
    bool has_organic = false;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->fragment_ids[i] == fragment_idx) {
            atom_count++;
            const atom_t* atom = &mol->atoms[i];

            if (atom->charge != 0) {
                charged_count++;
            }

            /* Check for organic atoms */
            if (atom->element == ELEM_C) {
                has_organic = true;
            }

            /* Single atom ions */
            if (atom_count == 1 && atom->num_neighbors == 0) {
                if (is_salt_cation(atom->element) && atom->charge > 0) {
                    return true;
                }
                if (is_salt_anion(atom->element) && atom->charge < 0) {
                    return true;
                }
            }
        }
    }

    /* Small charged inorganic fragments are likely salts */
    if (!has_organic && atom_count <= 5 && charged_count > 0) {
        return true;
    }

    return false;
}

bool molecule_fragment_is_organic(const molecule_t* mol, int fragment_idx) {
    return fragment_contains_carbon(mol, fragment_idx);
}

cchem_status_t molecule_get_fragment_info(molecule_t* mol,
                                          int* num_fragments,
                                          int* largest_fragment_idx,
                                          int* largest_fragment_atoms) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Ensure fragments are computed */
    if (!mol->fragment_ids) {
        cchem_status_t status = molecule_find_fragments(mol);
        if (status != CCHEM_OK) return status;
    }

    if (num_fragments) {
        *num_fragments = mol->num_fragments;
    }

    /* Find largest fragment */
    int max_atoms = 0;
    int max_idx = 0;

    for (int f = 0; f < mol->num_fragments; f++) {
        int count = count_fragment_heavy_atoms(mol, f);
        if (count > max_atoms) {
            max_atoms = count;
            max_idx = f;
        }
    }

    if (largest_fragment_idx) {
        *largest_fragment_idx = max_idx;
    }
    if (largest_fragment_atoms) {
        *largest_fragment_atoms = max_atoms;
    }

    return CCHEM_OK;
}

cchem_status_t molecule_keep_fragment(molecule_t* mol, int fragment_idx) {
    if (!mol || !mol->fragment_ids) return CCHEM_ERROR_INVALID_INPUT;
    if (fragment_idx < 0 || fragment_idx >= mol->num_fragments) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Mark atoms to remove (those not in the target fragment) */
    bool* remove = (bool*)calloc(mol->num_atoms, sizeof(bool));
    if (!remove) return CCHEM_ERROR_MEMORY;

    int remove_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->fragment_ids[i] != fragment_idx) {
            remove[i] = true;
            remove_count++;
        }
    }

    if (remove_count == 0) {
        free(remove);
        return CCHEM_OK;
    }

    /* Build new molecule with only the desired fragment */
    molecule_t* new_mol = molecule_create_with_capacity(
        mol->num_atoms - remove_count,
        mol->num_bonds);

    if (!new_mol) {
        free(remove);
        return CCHEM_ERROR_MEMORY;
    }

    /* Map old indices to new indices */
    int* index_map = (int*)calloc(mol->num_atoms, sizeof(int));
    if (!index_map) {
        free(remove);
        molecule_free(new_mol);
        return CCHEM_ERROR_MEMORY;
    }

    for (int i = 0; i < mol->num_atoms; i++) {
        index_map[i] = -1;
    }

    /* Copy atoms from the target fragment */
    int new_idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!remove[i]) {
            atom_t atom_copy;
            atom_copy = mol->atoms[i];

            /* Reset neighbor info (will be rebuilt) */
            atom_copy.num_neighbors = 0;
            atom_copy.num_stereo_neighbors = 0;

            molecule_add_atom_full(new_mol, &atom_copy);
            index_map[i] = new_idx++;
        }
    }

    /* Copy bonds within the fragment */
    for (int b = 0; b < mol->num_bonds; b++) {
        int a1 = mol->bonds[b].atom1;
        int a2 = mol->bonds[b].atom2;

        if (index_map[a1] >= 0 && index_map[a2] >= 0) {
            molecule_add_bond(new_mol, index_map[a1], index_map[a2],
                              mol->bonds[b].type);
        }
    }

    /* Update stereo neighbor indices for chirality */
    for (int i = 0; i < new_mol->num_atoms; i++) {
        atom_t* atom = &new_mol->atoms[i];
        int old_idx = -1;

        /* Find original atom index */
        for (int j = 0; j < mol->num_atoms; j++) {
            if (index_map[j] == i) {
                old_idx = j;
                break;
            }
        }

        if (old_idx >= 0 && mol->atoms[old_idx].num_stereo_neighbors > 0) {
            atom->num_stereo_neighbors = 0;
            for (int s = 0; s < mol->atoms[old_idx].num_stereo_neighbors; s++) {
                int old_neighbor = mol->atoms[old_idx].stereo_neighbors[s];
                if (index_map[old_neighbor] >= 0) {
                    atom->stereo_neighbors[atom->num_stereo_neighbors++] =
                        index_map[old_neighbor];
                }
            }
        }
    }

    /* Replace original molecule contents with new */
    /* Swap atoms */
    atom_t* old_atoms = mol->atoms;
    mol->atoms = new_mol->atoms;
    new_mol->atoms = old_atoms;
    mol->num_atoms = new_mol->num_atoms;
    int old_cap = mol->atoms_capacity;
    mol->atoms_capacity = new_mol->atoms_capacity;
    new_mol->atoms_capacity = old_cap;

    /* Swap bonds */
    bond_t* old_bonds = mol->bonds;
    mol->bonds = new_mol->bonds;
    new_mol->bonds = old_bonds;
    mol->num_bonds = new_mol->num_bonds;
    old_cap = mol->bonds_capacity;
    mol->bonds_capacity = new_mol->bonds_capacity;
    new_mol->bonds_capacity = old_cap;

    /* Reset fragment and ring info (needs recomputation) */
    if (mol->fragment_ids) {
        free(mol->fragment_ids);
        mol->fragment_ids = NULL;
    }
    mol->num_fragments = 0;
    mol->rings_computed = false;
    mol->is_canonical = false;

    free(remove);
    free(index_map);
    molecule_free(new_mol);

    return CCHEM_OK;
}

cchem_status_t molecule_unsalt(molecule_t* mol, const sanitize_options_t* options) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    const sanitize_options_t* opts = options ? options : &SANITIZE_OPTIONS_DEFAULT;

    /* Find fragments */
    cchem_status_t status = molecule_find_fragments(mol);
    if (status != CCHEM_OK) return status;

    /* If single fragment, nothing to do */
    if (mol->num_fragments <= 1) {
        return CCHEM_OK;
    }

    /* Calculate implicit hydrogens for water detection */
    molecule_calc_implicit_h(mol);

    /* Find best fragment to keep */
    int best_fragment = -1;
    int best_heavy_atoms = 0;

    for (int f = 0; f < mol->num_fragments; f++) {
        /* Skip salts */
        if (molecule_fragment_is_salt(mol, f)) {
            continue;
        }

        /* Skip water if requested */
        if (opts->remove_water && fragment_is_water(mol, f)) {
            continue;
        }

        /* Skip non-organic if not keeping all organic */
        if (!opts->keep_all_organic && !molecule_fragment_is_organic(mol, f)) {
            continue;
        }

        int heavy_atoms = count_fragment_heavy_atoms(mol, f);

        /* Skip small fragments */
        if (heavy_atoms < opts->min_fragment_atoms) {
            continue;
        }

        /* Track largest */
        if (heavy_atoms > best_heavy_atoms) {
            best_heavy_atoms = heavy_atoms;
            best_fragment = f;
        }
    }

    /* If no good fragment found, keep the largest one overall */
    if (best_fragment < 0) {
        for (int f = 0; f < mol->num_fragments; f++) {
            int heavy_atoms = count_fragment_heavy_atoms(mol, f);
            if (heavy_atoms > best_heavy_atoms) {
                best_heavy_atoms = heavy_atoms;
                best_fragment = f;
            }
        }
    }

    if (best_fragment < 0) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    return molecule_keep_fragment(mol, best_fragment);
}

cchem_status_t molecule_aromatize(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;
    return molecule_perceive_aromaticity(mol);
}

/* Kekulization state for backtracking */
typedef struct {
    molecule_t* mol;
    int* arom_bond_idx;
    int num_arom_bonds;
    bool* bond_is_double;
    bool* has_double;
} kekule_state_t;

/* Check if current assignment is valid for an atom */
static bool kekule_atom_valid(kekule_state_t* state, int atom_idx) {
    atom_t* atom = &state->mol->atoms[atom_idx];

    /* Count double bonds to this atom */
    int double_count = 0;
    if (state->has_double[atom_idx]) {
        double_count++;
    }

    for (int i = 0; i < state->num_arom_bonds; i++) {
        int b = state->arom_bond_idx[i];
        bond_t* bond = &state->mol->bonds[b];
        if ((bond->atom1 == atom_idx || bond->atom2 == atom_idx) &&
            state->bond_is_double[i]) {
            double_count++;
        }
    }

    /* Each sp2 carbon should have exactly one double bond */
    /* Heteroatoms may have zero (lone pair) or one */
    if (atom->element == ELEM_C) {
        return double_count <= 1;
    } else if (atom->element == ELEM_N || atom->element == ELEM_O ||
               atom->element == ELEM_S) {
        return double_count <= 1;
    }

    return true;
}

/* Check if assignment is complete and valid */
static bool kekule_is_complete(kekule_state_t* state) {
    /* Check all aromatic atoms have valid configuration */
    bool* checked = (bool*)calloc(state->mol->num_atoms, sizeof(bool));

    for (int i = 0; i < state->num_arom_bonds; i++) {
        int b = state->arom_bond_idx[i];
        bond_t* bond = &state->mol->bonds[b];

        if (!checked[bond->atom1]) {
            checked[bond->atom1] = true;
            int double_count = 0;
            if (state->has_double[bond->atom1]) double_count++;

            for (int j = 0; j < state->num_arom_bonds; j++) {
                int b2 = state->arom_bond_idx[j];
                bond_t* bond2 = &state->mol->bonds[b2];
                if ((bond2->atom1 == bond->atom1 || bond2->atom2 == bond->atom1) &&
                    state->bond_is_double[j]) {
                    double_count++;
                }
            }

            /* Carbon must have exactly one double bond */
            if (state->mol->atoms[bond->atom1].element == ELEM_C &&
                double_count != 1) {
                free(checked);
                return false;
            }
        }

        if (!checked[bond->atom2]) {
            checked[bond->atom2] = true;
            int double_count = 0;
            if (state->has_double[bond->atom2]) double_count++;

            for (int j = 0; j < state->num_arom_bonds; j++) {
                int b2 = state->arom_bond_idx[j];
                bond_t* bond2 = &state->mol->bonds[b2];
                if ((bond2->atom1 == bond->atom2 || bond2->atom2 == bond->atom2) &&
                    state->bond_is_double[j]) {
                    double_count++;
                }
            }

            if (state->mol->atoms[bond->atom2].element == ELEM_C &&
                double_count != 1) {
                free(checked);
                return false;
            }
        }
    }

    free(checked);
    return true;
}

/* Backtracking solver for Kekule assignment */
static bool kekule_solve(kekule_state_t* state, int bond_idx) {
    if (bond_idx >= state->num_arom_bonds) {
        return kekule_is_complete(state);
    }

    int b = state->arom_bond_idx[bond_idx];
    bond_t* bond = &state->mol->bonds[b];

    /* Try assigning as double bond */
    state->bond_is_double[bond_idx] = true;
    if (kekule_atom_valid(state, bond->atom1) &&
        kekule_atom_valid(state, bond->atom2)) {
        if (kekule_solve(state, bond_idx + 1)) {
            return true;
        }
    }

    /* Try assigning as single bond */
    state->bond_is_double[bond_idx] = false;
    if (kekule_solve(state, bond_idx + 1)) {
        return true;
    }

    return false;
}

cchem_status_t molecule_kekulize(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Collect aromatic bonds */
    kekule_state_t state;
    state.mol = mol;
    state.arom_bond_idx = (int*)malloc(mol->num_bonds * sizeof(int));
    if (!state.arom_bond_idx) return CCHEM_ERROR_MEMORY;

    state.num_arom_bonds = 0;

    /* Track which atoms are in aromatic system */
    bool* in_aromatic = (bool*)calloc(mol->num_atoms, sizeof(bool));
    if (!in_aromatic) {
        free(state.arom_bond_idx);
        return CCHEM_ERROR_MEMORY;
    }

    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->bonds[b].type == BOND_AROMATIC ||
            mol->bonds[b].type == BOND_RING_AROMATIC) {
            state.arom_bond_idx[state.num_arom_bonds++] = b;
            in_aromatic[mol->bonds[b].atom1] = true;
            in_aromatic[mol->bonds[b].atom2] = true;
        }
    }

    if (state.num_arom_bonds == 0) {
        free(state.arom_bond_idx);
        free(in_aromatic);
        return CCHEM_OK;
    }

    /* Mark atoms that have double bonds from NON-aromatic bonds */
    state.has_double = (bool*)calloc(mol->num_atoms, sizeof(bool));
    if (!state.has_double) {
        free(state.arom_bond_idx);
        free(in_aromatic);
        return CCHEM_ERROR_MEMORY;
    }

    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->bonds[b].type == BOND_DOUBLE) {
            int a1 = mol->bonds[b].atom1;
            int a2 = mol->bonds[b].atom2;

            if (!in_aromatic[a1]) {
                state.has_double[a1] = true;
            }
            if (!in_aromatic[a2]) {
                state.has_double[a2] = true;
            }
            /* Exocyclic double bond satisfies the aromatic atom */
            if (in_aromatic[a1] && !in_aromatic[a2]) {
                state.has_double[a1] = true;
            }
            if (in_aromatic[a2] && !in_aromatic[a1]) {
                state.has_double[a2] = true;
            }
        }
    }

    free(in_aromatic);

    /* Initialize bond assignment */
    state.bond_is_double = (bool*)calloc(state.num_arom_bonds, sizeof(bool));
    if (!state.bond_is_double) {
        free(state.arom_bond_idx);
        free(state.has_double);
        return CCHEM_ERROR_MEMORY;
    }

    /* Solve using backtracking */
    bool success = kekule_solve(&state, 0);

    if (success) {
        /* Apply the solution */
        for (int i = 0; i < state.num_arom_bonds; i++) {
            int b = state.arom_bond_idx[i];
            mol->bonds[b].type = state.bond_is_double[i] ? BOND_DOUBLE : BOND_SINGLE;
            mol->bonds[b].aromatic = false;
        }

        /* Clear aromatic flags on atoms */
        for (int i = 0; i < mol->num_atoms; i++) {
            mol->atoms[i].aromatic = false;
        }
    }

    free(state.bond_is_double);
    free(state.has_double);
    free(state.arom_bond_idx);

    return success ? CCHEM_OK : CCHEM_ERROR_INVALID_INPUT;
}

cchem_status_t molecule_neutralize(molecule_t* mol, const sanitize_options_t* options) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    const sanitize_options_t* opts = options ? options : &SANITIZE_OPTIONS_DEFAULT;

    for (int i = 0; i < mol->num_atoms; i++) {
        atom_t* atom = &mol->atoms[i];

        if (atom->charge == 0) continue;

        /* Skip quaternary nitrogen if requested */
        if (opts->preserve_quaternary_n &&
            atom->element == ELEM_N &&
            atom->charge > 0 &&
            atom->num_neighbors >= 4) {
            continue;
        }

        /* Neutralize positively charged atoms */
        if (atom->charge > 0 && opts->neutralize_bases) {
            /* R-NH3+ -> R-NH2: reduce charge, remove H */
            if (atom->element == ELEM_N) {
                /* First try to remove from explicit h_count */
                while (atom->charge > 0 && atom->h_count > 0) {
                    atom->charge--;
                    atom->h_count--;
                }
                /* Then from implicit */
                while (atom->charge > 0 && atom->implicit_h_count > 0) {
                    atom->charge--;
                    atom->implicit_h_count--;
                }
            }
            /* R-OH2+ -> R-OH */
            else if (atom->element == ELEM_O) {
                while (atom->charge > 0 && atom->h_count > 0) {
                    atom->charge--;
                    atom->h_count--;
                }
                while (atom->charge > 0 && atom->implicit_h_count > 0) {
                    atom->charge--;
                    atom->implicit_h_count--;
                }
            }
        }

        /* Neutralize negatively charged atoms */
        if (atom->charge < 0 && opts->neutralize_acids) {
            /* R-O- -> R-OH: add proton */
            if (atom->element == ELEM_O || atom->element == ELEM_S) {
                while (atom->charge < 0) {
                    atom->charge++;
                    /* For bracket atoms, use h_count; otherwise implicit */
                    if (atom->h_count >= 0) {
                        atom->h_count++;
                    } else {
                        atom->implicit_h_count++;
                    }
                }
            }
            /* R-N- -> R-NH: add proton (rare) */
            else if (atom->element == ELEM_N) {
                while (atom->charge < 0) {
                    atom->charge++;
                    if (atom->h_count >= 0) {
                        atom->h_count++;
                    } else {
                        atom->implicit_h_count++;
                    }
                }
            }
        }
    }

    return CCHEM_OK;
}

cchem_status_t molecule_normalize(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Normalize nitro groups: N(=O)=O or [N+](=O)[O-] -> canonical form */
    for (int i = 0; i < mol->num_atoms; i++) {
        atom_t* atom = &mol->atoms[i];

        if (atom->element != ELEM_N) continue;

        /* Check for nitro pattern: N connected to 2 oxygens */
        int oxygen_count = 0;
        int oxygen_neighbors[2] = {-1, -1};

        for (int j = 0; j < atom->num_neighbors; j++) {
            int neighbor_idx = atom->neighbors[j];
            atom_t* neighbor = &mol->atoms[neighbor_idx];

            if (neighbor->element == ELEM_O && neighbor->num_neighbors == 1) {
                if (oxygen_count < 2) {
                    oxygen_neighbors[oxygen_count] = neighbor_idx;
                }
                oxygen_count++;
            }
        }

        /* Found nitro group */
        if (oxygen_count == 2 && atom->num_neighbors == 3) {
            /* Normalize to [N+](=O)[O-] form */
            atom->charge = 1;

            /* One oxygen double bonded, one single with negative charge */
            for (int j = 0; j < 2; j++) {
                if (oxygen_neighbors[j] >= 0) {
                    int bond_idx = atom_get_bond_to(atom, oxygen_neighbors[j]);
                    if (j == 0) {
                        mol->bonds[bond_idx].type = BOND_DOUBLE;
                        mol->atoms[oxygen_neighbors[j]].charge = 0;
                    } else {
                        mol->bonds[bond_idx].type = BOND_SINGLE;
                        mol->atoms[oxygen_neighbors[j]].charge = -1;
                    }
                }
            }
        }
    }

    return CCHEM_OK;
}

cchem_status_t molecule_remove_explicit_h(molecule_t* mol,
                                          const sanitize_options_t* options) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    const sanitize_options_t* opts = options ? options : &SANITIZE_OPTIONS_DEFAULT;

    /* Find hydrogens to remove */
    bool* remove = (bool*)calloc(mol->num_atoms, sizeof(bool));
    if (!remove) return CCHEM_ERROR_MEMORY;

    int remove_count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        atom_t* atom = &mol->atoms[i];

        if (atom->element != ELEM_H) continue;

        /* Don't remove if it's the only atom */
        if (mol->num_atoms == 1) continue;

        /* Check if hydrogen should be preserved */
        bool preserve = false;

        if (!opts->remove_all_h) {
            /* Preserve H with isotope label */
            if (atom->isotope > 0) {
                preserve = true;
            }

            /* Preserve H involved in stereochemistry */
            if (atom->num_stereo_neighbors > 0) {
                preserve = true;
            }

            /* Check if attached atom has chirality */
            if (atom->num_neighbors > 0) {
                int parent_idx = atom->neighbors[0];
                if (mol->atoms[parent_idx].chirality != CHIRALITY_NONE) {
                    preserve = true;
                }
            }
        }

        if (!preserve) {
            remove[i] = true;
            remove_count++;

            /* Increment implicit H count on parent atom */
            if (atom->num_neighbors > 0) {
                int parent_idx = atom->neighbors[0];
                mol->atoms[parent_idx].implicit_h_count++;
            }
        }
    }

    if (remove_count == 0) {
        free(remove);
        return CCHEM_OK;
    }

    /* Remove the hydrogens by rebuilding molecule */
    molecule_t* new_mol = molecule_create_with_capacity(
        mol->num_atoms - remove_count,
        mol->num_bonds);

    if (!new_mol) {
        free(remove);
        return CCHEM_ERROR_MEMORY;
    }

    /* Map old indices to new */
    int* index_map = (int*)calloc(mol->num_atoms, sizeof(int));
    if (!index_map) {
        free(remove);
        molecule_free(new_mol);
        return CCHEM_ERROR_MEMORY;
    }

    for (int i = 0; i < mol->num_atoms; i++) {
        index_map[i] = -1;
    }

    /* Copy non-removed atoms */
    int new_idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!remove[i]) {
            atom_t copy = mol->atoms[i];
            copy.num_neighbors = 0;
            copy.num_stereo_neighbors = 0;

            molecule_add_atom_full(new_mol, &copy);
            index_map[i] = new_idx++;
        }
    }

    /* Copy bonds (only between remaining atoms) */
    for (int b = 0; b < mol->num_bonds; b++) {
        int a1 = mol->bonds[b].atom1;
        int a2 = mol->bonds[b].atom2;

        if (index_map[a1] >= 0 && index_map[a2] >= 0) {
            molecule_add_bond(new_mol, index_map[a1], index_map[a2],
                              mol->bonds[b].type);
        }
    }

    /* Swap contents */
    atom_t* old_atoms = mol->atoms;
    mol->atoms = new_mol->atoms;
    new_mol->atoms = old_atoms;
    mol->num_atoms = new_mol->num_atoms;
    int old_cap = mol->atoms_capacity;
    mol->atoms_capacity = new_mol->atoms_capacity;
    new_mol->atoms_capacity = old_cap;

    bond_t* old_bonds = mol->bonds;
    mol->bonds = new_mol->bonds;
    new_mol->bonds = old_bonds;
    mol->num_bonds = new_mol->num_bonds;
    old_cap = mol->bonds_capacity;
    mol->bonds_capacity = new_mol->bonds_capacity;
    new_mol->bonds_capacity = old_cap;

    /* Reset state */
    if (mol->fragment_ids) {
        free(mol->fragment_ids);
        mol->fragment_ids = NULL;
    }
    mol->num_fragments = 0;
    mol->rings_computed = false;
    mol->is_canonical = false;

    free(remove);
    free(index_map);
    molecule_free(new_mol);

    return CCHEM_OK;
}

cchem_status_t molecule_add_explicit_h(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Calculate implicit hydrogens if not done */
    molecule_calc_implicit_h(mol);

    /* Count total hydrogens to add */
    int total_h = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        total_h += mol->atoms[i].implicit_h_count;
    }

    if (total_h == 0) {
        return CCHEM_OK;
    }

    /* Add explicit hydrogens */
    for (int i = 0; i < mol->num_atoms; i++) {
        int h_count = mol->atoms[i].implicit_h_count;

        for (int h = 0; h < h_count; h++) {
            int h_idx = molecule_add_atom(mol, ELEM_H);
            if (h_idx < 0) {
                return CCHEM_ERROR_MEMORY;
            }

            if (molecule_add_bond(mol, i, h_idx, BOND_SINGLE) < 0) {
                return CCHEM_ERROR_MEMORY;
            }
        }

        mol->atoms[i].implicit_h_count = 0;
    }

    return CCHEM_OK;
}

cchem_status_t molecule_remove_stereo(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Clear atom chirality */
    for (int i = 0; i < mol->num_atoms; i++) {
        mol->atoms[i].chirality = CHIRALITY_NONE;
        mol->atoms[i].chirality_class = 0;
        mol->atoms[i].num_stereo_neighbors = 0;
    }

    /* Clear bond stereo */
    for (int b = 0; b < mol->num_bonds; b++) {
        mol->bonds[b].stereo = STEREO_NONE;
        mol->bonds[b].stereo_type = BOND_NONE;
        mol->bonds[b].stereo_atom = -1;
        if (mol->bonds[b].type == BOND_UP || mol->bonds[b].type == BOND_DOWN) {
            mol->bonds[b].type = BOND_SINGLE;
        }
    }

    return CCHEM_OK;
}

cchem_status_t molecule_remove_isotopes(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    for (int i = 0; i < mol->num_atoms; i++) {
        mol->atoms[i].isotope = 0;
    }

    return CCHEM_OK;
}

cchem_status_t molecule_sanitize(molecule_t* mol,
                                 const sanitize_options_t* options,
                                 char* error_buf, size_t error_buf_size) {
    if (!mol) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "NULL molecule");
        }
        return CCHEM_ERROR_INVALID_INPUT;
    }

    const sanitize_options_t* opts = options ? options : &SANITIZE_OPTIONS_DEFAULT;
    cchem_status_t status;

    /* 1. Validate structure first */
    if (opts->flags & SANITIZE_VALIDATE) {
        status = molecule_validate(mol);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Structure validation failed");
            }
            return status;
        }
    }

    /* 2. Remove salts/fragments */
    if (opts->flags & SANITIZE_UNSALT) {
        status = molecule_unsalt(mol, opts);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Salt removal failed");
            }
            return status;
        }
    }

    /* 3. Normalize functional groups */
    if (opts->flags & SANITIZE_NORMALIZE) {
        status = molecule_normalize(mol);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Normalization failed");
            }
            return status;
        }
    }

    /* 4. Neutralize charges */
    if (opts->flags & SANITIZE_NEUTRALIZE) {
        status = molecule_neutralize(mol, opts);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Neutralization failed");
            }
            return status;
        }
    }

    /* 5. Handle aromaticity (mutually exclusive) */
    if (opts->flags & SANITIZE_AROMATIZE) {
        status = molecule_aromatize(mol);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Aromatization failed");
            }
            return status;
        }
    } else if (opts->flags & SANITIZE_KEKULIZE) {
        status = molecule_kekulize(mol);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Kekulization failed");
            }
            return status;
        }
    }

    /* 6. Handle explicit hydrogens */
    if (opts->flags & SANITIZE_REMOVE_H) {
        status = molecule_remove_explicit_h(mol, opts);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Hydrogen removal failed");
            }
            return status;
        }
    }

    /* 7. Remove stereochemistry */
    if (opts->flags & SANITIZE_REMOVE_STEREO) {
        status = molecule_remove_stereo(mol);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Stereo removal failed");
            }
            return status;
        }
    }

    /* 8. Remove isotopes */
    if (opts->flags & SANITIZE_REMOVE_ISOTOPES) {
        status = molecule_remove_isotopes(mol);
        if (status != CCHEM_OK) {
            if (error_buf && error_buf_size > 0) {
                snprintf(error_buf, error_buf_size, "Isotope removal failed");
            }
            return status;
        }
    }

    return CCHEM_OK;
}

/* Tautomer enumeration - basic implementation */

/* Tautomeric sites: atoms that can participate in hydrogen shifts */
typedef struct {
    int atom_idx;
    bool is_donor;              /* Can donate H */
    bool is_acceptor;           /* Can accept H */
    int h_count;                /* Current H count */
} tautomer_site_t;

/* Find tautomeric sites in molecule */
static int find_tautomer_sites(const molecule_t* mol, tautomer_site_t* sites, int max_sites) {
    int num_sites = 0;

    for (int i = 0; i < mol->num_atoms && num_sites < max_sites; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Nitrogen with H */
        if (atom->element == ELEM_N) {
            int total_h = atom->implicit_h_count + atom->h_count;
            if (total_h > 0) {
                sites[num_sites].atom_idx = i;
                sites[num_sites].is_donor = true;
                sites[num_sites].is_acceptor = (atom->num_neighbors < 3);
                sites[num_sites].h_count = total_h;
                num_sites++;
            } else if (atom->num_neighbors < 3) {
                /* N without H but can accept */
                sites[num_sites].atom_idx = i;
                sites[num_sites].is_donor = false;
                sites[num_sites].is_acceptor = true;
                sites[num_sites].h_count = 0;
                num_sites++;
            }
        }
        /* Oxygen with H */
        else if (atom->element == ELEM_O) {
            int total_h = atom->implicit_h_count + atom->h_count;
            if (total_h > 0) {
                sites[num_sites].atom_idx = i;
                sites[num_sites].is_donor = true;
                sites[num_sites].is_acceptor = false;
                sites[num_sites].h_count = total_h;
                num_sites++;
            }
        }
        /* Carbon adjacent to C=O or C=N (keto-enol) */
        else if (atom->element == ELEM_C) {
            bool has_unsaturated_neighbor = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                int bond_idx = atom->neighbor_bonds[j];
                if (mol->bonds[bond_idx].type == BOND_DOUBLE) {
                    int neighbor = atom->neighbors[j];
                    if (mol->atoms[neighbor].element == ELEM_O ||
                        mol->atoms[neighbor].element == ELEM_N) {
                        has_unsaturated_neighbor = true;
                        break;
                    }
                }
            }

            if (has_unsaturated_neighbor) {
                int total_h = atom->implicit_h_count + atom->h_count;
                if (total_h > 0) {
                    sites[num_sites].atom_idx = i;
                    sites[num_sites].is_donor = true;
                    sites[num_sites].is_acceptor = true;
                    sites[num_sites].h_count = total_h;
                    num_sites++;
                }
            }
        }
    }

    return num_sites;
}

cchem_status_t tautomer_enumerate(const molecule_t* mol,
                                  const tautomer_options_t* options,
                                  tautomer_result_t* result) {
    if (!mol || !result) return CCHEM_ERROR_INVALID_INPUT;

    const tautomer_options_t* opts = options ? options : &TAUTOMER_OPTIONS_DEFAULT;

    memset(result, 0, sizeof(tautomer_result_t));

    /* Find tautomeric sites */
    tautomer_site_t sites[64];
    int num_sites = find_tautomer_sites(mol, sites, 64);

    if (num_sites < 2) {
        /* No tautomerism possible - return just the input */
        if (opts->include_input) {
            result->smiles = (char**)malloc(sizeof(char*));
            result->scores = (double*)malloc(sizeof(double));
            if (!result->smiles || !result->scores) {
                if (result->smiles) free(result->smiles);
                if (result->scores) free(result->scores);
                return CCHEM_ERROR_MEMORY;
            }

            result->smiles[0] = molecule_to_canonical_smiles(mol, NULL);
            result->scores[0] = 0.0;
            result->num_tautomers = 1;
            result->canonical_idx = 0;
        }
        return CCHEM_OK;
    }

    /* Allocate result arrays */
    int max_tautomers = opts->max_tautomers;
    result->smiles = (char**)calloc(max_tautomers, sizeof(char*));
    result->scores = (double*)calloc(max_tautomers, sizeof(double));

    if (!result->smiles || !result->scores) {
        if (result->smiles) free(result->smiles);
        if (result->scores) free(result->scores);
        return CCHEM_ERROR_MEMORY;
    }

    /* Add input tautomer */
    if (opts->include_input) {
        result->smiles[result->num_tautomers] = molecule_to_canonical_smiles(mol, NULL);
        result->scores[result->num_tautomers] = 0.0;
        result->num_tautomers++;
    }

    /* Generate tautomers by H shifts between donor/acceptor pairs */
    /* This is a simplified enumeration - full implementation would use
     * more sophisticated rules and scoring */
    for (int d = 0; d < num_sites && result->num_tautomers < max_tautomers; d++) {
        if (!sites[d].is_donor || sites[d].h_count == 0) continue;

        for (int a = 0; a < num_sites && result->num_tautomers < max_tautomers; a++) {
            if (d == a) continue;
            if (!sites[a].is_acceptor) continue;

            /* Clone molecule and perform H shift */
            molecule_t* taut = molecule_clone(mol);
            if (!taut) continue;

            int donor_idx = sites[d].atom_idx;
            int acceptor_idx = sites[a].atom_idx;

            /* Check if sites are connected (1,3 or 1,5 shift) */
            bool connected = false;
            for (int n = 0; n < taut->atoms[donor_idx].num_neighbors; n++) {
                int intermediate = taut->atoms[donor_idx].neighbors[n];
                for (int m = 0; m < taut->atoms[intermediate].num_neighbors; m++) {
                    if (taut->atoms[intermediate].neighbors[m] == acceptor_idx) {
                        connected = true;
                        break;
                    }
                }
                if (connected) break;
            }

            if (connected) {
                /* Perform H shift: donor loses H, acceptor gains H */
                taut->atoms[donor_idx].implicit_h_count--;
                taut->atoms[acceptor_idx].implicit_h_count++;

                /* Adjust bond orders for keto-enol */
                /* Find bond between donor and intermediate */
                for (int n = 0; n < taut->atoms[donor_idx].num_neighbors; n++) {
                    int intermediate = taut->atoms[donor_idx].neighbors[n];

                    /* Check if intermediate connects to acceptor */
                    int acc_bond = molecule_get_bond_index(taut, intermediate, acceptor_idx);
                    if (acc_bond >= 0) {
                        int don_bond = taut->atoms[donor_idx].neighbor_bonds[n];

                        /* Swap bond orders */
                        if (taut->bonds[acc_bond].type == BOND_DOUBLE) {
                            taut->bonds[acc_bond].type = BOND_SINGLE;
                            taut->bonds[don_bond].type = BOND_DOUBLE;
                        }
                        break;
                    }
                }

                /* Generate SMILES */
                char* smiles = molecule_to_canonical_smiles(taut, NULL);
                if (smiles) {
                    /* Check for duplicates */
                    bool duplicate = false;
                    for (int i = 0; i < result->num_tautomers; i++) {
                        if (strcmp(result->smiles[i], smiles) == 0) {
                            duplicate = true;
                            break;
                        }
                    }

                    if (!duplicate) {
                        result->smiles[result->num_tautomers] = smiles;
                        /* Simple score: prefer keto over enol */
                        result->scores[result->num_tautomers] = 1.0;
                        result->num_tautomers++;
                    } else {
                        free(smiles);
                    }
                }
            }

            molecule_free(taut);
        }
    }

    /* Find canonical tautomer (lowest score, or lexicographically first) */
    double min_score = result->scores[0];
    result->canonical_idx = 0;

    for (int i = 1; i < result->num_tautomers; i++) {
        if (result->scores[i] < min_score ||
            (result->scores[i] == min_score &&
             strcmp(result->smiles[i], result->smiles[result->canonical_idx]) < 0)) {
            min_score = result->scores[i];
            result->canonical_idx = i;
        }
    }

    return CCHEM_OK;
}

char* tautomer_canonical(const char* smiles,
                         const tautomer_options_t* options,
                         char* error_buf, size_t error_buf_size) {
    if (!smiles) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "NULL SMILES");
        }
        return NULL;
    }

    /* Parse molecule */
    molecule_t* mol = smiles_to_molecule(smiles, error_buf, error_buf_size);
    if (!mol) return NULL;

    /* Enumerate tautomers */
    tautomer_options_t opts = options ? *options : TAUTOMER_OPTIONS_DEFAULT;
    opts.canonical_only = true;

    tautomer_result_t result;
    cchem_status_t status = tautomer_enumerate(mol, &opts, &result);
    molecule_free(mol);

    if (status != CCHEM_OK || result.num_tautomers == 0) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "Tautomer enumeration failed");
        }
        tautomer_result_free(&result);
        return NULL;
    }

    /* Extract canonical tautomer */
    char* canonical = strdup(result.smiles[result.canonical_idx]);
    tautomer_result_free(&result);

    return canonical;
}

void tautomer_result_free(tautomer_result_t* result) {
    if (!result) return;

    if (result->smiles) {
        for (int i = 0; i < result->num_tautomers; i++) {
            if (result->smiles[i]) free(result->smiles[i]);
        }
        free(result->smiles);
    }

    if (result->scores) {
        free(result->scores);
    }

    memset(result, 0, sizeof(tautomer_result_t));
}

char* smiles_sanitize(const char* smiles, sanitize_flags_t flags,
                      char* error_buf, size_t error_buf_size) {
    if (!smiles) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "NULL SMILES");
        }
        return NULL;
    }

    /* Parse molecule */
    molecule_t* mol = smiles_to_molecule(smiles, error_buf, error_buf_size);
    if (!mol) return NULL;

    /* Set up options from flags */
    sanitize_options_t options = SANITIZE_OPTIONS_DEFAULT;
    options.flags = flags;

    /* Sanitize */
    cchem_status_t status = molecule_sanitize(mol, &options, error_buf, error_buf_size);
    if (status != CCHEM_OK) {
        molecule_free(mol);
        return NULL;
    }

    /* Generate canonical SMILES */
    char* result = molecule_to_canonical_smiles(mol, NULL);
    molecule_free(mol);

    if (!result && error_buf && error_buf_size > 0) {
        snprintf(error_buf, error_buf_size, "SMILES generation failed");
    }

    return result;
}

cchem_status_t sanitize_parse_flags(const char* str, sanitize_flags_t* flags) {
    if (!str || !flags) return CCHEM_ERROR_INVALID_INPUT;

    *flags = SANITIZE_NONE;

    /* Handle special preset names */
    if (strcasecmp(str, "complete") == 0) {
        *flags = SANITIZE_COMPLETE;
        return CCHEM_OK;
    }
    if (strcasecmp(str, "standard") == 0) {
        *flags = SANITIZE_STANDARD;
        return CCHEM_OK;
    }
    if (strcasecmp(str, "minimal") == 0) {
        *flags = SANITIZE_MINIMAL;
        return CCHEM_OK;
    }
    if (strcasecmp(str, "all") == 0) {
        *flags = SANITIZE_ALL;
        return CCHEM_OK;
    }
    if (strcasecmp(str, "none") == 0) {
        *flags = SANITIZE_NONE;
        return CCHEM_OK;
    }

    /* Parse comma-separated list */
    char* copy = strdup(str);
    if (!copy) return CCHEM_ERROR_MEMORY;

    char* token = strtok(copy, ",");
    while (token) {
        /* Trim whitespace */
        while (*token && isspace(*token)) token++;
        char* end = token + strlen(token) - 1;
        while (end > token && isspace(*end)) *end-- = '\0';

        if (strcasecmp(token, "unsalt") == 0) {
            *flags |= SANITIZE_UNSALT;
        } else if (strcasecmp(token, "aromatize") == 0) {
            *flags |= SANITIZE_AROMATIZE;
        } else if (strcasecmp(token, "kekulize") == 0) {
            *flags |= SANITIZE_KEKULIZE;
        } else if (strcasecmp(token, "neutralize") == 0) {
            *flags |= SANITIZE_NEUTRALIZE;
        } else if (strcasecmp(token, "remove-stereo") == 0 ||
                   strcasecmp(token, "removestereo") == 0) {
            *flags |= SANITIZE_REMOVE_STEREO;
        } else if (strcasecmp(token, "remove-isotopes") == 0 ||
                   strcasecmp(token, "removeisotopes") == 0) {
            *flags |= SANITIZE_REMOVE_ISOTOPES;
        } else if (strcasecmp(token, "remove-h") == 0 ||
                   strcasecmp(token, "removeh") == 0) {
            *flags |= SANITIZE_REMOVE_H;
        } else if (strcasecmp(token, "normalize") == 0) {
            *flags |= SANITIZE_NORMALIZE;
        } else if (strcasecmp(token, "validate") == 0) {
            *flags |= SANITIZE_VALIDATE;
        } else if (strcasecmp(token, "cleanup") == 0) {
            *flags |= SANITIZE_CLEANUP;
        } else {
            free(copy);
            return CCHEM_ERROR_INVALID_INPUT;
        }

        token = strtok(NULL, ",");
    }

    free(copy);
    return CCHEM_OK;
}

void sanitize_flags_to_string(sanitize_flags_t flags, char* buf, size_t buf_size) {
    if (!buf || buf_size == 0) return;

    buf[0] = '\0';
    size_t pos = 0;

    if (flags == SANITIZE_NONE) {
        snprintf(buf, buf_size, "none");
        return;
    }

    if (flags == SANITIZE_ALL) {
        snprintf(buf, buf_size, "all");
        return;
    }

    #define ADD_FLAG(flag, name) \
        if (flags & flag) { \
            if (pos > 0 && pos < buf_size - 1) buf[pos++] = ','; \
            int written = snprintf(buf + pos, buf_size - pos, "%s", name); \
            if (written > 0) pos += written; \
        }

    ADD_FLAG(SANITIZE_UNSALT, "unsalt")
    ADD_FLAG(SANITIZE_AROMATIZE, "aromatize")
    ADD_FLAG(SANITIZE_KEKULIZE, "kekulize")
    ADD_FLAG(SANITIZE_NEUTRALIZE, "neutralize")
    ADD_FLAG(SANITIZE_REMOVE_STEREO, "remove-stereo")
    ADD_FLAG(SANITIZE_REMOVE_ISOTOPES, "remove-isotopes")
    ADD_FLAG(SANITIZE_REMOVE_H, "remove-h")
    ADD_FLAG(SANITIZE_NORMALIZE, "normalize")
    ADD_FLAG(SANITIZE_VALIDATE, "validate")
    ADD_FLAG(SANITIZE_CLEANUP, "cleanup")

    #undef ADD_FLAG
}
