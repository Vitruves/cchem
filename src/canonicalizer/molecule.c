/**
 * @file molecule.c
 * @brief Molecular graph data structure implementation
 */

#include "cchem/compat.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/molecule.h"

#define INITIAL_ATOM_CAPACITY 64
#define INITIAL_BOND_CAPACITY 64
#define INITIAL_RING_CAPACITY 16

molecule_t* molecule_create(void) {
    return molecule_create_with_capacity(INITIAL_ATOM_CAPACITY, INITIAL_BOND_CAPACITY);
}

molecule_t* molecule_create_with_capacity(int atom_capacity, int bond_capacity) {
    molecule_t* mol = (molecule_t*)calloc(1, sizeof(molecule_t));
    if (!mol) return NULL;

    mol->atoms_capacity = atom_capacity > 0 ? atom_capacity : INITIAL_ATOM_CAPACITY;
    mol->bonds_capacity = bond_capacity > 0 ? bond_capacity : INITIAL_BOND_CAPACITY;
    mol->rings_capacity = INITIAL_RING_CAPACITY;

    mol->atoms = (atom_t*)calloc(mol->atoms_capacity, sizeof(atom_t));
    mol->bonds = (bond_t*)calloc(mol->bonds_capacity, sizeof(bond_t));
    mol->rings = (ring_t*)calloc(mol->rings_capacity, sizeof(ring_t));

    if (!mol->atoms || !mol->bonds || !mol->rings) {
        molecule_free(mol);
        return NULL;
    }

    mol->num_atoms = 0;
    mol->num_bonds = 0;
    mol->num_rings = 0;
    mol->num_fragments = 0;
    mol->fragment_ids = NULL;
    mol->canon_order = NULL;
    mol->is_canonical = false;
    mol->rings_computed = false;
    mol->original_smiles = NULL;
    mol->total_charge = 0;
    mol->molecular_weight = 0.0;
    mol->error_msg[0] = '\0';

    return mol;
}

void molecule_free(molecule_t* mol) {
    if (!mol) return;

    if (mol->atoms) free(mol->atoms);
    if (mol->bonds) free(mol->bonds);
    if (mol->rings) free(mol->rings);
    if (mol->fragment_ids) free(mol->fragment_ids);
    if (mol->canon_order) free(mol->canon_order);
    if (mol->original_smiles) free(mol->original_smiles);

    free(mol);
}

void molecule_reset(molecule_t* mol) {
    if (!mol) return;

    /* Clear atom data */
    for (int i = 0; i < mol->num_atoms; i++) {
        atom_reset(&mol->atoms[i]);
    }

    /* Clear bond data */
    for (int i = 0; i < mol->num_bonds; i++) {
        bond_reset(&mol->bonds[i]);
    }

    mol->num_atoms = 0;
    mol->num_bonds = 0;
    mol->num_rings = 0;
    mol->rings_computed = false;

    if (mol->fragment_ids) {
        free(mol->fragment_ids);
        mol->fragment_ids = NULL;
    }
    mol->num_fragments = 0;

    if (mol->canon_order) {
        free(mol->canon_order);
        mol->canon_order = NULL;
    }
    mol->is_canonical = false;

    if (mol->original_smiles) {
        free(mol->original_smiles);
        mol->original_smiles = NULL;
    }

    mol->total_charge = 0;
    mol->molecular_weight = 0.0;
    mol->error_msg[0] = '\0';
}

static int molecule_ensure_atom_capacity(molecule_t* mol, int needed) {
    if (mol->atoms_capacity >= needed) return 0;

    int new_capacity = mol->atoms_capacity * 2;
    while (new_capacity < needed) new_capacity *= 2;

    atom_t* new_atoms = (atom_t*)realloc(mol->atoms, new_capacity * sizeof(atom_t));
    if (!new_atoms) return -1;

    mol->atoms = new_atoms;
    mol->atoms_capacity = new_capacity;
    return 0;
}

static int molecule_ensure_bond_capacity(molecule_t* mol, int needed) {
    if (mol->bonds_capacity >= needed) return 0;

    int new_capacity = mol->bonds_capacity * 2;
    while (new_capacity < needed) new_capacity *= 2;

    bond_t* new_bonds = (bond_t*)realloc(mol->bonds, new_capacity * sizeof(bond_t));
    if (!new_bonds) return -1;

    mol->bonds = new_bonds;
    mol->bonds_capacity = new_capacity;
    return 0;
}

int molecule_add_atom(molecule_t* mol, element_t element) {
    if (!mol) return -1;

    if (molecule_ensure_atom_capacity(mol, mol->num_atoms + 1) < 0) {
        snprintf(mol->error_msg, sizeof(mol->error_msg), "Failed to allocate memory for atom");
        return -1;
    }

    int idx = mol->num_atoms;
    atom_init(&mol->atoms[idx], element);
    mol->atoms[idx].index = idx;
    mol->atoms[idx].input_order = idx;
    mol->num_atoms++;

    mol->is_canonical = false;
    return idx;
}

int molecule_add_atom_full(molecule_t* mol, const atom_t* atom) {
    if (!mol || !atom) return -1;

    if (molecule_ensure_atom_capacity(mol, mol->num_atoms + 1) < 0) {
        snprintf(mol->error_msg, sizeof(mol->error_msg), "Failed to allocate memory for atom");
        return -1;
    }

    int idx = mol->num_atoms;
    atom_copy(&mol->atoms[idx], atom);
    mol->atoms[idx].index = idx;
    if (mol->atoms[idx].input_order < 0) {
        mol->atoms[idx].input_order = idx;
    }
    mol->num_atoms++;

    mol->is_canonical = false;
    return idx;
}

int molecule_add_bond(molecule_t* mol, int atom1, int atom2, bond_type_t type) {
    if (!mol) return -1;
    if (atom1 < 0 || atom1 >= mol->num_atoms) return -1;
    if (atom2 < 0 || atom2 >= mol->num_atoms) return -1;
    if (atom1 == atom2) return -1;  /* No self-bonds */

    /* Check for existing bond */
    int existing = molecule_get_bond_index(mol, atom1, atom2);
    if (existing >= 0) {
        /* Update bond type if already exists */
        mol->bonds[existing].type = type;
        mol->bonds[existing].aromatic = (type == BOND_AROMATIC || type == BOND_RING_AROMATIC);
        /* If bond is aromatic, mark both atoms as aromatic */
        if (type == BOND_AROMATIC || type == BOND_RING_AROMATIC) {
            mol->atoms[atom1].aromatic = true;
            mol->atoms[atom2].aromatic = true;
        }
        return existing;
    }

    if (molecule_ensure_bond_capacity(mol, mol->num_bonds + 1) < 0) {
        snprintf(mol->error_msg, sizeof(mol->error_msg), "Failed to allocate memory for bond");
        return -1;
    }

    int idx = mol->num_bonds;
    bond_init(&mol->bonds[idx], atom1, atom2, type);
    mol->bonds[idx].index = idx;
    mol->num_bonds++;

    /* If bond is aromatic, mark both atoms as aromatic.
     * This handles explicit aromatic bond notation like C1:C:C:C:C:C:1 */
    if (type == BOND_AROMATIC || type == BOND_RING_AROMATIC) {
        mol->atoms[atom1].aromatic = true;
        mol->atoms[atom2].aromatic = true;
    }

    /* Update atom connectivity */
    atom_add_neighbor(&mol->atoms[atom1], atom2, idx);
    atom_add_neighbor(&mol->atoms[atom2], atom1, idx);

    /* Update bond order sum */
    int order = bond_get_int_order(&mol->bonds[idx]);
    mol->atoms[atom1].total_bond_order += order;
    mol->atoms[atom2].total_bond_order += order;

    mol->is_canonical = false;
    mol->rings_computed = false;
    return idx;
}

atom_t* molecule_get_atom(molecule_t* mol, int idx) {
    if (!mol || idx < 0 || idx >= mol->num_atoms) return NULL;
    return &mol->atoms[idx];
}

const atom_t* molecule_get_atom_const(const molecule_t* mol, int idx) {
    if (!mol || idx < 0 || idx >= mol->num_atoms) return NULL;
    return &mol->atoms[idx];
}

bond_t* molecule_get_bond(molecule_t* mol, int idx) {
    if (!mol || idx < 0 || idx >= mol->num_bonds) return NULL;
    return &mol->bonds[idx];
}

const bond_t* molecule_get_bond_const(const molecule_t* mol, int idx) {
    if (!mol || idx < 0 || idx >= mol->num_bonds) return NULL;
    return &mol->bonds[idx];
}

bond_t* molecule_get_bond_between(molecule_t* mol, int atom1, int atom2) {
    int idx = molecule_get_bond_index(mol, atom1, atom2);
    if (idx < 0) return NULL;
    return &mol->bonds[idx];
}

const bond_t* molecule_get_bond_between_const(const molecule_t* mol, int atom1, int atom2) {
    int idx = molecule_get_bond_index(mol, atom1, atom2);
    if (idx < 0) return NULL;
    return &mol->bonds[idx];
}

int molecule_get_bond_index(const molecule_t* mol, int atom1, int atom2) {
    if (!mol) return -1;

    /* Check atom1's neighbors */
    const atom_t* a1 = molecule_get_atom_const(mol, atom1);
    if (!a1) return -1;

    for (int i = 0; i < a1->num_neighbors; i++) {
        if (a1->neighbors[i] == atom2) {
            return a1->neighbor_bonds[i];
        }
    }

    return -1;
}

void molecule_calc_implicit_h(molecule_t* mol) {
    if (!mol) return;

    for (int i = 0; i < mol->num_atoms; i++) {
        mol->atoms[i].implicit_h_count = atom_calc_implicit_h(&mol->atoms[i], mol);
    }
}

/* Thread-local scratch for fragment finding BFS */
#define FRAG_FINDER_MAX_ATOMS 1024
static __thread int tl_frag_queue[FRAG_FINDER_MAX_ATOMS];

cchem_status_t molecule_find_fragments(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Allocate fragment IDs (only if not already allocated or too small) */
    if (!mol->fragment_ids || mol->num_atoms > mol->atoms_capacity) {
        if (mol->fragment_ids) free(mol->fragment_ids);
        mol->fragment_ids = (int*)calloc(mol->num_atoms, sizeof(int));
        if (!mol->fragment_ids) return CCHEM_ERROR_MEMORY;
    }

    /* Initialize all to -1 (unvisited) */
    for (int i = 0; i < mol->num_atoms; i++) {
        mol->fragment_ids[i] = -1;
    }

    mol->num_fragments = 0;

    /* Use thread-local queue for small molecules (common case) */
    bool use_tl = (mol->num_atoms <= FRAG_FINDER_MAX_ATOMS);
    int* queue;

    if (use_tl) {
        queue = tl_frag_queue;
    } else {
        queue = (int*)calloc(mol->num_atoms, sizeof(int));
        if (!queue) {
            free(mol->fragment_ids);
            mol->fragment_ids = NULL;
            return CCHEM_ERROR_MEMORY;
        }
    }

    for (int start = 0; start < mol->num_atoms; start++) {
        if (mol->fragment_ids[start] >= 0) continue;  /* Already visited */

        /* BFS from this atom */
        int front = 0, back = 0;
        queue[back++] = start;
        mol->fragment_ids[start] = mol->num_fragments;

        while (front < back) {
            int curr = queue[front++];
            atom_t* atom = &mol->atoms[curr];

            for (int i = 0; i < atom->num_neighbors; i++) {
                int neighbor = atom->neighbors[i];
                if (mol->fragment_ids[neighbor] < 0) {
                    mol->fragment_ids[neighbor] = mol->num_fragments;
                    queue[back++] = neighbor;
                }
            }
        }

        mol->num_fragments++;
    }

    if (!use_tl) {
        free(queue);
    }
    return CCHEM_OK;
}

double molecule_calc_weight(const molecule_t* mol) {
    if (!mol) return 0.0;

    double weight = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];

        /* Get atomic mass */
        double mass = element_atomic_mass(atom->element);
        if (atom->isotope > 0) {
            mass = (double)atom->isotope;
        }
        weight += mass;

        /* Add implicit hydrogens */
        weight += atom->implicit_h_count * 1.008;
    }

    return weight;
}

int molecule_calc_charge(const molecule_t* mol) {
    if (!mol) return 0;

    int charge = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        charge += mol->atoms[i].charge;
    }
    return charge;
}

cchem_status_t molecule_validate(const molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Check atom indices */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].index != i) {
            return CCHEM_ERROR_INVALID_INPUT;
        }

        /* Check neighbor validity */
        for (int j = 0; j < mol->atoms[i].num_neighbors; j++) {
            int neighbor = mol->atoms[i].neighbors[j];
            if (neighbor < 0 || neighbor >= mol->num_atoms) {
                return CCHEM_ERROR_INVALID_INPUT;
            }
        }
    }

    /* Check bond validity */
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->atom1 < 0 || bond->atom1 >= mol->num_atoms) {
            return CCHEM_ERROR_INVALID_INPUT;
        }
        if (bond->atom2 < 0 || bond->atom2 >= mol->num_atoms) {
            return CCHEM_ERROR_INVALID_INPUT;
        }
    }

    return CCHEM_OK;
}

bool molecule_is_connected(const molecule_t* mol) {
    if (!mol) return false;
    if (mol->num_atoms <= 1) return true;

    /* Quick check using fragment info if available */
    if (mol->fragment_ids && mol->num_fragments > 0) {
        return mol->num_fragments == 1;
    }

    /* BFS connectivity check */
    bool* visited = (bool*)calloc(mol->num_atoms, sizeof(bool));
    if (!visited) return false;

    int* queue = (int*)calloc(mol->num_atoms, sizeof(int));
    if (!queue) {
        free(visited);
        return false;
    }

    int front = 0, back = 0;
    queue[back++] = 0;
    visited[0] = true;
    int visited_count = 1;

    while (front < back) {
        int curr = queue[front++];
        const atom_t* atom = &mol->atoms[curr];

        for (int i = 0; i < atom->num_neighbors; i++) {
            int neighbor = atom->neighbors[i];
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue[back++] = neighbor;
                visited_count++;
            }
        }
    }

    free(visited);
    free(queue);

    return visited_count == mol->num_atoms;
}

int molecule_num_heavy_atoms(const molecule_t* mol) {
    if (!mol) return 0;

    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            count++;
        }
    }
    return count;
}

molecule_t* molecule_clone(const molecule_t* mol) {
    if (!mol) return NULL;

    molecule_t* clone = molecule_create_with_capacity(mol->num_atoms, mol->num_bonds);
    if (!clone) return NULL;

    /* Copy atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        atom_copy(&clone->atoms[i], &mol->atoms[i]);
    }
    clone->num_atoms = mol->num_atoms;

    /* Copy bonds */
    for (int i = 0; i < mol->num_bonds; i++) {
        bond_copy(&clone->bonds[i], &mol->bonds[i]);
    }
    clone->num_bonds = mol->num_bonds;

    /* Copy rings if computed */
    if (mol->rings_computed && mol->num_rings > 0) {
        if (mol->num_rings > clone->rings_capacity) {
            ring_t* new_rings = (ring_t*)realloc(clone->rings, mol->num_rings * sizeof(ring_t));
            if (new_rings) {
                clone->rings = new_rings;
                clone->rings_capacity = mol->num_rings;
            }
        }
        memcpy(clone->rings, mol->rings, mol->num_rings * sizeof(ring_t));
        clone->num_rings = mol->num_rings;
        clone->rings_computed = true;
    }

    /* Copy fragment info */
    if (mol->fragment_ids) {
        clone->fragment_ids = (int*)calloc(mol->num_atoms, sizeof(int));
        if (clone->fragment_ids) {
            memcpy(clone->fragment_ids, mol->fragment_ids, mol->num_atoms * sizeof(int));
        }
        clone->num_fragments = mol->num_fragments;
    }

    /* Copy canonical order */
    if (mol->canon_order) {
        clone->canon_order = (int*)calloc(mol->num_atoms, sizeof(int));
        if (clone->canon_order) {
            memcpy(clone->canon_order, mol->canon_order, mol->num_atoms * sizeof(int));
        }
        clone->is_canonical = mol->is_canonical;
    }

    /* Copy original SMILES */
    if (mol->original_smiles) {
        clone->original_smiles = strdup(mol->original_smiles);
    }

    clone->total_charge = mol->total_charge;
    clone->molecular_weight = mol->molecular_weight;

    return clone;
}

const char* molecule_get_error(const molecule_t* mol) {
    if (!mol) return "NULL molecule";
    return mol->error_msg;
}

/*
 * Check if an element can be aromatic
 */
static bool element_can_be_aromatic(element_t elem) {
    switch (elem) {
        case ELEM_C:
        case ELEM_N:
        case ELEM_O:
        case ELEM_S:
        case ELEM_P:
            return true;
        default:
            return false;
    }
}

/*
 * Count pi electrons contributed by an atom to an aromatic ring.
 * Uses simple Hückel-like counting:
 * - C with double bond: 1 electron
 * - C with single bonds only (carbocation-like): 0 electrons
 * - N with lone pair: 2 electrons
 * - N in pyridine-like (with double bond): 1 electron
 * - O/S with lone pair: 2 electrons
 */
static int count_pi_electrons(molecule_t* mol, int atom_idx, const ring_t* ring) {
    atom_t* atom = &mol->atoms[atom_idx];

    /* Check if this atom has a double or aromatic bond to another ring atom,
     * or if the atom is already marked aromatic (from fused ring system) */
    bool has_pi_character = atom->aromatic;  /* Already aromatic from another ring */

    if (!has_pi_character) {
        for (int i = 0; i < atom->num_neighbors; i++) {
            int neighbor = atom->neighbors[i];
            int bond_idx = atom->neighbor_bonds[i];

            /* Check if neighbor is in this ring */
            bool neighbor_in_ring = false;
            for (int j = 0; j < ring->size; j++) {
                if (ring->atoms[j] == neighbor) {
                    neighbor_in_ring = true;
                    break;
                }
            }

            if (neighbor_in_ring) {
                bond_type_t btype = mol->bonds[bond_idx].type;
                if (btype == BOND_DOUBLE || btype == BOND_AROMATIC) {
                    has_pi_character = true;
                    break;
                }
            }
        }
    }

    switch (atom->element) {
        case ELEM_C:
            /* Carbon contributes 1 electron if it has pi character */
            return has_pi_character ? 1 : 0;

        case ELEM_N:
            /* Nitrogen: pyridine-like (in conjugation) contributes 1,
               pyrrole-like (with lone pair) contributes 2 */
            if (has_pi_character) {
                return 1;  /* Pyridine-like */
            } else {
                return 2;  /* Pyrrole-like (lone pair) */
            }

        case ELEM_O:
        case ELEM_S:
            /* O/S contribute their lone pair (2 electrons) */
            return 2;

        default:
            return 0;
    }
}

/*
 * Check if a ring is aromatic using Hückel's rule (4n+2 pi electrons)
 */
static bool ring_is_aromatic(molecule_t* mol, const ring_t* ring) {
    /* Check all atoms can be aromatic */
    for (int i = 0; i < ring->size; i++) {
        if (!element_can_be_aromatic(mol->atoms[ring->atoms[i]].element)) {
            return false;
        }
    }

    /*
     * For aromaticity, EVERY atom in the ring must have pi character:
     * - Already marked aromatic, OR
     * - Has a double/aromatic bond to another ring atom, OR
     * - Is a heteroatom (N, O, S) that can donate a lone pair
     *
     * A ring with ANY saturated carbon (sp3) cannot be aromatic.
     */
    for (int i = 0; i < ring->size; i++) {
        int atom_idx = ring->atoms[i];
        atom_t* atom = &mol->atoms[atom_idx];

        /* Heteroatoms can always contribute (lone pair) */
        if (atom->element == ELEM_N || atom->element == ELEM_O || atom->element == ELEM_S) {
            continue;
        }

        /* Already aromatic is fine */
        if (atom->aromatic) {
            continue;
        }

        /* For carbon: must have double/aromatic bond to another ring atom */
        bool has_pi_bond = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int neighbor = atom->neighbors[j];
            int bond_idx = atom->neighbor_bonds[j];

            /* Check if neighbor is in this ring */
            bool neighbor_in_ring = false;
            for (int k = 0; k < ring->size; k++) {
                if (ring->atoms[k] == neighbor) {
                    neighbor_in_ring = true;
                    break;
                }
            }

            if (neighbor_in_ring) {
                bond_type_t btype = mol->bonds[bond_idx].type;
                if (btype == BOND_DOUBLE || btype == BOND_AROMATIC) {
                    has_pi_bond = true;
                    break;
                }
            }
        }

        if (!has_pi_bond) {
            return false;  /* This carbon is sp3, ring cannot be aromatic */
        }
    }

    /* Count pi electrons */
    int pi_electrons = 0;
    for (int i = 0; i < ring->size; i++) {
        pi_electrons += count_pi_electrons(mol, ring->atoms[i], ring);
    }

    /* Check Hückel's rule: 4n+2 for n=0,1,2,3... gives 2,6,10,14... */
    /* Common aromatic rings: benzene (6), naphthalene (10), pyridine (6), furan (6), etc. */
    if (pi_electrons >= 2 && (pi_electrons - 2) % 4 == 0) {
        return true;
    }

    return false;
}

cchem_status_t molecule_perceive_aromaticity(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Need rings to perceive aromaticity */
    if (!mol->rings_computed) {
        cchem_status_t status = molecule_find_rings(mol);
        if (status != CCHEM_OK) return status;
    }

    /* For fused ring systems, we need to iterate.
     * When one ring becomes aromatic, atoms shared with other rings
     * gain aromatic character, which may make those rings aromatic too.
     */
    bool changed = true;
    int max_iterations = mol->num_rings + 1;  /* Prevent infinite loop */

    while (changed && max_iterations-- > 0) {
        changed = false;

        for (int r = 0; r < mol->num_rings; r++) {
            ring_t* ring = &mol->rings[r];

            /* For already aromatic rings (from aromatic SMILES input),
             * we still need to compute pi_electrons for each atom */
            if (ring->aromatic) {
                for (int i = 0; i < ring->size; i++) {
                    int atom_idx = ring->atoms[i];
                    if (mol->atoms[atom_idx].pi_electrons == 0) {
                        mol->atoms[atom_idx].pi_electrons = count_pi_electrons(mol, atom_idx, ring);
                    }
                }
                continue;
            }

            if (ring_is_aromatic(mol, ring)) {
                ring->aromatic = true;
                changed = true;

                /* First pass: compute pi electrons BEFORE marking atoms as aromatic
                 * This is critical because count_pi_electrons uses atom->aromatic to determine
                 * has_pi_character, which affects whether N is pyridine-type (1) or pyrrole-type (2)
                 */
                for (int i = 0; i < ring->size; i++) {
                    int atom_idx = ring->atoms[i];
                    if (mol->atoms[atom_idx].pi_electrons == 0) {
                        mol->atoms[atom_idx].pi_electrons = count_pi_electrons(mol, atom_idx, ring);
                    }
                }

                /* Second pass: mark all atoms in ring as aromatic */
                for (int i = 0; i < ring->size; i++) {
                    int atom_idx = ring->atoms[i];
                    mol->atoms[atom_idx].aromatic = true;
                }

                /* Mark all bonds in ring as aromatic */
                for (int i = 0; i < ring->size; i++) {
                    int atom1 = ring->atoms[i];
                    int atom2 = ring->atoms[(i + 1) % ring->size];
                    int bond_idx = molecule_get_bond_index(mol, atom1, atom2);
                    if (bond_idx >= 0) {
                        mol->bonds[bond_idx].type = BOND_AROMATIC;
                        mol->bonds[bond_idx].in_ring = true;
                    }
                }
            }
        }
    }

    return CCHEM_OK;
}
