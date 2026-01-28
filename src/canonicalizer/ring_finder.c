/**
 * @file ring_finder.c
 * @brief Ring detection algorithms (SSSR)
 */

#include "cchem/compat.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/ring_finder.h"

/* Thread-local scratch memory to avoid malloc/free per molecule */
#define RING_FINDER_MAX_ATOMS 1024  /* Max atoms for thread-local buffers */

typedef struct {
    bool visited[RING_FINDER_MAX_ATOMS];
    int distance[RING_FINDER_MAX_ATOMS];
    int ring_path[MAX_RING_SIZE];
    int queue_data[RING_FINDER_MAX_ATOMS * 2];
    int queue_parent[RING_FINDER_MAX_ATOMS * 2];
    int queue_parent_bond[RING_FINDER_MAX_ATOMS * 2];
    bool initialized;
} ring_finder_scratch_t;

static __thread ring_finder_scratch_t tl_scratch = {0};

/* Queue for BFS (uses thread-local storage or heap for large molecules) */
typedef struct {
    int* data;
    int* parent;
    int* parent_bond;
    int front;
    int back;
    int capacity;
    bool heap_allocated;  /* True if using heap memory */
} bfs_queue_t;

static void bfs_queue_init_tl(bfs_queue_t* q, int capacity) {
    q->data = tl_scratch.queue_data;
    q->parent = tl_scratch.queue_parent;
    q->parent_bond = tl_scratch.queue_parent_bond;
    q->capacity = RING_FINDER_MAX_ATOMS * 2;
    q->front = 0;
    q->back = 0;
    q->heap_allocated = false;
    (void)capacity;
}

static bfs_queue_t* bfs_queue_create(int capacity) {
    bfs_queue_t* q = (bfs_queue_t*)calloc(1, sizeof(bfs_queue_t));
    if (!q) return NULL;

    q->data = (int*)calloc(capacity, sizeof(int));
    q->parent = (int*)calloc(capacity, sizeof(int));
    q->parent_bond = (int*)calloc(capacity, sizeof(int));

    if (!q->data || !q->parent || !q->parent_bond) {
        if (q->data) free(q->data);
        if (q->parent) free(q->parent);
        if (q->parent_bond) free(q->parent_bond);
        free(q);
        return NULL;
    }

    q->capacity = capacity;
    q->front = 0;
    q->back = 0;
    q->heap_allocated = true;
    return q;
}

static void bfs_queue_free(bfs_queue_t* q) {
    if (q && q->heap_allocated) {
        if (q->data) free(q->data);
        if (q->parent) free(q->parent);
        if (q->parent_bond) free(q->parent_bond);
        free(q);
    }
}

static void bfs_queue_push(bfs_queue_t* q, int atom, int parent, int parent_bond) {
    if (q->back < q->capacity) {
        q->data[q->back] = atom;
        q->parent[q->back] = parent;
        q->parent_bond[q->back] = parent_bond;
        q->back++;
    }
}

static bool bfs_queue_empty(bfs_queue_t* q) {
    return q->front >= q->back;
}

/* Ensure ring capacity */
static cchem_status_t molecule_ensure_ring_capacity(molecule_t* mol) {
    if (mol->num_rings < mol->rings_capacity) return CCHEM_OK;

    int new_cap = mol->rings_capacity * 2;
    ring_t* new_rings = (ring_t*)realloc(mol->rings, new_cap * sizeof(ring_t));
    if (!new_rings) return CCHEM_ERROR_MEMORY;

    mol->rings = new_rings;
    mol->rings_capacity = new_cap;
    return CCHEM_OK;
}

/* Add ring to molecule */
static cchem_status_t molecule_add_ring(molecule_t* mol, int* atoms, int size) {
    if (molecule_ensure_ring_capacity(mol) != CCHEM_OK) {
        return CCHEM_ERROR_MEMORY;
    }

    ring_t* ring = &mol->rings[mol->num_rings];
    memset(ring, 0, sizeof(ring_t));

    ring->size = size;
    ring->index = mol->num_rings;

    for (int i = 0; i < size; i++) {
        ring->atoms[i] = atoms[i];

        /* Update atom ring count */
        mol->atoms[atoms[i]].ring_count++;
    }

    /* Find bonds in ring */
    ring->aromatic = true;  /* Assume aromatic, check later */
    for (int i = 0; i < size; i++) {
        int a1 = atoms[i];
        int a2 = atoms[(i + 1) % size];
        int bond_idx = molecule_get_bond_index(mol, a1, a2);

        if (bond_idx >= 0) {
            ring->bonds[i] = bond_idx;
            mol->bonds[bond_idx].in_ring = true;

            /* Check aromaticity */
            if (!mol->bonds[bond_idx].aromatic) {
                ring->aromatic = false;
            }

            /* Add ring to bond */
            bond_t* bond = &mol->bonds[bond_idx];
            if (bond->num_rings < MAX_RINGS) {
                bond->ring_ids[bond->num_rings++] = mol->num_rings;
            }
        } else {
            ring->bonds[i] = -1;
        }
    }

    mol->num_rings++;
    return CCHEM_OK;
}

cchem_status_t ring_finder_sssr(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Reset ring data */
    mol->num_rings = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        mol->atoms[i].ring_count = 0;
    }
    for (int i = 0; i < mol->num_bonds; i++) {
        mol->bonds[i].in_ring = false;
        mol->bonds[i].num_rings = 0;
    }

    /* SSSR size = num_bonds - num_atoms + num_fragments */
    molecule_find_fragments(mol);
    int expected_rings = mol->num_bonds - mol->num_atoms + mol->num_fragments;

    if (expected_rings <= 0) {
        mol->rings_computed = true;
        return CCHEM_OK;
    }

    /* Use thread-local scratch memory for small molecules (common case) */
    bool use_tl = (mol->num_atoms <= RING_FINDER_MAX_ATOMS);

    bool* visited;
    int* distance;
    int* ring_path;
    bfs_queue_t queue_stack;
    bfs_queue_t* queue;

    if (use_tl) {
        /* Use pre-allocated thread-local buffers (zero malloc) */
        visited = tl_scratch.visited;
        distance = tl_scratch.distance;
        ring_path = tl_scratch.ring_path;
        bfs_queue_init_tl(&queue_stack, mol->num_atoms * 2);
        queue = &queue_stack;
    } else {
        /* Fall back to heap for large molecules */
        visited = (bool*)calloc(mol->num_atoms, sizeof(bool));
        distance = (int*)calloc(mol->num_atoms, sizeof(int));
        ring_path = (int*)calloc(MAX_RING_SIZE, sizeof(int));
        queue = bfs_queue_create(mol->num_atoms * 2);

        if (!visited || !distance || !ring_path || !queue) {
            if (visited) free(visited);
            if (distance) free(distance);
            if (ring_path) free(ring_path);
            if (queue) bfs_queue_free(queue);
            return CCHEM_ERROR_MEMORY;
        }
    }

    /* Process each bond that could be part of a ring */
    for (int bond_idx = 0; bond_idx < mol->num_bonds && mol->num_rings < expected_rings; bond_idx++) {
        bond_t* bond = &mol->bonds[bond_idx];

        /* Skip if already in a ring */
        if (bond->in_ring) continue;

        int start = bond->atom1;
        int end = bond->atom2;

        /* BFS to find shortest path excluding this bond */
        memset(visited, 0, mol->num_atoms * sizeof(bool));
        for (int i = 0; i < mol->num_atoms; i++) {
            distance[i] = -1;
        }

        queue->front = 0;
        queue->back = 0;

        bfs_queue_push(queue, start, -1, -1);
        visited[start] = true;
        distance[start] = 0;

        bool found_ring = false;

        while (!bfs_queue_empty(queue) && !found_ring) {
            int curr = queue->data[queue->front];
            int curr_dist = distance[curr];
            queue->front++;

            atom_t* atom = &mol->atoms[curr];

            for (int i = 0; i < atom->num_neighbors; i++) {
                int neighbor = atom->neighbors[i];
                int neighbor_bond = atom->neighbor_bonds[i];

                /* Skip the bond we're testing */
                if (neighbor_bond == bond_idx) continue;

                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    distance[neighbor] = curr_dist + 1;
                    bfs_queue_push(queue, neighbor, curr, neighbor_bond);

                    /* Found the end - we have a ring */
                    if (neighbor == end) {
                        /* Reconstruct ring path */
                        int ring_size = 0;
                        int temp = neighbor;

                        while (temp != -1 && ring_size < MAX_RING_SIZE) {
                            ring_path[ring_size++] = temp;

                            int parent = -1;
                            for (int j = 0; j < queue->back; j++) {
                                if (queue->data[j] == temp) {
                                    parent = queue->parent[j];
                                    break;
                                }
                            }
                            temp = parent;
                        }

                        /* Add ring to molecule */
                        if (ring_size >= 3) {
                            molecule_add_ring(mol, ring_path, ring_size);
                        }

                        found_ring = true;
                        break;
                    }
                }
            }
        }
    }

    /* Only free heap-allocated memory */
    if (!use_tl) {
        free(visited);
        free(distance);
        free(ring_path);
        bfs_queue_free(queue);
    }

    mol->rings_computed = true;
    return CCHEM_OK;
}

cchem_status_t molecule_find_rings(molecule_t* mol) {
    return ring_finder_sssr(mol);
}

int ring_finder_atom_ring_count(const molecule_t* mol, int atom_idx) {
    if (!mol || atom_idx < 0 || atom_idx >= mol->num_atoms) return 0;
    return mol->atoms[atom_idx].ring_count;
}

int ring_finder_bond_ring_count(const molecule_t* mol, int bond_idx) {
    if (!mol || bond_idx < 0 || bond_idx >= mol->num_bonds) return 0;
    return mol->bonds[bond_idx].num_rings;
}

/* Count pi electrons for aromaticity (Huckel's rule: 4n+2) */
int ring_finder_count_pi_electrons(const molecule_t* mol, const ring_t* ring) {
    if (!mol || !ring) return 0;

    int pi_electrons = 0;

    for (int i = 0; i < ring->size; i++) {
        int atom_idx = ring->atoms[i];
        const atom_t* atom = &mol->atoms[atom_idx];

        /* Carbon: 1 pi electron if part of aromatic system */
        if (atom->element == ELEM_C) {
            if (atom->aromatic) {
                pi_electrons += 1;
            }
        }
        /* Nitrogen: 1 or 2 pi electrons depending on lone pair */
        else if (atom->element == ELEM_N) {
            if (atom->charge == 1) {
                pi_electrons += 1;  /* Pyridinium-like */
            } else if (atom->num_neighbors == 2) {
                pi_electrons += 2;  /* Pyrrole-like (lone pair) */
            } else {
                pi_electrons += 1;  /* Pyridine-like */
            }
        }
        /* Oxygen: 2 pi electrons (lone pair) */
        else if (atom->element == ELEM_O) {
            pi_electrons += 2;
        }
        /* Sulfur: 2 pi electrons (lone pair) */
        else if (atom->element == ELEM_S) {
            pi_electrons += 2;
        }
    }

    return pi_electrons;
}

bool ring_finder_is_aromatic_ring(const molecule_t* mol, const ring_t* ring) {
    if (!mol || !ring) return false;

    /* Check Huckel's rule: 4n+2 pi electrons */
    int pi = ring_finder_count_pi_electrons(mol, ring);

    /* 4n+2 for n = 0, 1, 2, 3, ... gives 2, 6, 10, 14, ... */
    if (pi < 2) return false;

    int remainder = (pi - 2) % 4;
    return (remainder == 0);
}

cchem_status_t ring_finder_perceive_aromaticity(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* First find rings if not done */
    if (!mol->rings_computed) {
        cchem_status_t status = ring_finder_sssr(mol);
        if (status != CCHEM_OK) return status;
    }

    /* Check each ring for aromaticity */
    for (int r = 0; r < mol->num_rings; r++) {
        ring_t* ring = &mol->rings[r];
        ring->aromatic = ring_finder_is_aromatic_ring(mol, ring);
    }

    /* Update atom and bond aromaticity based on rings */
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];

        if (ring->aromatic) {
            /* Mark all atoms in ring as aromatic */
            for (int i = 0; i < ring->size; i++) {
                mol->atoms[ring->atoms[i]].aromatic = true;
            }

            /* Mark all bonds in ring as aromatic */
            for (int i = 0; i < ring->size; i++) {
                if (ring->bonds[i] >= 0) {
                    mol->bonds[ring->bonds[i]].aromatic = true;
                    mol->bonds[ring->bonds[i]].type = BOND_AROMATIC;
                }
            }
        }
    }

    return CCHEM_OK;
}
