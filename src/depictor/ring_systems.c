/**
 * @file ring_systems.c
 * @brief Ring system detection and classification for 2D layout
 *
 * Groups rings into fused/spiro/bridged systems and provides
 * ring traversal utilities for placement algorithms.
 */

#include "cchem/depictor/layout2d.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* ============== Union-Find for Ring Grouping ============== */

static int uf_find(int* parent, int i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]]; /* Path compression */
        i = parent[i];
    }
    return i;
}

static void uf_union(int* parent, int* rank, int a, int b) {
    int ra = uf_find(parent, a);
    int rb = uf_find(parent, b);
    if (ra == rb) return;

    /* Union by rank */
    if (rank[ra] < rank[rb]) {
        parent[ra] = rb;
    } else if (rank[ra] > rank[rb]) {
        parent[rb] = ra;
    } else {
        parent[rb] = ra;
        rank[ra]++;
    }
}

/* ============== Ring Connectivity ============== */

int layout_count_shared_atoms(const ring_t* r1, const ring_t* r2,
                              int* shared1_idx, int* shared2_idx, int max_shared) {
    int count = 0;
    for (int i = 0; i < r1->size && count < max_shared; i++) {
        for (int j = 0; j < r2->size; j++) {
            if (r1->atoms[i] == r2->atoms[j]) {
                if (shared1_idx) shared1_idx[count] = i;
                if (shared2_idx) shared2_idx[count] = j;
                count++;
                break;
            }
        }
    }
    return count;
}

bool layout_rings_are_fused(const ring_t* r1, const ring_t* r2) {
    return layout_count_shared_atoms(r1, r2, NULL, NULL, MAX_SHARED_ATOMS) >= 2;
}

/* ============== Ring Ordering ============== */

void layout_get_ring_order(const ring_t* ring, const molecule_t* mol,
                           int start_atom, int* ordered, int* num_ordered) {
    int n = ring->size;
    bool* in_ring = calloc(mol->num_atoms, sizeof(bool));
    bool* visited = calloc(mol->num_atoms, sizeof(bool));

    for (int i = 0; i < n; i++) {
        in_ring[ring->atoms[i]] = true;
    }

    *num_ordered = 0;
    int current = start_atom;

    while (*num_ordered < n) {
        ordered[(*num_ordered)++] = current;
        visited[current] = true;

        /* Find next ring atom that's bonded to current and not yet visited */
        const atom_t* atom = &mol->atoms[current];
        int next = -1;
        for (int i = 0; i < atom->num_neighbors; i++) {
            int nb = atom->neighbors[i];
            if (in_ring[nb] && !visited[nb]) {
                next = nb;
                break;
            }
        }

        if (next < 0) break;
        current = next;
    }

    free(in_ring);
    free(visited);
}

/* ============== Ring System Classification ============== */

static ring_system_type_t classify_ring_system(const molecule_t* mol,
                                               const int* ring_indices,
                                               int num_rings) {
    if (num_rings == 1) {
        return RING_SYSTEM_ISOLATED;
    }

    /* Check for spiro: any pair shares exactly 1 atom */
    bool has_spiro = false;
    bool has_fused = false;
    bool has_bridged = false;

    for (int i = 0; i < num_rings && !has_bridged; i++) {
        for (int j = i + 1; j < num_rings; j++) {
            const ring_t* r1 = &mol->rings[ring_indices[i]];
            const ring_t* r2 = &mol->rings[ring_indices[j]];
            int shared = layout_count_shared_atoms(r1, r2, NULL, NULL, MAX_SHARED_ATOMS);

            if (shared == 1) {
                has_spiro = true;
            } else if (shared == 2) {
                has_fused = true;
            } else if (shared > 2) {
                has_bridged = true;
                break;
            }
        }
    }

    if (has_bridged) return RING_SYSTEM_BRIDGED;
    if (has_fused) return RING_SYSTEM_FUSED;
    if (has_spiro) return RING_SYSTEM_SPIRO;
    return RING_SYSTEM_ISOLATED;
}

static void collect_shared_atoms(const molecule_t* mol,
                                 const int* ring_indices, int num_rings,
                                 int** shared_atoms, int* num_shared) {
    /* Use a bitmap to track which atoms are shared */
    bool* is_shared = calloc(mol->num_atoms, sizeof(bool));
    int* counts = calloc(mol->num_atoms, sizeof(int));

    /* Count occurrences of each atom across rings */
    for (int r = 0; r < num_rings; r++) {
        const ring_t* ring = &mol->rings[ring_indices[r]];
        for (int i = 0; i < ring->size; i++) {
            counts[ring->atoms[i]]++;
        }
    }

    /* Mark atoms that appear in multiple rings */
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (counts[i] > 1) {
            is_shared[i] = true;
            count++;
        }
    }

    /* Allocate and fill shared atoms array */
    *num_shared = count;
    if (count > 0) {
        *shared_atoms = malloc(count * sizeof(int));
        int idx = 0;
        for (int i = 0; i < mol->num_atoms; i++) {
            if (is_shared[i]) {
                (*shared_atoms)[idx++] = i;
            }
        }
    } else {
        *shared_atoms = NULL;
    }

    free(is_shared);
    free(counts);
}

static void collect_all_atoms(const molecule_t* mol,
                              const int* ring_indices, int num_rings,
                              int** all_atoms, int* num_atoms) {
    bool* in_system = calloc(mol->num_atoms, sizeof(bool));

    int count = 0;
    for (int r = 0; r < num_rings; r++) {
        const ring_t* ring = &mol->rings[ring_indices[r]];
        for (int i = 0; i < ring->size; i++) {
            if (!in_system[ring->atoms[i]]) {
                in_system[ring->atoms[i]] = true;
                count++;
            }
        }
    }

    *num_atoms = count;
    *all_atoms = malloc(count * sizeof(int));
    int idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (in_system[i]) {
            (*all_atoms)[idx++] = i;
        }
    }

    free(in_system);
}

/* ============== Main Ring System Detection ============== */

int layout_detect_ring_systems(layout_context_t* ctx) {
    if (!ctx || !ctx->mol) return -1;

    const molecule_t* mol = ctx->mol;

    /* Handle case with no rings */
    if (!mol->rings_computed || mol->num_rings == 0) {
        ctx->ring_systems = NULL;
        ctx->num_ring_systems = 0;
        return 0;
    }

    int nr = mol->num_rings;

    /* Initialize union-find */
    int* parent = malloc(nr * sizeof(int));
    int* rank = calloc(nr, sizeof(int));
    for (int i = 0; i < nr; i++) {
        parent[i] = i;
    }

    /* Group only FUSED rings (sharing 2+ atoms) into systems.
     * Spiro junctions (1 shared atom) are handled separately during placement. */
    for (int i = 0; i < nr; i++) {
        for (int j = i + 1; j < nr; j++) {
            if (layout_rings_are_fused(&mol->rings[i], &mol->rings[j])) {
                uf_union(parent, rank, i, j);
            }
        }
    }

    /* Count distinct ring systems */
    int* root_count = calloc(nr, sizeof(int));
    for (int i = 0; i < nr; i++) {
        root_count[uf_find(parent, i)]++;
    }

    int num_systems = 0;
    for (int i = 0; i < nr; i++) {
        if (root_count[i] > 0) num_systems++;
    }

    /* Allocate ring systems */
    ctx->ring_systems = calloc(num_systems, sizeof(ring_system_t));
    ctx->num_ring_systems = num_systems;

    /* Build ring systems */
    int sys_idx = 0;
    for (int root = 0; root < nr; root++) {
        if (root_count[root] == 0) continue;

        ring_system_t* sys = &ctx->ring_systems[sys_idx];

        /* Collect rings in this system */
        sys->num_rings = root_count[root];
        sys->ring_indices = malloc(sys->num_rings * sizeof(int));

        int idx = 0;
        for (int r = 0; r < nr; r++) {
            if (uf_find(parent, r) == root) {
                sys->ring_indices[idx++] = r;
            }
        }

        /* Classify the system */
        sys->type = classify_ring_system(mol, sys->ring_indices, sys->num_rings);

        /* Collect shared atoms */
        collect_shared_atoms(mol, sys->ring_indices, sys->num_rings,
                            &sys->shared_atoms, &sys->num_shared_atoms);

        /* Collect all atoms */
        collect_all_atoms(mol, sys->ring_indices, sys->num_rings,
                         &sys->all_atoms, &sys->num_atoms);

        sys->placed = false;

        /* Update atom-to-ring-system mapping */
        for (int i = 0; i < sys->num_atoms; i++) {
            ctx->atom_to_ring_system[sys->all_atoms[i]] = sys_idx;
        }

        sys_idx++;
    }

    free(parent);
    free(rank);
    free(root_count);

    return 0;
}
