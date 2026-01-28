/**
 * @file topology.c
 * @brief Fast topological descriptors for LogD prediction
 *
 * 30 original descriptors based on molecular graph topology:
 * - Kier-Hall connectivity indices (Chi0-Chi4)
 * - Zagreb indices (first, second, hyper)
 * - Wiener-type indices
 * - Balaban J index
 * - Path-based descriptors
 * - Neighborhood complexity measures
 * - Ionization-focused topology
 *
 * All descriptors are O(n) or O(n*m) for fast computation.
 */

#include "cchem/compat.h"
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* Maximum supported atoms for topology descriptors */
#define TOPO_MAX_ATOMS 512

/* Helper macro: Return default value for molecules too large */
#define CHECK_MOL_SIZE(mol, value, default_val) \
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { \
        value->d = default_val; \
        return CCHEM_OK; \
    }

/* ============================================================================
 * Helper: Get heavy atom degree (number of non-H neighbors)
 * ============================================================================ */

static int get_heavy_degree(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int degree = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            degree++;
        }
    }
    return degree;
}

/* ============================================================================
 * Helper: Get valence delta (for Kier-Hall indices)
 * delta_v = (Z_v - H) / (Z - Z_v - 1)
 * ============================================================================ */

static double get_valence_delta(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    if (atom->element == ELEM_H) return 0.0;

    int z_v = 0;  /* Valence electrons */
    int z = 0;    /* Atomic number */
    switch (atom->element) {
        case ELEM_C:  z_v = 4; z = 6; break;
        case ELEM_N:  z_v = 5; z = 7; break;
        case ELEM_O:  z_v = 6; z = 8; break;
        case ELEM_S:  z_v = 6; z = 16; break;
        case ELEM_P:  z_v = 5; z = 15; break;
        case ELEM_F:  z_v = 7; z = 9; break;
        case ELEM_Cl: z_v = 7; z = 17; break;
        case ELEM_Br: z_v = 7; z = 35; break;
        case ELEM_I:  z_v = 7; z = 53; break;
        case ELEM_Si: z_v = 4; z = 14; break;
        case ELEM_B:  z_v = 3; z = 5; break;
        default: z_v = 4; z = 6; break;
    }

    int h_count = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
    }

    int denom = z - z_v - 1;
    if (denom <= 0) denom = 1;
    return (double)(z_v - h_count) / denom;
}

/* ============================================================================
 * Helper: Is atom polar (N, O, S, P)
 * ============================================================================ */

static bool is_polar_atom(element_t elem) {
    return elem == ELEM_N || elem == ELEM_O || elem == ELEM_S || elem == ELEM_P;
}

/* ============================================================================
 * Helper: Is atom ionizable (basic N or acidic O)
 * ============================================================================ */

static bool is_ionizable(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    if (atom->element == ELEM_H) return false;

    /* Basic nitrogen (has lone pair, not in amide) */
    if (atom->element == ELEM_N) {
        int h_count = atom->implicit_h_count;
        for (int i = 0; i < atom->num_neighbors; i++) {
            if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
        }
        /* Amines with H */
        if (h_count > 0 && !atom->aromatic) return true;
        /* Aromatic N can be basic */
        if (atom->aromatic) return true;
    }

    /* Acidic oxygen (in COOH, phenol, etc) */
    if (atom->element == ELEM_O) {
        int h_count = atom->implicit_h_count;
        for (int i = 0; i < atom->num_neighbors; i++) {
            if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
        }
        if (h_count > 0) return true;  /* -OH groups */
    }

    return false;
}

/* ============================================================================
 * Helper: Is atom hydrophobic
 * ============================================================================ */

static bool is_hydrophobic(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    if (atom->element == ELEM_H) return false;
    if (atom->element == ELEM_C) return true;
    if (atom->element == ELEM_F || atom->element == ELEM_Cl ||
        atom->element == ELEM_Br || atom->element == ELEM_I) return true;
    return false;
}

/* ============================================================================
 * Kier-Hall Connectivity Indices
 * Chi_n = sum over subgraphs of 1/sqrt(product of deltas)
 * ============================================================================ */

/* Chi0: Zero-order connectivity (sum of 1/sqrt(delta) for each atom) */
static cchem_status_t desc_chi0(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double chi0 = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int deg = get_heavy_degree(mol, i);
        if (deg > 0) {
            chi0 += 1.0 / sqrt((double)deg);
        }
    }
    value->d = chi0;
    return CCHEM_OK;
}

/* Chi1: First-order connectivity (sum over bonds of 1/sqrt(d1*d2)) */
static cchem_status_t desc_chi1(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double chi1 = 0.0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        int d1 = get_heavy_degree(mol, bond->atom1);
        int d2 = get_heavy_degree(mol, bond->atom2);
        if (d1 > 0 && d2 > 0) {
            chi1 += 1.0 / sqrt((double)(d1 * d2));
        }
    }
    value->d = chi1;
    return CCHEM_OK;
}

/* Chi2Path: Second-order path connectivity (2-bond paths) */
static cchem_status_t desc_chi2_path(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double chi2 = 0.0;
    /* For each atom as middle of a 2-path */
    for (int mid = 0; mid < mol->num_atoms; mid++) {
        if (mol->atoms[mid].element == ELEM_H) continue;
        const atom_t* atom = &mol->atoms[mid];
        int d_mid = get_heavy_degree(mol, mid);
        if (d_mid < 2) continue;

        /* Find all pairs of heavy neighbors */
        int heavy_neighbors[16];
        int n_heavy = 0;
        for (int i = 0; i < atom->num_neighbors && n_heavy < 16; i++) {
            int nb = atom->neighbors[i];
            if (mol->atoms[nb].element != ELEM_H) {
                heavy_neighbors[n_heavy++] = nb;
            }
        }

        /* Sum over pairs */
        for (int i = 0; i < n_heavy; i++) {
            for (int j = i + 1; j < n_heavy; j++) {
                int d1 = get_heavy_degree(mol, heavy_neighbors[i]);
                int d2 = get_heavy_degree(mol, heavy_neighbors[j]);
                if (d1 > 0 && d2 > 0) {
                    chi2 += 1.0 / sqrt((double)(d1 * d_mid * d2));
                }
            }
        }
    }
    value->d = chi2;
    return CCHEM_OK;
}

/* Chi3Path: Third-order path connectivity (3-bond paths) */
static cchem_status_t desc_chi3_path(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double chi3 = 0.0;
    /* For each bond as the middle of a 3-path */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int a1 = bond->atom1, a2 = bond->atom2;
        if (mol->atoms[a1].element == ELEM_H || mol->atoms[a2].element == ELEM_H) continue;

        int d1 = get_heavy_degree(mol, a1);
        int d2 = get_heavy_degree(mol, a2);

        /* Find neighbors of a1 (not a2) and neighbors of a2 (not a1) */
        for (int i = 0; i < mol->atoms[a1].num_neighbors; i++) {
            int n1 = mol->atoms[a1].neighbors[i];
            if (n1 == a2 || mol->atoms[n1].element == ELEM_H) continue;
            int d_n1 = get_heavy_degree(mol, n1);

            for (int j = 0; j < mol->atoms[a2].num_neighbors; j++) {
                int n2 = mol->atoms[a2].neighbors[j];
                if (n2 == a1 || n2 == n1 || mol->atoms[n2].element == ELEM_H) continue;
                int d_n2 = get_heavy_degree(mol, n2);

                if (d_n1 > 0 && d1 > 0 && d2 > 0 && d_n2 > 0) {
                    chi3 += 1.0 / sqrt((double)(d_n1 * d1 * d2 * d_n2));
                }
            }
        }
    }
    value->d = chi3;
    return CCHEM_OK;
}

/* Chi3Cluster: Third-order cluster (star with 3 branches) */
static cchem_status_t desc_chi3_cluster(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double chi3c = 0.0;
    /* For each atom with degree >= 3 */
    for (int c = 0; c < mol->num_atoms; c++) {
        if (mol->atoms[c].element == ELEM_H) continue;
        const atom_t* center = &mol->atoms[c];
        int d_c = get_heavy_degree(mol, c);
        if (d_c < 3) continue;

        int heavy_neighbors[16];
        int n_heavy = 0;
        for (int i = 0; i < center->num_neighbors && n_heavy < 16; i++) {
            int nb = center->neighbors[i];
            if (mol->atoms[nb].element != ELEM_H) {
                heavy_neighbors[n_heavy++] = nb;
            }
        }

        /* Sum over all 3-combinations */
        for (int i = 0; i < n_heavy; i++) {
            for (int j = i + 1; j < n_heavy; j++) {
                for (int k = j + 1; k < n_heavy; k++) {
                    int d1 = get_heavy_degree(mol, heavy_neighbors[i]);
                    int d2 = get_heavy_degree(mol, heavy_neighbors[j]);
                    int d3 = get_heavy_degree(mol, heavy_neighbors[k]);
                    if (d1 > 0 && d2 > 0 && d3 > 0) {
                        chi3c += 1.0 / sqrt((double)(d1 * d2 * d3 * d_c));
                    }
                }
            }
        }
    }
    value->d = chi3c;
    return CCHEM_OK;
}

/* Chi0v: Zero-order valence connectivity */
static cchem_status_t desc_chi0v(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double chi0v = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        double delta_v = get_valence_delta(mol, i);
        if (delta_v > 0) {
            chi0v += 1.0 / sqrt(delta_v);
        }
    }
    value->d = chi0v;
    return CCHEM_OK;
}

/* Chi1v: First-order valence connectivity */
static cchem_status_t desc_chi1v(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double chi1v = 0.0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        double dv1 = get_valence_delta(mol, bond->atom1);
        double dv2 = get_valence_delta(mol, bond->atom2);
        if (dv1 > 0 && dv2 > 0) {
            chi1v += 1.0 / sqrt(dv1 * dv2);
        }
    }
    value->d = chi1v;
    return CCHEM_OK;
}

/* ============================================================================
 * Zagreb Indices
 * ============================================================================ */

/* Zagreb1: First Zagreb index = sum(degree^2) */
static cchem_status_t desc_zagreb1(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t z1 = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int d = get_heavy_degree(mol, i);
        z1 += d * d;
    }
    value->i = z1;
    return CCHEM_OK;
}

/* Zagreb2: Second Zagreb index = sum(d1*d2) over bonds */
static cchem_status_t desc_zagreb2(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t z2 = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        int d1 = get_heavy_degree(mol, bond->atom1);
        int d2 = get_heavy_degree(mol, bond->atom2);
        z2 += d1 * d2;
    }
    value->i = z2;
    return CCHEM_OK;
}

/* HyperZagreb: sum((d1+d2)^2) over bonds */
static cchem_status_t desc_hyper_zagreb(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t hz = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        int d1 = get_heavy_degree(mol, bond->atom1);
        int d2 = get_heavy_degree(mol, bond->atom2);
        int sum = d1 + d2;
        hz += sum * sum;
    }
    value->i = hz;
    return CCHEM_OK;
}

/* ============================================================================
 * Wiener-type Indices (using BFS for shortest paths)
 * ============================================================================ */

/* Helper: Compute sum of distances from one atom using BFS */
static int bfs_distance_sum(const molecule_t* mol, int start, int* heavy_map, int n_heavy __attribute__((unused))) {
    if (heavy_map[start] < 0) return 0;

    int dist[TOPO_MAX_ATOMS];
    int queue[TOPO_MAX_ATOMS];
    memset(dist, -1, sizeof(dist));

    int sum = 0;
    int head = 0, tail = 0;
    queue[tail++] = start;
    dist[start] = 0;

    while (head < tail) {
        int curr = queue[head++];
        const atom_t* atom = &mol->atoms[curr];

        for (int i = 0; i < atom->num_neighbors; i++) {
            int nb = atom->neighbors[i];
            if (mol->atoms[nb].element == ELEM_H) continue;
            if (dist[nb] < 0) {
                dist[nb] = dist[curr] + 1;
                queue[tail++] = nb;
                sum += dist[nb];
            }
        }
    }
    return sum;
}

/* WienerIndex: Sum of all shortest path distances */
static cchem_status_t desc_wiener(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->i = 0; return CCHEM_OK; }

    /* Map heavy atoms */
    int heavy_map[TOPO_MAX_ATOMS];
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms && i < TOPO_MAX_ATOMS; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_map[i] = n_heavy++;
        } else {
            heavy_map[i] = -1;
        }
    }

    int64_t wiener = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        wiener += bfs_distance_sum(mol, i, heavy_map, n_heavy);
    }
    /* Each pair counted twice */
    value->i = wiener / 2;
    return CCHEM_OK;
}

/* MeanPathLength: Average shortest path length */
static cchem_status_t desc_mean_path_length(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->d = 0.0; return CCHEM_OK; }

    int heavy_map[TOPO_MAX_ATOMS];
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms && i < TOPO_MAX_ATOMS; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_map[i] = n_heavy++;
        } else {
            heavy_map[i] = -1;
        }
    }

    if (n_heavy < 2) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    int64_t total = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        total += bfs_distance_sum(mol, i, heavy_map, n_heavy);
    }

    int64_t pairs = (int64_t)n_heavy * (n_heavy - 1);
    value->d = (double)total / pairs;
    return CCHEM_OK;
}

/* MaxPathLength: Molecular diameter (max shortest path) */
static cchem_status_t desc_max_path_length(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->i = 0; return CCHEM_OK; }

    int max_dist = 0;

    /* BFS from each heavy atom to find max distance */
    for (int start = 0; start < mol->num_atoms; start++) {
        if (mol->atoms[start].element == ELEM_H) continue;

        int dist[TOPO_MAX_ATOMS];
        int queue[TOPO_MAX_ATOMS];
        memset(dist, -1, sizeof(dist));

        int head = 0, tail = 0;
        queue[tail++] = start;
        dist[start] = 0;

        while (head < tail) {
            int curr = queue[head++];
            if (dist[curr] > max_dist) max_dist = dist[curr];

            const atom_t* atom = &mol->atoms[curr];
            for (int i = 0; i < atom->num_neighbors; i++) {
                int nb = atom->neighbors[i];
                if (mol->atoms[nb].element == ELEM_H) continue;
                if (dist[nb] < 0) {
                    dist[nb] = dist[curr] + 1;
                    queue[tail++] = nb;
                }
            }
        }
    }

    value->i = max_dist;
    return CCHEM_OK;
}

/* ============================================================================
 * Balaban J Index
 * J = m / (mu + 1) * sum(1/sqrt(Si*Sj)) over bonds
 * where Si = sum of distances from atom i
 * ============================================================================ */

static cchem_status_t desc_balaban_j(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->d = 0.0; return CCHEM_OK; }

    int heavy_map[TOPO_MAX_ATOMS];
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms && i < TOPO_MAX_ATOMS; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_map[i] = n_heavy++;
        } else {
            heavy_map[i] = -1;
        }
    }

    if (n_heavy < 2) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    /* Compute distance sums for each atom */
    int dist_sum[TOPO_MAX_ATOMS];
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) {
            dist_sum[i] = 0;
            continue;
        }
        dist_sum[i] = bfs_distance_sum(mol, i, heavy_map, n_heavy);
    }

    /* Count heavy bonds and cyclomatic number */
    int m = 0;  /* Number of heavy bonds */
    for (int i = 0; i < mol->num_bonds; i++) {
        if (mol->atoms[mol->bonds[i].atom1].element != ELEM_H &&
            mol->atoms[mol->bonds[i].atom2].element != ELEM_H) {
            m++;
        }
    }

    int mu = m - n_heavy + 1;  /* Cyclomatic number (simplified) */
    if (mu < 0) mu = 0;

    /* Calculate J */
    double sum = 0.0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        int s1 = dist_sum[bond->atom1];
        int s2 = dist_sum[bond->atom2];
        if (s1 > 0 && s2 > 0) {
            sum += 1.0 / sqrt((double)(s1 * s2));
        }
    }

    value->d = (double)m / (mu + 1) * sum;
    return CCHEM_OK;
}

/* ============================================================================
 * Neighborhood Complexity
 * ============================================================================ */

/* MeanNeighborDegree: Average degree of neighbors */
static cchem_status_t desc_mean_neighbor_degree(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        const atom_t* atom = &mol->atoms[i];

        for (int j = 0; j < atom->num_neighbors; j++) {
            int nb = atom->neighbors[j];
            if (mol->atoms[nb].element == ELEM_H) continue;
            sum += get_heavy_degree(mol, nb);
            count++;
        }
    }

    value->d = (count > 0) ? sum / count : 0.0;
    return CCHEM_OK;
}

/* TwoHopReachability: Average number of atoms reachable in 2 hops */
static cchem_status_t desc_two_hop_reach(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->d = 0.0; return CCHEM_OK; }

    int64_t total_reach = 0;
    int n_heavy = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        n_heavy++;

        bool visited[TOPO_MAX_ATOMS] = {false};
        int reach = 0;
        visited[i] = true;

        /* 1-hop neighbors */
        const atom_t* atom = &mol->atoms[i];
        for (int j = 0; j < atom->num_neighbors; j++) {
            int nb1 = atom->neighbors[j];
            if (mol->atoms[nb1].element == ELEM_H || visited[nb1]) continue;
            visited[nb1] = true;
            reach++;

            /* 2-hop neighbors */
            const atom_t* nb1_atom = &mol->atoms[nb1];
            for (int k = 0; k < nb1_atom->num_neighbors; k++) {
                int nb2 = nb1_atom->neighbors[k];
                if (mol->atoms[nb2].element == ELEM_H || visited[nb2]) continue;
                visited[nb2] = true;
                reach++;
            }
        }
        total_reach += reach;
    }

    value->d = (n_heavy > 0) ? (double)total_reach / n_heavy : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Ionization/Polarity Topology
 * ============================================================================ */

/* IonizableAtomDensity: Fraction of ionizable atoms */
static cchem_status_t desc_ionizable_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int ionizable = 0, heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        heavy++;
        if (is_ionizable(mol, i)) ionizable++;
    }

    value->d = (heavy > 0) ? (double)ionizable / heavy : 0.0;
    return CCHEM_OK;
}

/* PolarPeripheryRatio: Polar atoms with degree 1 / total polar */
static cchem_status_t desc_polar_periphery(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int polar_periph = 0, polar_total = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (!is_polar_atom(elem)) continue;
        polar_total++;
        if (get_heavy_degree(mol, i) <= 1) polar_periph++;
    }

    value->d = (polar_total > 0) ? (double)polar_periph / polar_total : 0.0;
    return CCHEM_OK;
}

/* HydrophobicCoreIndex: Hydrophobic atoms with high degree / total */
static cchem_status_t desc_hydrophobic_core(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int core = 0, total = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!is_hydrophobic(mol, i)) continue;
        total++;
        if (get_heavy_degree(mol, i) >= 3) core++;
    }

    value->d = (total > 0) ? (double)core / total : 0.0;
    return CCHEM_OK;
}

/* TopoPolarDegreeSum: Sum of degrees of polar atoms */
static cchem_status_t desc_topo_polar_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int sum = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (is_polar_atom(elem)) {
            sum += get_heavy_degree(mol, i);
        }
    }

    value->i = sum;
    return CCHEM_OK;
}

/* ============================================================================
 * Ring System Topology
 * ============================================================================ */

/* RingSystemCount: Number of separate ring systems */
static cchem_status_t desc_ring_system_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->i = 0; return CCHEM_OK; }

    if (mol->num_rings == 0) {
        value->i = 0;
        return CCHEM_OK;
    }

    /* Mark atoms in rings */
    bool in_ring[TOPO_MAX_ATOMS] = {false};
    for (int r = 0; r < mol->num_rings; r++) {
        for (int i = 0; i < mol->rings[r].size; i++) {
            in_ring[mol->rings[r].atoms[i]] = true;
        }
    }

    /* Count connected components of ring atoms using BFS */
    bool visited[TOPO_MAX_ATOMS] = {false};
    int systems = 0;

    for (int start = 0; start < mol->num_atoms; start++) {
        if (!in_ring[start] || visited[start]) continue;

        systems++;
        int queue[TOPO_MAX_ATOMS];
        int head = 0, tail = 0;
        queue[tail++] = start;
        visited[start] = true;

        while (head < tail) {
            int curr = queue[head++];
            const atom_t* atom = &mol->atoms[curr];

            for (int i = 0; i < atom->num_neighbors; i++) {
                int nb = atom->neighbors[i];
                if (in_ring[nb] && !visited[nb]) {
                    visited[nb] = true;
                    queue[tail++] = nb;
                }
            }
        }
    }

    value->i = systems;
    return CCHEM_OK;
}

/* FusedRingAtomRatio: Atoms in multiple rings / ring atoms */
static cchem_status_t desc_fused_ring_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->d = 0.0; return CCHEM_OK; }

    int ring_count[TOPO_MAX_ATOMS] = {0};
    int ring_atoms = 0, fused_atoms = 0;

    for (int r = 0; r < mol->num_rings; r++) {
        for (int i = 0; i < mol->rings[r].size; i++) {
            ring_count[mol->rings[r].atoms[i]]++;
        }
    }

    for (int i = 0; i < mol->num_atoms; i++) {
        if (ring_count[i] > 0) {
            ring_atoms++;
            if (ring_count[i] > 1) fused_atoms++;
        }
    }

    value->d = (ring_atoms > 0) ? (double)fused_atoms / ring_atoms : 0.0;
    return CCHEM_OK;
}

/* BridgeBondCount: Bonds connecting different ring systems */
static cchem_status_t desc_bridge_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int bridges = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        /* Bridge if connects ring to non-ring, or ring to ring but not in same ring */
        if ((a1->ring_count > 0) != (a2->ring_count > 0)) {
            bridges++;  /* Ring to chain */
        } else if (a1->ring_count > 0 && a2->ring_count > 0 && !bond->in_ring) {
            bridges++;  /* Between rings but not in ring */
        }
    }

    value->i = bridges;
    return CCHEM_OK;
}

/* ============================================================================
 * Additional LogD-relevant Descriptors
 * ============================================================================ */

/* ChargeablePathSum: Sum of shortest paths between ionizable atoms */
static cchem_status_t desc_chargeable_path_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->i = 0; return CCHEM_OK; }

    /* Find ionizable atoms */
    int ionizable[64];
    int n_ionizable = 0;
    for (int i = 0; i < mol->num_atoms && n_ionizable < 64; i++) {
        if (is_ionizable(mol, i)) {
            ionizable[n_ionizable++] = i;
        }
    }

    if (n_ionizable < 2) {
        value->i = 0;
        return CCHEM_OK;
    }

    /* BFS from first ionizable to find distances to others */
    int64_t total = 0;
    for (int s = 0; s < n_ionizable; s++) {
        int start = ionizable[s];
        int dist[TOPO_MAX_ATOMS];
        int queue[TOPO_MAX_ATOMS];
        memset(dist, -1, sizeof(dist));

        int head = 0, tail = 0;
        queue[tail++] = start;
        dist[start] = 0;

        while (head < tail) {
            int curr = queue[head++];
            const atom_t* atom = &mol->atoms[curr];
            for (int i = 0; i < atom->num_neighbors; i++) {
                int nb = atom->neighbors[i];
                if (mol->atoms[nb].element == ELEM_H) continue;
                if (dist[nb] < 0) {
                    dist[nb] = dist[curr] + 1;
                    queue[tail++] = nb;
                }
            }
        }

        for (int t = s + 1; t < n_ionizable; t++) {
            if (dist[ionizable[t]] > 0) {
                total += dist[ionizable[t]];
            }
        }
    }

    value->i = total;
    return CCHEM_OK;
}

/* LipophilicChainLength: Longest chain of hydrophobic atoms */
static cchem_status_t desc_lipophilic_chain(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->i = 0; return CCHEM_OK; }

    int max_chain = 0;

    /* DFS from each terminal hydrophobic atom */
    for (int start = 0; start < mol->num_atoms; start++) {
        if (!is_hydrophobic(mol, start)) continue;
        if (get_heavy_degree(mol, start) != 1) continue;

        /* BFS to find longest path through hydrophobic atoms */
        int dist[TOPO_MAX_ATOMS];
        int queue[TOPO_MAX_ATOMS];
        memset(dist, -1, sizeof(dist));

        int head = 0, tail = 0;
        queue[tail++] = start;
        dist[start] = 0;

        while (head < tail) {
            int curr = queue[head++];
            if (dist[curr] > max_chain) max_chain = dist[curr];

            const atom_t* atom = &mol->atoms[curr];
            for (int i = 0; i < atom->num_neighbors; i++) {
                int nb = atom->neighbors[i];
                if (!is_hydrophobic(mol, nb)) continue;
                if (dist[nb] < 0) {
                    dist[nb] = dist[curr] + 1;
                    queue[tail++] = nb;
                }
            }
        }
    }

    value->i = max_chain;
    return CCHEM_OK;
}

/* PolarClusterSize: Largest connected cluster of polar atoms */
static cchem_status_t desc_polar_cluster_size(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    if (mol->num_atoms >= TOPO_MAX_ATOMS) { value->i = 0; return CCHEM_OK; }

    bool visited[TOPO_MAX_ATOMS] = {false};
    int max_cluster = 0;

    for (int start = 0; start < mol->num_atoms; start++) {
        if (!is_polar_atom(mol->atoms[start].element)) continue;
        if (visited[start]) continue;

        /* BFS to find cluster size */
        int cluster_size = 0;
        int queue[TOPO_MAX_ATOMS];
        int head = 0, tail = 0;
        queue[tail++] = start;
        visited[start] = true;

        while (head < tail) {
            int curr = queue[head++];
            cluster_size++;

            const atom_t* atom = &mol->atoms[curr];
            for (int i = 0; i < atom->num_neighbors; i++) {
                int nb = atom->neighbors[i];
                if (!is_polar_atom(mol->atoms[nb].element)) continue;
                if (!visited[nb]) {
                    visited[nb] = true;
                    queue[tail++] = nb;
                }
            }
        }

        if (cluster_size > max_cluster) max_cluster = cluster_size;
    }

    value->i = max_cluster;
    return CCHEM_OK;
}

/* SumDegreeSquared: Sum of degree^2 / n (normalized Zagreb1) */
static cchem_status_t desc_norm_zagreb1(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t sum_sq = 0;
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        n_heavy++;
        int d = get_heavy_degree(mol, i);
        sum_sq += d * d;
    }

    value->d = (n_heavy > 0) ? (double)sum_sq / n_heavy : 0.0;
    return CCHEM_OK;
}

/* RandicIndex: Sum of 1/sqrt(d1*d2) - equivalent to Chi1 but traditional name */
static cchem_status_t desc_randic(const molecule_t* mol, descriptor_value_t* value) {
    /* Same as Chi1, included for completeness with standard name */
    return desc_chi1(mol, value);
}

/* ============================================================================
 * Registration
 * ============================================================================ */

void descriptors_register_topology(void) {
    descriptor_def_t def = {0};

    /* Kier-Hall Connectivity Indices */
    #define REG_TOPO(dname, ddesc, cat, vtype, fn) do { \
        memset(&def, 0, sizeof(def)); \
        strncpy(def.name, dname, MAX_DESCRIPTOR_NAME - 1); \
        strncpy(def.description, ddesc, sizeof(def.description) - 1); \
        def.category = cat; \
        def.value_type = vtype; \
        def.compute = fn; \
        descriptor_register(&def); \
    } while(0)

    REG_TOPO("Chi0", "Zero-order connectivity index", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_chi0);
    REG_TOPO("Chi1", "First-order connectivity index (path)", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_chi1);
    REG_TOPO("Chi2Path", "Second-order path connectivity", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_chi2_path);
    REG_TOPO("Chi3Path", "Third-order path connectivity", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_chi3_path);
    REG_TOPO("Chi3Cluster", "Third-order cluster connectivity", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_chi3_cluster);
    REG_TOPO("Chi0v", "Zero-order valence connectivity", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_chi0v);
    REG_TOPO("Chi1v", "First-order valence connectivity", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_chi1v);

    /* Zagreb Indices */
    REG_TOPO("Zagreb1", "First Zagreb index (sum degree^2)", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_zagreb1);
    REG_TOPO("Zagreb2", "Second Zagreb index (sum d1*d2)", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_zagreb2);
    REG_TOPO("HyperZagreb", "Hyper-Zagreb index (sum (d1+d2)^2)", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_hyper_zagreb);

    /* Wiener-type Indices */
    REG_TOPO("WienerIndex", "Sum of all shortest path distances", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_wiener);
    REG_TOPO("MeanPathLength", "Average shortest path length", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_mean_path_length);
    REG_TOPO("MolecularDiameter", "Maximum shortest path (diameter)", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_max_path_length);

    /* Balaban Index */
    REG_TOPO("BalabanJ", "Balaban J index", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_balaban_j);

    /* Neighborhood Complexity */
    REG_TOPO("MeanNeighborDegree", "Average neighbor degree", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_mean_neighbor_degree);
    REG_TOPO("TwoHopReachability", "Average atoms reachable in 2 hops", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_two_hop_reach);

    /* Ionization Topology */
    REG_TOPO("IonizableDensity", "Fraction of ionizable atoms", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_ionizable_density);
    REG_TOPO("PolarPeripheryRatio", "Polar atoms at periphery ratio", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_polar_periphery);
    REG_TOPO("HydrophobicCoreIndex", "Hydrophobic atoms in core ratio", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_hydrophobic_core);
    REG_TOPO("TopoPolarDegreeSum", "Sum of polar atom degrees", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_topo_polar_sum);

    /* Ring System Topology */
    REG_TOPO("RingSystemCount", "Number of separate ring systems", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_ring_system_count);
    REG_TOPO("FusedRingAtomRatio", "Atoms in fused rings / ring atoms", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_fused_ring_ratio);
    REG_TOPO("BridgeBondCount", "Bonds connecting ring systems", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_bridge_bond_count);

    /* LogD-relevant */
    REG_TOPO("ChargeablePathSum", "Sum of paths between ionizable atoms", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_chargeable_path_sum);
    REG_TOPO("LipophilicChainLength", "Longest hydrophobic chain", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_lipophilic_chain);
    REG_TOPO("PolarClusterSize", "Largest polar atom cluster", DESC_CATEGORY_PROPERTIES, DESC_VALUE_INT, desc_polar_cluster_size);
    REG_TOPO("NormZagreb1", "Normalized first Zagreb (Z1/n)", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_norm_zagreb1);
    REG_TOPO("RandicIndex", "Randic connectivity index", DESC_CATEGORY_PROPERTIES, DESC_VALUE_DOUBLE, desc_randic);

    #undef REG_TOPO
}
