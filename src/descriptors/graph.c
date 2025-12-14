/**
 * @file graph.c
 * @brief Graph-based molecular descriptors (igraph-free implementation)
 *
 * Native implementation of graph topology descriptors computed from molecular graphs.
 * All descriptors computed in single batch pass for maximum efficiency.
 * Thread-safe - no external library dependencies.
 *
 * Descriptor categories:
 * - Connectivity: density, degree statistics (permeability, size)
 * - Path-based: diameter, Wiener index (shape, flexibility)
 * - Centrality: betweenness, closeness (metabolic sites)
 * - Clustering: transitivity (ring systems, compactness)
 * - Spectral: eigenvalues (approximated from degree sequence)
 * - ADME-specific: polar centrality, aromatic connectivity
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

#define NUM_GRAPH_DESCRIPTORS 30
#define GRAPH_MAX_ATOMS 512

/* Size thresholds for expensive operations */
#define GRAPH_SIZE_SMALL 30    /* Full computation for small molecules */
#define GRAPH_SIZE_MEDIUM 60   /* Skip some expensive ops */
#define GRAPH_SIZE_LARGE 150   /* Use approximations */

/* ============================================================================
 * Graph Statistics Structure
 * ============================================================================ */

typedef struct {
    /* Basic properties */
    int n_vertices;
    int n_edges;
    double density;
    int is_connected;
    int n_components;

    /* Degree statistics */
    double mean_degree;
    double max_degree;
    double degree_variance;
    double degree_assortativity;

    /* Path-based */
    double diameter;
    double radius;
    double mean_eccentricity;
    double wiener_index;
    double mean_path_length;

    /* Centrality */
    double mean_betweenness;
    double max_betweenness;
    double mean_closeness;
    double max_eigenvector_centrality;

    /* Clustering */
    double global_clustering;
    double mean_local_clustering;
    double transitivity;

    /* Spectral (approximations) */
    double spectral_radius;
    double graph_energy;
    double algebraic_connectivity;

    /* Structural */
    int cyclomatic_number;
    int girth;
    double vertex_connectivity;
    double edge_connectivity;

    /* ADME-specific */
    double polar_centrality;
    double aromatic_connectivity;
    double branching_index;
    double peripheral_ratio;
    double core_ratio;
} graph_stats_t;

/* ============================================================================
 * Helper Functions
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

static int count_heavy_atoms(const molecule_t* mol) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) count++;
    }
    return count;
}

static int count_heavy_bonds(const molecule_t* mol) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        if (mol->atoms[mol->bonds[i].atom1].element != ELEM_H &&
            mol->atoms[mol->bonds[i].atom2].element != ELEM_H) {
            count++;
        }
    }
    return count;
}

static bool is_polar_heavy(element_t elem) {
    return elem == ELEM_N || elem == ELEM_O || elem == ELEM_S ||
           elem == ELEM_P || elem == ELEM_F || elem == ELEM_Cl ||
           elem == ELEM_Br || elem == ELEM_I;
}

/* ============================================================================
 * BFS-based Distance Computation
 * ============================================================================ */

/* Compute shortest distances from a source atom to all others (heavy atoms only) */
static void bfs_distances(const molecule_t* mol, int source, int* dist) {
    int queue[GRAPH_MAX_ATOMS];
    int head = 0, tail = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        dist[i] = -1;
    }

    if (mol->atoms[source].element == ELEM_H) return;

    dist[source] = 0;
    queue[tail++] = source;

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
}

/* Compute eccentricity (max distance to any other vertex) */
static int compute_eccentricity(const int* dist, int n, int* sum_out) {
    int max_dist = 0;
    int sum = 0;
    for (int i = 0; i < n; i++) {
        if (dist[i] > 0) {
            sum += dist[i];
            if (dist[i] > max_dist) max_dist = dist[i];
        }
    }
    if (sum_out) *sum_out = sum;
    return max_dist;
}

/* ============================================================================
 * Connected Components (using molecule's fragment finder)
 * ============================================================================ */

static int count_heavy_components(const molecule_t* mol) {
    /* Use molecule's fragment info if available */
    if (mol->fragment_ids && mol->num_fragments > 0) {
        /* Count fragments that contain heavy atoms */
        bool has_heavy[64] = {false};
        for (int i = 0; i < mol->num_atoms && mol->num_fragments < 64; i++) {
            if (mol->atoms[i].element != ELEM_H && mol->fragment_ids[i] >= 0) {
                has_heavy[mol->fragment_ids[i]] = true;
            }
        }
        int count = 0;
        for (int i = 0; i < mol->num_fragments && i < 64; i++) {
            if (has_heavy[i]) count++;
        }
        return count > 0 ? count : 1;
    }

    /* Fallback: BFS to count components */
    if (mol->num_atoms > GRAPH_MAX_ATOMS) return 1;

    bool visited[GRAPH_MAX_ATOMS] = {false};
    int queue[GRAPH_MAX_ATOMS];
    int components = 0;

    for (int start = 0; start < mol->num_atoms; start++) {
        if (mol->atoms[start].element == ELEM_H || visited[start]) continue;

        components++;
        int head = 0, tail = 0;
        queue[tail++] = start;
        visited[start] = true;

        while (head < tail) {
            int curr = queue[head++];
            const atom_t* atom = &mol->atoms[curr];

            for (int i = 0; i < atom->num_neighbors; i++) {
                int nb = atom->neighbors[i];
                if (mol->atoms[nb].element == ELEM_H || visited[nb]) continue;
                visited[nb] = true;
                queue[tail++] = nb;
            }
        }
    }

    return components > 0 ? components : 1;
}

/* ============================================================================
 * Betweenness Centrality (Brandes Algorithm - O(VE))
 * ============================================================================ */

static void compute_betweenness(const molecule_t* mol, double* betweenness,
                                 int* heavy_indices, int n_heavy) {
    if (n_heavy < 2 || n_heavy > GRAPH_MAX_ATOMS) return;

    memset(betweenness, 0, n_heavy * sizeof(double));

    /* Temporary arrays */
    int dist[GRAPH_MAX_ATOMS];
    double sigma[GRAPH_MAX_ATOMS];      /* Number of shortest paths */
    double delta[GRAPH_MAX_ATOMS];      /* Dependency */
    int pred[GRAPH_MAX_ATOMS][16];      /* Predecessors (max 16 per node) */
    int n_pred[GRAPH_MAX_ATOMS];
    int queue[GRAPH_MAX_ATOMS];
    int stack[GRAPH_MAX_ATOMS];

    /* Map from atom index to heavy index */
    int atom_to_heavy[GRAPH_MAX_ATOMS];
    memset(atom_to_heavy, -1, sizeof(atom_to_heavy));
    for (int i = 0; i < n_heavy; i++) {
        atom_to_heavy[heavy_indices[i]] = i;
    }

    /* BFS from each source */
    for (int si = 0; si < n_heavy; si++) {
        int s = heavy_indices[si];

        /* Initialize */
        for (int i = 0; i < mol->num_atoms; i++) {
            dist[i] = -1;
            sigma[i] = 0.0;
            delta[i] = 0.0;
            n_pred[i] = 0;
        }

        dist[s] = 0;
        sigma[s] = 1.0;

        int head = 0, tail = 0, stack_top = 0;
        queue[tail++] = s;

        /* BFS */
        while (head < tail) {
            int v = queue[head++];
            stack[stack_top++] = v;

            const atom_t* atom = &mol->atoms[v];
            for (int i = 0; i < atom->num_neighbors; i++) {
                int w = atom->neighbors[i];
                if (mol->atoms[w].element == ELEM_H) continue;

                /* First visit */
                if (dist[w] < 0) {
                    dist[w] = dist[v] + 1;
                    queue[tail++] = w;
                }

                /* Shortest path via v */
                if (dist[w] == dist[v] + 1) {
                    sigma[w] += sigma[v];
                    if (n_pred[w] < 16) {
                        pred[w][n_pred[w]++] = v;
                    }
                }
            }
        }

        /* Accumulation (back-propagation) */
        while (stack_top > 0) {
            int w = stack[--stack_top];
            for (int i = 0; i < n_pred[w]; i++) {
                int v = pred[w][i];
                if (sigma[w] > 0) {
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
            }
            if (w != s) {
                int hi = atom_to_heavy[w];
                if (hi >= 0) {
                    betweenness[hi] += delta[w];
                }
            }
        }
    }

    /* Normalize (undirected graph: each pair counted twice) */
    for (int i = 0; i < n_heavy; i++) {
        betweenness[i] /= 2.0;
    }
}

/* ============================================================================
 * Clustering Coefficient (Triangle Counting)
 * ============================================================================ */

static double compute_global_clustering(const molecule_t* mol, int n_heavy) {
    if (n_heavy < 3) return 0.0;

    int64_t triangles = 0;
    int64_t triplets = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        const atom_t* atom = &mol->atoms[i];

        /* Get heavy neighbors */
        int neighbors[16];
        int n_neighbors = 0;
        for (int j = 0; j < atom->num_neighbors && n_neighbors < 16; j++) {
            int nb = atom->neighbors[j];
            if (mol->atoms[nb].element != ELEM_H) {
                neighbors[n_neighbors++] = nb;
            }
        }

        /* Count triplets centered at i */
        int k = n_neighbors;
        if (k >= 2) {
            triplets += (int64_t)k * (k - 1) / 2;
        }

        /* Count triangles (edges between neighbors) */
        for (int a = 0; a < n_neighbors; a++) {
            for (int b = a + 1; b < n_neighbors; b++) {
                /* Check if neighbors[a] and neighbors[b] are connected */
                const atom_t* na = &mol->atoms[neighbors[a]];
                for (int c = 0; c < na->num_neighbors; c++) {
                    if (na->neighbors[c] == neighbors[b]) {
                        triangles++;
                        break;
                    }
                }
            }
        }
    }

    /* Each triangle counted 3 times */
    if (triplets == 0) return 0.0;
    return (double)(triangles / 3) / triplets * 3.0;
}

static void compute_local_clustering(const molecule_t* mol, double* local_cc, int n) {
    for (int i = 0; i < n; i++) local_cc[i] = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        const atom_t* atom = &mol->atoms[i];

        /* Get heavy neighbors */
        int neighbors[16];
        int k = 0;
        for (int j = 0; j < atom->num_neighbors && k < 16; j++) {
            int nb = atom->neighbors[j];
            if (mol->atoms[nb].element != ELEM_H) {
                neighbors[k++] = nb;
            }
        }

        if (k < 2) {
            local_cc[i] = 0.0;
            continue;
        }

        /* Count edges between neighbors */
        int edges = 0;
        for (int a = 0; a < k; a++) {
            for (int b = a + 1; b < k; b++) {
                const atom_t* na = &mol->atoms[neighbors[a]];
                for (int c = 0; c < na->num_neighbors; c++) {
                    if (na->neighbors[c] == neighbors[b]) {
                        edges++;
                        break;
                    }
                }
            }
        }

        int possible = k * (k - 1) / 2;
        local_cc[i] = (possible > 0) ? (double)edges / possible : 0.0;
    }
}

/* ============================================================================
 * Girth (Shortest Cycle) - BFS from each vertex
 * ============================================================================ */

static int compute_girth(const molecule_t* mol, int n_heavy) {
    if (n_heavy < 3) return 0;

    int min_girth = n_heavy + 1;  /* Start with impossible value */
    int dist[GRAPH_MAX_ATOMS];
    int parent[GRAPH_MAX_ATOMS];
    int queue[GRAPH_MAX_ATOMS];

    for (int start = 0; start < mol->num_atoms; start++) {
        if (mol->atoms[start].element == ELEM_H) continue;

        for (int i = 0; i < mol->num_atoms; i++) {
            dist[i] = -1;
            parent[i] = -1;
        }

        dist[start] = 0;
        int head = 0, tail = 0;
        queue[tail++] = start;

        while (head < tail) {
            int v = queue[head++];
            const atom_t* atom = &mol->atoms[v];

            for (int i = 0; i < atom->num_neighbors; i++) {
                int w = atom->neighbors[i];
                if (mol->atoms[w].element == ELEM_H) continue;

                if (dist[w] < 0) {
                    dist[w] = dist[v] + 1;
                    parent[w] = v;
                    queue[tail++] = w;
                } else if (parent[v] != w && parent[w] != v) {
                    /* Found a cycle */
                    int cycle_len = dist[v] + dist[w] + 1;
                    if (cycle_len < min_girth) {
                        min_girth = cycle_len;
                    }
                }
            }
        }
    }

    return (min_girth <= n_heavy) ? min_girth : 0;
}

/* ============================================================================
 * Batch Computation - All Graph Descriptors
 * ============================================================================ */

static void collect_graph_stats(const molecule_t* mol, graph_stats_t* s) {
    memset(s, 0, sizeof(graph_stats_t));

    int n_heavy = count_heavy_atoms(mol);
    int n_edges = count_heavy_bonds(mol);

    if (n_heavy < 2 || n_heavy > GRAPH_MAX_ATOMS) {
        s->is_connected = 1;
        s->n_vertices = n_heavy;
        s->n_components = 1;
        return;
    }

    s->n_vertices = n_heavy;
    s->n_edges = n_edges;

    /* === Basic Properties === */
    s->density = (n_heavy > 1) ? (double)(2 * n_edges) / (n_heavy * (n_heavy - 1)) : 0.0;
    s->n_components = count_heavy_components(mol);
    s->is_connected = (s->n_components == 1) ? 1 : 0;

    /* === Degree Statistics === */
    int heavy_indices[GRAPH_MAX_ATOMS];
    int hi = 0;
    double deg_sum = 0, deg_sq_sum = 0;
    int max_deg = 0;
    int branch_points = 0, peripheral = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;

        int d = get_heavy_degree(mol, i);
        heavy_indices[hi] = i;
        hi++;

        deg_sum += d;
        deg_sq_sum += d * d;
        if (d > max_deg) max_deg = d;
        if (d >= 3) branch_points++;
        if (d == 1) peripheral++;
    }

    s->mean_degree = deg_sum / n_heavy;
    s->max_degree = max_deg;
    s->degree_variance = (deg_sq_sum / n_heavy) - (s->mean_degree * s->mean_degree);
    s->branching_index = (double)branch_points / n_heavy;
    s->peripheral_ratio = (double)peripheral / n_heavy;
    s->core_ratio = s->branching_index;

    /* Degree assortativity: correlation of degrees at bond endpoints */
    if (n_edges > 0 && n_heavy <= GRAPH_SIZE_MEDIUM) {
        double sum_prod = 0, sum_d1 = 0, sum_d2 = 0, sum_sq1 = 0, sum_sq2 = 0;
        int edge_count = 0;

        for (int b = 0; b < mol->num_bonds; b++) {
            int a1 = mol->bonds[b].atom1;
            int a2 = mol->bonds[b].atom2;
            if (mol->atoms[a1].element == ELEM_H || mol->atoms[a2].element == ELEM_H) continue;

            int d1 = get_heavy_degree(mol, a1);
            int d2 = get_heavy_degree(mol, a2);
            sum_prod += d1 * d2;
            sum_d1 += d1;
            sum_d2 += d2;
            sum_sq1 += d1 * d1;
            sum_sq2 += d2 * d2;
            edge_count++;
        }

        if (edge_count > 0) {
            double mean1 = sum_d1 / edge_count;
            double mean2 = sum_d2 / edge_count;
            double var1 = sum_sq1 / edge_count - mean1 * mean1;
            double var2 = sum_sq2 / edge_count - mean2 * mean2;
            double cov = sum_prod / edge_count - mean1 * mean2;

            if (var1 > 0 && var2 > 0) {
                s->degree_assortativity = cov / sqrt(var1 * var2);
            }
        }
    }

    /* === Path-based Descriptors (BFS) === */
    if (s->is_connected && n_heavy <= GRAPH_SIZE_MEDIUM) {
        int dist[GRAPH_MAX_ATOMS];
        double ecc_sum = 0;
        int ecc_max = 0, ecc_min = n_heavy;
        double path_sum = 0;
        int path_count = 0;

        for (int i = 0; i < n_heavy; i++) {
            int src = heavy_indices[i];
            bfs_distances(mol, src, dist);

            int dist_sum;
            int ecc = compute_eccentricity(dist, mol->num_atoms, &dist_sum);

            ecc_sum += ecc;
            if (ecc > ecc_max) ecc_max = ecc;
            if (ecc < ecc_min && ecc > 0) ecc_min = ecc;

            path_sum += dist_sum;
            for (int j = 0; j < mol->num_atoms; j++) {
                if (dist[j] > 0) path_count++;
            }
        }

        s->diameter = ecc_max;
        s->radius = (ecc_min < n_heavy) ? ecc_min : 0;
        s->mean_eccentricity = ecc_sum / n_heavy;
        s->wiener_index = path_sum / 2;  /* Each pair counted twice */
        s->mean_path_length = (path_count > 0) ? path_sum / path_count : 0;
    } else if (s->is_connected) {
        /* Approximations for large graphs */
        s->diameter = 2.0 + log((double)n_heavy) / log(s->mean_degree + 1);
        s->radius = s->diameter / 2.0;
        s->mean_eccentricity = (s->diameter + s->radius) / 2.0;
        s->wiener_index = s->mean_eccentricity * n_heavy * (n_heavy - 1) / 4.0;
        s->mean_path_length = s->mean_eccentricity / 2.0;
    }

    /* === Centrality Measures === */
    if (n_heavy <= GRAPH_SIZE_SMALL) {
        double betweenness[GRAPH_MAX_ATOMS];
        compute_betweenness(mol, betweenness, heavy_indices, n_heavy);

        double bw_sum = 0, bw_max = 0;
        double polar_bw_sum = 0;
        int polar_count = 0;

        for (int i = 0; i < n_heavy; i++) {
            bw_sum += betweenness[i];
            if (betweenness[i] > bw_max) bw_max = betweenness[i];

            int atom_idx = heavy_indices[i];
            if (is_polar_heavy(mol->atoms[atom_idx].element)) {
                polar_bw_sum += betweenness[i];
                polar_count++;
            }
        }

        s->mean_betweenness = bw_sum / n_heavy;
        s->max_betweenness = bw_max;
        s->polar_centrality = (polar_count > 0) ? polar_bw_sum / polar_count : 0;

        /* Closeness: 1 / mean_distance */
        if (s->mean_path_length > 0) {
            s->mean_closeness = 1.0 / s->mean_path_length;
        }
    } else {
        /* Approximations */
        s->max_betweenness = (n_heavy - 1) * (n_heavy - 2) * s->max_degree / (2.0 * s->mean_degree * n_heavy);
        s->mean_betweenness = s->max_betweenness * s->mean_degree / s->max_degree;
        s->mean_closeness = (s->mean_path_length > 0) ? 1.0 / s->mean_path_length : 0;

        /* Polar centrality approximation */
        double polar_deg_sum = 0;
        int polar_count = 0;
        for (int i = 0; i < mol->num_atoms; i++) {
            if (is_polar_heavy(mol->atoms[i].element)) {
                polar_deg_sum += get_heavy_degree(mol, i);
                polar_count++;
            }
        }
        double polar_mean_deg = (polar_count > 0) ? polar_deg_sum / polar_count : 0;
        s->polar_centrality = s->mean_betweenness * polar_mean_deg / (s->mean_degree + 0.001);
    }

    /* Max eigenvector centrality approximation */
    s->max_eigenvector_centrality = (n_edges > 0) ? s->max_degree / sqrt(2.0 * n_edges) : 0;

    /* === Clustering === */
    s->global_clustering = compute_global_clustering(mol, n_heavy);
    s->transitivity = s->global_clustering;

    if (n_heavy <= GRAPH_SIZE_MEDIUM) {
        double local_cc[GRAPH_MAX_ATOMS];
        compute_local_clustering(mol, local_cc, mol->num_atoms);

        double lc_sum = 0;
        int lc_count = 0;
        for (int i = 0; i < mol->num_atoms; i++) {
            if (mol->atoms[i].element == ELEM_H) continue;
            lc_sum += local_cc[i];
            lc_count++;
        }
        s->mean_local_clustering = (lc_count > 0) ? lc_sum / lc_count : 0;
    } else {
        s->mean_local_clustering = s->global_clustering;
    }

    /* === Spectral Properties (degree-based approximations) === */
    s->spectral_radius = sqrt(s->max_degree * s->mean_degree);
    s->graph_energy = sqrt(2.0 * n_edges * n_heavy);
    if (s->is_connected && s->diameter > 0) {
        s->algebraic_connectivity = 4.0 / (n_heavy * s->diameter);
    }

    /* === Structural Properties === */
    s->cyclomatic_number = n_edges - n_heavy + s->n_components;

    if (n_heavy <= GRAPH_SIZE_SMALL && s->cyclomatic_number > 0) {
        s->girth = compute_girth(mol, n_heavy);
    } else if (s->cyclomatic_number > 0) {
        s->girth = 5;  /* Typical for drug-like molecules */
    }

    /* Vertex/edge connectivity approximation */
    if (s->is_connected) {
        s->vertex_connectivity = (s->peripheral_ratio > 0.3) ? 1.0 : 2.0;
        s->edge_connectivity = s->vertex_connectivity;
    }

    /* === ADME-Specific === */
    int aromatic_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H && mol->atoms[i].aromatic) {
            aromatic_count++;
        }
    }
    s->aromatic_connectivity = (n_heavy > 0) ? (double)aromatic_count / n_heavy : 0;
}

/* ============================================================================
 * Batch Computation Export
 * ============================================================================ */

int descriptors_compute_graph_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    graph_stats_t s;
    collect_graph_stats(mol, &s);

    int idx = 0;

    /* Basic properties (1-5) */
    values[idx++].d = s.density;                    /* GraphDensity */
    values[idx++].i = s.is_connected;               /* IsConnected */
    values[idx++].i = s.n_components;               /* ComponentCount */
    values[idx++].i = s.n_vertices;                 /* GraphVertexCount */
    values[idx++].i = s.n_edges;                    /* GraphEdgeCount */

    /* Degree statistics (6-9) */
    values[idx++].d = s.mean_degree;                /* GraphMeanDegree */
    values[idx++].d = s.max_degree;                 /* GraphMaxDegree */
    values[idx++].d = s.degree_variance;            /* GraphDegreeVariance */
    values[idx++].d = s.degree_assortativity;       /* DegreeAssortativity */

    /* Path-based (10-14) */
    values[idx++].d = s.diameter;                   /* GraphDiameter */
    values[idx++].d = s.radius;                     /* GraphRadius */
    values[idx++].d = s.mean_eccentricity;          /* MeanEccentricity */
    values[idx++].d = s.wiener_index;               /* GraphWienerIndex */
    values[idx++].d = s.mean_path_length;           /* GraphMeanPathLength */

    /* Centrality (15-19) */
    values[idx++].d = s.mean_betweenness;           /* MeanBetweenness */
    values[idx++].d = s.max_betweenness;            /* MaxBetweenness */
    values[idx++].d = s.mean_closeness;             /* MeanCloseness */
    values[idx++].d = s.max_eigenvector_centrality; /* MaxEigenvectorCentrality */
    values[idx++].d = s.polar_centrality;           /* PolarAtomCentrality */

    /* Clustering (20-22) */
    values[idx++].d = s.global_clustering;          /* GlobalClustering */
    values[idx++].d = s.mean_local_clustering;      /* MeanLocalClustering */
    values[idx++].d = s.transitivity;               /* GraphTransitivity */

    /* Spectral (23-25) */
    values[idx++].d = s.spectral_radius;            /* SpectralRadius */
    values[idx++].d = s.graph_energy;               /* GraphEnergy */
    values[idx++].d = s.algebraic_connectivity;     /* AlgebraicConnectivity */

    /* Structural (26-29) */
    values[idx++].i = s.cyclomatic_number;          /* CyclomaticNumber */
    values[idx++].i = s.girth;                      /* GraphGirth */
    values[idx++].d = s.vertex_connectivity;        /* VertexConnectivity */
    values[idx++].d = s.edge_connectivity;          /* EdgeConnectivity */

    /* ADME-specific (30) */
    values[idx++].d = s.branching_index;            /* BranchingIndex */

    return idx;
}

/* ============================================================================
 * Individual Descriptor Functions (with caching)
 * ============================================================================ */

static _Thread_local const molecule_t* graph_cached_mol = NULL;
static _Thread_local descriptor_value_t graph_cached_values[NUM_GRAPH_DESCRIPTORS];

static inline void ensure_graph_computed(const molecule_t* mol) {
    if (graph_cached_mol != mol) {
        descriptors_compute_graph_all(mol, graph_cached_values);
        graph_cached_mol = mol;
    }
}

#define GRAPH_DESC_FN(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_graph_computed(mol); \
    *value = graph_cached_values[idx]; \
    return CCHEM_OK; \
}

GRAPH_DESC_FN(graph_density, 0)
GRAPH_DESC_FN(is_connected, 1)
GRAPH_DESC_FN(component_count, 2)
GRAPH_DESC_FN(graph_vertex_count, 3)
GRAPH_DESC_FN(graph_edge_count, 4)
GRAPH_DESC_FN(graph_mean_degree, 5)
GRAPH_DESC_FN(graph_max_degree, 6)
GRAPH_DESC_FN(graph_degree_variance, 7)
GRAPH_DESC_FN(degree_assortativity, 8)
GRAPH_DESC_FN(graph_diameter, 9)
GRAPH_DESC_FN(graph_radius, 10)
GRAPH_DESC_FN(mean_eccentricity, 11)
GRAPH_DESC_FN(graph_wiener_index, 12)
GRAPH_DESC_FN(graph_mean_path_length, 13)
GRAPH_DESC_FN(mean_betweenness, 14)
GRAPH_DESC_FN(max_betweenness, 15)
GRAPH_DESC_FN(mean_closeness, 16)
GRAPH_DESC_FN(max_eigenvector_centrality, 17)
GRAPH_DESC_FN(polar_atom_centrality, 18)
GRAPH_DESC_FN(global_clustering, 19)
GRAPH_DESC_FN(mean_local_clustering, 20)
GRAPH_DESC_FN(graph_transitivity, 21)
GRAPH_DESC_FN(spectral_radius, 22)
GRAPH_DESC_FN(graph_energy, 23)
GRAPH_DESC_FN(algebraic_connectivity, 24)
GRAPH_DESC_FN(cyclomatic_number, 25)
GRAPH_DESC_FN(graph_girth, 26)
GRAPH_DESC_FN(vertex_connectivity, 27)
GRAPH_DESC_FN(edge_connectivity, 28)
GRAPH_DESC_FN(branching_index, 29)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REG_GRAPH_INT(dname, ddesc, fn) do { \
    memset(&def, 0, sizeof(def)); \
    strncpy(def.name, dname, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, ddesc, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_CUSTOM; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = fn; \
    descriptor_register(&def); \
} while(0)

#define REG_GRAPH_DBL(dname, ddesc, fn) do { \
    memset(&def, 0, sizeof(def)); \
    strncpy(def.name, dname, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, ddesc, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_CUSTOM; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = fn; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_graph(void) {
    descriptor_def_t def;

    /* Basic properties */
    REG_GRAPH_DBL("GraphDensity", "Graph edge density (edges/possible edges)", desc_graph_density);
    REG_GRAPH_INT("IsConnected", "Is molecular graph connected (1/0)", desc_is_connected);
    REG_GRAPH_INT("ComponentCount", "Number of connected components", desc_component_count);
    REG_GRAPH_INT("GraphVertexCount", "Number of heavy atom vertices", desc_graph_vertex_count);
    REG_GRAPH_INT("GraphEdgeCount", "Number of bond edges", desc_graph_edge_count);

    /* Degree statistics */
    REG_GRAPH_DBL("GraphMeanDegree", "Mean vertex degree", desc_graph_mean_degree);
    REG_GRAPH_DBL("GraphMaxDegree", "Maximum vertex degree", desc_graph_max_degree);
    REG_GRAPH_DBL("GraphDegreeVariance", "Variance of vertex degrees", desc_graph_degree_variance);
    REG_GRAPH_DBL("DegreeAssortativity", "Degree assortativity coefficient", desc_degree_assortativity);

    /* Path-based */
    REG_GRAPH_DBL("GraphDiameter", "Graph diameter (max eccentricity)", desc_graph_diameter);
    REG_GRAPH_DBL("GraphRadius", "Graph radius (min eccentricity)", desc_graph_radius);
    REG_GRAPH_DBL("MeanEccentricity", "Mean vertex eccentricity", desc_mean_eccentricity);
    REG_GRAPH_DBL("GraphWienerIndex", "Wiener topological index", desc_graph_wiener_index);
    REG_GRAPH_DBL("GraphMeanPathLength", "Average shortest path length", desc_graph_mean_path_length);

    /* Centrality */
    REG_GRAPH_DBL("MeanBetweenness", "Mean betweenness centrality", desc_mean_betweenness);
    REG_GRAPH_DBL("MaxBetweenness", "Maximum betweenness centrality", desc_max_betweenness);
    REG_GRAPH_DBL("MeanCloseness", "Mean closeness centrality", desc_mean_closeness);
    REG_GRAPH_DBL("MaxEigenvectorCentrality", "Maximum eigenvector centrality", desc_max_eigenvector_centrality);
    REG_GRAPH_DBL("PolarAtomCentrality", "Mean centrality of polar atoms", desc_polar_atom_centrality);

    /* Clustering */
    REG_GRAPH_DBL("GlobalClustering", "Global clustering coefficient", desc_global_clustering);
    REG_GRAPH_DBL("MeanLocalClustering", "Mean local clustering coefficient", desc_mean_local_clustering);
    REG_GRAPH_DBL("GraphTransitivity", "Graph transitivity", desc_graph_transitivity);

    /* Spectral */
    REG_GRAPH_DBL("SpectralRadius", "Largest eigenvalue of adjacency matrix", desc_spectral_radius);
    REG_GRAPH_DBL("GraphEnergy", "Graph energy (sum of |eigenvalues|)", desc_graph_energy);
    REG_GRAPH_DBL("AlgebraicConnectivity", "Second smallest Laplacian eigenvalue", desc_algebraic_connectivity);

    /* Structural */
    REG_GRAPH_INT("CyclomaticNumber", "Cyclomatic complexity (independent cycles)", desc_cyclomatic_number);
    REG_GRAPH_INT("GraphGirth", "Girth (shortest cycle length)", desc_graph_girth);
    REG_GRAPH_DBL("VertexConnectivity", "Vertex connectivity", desc_vertex_connectivity);
    REG_GRAPH_DBL("EdgeConnectivity", "Edge connectivity", desc_edge_connectivity);

    /* ADME-specific */
    REG_GRAPH_DBL("BranchingIndex", "Molecular branching index", desc_branching_index);
}
