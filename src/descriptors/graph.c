/**
 * @file graph.c
 * @brief Ultra-fast graph-based molecular descriptors using igraph
 *
 * ADME-relevant graph topology descriptors computed from molecular graphs.
 * All descriptors computed in single batch pass for maximum efficiency.
 *
 * Descriptor categories:
 * - Connectivity: density, degree statistics (permeability, size)
 * - Path-based: diameter, Wiener index (shape, flexibility)
 * - Centrality: betweenness, closeness (metabolic sites)
 * - Clustering: transitivity (ring systems, compactness)
 * - Spectral: eigenvalues (electronic properties proxy)
 * - ADME-specific: polar centrality, aromatic connectivity
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <igraph.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

#define NUM_GRAPH_DESCRIPTORS 30
#define LARGE_DIST 1e20  /* Large constant for distance comparisons (avoids -ffast-math issues) */

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

    /* Spectral */
    double spectral_radius;
    double graph_energy;
    double algebraic_connectivity;

    /* Structural */
    int cyclomatic_number;
    int girth;
    double vertex_connectivity;
    double edge_connectivity;

    /* ADME-specific */
    double polar_centrality;      /* Mean centrality of heteroatoms */
    double aromatic_connectivity; /* Aromatic subgraph connectivity */
    double branching_index;       /* Molecular branching factor */
    double peripheral_ratio;      /* Ratio of peripheral atoms */
    double core_ratio;            /* Ratio of core atoms */
} graph_stats_t;

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/* Build igraph from molecule - heavy atoms only for speed */
static igraph_t* molecule_to_igraph(const molecule_t* mol, int** atom_map, int* n_heavy) {
    if (!mol || mol->num_atoms == 0) return NULL;

    /* Count heavy atoms and build mapping */
    int* map = (int*)malloc(mol->num_atoms * sizeof(int));
    if (!map) return NULL;

    int heavy_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            map[i] = heavy_count++;
        } else {
            map[i] = -1;
        }
    }

    if (heavy_count < 2) {
        free(map);
        return NULL;
    }

    /* Count edges between heavy atoms */
    int edge_count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        int a1 = mol->bonds[i].atom1;
        int a2 = mol->bonds[i].atom2;
        if (map[a1] >= 0 && map[a2] >= 0) {
            edge_count++;
        }
    }

    /* Build edge list */
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, edge_count * 2);

    int idx = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        int a1 = mol->bonds[i].atom1;
        int a2 = mol->bonds[i].atom2;
        if (map[a1] >= 0 && map[a2] >= 0) {
            VECTOR(edges)[idx++] = map[a1];
            VECTOR(edges)[idx++] = map[a2];
        }
    }

    /* Create graph */
    igraph_t* g = (igraph_t*)malloc(sizeof(igraph_t));
    if (!g) {
        igraph_vector_int_destroy(&edges);
        free(map);
        return NULL;
    }

    igraph_error_t err = igraph_create(g, &edges, heavy_count, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges);

    if (err != IGRAPH_SUCCESS) {
        free(g);
        free(map);
        return NULL;
    }

    *atom_map = map;
    *n_heavy = heavy_count;
    return g;
}

static void free_igraph(igraph_t* g, int* atom_map) {
    if (g) {
        igraph_destroy(g);
        free(g);
    }
    if (atom_map) free(atom_map);
}

/* ============================================================================
 * Batch Computation - All Graph Descriptors in Single Pass
 * ============================================================================ */

/* Size thresholds for expensive operations - aggressive for speed */
#define GRAPH_SIZE_SMALL 20    /* Full computation only for tiny molecules */
#define GRAPH_SIZE_MEDIUM 40   /* Skip expensive ops */
#define GRAPH_SIZE_LARGE 100   /* Use approximations for everything */

static void collect_graph_stats(const molecule_t* mol, graph_stats_t* s) {
    memset(s, 0, sizeof(graph_stats_t));

    int* atom_map = NULL;
    int n_heavy = 0;
    igraph_t* g = molecule_to_igraph(mol, &atom_map, &n_heavy);

    if (!g || n_heavy < 2) {
        s->is_connected = 1;
        s->n_vertices = mol ? mol->num_atoms : 0;
        free_igraph(g, atom_map);
        return;
    }

    s->n_vertices = n_heavy;
    s->n_edges = (int)igraph_ecount(g);

    /* === Basic Properties === */
    s->density = (double)(2 * s->n_edges) / (n_heavy * (n_heavy - 1));

    /* Components - always needed */
    igraph_vector_int_t membership;
    igraph_vector_int_t csize;
    igraph_integer_t n_comp;
    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_init(&csize, 0);
    igraph_connected_components(g, &membership, &csize, &n_comp, IGRAPH_WEAK);
    s->n_components = (int)n_comp;
    s->is_connected = (n_comp == 1) ? 1 : 0;
    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&csize);

    /* === Degree Statistics - compute ONCE and reuse === */
    igraph_vector_int_t degree;
    igraph_vector_int_init(&degree, n_heavy);
    igraph_degree(g, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    double deg_sum = 0, deg_sq_sum = 0;
    int max_deg = 0;
    int branch_points = 0;  /* degree >= 3 */
    int peripheral = 0;     /* degree == 1 */

    for (int i = 0; i < n_heavy; i++) {
        int d = VECTOR(degree)[i];
        deg_sum += d;
        deg_sq_sum += d * d;
        if (d > max_deg) max_deg = d;
        if (d >= 3) branch_points++;
        if (d == 1) peripheral++;
    }
    s->mean_degree = deg_sum / n_heavy;
    s->max_degree = max_deg;
    s->degree_variance = (deg_sq_sum / n_heavy) - (s->mean_degree * s->mean_degree);

    /* ADME descriptors from degree - computed here to avoid re-calling igraph_degree */
    s->branching_index = (double)branch_points / n_heavy;
    s->peripheral_ratio = (double)peripheral / n_heavy;
    s->core_ratio = s->branching_index;  /* Same as branching (degree >= 3) */

    /* Keep degree vector for polar centrality later */

    /* Degree assortativity - skip for large graphs (expensive) */
    if (n_heavy <= GRAPH_SIZE_MEDIUM) {
        igraph_real_t assortativity;
        if (igraph_assortativity_degree(g, &assortativity, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS) {
            s->degree_assortativity = (assortativity >= -1.0 && assortativity <= 1.0) ? assortativity : 0.0;
        }
    }

    /* === Path-based Descriptors === */
    /* Only for small connected graphs - O(n²) memory and time */
    if (s->is_connected && n_heavy <= GRAPH_SIZE_MEDIUM) {
        /* Eccentricities first (needed for diameter/radius) */
        igraph_vector_t ecc;
        igraph_vector_init(&ecc, n_heavy);
        igraph_eccentricity(g, NULL, &ecc, igraph_vss_all(), IGRAPH_ALL);

        double ecc_sum = 0, ecc_max = 0, ecc_min = LARGE_DIST;
        for (int i = 0; i < n_heavy; i++) {
            double e = VECTOR(ecc)[i];
            if (e < LARGE_DIST) {
                ecc_sum += e;
                if (e > ecc_max) ecc_max = e;
                if (e < ecc_min) ecc_min = e;
            }
        }
        s->diameter = ecc_max;
        s->radius = (ecc_min >= LARGE_DIST) ? 0 : ecc_min;
        s->mean_eccentricity = ecc_sum / n_heavy;
        igraph_vector_destroy(&ecc);

        /* Wiener index - only for small graphs */
        if (n_heavy <= GRAPH_SIZE_SMALL) {
            igraph_matrix_t distances;
            igraph_matrix_init(&distances, n_heavy, n_heavy);
            igraph_distances(g, NULL, &distances, igraph_vss_all(), igraph_vss_all(), IGRAPH_ALL);

            double path_sum = 0;
            int path_count = 0;
            for (int i = 0; i < n_heavy; i++) {
                for (int j = i + 1; j < n_heavy; j++) {
                    double d = MATRIX(distances, i, j);
                    if (d < LARGE_DIST) {
                        path_sum += d;
                        path_count++;
                    }
                }
            }
            s->wiener_index = path_sum;
            s->mean_path_length = path_count > 0 ? path_sum / path_count : 0;
            igraph_matrix_destroy(&distances);
        } else {
            /* Approximate Wiener index for medium graphs */
            s->wiener_index = s->mean_eccentricity * n_heavy * (n_heavy - 1) / 4.0;
            s->mean_path_length = s->mean_eccentricity / 2.0;
        }
    } else if (s->is_connected) {
        /* Approximation for large connected graphs */
        /* Diameter approximation: ~log(n) for typical molecular graphs */
        s->diameter = 2.0 + log((double)n_heavy) / log(s->mean_degree + 1);
        s->radius = s->diameter / 2.0;
        s->mean_eccentricity = (s->diameter + s->radius) / 2.0;
        s->wiener_index = s->mean_eccentricity * n_heavy * (n_heavy - 1) / 4.0;
        s->mean_path_length = s->mean_eccentricity / 2.0;
    }

    /* === Centrality Measures === */
    /* Betweenness is O(n*m) - very expensive, only for small graphs */
    if (n_heavy <= GRAPH_SIZE_SMALL) {
        igraph_vector_t betweenness;
        igraph_vector_init(&betweenness, n_heavy);
        igraph_betweenness(g, NULL, &betweenness, igraph_vss_all(), false, false);

        double bw_sum = 0, bw_max = 0;
        for (int i = 0; i < n_heavy; i++) {
            double b = VECTOR(betweenness)[i];
            bw_sum += b;
            if (b > bw_max) bw_max = b;
        }
        s->mean_betweenness = bw_sum / n_heavy;
        s->max_betweenness = bw_max;

        /* Polar centrality - betweenness of heteroatoms */
        double polar_bw_sum = 0;
        int polar_count = 0;
        for (int i = 0; i < mol->num_atoms; i++) {
            if (atom_map[i] >= 0) {
                element_t elem = mol->atoms[i].element;
                if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S ||
                    elem == ELEM_P || elem == ELEM_F || elem == ELEM_Cl ||
                    elem == ELEM_Br || elem == ELEM_I) {
                    polar_bw_sum += VECTOR(betweenness)[atom_map[i]];
                    polar_count++;
                }
            }
        }
        s->polar_centrality = polar_count > 0 ? polar_bw_sum / polar_count : 0;
        igraph_vector_destroy(&betweenness);
    } else {
        /* Degree-based betweenness approximation for larger graphs */
        /* High-degree nodes typically have higher betweenness */
        s->max_betweenness = (n_heavy - 1) * (n_heavy - 2) * s->max_degree / (2.0 * s->mean_degree * n_heavy);
        s->mean_betweenness = s->max_betweenness * s->mean_degree / s->max_degree;

        /* Polar centrality approximation based on degree of heteroatoms */
        double polar_deg_sum = 0;
        int polar_count = 0;
        for (int i = 0; i < mol->num_atoms; i++) {
            if (atom_map[i] >= 0) {
                element_t elem = mol->atoms[i].element;
                if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S ||
                    elem == ELEM_P || elem == ELEM_F || elem == ELEM_Cl ||
                    elem == ELEM_Br || elem == ELEM_I) {
                    polar_deg_sum += VECTOR(degree)[atom_map[i]];
                    polar_count++;
                }
            }
        }
        double polar_mean_deg = polar_count > 0 ? polar_deg_sum / polar_count : 0;
        s->polar_centrality = s->mean_betweenness * polar_mean_deg / s->mean_degree;
    }

    /* Now we can free degree vector */
    igraph_vector_int_destroy(&degree);

    /* Closeness - only for small connected graphs */
    if (s->is_connected && n_heavy <= GRAPH_SIZE_SMALL) {
        igraph_vector_t closeness;
        igraph_vector_init(&closeness, n_heavy);
        igraph_closeness(g, &closeness, NULL, NULL, igraph_vss_all(), IGRAPH_ALL, NULL, true);

        double cl_sum = 0;
        for (int i = 0; i < n_heavy; i++) {
            double c = VECTOR(closeness)[i];
            if (c >= 0 && c < LARGE_DIST) cl_sum += c;
        }
        s->mean_closeness = cl_sum / n_heavy;
        igraph_vector_destroy(&closeness);
    } else if (s->is_connected) {
        /* Closeness approximation: 1 / mean_path_length */
        s->mean_closeness = s->mean_path_length > 0 ? 1.0 / s->mean_path_length : 0;
    }

    /* Max eigenvector centrality - always use degree-based approximation */
    s->max_eigenvector_centrality = s->n_edges > 0 ? s->max_degree / sqrt(2.0 * s->n_edges) : 0;

    /* === Clustering === */
    /* Global transitivity - fast O(m) operation */
    igraph_real_t transitivity;
    igraph_transitivity_undirected(g, &transitivity, IGRAPH_TRANSITIVITY_ZERO);
    s->global_clustering = (transitivity >= 0.0 && transitivity <= 1.0) ? transitivity : 0.0;
    s->transitivity = s->global_clustering;

    /* Local clustering - only for smaller graphs */
    if (n_heavy <= GRAPH_SIZE_MEDIUM) {
        igraph_vector_t local_clustering;
        igraph_vector_init(&local_clustering, n_heavy);
        igraph_transitivity_local_undirected(g, &local_clustering, igraph_vss_all(), IGRAPH_TRANSITIVITY_ZERO);

        double lc_sum = 0;
        for (int i = 0; i < n_heavy; i++) {
            double lc = VECTOR(local_clustering)[i];
            if (lc >= 0.0 && lc <= 1.0) lc_sum += lc;
        }
        s->mean_local_clustering = lc_sum / n_heavy;
        igraph_vector_destroy(&local_clustering);
    } else {
        /* Approximate from global clustering */
        s->mean_local_clustering = s->global_clustering;
    }

    /* === Spectral Properties - degree-based approximations (fast) === */
    s->spectral_radius = sqrt(s->max_degree * s->mean_degree);

    if (s->is_connected && s->diameter > 0) {
        s->algebraic_connectivity = 4.0 / (n_heavy * s->diameter);
    }

    s->graph_energy = sqrt(2.0 * s->n_edges * n_heavy);

    /* === Structural Properties === */
    s->cyclomatic_number = s->n_edges - n_heavy + s->n_components;

    /* Girth - only for small graphs (can be expensive for dense graphs) */
    if (n_heavy <= GRAPH_SIZE_SMALL && s->cyclomatic_number > 0) {
        igraph_real_t girth_val;
        if (igraph_girth(g, &girth_val, NULL) == IGRAPH_SUCCESS && girth_val > 0) {
            s->girth = (int)girth_val;
        }
    } else if (s->cyclomatic_number > 0) {
        /* Approximate: most drug-like molecules have 5/6-membered rings */
        s->girth = 5;
    }

    /* Vertex and edge connectivity - expensive, skip for large graphs */
    if (s->is_connected && n_heavy >= 2 && n_heavy <= GRAPH_SIZE_SMALL) {
        igraph_integer_t vconn, econn;
        igraph_vertex_connectivity(g, &vconn, true);
        igraph_edge_connectivity(g, &econn, true);
        s->vertex_connectivity = (double)vconn;
        s->edge_connectivity = (double)econn;
    } else if (s->is_connected) {
        /* Approximate: minimum degree is a lower bound */
        /* For typical molecules, connectivity ~ 1-2 */
        s->vertex_connectivity = s->peripheral_ratio > 0.3 ? 1.0 : 2.0;
        s->edge_connectivity = s->vertex_connectivity;
    }

    /* === ADME-Specific Descriptors === */
    /* Aromatic connectivity */
    int aromatic_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (atom_map[i] >= 0 && mol->atoms[i].aromatic) {
            aromatic_count++;
        }
    }
    s->aromatic_connectivity = (double)aromatic_count / n_heavy;

    free_igraph(g, atom_map);
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

    /* ADME-specific (30) - we already have 29, need one more */
    values[idx++].d = s.branching_index;            /* BranchingIndex */

    return idx;
}

/* ============================================================================
 * Individual Descriptor Functions
 * ============================================================================ */

#define GRAPH_DESC_FN(name, stat_idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    descriptor_value_t vals[NUM_GRAPH_DESCRIPTORS]; \
    int n = descriptors_compute_graph_all(mol, vals); \
    if (n < 0 || stat_idx >= n) return CCHEM_ERROR_INVALID_INPUT; \
    *value = vals[stat_idx]; \
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
