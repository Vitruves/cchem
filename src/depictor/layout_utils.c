/**
 * @file layout_utils.c
 * @brief Shared utilities for 2D molecular layout
 *
 * Provides layout context management and geometry helper functions
 * used by all layout modules.
 */

#include "cchem/depictor/layout2d.h"
#include "cchem/depictor/coords2d.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============== Layout Context Management ============== */

layout_context_t* layout_context_create(const molecule_t* mol,
                                        mol_coords_t* coords,
                                        const coords2d_options_t* options) {
    if (!mol || !coords) return NULL;

    layout_context_t* ctx = calloc(1, sizeof(layout_context_t));
    if (!ctx) return NULL;

    ctx->mol = mol;
    ctx->coords = coords;
    ctx->options = options;
    ctx->bond_length = options ? options->bond_length : 1.5;

    /* Allocate placement tracking */
    ctx->placed = calloc(mol->num_atoms, sizeof(bool));
    if (!ctx->placed) {
        free(ctx);
        return NULL;
    }

    /* Allocate chain direction tracking */
    ctx->chain_dir = calloc(mol->num_atoms, sizeof(int));
    if (!ctx->chain_dir) {
        free(ctx->placed);
        free(ctx);
        return NULL;
    }

    /* Allocate atom-to-ring-system mapping */
    ctx->atom_to_ring_system = malloc(mol->num_atoms * sizeof(int));
    if (!ctx->atom_to_ring_system) {
        free(ctx->chain_dir);
        free(ctx->placed);
        free(ctx);
        return NULL;
    }
    for (int i = 0; i < mol->num_atoms; i++) {
        ctx->atom_to_ring_system[i] = -1;
    }

    /* Ring systems will be populated by layout_detect_ring_systems */
    ctx->ring_systems = NULL;
    ctx->num_ring_systems = 0;

    return ctx;
}

void layout_context_free(layout_context_t* ctx) {
    if (!ctx) return;

    free(ctx->placed);
    free(ctx->chain_dir);
    free(ctx->atom_to_ring_system);

    /* Free ring systems */
    layout_free_ring_systems(ctx);

    free(ctx);
}

void layout_free_ring_systems(layout_context_t* ctx) {
    if (!ctx || !ctx->ring_systems) return;

    for (int i = 0; i < ctx->num_ring_systems; i++) {
        free(ctx->ring_systems[i].ring_indices);
        free(ctx->ring_systems[i].shared_atoms);
        free(ctx->ring_systems[i].all_atoms);
    }
    free(ctx->ring_systems);
    ctx->ring_systems = NULL;
    ctx->num_ring_systems = 0;
}

/* ============== Geometry Helpers ============== */

double layout_ring_circumradius(int num_sides, double edge_length) {
    if (num_sides < 3) return edge_length;
    return edge_length / (2.0 * sin(M_PI / num_sides));
}

point2d_t layout_centroid(const point2d_t* points, int num_points) {
    point2d_t sum = {0.0, 0.0};
    if (!points || num_points <= 0) return sum;

    for (int i = 0; i < num_points; i++) {
        sum.x += points[i].x;
        sum.y += points[i].y;
    }

    return (point2d_t){sum.x / num_points, sum.y / num_points};
}

point2d_t layout_atom_centroid(const layout_context_t* ctx,
                               const int* atom_indices, int num_atoms) {
    point2d_t sum = {0.0, 0.0};
    if (!ctx || !atom_indices || num_atoms <= 0) return sum;

    int count = 0;
    for (int i = 0; i < num_atoms; i++) {
        int idx = atom_indices[i];
        if (idx >= 0 && idx < ctx->coords->num_atoms && ctx->placed[idx]) {
            sum.x += ctx->coords->coords_2d[idx].x;
            sum.y += ctx->coords->coords_2d[idx].y;
            count++;
        }
    }

    if (count == 0) return sum;
    return (point2d_t){sum.x / count, sum.y / count};
}
