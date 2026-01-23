/**
 * @file coords2d.c
 * @brief 2D coordinate generation for molecular depiction
 *
 * Main coordinator for the modular 2D layout system. Delegates to specialized
 * modules for ring system detection, ring placement, chain placement, and
 * collision resolution.
 */

#include "cchem/depictor/coords2d.h"
#include "cchem/depictor/layout2d.h"
#include "cchem/depictor/force_directed.h"
#include "cchem/depictor/types.h"
#include <stdlib.h>
#include <math.h>

#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const coords2d_options_t DEFAULT_OPTIONS = COORDS2D_OPTIONS_DEFAULT;

/* ============== Main API ============== */

mol_coords_t* coords2d_generate(const molecule_t* mol, const coords2d_options_t* options) {
    if (!mol || mol->num_atoms == 0) return NULL;

    const coords2d_options_t* opts = options ? options : &DEFAULT_OPTIONS;

    /* Create coordinate storage */
    mol_coords_t* result = mol_coords_create(mol->num_atoms);
    if (!result) return NULL;
    result->has_2d = true;

    /* Create layout context */
    layout_context_t* ctx = layout_context_create(mol, result, opts);
    if (!ctx) {
        mol_coords_free(result);
        return NULL;
    }

    /* Detect and classify ring systems */
    layout_detect_ring_systems(ctx);

    if (opts->debug) {
        /* Print atom elements and bonds */
        fprintf(stderr, "DEBUG: Atom elements:\n");
        for (int i = 0; i < mol->num_atoms; i++) {
            fprintf(stderr, "  %2d: elem=%d ring=%d neighbors:", i, mol->atoms[i].element, mol->atoms[i].ring_count);
            for (int j = 0; j < mol->atoms[i].num_neighbors; j++) {
                fprintf(stderr, " %d", mol->atoms[i].neighbors[j]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "DEBUG: Double bonds:\n");
        for (int i = 0; i < mol->num_bonds; i++) {
            if (mol->bonds[i].type == BOND_DOUBLE) {
                fprintf(stderr, "  %d-%d\n", mol->bonds[i].atom1, mol->bonds[i].atom2);
            }
        }
        fprintf(stderr, "DEBUG: Detected %d ring systems, %d total rings\n",
                ctx->num_ring_systems, mol->num_rings);
        for (int s = 0; s < ctx->num_ring_systems; s++) {
            ring_system_t* sys = &ctx->ring_systems[s];
            fprintf(stderr, "  System %d: %d rings, type=%d, atoms: ",
                    s, sys->num_rings, sys->type);
            for (int i = 0; i < sys->num_atoms; i++) {
                fprintf(stderr, "%d ", sys->all_atoms[i]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "DEBUG: Ring contents:\n");
        for (int r = 0; r < mol->num_rings; r++) {
            fprintf(stderr, "  Ring %d (%d atoms): ", r, mol->rings[r].size);
            for (int i = 0; i < mol->rings[r].size; i++) {
                fprintf(stderr, "%d ", mol->rings[r].atoms[i]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "DEBUG: Connectivity for key atoms:\n");
        int key_atoms[] = {2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
        for (int k = 0; k < 26; k++) {
            int i = key_atoms[k];
            if (i < mol->num_atoms) {
                fprintf(stderr, "  Atom %d neighbors: ", i);
                for (int j = 0; j < mol->atoms[i].num_neighbors; j++) {
                    fprintf(stderr, "%d ", mol->atoms[i].neighbors[j]);
                }
                fprintf(stderr, "\n");
            }
        }
    }

    /* Place ring systems (includes chain placement between rings) */
    layout_place_ring_systems(ctx);

    if (opts->debug) {
        fprintf(stderr, "DEBUG: After ring placement, placed atoms:\n");
        int placed_count = 0;
        for (int i = 0; i < mol->num_atoms; i++) {
            if (ctx->placed[i]) {
                placed_count++;
                fprintf(stderr, "  Atom %d: (%.2f, %.2f)\n", i,
                        result->coords_2d[i].x, result->coords_2d[i].y);
            }
        }
        fprintf(stderr, "  Total placed: %d/%d\n", placed_count, mol->num_atoms);
    }

    /* Final chain placement for any remaining atoms */
    layout_place_chains(ctx);

    if (opts->debug) {
        fprintf(stderr, "DEBUG: After final chain placement:\n");
        int final_placed = 0;
        for (int i = 0; i < mol->num_atoms; i++) {
            if (ctx->placed[i]) final_placed++;
        }
        fprintf(stderr, "  Total placed: %d/%d\n", final_placed, mol->num_atoms);
        if (final_placed < mol->num_atoms) {
            fprintf(stderr, "  Unplaced atoms: ");
            for (int i = 0; i < mol->num_atoms; i++) {
                if (!ctx->placed[i]) fprintf(stderr, "%d ", i);
            }
            fprintf(stderr, "\n");
        }
    }

    /* Handle disconnected atoms */
    point2d_t* coords = result->coords_2d;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!ctx->placed[i]) {
            coords[i] = (point2d_t){i * opts->bond_length, 5 * opts->bond_length};
            ctx->placed[i] = true;
        }
    }

    /* Note: Force-directed refinement for bridged systems is disabled as it
     * tends to collapse the structure. Bridged ring systems like morphine
     * need template-based placement for good results. */
    (void)ctx; /* Suppress unused warning */

    /* Center the coordinates */
    coords2d_center(result);

    /* Clean up */
    layout_context_free(ctx);

    return result;
}

void coords2d_center(mol_coords_t* coords) {
    if (!coords || !coords->has_2d || coords->num_atoms == 0) return;

    point2d_t sum = {0, 0};
    for (int i = 0; i < coords->num_atoms; i++) {
        sum = point2d_add(sum, coords->coords_2d[i]);
    }
    point2d_t centroid = point2d_scale(sum, 1.0 / coords->num_atoms);
    for (int i = 0; i < coords->num_atoms; i++) {
        coords->coords_2d[i] = point2d_sub(coords->coords_2d[i], centroid);
    }
}

void coords2d_scale_to_fit(mol_coords_t* coords, double width, double height, double margin) {
    if (!coords || !coords->has_2d || coords->num_atoms == 0) return;

    point2d_t min = coords->coords_2d[0];
    point2d_t max = coords->coords_2d[0];

    for (int i = 1; i < coords->num_atoms; i++) {
        point2d_t p = coords->coords_2d[i];
        if (p.x < min.x) min.x = p.x;
        if (p.y < min.y) min.y = p.y;
        if (p.x > max.x) max.x = p.x;
        if (p.y > max.y) max.y = p.y;
    }

    double mol_w = fmax(max.x - min.x, 0.001);
    double mol_h = fmax(max.y - min.y, 0.001);
    double avail_w = width - 2 * margin;
    double avail_h = height - 2 * margin;
    double scale = fmin(avail_w / mol_w, avail_h / mol_h);

    double cx = width / 2.0;
    double cy = height / 2.0;
    double mol_cx = (min.x + max.x) / 2.0;
    double mol_cy = (min.y + max.y) / 2.0;

    for (int i = 0; i < coords->num_atoms; i++) {
        coords->coords_2d[i].x = cx + (coords->coords_2d[i].x - mol_cx) * scale;
        coords->coords_2d[i].y = cy + (coords->coords_2d[i].y - mol_cy) * scale;
    }
}

int coords2d_refine(mol_coords_t* coords, const molecule_t* mol,
                    const coords2d_options_t* options) {
    if (!coords || !mol || !coords->has_2d) return -1;

    const coords2d_options_t* opts = options ? options : &DEFAULT_OPTIONS;

    /* Create temporary context for refinement */
    layout_context_t* ctx = layout_context_create(mol, coords, opts);
    if (!ctx) return -1;

    /* Mark all atoms as placed */
    for (int i = 0; i < mol->num_atoms; i++) {
        ctx->placed[i] = true;
    }

    /* Run force-directed refinement with default or specified iterations */
    int iterations = opts->max_iterations > 0 ? opts->max_iterations : 100;
    int result = layout_force_directed_refine(ctx, iterations);

    layout_context_free(ctx);
    return result;
}
