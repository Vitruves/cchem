/**
 * @file coords3d.c
 * @brief 3D conformer generation using MMFF94 force field
 *
 * Uses 2D coordinates for initial x,y placement (which handles rings correctly),
 * then adds z-coordinates based on hybridization and minimizes using MMFF94.
 */

#include "cchem/depictor/coords3d.h"
#include "cchem/depictor/coords2d.h"
#include "cchem/depictor/colors.h"
#include "cchem/depictor/mmff94_types.h"
#include "cchem/depictor/mmff94_params.h"
#include "cchem/canonicalizer/bond.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846

/* ============================================================================
 * Ring and Geometry Helpers
 * ============================================================================ */

/* Check if ring is aromatic */
static bool is_ring_aromatic(const molecule_t* mol, int ring_idx) {
    const ring_t* ring = &mol->rings[ring_idx];
    for (int i = 0; i < ring->size; i++) {
        if (!mol->atoms[ring->atoms[i]].aromatic) return false;
    }
    return true;
}

/* Check if atom should be planar */
static bool is_planar(const molecule_t* mol, int idx) {
    const atom_t* atom = &mol->atoms[idx];
    if (atom->aromatic) return true;
    for (int i = 0; i < atom->num_neighbors; i++) {
        const bond_t* b = molecule_get_bond_between_const(mol, idx, atom->neighbors[i]);
        if (b && (b->type == BOND_DOUBLE || b->type == BOND_AROMATIC)) {
            return true;
        }
    }
    return false;
}

/* Get smallest ring containing atom */
static int get_smallest_ring(const molecule_t* mol, int idx) {
    int smallest = 0;
    for (int r = 0; r < mol->num_rings; r++) {
        for (int i = 0; i < mol->rings[r].size; i++) {
            if (mol->rings[r].atoms[i] == idx) {
                if (smallest == 0 || mol->rings[r].size < smallest) {
                    smallest = mol->rings[r].size;
                }
                break;
            }
        }
    }
    return smallest;
}

/* ============================================================================
 * Coordinate Scaling and Initialization
 * ============================================================================ */

/* Scale 2D coordinates to Angstrom-like units */
static void scale_to_angstroms(mol_coords_t* coords, const molecule_t* mol) {
    double sum = 0.0;
    int count = 0;

    for (int b = 0; b < mol->num_bonds; b++) {
        int i = mol->bonds[b].atom1;
        int j = mol->bonds[b].atom2;
        double dx = coords->coords_3d[j].x - coords->coords_3d[i].x;
        double dy = coords->coords_3d[j].y - coords->coords_3d[i].y;
        double dz = coords->coords_3d[j].z - coords->coords_3d[i].z;
        sum += sqrt(dx*dx + dy*dy + dz*dz);
        count++;
    }

    if (count == 0) return;

    double avg_current = sum / count;
    double target = 1.5;
    double scale = target / avg_current;

    for (int i = 0; i < mol->num_atoms; i++) {
        coords->coords_3d[i].x *= scale;
        coords->coords_3d[i].y *= scale;
        coords->coords_3d[i].z *= scale;
    }
}

/* Initialize z-coordinates for 3D effect */
static void initialize_z_coordinates(mol_coords_t* coords, const molecule_t* mol, int seed) {
    srand(seed);
    int n = mol->num_atoms;

    /* Start with all atoms at z=0 */
    for (int i = 0; i < n; i++) {
        coords->coords_3d[i].z = 0.0;
    }

    /* Process saturated rings - add chair/envelope conformations */
    for (int r = 0; r < mol->num_rings; r++) {
        if (is_ring_aromatic(mol, r)) continue;

        const ring_t* ring = &mol->rings[r];

        if (ring->size == 6) {
            /* Chair conformation for 6-membered saturated rings */
            double z_up = 0.5;
            double z_down = -0.5;
            double z_mid = 0.2;

            for (int i = 0; i < 6; i++) {
                int atom = ring->atoms[i];
                if (i == 0) coords->coords_3d[atom].z = z_up;
                else if (i == 1) coords->coords_3d[atom].z = z_mid;
                else if (i == 2) coords->coords_3d[atom].z = -z_mid;
                else if (i == 3) coords->coords_3d[atom].z = z_down;
                else if (i == 4) coords->coords_3d[atom].z = -z_mid;
                else if (i == 5) coords->coords_3d[atom].z = z_mid;
            }
        } else if (ring->size == 5) {
            /* Envelope conformation */
            for (int i = 0; i < 5; i++) {
                int atom = ring->atoms[i];
                coords->coords_3d[atom].z = (i == 0) ? 0.4 : 0.0;
            }
        } else if (ring->size == 7) {
            /* Twist-chair for 7-membered */
            for (int i = 0; i < 7; i++) {
                int atom = ring->atoms[i];
                double angle = 2.0 * PI * i / 7.0;
                coords->coords_3d[atom].z = 0.3 * sin(2.0 * angle);
            }
        }
    }

    /* Add z-displacement for sp3 atoms not in rings */
    for (int i = 0; i < n; i++) {
        if (is_planar(mol, i)) continue;
        if (get_smallest_ring(mol, i) > 0) continue;

        /* Random small z displacement */
        coords->coords_3d[i].z += 0.3 * (2.0 * rand() / RAND_MAX - 1.0);
    }
}

/* ============================================================================
 * MMFF94 Energy Minimization
 * ============================================================================ */

double coords3d_minimize(mol_coords_t* coords, const molecule_t* mol,
                         const coords3d_options_t* options) {
    if (!coords || !mol || !coords->has_3d || mol->num_atoms == 0) return 0.0;

    int n = mol->num_atoms;
    int max_iter = options ? options->max_iterations : 500;
    double tol = options ? options->gradient_tolerance : 1e-3;

    /* Allocate gradient arrays */
    double* gx = malloc(n * sizeof(double));
    double* gy = malloc(n * sizeof(double));
    double* gz = malloc(n * sizeof(double));

    if (!gx || !gy || !gz) {
        free(gx); free(gy); free(gz);
        return 0.0;
    }

    /* Create MMFF94 context */
    mmff94_context_t* mmff_ctx = mmff94_context_create(mol);
    if (!mmff_ctx) {
        free(gx); free(gy); free(gz);
        return 0.0;
    }

    /* Assign atom types and compute charges */
    cchem_status_t status = mmff94_assign_types(mol, mmff_ctx);
    if (status != CCHEM_OK) {
        mmff94_context_free(mmff_ctx);
        free(gx); free(gy); free(gz);
        return 0.0;
    }

    status = mmff94_compute_charges(mol, mmff_ctx);
    if (status != CCHEM_OK) {
        mmff94_context_free(mmff_ctx);
        free(gx); free(gy); free(gz);
        return 0.0;
    }

    /* Initial energy calculation */
    double energy = mmff94_calc_energy(mol, mmff_ctx, coords, gx, gy, gz);
    double step = 0.001;

    /* Allocate arrays to save gradients for backtracking */
    double* gx_save = malloc(n * sizeof(double));
    double* gy_save = malloc(n * sizeof(double));
    double* gz_save = malloc(n * sizeof(double));
    if (!gx_save || !gy_save || !gz_save) {
        free(gx); free(gy); free(gz);
        free(gx_save); free(gy_save); free(gz_save);
        mmff94_context_free(mmff_ctx);
        return 0.0;
    }

    /* Steepest descent minimization */
    for (int iter = 0; iter < max_iter; iter++) {
        /* Calculate gradient norm */
        double gnorm = 0.0;
        for (int i = 0; i < n; i++) {
            gnorm += gx[i]*gx[i] + gy[i]*gy[i] + gz[i]*gz[i];
        }
        gnorm = sqrt(gnorm);

        if (gnorm < tol) break;

        /* Check for NaN/inf gradients */
        if (isnan(gnorm) || isinf(gnorm)) {
            break;
        }

        /* Save current gradients for potential backtrack */
        memcpy(gx_save, gx, n * sizeof(double));
        memcpy(gy_save, gy, n * sizeof(double));
        memcpy(gz_save, gz, n * sizeof(double));

        /* Take step along negative gradient */
        for (int i = 0; i < n; i++) {
            coords->coords_3d[i].x -= step * gx[i];
            coords->coords_3d[i].y -= step * gy[i];
            coords->coords_3d[i].z -= step * gz[i];
        }

        /* Calculate new energy (and new gradients) */
        double new_energy = mmff94_calc_energy(mol, mmff_ctx, coords, gx, gy, gz);

        if (!isnan(new_energy) && new_energy < energy) {
            /* Accept step, try larger step next time */
            step *= 1.1;
            if (step > 0.01) step = 0.01;
            energy = new_energy;
        } else {
            /* Reject step, backtrack using SAVED gradients */
            for (int i = 0; i < n; i++) {
                coords->coords_3d[i].x += step * gx_save[i];
                coords->coords_3d[i].y += step * gy_save[i];
                coords->coords_3d[i].z += step * gz_save[i];
            }

            /* Restore saved gradients */
            memcpy(gx, gx_save, n * sizeof(double));
            memcpy(gy, gy_save, n * sizeof(double));
            memcpy(gz, gz_save, n * sizeof(double));

            step *= 0.5;
            if (step < 1e-8) break;
        }
    }

    free(gx_save);
    free(gy_save);
    free(gz_save);

    /* Cleanup */
    mmff94_context_free(mmff_ctx);
    free(gx);
    free(gy);
    free(gz);

    return energy;
}

/* ============================================================================
 * 3D Coordinate Generation
 * ============================================================================ */

mol_coords_t* coords3d_generate(const molecule_t* mol, const coords3d_options_t* options) {
    if (!mol || mol->num_atoms == 0) return NULL;

    coords3d_options_t opts = options ? *options : (coords3d_options_t)COORDS3D_OPTIONS_DEFAULT;

    /* First generate 2D coordinates (these handle rings correctly) */
    coords2d_options_t opts_2d = COORDS2D_OPTIONS_DEFAULT;
    mol_coords_t* coords = coords2d_generate(mol, &opts_2d);

    if (!coords) return NULL;

    coords->has_3d = true;

    /* Copy 2D to 3D as starting point */
    for (int i = 0; i < mol->num_atoms; i++) {
        coords->coords_3d[i].x = coords->coords_2d[i].x;
        coords->coords_3d[i].y = coords->coords_2d[i].y;
        coords->coords_3d[i].z = 0.0;
    }

    /* Scale to Angstrom-like units */
    scale_to_angstroms(coords, mol);

    /* Add z-coordinates for 3D structure */
    initialize_z_coordinates(coords, mol, opts.random_seed);

    /* Minimize with MMFF94 */
    coords3d_minimize(coords, mol, &opts);

    /* Post-process: restore aromatic ring geometry from 2D */
    /* The 2D coordinates have correct ring geometry; use those for x,y */
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];

        /* Check if ring is aromatic */
        bool is_aromatic = true;
        for (int i = 0; i < ring->size; i++) {
            if (!mol->atoms[ring->atoms[i]].aromatic) {
                is_aromatic = false;
                break;
            }
        }

        if (!is_aromatic) continue;

        /* Calculate current 3D centroid and average z */
        double cx3d = 0, cy3d = 0, cz = 0;
        for (int i = 0; i < ring->size; i++) {
            int idx = ring->atoms[i];
            cx3d += coords->coords_3d[idx].x;
            cy3d += coords->coords_3d[idx].y;
            cz += coords->coords_3d[idx].z;
        }
        cx3d /= ring->size;
        cy3d /= ring->size;
        cz /= ring->size;

        /* Calculate 2D centroid */
        double cx2d = 0, cy2d = 0;
        for (int i = 0; i < ring->size; i++) {
            int idx = ring->atoms[i];
            cx2d += coords->coords_2d[idx].x;
            cy2d += coords->coords_2d[idx].y;
        }
        cx2d /= ring->size;
        cy2d /= ring->size;

        /* Restore 2D ring geometry, centered at 3D centroid, flat at avg z */
        for (int i = 0; i < ring->size; i++) {
            int idx = ring->atoms[i];
            /* Use 2D relative positions + 3D center */
            coords->coords_3d[idx].x = cx3d + (coords->coords_2d[idx].x - cx2d);
            coords->coords_3d[idx].y = cy3d + (coords->coords_2d[idx].y - cy2d);
            coords->coords_3d[idx].z = cz;
        }
    }

    /* Update 2D from optimized 3D (orthographic projection) */
    for (int i = 0; i < mol->num_atoms; i++) {
        coords->coords_2d[i].x = coords->coords_3d[i].x;
        coords->coords_2d[i].y = coords->coords_3d[i].y;
    }

    return coords;
}

/* ============================================================================
 * Clash Detection
 * ============================================================================ */

int coords3d_count_clashes(const mol_coords_t* coords, const molecule_t* mol, double threshold) {
    if (!coords || !mol || !coords->has_3d) return 0;

    int clashes = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        for (int j = i + 1; j < mol->num_atoms; j++) {
            if (molecule_get_bond_between_const(mol, i, j)) continue;

            double dx = coords->coords_3d[j].x - coords->coords_3d[i].x;
            double dy = coords->coords_3d[j].y - coords->coords_3d[i].y;
            double dz = coords->coords_3d[j].z - coords->coords_3d[i].z;
            double d = sqrt(dx*dx + dy*dy + dz*dz);

            double r_min = (atom_get_vdw_radius(mol->atoms[i].element) +
                           atom_get_vdw_radius(mol->atoms[j].element)) * threshold;

            if (d < r_min) clashes++;
        }
    }

    return clashes;
}

/* ============================================================================
 * 3D Coordinate Generation with Energy Reporting
 * ============================================================================ */

mol_coords_t* coords3d_generate_with_energy(const molecule_t* mol, const coords3d_options_t* options,
                                            double* energy_initial, double* energy_final) {
    if (!mol || mol->num_atoms == 0) return NULL;

    coords3d_options_t opts = options ? *options : (coords3d_options_t)COORDS3D_OPTIONS_DEFAULT;

    /* First generate 2D coordinates (these handle rings correctly) */
    coords2d_options_t opts_2d = COORDS2D_OPTIONS_DEFAULT;
    mol_coords_t* coords = coords2d_generate(mol, &opts_2d);

    if (!coords) return NULL;

    coords->has_3d = true;

    /* Copy 2D to 3D as starting point */
    for (int i = 0; i < mol->num_atoms; i++) {
        coords->coords_3d[i].x = coords->coords_2d[i].x;
        coords->coords_3d[i].y = coords->coords_2d[i].y;
        coords->coords_3d[i].z = 0.0;
    }

    /* Scale to Angstrom-like units */
    scale_to_angstroms(coords, mol);

    /* Add z-coordinates for 3D structure */
    initialize_z_coordinates(coords, mol, opts.random_seed);

    /* Calculate initial energy before minimization */
    if (energy_initial) {
        mmff94_context_t* mmff_ctx = mmff94_context_create(mol);
        if (mmff_ctx) {
            mmff94_assign_types(mol, mmff_ctx);
            mmff94_compute_charges(mol, mmff_ctx);
            *energy_initial = mmff94_calc_energy(mol, mmff_ctx, coords, NULL, NULL, NULL);
            mmff94_context_free(mmff_ctx);
        } else {
            *energy_initial = 0.0;
        }
    }

    /* Minimize with MMFF94 */
    double final_e = coords3d_minimize(coords, mol, &opts);
    if (energy_final) {
        *energy_final = final_e;
    }

    /* Post-process: restore aromatic ring geometry from 2D */
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];

        bool is_aromatic = true;
        for (int i = 0; i < ring->size; i++) {
            if (!mol->atoms[ring->atoms[i]].aromatic) {
                is_aromatic = false;
                break;
            }
        }

        if (!is_aromatic) continue;

        double cx3d = 0, cy3d = 0, cz = 0;
        for (int i = 0; i < ring->size; i++) {
            int idx = ring->atoms[i];
            cx3d += coords->coords_3d[idx].x;
            cy3d += coords->coords_3d[idx].y;
            cz += coords->coords_3d[idx].z;
        }
        cx3d /= ring->size;
        cy3d /= ring->size;
        cz /= ring->size;

        double cx2d = 0, cy2d = 0;
        for (int i = 0; i < ring->size; i++) {
            int idx = ring->atoms[i];
            cx2d += coords->coords_2d[idx].x;
            cy2d += coords->coords_2d[idx].y;
        }
        cx2d /= ring->size;
        cy2d /= ring->size;

        for (int i = 0; i < ring->size; i++) {
            int idx = ring->atoms[i];
            coords->coords_3d[idx].x = cx3d + (coords->coords_2d[idx].x - cx2d);
            coords->coords_3d[idx].y = cy3d + (coords->coords_2d[idx].y - cy2d);
            coords->coords_3d[idx].z = cz;
        }
    }

    /* Update 2D from optimized 3D (orthographic projection) */
    for (int i = 0; i < mol->num_atoms; i++) {
        coords->coords_2d[i].x = coords->coords_3d[i].x;
        coords->coords_2d[i].y = coords->coords_3d[i].y;
    }

    return coords;
}
