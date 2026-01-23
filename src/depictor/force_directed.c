/**
 * @file force_directed.c
 * @brief Force-directed layout refinement for 2D molecular depiction
 *
 * Implements a force-directed algorithm using spring forces for bonds,
 * repulsion forces between non-bonded atoms, and angular forces for
 * ideal bond angles.
 */

#include "cchem/depictor/force_directed.h"
#include "cchem/depictor/layout2d.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const force_params_t DEFAULT_FORCE_PARAMS = FORCE_PARAMS_DEFAULT;

/* ============== Force Calculation ============== */

point2d_t force_spring(point2d_t p1, point2d_t p2,
                       double ideal_length, double strength) {
    point2d_t delta = point2d_sub(p2, p1);
    double dist = point2d_length(delta);

    if (dist < 0.001) {
        /* Avoid division by zero - apply small random force */
        return (point2d_t){strength * 0.1, 0};
    }

    double displacement = dist - ideal_length;
    double force_mag = strength * displacement;

    /* Force direction: towards p2 if stretched, away if compressed */
    point2d_t dir = point2d_scale(delta, 1.0 / dist);
    return point2d_scale(dir, force_mag);
}

point2d_t force_repulsion(point2d_t p1, point2d_t p2,
                          double min_distance, double strength) {
    point2d_t delta = point2d_sub(p1, p2);
    double dist = point2d_length(delta);

    if (dist > min_distance * 3.0) {
        /* Too far apart - no repulsion needed */
        return (point2d_t){0, 0};
    }

    if (dist < 0.001) {
        /* Overlap - apply strong random force */
        return (point2d_t){strength * 10.0, strength * 5.0};
    }

    /* Inverse square repulsion */
    double force_mag = strength * min_distance * min_distance / (dist * dist);

    /* Cap the maximum repulsion force */
    if (force_mag > strength * 10.0) {
        force_mag = strength * 10.0;
    }

    point2d_t dir = point2d_scale(delta, 1.0 / dist);
    return point2d_scale(dir, force_mag);
}

void force_angular(point2d_t center, point2d_t neighbor1, point2d_t neighbor2,
                   double ideal_angle, double strength,
                   point2d_t* out_force1, point2d_t* out_force2) {
    /* Calculate current angle */
    point2d_t v1 = point2d_sub(neighbor1, center);
    point2d_t v2 = point2d_sub(neighbor2, center);

    double len1 = point2d_length(v1);
    double len2 = point2d_length(v2);

    if (len1 < 0.001 || len2 < 0.001) {
        *out_force1 = (point2d_t){0, 0};
        *out_force2 = (point2d_t){0, 0};
        return;
    }

    double cos_angle = point2d_dot(v1, v2) / (len1 * len2);
    if (cos_angle > 1.0) cos_angle = 1.0;
    if (cos_angle < -1.0) cos_angle = -1.0;

    double current_angle = acos(cos_angle);
    double angle_diff = current_angle - ideal_angle;

    /* Small angles - no force needed */
    if (fabs(angle_diff) < 0.01) {
        *out_force1 = (point2d_t){0, 0};
        *out_force2 = (point2d_t){0, 0};
        return;
    }

    /* Calculate perpendicular directions for rotational force */
    point2d_t perp1 = {-v1.y / len1, v1.x / len1};
    point2d_t perp2 = {-v2.y / len2, v2.x / len2};

    /* Apply rotational force based on angle deviation */
    double force_mag = strength * angle_diff;

    /* Force on neighbor1 should rotate it towards ideal angle */
    double sign = (point2d_cross(v1, v2) > 0) ? 1.0 : -1.0;

    *out_force1 = point2d_scale(perp1, -sign * force_mag);
    *out_force2 = point2d_scale(perp2, sign * force_mag);
}

/* ============== Energy Calculation ============== */

double force_calculate_energy(const layout_context_t* ctx,
                              const force_params_t* params) {
    if (!ctx || !params) return 1e10;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;

    double energy = 0.0;

    /* Spring energy from bonds */
    for (int i = 0; i < mol->num_bonds; i++) {
        int a1 = mol->bonds[i].atom1;
        int a2 = mol->bonds[i].atom2;

        if (!placed[a1] || !placed[a2]) continue;

        double dist = point2d_distance(coords[a1], coords[a2]);
        double delta = dist - params->ideal_bond_length;
        energy += 0.5 * params->spring_strength * delta * delta;
    }

    /* Repulsion energy from non-bonded pairs */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!placed[i]) continue;

        for (int j = i + 1; j < mol->num_atoms; j++) {
            if (!placed[j]) continue;

            /* Skip bonded pairs */
            bool bonded = false;
            const atom_t* atom_i = &mol->atoms[i];
            for (int k = 0; k < atom_i->num_neighbors; k++) {
                if (atom_i->neighbors[k] == j) {
                    bonded = true;
                    break;
                }
            }
            if (bonded) continue;

            double dist = point2d_distance(coords[i], coords[j]);
            if (dist < params->min_atom_distance) {
                double delta = params->min_atom_distance - dist;
                energy += params->repulsion_strength * delta * delta * 10.0;
            }
        }
    }

    return energy;
}

/* ============== Force Computation ============== */

void force_calculate_forces(const layout_context_t* ctx,
                            const force_params_t* params,
                            point2d_t* forces) {
    if (!ctx || !params || !forces) return;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;

    /* Initialize forces to zero */
    for (int i = 0; i < mol->num_atoms; i++) {
        forces[i] = (point2d_t){0, 0};
    }

    /* Spring forces from bonds */
    for (int i = 0; i < mol->num_bonds; i++) {
        int a1 = mol->bonds[i].atom1;
        int a2 = mol->bonds[i].atom2;

        if (!placed[a1] || !placed[a2]) continue;

        point2d_t f = force_spring(coords[a1], coords[a2],
                                   params->ideal_bond_length,
                                   params->spring_strength);
        forces[a1] = point2d_add(forces[a1], f);
        forces[a2] = point2d_sub(forces[a2], f);
    }

    /* Repulsion forces from non-bonded pairs */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!placed[i]) continue;

        for (int j = i + 1; j < mol->num_atoms; j++) {
            if (!placed[j]) continue;

            /* Skip bonded pairs */
            bool bonded = false;
            const atom_t* atom_i = &mol->atoms[i];
            for (int k = 0; k < atom_i->num_neighbors; k++) {
                if (atom_i->neighbors[k] == j) {
                    bonded = true;
                    break;
                }
            }
            if (bonded) continue;

            point2d_t f = force_repulsion(coords[i], coords[j],
                                          params->min_atom_distance,
                                          params->repulsion_strength);
            forces[i] = point2d_add(forces[i], f);
            forces[j] = point2d_sub(forces[j], f);
        }
    }

    /* Angular forces */
    if (params->angular_strength > 0.001) {
        for (int i = 0; i < mol->num_atoms; i++) {
            if (!placed[i]) continue;

            const atom_t* atom = &mol->atoms[i];
            if (atom->num_neighbors < 2) continue;

            double ideal = force_ideal_angle(mol, i);

            /* Apply angular forces between consecutive neighbor pairs */
            for (int n1 = 0; n1 < atom->num_neighbors; n1++) {
                for (int n2 = n1 + 1; n2 < atom->num_neighbors; n2++) {
                    int nb1 = atom->neighbors[n1];
                    int nb2 = atom->neighbors[n2];

                    if (!placed[nb1] || !placed[nb2]) continue;

                    point2d_t f1, f2;
                    force_angular(coords[i], coords[nb1], coords[nb2],
                                  ideal, params->angular_strength,
                                  &f1, &f2);

                    forces[nb1] = point2d_add(forces[nb1], f1);
                    forces[nb2] = point2d_add(forces[nb2], f2);
                }
            }
        }
    }
}

/* ============== Ideal Angle Calculation ============== */

double force_ideal_angle(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Check hybridization based on bonds */
    bool has_triple = false;
    bool has_double = false;

    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (bond_idx >= 0) {
            if (mol->bonds[bond_idx].type == BOND_TRIPLE) {
                has_triple = true;
            } else if (mol->bonds[bond_idx].type == BOND_DOUBLE ||
                       mol->bonds[bond_idx].type == BOND_AROMATIC) {
                has_double = true;
            }
        }
    }

    /* sp hybridization (linear) */
    if (has_triple) {
        return M_PI;
    }

    /* sp2 hybridization (trigonal planar) */
    if (has_double || atom->aromatic) {
        return 2.0 * M_PI / 3.0; /* 120 degrees */
    }

    /* sp3 hybridization (tetrahedral) */
    /* Note: For 2D projection, we use ~109.5 degrees */
    return 109.5 * M_PI / 180.0;
}

void force_ideal_angles_for_neighbors(int num_neighbors, double* angles) {
    if (!angles || num_neighbors < 2) return;

    /* Distribute angles evenly around the atom */
    double step = 2.0 * M_PI / num_neighbors;
    for (int i = 0; i < num_neighbors - 1; i++) {
        angles[i] = step;
    }
}

/* ============== Integration ============== */

double force_step(layout_context_t* ctx, const force_params_t* params,
                  point2d_t* velocities) {
    if (!ctx || !params || !velocities) return 1e10;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;

    /* Calculate forces */
    point2d_t* forces = malloc(mol->num_atoms * sizeof(point2d_t));
    force_calculate_forces(ctx, params, forces);

    /* Update velocities and positions (Velocity Verlet-like) */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!placed[i]) continue;

        /* Skip fixed atoms */
        if (params->fix_ring_atoms && mol->atoms[i].ring_count > 0) continue;

        /* Update velocity with damping */
        velocities[i] = point2d_add(
            point2d_scale(velocities[i], params->damping),
            forces[i]
        );

        /* Limit velocity */
        double speed = point2d_length(velocities[i]);
        if (speed > params->max_displacement) {
            velocities[i] = point2d_scale(velocities[i],
                                          params->max_displacement / speed);
        }

        /* Update position */
        coords[i] = point2d_add(coords[i], velocities[i]);
    }

    /* Calculate new energy */
    double energy = force_calculate_energy(ctx, params);

    free(forces);
    return energy;
}

int force_refine(layout_context_t* ctx, const force_params_t* params,
                 int max_iterations) {
    if (!ctx) return 0;

    const force_params_t* p = params ? params : &DEFAULT_FORCE_PARAMS;
    const molecule_t* mol = ctx->mol;

    /* Initialize velocities to zero */
    point2d_t* velocities = calloc(mol->num_atoms, sizeof(point2d_t));
    if (!velocities) return 0;

    double prev_energy = force_calculate_energy(ctx, p);
    int iterations = 0;

    for (int iter = 0; iter < max_iterations; iter++) {
        double energy = force_step(ctx, p, velocities);
        iterations++;

        /* Check for convergence */
        double energy_change = fabs(prev_energy - energy);
        if (energy_change < p->convergence_threshold) {
            break;
        }

        prev_energy = energy;
    }

    free(velocities);
    return iterations;
}

/* ============== External Interface ============== */

int layout_force_directed_refine(layout_context_t* ctx, int iterations) {
    if (!ctx) return 0;

    force_params_t params = DEFAULT_FORCE_PARAMS;
    params.ideal_bond_length = ctx->bond_length;
    params.min_atom_distance = ctx->bond_length * 0.5;

    return force_refine(ctx, &params, iterations);
}
