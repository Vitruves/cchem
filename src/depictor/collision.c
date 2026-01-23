/**
 * @file collision.c
 * @brief Collision detection and resolution for 2D molecular layout
 *
 * Provides atom overlap detection, bond crossing detection, and
 * collision-free placement algorithms.
 */

#include "cchem/depictor/layout2d.h"
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============== Ring Collision Detection ============== */

bool layout_check_ring_collision(const layout_context_t* ctx,
                                 point2d_t center, double radius,
                                 const ring_t* ring, int exclude_atom) {
    if (!ctx || !ring) return false;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;

    for (int k = 0; k < mol->num_atoms; k++) {
        if (!placed[k]) continue;
        if (k == exclude_atom) continue;

        /* Check if atom is part of this ring */
        bool in_this_ring = false;
        for (int m = 0; m < ring->size; m++) {
            if (ring->atoms[m] == k) {
                in_this_ring = true;
                break;
            }
        }
        if (in_this_ring) continue;

        double dist = point2d_distance(center, coords[k]);
        /* If any atom is within ring radius of proposed center, collision */
        if (dist < radius * 1.2) {
            return true;
        }
    }

    return false;
}

bool layout_find_collision_free_angle(const layout_context_t* ctx,
                                      point2d_t neighbor_pos,
                                      const ring_t* ring,
                                      double base_angle, double radius,
                                      int exclude_atom,
                                      double* out_angle) {
    if (!ctx || !ring || !out_angle) return false;

    double bond_length = ctx->bond_length;

    /* Try different angles to find collision-free placement */
    double test_angles[] = {0, M_PI, 2*M_PI/3, -2*M_PI/3, M_PI/3, -M_PI/3, M_PI/2, -M_PI/2};
    int num_test = 8;

    for (int t = 0; t < num_test; t++) {
        double test_angle = base_angle + test_angles[t];
        point2d_t test_dir = {cos(test_angle), sin(test_angle)};
        point2d_t test_anchor = point2d_add(neighbor_pos,
                                            point2d_scale(test_dir, bond_length));
        point2d_t test_center = point2d_add(test_anchor,
                                            point2d_scale(test_dir, radius));

        if (!layout_check_ring_collision(ctx, test_center, radius, ring, exclude_atom)) {
            *out_angle = test_angle;
            return true;
        }
    }

    return false;
}

/* ============== Atom Overlap Detection ============== */

int layout_count_atom_overlaps(const layout_context_t* ctx, double threshold) {
    if (!ctx) return 0;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    int count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (!placed[i]) continue;

        for (int j = i + 1; j < mol->num_atoms; j++) {
            if (!placed[j]) continue;

            /* Skip bonded pairs - they're supposed to be close */
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
            if (dist < threshold) {
                count++;
            }
        }
    }

    return count;
}

/* ============== Bond Crossing Detection ============== */

/* Check if line segments (p1,p2) and (p3,p4) intersect */
static bool segments_intersect(point2d_t p1, point2d_t p2,
                               point2d_t p3, point2d_t p4) {
    /* Using cross product method */
    point2d_t d1 = point2d_sub(p2, p1);
    point2d_t d2 = point2d_sub(p4, p3);
    point2d_t d3 = point2d_sub(p3, p1);

    double cross1 = point2d_cross(d1, d2);

    /* Parallel or collinear */
    if (fabs(cross1) < 1e-10) return false;

    double t = point2d_cross(d3, d2) / cross1;
    double u = point2d_cross(d3, d1) / cross1;

    /* Check if intersection is within both segments (excluding endpoints) */
    double eps = 0.01;
    return (t > eps && t < 1.0 - eps && u > eps && u < 1.0 - eps);
}

int layout_count_bond_crossings(const layout_context_t* ctx) {
    if (!ctx) return 0;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    int count = 0;

    for (int i = 0; i < mol->num_bonds; i++) {
        int a1 = mol->bonds[i].atom1;
        int a2 = mol->bonds[i].atom2;

        if (!placed[a1] || !placed[a2]) continue;

        for (int j = i + 1; j < mol->num_bonds; j++) {
            int b1 = mol->bonds[j].atom1;
            int b2 = mol->bonds[j].atom2;

            if (!placed[b1] || !placed[b2]) continue;

            /* Skip if bonds share an atom */
            if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;

            if (segments_intersect(coords[a1], coords[a2],
                                   coords[b1], coords[b2])) {
                count++;
            }
        }
    }

    return count;
}

/* ============== Layout Quality Score ============== */

double layout_quality_score(const layout_context_t* ctx) {
    if (!ctx) return 1e10;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    double bond_length = ctx->bond_length;

    double score = 0.0;

    /* Penalize bond length deviations */
    for (int i = 0; i < mol->num_bonds; i++) {
        int a1 = mol->bonds[i].atom1;
        int a2 = mol->bonds[i].atom2;

        if (!placed[a1] || !placed[a2]) continue;

        double dist = point2d_distance(coords[a1], coords[a2]);
        double deviation = fabs(dist - bond_length) / bond_length;
        score += deviation * deviation * 10.0;
    }

    /* Penalize atom overlaps (non-bonded atoms too close) */
    double overlap_threshold = bond_length * 0.5;
    int overlaps = layout_count_atom_overlaps(ctx, overlap_threshold);
    score += overlaps * 100.0;

    /* Penalize bond crossings */
    int crossings = layout_count_bond_crossings(ctx);
    score += crossings * 50.0;

    return score;
}
