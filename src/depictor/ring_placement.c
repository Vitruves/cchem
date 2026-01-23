/**
 * @file ring_placement.c
 * @brief Ring placement algorithms for 2D molecular layout
 *
 * Provides algorithms for placing ring systems including template-based
 * placement, fused ring placement, spiro junctions, and substituent rings.
 */

#include "cchem/depictor/layout2d.h"
#include "cchem/depictor/templates.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* Debug output controlled by ctx->options->debug */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Forward declarations for collision detection */
extern bool layout_check_ring_collision(const layout_context_t* ctx,
                                        point2d_t center, double radius,
                                        const ring_t* ring, int exclude_atom);

extern bool layout_find_collision_free_angle(const layout_context_t* ctx,
                                             point2d_t neighbor_pos,
                                             const ring_t* ring,
                                             double base_angle, double radius,
                                             int exclude_atom,
                                             double* out_angle);

/* ============== Ring Polygon Placement ============== */

void layout_place_ring_polygon(layout_context_t* ctx, const ring_t* ring,
                               point2d_t center, double start_angle, int start_atom) {
    int n = ring->size;
    double radius = layout_ring_circumradius(n, ctx->bond_length);
    double angle_step = 2.0 * M_PI / n;

    /* Get ring atoms in bond-connected order starting from start_atom */
    int* ordered = malloc(n * sizeof(int));
    int num_ordered = 0;
    layout_get_ring_order(ring, ctx->mol, start_atom, ordered, &num_ordered);

    if (num_ordered < n) {
        /* Fallback: use ring->atoms[] directly */
        for (int i = 0; i < n; i++) {
            int atom = ring->atoms[i];
            double angle = start_angle - i * angle_step;
            ctx->coords->coords_2d[atom].x = center.x + radius * cos(angle);
            ctx->coords->coords_2d[atom].y = center.y + radius * sin(angle);
            ctx->placed[atom] = true;
        }
        free(ordered);
        return;
    }

    /* Place atoms in bond-connected order */
    for (int i = 0; i < n; i++) {
        int atom = ordered[i];
        double angle = start_angle - i * angle_step;
        ctx->coords->coords_2d[atom].x = center.x + radius * cos(angle);
        ctx->coords->coords_2d[atom].y = center.y + radius * sin(angle);
        ctx->placed[atom] = true;
    }

    free(ordered);
}

/* ============== Fused Ring Placement ============== */

bool layout_place_fused_ring(layout_context_t* ctx, const ring_t* ring) {
    int n = ring->size;
    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    double bond_length = ctx->bond_length;

    /* Get ring atoms in bond-connected order */
    int* ordered = malloc(n * sizeof(int));
    int num_ordered = 0;
    layout_get_ring_order(ring, mol, ring->atoms[0], ordered, &num_ordered);

    if (num_ordered < n) {
        free(ordered);
        return false;
    }

    /* Find placed atoms in the ordered ring */
    int placed_order_idx[20];
    int num_placed = 0;

    for (int i = 0; i < n && num_placed < 20; i++) {
        if (placed[ordered[i]]) {
            placed_order_idx[num_placed++] = i;
        }
    }

    if (num_placed < 2) {
        free(ordered);
        return false;
    }

    /* Find two adjacent placed atoms (consecutive in ordered list) to use as anchor */
    int anchor1_idx = -1, anchor2_idx = -1;

    for (int i = 0; i < num_placed; i++) {
        int curr_idx = placed_order_idx[i];
        for (int j = 0; j < num_placed; j++) {
            if (i == j) continue;
            int other_idx = placed_order_idx[j];

            /* Check if consecutive in ordered ring */
            if ((other_idx - curr_idx + n) % n == 1) {
                anchor1_idx = curr_idx;
                anchor2_idx = other_idx;
                break;
            }
        }
        if (anchor1_idx >= 0) break;
    }

    /* If no consecutive pair found, use first two */
    if (anchor1_idx < 0) {
        anchor1_idx = placed_order_idx[0];
        anchor2_idx = placed_order_idx[1];
    }

    int atom1 = ordered[anchor1_idx];
    int atom2 = ordered[anchor2_idx];
    point2d_t p1 = coords[atom1];
    point2d_t p2 = coords[atom2];

    /* Compute midpoint and perpendicular */
    point2d_t mid = point2d_scale(point2d_add(p1, p2), 0.5);
    point2d_t bond_vec = point2d_sub(p2, p1);
    double bond_len = point2d_length(bond_vec);
    if (bond_len < 0.001) {
        free(ordered);
        return false;
    }

    point2d_t bond_dir = point2d_scale(bond_vec, 1.0 / bond_len);
    point2d_t perp = (point2d_t){-bond_dir.y, bond_dir.x};

    /* Ring geometry - use ideal bond length for radius calculation */
    double radius = layout_ring_circumradius(n, bond_length);
    double half_chord = bond_len / 2.0;

    /* Clamp to avoid sqrt of negative */
    double apothem_sq = radius * radius - half_chord * half_chord;
    if (apothem_sq < 0.0001) apothem_sq = 0.0001;
    double apothem = sqrt(apothem_sq);

    /* Determine which side to place the new ring */
    point2d_t neighbor_sum = {0, 0};
    int neighbor_count = 0;

    for (int p = 0; p < num_placed; p++) {
        int atom = ordered[placed_order_idx[p]];
        const atom_t* a = &mol->atoms[atom];

        for (int j = 0; j < a->num_neighbors; j++) {
            int nb = a->neighbors[j];
            if (!placed[nb]) continue;

            /* Check if in this ring */
            bool in_ring = false;
            for (int k = 0; k < n; k++) {
                if (ordered[k] == nb) { in_ring = true; break; }
            }

            if (!in_ring) {
                neighbor_sum = point2d_add(neighbor_sum, coords[nb]);
                neighbor_count++;
            }
        }
    }

    double side = 1.0;
    if (neighbor_count > 0) {
        point2d_t neighbor_center = point2d_scale(neighbor_sum, 1.0 / neighbor_count);
        point2d_t to_neighbors = point2d_sub(neighbor_center, mid);
        side = (point2d_dot(to_neighbors, perp) > 0) ? -1.0 : 1.0;
    }

    point2d_t center = point2d_add(mid, point2d_scale(perp, side * apothem));

    /* Compute angle from center to anchor1 */
    point2d_t to_p1 = point2d_sub(p1, center);
    double angle_p1 = atan2(to_p1.y, to_p1.x);
    double angle_step = 2.0 * M_PI / n;

    /* Determine direction (CW or CCW) based on anchor positions */
    int steps_fwd = (anchor2_idx - anchor1_idx + n) % n;

    /* Try forward direction first */
    double test_angle = angle_p1 - steps_fwd * angle_step;
    point2d_t test_p2 = {
        center.x + radius * cos(test_angle),
        center.y + radius * sin(test_angle)
    };

    double dir = -1.0;  /* Default: counter-clockwise */
    if (point2d_distance(test_p2, p2) > bond_length * 0.5) {
        dir = 1.0;  /* Clockwise */
    }

    /* Place all unplaced atoms using ordered list */
    for (int i = 0; i < n; i++) {
        int atom = ordered[i];
        if (placed[atom]) continue;

        int steps = (i - anchor1_idx + n) % n;
        double angle = angle_p1 + dir * steps * angle_step;

        coords[atom].x = center.x + radius * cos(angle);
        coords[atom].y = center.y + radius * sin(angle);
        placed[atom] = true;
    }

    free(ordered);
    return true;
}

/* ============== Spiro Ring Placement ============== */

bool layout_place_spiro_ring(layout_context_t* ctx, const ring_t* ring, int anchor_atom) {
    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    double bond_length = ctx->bond_length;

    point2d_t anchor_pos = coords[anchor_atom];

    /* Find direction away from existing neighbors */
    point2d_t away_dir = {1, 0};
    const atom_t* anchor = &mol->atoms[anchor_atom];
    if (anchor->num_neighbors > 0) {
        point2d_t sum = {0, 0};
        int count = 0;
        for (int j = 0; j < anchor->num_neighbors; j++) {
            int nb = anchor->neighbors[j];
            if (placed[nb]) {
                sum = point2d_add(sum, point2d_sub(coords[nb], anchor_pos));
                count++;
            }
        }
        if (count > 0) {
            sum = point2d_scale(sum, 1.0 / count);
            double len = point2d_length(sum);
            if (len > 0.01) {
                away_dir = point2d_scale(sum, -1.0 / len);
            }
        }
    }

    /* Compute ring center position */
    int n = ring->size;
    double radius = layout_ring_circumradius(n, bond_length);
    point2d_t center = point2d_add(anchor_pos, point2d_scale(away_dir, radius));

    /* Place ring starting from anchor atom */
    double start_angle = atan2(anchor_pos.y - center.y, anchor_pos.x - center.x);
    layout_place_ring_polygon(ctx, ring, center, start_angle, anchor_atom);

    return true;
}

/* ============== Substituent Ring Placement ============== */

bool layout_place_substituent_ring(layout_context_t* ctx, const ring_t* ring,
                                   int anchor_atom, int placed_neighbor) {
    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    double bond_length = ctx->bond_length;

    point2d_t neighbor_pos = coords[placed_neighbor];

if (ctx->options->debug) {
    fprintf(stderr, "DEBUG substituent_ring: Placing ring with anchor=%d, neighbor=%d at (%.2f, %.2f)\n",
            anchor_atom, placed_neighbor, neighbor_pos.x, neighbor_pos.y);
}

    /* Find direction away from neighbor's other connections */
    point2d_t away_dir = {1, 0};
    const atom_t* nb_atom = &mol->atoms[placed_neighbor];
    point2d_t sum = {0, 0};
    int cnt = 0;

    for (int j = 0; j < nb_atom->num_neighbors; j++) {
        int other = nb_atom->neighbors[j];
        if (placed[other] && other != anchor_atom) {
            sum = point2d_add(sum, point2d_sub(coords[other], neighbor_pos));
            cnt++;
if (ctx->options->debug) {
            fprintf(stderr, "  neighbor %d's other connection: atom %d at (%.2f, %.2f)\n",
                    placed_neighbor, other, coords[other].x, coords[other].y);
}
        }
    }

    if (cnt > 0) {
        double len = point2d_length(sum);
        if (len > 0.01) {
            /* Calculate base angle pointing away from other connections */
            double base_angle = atan2(-sum.y, -sum.x);

            /* Only apply zigzag offset when the neighbor is a pure sp3 chain atom.
             * Don't apply zigzag for:
             * - Ring atoms (preserves ring-to-ring straight connections)
             * - sp2 atoms with double bonds (preserves planar amide/carbonyl geometry) */
            bool is_sp3_chain = (mol->atoms[placed_neighbor].ring_count == 0);

            if (is_sp3_chain) {
                /* Check if neighbor has any double/triple bonds (sp2/sp hybridization) */
                const atom_t* nb_atom_check = &mol->atoms[placed_neighbor];
                for (int k = 0; k < nb_atom_check->num_neighbors; k++) {
                    int bond_idx = nb_atom_check->neighbor_bonds[k];
                    if (bond_idx >= 0 &&
                        (mol->bonds[bond_idx].type == BOND_DOUBLE ||
                         mol->bonds[bond_idx].type == BOND_TRIPLE)) {
                        is_sp3_chain = false;
                        break;
                    }
                }
            }

            if (is_sp3_chain) {
                /* Neighbor is a pure sp3 chain atom - apply zigzag offset (60 degrees) */
                int neighbor_chain_dir = ctx->chain_dir[placed_neighbor];
                int zigzag_dir = -neighbor_chain_dir;
                if (zigzag_dir == 0) zigzag_dir = 1;

                base_angle = base_angle + zigzag_dir * M_PI / 3.0;
            }

            away_dir.x = cos(base_angle);
            away_dir.y = sin(base_angle);
        }
    }

if (ctx->options->debug) {
    fprintf(stderr, "  away_dir = (%.2f, %.2f), angle = %.1f deg\n",
            away_dir.x, away_dir.y, atan2(away_dir.y, away_dir.x) * 180 / M_PI);
}

    /* Compute ring center position */
    int n = ring->size;
    double radius = layout_ring_circumradius(n, bond_length);
    point2d_t anchor_pos = point2d_add(neighbor_pos, point2d_scale(away_dir, bond_length));
    point2d_t center = point2d_add(anchor_pos, point2d_scale(away_dir, radius));

if (ctx->options->debug) {
    fprintf(stderr, "  Initial: anchor_pos=(%.2f, %.2f), center=(%.2f, %.2f), radius=%.2f\n",
            anchor_pos.x, anchor_pos.y, center.x, center.y, radius);
}

    /* Check for collision with existing atoms */
    bool collision = layout_check_ring_collision(ctx, center, radius, ring, placed_neighbor);

if (ctx->options->debug) {
    fprintf(stderr, "  Collision check result: %s\n", collision ? "YES" : "NO");
}

    if (collision) {
        /* Try to find collision-free angle */
        double base_angle = atan2(away_dir.y, away_dir.x);
        double out_angle;
        if (layout_find_collision_free_angle(ctx, neighbor_pos, ring,
                                             base_angle, radius, placed_neighbor, &out_angle)) {
if (ctx->options->debug) {
            fprintf(stderr, "  Found collision-free angle: %.1f deg (was %.1f deg)\n",
                    out_angle * 180 / M_PI, base_angle * 180 / M_PI);
}
            away_dir.x = cos(out_angle);
            away_dir.y = sin(out_angle);
            anchor_pos = point2d_add(neighbor_pos, point2d_scale(away_dir, bond_length));
            center = point2d_add(anchor_pos, point2d_scale(away_dir, radius));
        } else {
if (ctx->options->debug) {
            fprintf(stderr, "  WARNING: No collision-free angle found!\n");
}
        }
    }

    /* Place anchor atom first */
    coords[anchor_atom] = anchor_pos;
    placed[anchor_atom] = true;

    /* Place ring starting from anchor atom */
    double start_angle = atan2(anchor_pos.y - center.y, anchor_pos.x - center.x);
    layout_place_ring_polygon(ctx, ring, center, start_angle, anchor_atom);

if (ctx->options->debug) {
    fprintf(stderr, "  Final positions:\n");
    for (int i = 0; i < ring->size; i++) {
        int atom = ring->atoms[i];
        fprintf(stderr, "    atom %d: (%.2f, %.2f)\n", atom, coords[atom].x, coords[atom].y);
    }
}

    return true;
}

/* ============== Find Best Starting Ring ============== */

static int find_most_connected_ring(const molecule_t* mol) {
    int nr = mol->num_rings;
    int first_ring = -1;
    int best_score = -1;

    /*
     * For complex polycyclic systems (especially bridged), we need to start with
     * a ring that has well-defined geometry. Priority:
     * 1. 6-membered aromatic rings (benzene) - most stable geometry
     * 2. 6-membered rings - good geometry
     * 3. 5-membered aromatic rings
     * 4. Other rings by connectivity
     *
     * Score: aromatic_bonus + size_bonus + connection_count
     */
    for (int i = 0; i < nr; i++) {
        const ring_t* ring = &mol->rings[i];
        int score = 0;

        /* Check if ring is aromatic (all atoms aromatic) */
        bool is_aromatic = true;
        for (int j = 0; j < ring->size; j++) {
            if (!mol->atoms[ring->atoms[j]].aromatic) {
                is_aromatic = false;
                break;
            }
        }

        /* Aromatic rings get high priority */
        if (is_aromatic) {
            score += 1000;
        }

        /* Prefer 6-membered rings (most common, best geometry) */
        if (ring->size == 6) {
            score += 500;
        } else if (ring->size == 5) {
            score += 200;
        }

        /* Add connectivity as tiebreaker */
        int connections = 0;
        for (int j = 0; j < ring->size; j++) {
            int atom = ring->atoms[j];
            for (int k = 0; k < mol->atoms[atom].num_neighbors; k++) {
                int nb = mol->atoms[atom].neighbors[k];
                bool in_ring = false;
                for (int m = 0; m < ring->size; m++) {
                    if (ring->atoms[m] == nb) { in_ring = true; break; }
                }
                if (!in_ring) connections++;
            }
        }
        score += connections;

        if (score > best_score) {
            best_score = score;
            first_ring = i;
        }
    }

    return first_ring;
}

/* ============== Main Ring System Placement ============== */

/* Forward declaration for chain placement (needed for interleaving) */
extern int layout_place_chains(layout_context_t* ctx);

int layout_place_ring_systems(layout_context_t* ctx) {
    if (!ctx || !ctx->mol) return -1;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    double bond_length = ctx->bond_length;

    /* Handle case with no rings */
    if (!mol->rings_computed || mol->num_rings == 0) {
        return 0;
    }

    int nr = mol->num_rings;
    bool* ring_placed = calloc(nr, sizeof(bool));

    /* Template matching disabled for now - templates don't handle inter-system connectivity.
     * TODO: Only use templates for the FIRST ring system with no placed atoms,
     * then use algorithmic placement for subsequent systems connected via chains.
     */
#if 0
    /* Try template-based placement for ring systems first */
    if (ctx->options && ctx->options->use_templates && ctx->ring_systems) {
if (ctx->options->debug) {
        fprintf(stderr, "DEBUG: Trying template matching for %d ring systems\n", ctx->num_ring_systems);
}
        for (int s = 0; s < ctx->num_ring_systems; s++) {
            ring_system_t* sys = &ctx->ring_systems[s];
            const ring_template_t* templ = template_find_match(ctx, sys);

            if (templ) {
if (ctx->options->debug) {
                fprintf(stderr, "DEBUG: Found template '%s' for system %d\n", templ->name, s);
}
                point2d_t center;
                double rotation;
                if (template_find_transformation(ctx, sys, templ, &center, &rotation)) {
                    if (template_apply(ctx, sys, templ, center, bond_length, rotation)) {
if (ctx->options->debug) {
                        fprintf(stderr, "DEBUG: Applied template to system %d\n", s);
}
                        sys->placed = true;
                        /* Mark all rings in system as placed */
                        for (int r = 0; r < sys->num_rings; r++) {
                            ring_placed[sys->ring_indices[r]] = true;
                        }
                    }
                }
            }
        }
    }
#endif

    /* Find most connected ring to start (if not already placed) */
    int first_ring = find_most_connected_ring(mol);

if (ctx->options->debug) {
    fprintf(stderr, "DEBUG ring_placement: first_ring=%d, ring_placed states: ", first_ring);
    for (int i = 0; i < nr; i++) {
        fprintf(stderr, "%d:%s ", i, ring_placed[i] ? "T" : "F");
    }
    fprintf(stderr, "\n");
}

    if (first_ring >= 0 && !ring_placed[first_ring]) {
        double start_angle = M_PI / 2.0 + M_PI / mol->rings[first_ring].size;
        layout_place_ring_polygon(ctx, &mol->rings[first_ring],
                                  (point2d_t){0, 0}, start_angle,
                                  mol->rings[first_ring].atoms[0]);
        ring_placed[first_ring] = true;
if (ctx->options->debug) {
        fprintf(stderr, "DEBUG: Placed first ring %d at origin\n", first_ring);
}
    }

    /* Iteratively place fused rings in the first ring system */
    bool progress = true;
    int max_iter = nr * nr;
    for (int iter = 0; iter < max_iter && progress; iter++) {
        progress = false;

        for (int r = 0; r < nr; r++) {
            if (ring_placed[r]) continue;

            /* Count placed atoms in this ring */
            int placed_count = 0;
            for (int i = 0; i < mol->rings[r].size; i++) {
                if (placed[mol->rings[r].atoms[i]]) {
                    placed_count++;
                }
            }

            if (placed_count >= 2) {
                /* Fused ring - use existing placement */
if (ctx->options->debug) {
                fprintf(stderr, "DEBUG: Placing fused ring %d (placed_count=%d)\n", r, placed_count);
}
                if (layout_place_fused_ring(ctx, &mol->rings[r])) {
                    ring_placed[r] = true;
                    progress = true;
                }
            }
        }
    }

if (ctx->options->debug) {
    fprintf(stderr, "DEBUG: After fused ring placement, calling layout_place_chains\n");
}

    /* CRITICAL: Place chain atoms connecting to first ring system before processing remaining rings
     * This allows rings connected via chains to find placed neighbors */
    layout_place_chains(ctx);

    /* Place remaining ring systems (connected via chains or disconnected) */
    bool remaining_progress = true;
    while (remaining_progress) {
        remaining_progress = false;

        for (int r = 0; r < nr; r++) {
            if (ring_placed[r]) continue;

            const ring_t* ring = &mol->rings[r];

            /* Count placed atoms in this ring */
            int placed_count = 0;
            int placed_ring_idx = -1;
            for (int i = 0; i < ring->size; i++) {
                if (placed[ring->atoms[i]]) {
                    placed_count++;
                    placed_ring_idx = i;
                }
            }

            if (placed_count >= 2) {
                /* Fused ring */
if (ctx->options->debug) {
                fprintf(stderr, "DEBUG remaining: Placing fused ring %d\n", r);
}
                if (layout_place_fused_ring(ctx, ring)) {
                    ring_placed[r] = true;
                    remaining_progress = true;
                }
            } else if (placed_count == 1) {
                /* Spiro ring - one atom placed, place ring around it */
                int anchor_atom = ring->atoms[placed_ring_idx];
if (ctx->options->debug) {
                fprintf(stderr, "DEBUG remaining: Placing spiro ring %d, anchor=%d\n", r, anchor_atom);
}
                if (layout_place_spiro_ring(ctx, ring, anchor_atom)) {
                    ring_placed[r] = true;
                    remaining_progress = true;
                }
            } else if (placed_count == 0) {
                /* Check if any ring atom has a placed neighbor (substituent ring) */
                int anchor_ring_idx = -1;
                int placed_neighbor = -1;

                for (int i = 0; i < ring->size && anchor_ring_idx < 0; i++) {
                    int atom = ring->atoms[i];
                    const atom_t* a = &mol->atoms[atom];
                    for (int j = 0; j < a->num_neighbors; j++) {
                        int nb = a->neighbors[j];
                        if (placed[nb]) {
                            anchor_ring_idx = i;
                            placed_neighbor = nb;
                            break;
                        }
                    }
                }

                if (anchor_ring_idx >= 0) {
                    int anchor_atom = ring->atoms[anchor_ring_idx];
if (ctx->options->debug) {
                    fprintf(stderr, "DEBUG remaining: Placing substituent ring %d, anchor=%d, neighbor=%d\n",
                            r, anchor_atom, placed_neighbor);
}
                    if (layout_place_substituent_ring(ctx, ring,
                                                      anchor_atom, placed_neighbor)) {
                        ring_placed[r] = true;
                        remaining_progress = true;
                    }
                } else {
if (ctx->options->debug) {
                    fprintf(stderr, "DEBUG remaining: Ring %d (atoms:", r);
                    for (int i = 0; i < ring->size; i++) {
                        fprintf(stderr, " %d", ring->atoms[i]);
                    }
                    fprintf(stderr, ") has no placed neighbor - skipping\n");
}
                }
                /* Skip disconnected rings for now - will place in final pass */
            }
        }

        /* Place chain atoms after each round so next iteration can find neighbors */
        if (remaining_progress) {
            layout_place_chains(ctx);
        }
    }

    /* Final pass: place any remaining disconnected rings */
    for (int r = 0; r < nr; r++) {
        if (!ring_placed[r]) {
            double offset_x = 0;
            for (int i = 0; i < mol->num_atoms; i++) {
                if (placed[i] && coords[i].x > offset_x) {
                    offset_x = coords[i].x;
                }
            }
            double start_angle = M_PI / 2.0 + M_PI / mol->rings[r].size;
            layout_place_ring_polygon(ctx, &mol->rings[r],
                                      (point2d_t){offset_x + 3 * bond_length, 0},
                                      start_angle, mol->rings[r].atoms[0]);
            ring_placed[r] = true;
        }
    }

    free(ring_placed);
    return 0;
}
