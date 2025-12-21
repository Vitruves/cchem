/**
 * @file coords2d.c
 * @brief 2D coordinate generation for molecular depiction
 *
 * Robust algorithm for fused ring systems including different-sized rings:
 * 1. Find all rings and group by shared atoms (not just edges)
 * 2. Place rings iteratively using placed atoms as anchors
 * 3. Handle spiro, fused, and bridged ring systems
 */

#include "cchem/depictor/coords2d.h"
#include "cchem/depictor/types.h"
#include "cchem/canonicalizer/bond.h"
#include "cchem/canonicalizer/element.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

/* Debug flag - set to 1 to enable debug output */
#define COORDS2D_DEBUG 0

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const coords2d_options_t DEFAULT_OPTIONS = COORDS2D_OPTIONS_DEFAULT;

/* ============== Ring Geometry ============== */

static double ring_circumradius(int size, double edge_length) {
    return edge_length / (2.0 * sin(M_PI / size));
}

/* ============== Ring Connectivity ============== */

/* Count shared atoms between two rings, optionally store them */
static int count_shared_atoms(const ring_t* r1, const ring_t* r2,
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

/* Check if rings share at least 2 atoms (fused) */
static bool rings_are_fused(const ring_t* r1, const ring_t* r2) {
    return count_shared_atoms(r1, r2, NULL, NULL, 10) >= 2;
}

/* ============== Union-Find for Ring Systems ============== */

static int uf_find(int* parent, int i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

static void uf_union(int* parent, int* rank, int a, int b) {
    int ra = uf_find(parent, a);
    int rb = uf_find(parent, b);
    if (ra == rb) return;
    if (rank[ra] < rank[rb]) parent[ra] = rb;
    else if (rank[ra] > rank[rb]) parent[rb] = ra;
    else { parent[rb] = ra; rank[ra]++; }
}

/* ============== Ring Placement ============== */

/* Get ring atoms in bond-connected order starting from a given atom */
static void get_ring_order(const ring_t* ring, const molecule_t* mol, int start_atom,
                           int* ordered, int* num_ordered) {
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

/* Place ring as regular polygon at given center, starting from specified atom (not index) */
static void place_ring_polygon_from(const ring_t* ring, const molecule_t* mol, point2d_t* coords, bool* placed,
                                    point2d_t center, double bond_length, double start_angle, int start_atom) {
    int n = ring->size;
    double radius = ring_circumradius(n, bond_length);
    double angle_step = 2.0 * M_PI / n;

    /* Get ring atoms in bond-connected order starting from start_atom */
    int* ordered = malloc(n * sizeof(int));
    int num_ordered = 0;
    get_ring_order(ring, mol, start_atom, ordered, &num_ordered);

    if (num_ordered < n) {
        /* Fallback: use ring->atoms[] directly */
        for (int i = 0; i < n; i++) {
            int atom = ring->atoms[i];
            double angle = start_angle - i * angle_step;
            coords[atom].x = center.x + radius * cos(angle);
            coords[atom].y = center.y + radius * sin(angle);
            placed[atom] = true;
        }
        free(ordered);
        return;
    }

    /* Place atoms in bond-connected order */
    for (int i = 0; i < n; i++) {
        int atom = ordered[i];
        double angle = start_angle - i * angle_step;
        coords[atom].x = center.x + radius * cos(angle);
        coords[atom].y = center.y + radius * sin(angle);
        placed[atom] = true;
    }

    free(ordered);
}

/* Place ring as regular polygon at given center */
static void place_ring_polygon(const ring_t* ring, const molecule_t* mol, point2d_t* coords, bool* placed,
                               point2d_t center, double bond_length, double start_angle) {
    place_ring_polygon_from(ring, mol, coords, placed, center, bond_length, start_angle, ring->atoms[0]);
}

/* Place a ring that shares atoms with already-placed atoms */
static bool place_fused_ring(const ring_t* ring, const molecule_t* mol,
                             point2d_t* coords, bool* placed, double bond_length) {
    int n = ring->size;

    /* Get ring atoms in bond-connected order */
    int* ordered = malloc(n * sizeof(int));
    int num_ordered = 0;
    get_ring_order(ring, mol, ring->atoms[0], ordered, &num_ordered);

    if (num_ordered < n) {
        free(ordered);
        return false;
    }

    /* Find placed atoms in the ordered ring */
    int placed_order_idx[20];  /* Index within ordered list */
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
    double radius = ring_circumradius(n, bond_length);
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

/* ============== Chain Placement ============== */

/* Trace a chain from a starting chain atom to find if it connects to another ring.
 * Returns the ring atom index at the other end, or -1 if chain is terminal.
 * Also returns the length of the chain in *chain_length if not NULL. */
static int trace_chain_to_ring(const molecule_t* mol, int start_chain_atom,
                                int from_ring_atom, int* chain_length) {
    if (chain_length) *chain_length = 0;

    int curr = start_chain_atom;
    int prev = from_ring_atom;
    int length = 1;

    while (true) {
        const atom_t* atom = &mol->atoms[curr];
        int next = -1;

        for (int i = 0; i < atom->num_neighbors; i++) {
            int nb = atom->neighbors[i];
            if (nb == prev) continue;

            if (mol->atoms[nb].ring_count > 0) {
                /* Found a ring at the other end */
                if (chain_length) *chain_length = length;
                return nb;
            }

            if (mol->atoms[nb].ring_count == 0) {
                /* Continue along chain */
                if (next < 0) {
                    next = nb;
                }
                /* If multiple chain neighbors, take the first non-prev one */
            }
        }

        if (next < 0) {
            /* Chain is terminal (no more chain atoms) */
            if (chain_length) *chain_length = length;
            return -1;
        }

        prev = curr;
        curr = next;
        length++;

        /* Safety: prevent infinite loop in case of cycles */
        if (length > mol->num_atoms) {
            if (chain_length) *chain_length = length;
            return -1;
        }
    }
}

static void place_chain_atoms(const molecule_t* mol, point2d_t* coords,
                               bool* placed, double bond_length) {
    int* queue = malloc(mol->num_atoms * sizeof(int));
    bool* in_queue = calloc(mol->num_atoms, sizeof(bool));
    int* chain_dir = calloc(mol->num_atoms, sizeof(int));
    int head = 0, tail = 0;

    /* Seed with placed atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (placed[i]) {
            queue[tail++] = i;
            in_queue[i] = true;
        }
    }

    /* If nothing placed, start with atom 0 */
    if (tail == 0 && mol->num_atoms > 0) {
        coords[0] = (point2d_t){0, 0};
        placed[0] = true;
        queue[tail++] = 0;
        in_queue[0] = true;
    }

    while (head < tail) {
        int curr = queue[head++];
        const atom_t* atom = &mol->atoms[curr];
        point2d_t curr_pos = coords[curr];

        /* Collect angles to placed neighbors */
        double angles[20];
        int n_angles = 0;

        for (int i = 0; i < atom->num_neighbors && n_angles < 20; i++) {
            int nb = atom->neighbors[i];
            if (placed[nb]) {
                point2d_t d = point2d_sub(coords[nb], curr_pos);
                angles[n_angles++] = atan2(d.y, d.x);
            }
        }

        /* Collect unplaced non-ring neighbors */
        int unplaced_nb[10];
        int n_unplaced = 0;
        for (int i = 0; i < atom->num_neighbors && n_unplaced < 10; i++) {
            int nb = atom->neighbors[i];
            if (in_queue[nb]) continue;
            if (mol->atoms[nb].ring_count > 0) continue;
            unplaced_nb[n_unplaced++] = nb;
        }

        /* For sp2 centers (carbonyl C), identify main chain vs side substituent */
        /* Main chain neighbor has more connections or is not terminal (O, N with H) */
        if (n_unplaced == 2 && n_angles == 1) {
            int main_nb = -1, side_nb = -1;
            for (int i = 0; i < n_unplaced; i++) {
                int nb = unplaced_nb[i];
                const atom_t* nb_atom = &mol->atoms[nb];
                /* Check if this is a terminal atom (O, S with no other heavy neighbors, or NH) */
                bool is_terminal = false;
                if (nb_atom->element == ELEM_O || nb_atom->element == ELEM_S) {
                    /* Count non-H neighbors */
                    int heavy_count = 0;
                    for (int j = 0; j < nb_atom->num_neighbors; j++) {
                        if (mol->atoms[nb_atom->neighbors[j]].element != ELEM_H) heavy_count++;
                    }
                    is_terminal = (heavy_count <= 1);
                }
                if (is_terminal) {
                    side_nb = nb;
                } else {
                    main_nb = nb;
                }
            }
            if (main_nb >= 0 && side_nb >= 0) {
                /* For sp2 centers (carbonyl), continue zigzag for main chain, place O perpendicular */
                double incoming = angles[0];
                int dir = -chain_dir[curr];
                if (dir == 0) {
                    /* When coming from ring or largest-gap placement, determine dir from geometry */
                    double norm_incoming = fmod(incoming + 2 * M_PI, 2 * M_PI);
                    dir = (norm_incoming > M_PI / 2 && norm_incoming < 3 * M_PI / 2) ? 1 : -1;
                }
                /* Main chain continues zigzag at 120° (same as normal chain placement) */
                double main_angle = incoming + M_PI + dir * M_PI / 3.0;
                /* Carbonyl O goes opposite side at 120° from main (perpendicular to chain) */
                double side_angle = incoming + M_PI - dir * M_PI / 3.0;

                coords[main_nb].x = curr_pos.x + bond_length * cos(main_angle);
                coords[main_nb].y = curr_pos.y + bond_length * sin(main_angle);
                placed[main_nb] = true;
                chain_dir[main_nb] = dir;
                queue[tail++] = main_nb;
                in_queue[main_nb] = true;

                coords[side_nb].x = curr_pos.x + bond_length * cos(side_angle);
                coords[side_nb].y = curr_pos.y + bond_length * sin(side_angle);
                placed[side_nb] = true;
                chain_dir[side_nb] = 0;
                queue[tail++] = side_nb;
                in_queue[side_nb] = true;

                continue;  /* Skip normal placement loop */
            }
        }

        /* Place unplaced neighbors (skip ring atoms - they'll be placed with ring geometry) */
        for (int i = 0; i < atom->num_neighbors; i++) {
            int nb = atom->neighbors[i];
            if (in_queue[nb]) continue;

            /* Skip atoms that are in rings - they should be placed with ring geometry */
            if (mol->atoms[nb].ring_count > 0) continue;

            /* If current atom is in a ring and nb is a chain atom, check if this chain
             * connects to another ring. If so, only the lower-indexed ring atom should
             * start placing the chain to ensure consistent zigzag direction. */
            if (mol->atoms[curr].ring_count > 0) {
                /* Trace the chain starting from nb to find if it connects to another ring */
                int other_ring_atom = trace_chain_to_ring(mol, nb, curr, NULL);

                if (other_ring_atom >= 0) {
                    /* Chain connects two rings. Only the lower-indexed ring atom
                     * should place the chain to ensure consistent zigzag. */
                    #if COORDS2D_DEBUG
                    fprintf(stderr, "DEBUG: Chain %d->%d connects to ring atom %d (curr=%d)\n",
                            curr, nb, other_ring_atom, curr);
                    #endif
                    if (curr > other_ring_atom) {
                        /* The other ring atom has lower index - let it handle this chain */
                        #if COORDS2D_DEBUG
                        fprintf(stderr, "DEBUG: Skipping chain %d->%d (deferring to atom %d)\n",
                                curr, nb, other_ring_atom);
                        #endif
                        continue;
                    }
                }
            }

            double angle;

            if (n_angles == 0) {
                angle = 0;
                chain_dir[nb] = 1;
            }
            else if (n_angles == 1) {
                /* Check if this is a triple bond (sp hybridization = linear 180°) */
                int bond_idx = atom->neighbor_bonds[i];
                bool is_triple = (bond_idx >= 0 && mol->bonds[bond_idx].type == BOND_TRIPLE);

                if (is_triple) {
                    /* Triple bond: place straight (180° from incoming) */
                    angle = angles[0] + M_PI;
                    chain_dir[nb] = 0;  /* No zigzag propagation */
                } else {
                    /* Zigzag: place at 120° internal angle (60° offset) from incoming */
                    int dir = -chain_dir[curr];
                    if (dir == 0) dir = 1;

                    /* Check if nb has an already-placed neighbor (e.g., connecting to ring).
                     * If so, we need to position nb so that BOTH the incoming bond (from curr)
                     * and the outgoing bond (to the placed neighbor) form a proper zigzag. */
                    const atom_t* nb_atom = &mol->atoms[nb];
                    int placed_neighbor = -1;
                    for (int k = 0; k < nb_atom->num_neighbors; k++) {
                        int nb_nb = nb_atom->neighbors[k];
                        if (nb_nb != curr && placed[nb_nb]) {
                            placed_neighbor = nb_nb;
                            break;
                        }
                    }

                    if (placed_neighbor >= 0) {
                        /* nb connects to an already-placed atom (likely a ring).
                         * Calculate standard zigzag angle, then check if it would
                         * create a near-parallel bond to the placed neighbor.
                         * If so, flip the direction to create visible zigzag. */
                        double standard_angle = angles[0] + M_PI + dir * M_PI / 3.0;

                        /* Calculate what the angle to placed_neighbor would be */
                        double nb_x = curr_pos.x + bond_length * cos(standard_angle);
                        double nb_y = curr_pos.y + bond_length * sin(standard_angle);
                        point2d_t placed_pos = coords[placed_neighbor];
                        double out_angle = atan2(placed_pos.y - nb_y, placed_pos.x - nb_x);

                        /* Normalize angles to [0, 2π] */
                        double in_norm = fmod(standard_angle + 4 * M_PI, 2 * M_PI);
                        double out_norm = fmod(out_angle + 4 * M_PI, 2 * M_PI);

                        /* Check if incoming and outgoing would be nearly parallel
                         * (within 30° of each other or 180° apart) */
                        double angle_diff = fabs(in_norm - out_norm);
                        if (angle_diff > M_PI) angle_diff = 2 * M_PI - angle_diff;

                        bool nearly_parallel = (angle_diff < M_PI / 6.0) ||
                                               (fabs(angle_diff - M_PI) < M_PI / 6.0);

                        if (nearly_parallel) {
                            /* Flip direction to create visible zigzag */
                            dir = -dir;
                            #if COORDS2D_DEBUG
                            fprintf(stderr, "DEBUG: ZigzagFlip %d->%d->%d: was parallel (diff=%.1f°), flipping dir to %d\n",
                                    curr, nb, placed_neighbor, angle_diff * 180.0 / M_PI, dir);
                            #endif
                        }

                        chain_dir[nb] = dir;
                        angle = angles[0] + M_PI + dir * M_PI / 3.0;
                        #if COORDS2D_DEBUG
                        fprintf(stderr, "DEBUG: Zigzag %d->%d: angles[0]=%.1f°, dir=%d, angle=%.1f°\n",
                                curr, nb, angles[0] * 180.0 / M_PI, dir, angle * 180.0 / M_PI);
                        #endif
                    } else {
                        /* Standard zigzag: no already-placed neighbor */
                        chain_dir[nb] = dir;
                        angle = angles[0] + M_PI + dir * M_PI / 3.0;
                        #if COORDS2D_DEBUG
                        fprintf(stderr, "DEBUG: Zigzag %d->%d: angles[0]=%.1f°, dir=%d, angle=%.1f°\n",
                                curr, nb, angles[0] * 180.0 / M_PI, dir, angle * 180.0 / M_PI);
                        #endif
                    }
                }
            }
            else {
                /* Find largest gap */
                for (int j = 0; j < n_angles - 1; j++) {
                    for (int k = j + 1; k < n_angles; k++) {
                        if (angles[j] > angles[k]) {
                            double tmp = angles[j];
                            angles[j] = angles[k];
                            angles[k] = tmp;
                        }
                    }
                }

                double max_gap = 0;
                int max_idx = 0;
                for (int j = 0; j < n_angles; j++) {
                    double next_angle = (j + 1 < n_angles) ? angles[j + 1] : angles[0] + 2 * M_PI;
                    double gap = next_angle - angles[j];
                    if (gap > max_gap) {
                        max_gap = gap;
                        max_idx = j;
                    }
                }

                angle = angles[max_idx] + max_gap / 2.0;
                /* Set chain_dir to create visible zigzag pattern.
                 * For best visual results, we want the zigzag to alternate symmetrically
                 * around the chain axis. The key is to start with a direction that
                 * will create clear left-right alternation.
                 *
                 * We use quadrant-based logic:
                 * - Q1 (0-90°): going up-right, zigzag should go down-right first (dir=-1)
                 * - Q2 (90-180°): going up-left, zigzag should go down-left first (dir=1)
                 * - Q3 (180-270°): going down-left, zigzag should go up-left first (dir=-1)
                 * - Q4 (270-360°): going down-right, zigzag should go up-right first (dir=1)
                 */
                double norm_angle = fmod(angle + 2 * M_PI, 2 * M_PI);

                /* Set chain_dir to create visible zigzag at chain start.
                 * For angles near vertical (up or down), we want the first
                 * zigzag to go more horizontal. For angles near horizontal,
                 * standard alternation works. */
                double vert = fabs(sin(norm_angle));  /* How vertical is this angle */

                if (vert > 0.5) {
                    /* Mostly vertical - zigzag should go more horizontal.
                     * For up-left (Q2) or down-right (Q4): use dir=1 to go right
                     * For up-right (Q1) or down-left (Q3): use dir=-1 to go left */
                    if ((norm_angle >= M_PI/2 && norm_angle < M_PI) ||
                        (norm_angle >= 3*M_PI/2)) {
                        chain_dir[nb] = -1;  /* Will create rightward zigzag */
                    } else {
                        chain_dir[nb] = 1;   /* Will create leftward zigzag */
                    }
                } else {
                    /* Mostly horizontal - standard alternation */
                    chain_dir[nb] = (norm_angle > M_PI) ? 1 : -1;
                }
                #if COORDS2D_DEBUG
                fprintf(stderr, "DEBUG: LargestGap %d->%d: angle=%.1f°, norm_angle=%.1f°, chain_dir=%d\n",
                        curr, nb, angle * 180.0 / M_PI, norm_angle * 180.0 / M_PI, chain_dir[nb]);
                #endif
            }

            coords[nb].x = curr_pos.x + bond_length * cos(angle);
            coords[nb].y = curr_pos.y + bond_length * sin(angle);
            placed[nb] = true;
            #if COORDS2D_DEBUG
            fprintf(stderr, "DEBUG: Placed atom %d (%s) at (%.2f, %.2f) from atom %d (%s) at (%.2f, %.2f), angle=%.1f°\n",
                    nb, element_to_symbol(mol->atoms[nb].element), coords[nb].x, coords[nb].y,
                    curr, element_to_symbol(mol->atoms[curr].element), curr_pos.x, curr_pos.y, angle * 180.0 / M_PI);
            #endif

            if (n_angles < 20) {
                angles[n_angles++] = angle;
            }

            queue[tail++] = nb;
            in_queue[nb] = true;
        }
    }

    free(queue);
    free(in_queue);
    free(chain_dir);
}

/* ============== Main API ============== */

mol_coords_t* coords2d_generate(const molecule_t* mol, const coords2d_options_t* options) {
    if (!mol || mol->num_atoms == 0) return NULL;

    coords2d_options_t opts = options ? *options : DEFAULT_OPTIONS;

    mol_coords_t* result = mol_coords_create(mol->num_atoms);
    if (!result) return NULL;
    result->has_2d = true;

    point2d_t* coords = result->coords_2d;
    bool* placed = calloc(mol->num_atoms, sizeof(bool));

    /* Handle rings */
    if (mol->rings_computed && mol->num_rings > 0) {
        int nr = mol->num_rings;

        /* Group rings into fused systems using union-find */
        int* parent = malloc(nr * sizeof(int));
        int* rank = calloc(nr, sizeof(int));
        for (int i = 0; i < nr; i++) parent[i] = i;

        for (int i = 0; i < nr; i++) {
            for (int j = i + 1; j < nr; j++) {
                if (rings_are_fused(&mol->rings[i], &mol->rings[j])) {
                    uf_union(parent, rank, i, j);
                }
            }
        }

        /* Find largest ring system */
        int* system_size = calloc(nr, sizeof(int));
        for (int i = 0; i < nr; i++) {
            int root = uf_find(parent, i);
            system_size[root] += mol->rings[i].size;
        }

        int largest_root = 0;
        for (int i = 0; i < nr; i++) {
            if (system_size[i] > system_size[largest_root]) {
                largest_root = i;
            }
        }

        /* Place rings in largest system first */
        bool* ring_placed = calloc(nr, sizeof(bool));

        /* Find most connected ring to start (most substituents) - consider ALL rings */
        int first_ring = -1;
        int max_connections = -1;
        for (int i = 0; i < nr; i++) {
            /* Count atoms in ring that have neighbors outside the ring */
            int connections = 0;
            for (int j = 0; j < mol->rings[i].size; j++) {
                int atom = mol->rings[i].atoms[j];
                for (int k = 0; k < mol->atoms[atom].num_neighbors; k++) {
                    int nb = mol->atoms[atom].neighbors[k];
                    /* Check if neighbor is not in this ring */
                    bool in_ring = false;
                    for (int m = 0; m < mol->rings[i].size; m++) {
                        if (mol->rings[i].atoms[m] == nb) { in_ring = true; break; }
                    }
                    if (!in_ring) connections++;
                }
            }
            if (connections > max_connections) {
                max_connections = connections;
                first_ring = i;
            }
        }

        if (first_ring >= 0) {
            place_ring_polygon(&mol->rings[first_ring], mol, coords, placed,
                              (point2d_t){0, 0}, opts.bond_length, M_PI / 2.0 + M_PI / mol->rings[first_ring].size);
            ring_placed[first_ring] = true;

            /* Iteratively place connected rings */
            bool progress = true;
            int max_iter = nr * nr;
            for (int iter = 0; iter < max_iter && progress; iter++) {
                progress = false;

                for (int r = 0; r < nr; r++) {
                    if (ring_placed[r]) continue;

                    /* Count placed atoms in this ring */
                    int placed_count = 0;
                    for (int i = 0; i < mol->rings[r].size; i++) {
                        if (placed[mol->rings[r].atoms[i]]) placed_count++;
                    }

                    if (placed_count >= 2) {
                        if (place_fused_ring(&mol->rings[r], mol, coords, placed, opts.bond_length)) {
                            ring_placed[r] = true;
                            progress = true;
                        }
                    }
                }
            }
        }

        /* Place chain atoms connecting to first ring system before processing remaining rings */
        place_chain_atoms(mol, coords, placed, opts.bond_length);

        /* Place remaining ring systems (connected via chains or disconnected) */
        bool remaining_progress = true;
        while (remaining_progress) {
            remaining_progress = false;
            for (int r = 0; r < nr; r++) {
                if (!ring_placed[r]) {
                    /* Check if any atoms already placed */
                    int placed_count = 0;
                    int placed_ring_idx = -1;
                    for (int i = 0; i < mol->rings[r].size; i++) {
                        if (placed[mol->rings[r].atoms[i]]) {
                            placed_count++;
                            placed_ring_idx = i;
                        }
                    }

                    if (placed_count >= 2) {
                        place_fused_ring(&mol->rings[r], mol, coords, placed, opts.bond_length);
                        ring_placed[r] = true;
                        remaining_progress = true;
                    } else if (placed_count == 1) {
                        /* Substituent ring - one atom placed (spiro), place ring around it */
                        const ring_t* ring = &mol->rings[r];
                        int anchor_atom = ring->atoms[placed_ring_idx];
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
                        double radius = opts.bond_length / (2.0 * sin(M_PI / n));
                        point2d_t center = point2d_add(anchor_pos, point2d_scale(away_dir, radius));

                        /* Place ring starting from anchor atom */
                        double start_angle = atan2(anchor_pos.y - center.y, anchor_pos.x - center.x);
                        place_ring_polygon_from(ring, mol, coords, placed, center, opts.bond_length, start_angle, anchor_atom);
                        ring_placed[r] = true;
                        remaining_progress = true;
                    } else if (placed_count == 0) {
                        /* Check if any ring atom has a placed neighbor (substituent ring) */
                        const ring_t* ring = &mol->rings[r];
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
                            /* Place ring as substituent connected to placed_neighbor */
                            int anchor_atom = ring->atoms[anchor_ring_idx];
                            point2d_t neighbor_pos = coords[placed_neighbor];

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
                                }
                            }
                            if (cnt > 0) {
                                double len = point2d_length(sum);
                                if (len > 0.01) {
                                    double base_angle = atan2(-sum.y, -sum.x);
                                    /* Only add zigzag offset when this ring is at the end of a
                                     * chain that connects to ANOTHER ring. Terminal substituent
                                     * rings (like cyclopropyl) should be placed straight. */
                                    if (mol->atoms[placed_neighbor].ring_count == 0) {
                                        /* Check if chain from placed_neighbor leads to another ring
                                         * through a proper chain (at least 2 atoms). Direct ring
                                         * attachments (like C=O to ring) should not trigger zigzag. */
                                        bool connects_via_chain = false;
                                        const atom_t* pn_atom = &mol->atoms[placed_neighbor];
                                        for (int k = 0; k < pn_atom->num_neighbors; k++) {
                                            int pn_nb = pn_atom->neighbors[k];
                                            if (pn_nb == anchor_atom) continue;
                                            if (mol->atoms[pn_nb].ring_count > 0) continue; /* Skip direct ring connections */
                                            /* Trace from this chain neighbor to see if it leads to a ring */
                                            int chain_len = 0;
                                            int other_ring = trace_chain_to_ring(mol, pn_nb, placed_neighbor, &chain_len);
                                            if (other_ring >= 0 && chain_len >= 2) {
                                                connects_via_chain = true;
                                                break;
                                            }
                                        }
                                        if (connects_via_chain) {
                                            double offset = M_PI / 3.0;  /* 60° zigzag offset */
                                            /* Alternate offset direction based on anchor atom index */
                                            if (anchor_atom % 2 == 0) offset = -offset;
                                            base_angle += offset;
                                        }
                                    }
                                    away_dir.x = cos(base_angle);
                                    away_dir.y = sin(base_angle);
                                }
                            }

                            /* Compute ring center position */
                            int n = ring->size;
                            double radius = opts.bond_length / (2.0 * sin(M_PI / n));
                            point2d_t anchor_pos = point2d_add(neighbor_pos,
                                                              point2d_scale(away_dir, opts.bond_length));
                            point2d_t center = point2d_add(anchor_pos, point2d_scale(away_dir, radius));

                            /* Check for collision with existing atoms - if any placed atom is within ring radius, flip */
                            bool collision = false;
                            for (int k = 0; k < mol->num_atoms; k++) {
                                if (!placed[k]) continue;
                                /* Skip the neighbor we're attaching to */
                                if (k == placed_neighbor) continue;
                                /* Check if not part of this ring */
                                bool in_this_ring = false;
                                for (int m = 0; m < ring->size; m++) {
                                    if (ring->atoms[m] == k) { in_this_ring = true; break; }
                                }
                                if (in_this_ring) continue;
                                double dist = point2d_distance(center, coords[k]);
                                /* If any atom is within ring radius of proposed center, collision */
                                if (dist < radius * 1.2) {
                                    collision = true;
                                    break;
                                }
                            }
                            if (collision) {
                                /* Try different angles to find collision-free placement */
                                double base_angle = atan2(away_dir.y, away_dir.x);
                                double test_angles[] = {M_PI, 2*M_PI/3, -2*M_PI/3, M_PI/3, -M_PI/3, M_PI/2, -M_PI/2};
                                int num_test = 7;
                                bool found = false;

                                for (int t = 0; t < num_test && !found; t++) {
                                    double test_angle = base_angle + test_angles[t];
                                    point2d_t test_dir = {cos(test_angle), sin(test_angle)};
                                    point2d_t test_anchor = point2d_add(neighbor_pos,
                                                                        point2d_scale(test_dir, opts.bond_length));
                                    point2d_t test_center = point2d_add(test_anchor, point2d_scale(test_dir, radius));

                                    bool test_collision = false;
                                    for (int k = 0; k < mol->num_atoms; k++) {
                                        if (!placed[k]) continue;
                                        if (k == placed_neighbor) continue;
                                        bool in_this_ring = false;
                                        for (int m = 0; m < ring->size; m++) {
                                            if (ring->atoms[m] == k) { in_this_ring = true; break; }
                                        }
                                        if (in_this_ring) continue;
                                        double dist = point2d_distance(test_center, coords[k]);
                                        if (dist < radius * 1.2) {
                                            test_collision = true;
                                            break;
                                        }
                                    }
                                    if (!test_collision) {
                                        away_dir = test_dir;
                                        anchor_pos = test_anchor;
                                        center = test_center;
                                        found = true;
                                    }
                                }
                            }

                            /* Place ring starting from anchor atom */
                            double start_angle = atan2(anchor_pos.y - center.y, anchor_pos.x - center.x);
                            place_ring_polygon_from(ring, mol, coords, placed, center, opts.bond_length, start_angle, anchor_atom);
                            ring_placed[r] = true;
                            remaining_progress = true;
                        }
                        /* Skip disconnected rings for now - will place in final pass */
                    }
                }
            }
            /* Place chain atoms after each round so next iteration can find neighbors */
            if (remaining_progress) {
                place_chain_atoms(mol, coords, placed, opts.bond_length);
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
                place_ring_polygon(&mol->rings[r], mol, coords, placed,
                                  (point2d_t){offset_x + 3 * opts.bond_length, 0},
                                  opts.bond_length, M_PI / 2.0);
                ring_placed[r] = true;
            }
        }

        free(parent);
        free(rank);
        free(system_size);
        free(ring_placed);
    }

    /* Place chain atoms */
    place_chain_atoms(mol, coords, placed, opts.bond_length);

    /* Handle disconnected atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!placed[i]) {
            coords[i] = (point2d_t){i * opts.bond_length, 5 * opts.bond_length};
            placed[i] = true;
        }
    }

    coords2d_center(result);
    free(placed);
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
    (void)options;
    if (!coords || !mol || !coords->has_2d) return -1;
    return 0;
}
