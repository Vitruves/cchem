/**
 * @file chain_placement.c
 * @brief Chain and substituent placement for 2D molecular layout
 *
 * Handles placement of non-ring atoms using BFS traversal with
 * zigzag angle calculation for optimal visual appearance.
 */

#include "cchem/depictor/layout2d.h"
#include "cchem/canonicalizer/bond.h"
#include "cchem/canonicalizer/element.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* Debug output controlled by ctx->options->debug */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============== Chain Tracing ============== */

int layout_trace_chain_to_ring(const molecule_t* mol, int start_chain_atom,
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

/* ============== Carbonyl/Functional Group Detection ============== */

static bool is_terminal_atom(const molecule_t* mol, int atom_idx, int from_atom) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Check if this is a terminal atom (O, S with no other heavy neighbors, or NH) */
    if (atom->element == ELEM_O || atom->element == ELEM_S) {
        /* Count non-H neighbors */
        int heavy_count = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int nb = atom->neighbors[j];
            if (nb == from_atom) continue;
            if (mol->atoms[nb].element != ELEM_H) heavy_count++;
        }
        return (heavy_count == 0);
    }

    return false;
}

/* ============== Angle Calculation ============== */

static double calculate_zigzag_angle(const layout_context_t* ctx, int curr, int nb,
                                     double incoming_angle, int* chain_dir) {
    const molecule_t* mol = ctx->mol;
    int curr_dir = ctx->chain_dir[curr];

    /* Check if this is a triple bond (sp hybridization = linear 180 deg) */
    const atom_t* curr_atom = &mol->atoms[curr];
    for (int i = 0; i < curr_atom->num_neighbors; i++) {
        if (curr_atom->neighbors[i] == nb) {
            int bond_idx = curr_atom->neighbor_bonds[i];
            if (bond_idx >= 0 && mol->bonds[bond_idx].type == BOND_TRIPLE) {
                /* Triple bond: place straight (180 deg from incoming) */
                *chain_dir = 0;
                return incoming_angle + M_PI;
            }
            break;
        }
    }

    /* Standard zigzag: place at 120 deg internal angle */
    int dir = -curr_dir;
    if (dir == 0) dir = 1;

    /* Check if nb has an already-placed neighbor (e.g., connecting to ring) */
    const atom_t* nb_atom = &mol->atoms[nb];
    int placed_neighbor = -1;
    for (int k = 0; k < nb_atom->num_neighbors; k++) {
        int nb_nb = nb_atom->neighbors[k];
        if (nb_nb != curr && ctx->placed[nb_nb]) {
            placed_neighbor = nb_nb;
            break;
        }
    }

    if (placed_neighbor >= 0) {
        /* nb connects to an already-placed atom (likely a ring) */
        double standard_angle = incoming_angle + M_PI + dir * M_PI / 3.0;

        /* Calculate what the angle to placed_neighbor would be */
        double nb_x = ctx->coords->coords_2d[curr].x + ctx->bond_length * cos(standard_angle);
        double nb_y = ctx->coords->coords_2d[curr].y + ctx->bond_length * sin(standard_angle);
        point2d_t placed_pos = ctx->coords->coords_2d[placed_neighbor];
        double out_angle = atan2(placed_pos.y - nb_y, placed_pos.x - nb_x);

        /* Normalize angles to [0, 2pi] */
        double in_norm = fmod(standard_angle + 4 * M_PI, 2 * M_PI);
        double out_norm = fmod(out_angle + 4 * M_PI, 2 * M_PI);

        /* Check if incoming and outgoing would be nearly parallel */
        double angle_diff = fabs(in_norm - out_norm);
        if (angle_diff > M_PI) angle_diff = 2 * M_PI - angle_diff;

        bool nearly_parallel = (angle_diff < M_PI / 6.0) ||
                               (fabs(angle_diff - M_PI) < M_PI / 6.0);

        if (nearly_parallel) {
            dir = -dir;
        }
    }

    *chain_dir = dir;
    return incoming_angle + M_PI + dir * M_PI / 3.0;
}

static double calculate_largest_gap_angle(const double* angles, int n_angles,
                                          int* out_chain_dir) {
    /* Sort angles */
    double sorted[20];
    for (int i = 0; i < n_angles && i < 20; i++) {
        sorted[i] = angles[i];
    }
    for (int i = 0; i < n_angles - 1; i++) {
        for (int j = i + 1; j < n_angles; j++) {
            if (sorted[i] > sorted[j]) {
                double tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }

    /* Find largest gap */
    double max_gap = 0;
    int max_idx = 0;
    for (int j = 0; j < n_angles; j++) {
        double next_angle = (j + 1 < n_angles) ? sorted[j + 1] : sorted[0] + 2 * M_PI;
        double gap = next_angle - sorted[j];
        if (gap > max_gap) {
            max_gap = gap;
            max_idx = j;
        }
    }

    double angle = sorted[max_idx] + max_gap / 2.0;

    /* Set chain direction based on angle quadrant */
    double norm_angle = fmod(angle + 2 * M_PI, 2 * M_PI);
    double vert = fabs(sin(norm_angle));

    if (vert > 0.5) {
        /* Mostly vertical */
        if ((norm_angle >= M_PI/2 && norm_angle < M_PI) ||
            (norm_angle >= 3*M_PI/2)) {
            *out_chain_dir = -1;
        } else {
            *out_chain_dir = 1;
        }
    } else {
        /* Mostly horizontal */
        *out_chain_dir = (norm_angle > M_PI) ? 1 : -1;
    }

    return angle;
}

/* ============== Main Chain Placement ============== */

int layout_place_chains(layout_context_t* ctx) {
    if (!ctx || !ctx->mol) return -1;

    const molecule_t* mol = ctx->mol;
    point2d_t* coords = ctx->coords->coords_2d;
    bool* placed = ctx->placed;
    int* chain_dir = ctx->chain_dir;
    double bond_length = ctx->bond_length;

    int* queue = malloc(mol->num_atoms * sizeof(int));
    bool* in_queue = calloc(mol->num_atoms, sizeof(bool));
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

        /* Also collect unplaced ring neighbors (for sp2 handling of carbonyls connected to rings) */
        int unplaced_ring_nb[10];
        int n_unplaced_ring = 0;
        for (int i = 0; i < atom->num_neighbors && n_unplaced_ring < 10; i++) {
            int nb = atom->neighbors[i];
            if (in_queue[nb]) continue;
            if (mol->atoms[nb].ring_count == 0) continue;
            unplaced_ring_nb[n_unplaced_ring++] = nb;
        }

        /* Handle sp2 centers (carbonyl C) - identify main chain vs side substituent
         * This handles cases where:
         * 1. Two chain neighbors: one terminal (O/S), one main chain
         * 2. One chain neighbor (terminal O/S) + one ring neighbor (N in ring) */
        if (n_unplaced == 2 && n_angles == 1) {
            int main_nb = -1, side_nb = -1;
            for (int i = 0; i < n_unplaced; i++) {
                int nb = unplaced_nb[i];
                if (is_terminal_atom(mol, nb, curr)) {
                    side_nb = nb;
                } else {
                    main_nb = nb;
                }
            }
            if (main_nb >= 0 && side_nb >= 0) {
                double incoming = angles[0];
                int dir = -chain_dir[curr];
                if (dir == 0) {
                    double norm_incoming = fmod(incoming + 2 * M_PI, 2 * M_PI);
                    dir = (norm_incoming > M_PI / 2 && norm_incoming < 3 * M_PI / 2) ? 1 : -1;
                }

                /* Main chain continues zigzag */
                double main_angle = incoming + M_PI + dir * M_PI / 3.0;
                /* Carbonyl O goes opposite side */
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

                continue;
            }
        }

        /* Handle sp2 centers connected to ring (e.g., amide C=O connected to ring N)
         * Case: 1 chain neighbor (terminal O/S) + 1 unplaced ring neighbor */
        if (n_unplaced == 1 && n_unplaced_ring == 1 && n_angles == 1) {
            int chain_nb = unplaced_nb[0];
            if (is_terminal_atom(mol, chain_nb, curr)) {
                /* This is a carbonyl connected to a ring atom */
                double incoming = angles[0];
                int dir = -chain_dir[curr];
                if (dir == 0) {
                    double norm_incoming = fmod(incoming + 2 * M_PI, 2 * M_PI);
                    dir = (norm_incoming > M_PI / 2 && norm_incoming < 3 * M_PI / 2) ? 1 : -1;
                }

                /* Carbonyl O goes to zigzag side (will be visible) */
                double side_angle = incoming + M_PI - dir * M_PI / 3.0;

                coords[chain_nb].x = curr_pos.x + bond_length * cos(side_angle);
                coords[chain_nb].y = curr_pos.y + bond_length * sin(side_angle);
                placed[chain_nb] = true;
                chain_dir[chain_nb] = 0;
                queue[tail++] = chain_nb;
                in_queue[chain_nb] = true;

                /* Store the "main chain" direction for the ring to use later */
                /* The ring will be placed by ring_placement, but we set chain_dir */
                chain_dir[curr] = dir;

                continue;
            }
        }

        /* Place unplaced neighbors */
        for (int i = 0; i < atom->num_neighbors; i++) {
            int nb = atom->neighbors[i];
            if (ctx->options->debug && (nb == 6 || curr == 7)) {
                fprintf(stderr, "DEBUG chain: curr=%d checking neighbor %d: in_queue=%d ring_count=%d placed=%d\n",
                        curr, nb, in_queue[nb], mol->atoms[nb].ring_count, placed[nb]);
            }
            if (in_queue[nb]) continue;
            if (mol->atoms[nb].ring_count > 0) continue;

            /* Check if chain connects to another ring via lower-indexed atom.
             * Only defer if the other ring atom is ALREADY PLACED - otherwise
             * we need to place the chain so the other ring can detect a placed neighbor. */
            if (mol->atoms[curr].ring_count > 0) {
                int other_ring_atom = layout_trace_chain_to_ring(mol, nb, curr, NULL);
                if (other_ring_atom >= 0 && curr > other_ring_atom && placed[other_ring_atom]) {
                    continue; /* Let other ring atom handle this chain */
                }
            }

            double angle;
            int new_chain_dir;

            if (n_angles == 0) {
                angle = 0;
                new_chain_dir = 1;
            } else if (n_angles == 1) {
                angle = calculate_zigzag_angle(ctx, curr, nb, angles[0], &new_chain_dir);
            } else {
                angle = calculate_largest_gap_angle(angles, n_angles, &new_chain_dir);
            }

            coords[nb].x = curr_pos.x + bond_length * cos(angle);
            coords[nb].y = curr_pos.y + bond_length * sin(angle);
            placed[nb] = true;
            chain_dir[nb] = new_chain_dir;

            if (n_angles < 20) {
                angles[n_angles++] = angle;
            }

            queue[tail++] = nb;
            in_queue[nb] = true;
        }
    }

    free(queue);
    free(in_queue);

    return 0;
}
