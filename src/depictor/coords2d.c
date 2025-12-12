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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

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

/* Place ring as regular polygon at given center */
static void place_ring_polygon(const ring_t* ring, point2d_t* coords, bool* placed,
                               point2d_t center, double bond_length, double start_angle) {
    int n = ring->size;
    double radius = ring_circumradius(n, bond_length);
    double angle_step = 2.0 * M_PI / n;

    for (int i = 0; i < n; i++) {
        double angle = start_angle - i * angle_step;
        coords[ring->atoms[i]].x = center.x + radius * cos(angle);
        coords[ring->atoms[i]].y = center.y + radius * sin(angle);
        placed[ring->atoms[i]] = true;
    }
}

/* Place a ring that shares atoms with already-placed atoms */
static bool place_fused_ring(const ring_t* ring, const molecule_t* mol,
                             point2d_t* coords, bool* placed, double bond_length) {
    int n = ring->size;

    /* Find placed atoms in this ring */
    int placed_ring_idx[20];  /* Index within ring */
    int placed_atom_idx[20];  /* Atom index */
    int num_placed = 0;

    for (int i = 0; i < n && num_placed < 20; i++) {
        if (placed[ring->atoms[i]]) {
            placed_ring_idx[num_placed] = i;
            placed_atom_idx[num_placed] = ring->atoms[i];
            num_placed++;
        }
    }

    if (num_placed < 2) return false;

    /* Find two adjacent placed atoms to use as anchor */
    int anchor1_ring = -1, anchor2_ring = -1;

    for (int i = 0; i < num_placed; i++) {
        int curr_ring = placed_ring_idx[i];
        for (int j = 0; j < num_placed; j++) {
            if (i == j) continue;
            int other_ring = placed_ring_idx[j];

            /* Check if consecutive in ring */
            if ((other_ring - curr_ring + n) % n == 1) {
                anchor1_ring = curr_ring;
                anchor2_ring = other_ring;
                break;
            }
        }
        if (anchor1_ring >= 0) break;
    }

    /* If no consecutive pair, use first two */
    if (anchor1_ring < 0) {
        anchor1_ring = placed_ring_idx[0];
        anchor2_ring = placed_ring_idx[1];
    }

    int atom1 = ring->atoms[anchor1_ring];
    int atom2 = ring->atoms[anchor2_ring];
    point2d_t p1 = coords[atom1];
    point2d_t p2 = coords[atom2];

    /* Compute midpoint and perpendicular */
    point2d_t mid = point2d_scale(point2d_add(p1, p2), 0.5);
    point2d_t bond_vec = point2d_sub(p2, p1);
    double bond_len = point2d_length(bond_vec);
    if (bond_len < 0.001) return false;

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
    /* Look at neighbors of anchor atoms that are placed but not in this ring */
    point2d_t neighbor_sum = {0, 0};
    int neighbor_count = 0;

    for (int p = 0; p < num_placed; p++) {
        int atom = placed_atom_idx[p];
        const atom_t* a = &mol->atoms[atom];

        for (int j = 0; j < a->num_neighbors; j++) {
            int nb = a->neighbors[j];
            if (!placed[nb]) continue;

            /* Check if in this ring */
            bool in_ring = false;
            for (int k = 0; k < n; k++) {
                if (ring->atoms[k] == nb) { in_ring = true; break; }
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

    /* Determine direction (CW or CCW) */
    /* From anchor1 to anchor2 in ring indices */
    int steps_fwd = (anchor2_ring - anchor1_ring + n) % n;

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

    /* Place all unplaced atoms */
    for (int i = 0; i < n; i++) {
        int atom = ring->atoms[i];
        if (placed[atom]) continue;

        int steps = (i - anchor1_ring + n) % n;
        double angle = angle_p1 + dir * steps * angle_step;

        coords[atom].x = center.x + radius * cos(angle);
        coords[atom].y = center.y + radius * sin(angle);
        placed[atom] = true;
    }

    return true;
}

/* ============== Chain Placement ============== */

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

        /* Place unplaced neighbors */
        for (int i = 0; i < atom->num_neighbors; i++) {
            int nb = atom->neighbors[i];
            if (in_queue[nb]) continue;

            double angle;

            if (n_angles == 0) {
                angle = 0;
                chain_dir[nb] = 1;
            }
            else if (n_angles == 1) {
                /* Zigzag */
                int dir = -chain_dir[curr];
                if (dir == 0) dir = 1;
                chain_dir[nb] = dir;
                angle = angles[0] + M_PI + dir * M_PI / 3.0;
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
            }

            coords[nb].x = curr_pos.x + bond_length * cos(angle);
            coords[nb].y = curr_pos.y + bond_length * sin(angle);
            placed[nb] = true;

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

        /* Find largest ring in this system to start */
        int first_ring = -1;
        int max_size = 0;
        for (int i = 0; i < nr; i++) {
            if (uf_find(parent, i) == largest_root) {
                if (mol->rings[i].size > max_size) {
                    max_size = mol->rings[i].size;
                    first_ring = i;
                }
            }
        }

        if (first_ring >= 0) {
            place_ring_polygon(&mol->rings[first_ring], coords, placed,
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

        /* Place remaining ring systems (disconnected from main) */
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
                } else if (placed_count == 1) {
                    /* Substituent ring - one atom placed, place ring around it */
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

                    /* Place ring with anchor at edge, extending in away direction */
                    int n = ring->size;
                    double radius = opts.bond_length / (2.0 * sin(M_PI / n));
                    point2d_t center = point2d_add(anchor_pos, point2d_scale(away_dir, radius));

                    double anchor_angle = atan2(anchor_pos.y - center.y, anchor_pos.x - center.x);
                    double angle_step = 2.0 * M_PI / n;

                    for (int i = 0; i < n; i++) {
                        int atom = ring->atoms[i];
                        if (placed[atom]) continue;

                        int steps = (i - placed_ring_idx + n) % n;
                        double angle = anchor_angle - steps * angle_step;

                        coords[atom].x = center.x + radius * cos(angle);
                        coords[atom].y = center.y + radius * sin(angle);
                        placed[atom] = true;
                    }
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
                                away_dir = point2d_scale(sum, -1.0 / len);
                            }
                        }

                        /* Position anchor atom at bond length from neighbor */
                        point2d_t anchor_pos = point2d_add(neighbor_pos,
                                                          point2d_scale(away_dir, opts.bond_length));
                        coords[anchor_atom] = anchor_pos;
                        placed[anchor_atom] = true;

                        /* Place rest of ring around anchor */
                        int n = ring->size;
                        double radius = opts.bond_length / (2.0 * sin(M_PI / n));
                        point2d_t center = point2d_add(anchor_pos, point2d_scale(away_dir, radius));

                        double anchor_angle = atan2(anchor_pos.y - center.y, anchor_pos.x - center.x);
                        double angle_step = 2.0 * M_PI / n;

                        for (int i = 0; i < n; i++) {
                            int atom = ring->atoms[i];
                            if (placed[atom]) continue;

                            int steps = (i - anchor_ring_idx + n) % n;
                            double angle = anchor_angle - steps * angle_step;

                            coords[atom].x = center.x + radius * cos(angle);
                            coords[atom].y = center.y + radius * sin(angle);
                            placed[atom] = true;
                        }
                    } else {
                        /* Completely disconnected - place at offset */
                        double offset_x = 0;
                        for (int i = 0; i < mol->num_atoms; i++) {
                            if (placed[i] && coords[i].x > offset_x) {
                                offset_x = coords[i].x;
                            }
                        }
                        place_ring_polygon(&mol->rings[r], coords, placed,
                                          (point2d_t){offset_x + 3 * opts.bond_length, 0},
                                          opts.bond_length, M_PI / 2.0);
                    }
                }
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
