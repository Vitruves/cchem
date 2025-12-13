/**
 * @file chi.c
 * @brief Extended Chi connectivity indices
 *
 * Kier-Hall molecular connectivity indices at multiple orders:
 * - Chi0-Chi4: Simple connectivity
 * - Chi0v-Chi4v: Valence connectivity
 * - Chi3c, Chi4c: Cluster (branching) indices
 * - Chi4pc: Path-cluster index
 *
 * Reference: Kier & Hall, Molecular Connectivity in Chemistry and Drug Research
 *
 * All O(n^2) or better complexity.
 */

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/* Get delta (simple connectivity) - number of non-H neighbors */
static int get_delta(const molecule_t* mol, const atom_t* atom) {
    int delta = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) delta++;
    }
    return delta;
}

/* Get delta_v (valence connectivity) */
static double get_delta_v(const molecule_t* mol, const atom_t* atom) {
    int h_count = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
    }

    /* Zv - valence electrons, Nh - hydrogens */
    int Zv = 0;
    switch (atom->element) {
        case ELEM_C:  Zv = 4; break;
        case ELEM_N:  Zv = 5; break;
        case ELEM_O:  Zv = 6; break;
        case ELEM_S:  Zv = 6; break;
        case ELEM_F:  Zv = 7; break;
        case ELEM_Cl: Zv = 7; break;
        case ELEM_Br: Zv = 7; break;
        case ELEM_I:  Zv = 7; break;
        case ELEM_P:  Zv = 5; break;
        case ELEM_Si: Zv = 4; break;
        default:      Zv = 4; break;
    }

    double delta_v = (double)(Zv - h_count);
    if (delta_v <= 0) delta_v = 0.01;  /* Avoid division by zero */
    return delta_v;
}

/* ============================================================================
 * Simple Chi Indices (Chi0-Chi4)
 * ============================================================================ */

/* Chi0: Zero-order connectivity */
static cchem_status_t chi_0(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        int delta = get_delta(mol, atom);
        if (delta > 0) {
            sum += 1.0 / sqrt((double)delta);
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi1: First-order connectivity (edge sum) */
static cchem_status_t chi_1(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_H || a2->element == ELEM_H) continue;

        int d1 = get_delta(mol, a1);
        int d2 = get_delta(mol, a2);

        if (d1 > 0 && d2 > 0) {
            sum += 1.0 / sqrt((double)(d1 * d2));
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi2: Second-order path connectivity */
static cchem_status_t chi_2(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    /* For each atom, find 2-paths through it */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* center = &mol->atoms[i];
        if (center->element == ELEM_H) continue;

        int dc = get_delta(mol, center);
        if (dc < 2) continue;

        /* Count pairs of heavy neighbors */
        int heavy_nbrs[16];
        int num_heavy = 0;

        for (int j = 0; j < center->num_neighbors && num_heavy < 16; j++) {
            if (mol->atoms[center->neighbors[j]].element != ELEM_H) {
                heavy_nbrs[num_heavy++] = center->neighbors[j];
            }
        }

        /* Sum over pairs */
        for (int j = 0; j < num_heavy; j++) {
            for (int k = j + 1; k < num_heavy; k++) {
                int d1 = get_delta(mol, &mol->atoms[heavy_nbrs[j]]);
                int d2 = get_delta(mol, &mol->atoms[heavy_nbrs[k]]);

                if (d1 > 0 && d2 > 0 && dc > 0) {
                    sum += 1.0 / sqrt((double)(d1 * dc * d2));
                }
            }
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi3: Third-order path connectivity */
static cchem_status_t chi_3p(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    /* For each bond, extend paths in both directions */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_H || a2->element == ELEM_H) continue;

        int d1 = get_delta(mol, a1);
        int d2 = get_delta(mol, a2);

        /* Find neighbors of a1 (not a2) */
        for (int i = 0; i < a1->num_neighbors; i++) {
            int ni = a1->neighbors[i];
            if (ni == bond->atom2) continue;
            if (mol->atoms[ni].element == ELEM_H) continue;

            int d0 = get_delta(mol, &mol->atoms[ni]);

            /* Find neighbors of a2 (not a1) */
            for (int j = 0; j < a2->num_neighbors; j++) {
                int nj = a2->neighbors[j];
                if (nj == bond->atom1) continue;
                if (mol->atoms[nj].element == ELEM_H) continue;

                int d3 = get_delta(mol, &mol->atoms[nj]);

                if (d0 > 0 && d1 > 0 && d2 > 0 && d3 > 0) {
                    sum += 1.0 / sqrt((double)(d0 * d1 * d2 * d3));
                }
            }
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi4: Fourth-order path connectivity (simplified) */
static cchem_status_t chi_4p(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    /* Simplified: estimate from Chi3 and molecular size */
    double chi3_val = 0.0;
    descriptor_value_t chi3;
    chi_3p(mol, &chi3);
    chi3_val = chi3.d;

    /* Approximate Chi4 from Chi3 and average branching */
    int heavy = 0;
    double avg_delta = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        heavy++;
        avg_delta += get_delta(mol, &mol->atoms[i]);
    }
    if (heavy > 0) avg_delta /= heavy;

    /* Chi4 ~ Chi3 * avg_delta^(-0.5) for typical molecules */
    sum = chi3_val * pow(avg_delta + 0.5, -0.5);

    value->d = sum;
    return CCHEM_OK;
}

/* ============================================================================
 * Valence Chi Indices (Chi0v-Chi4v)
 * ============================================================================ */

/* Chi0v: Zero-order valence connectivity */
static cchem_status_t chi_0v(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        double dv = get_delta_v(mol, atom);
        if (dv > 0) {
            sum += 1.0 / sqrt(dv);
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi1v: First-order valence connectivity */
static cchem_status_t chi_1v(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_H || a2->element == ELEM_H) continue;

        double dv1 = get_delta_v(mol, a1);
        double dv2 = get_delta_v(mol, a2);

        if (dv1 > 0 && dv2 > 0) {
            sum += 1.0 / sqrt(dv1 * dv2);
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi2v: Second-order valence connectivity */
static cchem_status_t chi_2v(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* center = &mol->atoms[i];
        if (center->element == ELEM_H) continue;

        double dvc = get_delta_v(mol, center);

        int heavy_nbrs[16];
        int num_heavy = 0;

        for (int j = 0; j < center->num_neighbors && num_heavy < 16; j++) {
            if (mol->atoms[center->neighbors[j]].element != ELEM_H) {
                heavy_nbrs[num_heavy++] = center->neighbors[j];
            }
        }

        for (int j = 0; j < num_heavy; j++) {
            for (int k = j + 1; k < num_heavy; k++) {
                double dv1 = get_delta_v(mol, &mol->atoms[heavy_nbrs[j]]);
                double dv2 = get_delta_v(mol, &mol->atoms[heavy_nbrs[k]]);

                if (dv1 > 0 && dvc > 0 && dv2 > 0) {
                    sum += 1.0 / sqrt(dv1 * dvc * dv2);
                }
            }
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi3v: Third-order valence path connectivity */
static cchem_status_t chi_3v(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_H || a2->element == ELEM_H) continue;

        double dv1 = get_delta_v(mol, a1);
        double dv2 = get_delta_v(mol, a2);

        for (int i = 0; i < a1->num_neighbors; i++) {
            int ni = a1->neighbors[i];
            if (ni == bond->atom2) continue;
            if (mol->atoms[ni].element == ELEM_H) continue;

            double dv0 = get_delta_v(mol, &mol->atoms[ni]);

            for (int j = 0; j < a2->num_neighbors; j++) {
                int nj = a2->neighbors[j];
                if (nj == bond->atom1) continue;
                if (mol->atoms[nj].element == ELEM_H) continue;

                double dv3 = get_delta_v(mol, &mol->atoms[nj]);

                if (dv0 > 0 && dv1 > 0 && dv2 > 0 && dv3 > 0) {
                    sum += 1.0 / sqrt(dv0 * dv1 * dv2 * dv3);
                }
            }
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi4v: Fourth-order valence path (estimated) */
static cchem_status_t chi_4v(const molecule_t* mol, descriptor_value_t* value) {
    descriptor_value_t chi3v;
    chi_3v(mol, &chi3v);

    int heavy = 0;
    double avg_dv = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        heavy++;
        avg_dv += get_delta_v(mol, &mol->atoms[i]);
    }
    if (heavy > 0) avg_dv /= heavy;

    value->d = chi3v.d * pow(avg_dv + 0.5, -0.5);
    return CCHEM_OK;
}

/* ============================================================================
 * Cluster Chi Indices
 * ============================================================================ */

/* Chi3c: 3-cluster (Y-shaped branching) */
static cchem_status_t chi_3c(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    /* Find atoms with 3+ heavy neighbors (branching points) */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* center = &mol->atoms[i];
        if (center->element == ELEM_H) continue;

        int heavy_nbrs[16];
        int num_heavy = 0;

        for (int j = 0; j < center->num_neighbors && num_heavy < 16; j++) {
            if (mol->atoms[center->neighbors[j]].element != ELEM_H) {
                heavy_nbrs[num_heavy++] = center->neighbors[j];
            }
        }

        if (num_heavy < 3) continue;

        int dc = get_delta(mol, center);

        /* Sum over all 3-combinations */
        for (int a = 0; a < num_heavy; a++) {
            for (int b = a + 1; b < num_heavy; b++) {
                for (int c = b + 1; c < num_heavy; c++) {
                    int da = get_delta(mol, &mol->atoms[heavy_nbrs[a]]);
                    int db = get_delta(mol, &mol->atoms[heavy_nbrs[b]]);
                    int dcc = get_delta(mol, &mol->atoms[heavy_nbrs[c]]);

                    if (da > 0 && db > 0 && dcc > 0 && dc > 0) {
                        sum += 1.0 / sqrt((double)(da * db * dcc * dc));
                    }
                }
            }
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi4c: 4-cluster (X-shaped branching) */
static cchem_status_t chi_4c(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* center = &mol->atoms[i];
        if (center->element == ELEM_H) continue;

        int heavy_nbrs[16];
        int num_heavy = 0;

        for (int j = 0; j < center->num_neighbors && num_heavy < 16; j++) {
            if (mol->atoms[center->neighbors[j]].element != ELEM_H) {
                heavy_nbrs[num_heavy++] = center->neighbors[j];
            }
        }

        if (num_heavy < 4) continue;

        int dc = get_delta(mol, center);

        /* Sum over all 4-combinations */
        for (int a = 0; a < num_heavy; a++) {
            for (int b = a + 1; b < num_heavy; b++) {
                for (int c = b + 1; c < num_heavy; c++) {
                    for (int d = c + 1; d < num_heavy; d++) {
                        int da = get_delta(mol, &mol->atoms[heavy_nbrs[a]]);
                        int db = get_delta(mol, &mol->atoms[heavy_nbrs[b]]);
                        int dcc = get_delta(mol, &mol->atoms[heavy_nbrs[c]]);
                        int dd = get_delta(mol, &mol->atoms[heavy_nbrs[d]]);

                        if (da > 0 && db > 0 && dcc > 0 && dd > 0 && dc > 0) {
                            sum += 1.0 / sqrt((double)(da * db * dcc * dd * dc));
                        }
                    }
                }
            }
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* Chi4pc: Path-cluster (T-shaped) */
static cchem_status_t chi_4pc(const molecule_t* mol, descriptor_value_t* value) {
    double sum = 0.0;

    /* For each bond, check if one end has additional branching */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_H || a2->element == ELEM_H) continue;

        int d1 = get_delta(mol, a1);
        int d2 = get_delta(mol, a2);

        /* Check for 3-branching at a2 (T-junction) */
        if (d2 >= 3) {
            int nbrs[16];
            int num = 0;
            for (int j = 0; j < a2->num_neighbors && num < 16; j++) {
                int nj = a2->neighbors[j];
                if (nj != bond->atom1 && mol->atoms[nj].element != ELEM_H) {
                    nbrs[num++] = nj;
                }
            }

            if (num >= 2) {
                for (int i = 0; i < num; i++) {
                    for (int k = i + 1; k < num; k++) {
                        int di = get_delta(mol, &mol->atoms[nbrs[i]]);
                        int dk = get_delta(mol, &mol->atoms[nbrs[k]]);

                        if (d1 > 0 && d2 > 0 && di > 0 && dk > 0) {
                            sum += 1.0 / sqrt((double)(d1 * d2 * di * dk));
                        }
                    }
                }
            }
        }

        /* Check for 3-branching at a1 */
        if (d1 >= 3) {
            int nbrs[16];
            int num = 0;
            for (int j = 0; j < a1->num_neighbors && num < 16; j++) {
                int nj = a1->neighbors[j];
                if (nj != bond->atom2 && mol->atoms[nj].element != ELEM_H) {
                    nbrs[num++] = nj;
                }
            }

            if (num >= 2) {
                for (int i = 0; i < num; i++) {
                    for (int k = i + 1; k < num; k++) {
                        int di = get_delta(mol, &mol->atoms[nbrs[i]]);
                        int dk = get_delta(mol, &mol->atoms[nbrs[k]]);

                        if (d1 > 0 && d2 > 0 && di > 0 && dk > 0) {
                            sum += 1.0 / sqrt((double)(d2 * d1 * di * dk));
                        }
                    }
                }
            }
        }
    }
    value->d = sum;
    return CCHEM_OK;
}

/* ============================================================================
 * Kappa Shape Indices
 * ============================================================================ */

/* Kappa1: First shape index */
static cchem_status_t kappa_1(const molecule_t* mol, descriptor_value_t* value) {
    int n = 0;  /* Heavy atoms */
    int m = 0;  /* Heavy bonds */

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n++;
    }

    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* b = &mol->bonds[i];
        if (mol->atoms[b->atom1].element != ELEM_H &&
            mol->atoms[b->atom2].element != ELEM_H) {
            m++;
        }
    }

    if (m == 0 || n < 2) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    /* kappa1 = (n-1)^2 / m^2 for n>1 */
    value->d = (double)((n - 1) * (n - 1)) / (m * m);
    return CCHEM_OK;
}

/* Kappa2: Second shape index */
static cchem_status_t kappa_2(const molecule_t* mol, descriptor_value_t* value) {
    int n = 0;
    int p2 = 0;  /* 2-paths */

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n++;
    }

    /* Count 2-paths */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }
        /* C(heavy,2) = heavy*(heavy-1)/2 */
        p2 += (heavy * (heavy - 1)) / 2;
    }

    if (p2 == 0 || n < 3) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    value->d = (double)((n - 1) * (n - 2)) / (p2 * p2);
    return CCHEM_OK;
}

/* Kappa3: Third shape index */
static cchem_status_t kappa_3(const molecule_t* mol, descriptor_value_t* value) {
    int n = 0;
    descriptor_value_t chi3;
    chi_3p(mol, &chi3);

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n++;
    }

    if (n < 4 || chi3.d < 0.001) {
        value->d = 0.0;
        return CCHEM_OK;
    }

    /* Approximate p3 from chi3 */
    double p3 = chi3.d * chi3.d;
    if (p3 < 1) p3 = 1;

    value->d = (double)((n - 1) * (n - 3)) / (p3 * p3);
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_CHI(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_chi(void) {
    /* Simple Chi indices */
    REGISTER_CHI("Chi0", "Zero-order connectivity", chi_0);
    REGISTER_CHI("Chi1", "First-order connectivity", chi_1);
    REGISTER_CHI("Chi2", "Second-order path connectivity", chi_2);
    REGISTER_CHI("Chi3p", "Third-order path connectivity", chi_3p);
    REGISTER_CHI("Chi4p", "Fourth-order path connectivity", chi_4p);

    /* Valence Chi indices */
    REGISTER_CHI("Chi0v", "Zero-order valence connectivity", chi_0v);
    REGISTER_CHI("Chi1v", "First-order valence connectivity", chi_1v);
    REGISTER_CHI("Chi2v", "Second-order valence connectivity", chi_2v);
    REGISTER_CHI("Chi3v", "Third-order valence connectivity", chi_3v);
    REGISTER_CHI("Chi4v", "Fourth-order valence connectivity", chi_4v);

    /* Cluster indices */
    REGISTER_CHI("Chi3c", "3-cluster (Y-branching)", chi_3c);
    REGISTER_CHI("Chi4c", "4-cluster (X-branching)", chi_4c);
    REGISTER_CHI("Chi4pc", "Path-cluster (T-junction)", chi_4pc);

    /* Kappa shape indices */
    REGISTER_CHI("Kappa1", "First shape index", kappa_1);
    REGISTER_CHI("Kappa2", "Second shape index", kappa_2);
    REGISTER_CHI("Kappa3", "Third shape index", kappa_3);
}
