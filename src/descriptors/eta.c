/**
 * @file eta.c
 * @brief Extended Topological Atom (ETA) Descriptors
 *
 * Based on Roy & Chatterjee extended topological atom indices:
 * - ETA_Shape: Molecular shape/branching indices
 * - ETA_Core: Core count descriptors
 * - ETA_Epsilon: Electronic/functional indices
 * - ETA_Composite: Combined indices
 *
 * Total: 24 descriptors
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

#define NUM_ETA_DESCRIPTORS 24
#define MAX_ETA_ATOMS 512

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static int get_heavy_degree(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int degree = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            degree++;
        }
    }
    return degree;
}

static int count_heavy_atoms(const molecule_t* mol) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) count++;
    }
    return count;
}

/* Get valence electrons for element */
static int get_valence_electrons(element_t e) {
    switch (e) {
        case ELEM_C:  return 4;
        case ELEM_N:  return 5;
        case ELEM_O:  return 6;
        case ELEM_S:  return 6;
        case ELEM_P:  return 5;
        case ELEM_F:  return 7;
        case ELEM_Cl: return 7;
        case ELEM_Br: return 7;
        case ELEM_I:  return 7;
        case ELEM_Si: return 4;
        case ELEM_B:  return 3;
        default:      return 4;
    }
}

/* Get principal quantum number for element */
static int get_principal_quantum(element_t e) {
    switch (e) {
        case ELEM_H:  return 1;
        case ELEM_C: case ELEM_N: case ELEM_O: case ELEM_F: return 2;
        case ELEM_Si: case ELEM_P: case ELEM_S: case ELEM_Cl: return 3;
        case ELEM_Br: return 4;
        case ELEM_I:  return 5;
        default:      return 2;
    }
}

/* Compute sigma electrons for atom (bonded electrons) */
static int get_sigma_electrons(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int sigma = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        sigma++;  /* Each bond contributes 1 sigma electron */
    }
    sigma += atom->implicit_h_count;  /* Implicit H bonds */
    return sigma;
}

/* Compute pi electrons for atom */
static int get_pi_electrons(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int pi = 0;

    for (int i = 0; i < atom->num_neighbors; i++) {
        int neighbor_idx = atom->neighbors[i];
        /* Find bond to this neighbor */
        for (int b = 0; b < mol->num_bonds; b++) {
            const bond_t* bond = &mol->bonds[b];
            if ((bond->atom1 == atom_idx && bond->atom2 == neighbor_idx) ||
                (bond->atom2 == atom_idx && bond->atom1 == neighbor_idx)) {
                if (bond->aromatic) {
                    pi += 1;  /* Aromatic contributes ~0.5 each side, we count 1 */
                } else if (bond->type == BOND_DOUBLE) {
                    pi += 1;
                } else if (bond->type == BOND_TRIPLE) {
                    pi += 2;
                }
                break;
            }
        }
    }
    return pi;
}

/* Compute lone pair electrons */
static int get_lone_pairs(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int valence = get_valence_electrons(atom->element);
    int sigma = get_sigma_electrons(mol, atom_idx);
    int pi = get_pi_electrons(mol, atom_idx);

    /* lone pairs = (valence - bonding electrons) / 2 */
    int bonding = sigma + pi;
    int non_bonding = valence - bonding;
    return (non_bonding > 0) ? non_bonding / 2 : 0;
}

/* ============================================================================
 * ETA Alpha Values (vertex contribution)
 * alpha = (Zv - h) / Zv where Zv = valence electrons, h = hydrogens
 * ============================================================================ */

static double compute_eta_alpha(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    if (atom->element == ELEM_H) return 0.0;

    int zv = get_valence_electrons(atom->element);
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h++;
    }

    return (zv > 0) ? (double)(zv - h) / zv : 0.0;
}

/* Ref alpha (reference state contribution) */
static double compute_eta_alpha_ref(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    if (atom->element == ELEM_H) return 0.0;

    int pq = get_principal_quantum(atom->element);
    int zv = get_valence_electrons(atom->element);

    /* Reference alpha incorporates quantum number */
    return (zv > 0) ? (double)pq / zv : 0.0;
}

/* ============================================================================
 * ETA Shape Indices
 * ============================================================================ */

/* ETA_Shape_Y: Branching index */
static double compute_eta_shape_y(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy <= 1) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int d = get_heavy_degree(mol, i);
        if (d >= 3) {
            sum += (d - 2);  /* Branching contribution */
        }
    }

    return sum / (n_heavy - 1);
}

/* ETA_Shape_X: Extension index */
static double compute_eta_shape_x(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy <= 1) return 0.0;

    int n_terminal = 0;  /* Degree 1 atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        if (get_heavy_degree(mol, i) == 1) n_terminal++;
    }

    return (double)n_terminal / n_heavy;
}

/* ETA_Shape_P: Path index (based on longest path) */
static double compute_eta_shape_p(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy <= 1) return 0.0;

    /* Find longest path using BFS from each terminal */
    int max_path = 0;

    for (int start = 0; start < mol->num_atoms; start++) {
        if (mol->atoms[start].element == ELEM_H) continue;
        if (get_heavy_degree(mol, start) != 1) continue;  /* Start from terminal */

        int dist[MAX_ETA_ATOMS];
        int queue[MAX_ETA_ATOMS];
        memset(dist, -1, sizeof(int) * mol->num_atoms);

        int head = 0, tail = 0;
        queue[tail++] = start;
        dist[start] = 0;

        while (head < tail) {
            int curr = queue[head++];
            if (dist[curr] > max_path) max_path = dist[curr];

            const atom_t* atom = &mol->atoms[curr];
            for (int i = 0; i < atom->num_neighbors; i++) {
                int nb = atom->neighbors[i];
                if (mol->atoms[nb].element == ELEM_H) continue;
                if (dist[nb] < 0) {
                    dist[nb] = dist[curr] + 1;
                    queue[tail++] = nb;
                }
            }
        }
    }

    return (n_heavy > 1) ? (double)max_path / (n_heavy - 1) : 0.0;
}

/* ============================================================================
 * ETA Core Count Indices
 * ============================================================================ */

/* ETA_dAlpha_A: Sum of alpha differences from reference */
static double compute_eta_dalpha_a(const molecule_t* mol) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        double alpha = compute_eta_alpha(mol, i);
        double alpha_ref = compute_eta_alpha_ref(mol, i);
        sum += fabs(alpha - alpha_ref);
    }
    return sum;
}

/* ETA_dAlpha_B: Bond-weighted alpha sum */
static double compute_eta_dalpha_b(const molecule_t* mol) {
    double sum = 0.0;
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        double alpha1 = compute_eta_alpha(mol, bond->atom1);
        double alpha2 = compute_eta_alpha(mol, bond->atom2);

        double bond_order = 1.0;
        if (bond->aromatic) bond_order = 1.5;
        else if (bond->type == BOND_DOUBLE) bond_order = 2.0;
        else if (bond->type == BOND_TRIPLE) bond_order = 3.0;

        sum += sqrt(alpha1 * alpha2) * bond_order;
    }
    return sum;
}

/* ============================================================================
 * ETA Epsilon Indices (electronic)
 * ============================================================================ */

/* ETA_Epsilon_1: First epsilon index (sigma electron contribution) */
static double compute_eta_epsilon_1(const molecule_t* mol) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int sigma = get_sigma_electrons(mol, i);
        int zv = get_valence_electrons(mol->atoms[i].element);
        if (zv > 0) {
            sum += (double)sigma / zv;
        }
    }
    return sum;
}

/* ETA_Epsilon_2: Second epsilon index (pi electron contribution) */
static double compute_eta_epsilon_2(const molecule_t* mol) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int pi = get_pi_electrons(mol, i);
        int zv = get_valence_electrons(mol->atoms[i].element);
        if (zv > 0) {
            sum += (double)pi / zv;
        }
    }
    return sum;
}

/* ETA_Epsilon_3: Third epsilon index (lone pair contribution) */
static double compute_eta_epsilon_3(const molecule_t* mol) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int lp = get_lone_pairs(mol, i);
        sum += (double)lp;
    }
    return sum;
}

/* ETA_Epsilon_4: Fourth epsilon index (combined electronic) */
static double compute_eta_epsilon_4(const molecule_t* mol) {
    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int sigma = get_sigma_electrons(mol, i);
        int pi = get_pi_electrons(mol, i);
        int lp = get_lone_pairs(mol, i);
        int zv = get_valence_electrons(mol->atoms[i].element);
        if (zv > 0) {
            sum += (double)(sigma + 2*pi + 3*lp) / zv;
        }
    }
    return sum;
}

/* Pauling electronegativity lookup */
static double get_pauling_en(element_t e) {
    switch (e) {
        case ELEM_H: return 2.20;
        case ELEM_C: return 2.55;
        case ELEM_N: return 3.04;
        case ELEM_O: return 3.44;
        case ELEM_F: return 3.98;
        case ELEM_S: return 2.58;
        case ELEM_P: return 2.19;
        case ELEM_Cl: return 3.16;
        case ELEM_Br: return 2.96;
        case ELEM_I: return 2.66;
        default: return 2.55;
    }
}

/* ETA_Epsilon_5: Fifth epsilon index (electronegativity-weighted) */
static double compute_eta_epsilon_5(const molecule_t* mol) {
    double sum = 0.0;
    double en_c = 2.55;  /* Carbon reference */

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        double en = get_pauling_en(mol->atoms[i].element);
        double alpha = compute_eta_alpha(mol, i);
        sum += alpha * (en / en_c);
    }
    return sum;
}

/* ============================================================================
 * ETA Composite Indices
 * ============================================================================ */

/* ETA_BetaS: Sum of local beta indices (sigma contribution) */
static double compute_eta_beta_s(const molecule_t* mol) {
    double sum = 0.0;
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        double alpha1 = compute_eta_alpha(mol, bond->atom1);
        double alpha2 = compute_eta_alpha(mol, bond->atom2);

        if (alpha1 > 0 && alpha2 > 0) {
            sum += 1.0 / sqrt(alpha1 * alpha2);
        }
    }
    return sum;
}

/* ETA_BetaNS: Normalized beta sigma */
static double compute_eta_beta_ns(const molecule_t* mol) {
    int n_heavy_bonds = 0;
    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->atoms[mol->bonds[b].atom1].element != ELEM_H &&
            mol->atoms[mol->bonds[b].atom2].element != ELEM_H) {
            n_heavy_bonds++;
        }
    }
    if (n_heavy_bonds == 0) return 0.0;

    return compute_eta_beta_s(mol) / n_heavy_bonds;
}

/* ETA_Psi_1: First psi index (size-dependent complexity) */
static double compute_eta_psi_1(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy == 0) return 0.0;

    double alpha_sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        alpha_sum += compute_eta_alpha(mol, i);
    }

    return alpha_sum / sqrt((double)n_heavy);
}

/* ETA_dBeta: Bond-type differentiation */
static double compute_eta_dbeta(const molecule_t* mol) {
    double single_sum = 0.0, multiple_sum = 0.0;
    int n_single = 0, n_multiple = 0;

    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        double alpha1 = compute_eta_alpha(mol, bond->atom1);
        double alpha2 = compute_eta_alpha(mol, bond->atom2);
        double contrib = (alpha1 > 0 && alpha2 > 0) ? 1.0 / sqrt(alpha1 * alpha2) : 0.0;

        if (bond->aromatic || bond->type == BOND_DOUBLE || bond->type == BOND_TRIPLE) {
            multiple_sum += contrib;
            n_multiple++;
        } else {
            single_sum += contrib;
            n_single++;
        }
    }

    double avg_single = (n_single > 0) ? single_sum / n_single : 0.0;
    double avg_multiple = (n_multiple > 0) ? multiple_sum / n_multiple : 0.0;

    return avg_multiple - avg_single;
}

/* ETA_dEpsilon: Electronic differentiation */
static double compute_eta_deps(const molecule_t* mol) {
    double eps1 = compute_eta_epsilon_1(mol);
    double eps2 = compute_eta_epsilon_2(mol);
    double eps3 = compute_eta_epsilon_3(mol);

    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy == 0) return 0.0;

    return (eps2 + eps3) / (eps1 + 0.001);  /* Pi + LP contribution relative to sigma */
}

/* ETA_Eta: Main eta index (composite) */
static double compute_eta_eta(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy == 0) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int d = get_heavy_degree(mol, i);
        double alpha = compute_eta_alpha(mol, i);
        sum += d * alpha;
    }

    return sum / n_heavy;
}

/* ETA_EtaR: Ring-weighted eta */
static double compute_eta_eta_ring(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy == 0) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        const atom_t* atom = &mol->atoms[i];
        int d = get_heavy_degree(mol, i);
        double alpha = compute_eta_alpha(mol, i);

        /* Weight by ring membership */
        double ring_weight = (atom->ring_count > 0) ? 1.0 + 0.5 * atom->ring_count : 1.0;
        sum += d * alpha * ring_weight;
    }

    return sum / n_heavy;
}

/* ETA_EtaP: Path-based eta */
static double compute_eta_eta_path(const molecule_t* mol) {
    double shape_p = compute_eta_shape_p(mol);
    double eta = compute_eta_eta(mol);
    return eta * (1.0 + shape_p);
}

/* ETA_EtaF: Functional eta (heteroatom weighted) */
static double compute_eta_eta_f(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy == 0) return 0.0;

    double sum = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (e == ELEM_H) continue;

        int d = get_heavy_degree(mol, i);
        double alpha = compute_eta_alpha(mol, i);

        if (e != ELEM_C) {
            sum += d * alpha * 2.0;  /* Weight heteroatoms */
        } else {
            sum += d * alpha;
        }
    }

    return sum / n_heavy;
}

/* ============================================================================
 * Batch Computation
 * ============================================================================ */

int descriptors_compute_eta_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_ETA_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    if (mol->num_atoms == 0 || mol->num_atoms > MAX_ETA_ATOMS) {
        return NUM_ETA_DESCRIPTORS;
    }

    int idx = 0;

    /* Shape indices */
    values[idx++].d = compute_eta_shape_y(mol);     /* 0: ETA_Shape_Y */
    values[idx++].d = compute_eta_shape_x(mol);     /* 1: ETA_Shape_X */
    values[idx++].d = compute_eta_shape_p(mol);     /* 2: ETA_Shape_P */

    /* Core indices */
    values[idx++].d = compute_eta_dalpha_a(mol);    /* 3: ETA_dAlpha_A */
    values[idx++].d = compute_eta_dalpha_b(mol);    /* 4: ETA_dAlpha_B */

    /* Epsilon indices */
    values[idx++].d = compute_eta_epsilon_1(mol);   /* 5: ETA_Epsilon_1 */
    values[idx++].d = compute_eta_epsilon_2(mol);   /* 6: ETA_Epsilon_2 */
    values[idx++].d = compute_eta_epsilon_3(mol);   /* 7: ETA_Epsilon_3 */
    values[idx++].d = compute_eta_epsilon_4(mol);   /* 8: ETA_Epsilon_4 */
    values[idx++].d = compute_eta_epsilon_5(mol);   /* 9: ETA_Epsilon_5 */

    /* Composite indices */
    values[idx++].d = compute_eta_beta_s(mol);      /* 10: ETA_BetaS */
    values[idx++].d = compute_eta_beta_ns(mol);     /* 11: ETA_BetaNS */
    values[idx++].d = compute_eta_psi_1(mol);       /* 12: ETA_Psi_1 */
    values[idx++].d = compute_eta_dbeta(mol);       /* 13: ETA_dBeta */
    values[idx++].d = compute_eta_deps(mol);        /* 14: ETA_dEpsilon */
    values[idx++].d = compute_eta_eta(mol);         /* 15: ETA_Eta */
    values[idx++].d = compute_eta_eta_ring(mol);    /* 16: ETA_EtaR */
    values[idx++].d = compute_eta_eta_path(mol);    /* 17: ETA_EtaP */
    values[idx++].d = compute_eta_eta_f(mol);       /* 18: ETA_EtaF */

    /* Additional summary indices */
    int n_heavy = count_heavy_atoms(mol);
    double alpha_sum = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            alpha_sum += compute_eta_alpha(mol, i);
        }
    }

    values[idx++].d = alpha_sum;                              /* 19: ETA_Alpha */
    values[idx++].d = (n_heavy > 0) ? alpha_sum / n_heavy : 0.0;  /* 20: ETA_AlphaMean */
    values[idx++].d = values[3].d + values[4].d;              /* 21: ETA_dAlpha */
    values[idx++].d = values[5].d + values[6].d + values[7].d; /* 22: ETA_EpsilonSum */
    values[idx++].d = values[15].d * n_heavy;                 /* 23: ETA_EtaTotal */

    return NUM_ETA_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* eta_cached_mol = NULL;
static _Thread_local descriptor_value_t eta_cached_values[NUM_ETA_DESCRIPTORS];

static inline void ensure_eta_computed(const molecule_t* mol) {
    if (eta_cached_mol != mol) {
        descriptors_compute_eta_all(mol, eta_cached_values);
        eta_cached_mol = mol;
    }
}

#define DEFINE_ETA_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_eta_computed(mol); \
    value->d = eta_cached_values[idx].d; \
    return CCHEM_OK; \
}

DEFINE_ETA_FUNC(eta_shape_y, 0)
DEFINE_ETA_FUNC(eta_shape_x, 1)
DEFINE_ETA_FUNC(eta_shape_p, 2)
DEFINE_ETA_FUNC(eta_dalpha_a, 3)
DEFINE_ETA_FUNC(eta_dalpha_b, 4)
DEFINE_ETA_FUNC(eta_epsilon_1, 5)
DEFINE_ETA_FUNC(eta_epsilon_2, 6)
DEFINE_ETA_FUNC(eta_epsilon_3, 7)
DEFINE_ETA_FUNC(eta_epsilon_4, 8)
DEFINE_ETA_FUNC(eta_epsilon_5, 9)
DEFINE_ETA_FUNC(eta_beta_s, 10)
DEFINE_ETA_FUNC(eta_beta_ns, 11)
DEFINE_ETA_FUNC(eta_psi_1, 12)
DEFINE_ETA_FUNC(eta_dbeta, 13)
DEFINE_ETA_FUNC(eta_deps, 14)
DEFINE_ETA_FUNC(eta_eta, 15)
DEFINE_ETA_FUNC(eta_eta_ring, 16)
DEFINE_ETA_FUNC(eta_eta_path, 17)
DEFINE_ETA_FUNC(eta_eta_f, 18)
DEFINE_ETA_FUNC(eta_alpha, 19)
DEFINE_ETA_FUNC(eta_alpha_mean, 20)
DEFINE_ETA_FUNC(eta_dalpha, 21)
DEFINE_ETA_FUNC(eta_epsilon_sum, 22)
DEFINE_ETA_FUNC(eta_eta_total, 23)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_ETA(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_eta(void) {
    /* Shape indices */
    REGISTER_ETA("ETA_Shape_Y", "ETA branching index", desc_eta_shape_y);
    REGISTER_ETA("ETA_Shape_X", "ETA extension index", desc_eta_shape_x);
    REGISTER_ETA("ETA_Shape_P", "ETA path index", desc_eta_shape_p);

    /* Core indices */
    REGISTER_ETA("ETA_dAlpha_A", "ETA alpha difference A", desc_eta_dalpha_a);
    REGISTER_ETA("ETA_dAlpha_B", "ETA bond-weighted alpha", desc_eta_dalpha_b);

    /* Epsilon indices */
    REGISTER_ETA("ETA_Epsilon_1", "ETA sigma electron index", desc_eta_epsilon_1);
    REGISTER_ETA("ETA_Epsilon_2", "ETA pi electron index", desc_eta_epsilon_2);
    REGISTER_ETA("ETA_Epsilon_3", "ETA lone pair index", desc_eta_epsilon_3);
    REGISTER_ETA("ETA_Epsilon_4", "ETA combined electronic index", desc_eta_epsilon_4);
    REGISTER_ETA("ETA_Epsilon_5", "ETA EN-weighted index", desc_eta_epsilon_5);

    /* Composite indices */
    REGISTER_ETA("ETA_BetaS", "ETA sigma beta sum", desc_eta_beta_s);
    REGISTER_ETA("ETA_BetaNS", "ETA normalized beta", desc_eta_beta_ns);
    REGISTER_ETA("ETA_Psi_1", "ETA psi complexity index", desc_eta_psi_1);
    REGISTER_ETA("ETA_dBeta", "ETA bond differentiation", desc_eta_dbeta);
    REGISTER_ETA("ETA_dEpsilon", "ETA electronic differentiation", desc_eta_deps);
    REGISTER_ETA("ETA_Eta", "ETA main index", desc_eta_eta);
    REGISTER_ETA("ETA_EtaR", "ETA ring-weighted index", desc_eta_eta_ring);
    REGISTER_ETA("ETA_EtaP", "ETA path-based index", desc_eta_eta_path);
    REGISTER_ETA("ETA_EtaF", "ETA functional index", desc_eta_eta_f);

    /* Summary indices */
    REGISTER_ETA("ETA_Alpha", "ETA total alpha", desc_eta_alpha);
    REGISTER_ETA("ETA_AlphaMean", "ETA mean alpha", desc_eta_alpha_mean);
    REGISTER_ETA("ETA_dAlpha", "ETA combined dAlpha", desc_eta_dalpha);
    REGISTER_ETA("ETA_EpsilonSum", "ETA epsilon sum", desc_eta_epsilon_sum);
    REGISTER_ETA("ETA_EtaTotal", "ETA total index", desc_eta_eta_total);
}
