/**
 * @file ratios.c
 * @brief Ratio-based molecular descriptors
 *
 * Ultra-fast ratio descriptors computed in a single pass over molecule atoms.
 * All ratios are dimensionless and normalized for cross-molecule comparison.
 *
 * Group A: Elemental Stoichiometry
 * Group B: Electronic Hybridization
 * Group C: Electronegativity & Hardness
 * Group D: Bond Dynamics
 * Group E: Charge & Ionization
 * Group F: Topological Electronics
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Physical Constants
 * ============================================================================ */

/* Pauling electronegativities */
static const double ELECTRONEGATIVITY[] = {
    [ELEM_H]  = 2.20,
    [ELEM_C]  = 2.55,
    [ELEM_N]  = 3.04,
    [ELEM_O]  = 3.44,
    [ELEM_F]  = 3.98,
    [ELEM_P]  = 2.19,
    [ELEM_S]  = 2.58,
    [ELEM_Cl] = 3.16,
    [ELEM_Br] = 2.96,
    [ELEM_I]  = 2.66,
    [ELEM_Si] = 1.90,
    [ELEM_B]  = 2.04,
    [ELEM_Se] = 2.55,
    [ELEM_As] = 2.18,
};

/* Atomic masses */
static const double ATOMIC_MASS[] = {
    [ELEM_H]  = 1.008,
    [ELEM_C]  = 12.011,
    [ELEM_N]  = 14.007,
    [ELEM_O]  = 15.999,
    [ELEM_F]  = 18.998,
    [ELEM_P]  = 30.974,
    [ELEM_S]  = 32.065,
    [ELEM_Cl] = 35.453,
    [ELEM_Br] = 79.904,
    [ELEM_I]  = 126.904,
    [ELEM_Si] = 28.086,
    [ELEM_B]  = 10.811,
    [ELEM_Se] = 78.96,
    [ELEM_As] = 74.922,
};

/* Valence electrons */
static const int VALENCE_ELECTRONS[] = {
    [ELEM_H]  = 1,
    [ELEM_C]  = 4,
    [ELEM_N]  = 5,
    [ELEM_O]  = 6,
    [ELEM_F]  = 7,
    [ELEM_P]  = 5,
    [ELEM_S]  = 6,
    [ELEM_Cl] = 7,
    [ELEM_Br] = 7,
    [ELEM_I]  = 7,
    [ELEM_Si] = 4,
    [ELEM_B]  = 3,
    [ELEM_Se] = 6,
    [ELEM_As] = 5,
};

/* First ionization potential in eV */
static const double IONIZATION_POTENTIAL[] = {
    [ELEM_H]  = 13.598,
    [ELEM_C]  = 11.260,
    [ELEM_N]  = 14.534,
    [ELEM_O]  = 13.618,
    [ELEM_F]  = 17.422,
    [ELEM_P]  = 10.486,
    [ELEM_S]  = 10.360,
    [ELEM_Cl] = 12.967,
    [ELEM_Br] = 11.814,
    [ELEM_I]  = 10.451,
    [ELEM_Si] = 8.151,
    [ELEM_B]  = 8.298,
    [ELEM_Se] = 9.752,
    [ELEM_As] = 9.815,
};

/* Van der Waals radii in Angstroms (Bondi) */
static const double VDW_RADII[] = {
    [ELEM_H]  = 1.20,
    [ELEM_C]  = 1.70,
    [ELEM_N]  = 1.55,
    [ELEM_O]  = 1.52,
    [ELEM_F]  = 1.47,
    [ELEM_P]  = 1.80,
    [ELEM_S]  = 1.80,
    [ELEM_Cl] = 1.75,
    [ELEM_Br] = 1.85,
    [ELEM_I]  = 1.98,
    [ELEM_Si] = 2.10,
    [ELEM_B]  = 1.92,
    [ELEM_Se] = 1.90,
    [ELEM_As] = 1.85,
};

/* TPSA contributions per atom type (simplified Ertl model) */
/* For N/O based on hybridization and H-count */
static const double TPSA_N_SP3_H2 = 26.03;   /* -NH2 */
static const double TPSA_N_SP3_H1 = 12.03;   /* -NH- */
static const double TPSA_N_SP3_H0 = 3.24;    /* -N< */
static const double TPSA_N_SP2_H1 = 12.36;   /* =NH */
static const double TPSA_N_SP2_H0 = 3.01;    /* =N- */
static const double TPSA_N_AROM_H1 = 15.79;  /* aromatic NH */
static const double TPSA_N_AROM_H0 = 12.89;  /* aromatic N */
static const double TPSA_O_SP3_H1 = 20.23;   /* -OH */
static const double TPSA_O_SP3_H0 = 9.23;    /* -O- */
static const double TPSA_O_SP2 = 17.07;      /* =O */
static const double TPSA_S_H1 = 25.30;       /* -SH */
static const double TPSA_S_H0 = 25.30;       /* -S- (simplified) */

/* Lone pairs by element (for Lewis basicity estimate) */
static const int LONE_PAIRS[] = {
    [ELEM_H]  = 0,
    [ELEM_C]  = 0,
    [ELEM_N]  = 1,  /* sp3: 1 LP, sp2: 1 LP */
    [ELEM_O]  = 2,  /* sp3: 2 LP, sp2: 2 LP */
    [ELEM_F]  = 3,
    [ELEM_P]  = 1,
    [ELEM_S]  = 2,
    [ELEM_Cl] = 3,
    [ELEM_Br] = 3,
    [ELEM_I]  = 3,
    [ELEM_Si] = 0,
    [ELEM_B]  = 0,
    [ELEM_Se] = 2,
    [ELEM_As] = 1,
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_electronegativity(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ELECTRONEGATIVITY)/sizeof(ELECTRONEGATIVITY[0])))
        return 2.55;
    double chi = ELECTRONEGATIVITY[elem];
    return (chi > 0.0) ? chi : 2.55;
}

static inline double get_atomic_mass(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ATOMIC_MASS)/sizeof(ATOMIC_MASS[0])))
        return 12.011;
    double m = ATOMIC_MASS[elem];
    return (m > 0.0) ? m : 12.011;
}

static inline int get_valence_electrons(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(VALENCE_ELECTRONS)/sizeof(VALENCE_ELECTRONS[0])))
        return 4;
    int v = VALENCE_ELECTRONS[elem];
    return (v > 0) ? v : 4;
}

static inline double get_ionization_potential(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(IONIZATION_POTENTIAL)/sizeof(IONIZATION_POTENTIAL[0])))
        return 11.26;
    double ip = IONIZATION_POTENTIAL[elem];
    return (ip > 0.0) ? ip : 11.26;
}

static inline double get_vdw_radius(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(VDW_RADII)/sizeof(VDW_RADII[0])))
        return 1.70;
    double r = VDW_RADII[elem];
    return (r > 0.0) ? r : 1.70;
}

static inline int get_lone_pairs(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(LONE_PAIRS)/sizeof(LONE_PAIRS[0])))
        return 0;
    return LONE_PAIRS[elem];
}

static inline double calc_vdw_volume(double radius) {
    /* V = (4/3) * pi * r^3 */
    return (4.0 / 3.0) * 3.14159265359 * radius * radius * radius;
}

static inline bool is_halogen(element_t e) {
    return (e == ELEM_F || e == ELEM_Cl || e == ELEM_Br || e == ELEM_I);
}

static inline bool is_polar_atom(element_t e) {
    return (e == ELEM_N || e == ELEM_O || e == ELEM_S || e == ELEM_P);
}

/* ============================================================================
 * Statistics Collection Structure
 * ============================================================================ */

typedef struct {
    /* Element counts */
    int n_total;
    int n_heavy;
    int n_carbon;
    int n_hydrogen;
    int n_nitrogen;
    int n_oxygen;
    int n_halogen;
    int n_heteroatom;

    /* Individual halogen counts (Group G) */
    int n_fluorine;
    int n_chlorine;
    int n_bromine;
    int n_iodine;
    int n_sulfur;
    int n_phosphorus;

    /* Hybridization counts */
    int n_sp3;
    int n_sp2;
    int n_sp;
    int n_aromatic;

    /* Bond counts */
    int n_bonds;
    int n_single;
    int n_double;
    int n_triple;
    int n_rotatable;
    int n_ring_bonds;
    int n_aromatic_bonds;  /* Group H */

    /* Ring counts */
    int n_ring_atoms;
    int n_aromatic_atoms;
    int n_rings;
    int n_aromatic_rings;
    int n_saturated_rings;

    /* Charge/polarity */
    int n_charged;
    int n_positive;
    int n_negative;
    int n_donors;
    int n_acceptors;

    /* Sums */
    double sum_chi;
    double sum_chi_sq;
    double sum_mass;
    double sum_valence;
    double sum_ip;
    double max_chi;
    double min_chi;

    /* Mass-based (Group I) */
    double sum_halogen_mass;
    double sum_hydrogen_mass;
    double sum_polar_mass;  /* N + O + S */
    double sum_vdw_volume;

    /* Connectivity (Group J) */
    int n_terminal;
    int n_branch;
    int n_degree2;  /* linker atoms */
    int n_degree3plus;  /* core atoms */
    int n_chiral_carbons;

    /* Ring environment (Group K) */
    int n_fused_ring_atoms;  /* atoms in 2+ rings */
    int n_hetero_in_rings;
    int total_ring_atoms;  /* for ring fractions */

    /* Electronic/solvation (Group L) */
    int n_conjugated;
    int n_lone_pairs;
    int n_polar_buried;  /* polar with degree > 2 */
    int n_lipophilic_buried;  /* C/halogen with degree > 2 */
    int n_polar_total;
    int n_lipophilic_total;  /* C + halogens */
    double sum_tpsa;
    double sum_pi_electrons;
} mol_stats_t;

static void collect_mol_stats(const molecule_t* mol, mol_stats_t* stats) {
    memset(stats, 0, sizeof(mol_stats_t));
    stats->min_chi = 999.0;

    /* === SINGLE PASS THROUGH ATOMS === */
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        element_t elem = atom->element;

        double chi = get_electronegativity(elem);
        double mass = get_atomic_mass(elem);
        int valence = get_valence_electrons(elem);
        double ip = get_ionization_potential(elem);
        double vdw_r = get_vdw_radius(elem);

        stats->n_total++;
        stats->sum_chi += chi;
        stats->sum_chi_sq += chi * chi;
        stats->sum_mass += mass;
        stats->sum_valence += valence;
        stats->sum_ip += ip;
        stats->sum_vdw_volume += calc_vdw_volume(vdw_r);

        if (chi > stats->max_chi) stats->max_chi = chi;
        if (chi < stats->min_chi) stats->min_chi = chi;

        if (elem == ELEM_H) {
            stats->n_hydrogen++;
            stats->sum_hydrogen_mass += mass;
            continue;
        }

        stats->n_heavy++;

        /* Element classification with individual counts */
        switch (elem) {
            case ELEM_C:
                stats->n_carbon++;
                stats->n_lipophilic_total++;
                break;
            case ELEM_N:
                stats->n_nitrogen++;
                stats->n_heteroatom++;
                stats->sum_polar_mass += mass;
                stats->n_polar_total++;
                stats->n_lone_pairs += get_lone_pairs(elem);
                break;
            case ELEM_O:
                stats->n_oxygen++;
                stats->n_heteroatom++;
                stats->sum_polar_mass += mass;
                stats->n_polar_total++;
                stats->n_lone_pairs += get_lone_pairs(elem);
                break;
            case ELEM_S:
                stats->n_sulfur++;
                stats->n_heteroatom++;
                stats->sum_polar_mass += mass;
                stats->n_polar_total++;
                stats->n_lone_pairs += get_lone_pairs(elem);
                break;
            case ELEM_P:
                stats->n_phosphorus++;
                stats->n_heteroatom++;
                stats->n_polar_total++;
                stats->n_lone_pairs += get_lone_pairs(elem);
                break;
            case ELEM_F:
                stats->n_fluorine++;
                stats->n_halogen++;
                stats->n_heteroatom++;
                stats->sum_halogen_mass += mass;
                stats->n_lipophilic_total++;
                break;
            case ELEM_Cl:
                stats->n_chlorine++;
                stats->n_halogen++;
                stats->n_heteroatom++;
                stats->sum_halogen_mass += mass;
                stats->n_lipophilic_total++;
                break;
            case ELEM_Br:
                stats->n_bromine++;
                stats->n_halogen++;
                stats->n_heteroatom++;
                stats->sum_halogen_mass += mass;
                stats->n_lipophilic_total++;
                break;
            case ELEM_I:
                stats->n_iodine++;
                stats->n_halogen++;
                stats->n_heteroatom++;
                stats->sum_halogen_mass += mass;
                stats->n_lipophilic_total++;
                break;
            default:
                stats->n_heteroatom++;
                break;
        }

        /* Hybridization (approximate from connectivity) */
        int degree = atom->num_neighbors;
        int h_count = atom->implicit_h_count;
        int total_bonds = degree + h_count;
        int total_h = h_count;

        /* Check for double/triple bonds and count explicit H */
        int double_bonds = 0, triple_bonds = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt == BOND_DOUBLE) double_bonds++;
            else if (bt == BOND_TRIPLE) triple_bonds++;
            if (mol->atoms[atom->neighbors[j]].element == ELEM_H) total_h++;
        }

        /* Hybridization classification */
        bool is_sp2 = false;
        if (atom->aromatic) {
            stats->n_aromatic++;
            stats->n_sp2++;
            is_sp2 = true;
        } else if (triple_bonds > 0 || (elem == ELEM_C && total_bonds == 2)) {
            stats->n_sp++;
        } else if (double_bonds > 0 || (elem == ELEM_C && total_bonds == 3)) {
            stats->n_sp2++;
            is_sp2 = true;
        } else {
            stats->n_sp3++;
        }

        /* Ring membership */
        if (atom->ring_count > 0) {
            stats->n_ring_atoms++;
            stats->total_ring_atoms++;
            if (atom->aromatic) stats->n_aromatic_atoms++;
            /* Fused ring atoms (in 2+ rings) */
            if (atom->ring_count >= 2) stats->n_fused_ring_atoms++;
            /* Heteroatoms in rings */
            if (elem != ELEM_C) stats->n_hetero_in_rings++;
        }

        /* Charge */
        if (atom->charge != 0) {
            stats->n_charged++;
            if (atom->charge > 0) stats->n_positive++;
            else stats->n_negative++;
        }

        /* H-bond donors/acceptors and TPSA */
        if (elem == ELEM_N) {
            stats->n_acceptors++;
            if (total_h > 0) stats->n_donors++;
            /* TPSA for N */
            if (atom->aromatic) {
                stats->sum_tpsa += (total_h > 0) ? TPSA_N_AROM_H1 : TPSA_N_AROM_H0;
            } else if (is_sp2 || double_bonds > 0) {
                stats->sum_tpsa += (total_h > 0) ? TPSA_N_SP2_H1 : TPSA_N_SP2_H0;
            } else {
                if (total_h >= 2) stats->sum_tpsa += TPSA_N_SP3_H2;
                else if (total_h == 1) stats->sum_tpsa += TPSA_N_SP3_H1;
                else stats->sum_tpsa += TPSA_N_SP3_H0;
            }
        } else if (elem == ELEM_O) {
            stats->n_acceptors++;
            if (total_h > 0) stats->n_donors++;
            /* TPSA for O */
            if (double_bonds > 0) {
                stats->sum_tpsa += TPSA_O_SP2;
            } else {
                stats->sum_tpsa += (total_h > 0) ? TPSA_O_SP3_H1 : TPSA_O_SP3_H0;
            }
        } else if (elem == ELEM_S) {
            /* TPSA for S */
            stats->sum_tpsa += (total_h > 0) ? TPSA_S_H1 : TPSA_S_H0;
        }

        /* Connectivity - count heavy neighbors */
        int heavy_neighbors = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy_neighbors++;
        }

        if (heavy_neighbors <= 1) stats->n_terminal++;
        if (heavy_neighbors >= 3) stats->n_branch++;
        if (heavy_neighbors == 2) stats->n_degree2++;
        if (heavy_neighbors >= 3) stats->n_degree3plus++;

        /* Buried atoms (degree > 2) */
        if (heavy_neighbors > 2) {
            if (is_polar_atom(elem)) stats->n_polar_buried++;
            if (elem == ELEM_C || is_halogen(elem)) stats->n_lipophilic_buried++;
        }

        /* Chirality */
        if (elem == ELEM_C && atom->chirality != CHIRALITY_NONE) {
            stats->n_chiral_carbons++;
        }

        /* Conjugated atoms (aromatic or has both single and double bonds) */
        bool has_single = false, has_double = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) has_single = true;
            if (bt == BOND_DOUBLE) has_double = true;
        }
        if (atom->aromatic || (has_single && has_double)) stats->n_conjugated++;

        /* Pi electrons: double bonds contribute 2, aromatic ~1.5, triple 2 (linear) */
        if (atom->aromatic) {
            stats->sum_pi_electrons += 1.0;  /* each aromatic atom contributes ~1 */
        }

        /* Add implicit H contributions */
        stats->n_total += h_count;
        stats->n_hydrogen += h_count;
        stats->sum_hydrogen_mass += h_count * get_atomic_mass(ELEM_H);
        stats->sum_chi += h_count * get_electronegativity(ELEM_H);
        stats->sum_chi_sq += h_count * get_electronegativity(ELEM_H) * get_electronegativity(ELEM_H);
        stats->sum_mass += h_count * get_atomic_mass(ELEM_H);
        stats->sum_valence += h_count * get_valence_electrons(ELEM_H);
        stats->sum_ip += h_count * get_ionization_potential(ELEM_H);
        stats->sum_vdw_volume += h_count * calc_vdw_volume(get_vdw_radius(ELEM_H));
    }

    /* === SINGLE PASS THROUGH BONDS === */
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];

        stats->n_bonds++;

        switch (bond->type) {
            case BOND_SINGLE:
            case BOND_UP:
            case BOND_DOWN:
                stats->n_single++;
                break;
            case BOND_DOUBLE:
                stats->n_double++;
                stats->sum_pi_electrons += 2.0;  /* double bond = 2 pi electrons */
                break;
            case BOND_TRIPLE:
                stats->n_triple++;
                stats->sum_pi_electrons += 2.0;  /* triple bond = 2 pi electrons (one pi pair counted) */
                break;
            case BOND_AROMATIC:
                stats->n_aromatic_bonds++;
                /* pi electrons already counted via aromatic atoms */
                break;
            default:
                stats->n_single++;
                break;
        }

        /* Check aromatic flag */
        if (bond->aromatic) stats->n_aromatic_bonds++;

        if (bond->in_ring) stats->n_ring_bonds++;

        /* Rotatable: single, non-ring, non-terminal */
        if ((bond->type == BOND_SINGLE || bond->type == BOND_UP || bond->type == BOND_DOWN) && !bond->in_ring) {
            const atom_t* a1 = &mol->atoms[bond->atom1];
            const atom_t* a2 = &mol->atoms[bond->atom2];
            if (a1->element != ELEM_H && a2->element != ELEM_H) {
                int h1 = 0, h2 = 0;
                for (int j = 0; j < a1->num_neighbors; j++) {
                    if (mol->atoms[a1->neighbors[j]].element != ELEM_H) h1++;
                }
                for (int j = 0; j < a2->num_neighbors; j++) {
                    if (mol->atoms[a2->neighbors[j]].element != ELEM_H) h2++;
                }
                if (h1 > 1 && h2 > 1) stats->n_rotatable++;
            }
        }
    }

    /* === RING STATISTICS === */
    stats->n_rings = mol->num_rings;
    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        if (ring->aromatic) {
            stats->n_aromatic_rings++;
        } else {
            /* Check if saturated (no double bonds in ring) */
            bool all_single = true;
            for (int j = 0; j < ring->size && all_single; j++) {
                if (j < ring->size) {
                    int bond_idx = ring->bonds[j];
                    if (bond_idx >= 0 && bond_idx < mol->num_bonds) {
                        bond_type_t bt = mol->bonds[bond_idx].type;
                        if (bt == BOND_DOUBLE || bt == BOND_TRIPLE) all_single = false;
                    }
                }
            }
            if (all_single) stats->n_saturated_rings++;
        }
    }
}

/* ============================================================================
 * Batch Computation (thread-safe: stack-allocated stats)
 * ============================================================================ */

#define NUM_RATIOS_DESCRIPTORS 60

/* Batch compute ALL ratio descriptors - collects stats once */
int descriptors_compute_ratios_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    mol_stats_t stats;
    collect_mol_stats(mol, &stats);

    int idx = 0;

    /* Group A: Elemental Stoichiometry (5) */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_carbon / stats.n_heavy : 0.0;           /* CarbonHeavyRatio */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_heteroatom / stats.n_heavy : 0.0;       /* HeteroHeavyRatio */
    values[idx++].d = (stats.n_carbon > 0) ? (double)stats.n_hydrogen / stats.n_carbon : 0.0;       /* HydrogenCarbonRatio */
    values[idx++].d = (stats.n_total > 0) ? (double)stats.n_halogen / stats.n_total : 0.0;          /* HalogenTotalRatio */
    values[idx++].d = (stats.n_heavy > 0) ? (double)(stats.n_nitrogen + stats.n_oxygen) / stats.n_heavy : 0.0; /* PolarAtomRatio */

    /* Group B: Electronic Hybridization (5) */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_sp3 / stats.n_heavy : 0.0;              /* Sp3Fraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_sp2 / stats.n_heavy : 0.0;              /* Sp2Fraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_sp / stats.n_heavy : 0.0;               /* SpFraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_aromatic / stats.n_heavy : 0.0;         /* AromaticFraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_ring_atoms / stats.n_heavy : 0.0;       /* RingAtomFraction */

    /* Group C: Electronegativity & Hardness (5) */
    double mean_chi = (stats.n_total > 0) ? stats.sum_chi / stats.n_total : 0.0;
    double var_chi = (stats.n_total > 0) ? (stats.sum_chi_sq / stats.n_total) - (mean_chi * mean_chi) : 0.0;
    if (var_chi < 0) var_chi = 0.0;

    values[idx++].d = mean_chi;                                                                      /* MeanElectronegativity */
    values[idx++].d = stats.max_chi - stats.min_chi;                                                /* ChiRange */
    values[idx++].d = var_chi;                                                                       /* ChiVariance */
    values[idx++].d = (mean_chi > 0.001) ? (stats.max_chi - stats.min_chi) / mean_chi : 0.0;       /* HardnessIndex */
    values[idx++].d = (stats.sum_mass > 0) ? stats.sum_chi / stats.sum_mass : 0.0;                 /* ChiMassRatio */

    /* Group D: Bond Dynamics (5) */
    values[idx++].d = (stats.n_bonds > 0) ? (double)stats.n_double / stats.n_bonds : 0.0;           /* DoubleBondFraction */
    values[idx++].d = (stats.n_bonds > 0) ? (double)stats.n_triple / stats.n_bonds : 0.0;           /* TripleBondFraction */
    values[idx++].d = (stats.n_bonds > 0) ? (double)stats.n_rotatable / stats.n_bonds : 0.0;        /* RotatableFraction */
    values[idx++].d = (stats.n_bonds > 0) ? (double)stats.n_ring_bonds / stats.n_bonds : 0.0;       /* RingBondFraction */
    int weighted_unsat = 2 * stats.n_double + 3 * stats.n_triple;
    values[idx++].d = (stats.n_bonds > 0) ? (double)weighted_unsat / stats.n_bonds : 0.0;           /* UnsaturationIndex */

    /* Group E: Charge & Ionization (4) */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_charged / stats.n_heavy : 0.0;          /* ChargedFraction */
    values[idx++].d = (stats.n_charged > 0) ? (double)(stats.n_positive - stats.n_negative) / stats.n_charged : 0.0; /* ChargeBalance */
    values[idx++].d = (stats.n_total > 0) ? stats.sum_ip / stats.n_total : 0.0;                    /* MeanIonizationPotential */
    values[idx++].d = (stats.n_heavy > 0) ? stats.sum_valence / stats.n_heavy : 0.0;               /* ValenceDensity */

    /* Group F: Topological Electronics (6) */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_donors / stats.n_heavy : 0.0;           /* DonorFraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_acceptors / stats.n_heavy : 0.0;        /* AcceptorFraction */
    values[idx++].d = (stats.n_acceptors > 0) ? (double)stats.n_donors / stats.n_acceptors : 0.0;   /* DonorAcceptorRatio */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_terminal / stats.n_heavy : 0.0;         /* TerminalFraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_branch / stats.n_heavy : 0.0;           /* BranchFraction */
    int cyclomatic = stats.n_bonds - stats.n_heavy + 1;
    values[idx++].d = (stats.n_heavy > 0) ? (double)cyclomatic / stats.n_heavy : 0.0;              /* TopoComplexity */

    /* Group G: Elemental Balance Ratios (5) */
    int n_plus_o = stats.n_nitrogen + stats.n_oxygen;
    values[idx++].d = (n_plus_o > 0) ? (double)stats.n_nitrogen / n_plus_o : 0.0;                  /* NitrogenOxygenBalance */
    values[idx++].d = (stats.n_heteroatom > 0) ? (double)stats.n_sulfur / stats.n_heteroatom : 0.0; /* SulfurHeteroFraction */
    values[idx++].d = (stats.n_halogen > 0) ? (double)stats.n_fluorine / stats.n_halogen : 0.0;    /* FluoroHalogenFraction */
    values[idx++].d = (stats.n_halogen > 0) ? (double)stats.n_chlorine / stats.n_halogen : 0.0;    /* ChlorineHalogenFraction */
    int hydrophobic = stats.n_carbon + stats.n_halogen;
    int hydrophilic = stats.n_nitrogen + stats.n_oxygen + stats.n_sulfur + stats.n_phosphorus;
    values[idx++].d = (hydrophilic > 0) ? (double)hydrophobic / hydrophilic : (hydrophobic > 0 ? 999.0 : 0.0); /* HydrophobicHydrophilicRatio */

    /* Group H: Bond Distribution Ratios (5) */
    values[idx++].d = (stats.n_bonds > 0) ? (double)stats.n_aromatic_bonds / stats.n_bonds : 0.0;  /* AromaticBondFraction */
    int acyclic_bonds = stats.n_bonds - stats.n_ring_bonds;
    values[idx++].d = (stats.n_bonds > 0) ? (double)acyclic_bonds / stats.n_bonds : 0.0;           /* AcyclicBondFraction */
    values[idx++].d = (stats.n_bonds > 0) ? (double)stats.n_single / stats.n_bonds : 0.0;          /* SingleBondFraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_rotatable / stats.n_heavy : 0.0;       /* RotatableDensity */
    /* Pi electron density: 2*double + 1.5*aromatic + 2*triple per heavy atom */
    double pi_density = stats.sum_pi_electrons;
    values[idx++].d = (stats.n_heavy > 0) ? pi_density / stats.n_heavy : 0.0;                      /* PiElectronDensity */

    /* Group I: Mass-Based Fractions (5) */
    values[idx++].d = (stats.sum_mass > 0) ? stats.sum_halogen_mass / stats.sum_mass : 0.0;        /* HalogenMassFraction */
    values[idx++].d = (stats.sum_mass > 0) ? stats.sum_hydrogen_mass / stats.sum_mass : 0.0;       /* HydrogenMassFraction */
    values[idx++].d = (stats.sum_mass > 0) ? stats.sum_polar_mass / stats.sum_mass : 0.0;          /* PolarMassFraction */
    values[idx++].d = (stats.sum_mass > 0) ? (double)stats.n_heavy / stats.sum_mass : 0.0;         /* HeavyAtomDensity */
    values[idx++].d = (stats.sum_mass > 0) ? stats.sum_vdw_volume / stats.sum_mass : 0.0;          /* SpecificVolume */

    /* Group J: Topological Connectivity Ratios (5) */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_degree2 / stats.n_heavy : 0.0;         /* LinkerAtomFraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_degree3plus / stats.n_heavy : 0.0;     /* CoreAtomFraction */
    values[idx++].d = (stats.n_ring_atoms + 1.0 > 0) ? (double)stats.n_terminal / (stats.n_ring_atoms + 1.0) : 0.0; /* ShapeIndex1 */
    int graph_cyclomatic = stats.n_bonds - stats.n_heavy + 1;
    values[idx++].d = (stats.n_heavy > 0) ? (double)graph_cyclomatic / stats.n_heavy : 0.0;        /* GraphCompactness */
    values[idx++].d = (stats.n_carbon > 0) ? (double)stats.n_chiral_carbons / stats.n_carbon : 0.0; /* ChiralCarbonFraction */

    /* Group K: Ring Environment Ratios (4) */
    values[idx++].d = (stats.total_ring_atoms > 0) ? (double)stats.n_fused_ring_atoms / stats.total_ring_atoms : 0.0; /* FusedRingAtomFraction */
    values[idx++].d = (stats.total_ring_atoms > 0) ? (double)stats.n_hetero_in_rings / stats.total_ring_atoms : 0.0;  /* HeteroInRingFraction */
    values[idx++].d = (stats.n_rings > 0) ? (double)stats.n_aromatic_rings / stats.n_rings : 0.0;  /* AromaticRingFraction */
    values[idx++].d = (stats.n_rings > 0) ? (double)stats.n_saturated_rings / stats.n_rings : 0.0; /* SaturatedRingFraction */

    /* Group L: Electronic & Solvation Environment (6) */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_conjugated / stats.n_heavy : 0.0;      /* ConjugatedAtomFraction */
    values[idx++].d = (stats.n_heavy > 0) ? (double)stats.n_lone_pairs / stats.n_heavy : 0.0;      /* LonePairDensity */
    values[idx++].d = (stats.n_polar_total > 0) ? (double)stats.n_polar_buried / stats.n_polar_total : 0.0; /* BuriedPolarFraction */
    values[idx++].d = (stats.n_lipophilic_total > 0) ? (double)stats.n_lipophilic_buried / stats.n_lipophilic_total : 0.0; /* BuriedLipophilicFraction */
    values[idx++].d = (stats.sum_vdw_volume > 0) ? stats.sum_tpsa / stats.sum_vdw_volume : 0.0;    /* PolarSurfaceDensity */
    double unsat_index = (stats.n_bonds > 0) ? (double)(2 * stats.n_double + 3 * stats.n_triple) / stats.n_bonds : 0.0;
    values[idx++].d = (stats.sum_mass > 0) ? unsat_index / stats.sum_mass * 100.0 : 0.0;           /* UnsaturationDensity */

    return NUM_RATIOS_DESCRIPTORS;
}

/* ============================================================================
 * Group A: Elemental Stoichiometry (5 descriptors)
 * ============================================================================ */

/* 1. Carbon/Heavy Atom Ratio */
static cchem_status_t desc_carbon_heavy_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_carbon / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 2. Heteroatom/Heavy Atom Ratio */
static cchem_status_t desc_hetero_heavy_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_heteroatom / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 3. Hydrogen/Carbon Ratio */
static cchem_status_t desc_hydrogen_carbon_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_carbon > 0) ? (double)stats.n_hydrogen / stats.n_carbon : 0.0;
    return CCHEM_OK;
}

/* 4. Halogen/Total Ratio */
static cchem_status_t desc_halogen_total_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_total > 0) ? (double)stats.n_halogen / stats.n_total : 0.0;
    return CCHEM_OK;
}

/* 5. Nitrogen+Oxygen/Heavy Ratio (Polar Atom Fraction) */
static cchem_status_t desc_polar_atom_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    int polar = stats.n_nitrogen + stats.n_oxygen;
    value->d = (stats.n_heavy > 0) ? (double)polar / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group B: Electronic Hybridization (5 descriptors)
 * ============================================================================ */

/* 6. sp3 Fraction */
static cchem_status_t desc_sp3_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_sp3 / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 7. sp2 Fraction */
static cchem_status_t desc_sp2_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_sp2 / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 8. sp Fraction */
static cchem_status_t desc_sp_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_sp / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 9. Aromatic Fraction */
static cchem_status_t desc_aromatic_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_aromatic / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 10. Ring Atom Fraction */
static cchem_status_t desc_ring_atom_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_ring_atoms / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group C: Electronegativity & Hardness (5 descriptors)
 * ============================================================================ */

/* 11. Mean Electronegativity */
static cchem_status_t desc_mean_electronegativity(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_total > 0) ? stats.sum_chi / stats.n_total : 0.0;
    return CCHEM_OK;
}

/* 12. Electronegativity Range */
static cchem_status_t desc_chi_range(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = stats.max_chi - stats.min_chi;
    return CCHEM_OK;
}

/* 13. Electronegativity Variance */
static cchem_status_t desc_chi_variance(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    if (stats.n_total > 0) {
        double mean = stats.sum_chi / stats.n_total;
        double variance = (stats.sum_chi_sq / stats.n_total) - (mean * mean);
        value->d = (variance > 0) ? variance : 0.0;
    } else {
        value->d = 0.0;
    }
    return CCHEM_OK;
}

/* 14. Hardness Index: (chi_max - chi_min) / chi_mean */
static cchem_status_t desc_hardness_index(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    double mean = (stats.n_total > 0) ? stats.sum_chi / stats.n_total : 0.0;
    value->d = (mean > 0.001) ? (stats.max_chi - stats.min_chi) / mean : 0.0;
    return CCHEM_OK;
}

/* 15. Electronegativity/Mass Ratio */
static cchem_status_t desc_chi_mass_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.sum_mass > 0) ? stats.sum_chi / stats.sum_mass : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group D: Bond Dynamics (5 descriptors)
 * ============================================================================ */

/* 16. Double Bond Fraction */
static cchem_status_t desc_double_bond_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_bonds > 0) ? (double)stats.n_double / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* 17. Triple Bond Fraction */
static cchem_status_t desc_triple_bond_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_bonds > 0) ? (double)stats.n_triple / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* 18. Rotatable Bond Fraction */
static cchem_status_t desc_rotatable_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_bonds > 0) ? (double)stats.n_rotatable / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* 19. Ring Bond Fraction */
static cchem_status_t desc_ring_bond_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_bonds > 0) ? (double)stats.n_ring_bonds / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* 20. Unsaturation Index: (2*double + 3*triple) / n_bonds */
static cchem_status_t desc_unsaturation_index(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    int weighted = 2 * stats.n_double + 3 * stats.n_triple;
    value->d = (stats.n_bonds > 0) ? (double)weighted / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group E: Charge & Ionization (4 descriptors)
 * ============================================================================ */

/* 21. Charged Atom Fraction */
static cchem_status_t desc_charged_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_charged / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 22. Charge Balance: (positive - negative) / total_charged */
static cchem_status_t desc_charge_balance(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    if (stats.n_charged > 0) {
        value->d = (double)(stats.n_positive - stats.n_negative) / stats.n_charged;
    } else {
        value->d = 0.0;
    }
    return CCHEM_OK;
}

/* 23. Mean Ionization Potential */
static cchem_status_t desc_mean_ip(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_total > 0) ? stats.sum_ip / stats.n_total : 0.0;
    return CCHEM_OK;
}

/* 24. Valence Electron Density: sum_valence / n_heavy */
static cchem_status_t desc_valence_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? stats.sum_valence / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group F: Topological Electronics (6 descriptors)
 * ============================================================================ */

/* 25. H-Bond Donor Fraction */
static cchem_status_t desc_donor_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_donors / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 26. H-Bond Acceptor Fraction */
static cchem_status_t desc_acceptor_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_acceptors / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 27. Donor/Acceptor Ratio */
static cchem_status_t desc_donor_acceptor_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_acceptors > 0) ? (double)stats.n_donors / stats.n_acceptors : 0.0;
    return CCHEM_OK;
}

/* 28. Terminal Atom Fraction */
static cchem_status_t desc_terminal_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_terminal / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 29. Branch Point Fraction */
static cchem_status_t desc_branch_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_branch / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 30. Topological Complexity: (n_bonds - n_atoms + 1) / n_atoms (cyclomatic density) */
static cchem_status_t desc_topo_complexity(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    /* Cyclomatic complexity = edges - nodes + connected_components */
    /* For a single molecule, connected_components = 1 */
    int cyclomatic = stats.n_bonds - stats.n_heavy + 1;
    value->d = (stats.n_heavy > 0) ? (double)cyclomatic / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group G: Elemental Balance Ratios (5 descriptors)
 * ============================================================================ */

/* 31. Nitrogen/Oxygen Balance */
static cchem_status_t desc_nitrogen_oxygen_balance(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    int n_plus_o = stats.n_nitrogen + stats.n_oxygen;
    value->d = (n_plus_o > 0) ? (double)stats.n_nitrogen / n_plus_o : 0.0;
    return CCHEM_OK;
}

/* 32. Sulfur/Heteroatom Fraction */
static cchem_status_t desc_sulfur_hetero_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heteroatom > 0) ? (double)stats.n_sulfur / stats.n_heteroatom : 0.0;
    return CCHEM_OK;
}

/* 33. Fluorine/Halogen Fraction */
static cchem_status_t desc_fluoro_halogen_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_halogen > 0) ? (double)stats.n_fluorine / stats.n_halogen : 0.0;
    return CCHEM_OK;
}

/* 34. Chlorine/Halogen Fraction */
static cchem_status_t desc_chlorine_halogen_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_halogen > 0) ? (double)stats.n_chlorine / stats.n_halogen : 0.0;
    return CCHEM_OK;
}

/* 35. Hydrophobic/Hydrophilic Ratio */
static cchem_status_t desc_hydrophobic_hydrophilic_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    int hydrophobic = stats.n_carbon + stats.n_halogen;
    int hydrophilic = stats.n_nitrogen + stats.n_oxygen + stats.n_sulfur + stats.n_phosphorus;
    value->d = (hydrophilic > 0) ? (double)hydrophobic / hydrophilic : (hydrophobic > 0 ? 999.0 : 0.0);
    return CCHEM_OK;
}

/* ============================================================================
 * Group H: Bond Distribution Ratios (5 descriptors)
 * ============================================================================ */

/* 36. Aromatic Bond Fraction */
static cchem_status_t desc_aromatic_bond_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_bonds > 0) ? (double)stats.n_aromatic_bonds / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* 37. Acyclic Bond Fraction */
static cchem_status_t desc_acyclic_bond_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    int acyclic = stats.n_bonds - stats.n_ring_bonds;
    value->d = (stats.n_bonds > 0) ? (double)acyclic / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* 38. Single Bond Fraction */
static cchem_status_t desc_single_bond_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_bonds > 0) ? (double)stats.n_single / stats.n_bonds : 0.0;
    return CCHEM_OK;
}

/* 39. Rotatable Bond Density */
static cchem_status_t desc_rotatable_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_rotatable / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 40. Pi Electron Density */
static cchem_status_t desc_pi_electron_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? stats.sum_pi_electrons / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group I: Mass-Based Fractions (5 descriptors)
 * ============================================================================ */

/* 41. Halogen Mass Fraction */
static cchem_status_t desc_halogen_mass_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.sum_mass > 0) ? stats.sum_halogen_mass / stats.sum_mass : 0.0;
    return CCHEM_OK;
}

/* 42. Hydrogen Mass Fraction */
static cchem_status_t desc_hydrogen_mass_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.sum_mass > 0) ? stats.sum_hydrogen_mass / stats.sum_mass : 0.0;
    return CCHEM_OK;
}

/* 43. Polar Mass Fraction */
static cchem_status_t desc_polar_mass_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.sum_mass > 0) ? stats.sum_polar_mass / stats.sum_mass : 0.0;
    return CCHEM_OK;
}

/* 44. Heavy Atom Density */
static cchem_status_t desc_heavy_atom_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.sum_mass > 0) ? (double)stats.n_heavy / stats.sum_mass : 0.0;
    return CCHEM_OK;
}

/* 45. Specific Volume (VdW Volume / MW) */
static cchem_status_t desc_specific_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.sum_mass > 0) ? stats.sum_vdw_volume / stats.sum_mass : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group J: Topological Connectivity Ratios (5 descriptors)
 * ============================================================================ */

/* 46. Linker Atom Fraction (degree-2 atoms) */
static cchem_status_t desc_linker_atom_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_degree2 / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 47. Core Atom Fraction (degree-3+ atoms) */
static cchem_status_t desc_core_atom_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_degree3plus / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 48. Shape Index 1: terminal / (ring_atoms + 1) */
static cchem_status_t desc_shape_index1(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (double)stats.n_terminal / (stats.n_ring_atoms + 1.0);
    return CCHEM_OK;
}

/* 49. Graph Compactness (cyclomatic complexity / heavy atoms) */
static cchem_status_t desc_graph_compactness(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    int cyclomatic = stats.n_bonds - stats.n_heavy + 1;
    value->d = (stats.n_heavy > 0) ? (double)cyclomatic / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 50. Chiral Carbon Fraction */
static cchem_status_t desc_chiral_carbon_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_carbon > 0) ? (double)stats.n_chiral_carbons / stats.n_carbon : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group K: Ring Environment Ratios (4 descriptors)
 * ============================================================================ */

/* 51. Fused Ring Atom Fraction */
static cchem_status_t desc_fused_ring_atom_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.total_ring_atoms > 0) ? (double)stats.n_fused_ring_atoms / stats.total_ring_atoms : 0.0;
    return CCHEM_OK;
}

/* 52. Heteroatom in Ring Fraction */
static cchem_status_t desc_hetero_in_ring_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.total_ring_atoms > 0) ? (double)stats.n_hetero_in_rings / stats.total_ring_atoms : 0.0;
    return CCHEM_OK;
}

/* 53. Aromatic Ring Fraction */
static cchem_status_t desc_aromatic_ring_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_rings > 0) ? (double)stats.n_aromatic_rings / stats.n_rings : 0.0;
    return CCHEM_OK;
}

/* 54. Saturated Ring Fraction */
static cchem_status_t desc_saturated_ring_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_rings > 0) ? (double)stats.n_saturated_rings / stats.n_rings : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group L: Electronic & Solvation Environment (6 descriptors)
 * ============================================================================ */

/* 55. Conjugated Atom Fraction */
static cchem_status_t desc_conjugated_atom_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_conjugated / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 56. Lone Pair Density */
static cchem_status_t desc_lone_pair_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_heavy > 0) ? (double)stats.n_lone_pairs / stats.n_heavy : 0.0;
    return CCHEM_OK;
}

/* 57. Buried Polar Fraction */
static cchem_status_t desc_buried_polar_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_polar_total > 0) ? (double)stats.n_polar_buried / stats.n_polar_total : 0.0;
    return CCHEM_OK;
}

/* 58. Buried Lipophilic Fraction */
static cchem_status_t desc_buried_lipophilic_fraction(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.n_lipophilic_total > 0) ? (double)stats.n_lipophilic_buried / stats.n_lipophilic_total : 0.0;
    return CCHEM_OK;
}

/* 59. Polar Surface Density (TPSA / VdW Volume) */
static cchem_status_t desc_polar_surface_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    value->d = (stats.sum_vdw_volume > 0) ? stats.sum_tpsa / stats.sum_vdw_volume : 0.0;
    return CCHEM_OK;
}

/* 60. Unsaturation Density */
static cchem_status_t desc_unsaturation_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    mol_stats_t stats;
    collect_mol_stats(mol, &stats);
    double unsat_index = (stats.n_bonds > 0) ? (double)(2 * stats.n_double + 3 * stats.n_triple) / stats.n_bonds : 0.0;
    value->d = (stats.sum_mass > 0) ? unsat_index / stats.sum_mass * 100.0 : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_RATIO_DESC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_RATIOS; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_ratios(void) {
    /* Group A: Elemental Stoichiometry */
    REGISTER_RATIO_DESC("CarbonHeavyRatio", "Carbon atoms / heavy atoms", desc_carbon_heavy_ratio);
    REGISTER_RATIO_DESC("HeteroHeavyRatio", "Heteroatoms / heavy atoms", desc_hetero_heavy_ratio);
    REGISTER_RATIO_DESC("HydrogenCarbonRatio", "Hydrogen / carbon atoms", desc_hydrogen_carbon_ratio);
    REGISTER_RATIO_DESC("HalogenTotalRatio", "Halogens / total atoms", desc_halogen_total_ratio);
    REGISTER_RATIO_DESC("PolarAtomRatio", "(N+O) / heavy atoms", desc_polar_atom_ratio);

    /* Group B: Electronic Hybridization */
    REGISTER_RATIO_DESC("Sp3Fraction", "sp3 carbons / heavy atoms", desc_sp3_fraction);
    REGISTER_RATIO_DESC("Sp2Fraction", "sp2 carbons / heavy atoms", desc_sp2_fraction);
    REGISTER_RATIO_DESC("SpFraction", "sp carbons / heavy atoms", desc_sp_fraction);
    REGISTER_RATIO_DESC("AromaticFraction", "Aromatic atoms / heavy atoms", desc_aromatic_fraction);
    REGISTER_RATIO_DESC("RingAtomFraction", "Ring atoms / heavy atoms", desc_ring_atom_fraction);

    /* Group C: Electronegativity & Hardness */
    REGISTER_RATIO_DESC("MeanElectronegativity", "Average Pauling electronegativity", desc_mean_electronegativity);
    REGISTER_RATIO_DESC("ChiRange", "Max - min electronegativity", desc_chi_range);
    REGISTER_RATIO_DESC("ChiVariance", "Electronegativity variance", desc_chi_variance);
    REGISTER_RATIO_DESC("HardnessIndex", "Chi range / mean chi", desc_hardness_index);
    REGISTER_RATIO_DESC("ChiMassRatio", "Total chi / molecular weight", desc_chi_mass_ratio);

    /* Group D: Bond Dynamics */
    REGISTER_RATIO_DESC("DoubleBondFraction", "Double bonds / total bonds", desc_double_bond_fraction);
    REGISTER_RATIO_DESC("TripleBondFraction", "Triple bonds / total bonds", desc_triple_bond_fraction);
    REGISTER_RATIO_DESC("RotatableFraction", "Rotatable bonds / total bonds", desc_rotatable_fraction);
    REGISTER_RATIO_DESC("RingBondFraction", "Ring bonds / total bonds", desc_ring_bond_fraction);
    REGISTER_RATIO_DESC("UnsaturationIndex", "Weighted unsaturation", desc_unsaturation_index);

    /* Group E: Charge & Ionization */
    REGISTER_RATIO_DESC("ChargedFraction", "Charged atoms / heavy atoms", desc_charged_fraction);
    REGISTER_RATIO_DESC("ChargeBalance", "Charge asymmetry index", desc_charge_balance);
    REGISTER_RATIO_DESC("MeanIonizationPotential", "Average first IP (eV)", desc_mean_ip);
    REGISTER_RATIO_DESC("ValenceDensity", "Valence electrons / heavy atoms", desc_valence_density);

    /* Group F: Topological Electronics */
    REGISTER_RATIO_DESC("DonorFraction", "H-bond donors / heavy atoms", desc_donor_fraction);
    REGISTER_RATIO_DESC("AcceptorFraction", "H-bond acceptors / heavy atoms", desc_acceptor_fraction);
    REGISTER_RATIO_DESC("DonorAcceptorRatio", "Donors / acceptors", desc_donor_acceptor_ratio);
    REGISTER_RATIO_DESC("TerminalFraction", "Terminal atoms / heavy atoms", desc_terminal_fraction);
    REGISTER_RATIO_DESC("BranchFraction", "Branch points / heavy atoms", desc_branch_fraction);
    REGISTER_RATIO_DESC("TopoComplexity", "Cyclomatic complexity density", desc_topo_complexity);

    /* Group G: Elemental Balance Ratios */
    REGISTER_RATIO_DESC("NitrogenOxygenBalance", "N / (N+O) balance", desc_nitrogen_oxygen_balance);
    REGISTER_RATIO_DESC("SulfurHeteroFraction", "Sulfur / heteroatoms", desc_sulfur_hetero_fraction);
    REGISTER_RATIO_DESC("FluoroHalogenFraction", "Fluorine / halogens", desc_fluoro_halogen_fraction);
    REGISTER_RATIO_DESC("ChlorineHalogenFraction", "Chlorine / halogens", desc_chlorine_halogen_fraction);
    REGISTER_RATIO_DESC("HydrophobicHydrophilicRatio", "(C+Hal) / (N+O+S+P)", desc_hydrophobic_hydrophilic_ratio);

    /* Group H: Bond Distribution Ratios */
    REGISTER_RATIO_DESC("AromaticBondFraction", "Aromatic bonds / total bonds", desc_aromatic_bond_fraction);
    REGISTER_RATIO_DESC("AcyclicBondFraction", "Non-ring bonds / total bonds", desc_acyclic_bond_fraction);
    REGISTER_RATIO_DESC("SingleBondFraction", "Single bonds / total bonds", desc_single_bond_fraction);
    REGISTER_RATIO_DESC("RotatableDensity", "Rotatable bonds / heavy atoms", desc_rotatable_density);
    REGISTER_RATIO_DESC("PiElectronDensity", "Pi electrons / heavy atoms", desc_pi_electron_density);

    /* Group I: Mass-Based Fractions */
    REGISTER_RATIO_DESC("HalogenMassFraction", "Halogen mass / MW", desc_halogen_mass_fraction);
    REGISTER_RATIO_DESC("HydrogenMassFraction", "Hydrogen mass / MW", desc_hydrogen_mass_fraction);
    REGISTER_RATIO_DESC("PolarMassFraction", "(N+O+S) mass / MW", desc_polar_mass_fraction);
    REGISTER_RATIO_DESC("HeavyAtomDensity", "Heavy atoms / MW", desc_heavy_atom_density);
    REGISTER_RATIO_DESC("SpecificVolume", "VdW volume / MW", desc_specific_volume);

    /* Group J: Topological Connectivity Ratios */
    REGISTER_RATIO_DESC("LinkerAtomFraction", "Degree-2 atoms / heavy atoms", desc_linker_atom_fraction);
    REGISTER_RATIO_DESC("CoreAtomFraction", "Degree-3+ atoms / heavy atoms", desc_core_atom_fraction);
    REGISTER_RATIO_DESC("ShapeIndex1", "Terminal / (ring atoms + 1)", desc_shape_index1);
    REGISTER_RATIO_DESC("GraphCompactness", "Cyclomatic / heavy atoms", desc_graph_compactness);
    REGISTER_RATIO_DESC("ChiralCarbonFraction", "Chiral C / total C", desc_chiral_carbon_fraction);

    /* Group K: Ring Environment Ratios */
    REGISTER_RATIO_DESC("FusedRingAtomFraction", "Fused ring atoms / ring atoms", desc_fused_ring_atom_fraction);
    REGISTER_RATIO_DESC("HeteroInRingFraction", "Hetero in rings / ring atoms", desc_hetero_in_ring_fraction);
    REGISTER_RATIO_DESC("AromaticRingFraction", "Aromatic rings / total rings", desc_aromatic_ring_fraction);
    REGISTER_RATIO_DESC("SaturatedRingFraction", "Saturated rings / total rings", desc_saturated_ring_fraction);

    /* Group L: Electronic & Solvation Environment */
    REGISTER_RATIO_DESC("ConjugatedAtomFraction", "Conjugated atoms / heavy atoms", desc_conjugated_atom_fraction);
    REGISTER_RATIO_DESC("LonePairDensity", "Lone pairs / heavy atoms", desc_lone_pair_density);
    REGISTER_RATIO_DESC("BuriedPolarFraction", "Buried polar / total polar", desc_buried_polar_fraction);
    REGISTER_RATIO_DESC("BuriedLipophilicFraction", "Buried lipophilic / total lipophilic", desc_buried_lipophilic_fraction);
    REGISTER_RATIO_DESC("PolarSurfaceDensity", "TPSA / VdW volume", desc_polar_surface_density);
    REGISTER_RATIO_DESC("UnsaturationDensity", "Unsaturation / MW", desc_unsaturation_density);
}
