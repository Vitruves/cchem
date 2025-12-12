/**
 * @file fractional.c
 * @brief Fractional molecular descriptors
 *
 * Ultra-fast fractional descriptors representing proportions and ratios
 * of molecular properties. All descriptors are computed in O(n) time
 * with a single pass through atoms and bonds.
 *
 * Groups:
 * - Element MW fractions (FcC, FcN, FcO, FcS, FcF, FcCl, FcBr, FcI, FcHalo, FcHetero)
 * - Electronegativity fractions (FcPolar, FcApolar, FcENAboveAvg, FcENHigh, etc.)
 * - Bond type fractions (FcCSp3, FcCSp2, FcPol, FcUnpol, FcBondC, FcBondN, etc.)
 * - Structural fractions (FcAromaticAtoms, FcRingAtoms, FcChargedAtoms, etc.)
 * - Physical property fractions (FcSmallR, FcLargeR, FcHighPolz, FcSmallVdW, etc.)
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Physical Constants
 * ============================================================================ */

/* Atomic weights (g/mol) */
static const double ATOMIC_WEIGHT[] = {
    [ELEM_H]  = 1.008,
    [ELEM_C]  = 12.011,
    [ELEM_N]  = 14.007,
    [ELEM_O]  = 15.999,
    [ELEM_F]  = 18.998,
    [ELEM_P]  = 30.974,
    [ELEM_S]  = 32.065,
    [ELEM_Cl] = 35.453,
    [ELEM_Br] = 79.904,
    [ELEM_I]  = 126.90,
    [ELEM_Si] = 28.086,
    [ELEM_B]  = 10.811,
    [ELEM_Se] = 78.96,
    [ELEM_As] = 74.922,
};

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

/* Covalent radii in Angstroms */
static const double COVALENT_RADIUS[] = {
    [ELEM_H]  = 0.31,
    [ELEM_C]  = 0.77,
    [ELEM_N]  = 0.71,
    [ELEM_O]  = 0.66,
    [ELEM_F]  = 0.57,
    [ELEM_P]  = 1.07,
    [ELEM_S]  = 1.05,
    [ELEM_Cl] = 1.02,
    [ELEM_Br] = 1.20,
    [ELEM_I]  = 1.39,
    [ELEM_Si] = 1.11,
    [ELEM_B]  = 0.84,
    [ELEM_Se] = 1.20,
    [ELEM_As] = 1.19,
};

/* VdW volumes in Angstrom^3 (approximate spherical) */
static const double VDW_VOLUME[] = {
    [ELEM_H]  = 7.24,
    [ELEM_C]  = 20.58,
    [ELEM_N]  = 15.60,
    [ELEM_O]  = 14.71,
    [ELEM_F]  = 13.31,
    [ELEM_P]  = 24.43,
    [ELEM_S]  = 24.43,
    [ELEM_Cl] = 22.45,
    [ELEM_Br] = 26.52,
    [ELEM_I]  = 32.52,
    [ELEM_Si] = 38.79,
    [ELEM_B]  = 29.82,
    [ELEM_Se] = 28.73,
    [ELEM_As] = 26.52,
};

/* Atomic polarizabilities in Angstrom^3 */
static const double POLARIZABILITY[] = {
    [ELEM_H]  = 0.667,
    [ELEM_C]  = 1.76,
    [ELEM_N]  = 1.10,
    [ELEM_O]  = 0.802,
    [ELEM_F]  = 0.557,
    [ELEM_P]  = 3.63,
    [ELEM_S]  = 2.90,
    [ELEM_Cl] = 2.18,
    [ELEM_Br] = 3.05,
    [ELEM_I]  = 4.70,
    [ELEM_Si] = 5.38,
    [ELEM_B]  = 3.03,
    [ELEM_Se] = 3.77,
    [ELEM_As] = 4.31,
};

/* First ionization energy in eV */
static const double IONIZATION_ENERGY[] = {
    [ELEM_H]  = 13.60,
    [ELEM_C]  = 11.26,
    [ELEM_N]  = 14.53,
    [ELEM_O]  = 13.62,
    [ELEM_F]  = 17.42,
    [ELEM_P]  = 10.49,
    [ELEM_S]  = 10.36,
    [ELEM_Cl] = 12.97,
    [ELEM_Br] = 11.81,
    [ELEM_I]  = 10.45,
    [ELEM_Si] = 8.15,
    [ELEM_B]  = 8.30,
    [ELEM_Se] = 9.75,
    [ELEM_As] = 9.79,
};

/* Electron affinity in eV */
static const double ELECTRON_AFFINITY[] = {
    [ELEM_H]  = 0.75,
    [ELEM_C]  = 1.26,
    [ELEM_N]  = -0.07,  /* N has negative EA */
    [ELEM_O]  = 1.46,
    [ELEM_F]  = 3.40,
    [ELEM_P]  = 0.75,
    [ELEM_S]  = 2.08,
    [ELEM_Cl] = 3.61,
    [ELEM_Br] = 3.36,
    [ELEM_I]  = 3.06,
    [ELEM_Si] = 1.39,
    [ELEM_B]  = 0.28,
    [ELEM_Se] = 2.02,
    [ELEM_As] = 0.81,
};

/* Common maximum oxidation states (approximate) */
static const int MAX_OXIDATION[] = {
    [ELEM_H]  = 1,
    [ELEM_C]  = 4,
    [ELEM_N]  = 5,
    [ELEM_O]  = 2,
    [ELEM_F]  = 1,
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

/* Periodic table group for group 16 check */
static const int PERIODIC_GROUP[] = {
    [ELEM_H]  = 1,
    [ELEM_C]  = 14,
    [ELEM_N]  = 15,
    [ELEM_O]  = 16,
    [ELEM_F]  = 17,
    [ELEM_P]  = 15,
    [ELEM_S]  = 16,
    [ELEM_Cl] = 17,
    [ELEM_Br] = 17,
    [ELEM_I]  = 17,
    [ELEM_Si] = 14,
    [ELEM_B]  = 13,
    [ELEM_Se] = 16,
    [ELEM_As] = 15,
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

/* Metalloid check */
static inline bool is_metalloid(element_t elem) {
    return elem == ELEM_B || elem == ELEM_Si || elem == ELEM_As || elem == ELEM_Se;
}

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_atomic_weight(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ATOMIC_WEIGHT)/sizeof(ATOMIC_WEIGHT[0])))
        return 12.011;
    double w = ATOMIC_WEIGHT[elem];
    return (w > 0.0) ? w : 12.011;
}

static inline double get_electronegativity(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ELECTRONEGATIVITY)/sizeof(ELECTRONEGATIVITY[0])))
        return 2.55;
    double chi = ELECTRONEGATIVITY[elem];
    return (chi > 0.0) ? chi : 2.55;
}

static inline double get_covalent_radius(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(COVALENT_RADIUS)/sizeof(COVALENT_RADIUS[0])))
        return 0.77;
    double r = COVALENT_RADIUS[elem];
    return (r > 0.0) ? r : 0.77;
}

static inline double get_vdw_volume(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(VDW_VOLUME)/sizeof(VDW_VOLUME[0])))
        return 20.58;
    double v = VDW_VOLUME[elem];
    return (v > 0.0) ? v : 20.58;
}

static inline double get_polarizability(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(POLARIZABILITY)/sizeof(POLARIZABILITY[0])))
        return 1.76;
    double p = POLARIZABILITY[elem];
    return (p > 0.0) ? p : 1.76;
}

static inline double get_ionization_energy(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(IONIZATION_ENERGY)/sizeof(IONIZATION_ENERGY[0])))
        return 11.26;
    double ie = IONIZATION_ENERGY[elem];
    return (ie > 0.0) ? ie : 11.26;
}

static inline double get_electron_affinity(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ELECTRON_AFFINITY)/sizeof(ELECTRON_AFFINITY[0])))
        return 1.26;
    return ELECTRON_AFFINITY[elem];
}

static inline int get_max_oxidation(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(MAX_OXIDATION)/sizeof(MAX_OXIDATION[0])))
        return 4;
    return MAX_OXIDATION[elem];
}

static inline int get_periodic_group(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(PERIODIC_GROUP)/sizeof(PERIODIC_GROUP[0])))
        return 14;
    return PERIODIC_GROUP[elem];
}

static inline int get_valence_electrons(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(VALENCE_ELECTRONS)/sizeof(VALENCE_ELECTRONS[0])))
        return 4;
    return VALENCE_ELECTRONS[elem];
}

/* Check if atom is a halogen */
static inline bool is_halogen(element_t elem) {
    return elem == ELEM_F || elem == ELEM_Cl || elem == ELEM_Br || elem == ELEM_I;
}

/* Check if atom is heteroatom (non-C, non-H) */
static inline bool is_heteroatom(element_t elem) {
    return elem != ELEM_C && elem != ELEM_H;
}

/* Check if atom is heavy (Z > 10) */
static inline bool is_heavy_z(element_t elem) {
    return elem == ELEM_P || elem == ELEM_S || elem == ELEM_Cl ||
           elem == ELEM_Br || elem == ELEM_I || elem == ELEM_Si ||
           elem == ELEM_Se || elem == ELEM_As;
}

/* Get hybridization from connectivity */
static inline int get_hybridization(const molecule_t* mol, int atom_idx) {
    const atom_t* a = &mol->atoms[atom_idx];
    if (a->element == ELEM_H) return 0;

    int double_bonds = 0, triple_bonds = 0;
    for (int i = 0; i < a->num_neighbors; i++) {
        bond_type_t bt = mol->bonds[a->neighbor_bonds[i]].type;
        if (bt == BOND_DOUBLE) double_bonds++;
        else if (bt == BOND_TRIPLE) triple_bonds++;
    }

    if (triple_bonds > 0 || double_bonds >= 2) return 1;  /* sp */
    if (double_bonds == 1 || a->aromatic) return 2;        /* sp2 */
    return 3;  /* sp3 */
}

/* Check if atom is H-bond donor */
static inline bool is_hb_donor(const molecule_t* mol, int atom_idx) {
    const atom_t* a = &mol->atoms[atom_idx];
    element_t e = a->element;
    if (e != ELEM_N && e != ELEM_O && e != ELEM_S) return false;

    /* Check for attached H */
    if (a->implicit_h_count > 0) return true;
    for (int i = 0; i < a->num_neighbors; i++) {
        if (mol->atoms[a->neighbors[i]].element == ELEM_H) return true;
    }
    return false;
}

/* Check if atom is H-bond acceptor */
static inline bool is_hb_acceptor(const atom_t* a) {
    element_t e = a->element;
    return e == ELEM_N || e == ELEM_O || e == ELEM_F || e == ELEM_S;
}

/* ============================================================================
 * Batch Computation - All fractional descriptors in single pass
 * ============================================================================ */

#define NUM_FRAC_DESCRIPTORS 61

typedef struct {
    /* MW components */
    double total_mw;
    double mw_C, mw_N, mw_O, mw_S, mw_F, mw_Cl, mw_Br, mw_I;
    double mw_halo, mw_hetero;

    /* Atom counts */
    int n_atoms, n_heavy;
    int n_C, n_N, n_O, n_S, n_P;
    int n_polar, n_apolar;
    int n_aromatic, n_ring;
    int n_charged, n_pos, n_neg;
    int n_hbd, n_hba;
    int n_bridge;
    int n_sp, n_sp2, n_sp3;
    int n_hetero, n_halo;
    int n_metalloid;

    /* EN-based counts */
    int n_en_high;        /* EN > 3.5 */
    int n_en_above_avg;   /* EN > molecule average */
    int n_en_below_avg;

    /* Size-based counts */
    int n_small_r;        /* covalent radius < 0.77 */
    int n_large_r;        /* covalent radius > 1.1 */
    int n_small_vdw;      /* VdW volume < 15 */
    int n_large_vdw;      /* VdW volume > 25 */

    /* Polarizability counts */
    int n_low_polz;       /* polarizability < 1.0 (low threshold for atoms) */
    int n_high_polz;      /* polarizability > 3.5 */

    /* EA/IE counts */
    int n_high_ea;        /* EA > 2 */
    int n_low_ea;         /* EA < 0.5 */
    int n_ie_odd;         /* int(IE) is odd */
    int n_low_ie_heavy;   /* heavy atoms with IE < 11 */

    /* Valence counts */
    int n_even_valence;

    /* Oxidation state proxies */
    int n_high_ox;        /* max ox state >= 4 */
    int n_ox_en_above;    /* abs(ox) * EN > 10 */

    /* Mixed property counts */
    int n_hetero_polar;   /* heteroatoms that are polar */
    int n_halo_polar_bond;
    int n_heavy_polar_bond;
    int n_sp3_heavy;
    int n_sp2_en_above;
    int n_nonhetero_high_ea;
    int n_hetero_low_polz;
    int n_group16;
    int n_sp_high_en;
    int n_ring_high_ox;
    int n_heavy_formal;
    int n_sp3_high_polz;
    int n_heavy_ox_neg;
    int n_radius_mw_above;

    /* Bond counts */
    int n_bonds;
    int n_polar_bonds;
    int n_unpolar_bonds;
    int n_cc_sp3, n_cc_sp2, n_cc_total;
    int n_C_nonsingle, n_C_total_bonds;
    int n_N_nonsingle, n_N_total_bonds;
    int n_O_nonsingle, n_O_total_bonds;
    int n_S_nonsingle, n_S_total_bonds;
    int n_P_nonsingle, n_P_total_bonds;
    int n_en_bonded;      /* bonds where both atoms have EN > 3.0 */
    int total_bond_order;

    /* Sums for weighted descriptors */
    double sum_en;
    double sum_en_mw;     /* sum(EN * atom_MW) */
    double sum_vdw_mw;
    double sum_rcov_mw;
    double avg_en;        /* for computing above/below average */
} frac_stats_t;

static void collect_frac_stats(const molecule_t* mol, frac_stats_t* s) {
    memset(s, 0, sizeof(frac_stats_t));

    const int n_atoms = mol->num_atoms;
    const int n_bonds = mol->num_bonds;

    /* First pass: collect basic counts and sums */
    for (int i = 0; i < n_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        element_t e = a->element;

        if (e == ELEM_H) continue;  /* Skip explicit H for most calculations */

        s->n_atoms++;
        s->n_heavy++;

        double mw = get_atomic_weight(e);
        double en = get_electronegativity(e);
        double rcov = get_covalent_radius(e);
        double vdw = get_vdw_volume(e);
        double polz = get_polarizability(e);
        double ie = get_ionization_energy(e);
        double ea = get_electron_affinity(e);
        int ox = get_max_oxidation(e);
        int val = get_valence_electrons(e);
        int group = get_periodic_group(e);

        s->total_mw += mw;
        s->sum_en += en;
        s->sum_en_mw += en * mw;
        s->sum_vdw_mw += vdw * mw;
        s->sum_rcov_mw += rcov * mw;

        /* Element-specific MW */
        switch (e) {
            case ELEM_C:  s->mw_C += mw; s->n_C++; break;
            case ELEM_N:  s->mw_N += mw; s->n_N++; s->mw_hetero += mw; s->n_hetero++; break;
            case ELEM_O:  s->mw_O += mw; s->n_O++; s->mw_hetero += mw; s->n_hetero++; break;
            case ELEM_S:  s->mw_S += mw; s->n_S++; s->mw_hetero += mw; s->n_hetero++; break;
            case ELEM_P:  s->n_P++; s->mw_hetero += mw; s->n_hetero++; break;
            case ELEM_F:  s->mw_F += mw; s->mw_halo += mw; s->mw_hetero += mw; s->n_halo++; s->n_hetero++; break;
            case ELEM_Cl: s->mw_Cl += mw; s->mw_halo += mw; s->mw_hetero += mw; s->n_halo++; s->n_hetero++; break;
            case ELEM_Br: s->mw_Br += mw; s->mw_halo += mw; s->mw_hetero += mw; s->n_halo++; s->n_hetero++; break;
            case ELEM_I:  s->mw_I += mw; s->mw_halo += mw; s->mw_hetero += mw; s->n_halo++; s->n_hetero++; break;
            default:      if (is_heteroatom(e)) { s->mw_hetero += mw; s->n_hetero++; } break;
        }

        /* EN-based: polar = EN > C (2.55) */
        if (en > 2.55) s->n_polar++;
        else s->n_apolar++;

        if (en > 3.5) s->n_en_high++;

        /* Size-based */
        if (rcov < 0.77) s->n_small_r++;
        if (rcov > 1.1) s->n_large_r++;
        if (vdw < 15.0) s->n_small_vdw++;
        if (vdw > 25.0) s->n_large_vdw++;

        /* Polarizability */
        if (polz < 1.0) s->n_low_polz++;
        if (polz > 3.5) s->n_high_polz++;

        /* EA/IE */
        if (ea > 2.0) s->n_high_ea++;
        if (ea < 0.5) s->n_low_ea++;
        if (((int)ie) % 2 == 1) s->n_ie_odd++;
        if (is_heavy_z(e) && ie < 11.0) s->n_low_ie_heavy++;

        /* Valence */
        if (val % 2 == 0) s->n_even_valence++;

        /* Oxidation */
        if (ox >= 4) s->n_high_ox++;
        if (abs(ox) * en > 10.0) s->n_ox_en_above++;

        /* Structural */
        if (a->aromatic) s->n_aromatic++;
        if (a->ring_count > 0) s->n_ring++;
        if (a->charge != 0) {
            s->n_charged++;
            if (a->charge > 0) s->n_pos++;
            else s->n_neg++;
        }

        /* H-bonding */
        if (is_hb_donor(mol, i)) s->n_hbd++;
        if (is_hb_acceptor(a)) s->n_hba++;

        /* Bridgehead: in ring with degree > 2 */
        if (a->ring_count > 0 && a->num_neighbors > 2) {
            int ring_neighbors = 0;
            for (int j = 0; j < a->num_neighbors; j++) {
                if (mol->atoms[a->neighbors[j]].ring_count > 0) ring_neighbors++;
            }
            if (ring_neighbors > 2) s->n_bridge++;
        }

        /* Hybridization */
        int hyb = get_hybridization(mol, i);
        if (hyb == 1) s->n_sp++;
        else if (hyb == 2) s->n_sp2++;
        else if (hyb == 3) s->n_sp3++;

        /* Metalloid */
        if (is_metalloid(e)) s->n_metalloid++;

        /* Group 16 */
        if (group == 16) s->n_group16++;

        /* Mixed property counts */
        if (is_heteroatom(e) && en > 2.55) s->n_hetero_polar++;
        if (is_heavy_z(e) && hyb == 3) s->n_sp3_heavy++;
        if (hyb == 2 && en > 2.55) s->n_sp2_en_above++;
        if (!is_heteroatom(e) && ea > 1.0) s->n_nonhetero_high_ea++;
        if (is_heteroatom(e) && polz < 1.5) s->n_hetero_low_polz++;
        if (hyb == 1 && en > 2.5) s->n_sp_high_en++;
        if (a->ring_count > 0 && ox > 2) s->n_ring_high_ox++;
        if (is_heavy_z(e) && a->charge != 0) s->n_heavy_formal++;
        if (hyb == 3 && polz > 3.5) s->n_sp3_high_polz++;
        if (is_heavy_z(e) && ox < 0) s->n_heavy_ox_neg++;
        if (rcov / mw > 0.07) s->n_radius_mw_above++;
    }

    /* Add implicit H to total MW */
    for (int i = 0; i < n_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        s->total_mw += a->implicit_h_count * 1.008;
        s->n_atoms += a->implicit_h_count;
    }

    /* Compute average EN for above/below calculations */
    s->avg_en = (s->n_heavy > 0) ? s->sum_en / s->n_heavy : 2.55;

    /* Second pass: EN above/below average */
    for (int i = 0; i < n_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        if (a->element == ELEM_H) continue;
        double en = get_electronegativity(a->element);
        if (en > s->avg_en) s->n_en_above_avg++;
        else s->n_en_below_avg++;
    }

    /* Bond pass */
    s->n_bonds = n_bonds;
    for (int i = 0; i < n_bonds; i++) {
        const bond_t* b = &mol->bonds[i];
        element_t e1 = mol->atoms[b->atom1].element;
        element_t e2 = mol->atoms[b->atom2].element;

        if (e1 == ELEM_H || e2 == ELEM_H) continue;

        double en1 = get_electronegativity(e1);
        double en2 = get_electronegativity(e2);

        /* Polar vs unpolar bonds */
        if (fabs(en1 - en2) < 0.4) s->n_unpolar_bonds++;
        else s->n_polar_bonds++;

        /* Both atoms EN > 3.0 */
        if (en1 > 3.0 && en2 > 3.0) s->n_en_bonded++;

        /* C-C bond types */
        if (e1 == ELEM_C && e2 == ELEM_C) {
            s->n_cc_total++;
            int hyb1 = get_hybridization(mol, b->atom1);
            int hyb2 = get_hybridization(mol, b->atom2);
            if (hyb1 == 3 && hyb2 == 3) s->n_cc_sp3++;
            else if (hyb1 == 2 || hyb2 == 2) s->n_cc_sp2++;
        }

        /* Bond order */
        int order = 1;
        if (b->type == BOND_DOUBLE) order = 2;
        else if (b->type == BOND_TRIPLE) order = 3;
        else if (b->type == BOND_AROMATIC) order = 1;  /* count as 1.5 later */
        s->total_bond_order += order;

        bool is_nonsingle = (b->type != BOND_SINGLE);

        /* Element-specific non-single bonds */
        if (e1 == ELEM_C || e2 == ELEM_C) {
            s->n_C_total_bonds++;
            if (is_nonsingle) s->n_C_nonsingle++;
        }
        if (e1 == ELEM_N || e2 == ELEM_N) {
            s->n_N_total_bonds++;
            if (is_nonsingle) s->n_N_nonsingle++;
        }
        if (e1 == ELEM_O || e2 == ELEM_O) {
            s->n_O_total_bonds++;
            if (is_nonsingle) s->n_O_nonsingle++;
        }
        if (e1 == ELEM_S || e2 == ELEM_S) {
            s->n_S_total_bonds++;
            if (is_nonsingle) s->n_S_nonsingle++;
        }
        if (e1 == ELEM_P || e2 == ELEM_P) {
            s->n_P_total_bonds++;
            if (is_nonsingle) s->n_P_nonsingle++;
        }

        /* Heavy atom polar bonds */
        if ((is_heavy_z(e1) || is_heavy_z(e2)) && fabs(en1 - en2) >= 0.4) {
            s->n_heavy_polar_bond++;
        }

        /* Halogen polar bonds */
        if ((is_halogen(e1) || is_halogen(e2)) && fabs(en1 - en2) >= 0.4) {
            s->n_halo_polar_bond++;
        }
    }
}

int descriptors_compute_fractional_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    frac_stats_t s;
    collect_frac_stats(mol, &s);

    int idx = 0;
    double safe_mw = (s.total_mw > 0) ? s.total_mw : 1.0;
    double safe_atoms = (s.n_atoms > 0) ? (double)s.n_atoms : 1.0;
    double safe_heavy = (s.n_heavy > 0) ? (double)s.n_heavy : 1.0;
    double safe_bonds = (s.n_bonds > 0) ? (double)s.n_bonds : 1.0;

    /* MW fractions */
    values[idx++].d = s.mw_C / safe_mw;          /* FcC */
    values[idx++].d = s.mw_N / safe_mw;          /* FcN */
    values[idx++].d = s.mw_O / safe_mw;          /* FcO */
    values[idx++].d = s.mw_S / safe_mw;          /* FcS */
    values[idx++].d = s.mw_F / safe_mw;          /* FcF */
    values[idx++].d = s.mw_Cl / safe_mw;         /* FcCl */
    values[idx++].d = s.mw_Br / safe_mw;         /* FcBr */
    values[idx++].d = s.mw_I / safe_mw;          /* FcI */
    values[idx++].d = s.mw_halo / safe_mw;       /* FcHalo */
    values[idx++].d = s.mw_hetero / safe_mw;     /* FcHetero */

    /* EN-based fractions */
    values[idx++].d = s.n_polar / safe_heavy;         /* FcPolar */
    values[idx++].d = s.n_apolar / safe_heavy;        /* FcApolar */
    values[idx++].d = s.n_en_above_avg / safe_heavy;  /* FcENAboveAvg */
    values[idx++].d = s.n_en_below_avg / safe_heavy;  /* FcENBelowAvg */
    values[idx++].d = s.n_en_high / safe_heavy;       /* FcENHigh */

    /* Bond type fractions */
    double safe_cc = (s.n_cc_total > 0) ? (double)s.n_cc_total : 1.0;
    values[idx++].d = s.n_cc_sp3 / safe_cc;           /* FcCSp3 */
    values[idx++].d = s.n_cc_sp2 / safe_cc;           /* FcCSp2 */
    values[idx++].d = s.n_unpolar_bonds / safe_bonds; /* FcUnpol */
    values[idx++].d = s.n_polar_bonds / safe_bonds;   /* FcPol */

    /* Electronegativity averages */
    values[idx++].d = s.sum_en / safe_heavy;           /* FcSumPolAt */
    values[idx++].d = s.sum_en / safe_mw;              /* FcSumPolMW */

    /* Bond density */
    values[idx++].d = s.total_bond_order / safe_atoms; /* FcBondAt */

    /* Element-specific bond fractions */
    values[idx++].d = (s.n_N_total_bonds > 0) ? (double)s.n_N_nonsingle / s.n_N_total_bonds : 0.0;  /* FcBondN */
    values[idx++].d = (s.n_O_total_bonds > 0) ? (double)s.n_O_nonsingle / s.n_O_total_bonds : 0.0;  /* FcBondO */
    values[idx++].d = (s.n_C_total_bonds > 0) ? (double)s.n_C_nonsingle / s.n_C_total_bonds : 0.0;  /* FcBondC */
    values[idx++].d = (s.n_S_total_bonds > 0) ? (double)s.n_S_nonsingle / s.n_S_total_bonds : 0.0;  /* FcBondS */
    values[idx++].d = (s.n_P_total_bonds > 0) ? (double)s.n_P_nonsingle / s.n_P_total_bonds : 0.0;  /* FcBondP */

    /* Structural fractions */
    values[idx++].d = s.n_hbd / safe_heavy;       /* FcHBDonors */
    values[idx++].d = s.n_hba / safe_heavy;       /* FcHBAcceptors */
    values[idx++].d = s.n_aromatic / safe_heavy;  /* FcAromaticAtoms */
    values[idx++].d = s.n_ring / safe_heavy;      /* FcRingAtoms */
    values[idx++].d = s.n_bridge / safe_heavy;    /* FcBridgeAtoms */
    values[idx++].d = s.n_charged / safe_heavy;   /* FcChargedAtoms */

    /* Physical property fractions */
    values[idx++].d = s.n_small_r / safe_heavy;   /* FcSmallR */
    values[idx++].d = s.n_large_r / safe_heavy;   /* FcLargeR */
    values[idx++].d = s.n_low_polz / safe_heavy;  /* FcLowPolz */
    values[idx++].d = s.n_high_polz / safe_heavy; /* FcHighPolz */
    values[idx++].d = s.n_high_ea / safe_heavy;   /* FcHighEA */
    values[idx++].d = s.n_low_ea / safe_heavy;    /* FcLowEA */
    values[idx++].d = s.n_small_vdw / safe_heavy; /* FcSmallVdW */
    values[idx++].d = s.n_large_vdw / safe_heavy; /* FcLargeVdW */

    /* MW-weighted descriptors */
    values[idx++].d = s.sum_en_mw / safe_mw;      /* FcENMW */
    values[idx++].d = s.n_en_bonded / safe_bonds; /* FcENBonded */
    values[idx++].d = s.sum_vdw_mw / safe_mw;     /* FcVdWMW */
    values[idx++].d = s.sum_rcov_mw / safe_mw;    /* FcRcovMW */

    /* Oxidation and charge fractions */
    values[idx++].d = s.n_high_ox / safe_heavy;           /* FcHighOxState */
    values[idx++].d = s.n_metalloid / safe_heavy;         /* FcMetalloid */
    values[idx++].d = s.n_hetero_polar / safe_heavy;      /* FcHETpol */
    values[idx++].d = (s.n_halo > 0) ? (double)s.n_halo_polar_bond / s.n_halo : 0.0;  /* FcHALpol */

    /* Heavy atom fractions */
    int n_heavy_z = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (is_heavy_z(mol->atoms[i].element)) n_heavy_z++;
    }
    double safe_heavy_z = (n_heavy_z > 0) ? (double)n_heavy_z : 1.0;
    values[idx++].d = (n_heavy_z > 0) ? (double)s.n_heavy_polar_bond / n_heavy_z : 0.0;  /* FcHeavyPol */
    values[idx++].d = s.n_sp3_heavy / safe_heavy_z;       /* FcSp3HeavyAtoms */

    /* Mixed fractions */
    double safe_sp2 = (s.n_sp2 > 0) ? (double)s.n_sp2 : 1.0;
    values[idx++].d = s.n_sp2_en_above / safe_sp2;        /* FcSp2ENAboveAvg */
    values[idx++].d = s.n_ie_odd / safe_heavy;            /* FcIEOdd */
    values[idx++].d = s.n_even_valence / safe_heavy;      /* FcEvenValenceAtoms */
    values[idx++].d = s.n_nonhetero_high_ea / safe_heavy; /* FcNonHeteroHighEA */
    values[idx++].d = (n_heavy_z > 0) ? (double)s.n_low_ie_heavy / n_heavy_z : 0.0;  /* FcHeavyLowIE */

    double safe_hetero = (s.n_hetero > 0) ? (double)s.n_hetero : 1.0;
    values[idx++].d = s.n_hetero_low_polz / safe_hetero;  /* FcHeteroLowPolz */
    values[idx++].d = s.n_ox_en_above / safe_heavy;       /* FcOxENAboveThreshold */

    /* Formal charge fractions */
    values[idx++].d = s.n_charged / safe_heavy;           /* FcFormalChargeNonZero */
    values[idx++].d = s.n_pos / safe_heavy;               /* FcFormalChargePositive */
    values[idx++].d = s.n_neg / safe_heavy;               /* FcFormalChargeNegative */

    return idx;
}

/* ============================================================================
 * Individual Descriptor Functions
 * ============================================================================ */

/* For registration, we create individual compute functions that extract
 * specific values. For efficiency in batch mode, use compute_fractional_all. */

#define FRAC_DESC_FN(name, stat_idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    descriptor_value_t vals[NUM_FRAC_DESCRIPTORS]; \
    int n = descriptors_compute_fractional_all(mol, vals); \
    if (n < 0 || stat_idx >= n) return CCHEM_ERROR_INVALID_INPUT; \
    *value = vals[stat_idx]; \
    return CCHEM_OK; \
}

FRAC_DESC_FN(fc_c, 0)
FRAC_DESC_FN(fc_n, 1)
FRAC_DESC_FN(fc_o, 2)
FRAC_DESC_FN(fc_s, 3)
FRAC_DESC_FN(fc_f, 4)
FRAC_DESC_FN(fc_cl, 5)
FRAC_DESC_FN(fc_br, 6)
FRAC_DESC_FN(fc_i, 7)
FRAC_DESC_FN(fc_halo, 8)
FRAC_DESC_FN(fc_hetero, 9)
FRAC_DESC_FN(fc_polar, 10)
FRAC_DESC_FN(fc_apolar, 11)
FRAC_DESC_FN(fc_en_above_avg, 12)
FRAC_DESC_FN(fc_en_below_avg, 13)
FRAC_DESC_FN(fc_en_high, 14)
FRAC_DESC_FN(fc_csp3, 15)
FRAC_DESC_FN(fc_csp2, 16)
FRAC_DESC_FN(fc_unpol, 17)
FRAC_DESC_FN(fc_pol, 18)
FRAC_DESC_FN(fc_sum_pol_at, 19)
FRAC_DESC_FN(fc_sum_pol_mw, 20)
FRAC_DESC_FN(fc_bond_at, 21)
FRAC_DESC_FN(fc_bond_n, 22)
FRAC_DESC_FN(fc_bond_o, 23)
FRAC_DESC_FN(fc_bond_c, 24)
FRAC_DESC_FN(fc_bond_s, 25)
FRAC_DESC_FN(fc_bond_p, 26)
FRAC_DESC_FN(fc_hb_donors, 27)
FRAC_DESC_FN(fc_hb_acceptors, 28)
FRAC_DESC_FN(fc_aromatic_atoms, 29)
FRAC_DESC_FN(fc_ring_atoms, 30)
FRAC_DESC_FN(fc_bridge_atoms, 31)
FRAC_DESC_FN(fc_charged_atoms, 32)
FRAC_DESC_FN(fc_small_r, 33)
FRAC_DESC_FN(fc_large_r, 34)
FRAC_DESC_FN(fc_low_polz, 35)
FRAC_DESC_FN(fc_high_polz, 36)
FRAC_DESC_FN(fc_high_ea, 37)
FRAC_DESC_FN(fc_low_ea, 38)
FRAC_DESC_FN(fc_small_vdw, 39)
FRAC_DESC_FN(fc_large_vdw, 40)
FRAC_DESC_FN(fc_en_mw, 41)
FRAC_DESC_FN(fc_en_bonded, 42)
FRAC_DESC_FN(fc_vdw_mw, 43)
FRAC_DESC_FN(fc_rcov_mw, 44)
FRAC_DESC_FN(fc_high_ox_state, 45)
FRAC_DESC_FN(fc_metalloid, 46)
FRAC_DESC_FN(fc_het_pol, 47)
FRAC_DESC_FN(fc_hal_pol, 48)
FRAC_DESC_FN(fc_heavy_pol, 49)
FRAC_DESC_FN(fc_sp3_heavy, 50)
FRAC_DESC_FN(fc_sp2_en_above, 51)
FRAC_DESC_FN(fc_ie_odd, 52)
FRAC_DESC_FN(fc_even_valence, 53)
FRAC_DESC_FN(fc_nonhetero_high_ea, 54)
FRAC_DESC_FN(fc_heavy_low_ie, 55)
FRAC_DESC_FN(fc_hetero_low_polz, 56)
FRAC_DESC_FN(fc_ox_en_above, 57)
FRAC_DESC_FN(fc_formal_nonzero, 58)
FRAC_DESC_FN(fc_formal_pos, 59)
FRAC_DESC_FN(fc_formal_neg, 60)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REG_FRAC(dname, ddesc, fn) do { \
    memset(&def, 0, sizeof(def)); \
    strncpy(def.name, dname, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, ddesc, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = fn; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_fractional(void) {
    descriptor_def_t def;

    /* MW fractions */
    REG_FRAC("FcC", "Carbon MW fraction", desc_fc_c);
    REG_FRAC("FcN", "Nitrogen MW fraction", desc_fc_n);
    REG_FRAC("FcO", "Oxygen MW fraction", desc_fc_o);
    REG_FRAC("FcS", "Sulfur MW fraction", desc_fc_s);
    REG_FRAC("FcF", "Fluorine MW fraction", desc_fc_f);
    REG_FRAC("FcCl", "Chlorine MW fraction", desc_fc_cl);
    REG_FRAC("FcBr", "Bromine MW fraction", desc_fc_br);
    REG_FRAC("FcI", "Iodine MW fraction", desc_fc_i);
    REG_FRAC("FcHalo", "Combined halogen MW fraction", desc_fc_halo);
    REG_FRAC("FcHetero", "Heteroatom MW fraction", desc_fc_hetero);

    /* EN-based fractions */
    REG_FRAC("FcPolar", "Fraction of atoms with EN > C", desc_fc_polar);
    REG_FRAC("FcApolar", "Fraction of atoms with EN <= C", desc_fc_apolar);
    REG_FRAC("FcENAboveAvg", "Atoms with EN above molecule avg", desc_fc_en_above_avg);
    REG_FRAC("FcENBelowAvg", "Atoms with EN below molecule avg", desc_fc_en_below_avg);
    REG_FRAC("FcENHigh", "Fraction with EN > 3.5", desc_fc_en_high);

    /* Bond fractions */
    REG_FRAC("FcCSp3", "Fraction of C-C sp3 bonds", desc_fc_csp3);
    REG_FRAC("FcCSp2", "Fraction of C-C sp2 bonds", desc_fc_csp2);
    REG_FRAC("FcUnpol", "Fraction of unpolar bonds", desc_fc_unpol);
    REG_FRAC("FcPol", "Fraction of polar bonds", desc_fc_pol);

    /* EN averages */
    REG_FRAC("FcSumPolAt", "Mean electronegativity per atom", desc_fc_sum_pol_at);
    REG_FRAC("FcSumPolMW", "Electronegativity per MW", desc_fc_sum_pol_mw);

    /* Bond density */
    REG_FRAC("FcBondAt", "Total bond order per atom", desc_fc_bond_at);

    /* Element bond fractions */
    REG_FRAC("FcBondN", "Non-single bonds on N", desc_fc_bond_n);
    REG_FRAC("FcBondO", "Non-single bonds on O", desc_fc_bond_o);
    REG_FRAC("FcBondC", "Non-single bonds on C", desc_fc_bond_c);
    REG_FRAC("FcBondS", "Non-single bonds on S", desc_fc_bond_s);
    REG_FRAC("FcBondP", "Non-single bonds on P", desc_fc_bond_p);

    /* Structural fractions */
    REG_FRAC("FcHBDonors", "Fraction H-bond donors", desc_fc_hb_donors);
    REG_FRAC("FcHBAcceptors", "Fraction H-bond acceptors", desc_fc_hb_acceptors);
    REG_FRAC("FcAromaticAtoms", "Fraction aromatic atoms", desc_fc_aromatic_atoms);
    REG_FRAC("FcRingAtoms", "Fraction ring atoms", desc_fc_ring_atoms);
    REG_FRAC("FcBridgeAtoms", "Fraction bridge atoms", desc_fc_bridge_atoms);
    REG_FRAC("FcChargedAtoms", "Fraction charged atoms", desc_fc_charged_atoms);

    /* Size fractions */
    REG_FRAC("FcSmallR", "Fraction covalent radius < 0.77", desc_fc_small_r);
    REG_FRAC("FcLargeR", "Fraction covalent radius > 1.1", desc_fc_large_r);
    REG_FRAC("FcLowPolz", "Fraction polarizability < 1.0", desc_fc_low_polz);
    REG_FRAC("FcHighPolz", "Fraction polarizability > 3.5", desc_fc_high_polz);
    REG_FRAC("FcHighEA", "Fraction electron affinity > 2", desc_fc_high_ea);
    REG_FRAC("FcLowEA", "Fraction electron affinity < 0.5", desc_fc_low_ea);
    REG_FRAC("FcSmallVdW", "Fraction VdW volume < 15", desc_fc_small_vdw);
    REG_FRAC("FcLargeVdW", "Fraction VdW volume > 25", desc_fc_large_vdw);

    /* MW-weighted */
    REG_FRAC("FcENMW", "EN-weighted MW fraction", desc_fc_en_mw);
    REG_FRAC("FcENBonded", "Bonds with both EN > 3.0", desc_fc_en_bonded);
    REG_FRAC("FcVdWMW", "VdW-weighted MW fraction", desc_fc_vdw_mw);
    REG_FRAC("FcRcovMW", "Covalent radius-weighted MW", desc_fc_rcov_mw);

    /* Oxidation and special */
    REG_FRAC("FcHighOxState", "Fraction max ox state >= 4", desc_fc_high_ox_state);
    REG_FRAC("FcMetalloid", "Fraction metalloid atoms", desc_fc_metalloid);
    REG_FRAC("FcHETpol", "Heteroatoms that are polar", desc_fc_het_pol);
    REG_FRAC("FcHALpol", "Halogens with polar bonds", desc_fc_hal_pol);
    REG_FRAC("FcHeavyPol", "Heavy atoms in polar bonds", desc_fc_heavy_pol);
    REG_FRAC("FcSp3HeavyAtoms", "sp3 heavy atoms fraction", desc_fc_sp3_heavy);
    REG_FRAC("FcSp2ENAboveAvg", "sp2 atoms with EN > avg", desc_fc_sp2_en_above);
    REG_FRAC("FcIEOdd", "Atoms with odd int(IE)", desc_fc_ie_odd);
    REG_FRAC("FcEvenValence", "Even valence atoms", desc_fc_even_valence);
    REG_FRAC("FcNonHeteroHighEA", "Non-hetero with EA > 1.0", desc_fc_nonhetero_high_ea);
    REG_FRAC("FcHeavyLowIE", "Heavy atoms with IE < 11", desc_fc_heavy_low_ie);
    REG_FRAC("FcHeteroLowPolz", "Heteroatoms polz < 1.5", desc_fc_hetero_low_polz);
    REG_FRAC("FcOxENAbove", "Ox state * EN > 10", desc_fc_ox_en_above);

    /* Formal charge fractions */
    REG_FRAC("FcFormalNonZero", "Non-zero formal charge", desc_fc_formal_nonzero);
    REG_FRAC("FcFormalPos", "Positive formal charge", desc_fc_formal_pos);
    REG_FRAC("FcFormalNeg", "Negative formal charge", desc_fc_formal_neg);
}
