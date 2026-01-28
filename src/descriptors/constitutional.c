/**
 * @file constitutional.c
 * @brief Constitutional Extension Descriptors - Extended atom and bond counts
 *
 * Categories:
 * 1. Element Combination Counts (12 descriptors)
 *    - Combined element groups (CNO, CNS, CNOS, etc.)
 *    - Period-based counts, electronegativity-based counts
 *
 * 2. Bond Environment Counts (12 descriptors)
 *    - Polar/nonpolar bonds, heterobonds, carbon-only bonds
 *    - Ring double bonds, conjugated bonds, bridging bonds
 *
 * 3. Hybrid Counts (10 descriptors)
 *    - sp3/sp2 carbons bonded to heteroatoms/halogens
 *    - N/O/S in various environments
 *
 * Total: 34 descriptors
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* Total constitutional extension descriptors */
#define NUM_CONST_DESCRIPTORS 34

/* Electronegativity values */
static const double CONST_EN[] = {
    [ELEM_H]  = 2.20,  [ELEM_C]  = 2.55,  [ELEM_N]  = 3.04,
    [ELEM_O]  = 3.44,  [ELEM_F]  = 3.98,  [ELEM_P]  = 2.19,
    [ELEM_S]  = 2.58,  [ELEM_Cl] = 3.16,  [ELEM_Br] = 2.96,
    [ELEM_I]  = 2.66,  [ELEM_Si] = 1.90,  [ELEM_B]  = 2.04,
    [ELEM_Se] = 2.55,  [ELEM_As] = 2.18,  [ELEM_Li] = 0.98,
    [ELEM_Na] = 0.93,  [ELEM_K]  = 0.82,  [ELEM_Mg] = 1.31,
    [ELEM_Ca] = 1.00,
};

/* Ionization potential (eV) */
static const double CONST_IP[] = {
    [ELEM_H]  = 13.60, [ELEM_C]  = 11.26, [ELEM_N]  = 14.53,
    [ELEM_O]  = 13.62, [ELEM_F]  = 17.42, [ELEM_P]  = 10.49,
    [ELEM_S]  = 10.36, [ELEM_Cl] = 12.97, [ELEM_Br] = 11.81,
    [ELEM_I]  = 10.45, [ELEM_Si] = 8.15,  [ELEM_B]  = 8.30,
    [ELEM_Se] = 9.75,  [ELEM_As] = 9.82,
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_en(element_t elem) {
    if (elem <= 0 || elem >= 128) return 2.55;
    double v = CONST_EN[elem];
    return (v > 0) ? v : 2.55;
}

static inline double get_ip(element_t elem) {
    if (elem <= 0 || elem >= 128) return 11.26;
    double v = CONST_IP[elem];
    return (v > 0) ? v : 11.26;
}

/* Check if element is a halogen */
static inline bool is_halogen(element_t elem) {
    return elem == ELEM_F || elem == ELEM_Cl || elem == ELEM_Br || elem == ELEM_I;
}

/* Check if element is a heteroatom (not C, H) */
static inline bool is_heteroatom(element_t elem) {
    return elem != ELEM_C && elem != ELEM_H;
}

/* Check if atom is in any ring */
static bool is_ring_atom(const molecule_t* mol, int atom_idx) {
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        for (int i = 0; i < ring->size; i++) {
            if (ring->atoms[i] == atom_idx) return true;
        }
    }
    return false;
}

/* Check if bond is in any ring */
static bool is_ring_bond(const molecule_t* mol, int atom1, int atom2) {
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        bool found1 = false, found2 = false;
        for (int i = 0; i < ring->size; i++) {
            if (ring->atoms[i] == atom1) found1 = true;
            if (ring->atoms[i] == atom2) found2 = true;
        }
        if (found1 && found2) return true;
    }
    return false;
}

/* Get carbon hybridization */
static int get_carbon_hybridization(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    if (atom->element != ELEM_C) return 0;

    /* Count multiple bonds */
    int triple = 0, double_bonds = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int nb = atom->neighbors[i];
        /* Find bond type */
        for (int b = 0; b < mol->num_bonds; b++) {
            const bond_t* bond = &mol->bonds[b];
            if ((bond->atom1 == atom_idx && bond->atom2 == nb) ||
                (bond->atom2 == atom_idx && bond->atom1 == nb)) {
                if (bond->type == BOND_TRIPLE) triple++;
                else if (bond->type == BOND_DOUBLE) double_bonds++;
                break;
            }
        }
    }

    if (triple > 0) return 1;  /* sp */
    if (double_bonds > 0) return 2;  /* sp2 */
    return 3;  /* sp3 */
}

/* Check if atom is aromatic */
static bool is_aromatic_atom(const molecule_t* mol, int atom_idx) {
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if ((bond->atom1 == atom_idx || bond->atom2 == atom_idx) && bond->aromatic) {
            return true;
        }
    }
    return false;
}

/* ============================================================================
 * Main Computation
 * ============================================================================ */

int descriptors_compute_constitutional_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_CONST_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    /* ====================================================================
     * Section 1: Element Combination Counts (12 descriptors, idx 0-11)
     * ==================================================================== */

    int cnt_cno = 0, cnt_cns = 0, cnt_cnos = 0;
    int cnt_heterohal = 0, cnt_polar_heavy = 0, cnt_nonpolar_heavy = 0;
    int cnt_period2 = 0, cnt_period3 = 0;
    int cnt_high_en = 0, cnt_low_en = 0;
    int cnt_high_ip = 0, cnt_low_ip = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_H) continue;

        /* CNO, CNS, CNOS */
        if (elem == ELEM_C || elem == ELEM_N || elem == ELEM_O) cnt_cno++;
        if (elem == ELEM_C || elem == ELEM_N || elem == ELEM_S) cnt_cns++;
        if (elem == ELEM_C || elem == ELEM_N || elem == ELEM_O || elem == ELEM_S) cnt_cnos++;

        /* Hetero + halogens */
        if (is_heteroatom(elem) || is_halogen(elem)) cnt_heterohal++;

        /* Polar heavy (N, O, S, P) */
        if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S || elem == ELEM_P) {
            cnt_polar_heavy++;
        }

        /* Non-polar heavy (C, Si, halogens) */
        if (elem == ELEM_C || elem == ELEM_Si || is_halogen(elem)) {
            cnt_nonpolar_heavy++;
        }

        /* Period 2 (C, N, O, F) */
        if (elem == ELEM_C || elem == ELEM_N || elem == ELEM_O || elem == ELEM_F) {
            cnt_period2++;
        }

        /* Period 3 (Si, P, S, Cl) */
        if (elem == ELEM_Si || elem == ELEM_P || elem == ELEM_S || elem == ELEM_Cl) {
            cnt_period3++;
        }

        /* EN thresholds */
        double en = get_en(elem);
        if (en > 3.0) cnt_high_en++;
        if (en < 2.5) cnt_low_en++;

        /* IP thresholds */
        double ip = get_ip(elem);
        if (ip > 12.0) cnt_high_ip++;
        if (ip < 10.0) cnt_low_ip++;
    }

    int idx = 0;
    values[idx++].d = cnt_cno;
    values[idx++].d = cnt_cns;
    values[idx++].d = cnt_cnos;
    values[idx++].d = cnt_heterohal;
    values[idx++].d = cnt_polar_heavy;
    values[idx++].d = cnt_nonpolar_heavy;
    values[idx++].d = cnt_period2;
    values[idx++].d = cnt_period3;
    values[idx++].d = cnt_high_en;
    values[idx++].d = cnt_low_en;
    values[idx++].d = cnt_high_ip;
    values[idx++].d = cnt_low_ip;

    /* ====================================================================
     * Section 2: Bond Environment Counts (12 descriptors, idx 12-23)
     * ==================================================================== */

    int polar_bonds = 0, nonpolar_bonds = 0;
    int hetero_bonds = 0, carbon_only_bonds = 0;
    int ring_double_bonds = 0, chain_double_bonds = 0;
    int conjugated_bonds = 0, rot_in_ring = 0;
    int multi_ring_bonds = 0, bridging_bonds = 0;
    int terminal_bonds = 0, branch_bonds = 0;

    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        element_t e1 = mol->atoms[bond->atom1].element;
        element_t e2 = mol->atoms[bond->atom2].element;

        if (e1 == ELEM_H || e2 == ELEM_H) continue;

        /* Polar/nonpolar based on EN difference */
        double en_diff = fabs(get_en(e1) - get_en(e2));
        if (en_diff > 0.5) polar_bonds++;
        else nonpolar_bonds++;

        /* Hetero bonds */
        if (is_heteroatom(e1) || is_heteroatom(e2)) {
            hetero_bonds++;
        }

        /* Carbon-only bonds */
        if (e1 == ELEM_C && e2 == ELEM_C) {
            carbon_only_bonds++;
        }

        /* Ring vs chain double bonds */
        if (bond->type == BOND_DOUBLE) {
            if (is_ring_bond(mol, bond->atom1, bond->atom2)) {
                ring_double_bonds++;
            } else {
                chain_double_bonds++;
            }
        }

        /* Conjugated bonds (aromatic or double adjacent to double) */
        if (bond->aromatic) {
            conjugated_bonds++;
        }

        /* Rotatable bonds in ring (single, not aromatic, in ring) */
        if (bond->type == BOND_SINGLE && !bond->aromatic &&
            is_ring_bond(mol, bond->atom1, bond->atom2)) {
            rot_in_ring++;
        }

        /* Multi-ring bonds (in more than one ring) */
        int ring_count = 0;
        for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
            const ring_t* ring = &mol->rings[r];
            bool found1 = false, found2 = false;
            for (int i = 0; i < ring->size; i++) {
                if (ring->atoms[i] == bond->atom1) found1 = true;
                if (ring->atoms[i] == bond->atom2) found2 = true;
            }
            if (found1 && found2) ring_count++;
        }
        if (ring_count > 1) multi_ring_bonds++;

        /* Bridging bonds (connects ring systems) */
        bool a1_ring = is_ring_atom(mol, bond->atom1);
        bool a2_ring = is_ring_atom(mol, bond->atom2);
        if (a1_ring && a2_ring && !is_ring_bond(mol, bond->atom1, bond->atom2)) {
            bridging_bonds++;
        }

        /* Terminal bonds (to degree-1 atom) */
        int deg1 = 0, deg2 = 0;
        for (int n = 0; n < mol->atoms[bond->atom1].num_neighbors; n++) {
            if (mol->atoms[mol->atoms[bond->atom1].neighbors[n]].element != ELEM_H) deg1++;
        }
        for (int n = 0; n < mol->atoms[bond->atom2].num_neighbors; n++) {
            if (mol->atoms[mol->atoms[bond->atom2].neighbors[n]].element != ELEM_H) deg2++;
        }
        if (deg1 == 1 || deg2 == 1) terminal_bonds++;
        if (deg1 >= 3 || deg2 >= 3) branch_bonds++;
    }

    values[idx++].d = polar_bonds;
    values[idx++].d = nonpolar_bonds;
    values[idx++].d = hetero_bonds;
    values[idx++].d = carbon_only_bonds;
    values[idx++].d = ring_double_bonds;
    values[idx++].d = chain_double_bonds;
    values[idx++].d = conjugated_bonds;
    values[idx++].d = rot_in_ring;
    values[idx++].d = multi_ring_bonds;
    values[idx++].d = bridging_bonds;
    values[idx++].d = terminal_bonds;
    values[idx++].d = branch_bonds;

    /* ====================================================================
     * Section 3: Hybrid Counts (10 descriptors, idx 24-33)
     * ==================================================================== */

    int csp3_hetero = 0, csp2_hetero = 0;
    int csp3_hal = 0, csp2_hal = 0;
    int n_aromatic = 0, n_aliphatic = 0;
    int o_carbonyl = 0, o_ether = 0;
    int s_thioether = 0, s_oxidized = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        const atom_t* atom = &mol->atoms[i];

        if (elem == ELEM_C) {
            int hyb = get_carbon_hybridization(mol, i);
            bool has_hetero = false, has_hal = false;

            for (int n = 0; n < atom->num_neighbors; n++) {
                element_t nb_elem = mol->atoms[atom->neighbors[n]].element;
                if (is_heteroatom(nb_elem)) has_hetero = true;
                if (is_halogen(nb_elem)) has_hal = true;
            }

            if (hyb == 3 && has_hetero) csp3_hetero++;
            if (hyb == 2 && has_hetero) csp2_hetero++;
            if (hyb == 3 && has_hal) csp3_hal++;
            if (hyb == 2 && has_hal) csp2_hal++;
        }
        else if (elem == ELEM_N) {
            if (is_aromatic_atom(mol, i)) n_aromatic++;
            else n_aliphatic++;
        }
        else if (elem == ELEM_O) {
            /* Check for carbonyl vs ether */
            bool is_carbonyl = false;
            for (int n = 0; n < atom->num_neighbors; n++) {
                int nb = atom->neighbors[n];
                if (mol->atoms[nb].element == ELEM_C) {
                    /* Check bond type */
                    for (int b = 0; b < mol->num_bonds; b++) {
                        const bond_t* bond = &mol->bonds[b];
                        if ((bond->atom1 == i && bond->atom2 == nb) ||
                            (bond->atom2 == i && bond->atom1 == nb)) {
                            if (bond->type == BOND_DOUBLE) is_carbonyl = true;
                            break;
                        }
                    }
                }
            }
            if (is_carbonyl) o_carbonyl++;
            else if (atom->num_neighbors == 2) o_ether++;
        }
        else if (elem == ELEM_S) {
            /* Check for thioether vs oxidized */
            int o_neighbors = 0;
            for (int n = 0; n < atom->num_neighbors; n++) {
                if (mol->atoms[atom->neighbors[n]].element == ELEM_O) o_neighbors++;
            }
            if (o_neighbors >= 1) s_oxidized++;
            else if (atom->num_neighbors == 2) s_thioether++;
        }
    }

    values[idx++].d = csp3_hetero;
    values[idx++].d = csp2_hetero;
    values[idx++].d = csp3_hal;
    values[idx++].d = csp2_hal;
    values[idx++].d = n_aromatic;
    values[idx++].d = n_aliphatic;
    values[idx++].d = o_carbonyl;
    values[idx++].d = o_ether;
    values[idx++].d = s_thioether;
    values[idx++].d = s_oxidized;

    return NUM_CONST_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* const_cached_mol = NULL;
static _Thread_local uint64_t const_cached_gen = 0;
static _Thread_local descriptor_value_t const_cached_values[NUM_CONST_DESCRIPTORS];

static inline void ensure_constitutional_computed(const molecule_t* mol) {
    uint64_t current_gen = descriptor_cache_generation();
    if (const_cached_mol != mol || const_cached_gen != current_gen) {
        descriptors_compute_constitutional_all(mol, const_cached_values);
        const_cached_mol = mol;
        const_cached_gen = current_gen;
    }
}

#define DEFINE_CONST_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_constitutional_computed(mol); \
    value->d = const_cached_values[idx].d; \
    return CCHEM_OK; \
}

/* Element Combination Counts (12) */
DEFINE_CONST_FUNC(const_cno, 0)
DEFINE_CONST_FUNC(const_cns, 1)
DEFINE_CONST_FUNC(const_cnos, 2)
DEFINE_CONST_FUNC(const_heterohal, 3)
DEFINE_CONST_FUNC(const_polarheavy, 4)
DEFINE_CONST_FUNC(const_nonpolarheavy, 5)
DEFINE_CONST_FUNC(const_period2, 6)
DEFINE_CONST_FUNC(const_period3, 7)
DEFINE_CONST_FUNC(const_highen, 8)
DEFINE_CONST_FUNC(const_lowen, 9)
DEFINE_CONST_FUNC(const_highip, 10)
DEFINE_CONST_FUNC(const_lowip, 11)

/* Bond Environment Counts (12) */
DEFINE_CONST_FUNC(const_polarbonds, 12)
DEFINE_CONST_FUNC(const_nonpolarbonds, 13)
DEFINE_CONST_FUNC(const_heterobonds, 14)
DEFINE_CONST_FUNC(const_carbonbonds, 15)
DEFINE_CONST_FUNC(const_ringdouble, 16)
DEFINE_CONST_FUNC(const_chaindouble, 17)
DEFINE_CONST_FUNC(const_conjugated, 18)
DEFINE_CONST_FUNC(const_rotinring, 19)
DEFINE_CONST_FUNC(const_multiringbonds, 20)
DEFINE_CONST_FUNC(const_bridging, 21)
DEFINE_CONST_FUNC(const_terminalbonds, 22)
DEFINE_CONST_FUNC(const_branchbonds, 23)

/* Hybrid Counts (10) */
DEFINE_CONST_FUNC(const_csp3hetero, 24)
DEFINE_CONST_FUNC(const_csp2hetero, 25)
DEFINE_CONST_FUNC(const_csp3hal, 26)
DEFINE_CONST_FUNC(const_csp2hal, 27)
DEFINE_CONST_FUNC(const_naromatic, 28)
DEFINE_CONST_FUNC(const_naliphatic, 29)
DEFINE_CONST_FUNC(const_ocarbonyl, 30)
DEFINE_CONST_FUNC(const_oether, 31)
DEFINE_CONST_FUNC(const_sthioether, 32)
DEFINE_CONST_FUNC(const_soxidized, 33)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_CONST(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_CONSTITUTIONAL; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_constitutional(void) {
    /* Element Combination Counts */
    REGISTER_CONST("Const_CNO", "Atoms that are C, N, or O", desc_const_cno);
    REGISTER_CONST("Const_CNS", "Atoms that are C, N, or S", desc_const_cns);
    REGISTER_CONST("Const_CNOS", "Atoms that are C, N, O, or S", desc_const_cnos);
    REGISTER_CONST("Const_HeteroHal", "Heteroatoms + halogens", desc_const_heterohal);
    REGISTER_CONST("Const_PolarHeavy", "Polar heavy atoms (N, O, S, P)", desc_const_polarheavy);
    REGISTER_CONST("Const_NonPolarHeavy", "Non-polar heavy atoms", desc_const_nonpolarheavy);
    REGISTER_CONST("Const_Period2", "Period 2 atoms (C, N, O, F)", desc_const_period2);
    REGISTER_CONST("Const_Period3", "Period 3 atoms (Si, P, S, Cl)", desc_const_period3);
    REGISTER_CONST("Const_HighEN", "Atoms with EN > 3.0", desc_const_highen);
    REGISTER_CONST("Const_LowEN", "Atoms with EN < 2.5", desc_const_lowen);
    REGISTER_CONST("Const_HighIP", "Atoms with IP > 12 eV", desc_const_highip);
    REGISTER_CONST("Const_LowIP", "Atoms with IP < 10 eV", desc_const_lowip);

    /* Bond Environment Counts */
    REGISTER_CONST("Const_PolarBonds", "Bonds with |EN diff| > 0.5", desc_const_polarbonds);
    REGISTER_CONST("Const_NonpolarBonds", "Bonds with |EN diff| < 0.5", desc_const_nonpolarbonds);
    REGISTER_CONST("Const_HeteroBonds", "Bonds involving heteroatoms", desc_const_heterobonds);
    REGISTER_CONST("Const_CarbonOnlyBonds", "C-C bonds", desc_const_carbonbonds);
    REGISTER_CONST("Const_RingDoubleBonds", "Double bonds in rings", desc_const_ringdouble);
    REGISTER_CONST("Const_ChainDoubleBonds", "Double bonds in chains", desc_const_chaindouble);
    REGISTER_CONST("Const_ConjugatedBonds", "Bonds in conjugated systems", desc_const_conjugated);
    REGISTER_CONST("Const_RotInRing", "Would-be-rotatable in ring", desc_const_rotinring);
    REGISTER_CONST("Const_MultiRingBonds", "Bonds in multiple rings", desc_const_multiringbonds);
    REGISTER_CONST("Const_BridgingBonds", "Bonds connecting ring systems", desc_const_bridging);
    REGISTER_CONST("Const_TerminalBonds", "Bonds to terminal atoms", desc_const_terminalbonds);
    REGISTER_CONST("Const_BranchBonds", "Bonds at branch points", desc_const_branchbonds);

    /* Hybrid Counts */
    REGISTER_CONST("Const_Csp3_Hetero", "sp3 C bonded to heteroatom", desc_const_csp3hetero);
    REGISTER_CONST("Const_Csp2_Hetero", "sp2 C bonded to heteroatom", desc_const_csp2hetero);
    REGISTER_CONST("Const_Csp3_Hal", "sp3 C bonded to halogen", desc_const_csp3hal);
    REGISTER_CONST("Const_Csp2_Hal", "sp2 C bonded to halogen", desc_const_csp2hal);
    REGISTER_CONST("Const_N_Aromatic", "N in aromatic environment", desc_const_naromatic);
    REGISTER_CONST("Const_N_Aliphatic", "N in aliphatic environment", desc_const_naliphatic);
    REGISTER_CONST("Const_O_Carbonyl", "O in carbonyl", desc_const_ocarbonyl);
    REGISTER_CONST("Const_O_Ether", "O in ether", desc_const_oether);
    REGISTER_CONST("Const_S_Thioether", "S in thioether", desc_const_sthioether);
    REGISTER_CONST("Const_S_Oxidized", "S in sulfoxide/sulfone", desc_const_soxidized);
}
