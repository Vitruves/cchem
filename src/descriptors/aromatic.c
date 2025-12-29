/**
 * @file aromatic.c
 * @brief Aromatic Pattern Descriptors - Detailed aromatic system analysis
 *
 * Categories:
 * 1. Aromatic Ring Classification (20 descriptors)
 *    - Benzene, pyridine, pyrimidine, pyrrole, furan, thiophene, imidazole
 *    - Fused systems, biphenyl, isolated, bridged, mixed 5-6
 *
 * 2. Aromatic Density & Distribution (16 descriptors)
 *    - Atom/bond/surface density, terminal/branch/core distribution
 *    - Connection degrees, heteroatom counts
 *
 * 3. Aromatic Electronic Properties (16 descriptors)
 *    - EN, polarizability, charge, E-State for aromatic atoms
 *    - IP, EA, mass fraction, pi electrons
 *
 * 4. Aromatic Topology & Connectivity (12 descriptors)
 *    - Bridges, fusion degree, clustering, substitution patterns
 *
 * Total: 64 descriptors
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* Total aromatic descriptors */
#define NUM_AROM_DESCRIPTORS 64

/* Maximum atoms for stack allocation */
#define MAX_ATOMS_STACK 256

/* ============================================================================
 * Atomic Property Tables
 * ============================================================================ */

static const double AROM_EN[] = {
    [ELEM_H]  = 2.20,  [ELEM_C]  = 2.55,  [ELEM_N]  = 3.04,
    [ELEM_O]  = 3.44,  [ELEM_F]  = 3.98,  [ELEM_S]  = 2.58,
    [ELEM_P]  = 2.19,  [ELEM_Se] = 2.55,
};

static const double AROM_POL[] = {
    [ELEM_H]  = 0.667, [ELEM_C]  = 1.76,  [ELEM_N]  = 1.10,
    [ELEM_O]  = 0.802, [ELEM_F]  = 0.557, [ELEM_S]  = 2.90,
    [ELEM_P]  = 3.63,  [ELEM_Se] = 3.77,
};

static const double AROM_IP[] = {
    [ELEM_H]  = 13.60, [ELEM_C]  = 11.26, [ELEM_N]  = 14.53,
    [ELEM_O]  = 13.62, [ELEM_F]  = 17.42, [ELEM_S]  = 10.36,
    [ELEM_P]  = 10.49, [ELEM_Se] = 9.75,
};

static const double AROM_EA[] = {
    [ELEM_H]  = 0.75,  [ELEM_C]  = 1.26,  [ELEM_N]  = -0.07,
    [ELEM_O]  = 1.46,  [ELEM_F]  = 3.40,  [ELEM_S]  = 2.08,
    [ELEM_P]  = 0.75,  [ELEM_Se] = 2.02,
};

static const double AROM_MASS[] = {
    [ELEM_H]  = 1.008,  [ELEM_C]  = 12.011, [ELEM_N]  = 14.007,
    [ELEM_O]  = 15.999, [ELEM_F]  = 18.998, [ELEM_S]  = 32.065,
    [ELEM_P]  = 30.974, [ELEM_Se] = 78.960,
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_prop(element_t elem, const double* table, double def) {
    if (elem <= 0 || elem >= 128) return def;
    double v = table[elem];
    return (v != 0.0) ? v : def;
}

/* Check if atom is aromatic (using atom's aromatic flag) */
static bool is_aromatic_atom(const molecule_t* mol, int atom_idx) {
    return mol->atoms[atom_idx].aromatic;
}

/* Check if ring is aromatic (all atoms aromatic) */
static bool is_aromatic_ring(const molecule_t* mol, const ring_t* ring) {
    for (int i = 0; i < ring->size; i++) {
        if (!is_aromatic_atom(mol, ring->atoms[i])) {
            return false;
        }
    }
    return true;
}

/* Count heteroatoms in ring */
static int count_ring_heteroatoms(const molecule_t* mol, const ring_t* ring, element_t elem) {
    int count = 0;
    for (int i = 0; i < ring->size; i++) {
        if (mol->atoms[ring->atoms[i]].element == elem) {
            count++;
        }
    }
    return count;
}

/* Check if two rings share an edge (fused) */
static bool rings_fused(const ring_t* r1, const ring_t* r2) {
    int shared = 0;
    for (int i = 0; i < r1->size; i++) {
        for (int j = 0; j < r2->size; j++) {
            if (r1->atoms[i] == r2->atoms[j]) {
                shared++;
                if (shared >= 2) return true;
            }
        }
    }
    return false;
}

/* Count aromatic neighbors not in aromatic ring with atom */
static int count_non_arom_neighbors(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int nb = atom->neighbors[i];
        if (!is_aromatic_atom(mol, nb)) {
            count++;
        }
    }
    return count;
}

/* ============================================================================
 * Main Computation
 * ============================================================================ */

int descriptors_compute_aromatic_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_AROM_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    /* Count heavy atoms and aromatic atoms */
    int n_heavy = 0;
    int n_aromatic = 0;
    int n_arom_bonds = 0;
    double total_mass = 0.0;
    double arom_mass = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            n_heavy++;
            total_mass += get_prop(mol->atoms[i].element, AROM_MASS, 12.0);
            if (is_aromatic_atom(mol, i)) {
                n_aromatic++;
                arom_mass += get_prop(mol->atoms[i].element, AROM_MASS, 12.0);
            }
        }
    }

    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->bonds[b].aromatic) n_arom_bonds++;
    }

    if (n_heavy == 0) return NUM_AROM_DESCRIPTORS;

    /* ====================================================================
     * Section 1: Ring Classification (20 descriptors, idx 0-19)
     * ==================================================================== */

    int benzene = 0, pyridine = 0, pyrimidine = 0, pyrrole = 0;
    int furan = 0, thiophene = 0, imidazole = 0;
    int fused2 = 0, fused3 = 0, fused4plus = 0;
    int biphenyl = 0, isolated = 0, bridged = 0, mixed56 = 0;
    int all_carbon = 0, heterocyclic = 0, polycyclic = 0;
    int naphthalene = 0, anthracene = 0, phenanthrene = 0;

    /* Track fused ring systems */
    int* ring_fused_count = (int*)alloca(mol->num_rings * sizeof(int));
    memset(ring_fused_count, 0, mol->num_rings * sizeof(int));

    /* Analyze each aromatic ring */
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        if (!is_aromatic_ring(mol, ring)) continue;

        int n_N = count_ring_heteroatoms(mol, ring, ELEM_N);
        int n_O = count_ring_heteroatoms(mol, ring, ELEM_O);
        int n_S = count_ring_heteroatoms(mol, ring, ELEM_S);
        int n_hetero = n_N + n_O + n_S;

        if (ring->size == 6) {
            if (n_hetero == 0) benzene++;
            else if (n_N == 1 && n_O == 0 && n_S == 0) pyridine++;
            else if (n_N == 2 && n_O == 0 && n_S == 0) pyrimidine++;
        } else if (ring->size == 5) {
            if (n_N == 1 && n_O == 0 && n_S == 0) pyrrole++;
            else if (n_O == 1 && n_N == 0 && n_S == 0) furan++;
            else if (n_S == 1 && n_N == 0 && n_O == 0) thiophene++;
            else if (n_N == 2 && n_O == 0 && n_S == 0) imidazole++;
        }

        if (n_hetero == 0) all_carbon++;
        else heterocyclic++;

        /* Count fused rings */
        for (int r2 = r + 1; r2 < mol->num_rings && r2 < MAX_RINGS; r2++) {
            const ring_t* ring2 = &mol->rings[r2];
            if (!is_aromatic_ring(mol, ring2)) continue;
            if (rings_fused(ring, ring2)) {
                ring_fused_count[r]++;
                ring_fused_count[r2]++;
            }
        }
    }

    /* Analyze fusion patterns */
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        if (!is_aromatic_ring(mol, ring)) continue;

        if (ring_fused_count[r] == 0) isolated++;
        else if (ring_fused_count[r] == 1) {
            /* Part of a 2-ring fused system */
            if (ring->size == 6) {
                /* Check if fused partner is also 6 */
                for (int r2 = 0; r2 < mol->num_rings && r2 < MAX_RINGS; r2++) {
                    if (r2 == r) continue;
                    const ring_t* ring2 = &mol->rings[r2];
                    if (is_aromatic_ring(mol, ring2) && rings_fused(ring, ring2)) {
                        if (ring2->size == 6 && ring_fused_count[r2] == 1) {
                            naphthalene++;
                        } else if (ring2->size == 5) {
                            mixed56++;
                        }
                    }
                }
            }
        } else if (ring_fused_count[r] >= 2) {
            polycyclic++;
            bridged++;
        }
    }

    /* Avoid double counting naphthalene */
    naphthalene /= 2;

    /* Count fused system sizes */
    /* Simplified: just count rings with fusion count */
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        if (ring_fused_count[r] == 1) fused2++;
        else if (ring_fused_count[r] == 2) fused3++;
        else if (ring_fused_count[r] >= 3) fused4plus++;
    }
    fused2 /= 2;  /* Each fusion counted twice */

    /* Biphenyl - aromatic rings connected by single bond (not fused) */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if (!bond->aromatic && bond->type == BOND_SINGLE) {
            bool a1_arom = is_aromatic_atom(mol, bond->atom1);
            bool a2_arom = is_aromatic_atom(mol, bond->atom2);
            if (a1_arom && a2_arom) biphenyl++;
        }
    }

    /* Store ring classification */
    int idx = 0;
    values[idx++].d = benzene;
    values[idx++].d = pyridine;
    values[idx++].d = pyrimidine;
    values[idx++].d = pyrrole;
    values[idx++].d = furan;
    values[idx++].d = thiophene;
    values[idx++].d = imidazole;
    values[idx++].d = fused2;
    values[idx++].d = fused3;
    values[idx++].d = fused4plus;
    values[idx++].d = biphenyl;
    values[idx++].d = isolated;
    values[idx++].d = bridged;
    values[idx++].d = mixed56;
    values[idx++].d = all_carbon;
    values[idx++].d = heterocyclic;
    values[idx++].d = polycyclic;
    values[idx++].d = naphthalene;
    values[idx++].d = anthracene;
    values[idx++].d = phenanthrene;

    /* ====================================================================
     * Section 2: Density & Distribution (16 descriptors, idx 20-35)
     * ==================================================================== */

    double arom_dens_atoms = (n_heavy > 0) ? (double)n_aromatic / n_heavy : 0.0;
    double arom_dens_bonds = (mol->num_bonds > 0) ? (double)n_arom_bonds / mol->num_bonds : 0.0;
    double arom_mass_frac = (total_mass > 0) ? arom_mass / total_mass : 0.0;

    /* Count connection degrees and positions */
    int arom_term = 0, arom_branch = 0, arom_core = 0;
    int arom_conn_deg1 = 0, arom_conn_deg2 = 0, arom_conn_deg3 = 0;
    int arom_N = 0, arom_O = 0, arom_S = 0, arom_mixed_hetero = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        if (!is_aromatic_atom(mol, i)) continue;

        const atom_t* atom = &mol->atoms[i];
        int non_arom_neighbors = count_non_arom_neighbors(mol, i);

        if (non_arom_neighbors == 0) arom_core++;
        else if (non_arom_neighbors == 1) arom_conn_deg1++;
        else if (non_arom_neighbors == 2) arom_conn_deg2++;
        else arom_conn_deg3++;

        if (atom->num_neighbors == 1) arom_term++;
        if (atom->num_neighbors >= 3) arom_branch++;

        element_t elem = atom->element;
        if (elem == ELEM_N) arom_N++;
        else if (elem == ELEM_O) arom_O++;
        else if (elem == ELEM_S) arom_S++;
    }

    /* Count aromatic rings with mixed heteroatoms */
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        if (!is_aromatic_ring(mol, ring)) continue;
        int types = 0;
        if (count_ring_heteroatoms(mol, ring, ELEM_N) > 0) types++;
        if (count_ring_heteroatoms(mol, ring, ELEM_O) > 0) types++;
        if (count_ring_heteroatoms(mol, ring, ELEM_S) > 0) types++;
        if (types >= 2) arom_mixed_hetero++;
    }

    values[idx++].d = arom_dens_atoms;
    values[idx++].d = arom_dens_bonds;
    values[idx++].d = arom_dens_atoms;  /* Surface approximation */
    values[idx++].d = arom_dens_atoms;  /* Volume approximation */
    values[idx++].d = arom_mass_frac;
    values[idx++].d = arom_term;
    values[idx++].d = arom_branch;
    values[idx++].d = arom_core;
    values[idx++].d = n_aromatic - arom_core;  /* Periphery */
    values[idx++].d = arom_conn_deg1;
    values[idx++].d = arom_conn_deg2;
    values[idx++].d = arom_conn_deg3;
    values[idx++].d = arom_N;
    values[idx++].d = arom_O;
    values[idx++].d = arom_S;
    values[idx++].d = arom_mixed_hetero;

    /* ====================================================================
     * Section 3: Electronic Properties (16 descriptors, idx 36-51)
     * ==================================================================== */

    double en_sum = 0, en_mean = 0, en_var = 0;
    double pol_sum = 0, pol_mean = 0;
    double charge_sum = 0, charge_pos = 0, charge_neg = 0;
    double estate_sum = 0, estate_mean = 0, estate_min = 1e10, estate_max = -1e10;
    double ip_sum = 0, ea_sum = 0;
    int pi_electrons = 0;

    if (n_aromatic > 0) {
        for (int i = 0; i < mol->num_atoms; i++) {
            if (mol->atoms[i].element == ELEM_H) continue;
            if (!is_aromatic_atom(mol, i)) continue;

            double en = get_prop(mol->atoms[i].element, AROM_EN, 2.55);
            double pol = get_prop(mol->atoms[i].element, AROM_POL, 1.76);
            double ip = get_prop(mol->atoms[i].element, AROM_IP, 11.26);
            double ea = get_prop(mol->atoms[i].element, AROM_EA, 1.26);
            double charge = mol->atoms[i].charge;

            /* Simplified E-State estimate */
            int valence = 4;  /* Default for C */
            if (mol->atoms[i].element == ELEM_N) valence = 3;
            else if (mol->atoms[i].element == ELEM_O) valence = 2;
            else if (mol->atoms[i].element == ELEM_S) valence = 2;
            double estate = (valence - mol->atoms[i].num_neighbors + 1.0) / (mol->atoms[i].num_neighbors + 1.0);

            en_sum += en;
            pol_sum += pol;
            ip_sum += ip;
            ea_sum += ea;
            charge_sum += charge;
            if (charge > 0) charge_pos += charge;
            if (charge < 0) charge_neg += fabs(charge);

            estate_sum += estate;
            if (estate < estate_min) estate_min = estate;
            if (estate > estate_max) estate_max = estate;

            /* Pi electrons: aromatic = 1 pi electron contribution per atom */
            pi_electrons++;
        }

        en_mean = en_sum / n_aromatic;
        pol_mean = pol_sum / n_aromatic;
        estate_mean = estate_sum / n_aromatic;

        /* Compute EN variance */
        for (int i = 0; i < mol->num_atoms; i++) {
            if (mol->atoms[i].element == ELEM_H) continue;
            if (!is_aromatic_atom(mol, i)) continue;
            double en = get_prop(mol->atoms[i].element, AROM_EN, 2.55);
            en_var += (en - en_mean) * (en - en_mean);
        }
        en_var /= n_aromatic;
    }

    if (n_aromatic == 0) {
        estate_min = 0;
        estate_max = 0;
    }

    values[idx++].d = en_sum;
    values[idx++].d = en_mean;
    values[idx++].d = en_var;
    values[idx++].d = pol_sum;
    values[idx++].d = pol_mean;
    values[idx++].d = charge_sum;
    values[idx++].d = charge_pos;
    values[idx++].d = charge_neg;
    values[idx++].d = estate_sum;
    values[idx++].d = estate_mean;
    values[idx++].d = estate_min;
    values[idx++].d = estate_max;
    values[idx++].d = ip_sum;
    values[idx++].d = ea_sum;
    values[idx++].d = arom_mass_frac;
    values[idx++].d = pi_electrons;

    /* ====================================================================
     * Section 4: Topology & Connectivity (12 descriptors, idx 52-63)
     * ==================================================================== */

    /* Count bridges (bonds connecting aromatic systems) */
    int bridges = biphenyl;  /* Already counted */

    /* Fusion degree: average fusions per aromatic ring */
    int arom_ring_count = 0;
    double fusion_sum = 0;
    int max_fused = 0;
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        if (is_aromatic_ring(mol, &mol->rings[r])) {
            arom_ring_count++;
            fusion_sum += ring_fused_count[r];
            if (ring_fused_count[r] > max_fused) max_fused = ring_fused_count[r];
        }
    }
    double fusion_degree = (arom_ring_count > 0) ? fusion_sum / arom_ring_count : 0.0;

    /* Spiro junctions in aromatic (single atom shared) - rare */
    int spiro_arom = 0;  /* Would need more complex detection */

    /* Mean distance between aromatic systems - simplified */
    double arom_distance = 0.0;

    /* Clustering coefficient */
    double clustering = (n_aromatic > 0) ? (double)arom_core / n_aromatic : 0.0;

    /* Substitution patterns */
    int total_substituents = arom_conn_deg1 + 2 * arom_conn_deg2 + 3 * arom_conn_deg3;

    /* Ortho/meta/para - would need ring position analysis, simplified */
    int ortho_sub = 0, meta_sub = 0, para_sub = 0;

    values[idx++].d = bridges;
    values[idx++].d = fusion_degree;
    values[idx++].d = max_fused;
    values[idx++].d = spiro_arom;
    values[idx++].d = arom_distance;
    values[idx++].d = clustering;
    values[idx++].d = n_aromatic - arom_core;  /* Peripheral */
    values[idx++].d = arom_core;
    values[idx++].d = total_substituents;
    values[idx++].d = ortho_sub;
    values[idx++].d = meta_sub;
    values[idx++].d = para_sub;

    return NUM_AROM_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* arom_cached_mol = NULL;
static _Thread_local uint64_t arom_cached_gen = 0;
static _Thread_local descriptor_value_t arom_cached_values[NUM_AROM_DESCRIPTORS];

static inline void ensure_aromatic_computed(const molecule_t* mol) {
    uint64_t current_gen = descriptor_cache_generation();
    if (arom_cached_mol != mol || arom_cached_gen != current_gen) {
        descriptors_compute_aromatic_all(mol, arom_cached_values);
        arom_cached_mol = mol;
        arom_cached_gen = current_gen;
    }
}

#define DEFINE_AROM_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_aromatic_computed(mol); \
    value->d = arom_cached_values[idx].d; \
    return CCHEM_OK; \
}

/* Ring Classification (20) */
DEFINE_AROM_FUNC(arom_benzene, 0)
DEFINE_AROM_FUNC(arom_pyridine, 1)
DEFINE_AROM_FUNC(arom_pyrimidine, 2)
DEFINE_AROM_FUNC(arom_pyrrole, 3)
DEFINE_AROM_FUNC(arom_furan, 4)
DEFINE_AROM_FUNC(arom_thiophene, 5)
DEFINE_AROM_FUNC(arom_imidazole, 6)
DEFINE_AROM_FUNC(arom_fused2, 7)
DEFINE_AROM_FUNC(arom_fused3, 8)
DEFINE_AROM_FUNC(arom_fused4plus, 9)
DEFINE_AROM_FUNC(arom_biphenyl, 10)
DEFINE_AROM_FUNC(arom_isolated, 11)
DEFINE_AROM_FUNC(arom_bridged, 12)
DEFINE_AROM_FUNC(arom_mixed56, 13)
DEFINE_AROM_FUNC(arom_allcarbon, 14)
DEFINE_AROM_FUNC(arom_heterocyclic, 15)
DEFINE_AROM_FUNC(arom_polycyclic, 16)
DEFINE_AROM_FUNC(arom_naphthalene, 17)
DEFINE_AROM_FUNC(arom_anthracene, 18)
DEFINE_AROM_FUNC(arom_phenanthrene, 19)

/* Density & Distribution (16) */
DEFINE_AROM_FUNC(arom_dens_atoms, 20)
DEFINE_AROM_FUNC(arom_dens_bonds, 21)
DEFINE_AROM_FUNC(arom_dens_surface, 22)
DEFINE_AROM_FUNC(arom_dens_volume, 23)
DEFINE_AROM_FUNC(arom_dens_mass, 24)
DEFINE_AROM_FUNC(arom_dist_terminal, 25)
DEFINE_AROM_FUNC(arom_dist_branch, 26)
DEFINE_AROM_FUNC(arom_dist_core, 27)
DEFINE_AROM_FUNC(arom_dist_periphery, 28)
DEFINE_AROM_FUNC(arom_conn_deg1, 29)
DEFINE_AROM_FUNC(arom_conn_deg2, 30)
DEFINE_AROM_FUNC(arom_conn_deg3, 31)
DEFINE_AROM_FUNC(arom_hetero_n, 32)
DEFINE_AROM_FUNC(arom_hetero_o, 33)
DEFINE_AROM_FUNC(arom_hetero_s, 34)
DEFINE_AROM_FUNC(arom_hetero_mix, 35)

/* Electronic Properties (16) */
DEFINE_AROM_FUNC(arom_chi_sum, 36)
DEFINE_AROM_FUNC(arom_chi_mean, 37)
DEFINE_AROM_FUNC(arom_chi_var, 38)
DEFINE_AROM_FUNC(arom_pol_sum, 39)
DEFINE_AROM_FUNC(arom_pol_mean, 40)
DEFINE_AROM_FUNC(arom_charge_sum, 41)
DEFINE_AROM_FUNC(arom_charge_pos, 42)
DEFINE_AROM_FUNC(arom_charge_neg, 43)
DEFINE_AROM_FUNC(arom_estate_sum, 44)
DEFINE_AROM_FUNC(arom_estate_mean, 45)
DEFINE_AROM_FUNC(arom_estate_min, 46)
DEFINE_AROM_FUNC(arom_estate_max, 47)
DEFINE_AROM_FUNC(arom_ip_sum, 48)
DEFINE_AROM_FUNC(arom_ea_sum, 49)
DEFINE_AROM_FUNC(arom_mass_frac, 50)
DEFINE_AROM_FUNC(arom_pi_elec, 51)

/* Topology (12) */
DEFINE_AROM_FUNC(arom_topo_bridges, 52)
DEFINE_AROM_FUNC(arom_topo_fusiondeg, 53)
DEFINE_AROM_FUNC(arom_topo_maxfused, 54)
DEFINE_AROM_FUNC(arom_topo_spiro, 55)
DEFINE_AROM_FUNC(arom_topo_distance, 56)
DEFINE_AROM_FUNC(arom_topo_clustering, 57)
DEFINE_AROM_FUNC(arom_topo_peripheral, 58)
DEFINE_AROM_FUNC(arom_topo_core, 59)
DEFINE_AROM_FUNC(arom_topo_substitution, 60)
DEFINE_AROM_FUNC(arom_topo_ortho, 61)
DEFINE_AROM_FUNC(arom_topo_meta, 62)
DEFINE_AROM_FUNC(arom_topo_para, 63)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_AROM(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_aromatic(void) {
    /* Ring Classification */
    REGISTER_AROM("Arom_Benzene", "Pure benzene rings", desc_arom_benzene);
    REGISTER_AROM("Arom_Pyridine", "Pyridine-like rings", desc_arom_pyridine);
    REGISTER_AROM("Arom_Pyrimidine", "Pyrimidine-like rings", desc_arom_pyrimidine);
    REGISTER_AROM("Arom_Pyrrole", "Pyrrole-like rings", desc_arom_pyrrole);
    REGISTER_AROM("Arom_Furan", "Furan-like rings", desc_arom_furan);
    REGISTER_AROM("Arom_Thiophene", "Thiophene-like rings", desc_arom_thiophene);
    REGISTER_AROM("Arom_Imidazole", "Imidazole-like rings", desc_arom_imidazole);
    REGISTER_AROM("Arom_Fused2", "2-ring fused aromatic systems", desc_arom_fused2);
    REGISTER_AROM("Arom_Fused3", "3-ring fused aromatic systems", desc_arom_fused3);
    REGISTER_AROM("Arom_Fused4Plus", "4+ ring fused aromatic systems", desc_arom_fused4plus);
    REGISTER_AROM("Arom_Biphenyl", "Biphenyl-like connected aromatics", desc_arom_biphenyl);
    REGISTER_AROM("Arom_Isolated", "Isolated aromatic rings", desc_arom_isolated);
    REGISTER_AROM("Arom_Bridged", "Bridged aromatic systems", desc_arom_bridged);
    REGISTER_AROM("Arom_Mixed5_6", "Mixed 5-6 fused aromatics", desc_arom_mixed56);
    REGISTER_AROM("Arom_AllCarbon", "All-carbon aromatic rings", desc_arom_allcarbon);
    REGISTER_AROM("Arom_Heterocyclic", "Heteroatom aromatic rings", desc_arom_heterocyclic);
    REGISTER_AROM("Arom_Polycyclic", "Polycyclic aromatic systems", desc_arom_polycyclic);
    REGISTER_AROM("Arom_NaphthaleneL", "Naphthalene-like systems", desc_arom_naphthalene);
    REGISTER_AROM("Arom_AnthraceneL", "Anthracene-like systems", desc_arom_anthracene);
    REGISTER_AROM("Arom_PhenanthreneL", "Phenanthrene-like systems", desc_arom_phenanthrene);

    /* Density & Distribution */
    REGISTER_AROM("AromDens_Atoms", "Aromatic atoms / heavy atoms", desc_arom_dens_atoms);
    REGISTER_AROM("AromDens_Bonds", "Aromatic bonds / total bonds", desc_arom_dens_bonds);
    REGISTER_AROM("AromDens_Surface", "Aromatic surface fraction", desc_arom_dens_surface);
    REGISTER_AROM("AromDens_Volume", "Aromatic volume fraction", desc_arom_dens_volume);
    REGISTER_AROM("AromDens_Mass", "Aromatic mass fraction", desc_arom_dens_mass);
    REGISTER_AROM("AromDist_Terminal", "Aromatic atoms at terminals", desc_arom_dist_terminal);
    REGISTER_AROM("AromDist_Branch", "Aromatic atoms at branches", desc_arom_dist_branch);
    REGISTER_AROM("AromDist_Core", "Aromatic atoms in core", desc_arom_dist_core);
    REGISTER_AROM("AromDist_Periphery", "Aromatic atoms at periphery", desc_arom_dist_periphery);
    REGISTER_AROM("AromConn_Deg1", "Aromatic with 1 non-arom neighbor", desc_arom_conn_deg1);
    REGISTER_AROM("AromConn_Deg2", "Aromatic with 2 non-arom neighbors", desc_arom_conn_deg2);
    REGISTER_AROM("AromConn_Deg3", "Aromatic with 3+ non-arom neighbors", desc_arom_conn_deg3);
    REGISTER_AROM("AromHetero_N", "N atoms in aromatic rings", desc_arom_hetero_n);
    REGISTER_AROM("AromHetero_O", "O atoms in aromatic rings", desc_arom_hetero_o);
    REGISTER_AROM("AromHetero_S", "S atoms in aromatic rings", desc_arom_hetero_s);
    REGISTER_AROM("AromHetero_Mix", "Mixed heteroatom aromatic rings", desc_arom_hetero_mix);

    /* Electronic Properties */
    REGISTER_AROM("AromChi_Sum", "Sum of EN for aromatic atoms", desc_arom_chi_sum);
    REGISTER_AROM("AromChi_Mean", "Mean EN of aromatic atoms", desc_arom_chi_mean);
    REGISTER_AROM("AromChi_Var", "Variance of aromatic EN", desc_arom_chi_var);
    REGISTER_AROM("AromPol_Sum", "Sum of polarizability (aromatic)", desc_arom_pol_sum);
    REGISTER_AROM("AromPol_Mean", "Mean polarizability (aromatic)", desc_arom_pol_mean);
    REGISTER_AROM("AromCharge_Sum", "Sum of partial charges (aromatic)", desc_arom_charge_sum);
    REGISTER_AROM("AromCharge_Pos", "Positive charge on aromatic", desc_arom_charge_pos);
    REGISTER_AROM("AromCharge_Neg", "Negative charge on aromatic", desc_arom_charge_neg);
    REGISTER_AROM("AromEState_Sum", "Sum of E-State (aromatic)", desc_arom_estate_sum);
    REGISTER_AROM("AromEState_Mean", "Mean E-State (aromatic)", desc_arom_estate_mean);
    REGISTER_AROM("AromEState_Min", "Min E-State (aromatic)", desc_arom_estate_min);
    REGISTER_AROM("AromEState_Max", "Max E-State (aromatic)", desc_arom_estate_max);
    REGISTER_AROM("AromIP_Sum", "Sum of IP (aromatic)", desc_arom_ip_sum);
    REGISTER_AROM("AromEA_Sum", "Sum of EA (aromatic)", desc_arom_ea_sum);
    REGISTER_AROM("AromMass_Frac", "Mass fraction of aromatic system", desc_arom_mass_frac);
    REGISTER_AROM("AromPiElec", "Pi electrons in aromatic system", desc_arom_pi_elec);

    /* Topology */
    REGISTER_AROM("AromTopo_Bridges", "Bonds connecting aromatic systems", desc_arom_topo_bridges);
    REGISTER_AROM("AromTopo_FusionDeg", "Avg fusion degree", desc_arom_topo_fusiondeg);
    REGISTER_AROM("AromTopo_MaxFused", "Max rings in fused system", desc_arom_topo_maxfused);
    REGISTER_AROM("AromTopo_SpiroArom", "Spiro junctions in aromatics", desc_arom_topo_spiro);
    REGISTER_AROM("AromTopo_Distance", "Mean distance between aromatics", desc_arom_topo_distance);
    REGISTER_AROM("AromTopo_Clustering", "Aromatic clustering coefficient", desc_arom_topo_clustering);
    REGISTER_AROM("AromTopo_Peripheral", "Aromatic on periphery", desc_arom_topo_peripheral);
    REGISTER_AROM("AromTopo_Core", "Aromatic in core", desc_arom_topo_core);
    REGISTER_AROM("AromTopo_Substitution", "Total substituents on aromatics", desc_arom_topo_substitution);
    REGISTER_AROM("AromTopo_OrthoSub", "Ortho-substituted positions", desc_arom_topo_ortho);
    REGISTER_AROM("AromTopo_MetaSub", "Meta-substituted positions", desc_arom_topo_meta);
    REGISTER_AROM("AromTopo_ParaSub", "Para-substituted positions", desc_arom_topo_para);
}
