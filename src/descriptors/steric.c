/**
 * @file steric.c
 * @brief Steric (volume and surface) molecular descriptors
 *
 * Ultra-fast steric descriptors using additive atomic contributions.
 * All volumes use McGowan VdW volumes, surface areas use Ertl TPSA fragments.
 *
 * Group A: Additive Volume & Refractivity
 * Group B: Surface Areas (TPSA Family)
 * Group C: Volume Partitioning
 * Group D: Shape & Compactness Topology
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Physical Constants: Atomic Volumes and Radii
 * ============================================================================ */

/* McGowan characteristic volumes in cm³/mol (divided by 100 for convenience) */
static const double MCGOWAN_VOLUME[] = {
    [ELEM_H]  = 0.0871,
    [ELEM_C]  = 0.1635,
    [ELEM_N]  = 0.1402,
    [ELEM_O]  = 0.1208,
    [ELEM_F]  = 0.1048,
    [ELEM_P]  = 0.2456,
    [ELEM_S]  = 0.2265,
    [ELEM_Cl] = 0.2053,
    [ELEM_Br] = 0.2571,
    [ELEM_I]  = 0.3415,
    [ELEM_Si] = 0.2662,
    [ELEM_B]  = 0.1853,
    [ELEM_Se] = 0.2692,
    [ELEM_As] = 0.2584,
};

/* Atomic molar refractivity (Å³) */
static const double ATOMIC_MR[] = {
    [ELEM_H]  = 1.057,
    [ELEM_C]  = 2.503,
    [ELEM_N]  = 2.262,
    [ELEM_O]  = 1.607,
    [ELEM_F]  = 0.997,
    [ELEM_P]  = 6.920,
    [ELEM_S]  = 7.365,
    [ELEM_Cl] = 5.853,
    [ELEM_Br] = 8.927,
    [ELEM_I]  = 13.900,
    [ELEM_Si] = 6.840,
    [ELEM_B]  = 3.710,
    [ELEM_Se] = 7.990,
};

/* Approximate atomic surface area (4*pi*r²) in Ų */
static const double ATOMIC_SURFACE[] = {
    [ELEM_H]  = 18.10,
    [ELEM_C]  = 36.32,
    [ELEM_N]  = 30.19,
    [ELEM_O]  = 29.03,
    [ELEM_F]  = 27.14,
    [ELEM_P]  = 40.72,
    [ELEM_S]  = 40.72,
    [ELEM_Cl] = 38.48,
    [ELEM_Br] = 43.01,
    [ELEM_I]  = 49.26,
    [ELEM_Si] = 55.42,
    [ELEM_B]  = 46.33,
};

static inline double get_mcgowan_volume(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(MCGOWAN_VOLUME)/sizeof(MCGOWAN_VOLUME[0])))
        return 0.1635;
    double v = MCGOWAN_VOLUME[elem];
    return (v > 0.0) ? v : 0.1635;
}

static inline double get_atomic_mr(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ATOMIC_MR)/sizeof(ATOMIC_MR[0])))
        return 2.503;
    double mr = ATOMIC_MR[elem];
    return (mr > 0.0) ? mr : 2.503;
}

static inline double get_atomic_surface(element_t elem) {
    if (elem <= 0 || elem >= (int)(sizeof(ATOMIC_SURFACE)/sizeof(ATOMIC_SURFACE[0])))
        return 36.32;
    double s = ATOMIC_SURFACE[elem];
    return (s > 0.0) ? s : 36.32;
}

/* ============================================================================
 * TPSA Fragment Contributions (Ertl method)
 * ============================================================================ */

static double get_tpsa_contribution(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    element_t elem = atom->element;

    if (elem != ELEM_N && elem != ELEM_O && elem != ELEM_S && elem != ELEM_P) {
        return 0.0;
    }

    int h_count = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
    }

    int double_bonds = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE) double_bonds++;
    }

    if (elem == ELEM_N) {
        if (atom->aromatic) {
            return (h_count == 1) ? 15.79 : 12.89;
        }
        if (atom->charge > 0) {
            if (h_count == 3) return 27.64;
            else if (h_count == 2) return 25.59;
            else if (h_count == 1) return 23.79;
            else return 16.61;
        }
        if (double_bonds > 0) {
            return (h_count == 1) ? 23.85 : 12.36;
        }
        if (h_count == 2) return 26.03;
        else if (h_count == 1) return 12.03;
        else return 3.24;
    }

    if (elem == ELEM_O) {
        if (atom->aromatic) return 13.14;
        if (atom->charge < 0) return 23.06;
        if (double_bonds > 0) return 17.07;
        if (h_count >= 1) return 20.23;
        return 9.23;
    }

    if (elem == ELEM_S) {
        if (atom->aromatic) return 28.24;
        if (double_bonds == 2) return 36.28;
        if (double_bonds == 1) return 32.09;
        if (h_count >= 1) return 38.80;
        return 25.30;
    }

    if (elem == ELEM_P) {
        return 13.59;
    }

    return 0.0;
}

/* ============================================================================
 * Helper: Collect steric statistics
 * ============================================================================ */

typedef struct {
    double total_volume;
    double heavy_volume;
    double total_mr;
    double total_surface;
    double tpsa;
    double aromatic_volume;
    double aliphatic_volume;
    double halogen_volume;
    double terminal_volume;
    double branch_volume;
    double stereo_volume;
    int n_total;
    int n_bonds;
    int n_ring_atoms;
} steric_stats_t;

static void collect_steric_stats(const molecule_t* mol, steric_stats_t* stats) {
    memset(stats, 0, sizeof(steric_stats_t));

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        element_t elem = atom->element;

        double vol = get_mcgowan_volume(elem);
        double mr = get_atomic_mr(elem);
        double surf = get_atomic_surface(elem);

        if (elem == ELEM_H) {
            stats->total_volume += vol;
            stats->total_mr += mr;
            stats->total_surface += surf * 0.5;
            continue;
        }

        stats->n_total++;
        stats->total_volume += vol;
        stats->heavy_volume += vol;
        stats->total_mr += mr;
        stats->total_surface += surf;

        int h_implicit = atom->implicit_h_count;
        stats->total_volume += h_implicit * get_mcgowan_volume(ELEM_H);
        stats->total_mr += h_implicit * get_atomic_mr(ELEM_H);
        stats->total_surface += h_implicit * get_atomic_surface(ELEM_H) * 0.5;

        stats->tpsa += get_tpsa_contribution(mol, i);

        if (atom->aromatic) {
            stats->aromatic_volume += vol;
        } else {
            stats->aliphatic_volume += vol;
        }

        if (elem == ELEM_F || elem == ELEM_Cl || elem == ELEM_Br || elem == ELEM_I) {
            stats->halogen_volume += vol;
        }

        if (atom->num_neighbors == 1) {
            stats->terminal_volume += vol;
        }

        if (atom->num_neighbors >= 3) {
            stats->branch_volume += vol;
        }

        if (atom->chirality != CHIRALITY_NONE) {
            stats->stereo_volume += vol;
        }

        if (atom->ring_count > 0) {
            stats->n_ring_atoms++;
        }
    }

    stats->n_bonds = mol->num_bonds;
}

/* ============================================================================
 * Batch Computation (thread-safe: stack-allocated stats)
 * ============================================================================ */

#define NUM_STERIC_DESCRIPTORS 20

/* Batch compute ALL steric descriptors - collects stats once */
int descriptors_compute_steric_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    steric_stats_t stats;
    collect_steric_stats(mol, &stats);

    int idx = 0;

    /* Group A: Additive Volume & Refractivity (5) */
    values[idx++].d = stats.total_volume;                                                       /* TotalVdWVolume */
    values[idx++].d = stats.heavy_volume;                                                       /* HeavyAtomVolume */
    values[idx++].d = stats.total_mr;                                                           /* MolarRefractivityS */
    values[idx++].d = (stats.n_total > 0) ? stats.heavy_volume / stats.n_total : 0.0;          /* MeanAtomicVolume */
    values[idx++].d = (stats.total_volume > 0) ? cbrt((3.0 * stats.total_volume) / (4.0 * M_PI)) : 0.0; /* EffectiveRadius */

    /* Group B: Surface Areas (TPSA Family) (5) */
    values[idx++].d = stats.tpsa;                                                               /* ClassicalTPSA */
    values[idx++].d = stats.total_surface - stats.tpsa;                                         /* NonPolarSurface */
    values[idx++].d = stats.total_surface;                                                      /* TotalSurfaceArea */
    values[idx++].d = (stats.total_volume > 0) ? stats.tpsa / stats.total_volume : 0.0;        /* TPSADensityS */
    values[idx++].d = (stats.total_surface > 0) ? stats.tpsa / stats.total_surface : 0.0;      /* RelativePolarSurface */

    /* Group C: Volume Partitioning (5) */
    values[idx++].d = stats.aromatic_volume;                                                    /* AromaticVolume */
    values[idx++].d = stats.aliphatic_volume;                                                   /* AliphaticVolume */
    values[idx++].d = (stats.heavy_volume > 0) ? stats.aromatic_volume / stats.heavy_volume : 0.0; /* RingVolumeRatio */
    values[idx++].d = (stats.heavy_volume > 0) ? stats.aliphatic_volume / stats.heavy_volume : 0.0; /* SideChainVolumeRatio */
    values[idx++].d = stats.halogen_volume;                                                     /* HalogenBulk */

    /* Group D: Shape & Compactness Topology (5) */
    /* Hall-Kier Alpha */
    int n_heavy = 0;
    double hetero_correction = 0.0;
    static const double ALPHA_CORRECTION[] = {
        [ELEM_N]  = -0.20, [ELEM_O]  = -0.40, [ELEM_S]  = 0.00,
        [ELEM_F]  = -0.35, [ELEM_Cl] = 0.10,  [ELEM_Br] = 0.30,
        [ELEM_I]  = 0.50,  [ELEM_P]  = -0.10,
    };
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_H) continue;
        n_heavy++;
        if (elem < (int)(sizeof(ALPHA_CORRECTION)/sizeof(ALPHA_CORRECTION[0]))) {
            double corr = ALPHA_CORRECTION[elem];
            if (corr != 0.0) hetero_correction += corr;
        }
    }
    values[idx++].d = (double)(n_heavy - mol->num_rings - 1) + hetero_correction;              /* HallKierAlpha */

    /* Globularity */
    values[idx++].d = (stats.total_volume > 0) ? pow(stats.total_surface, 1.5) / stats.total_volume : 0.0; /* GlobularityProxy */
    values[idx++].d = stats.terminal_volume;                                                    /* TerminalVolume */
    values[idx++].d = stats.branch_volume;                                                      /* BranchingVolume */
    values[idx++].d = stats.stereo_volume;                                                      /* StereoVolume */

    return NUM_STERIC_DESCRIPTORS;
}

/* ============================================================================
 * Group A: Additive Volume & Refractivity (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_total_vdw_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.total_volume;
    return CCHEM_OK;
}

static cchem_status_t desc_heavy_atom_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.heavy_volume;
    return CCHEM_OK;
}

static cchem_status_t desc_molar_refractivity_steric(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.total_mr;
    return CCHEM_OK;
}

static cchem_status_t desc_mean_atomic_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = (stats.n_total > 0) ? stats.heavy_volume / stats.n_total : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_effective_radius(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    if (stats.total_volume > 0) {
        value->d = cbrt((3.0 * stats.total_volume) / (4.0 * M_PI));
    } else {
        value->d = 0.0;
    }
    return CCHEM_OK;
}

/* ============================================================================
 * Group B: Surface Areas (TPSA Family) (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_classical_tpsa(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.tpsa;
    return CCHEM_OK;
}

static cchem_status_t desc_non_polar_surface(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.total_surface - stats.tpsa;
    return CCHEM_OK;
}

static cchem_status_t desc_total_surface_area(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.total_surface;
    return CCHEM_OK;
}

static cchem_status_t desc_tpsa_density_steric(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = (stats.total_volume > 0) ? stats.tpsa / stats.total_volume : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_relative_polar_surface(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = (stats.total_surface > 0) ? stats.tpsa / stats.total_surface : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group C: Volume Partitioning (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_aromatic_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.aromatic_volume;
    return CCHEM_OK;
}

static cchem_status_t desc_aliphatic_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.aliphatic_volume;
    return CCHEM_OK;
}

static cchem_status_t desc_ring_volume_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = (stats.heavy_volume > 0) ? stats.aromatic_volume / stats.heavy_volume : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_sidechain_volume_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = (stats.heavy_volume > 0) ? stats.aliphatic_volume / stats.heavy_volume : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_halogen_bulk(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.halogen_volume;
    return CCHEM_OK;
}

/* ============================================================================
 * Group D: Shape & Compactness Topology (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_hall_kier_alpha(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int n_heavy = 0;
    int n_rings = mol->num_rings;
    double hetero_correction = 0.0;

    static const double ALPHA_CORRECTION[] = {
        [ELEM_N]  = -0.20,
        [ELEM_O]  = -0.40,
        [ELEM_S]  = 0.00,
        [ELEM_F]  = -0.35,
        [ELEM_Cl] = 0.10,
        [ELEM_Br] = 0.30,
        [ELEM_I]  = 0.50,
        [ELEM_P]  = -0.10,
    };

    for (int i = 0; i < mol->num_atoms; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_H) continue;
        n_heavy++;

        if (elem < (int)(sizeof(ALPHA_CORRECTION)/sizeof(ALPHA_CORRECTION[0]))) {
            double corr = ALPHA_CORRECTION[elem];
            if (corr != 0.0) {
                hetero_correction += corr;
            }
        }
    }

    double alpha = (double)(n_heavy - n_rings - 1) + hetero_correction;
    value->d = alpha;
    return CCHEM_OK;
}

static cchem_status_t desc_globularity_proxy(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);

    if (stats.total_volume > 0) {
        double surf_15 = pow(stats.total_surface, 1.5);
        value->d = surf_15 / stats.total_volume;
    } else {
        value->d = 0.0;
    }
    return CCHEM_OK;
}

static cchem_status_t desc_terminal_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.terminal_volume;
    return CCHEM_OK;
}

static cchem_status_t desc_branching_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.branch_volume;
    return CCHEM_OK;
}

static cchem_status_t desc_stereo_volume(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    steric_stats_t stats;
    collect_steric_stats(mol, &stats);
    value->d = stats.stereo_volume;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_STERIC_DESC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_STERIC; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_steric(void) {
    /* Group A: Additive Volume & Refractivity */
    REGISTER_STERIC_DESC("TotalVdWVolume", "Total van der Waals volume (McGowan)", desc_total_vdw_volume);
    REGISTER_STERIC_DESC("HeavyAtomVolume", "VdW volume of heavy atoms only", desc_heavy_atom_volume);
    REGISTER_STERIC_DESC("MolarRefractivityS", "Sum of atomic molar refractivity", desc_molar_refractivity_steric);
    REGISTER_STERIC_DESC("MeanAtomicVolume", "Average volume per heavy atom", desc_mean_atomic_volume);
    REGISTER_STERIC_DESC("EffectiveRadius", "Radius of equivalent sphere", desc_effective_radius);

    /* Group B: Surface Areas (TPSA Family) */
    REGISTER_STERIC_DESC("ClassicalTPSA", "Ertl polar surface area (N,O,S,P)", desc_classical_tpsa);
    REGISTER_STERIC_DESC("NonPolarSurface", "Total surface - TPSA", desc_non_polar_surface);
    REGISTER_STERIC_DESC("TotalSurfaceArea", "Approximate total surface area", desc_total_surface_area);
    REGISTER_STERIC_DESC("TPSADensityS", "TPSA / Volume (BBB predictor)", desc_tpsa_density_steric);
    REGISTER_STERIC_DESC("RelativePolarSurface", "TPSA / Total surface", desc_relative_polar_surface);

    /* Group C: Volume Partitioning */
    REGISTER_STERIC_DESC("AromaticVolume", "Volume of aromatic atoms", desc_aromatic_volume);
    REGISTER_STERIC_DESC("AliphaticVolume", "Volume of aliphatic atoms", desc_aliphatic_volume);
    REGISTER_STERIC_DESC("RingVolumeRatio", "Aromatic volume / total volume", desc_ring_volume_ratio);
    REGISTER_STERIC_DESC("SideChainVolumeRatio", "Aliphatic volume / total volume", desc_sidechain_volume_ratio);
    REGISTER_STERIC_DESC("HalogenBulk", "Volume of halogen atoms", desc_halogen_bulk);

    /* Group D: Shape & Compactness Topology */
    REGISTER_STERIC_DESC("HallKierAlpha", "Hall-Kier molecular flexibility index", desc_hall_kier_alpha);
    REGISTER_STERIC_DESC("GlobularityProxy", "Surface^1.5 / Volume (shape index)", desc_globularity_proxy);
    REGISTER_STERIC_DESC("TerminalVolume", "Volume of terminal atoms", desc_terminal_volume);
    REGISTER_STERIC_DESC("BranchingVolume", "Volume of branch point atoms", desc_branching_volume);
    REGISTER_STERIC_DESC("StereoVolume", "Volume of chiral centers", desc_stereo_volume);
}
