/**
 * @file energetics.c
 * @brief Energetic molecular descriptors
 *
 * Ultra-fast energetic descriptors based on thermodynamic and quantum mechanical
 * concepts. Includes Born solvation proxies, FMO-based reactivity indices,
 * Hansen solubility parameters, and entropy-related descriptors.
 *
 * Group A: Born Solvation Proxies
 * Group B: Electronic Hardness & Reactivity (FMO Theory)
 * Group C: Hydrophobic Effect & Cavity Formation
 * Group D: Hansen Solubility Parameter Components
 * Group E: Linear Moments of Distribution
 * Group F: Rotational & Vibrational Entropy
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

static const double ELECTRONEGATIVITY[] = {
    [ELEM_H]  = 2.20, [ELEM_C]  = 2.55, [ELEM_N]  = 3.04, [ELEM_O]  = 3.44,
    [ELEM_F]  = 3.98, [ELEM_P]  = 2.19, [ELEM_S]  = 2.58, [ELEM_Cl] = 3.16,
    [ELEM_Br] = 2.96, [ELEM_I]  = 2.66, [ELEM_Si] = 1.90, [ELEM_B]  = 2.04,
    [ELEM_Se] = 2.55, [ELEM_As] = 2.18,
};

static const double VDW_RADIUS[] = {
    [ELEM_H]  = 1.20, [ELEM_C]  = 1.70, [ELEM_N]  = 1.55, [ELEM_O]  = 1.52,
    [ELEM_F]  = 1.47, [ELEM_P]  = 1.80, [ELEM_S]  = 1.80, [ELEM_Cl] = 1.75,
    [ELEM_Br] = 1.85, [ELEM_I]  = 1.98, [ELEM_Si] = 2.10, [ELEM_B]  = 1.92,
    [ELEM_Se] = 1.90, [ELEM_As] = 1.85,
};

static const double POLARIZABILITY[] = {
    [ELEM_H]  = 0.667, [ELEM_C]  = 1.76, [ELEM_N]  = 1.10, [ELEM_O]  = 0.802,
    [ELEM_F]  = 0.557, [ELEM_P]  = 3.63, [ELEM_S]  = 2.90, [ELEM_Cl] = 2.18,
    [ELEM_Br] = 3.05, [ELEM_I]  = 4.70, [ELEM_Si] = 5.38, [ELEM_B]  = 3.03,
    [ELEM_Se] = 3.77, [ELEM_As] = 4.31,
};

static const double MOLAR_REFRACTIVITY[] = {
    [ELEM_H]  = 1.057, [ELEM_C]  = 2.503, [ELEM_N]  = 2.262, [ELEM_O]  = 1.607,
    [ELEM_F]  = 0.997, [ELEM_P]  = 6.920, [ELEM_S]  = 7.365, [ELEM_Cl] = 5.853,
    [ELEM_Br] = 8.927, [ELEM_I]  = 13.900, [ELEM_Si] = 6.840, [ELEM_B]  = 3.710,
    [ELEM_Se] = 7.990, [ELEM_As] = 6.280,
};

static const double IONIZATION_POTENTIAL[] = {
    [ELEM_H]  = 13.598, [ELEM_C]  = 11.260, [ELEM_N]  = 14.534, [ELEM_O]  = 13.618,
    [ELEM_F]  = 17.422, [ELEM_P]  = 10.486, [ELEM_S]  = 10.360, [ELEM_Cl] = 12.967,
    [ELEM_Br] = 11.814, [ELEM_I]  = 10.451, [ELEM_Si] = 8.151, [ELEM_B]  = 8.298,
    [ELEM_Se] = 9.752, [ELEM_As] = 9.815,
};

static const double ELECTRON_AFFINITY[] = {
    [ELEM_H]  = 0.754, [ELEM_C]  = 1.262, [ELEM_N]  = -0.07, [ELEM_O]  = 1.461,
    [ELEM_F]  = 3.401, [ELEM_P]  = 0.746, [ELEM_S]  = 2.077, [ELEM_Cl] = 3.612,
    [ELEM_Br] = 3.364, [ELEM_I]  = 3.059, [ELEM_Si] = 1.389, [ELEM_B]  = 0.277,
    [ELEM_Se] = 2.020, [ELEM_As] = 0.804,
};

static const double HYDRATION_ENERGY[] = {
    [ELEM_H]  = 0.0, [ELEM_C]  = 8.4, [ELEM_N]  = -20.0, [ELEM_O]  = -25.0,
    [ELEM_F]  = -15.0, [ELEM_P]  = -10.0, [ELEM_S]  = -5.0, [ELEM_Cl] = 5.0,
    [ELEM_Br] = 7.0, [ELEM_I]  = 10.0, [ELEM_Si] = 5.0, [ELEM_B]  = -5.0,
    [ELEM_Se] = 3.0, [ELEM_As] = 5.0,
};

static const double ATOMIC_MASS[] = {
    [ELEM_H]  = 1.008, [ELEM_C]  = 12.011, [ELEM_N]  = 14.007, [ELEM_O]  = 15.999,
    [ELEM_F]  = 18.998, [ELEM_P]  = 30.974, [ELEM_S]  = 32.065, [ELEM_Cl] = 35.453,
    [ELEM_Br] = 79.904, [ELEM_I]  = 126.904, [ELEM_Si] = 28.086, [ELEM_B]  = 10.811,
    [ELEM_Se] = 78.96, [ELEM_As] = 74.922,
};

/* Helper functions */
#define GET_PROP(arr, elem, def) \
    ((elem > 0 && elem < (int)(sizeof(arr)/sizeof(arr[0])) && arr[elem] != 0.0) ? arr[elem] : def)

static inline double get_chi(element_t e) { return GET_PROP(ELECTRONEGATIVITY, e, 2.55); }
static inline double get_radius(element_t e) { return GET_PROP(VDW_RADIUS, e, 1.70); }
static inline double get_pol(element_t e) { return GET_PROP(POLARIZABILITY, e, 1.76); }
static inline double get_mr(element_t e) { return GET_PROP(MOLAR_REFRACTIVITY, e, 2.503); }
static inline double get_ip(element_t e) { return GET_PROP(IONIZATION_POTENTIAL, e, 11.26); }
static inline double get_ea(element_t e) { return ELECTRON_AFFINITY[e]; }
static inline double get_hydration(element_t e) { return HYDRATION_ENERGY[e]; }
static inline double get_mass(element_t e) { return GET_PROP(ATOMIC_MASS, e, 12.011); }

static inline int is_hydrophobic(element_t e) {
    return (e == ELEM_C || e == ELEM_Cl || e == ELEM_Br || e == ELEM_I);
}

static inline int is_polar(element_t e) {
    return (e == ELEM_N || e == ELEM_O || e == ELEM_F || e == ELEM_S);
}

static int is_hbond_donor(const molecule_t* mol, int idx) {
    const atom_t* a = &mol->atoms[idx];
    if (a->element != ELEM_N && a->element != ELEM_O) return 0;
    if (a->implicit_h_count > 0) return 1;
    for (int i = 0; i < a->num_neighbors; i++) {
        if (mol->atoms[a->neighbors[i]].element == ELEM_H) return 1;
    }
    return 0;
}

static inline int is_hbond_acceptor(element_t e) {
    return (e == ELEM_N || e == ELEM_O || e == ELEM_F);
}

/* Statistics structure */
typedef struct {
    int n_heavy, n_total, n_rotatable, n_ring_atoms, n_double_bonds;
    int n_polar, n_hydrophobic, n_hbond_donor, n_hbond_acceptor;
    double sum_chi, sum_chi_sq_over_r, sum_radius, sum_pol, sum_mr;
    double sum_ip, sum_ea, sum_hydration, sum_mass, sum_polar_mass;
    double local_dipole;
} energetic_stats_t;

static void collect_stats(const molecule_t* mol, energetic_stats_t* s) {
    memset(s, 0, sizeof(energetic_stats_t));

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        element_t e = a->element;

        double chi = get_chi(e), r = get_radius(e), pol = get_pol(e);
        double mr = get_mr(e), ip = get_ip(e), ea = get_ea(e);
        double mass = get_mass(e);

        s->n_total++;
        s->sum_chi += chi;
        s->sum_radius += r;
        s->sum_pol += pol;
        s->sum_mr += mr;
        s->sum_ip += ip;
        s->sum_ea += ea;
        s->sum_mass += mass;

        if (e == ELEM_H) continue;

        s->n_heavy++;
        s->sum_chi_sq_over_r += (chi * chi) / r;
        s->sum_hydration += get_hydration(e);

        int h_impl = a->implicit_h_count;
        s->n_total += h_impl;
        s->sum_chi += h_impl * get_chi(ELEM_H);
        s->sum_radius += h_impl * get_radius(ELEM_H);
        s->sum_pol += h_impl * get_pol(ELEM_H);
        s->sum_mr += h_impl * get_mr(ELEM_H);
        s->sum_ip += h_impl * get_ip(ELEM_H);
        s->sum_ea += h_impl * get_ea(ELEM_H);
        s->sum_mass += h_impl * get_mass(ELEM_H);

        if (is_polar(e)) { s->n_polar++; s->sum_polar_mass += mass; }
        if (is_hydrophobic(e)) s->n_hydrophobic++;
        if (is_hbond_donor(mol, i)) s->n_hbond_donor++;
        if (is_hbond_acceptor(e)) s->n_hbond_acceptor++;
        if (a->ring_count > 0) s->n_ring_atoms++;

        for (int j = 0; j < a->num_neighbors; j++) {
            int ni = a->neighbors[j];
            if (ni > i) {
                s->local_dipole += fabs(chi - get_chi(mol->atoms[ni].element));
            }
        }
    }

    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* b = &mol->bonds[i];
        if (b->type == BOND_DOUBLE) s->n_double_bonds++;
        if (b->type == BOND_SINGLE && !b->in_ring) {
            const atom_t* a1 = &mol->atoms[b->atom1];
            const atom_t* a2 = &mol->atoms[b->atom2];
            if (a1->element != ELEM_H && a2->element != ELEM_H) {
                int h1 = 0, h2 = 0;
                for (int j = 0; j < a1->num_neighbors; j++)
                    if (mol->atoms[a1->neighbors[j]].element != ELEM_H) h1++;
                for (int j = 0; j < a2->num_neighbors; j++)
                    if (mol->atoms[a2->neighbors[j]].element != ELEM_H) h2++;
                if (h1 > 1 && h2 > 1) s->n_rotatable++;
            }
        }
    }
}

/* ============================================================================
 * Batch Computation (thread-safe: stack-allocated stats)
 * ============================================================================ */

#define NUM_ENERGETIC_DESCRIPTORS 30

/* Batch compute ALL energetic descriptors - collects stats once */
int descriptors_compute_energetic_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    energetic_stats_t s;
    collect_stats(mol, &s);

    int idx = 0;

    /* Group A: Born Solvation Proxies (5) */
    values[idx++].d = s.sum_chi_sq_over_r;                                                      /* BornEnergyPotential */
    values[idx++].d = s.local_dipole;                                                           /* LocalDipolePotential */
    values[idx++].d = (s.sum_radius > 0) ? (double)s.n_polar / s.sum_radius : 0.0;             /* DielectricSaturation */
    values[idx++].d = s.sum_pol;                                                                /* PolarizabilityVolume */
    values[idx++].d = s.sum_mr;                                                                 /* MolarRefractivityEnergy */

    /* Group B: Electronic Hardness & Reactivity (5) */
    values[idx++].d = s.sum_ip;                                                                 /* TotalIonizationPotential */
    values[idx++].d = s.sum_ea;                                                                 /* TotalElectronAffinity */
    values[idx++].d = (s.sum_ip - s.sum_ea) / 2.0;                                             /* ChemicalHardness */
    double chi_mol = (s.sum_ip + s.sum_ea) / 2.0;
    double eta_mol = (s.sum_ip - s.sum_ea) / 2.0;
    values[idx++].d = (eta_mol > 0.001) ? (chi_mol * chi_mol) / (2.0 * eta_mol) : 0.0;        /* ElectrophilicityIndex */
    values[idx++].d = chi_mol;                                                                  /* AbsoluteElectronegativity */

    /* Group C: Hydrophobic Effect & Cavity Formation (5) */
    values[idx++].d = s.sum_radius * s.sum_radius * 0.072;                                     /* CavityEnergyProxy */
    values[idx++].d = (s.sum_radius > 0) ? (double)s.n_hydrophobic / s.sum_radius : 0.0;      /* HydrophobicSurfaceDensity */
    values[idx++].d = s.sum_hydration;                                                          /* HydrationEnergyIndex */

    /* AmphiphilicMoment - need to compute separately */
    double hc = 0, pc = 0; int nh = 0, np = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (is_hydrophobic(e)) { hc += i; nh++; }
        if (is_polar(e)) { pc += i; np++; }
    }
    if (nh > 0) hc /= nh;
    if (np > 0) pc /= np;
    values[idx++].d = fabs(hc - pc);                                                            /* AmphiphilicMoment */

    /* CarbonScaleEntropy */
    double entropy = 0.0;
    double ln_c = log(get_mass(ELEM_C));
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_C) entropy += ln_c;
    }
    values[idx++].d = entropy;                                                                  /* CarbonScaleEntropy */

    /* Group D: Hansen Solubility Parameters (5) */
    double hansen_d = 0.0, hansen_p = 0.0, hansen_h = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (e != ELEM_H) hansen_d += get_pol(e);
        if (e != ELEM_H && e != ELEM_C) hansen_p += get_chi(e) - 2.55;
        if (is_hbond_donor(mol, i)) hansen_h += 5.0;
        else if (is_hbond_acceptor(e)) hansen_h += 3.0;
    }
    values[idx++].d = hansen_d;                                                                 /* HansenDispersion */
    values[idx++].d = hansen_p;                                                                 /* HansenPolar */
    values[idx++].d = hansen_h;                                                                 /* HansenHBond */
    values[idx++].d = hansen_d*hansen_d + hansen_p*hansen_p + hansen_h*hansen_h;              /* CohesiveEnergyDensity */
    values[idx++].d = (s.sum_mass > 0) ? 20.0 * s.sum_polar_mass / s.sum_mass : 0.0;          /* HLB */

    /* Group E: Linear Moments (5) */
    double wsum = 0, msum = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double m = get_mass(mol->atoms[i].element);
        wsum += i * m; msum += m;
    }
    double com = (msum > 0) ? wsum / msum : 0.0;
    values[idx++].d = com;                                                                      /* LinearCenterOfMass */

    wsum = 0; double csum = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double c = get_chi(mol->atoms[i].element);
        wsum += i * c; csum += c;
    }
    values[idx++].d = (csum > 0) ? wsum / csum : 0.0;                                          /* LinearCenterOfPolarity */

    /* Moment of inertia */
    double moi = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double diff = i - com;
        moi += get_mass(mol->atoms[i].element) * diff * diff;
    }
    values[idx++].d = moi;                                                                      /* LinearMomentOfInertia */

    /* Polar spread */
    double p_sum = 0, p_sum_sq = 0; int p_n = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (is_polar(mol->atoms[i].element)) {
            p_sum += i; p_sum_sq += i * i; p_n++;
        }
    }
    if (p_n > 1) {
        double p_mean = p_sum / p_n;
        double p_var = (p_sum_sq / p_n) - (p_mean * p_mean);
        values[idx++].d = (p_var > 0) ? sqrt(p_var) : 0.0;
    } else {
        values[idx++].d = 0.0;
    }                                                                                           /* PolarDistributionSpread */

    /* Mass asymmetry */
    double m2 = 0, m3 = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double m = get_mass(mol->atoms[i].element);
        double diff = i - com;
        m2 += m * diff * diff;
        m3 += m * diff * diff * diff;
    }
    if (msum > 0) { m2 /= msum; m3 /= msum; }
    double std_m = sqrt(m2);
    values[idx++].d = (std_m > 0.001) ? m3 / (std_m * std_m * std_m) : 0.0;                   /* MassDistributionAsymmetry */

    /* Group F: Rotational & Vibrational Entropy (5) */
    values[idx++].d = (s.sum_mass > 0) ? (s.n_rotatable * 2.479) / s.sum_mass : 0.0;          /* RotationalFreedomCost */
    int rigid = s.n_ring_atoms + s.n_double_bonds;
    values[idx++].d = (s.n_total > 0) ? (double)rigid / s.n_total : 0.0;                      /* RigidityPenalty */

    /* Kier Flexibility */
    int n_heavy = 0, n_bonds = 0;
    for (int i = 0; i < mol->num_atoms; i++)
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* b = &mol->bonds[i];
        if (mol->atoms[b->atom1].element != ELEM_H && mol->atoms[b->atom2].element != ELEM_H)
            n_bonds++;
    }
    double k1 = 0, k2 = 0;
    if (n_bonds > 0) k1 = (double)n_heavy * (n_heavy - 1) * (n_heavy - 1) / ((double)n_bonds * n_bonds);
    if (n_heavy > 2 && n_bonds > 1) {
        double p2 = n_bonds * 0.8;
        if (p2 > 0) k2 = (double)(n_heavy-1)*(n_heavy-1)*(n_heavy-2)*(n_heavy-2) / (p2*p2);
    }
    values[idx++].d = k1 * k2;                                                                  /* KierFlexibility */

    /* Terminal mass fraction */
    double term = 0, total = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        double m = get_mass(a->element);
        total += m;
        if (a->element != ELEM_H) {
            int hn = 0;
            for (int j = 0; j < a->num_neighbors; j++)
                if (mol->atoms[a->neighbors[j]].element != ELEM_H) hn++;
            if (hn <= 1) term += m;
        }
    }
    values[idx++].d = (total > 0) ? term / total : 0.0;                                        /* TerminalMassFraction */

    /* Folding potential */
    int n = mol->num_atoms;
    if (n >= 4) {
        double q1 = n * 0.25, q3 = n * 0.75;
        int mid = 0, outer = 0;
        for (int i = 0; i < n; i++) {
            if (i >= q1 && i <= q3) mid++;
            else outer++;
        }
        values[idx++].d = (outer > 0) ? (double)mid / outer : (double)mid;
    } else {
        values[idx++].d = 1.0;
    }                                                                                           /* FoldingPotential */

    return NUM_ENERGETIC_DESCRIPTORS;
}

/* ============================================================================
 * Group A: Born Solvation Proxies (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_born_energy(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.sum_chi_sq_over_r;
    return CCHEM_OK;
}

static cchem_status_t desc_local_dipole(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.local_dipole;
    return CCHEM_OK;
}

static cchem_status_t desc_dielectric_sat(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = (s.sum_radius > 0) ? (double)s.n_polar / s.sum_radius : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_pol_volume(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.sum_pol;
    return CCHEM_OK;
}

static cchem_status_t desc_mr_energy(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.sum_mr;
    return CCHEM_OK;
}

/* ============================================================================
 * Group B: Electronic Hardness & Reactivity (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_total_ip(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.sum_ip;
    return CCHEM_OK;
}

static cchem_status_t desc_total_ea(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.sum_ea;
    return CCHEM_OK;
}

static cchem_status_t desc_chem_hardness(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = (s.sum_ip - s.sum_ea) / 2.0;
    return CCHEM_OK;
}

static cchem_status_t desc_electrophilicity(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    double chi = (s.sum_ip + s.sum_ea) / 2.0;
    double eta = (s.sum_ip - s.sum_ea) / 2.0;
    v->d = (eta > 0.001) ? (chi * chi) / (2.0 * eta) : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_abs_chi(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = (s.sum_ip + s.sum_ea) / 2.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group C: Hydrophobic Effect & Cavity Formation (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_cavity_energy(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.sum_radius * s.sum_radius * 0.072;
    return CCHEM_OK;
}

static cchem_status_t desc_hydrophobic_density(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = (s.sum_radius > 0) ? (double)s.n_hydrophobic / s.sum_radius : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_hydration_energy(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = s.sum_hydration;
    return CCHEM_OK;
}

static cchem_status_t desc_amphiphilic(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double hc = 0, pc = 0; int nh = 0, np = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (is_hydrophobic(e)) { hc += i; nh++; }
        if (is_polar(e)) { pc += i; np++; }
    }
    if (nh > 0) hc /= nh;
    if (np > 0) pc /= np;
    v->d = fabs(hc - pc);
    return CCHEM_OK;
}

static cchem_status_t desc_carbon_entropy(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double entropy = 0.0;
    double ln_c = log(get_mass(ELEM_C));
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_C) entropy += ln_c;
    }
    v->d = entropy;
    return CCHEM_OK;
}

/* ============================================================================
 * Group D: Hansen Solubility Parameters (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_hansen_d(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double d = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (e != ELEM_H) d += get_pol(e);
    }
    v->d = d;
    return CCHEM_OK;
}

static cchem_status_t desc_hansen_p(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double p = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        if (e != ELEM_H && e != ELEM_C) p += get_chi(e) - 2.55;
    }
    v->d = p;
    return CCHEM_OK;
}

static cchem_status_t desc_hansen_h(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double h = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (is_hbond_donor(mol, i)) h += 5.0;
        else if (is_hbond_acceptor(mol->atoms[i].element)) h += 3.0;
    }
    v->d = h;
    return CCHEM_OK;
}

static cchem_status_t desc_ced(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    descriptor_value_t d, p, h;
    desc_hansen_d(mol, &d);
    desc_hansen_p(mol, &p);
    desc_hansen_h(mol, &h);
    v->d = d.d*d.d + p.d*p.d + h.d*h.d;
    return CCHEM_OK;
}

static cchem_status_t desc_hlb(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = (s.sum_mass > 0) ? 20.0 * s.sum_polar_mass / s.sum_mass : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group E: Linear Moments (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_linear_com(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double wsum = 0, msum = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double m = get_mass(mol->atoms[i].element);
        wsum += i * m; msum += m;
    }
    v->d = (msum > 0) ? wsum / msum : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_linear_cop(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double wsum = 0, csum = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double c = get_chi(mol->atoms[i].element);
        wsum += i * c; csum += c;
    }
    v->d = (csum > 0) ? wsum / csum : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_linear_moi(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double wsum = 0, msum = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double m = get_mass(mol->atoms[i].element);
        wsum += i * m; msum += m;
    }
    double center = (msum > 0) ? wsum / msum : 0.0;
    double moi = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double diff = i - center;
        moi += get_mass(mol->atoms[i].element) * diff * diff;
    }
    v->d = moi;
    return CCHEM_OK;
}

static cchem_status_t desc_polar_spread(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double sum = 0, sum_sq = 0; int n = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (is_polar(mol->atoms[i].element)) {
            sum += i; sum_sq += i * i; n++;
        }
    }
    if (n > 1) {
        double mean = sum / n;
        double var = (sum_sq / n) - (mean * mean);
        v->d = (var > 0) ? sqrt(var) : 0.0;
    } else {
        v->d = 0.0;
    }
    return CCHEM_OK;
}

static cchem_status_t desc_mass_asym(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double sum = 0, mtot = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double m = get_mass(mol->atoms[i].element);
        sum += i * m; mtot += m;
    }
    double mean = (mtot > 0) ? sum / mtot : 0.0;
    double m2 = 0, m3 = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double m = get_mass(mol->atoms[i].element);
        double diff = i - mean;
        m2 += m * diff * diff;
        m3 += m * diff * diff * diff;
    }
    m2 /= mtot; m3 /= mtot;
    double std = sqrt(m2);
    v->d = (std > 0.001) ? m3 / (std * std * std) : 0.0;
    return CCHEM_OK;
}

/* ============================================================================
 * Group F: Rotational & Vibrational Entropy (5 descriptors)
 * ============================================================================ */

static cchem_status_t desc_rot_freedom(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    v->d = (s.sum_mass > 0) ? (s.n_rotatable * 2.479) / s.sum_mass : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_rigidity(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    energetic_stats_t s; collect_stats(mol, &s);
    int rigid = s.n_ring_atoms + s.n_double_bonds;
    v->d = (s.n_total > 0) ? (double)rigid / s.n_total : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_kier_flex(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    int n_heavy = 0, n_bonds = 0;
    for (int i = 0; i < mol->num_atoms; i++)
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* b = &mol->bonds[i];
        if (mol->atoms[b->atom1].element != ELEM_H && mol->atoms[b->atom2].element != ELEM_H)
            n_bonds++;
    }
    double k1 = 0, k2 = 0;
    if (n_bonds > 0) {
        k1 = (double)n_heavy * (n_heavy - 1) * (n_heavy - 1) / ((double)n_bonds * n_bonds);
    }
    if (n_heavy > 2 && n_bonds > 1) {
        double p2 = n_bonds * 0.8;
        if (p2 > 0) k2 = (double)(n_heavy-1)*(n_heavy-1)*(n_heavy-2)*(n_heavy-2) / (p2*p2);
    }
    v->d = k1 * k2;
    return CCHEM_OK;
}

static cchem_status_t desc_terminal_mass(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    double term = 0, total = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        double m = get_mass(a->element);
        total += m;
        if (a->element != ELEM_H) {
            int hn = 0;
            for (int j = 0; j < a->num_neighbors; j++)
                if (mol->atoms[a->neighbors[j]].element != ELEM_H) hn++;
            if (hn <= 1) term += m;
        }
    }
    v->d = (total > 0) ? term / total : 0.0;
    return CCHEM_OK;
}

static cchem_status_t desc_folding(const molecule_t* mol, descriptor_value_t* v) {
    if (!mol || !v) return CCHEM_ERROR_INVALID_INPUT;
    int n = mol->num_atoms;
    if (n < 4) { v->d = 1.0; return CCHEM_OK; }
    double q1 = n * 0.25, q3 = n * 0.75;
    int mid = 0, outer = 0;
    for (int i = 0; i < n; i++) {
        if (i >= q1 && i <= q3) mid++;
        else outer++;
    }
    v->d = (outer > 0) ? (double)mid / outer : (double)mid;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REG(n, d, f) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, n, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, d, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_ENERGETIC; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = f; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_energetic(void) {
    /* Group A: Born Solvation Proxies */
    REG("BornEnergyPotential", "Total Born solvation energy proxy (chi^2/r)", desc_born_energy);
    REG("LocalDipolePotential", "Sum of bond electronegativity differences", desc_local_dipole);
    REG("DielectricSaturation", "Polar atom density (n_polar/volume)", desc_dielectric_sat);
    REG("PolarizabilityVolume", "Sum of atomic polarizabilities", desc_pol_volume);
    REG("MolarRefractivityEnergy", "Sum of atomic molar refractivity", desc_mr_energy);

    /* Group B: Electronic Hardness & Reactivity */
    REG("TotalIonizationPotential", "Sum of first ionization potentials", desc_total_ip);
    REG("TotalElectronAffinity", "Sum of electron affinities", desc_total_ea);
    REG("ChemicalHardness", "Global chemical hardness (IP-EA)/2", desc_chem_hardness);
    REG("ElectrophilicityIndex", "Global electrophilicity chi^2/2eta", desc_electrophilicity);
    REG("AbsoluteElectronegativity", "Mulliken electronegativity (IP+EA)/2", desc_abs_chi);

    /* Group C: Hydrophobic Effect & Cavity Formation */
    REG("CavityEnergyProxy", "Cavity formation cost (r^2 * gamma)", desc_cavity_energy);
    REG("HydrophobicSurfaceDensity", "Hydrophobic atoms / total radius", desc_hydrophobic_density);
    REG("HydrationEnergyIndex", "Sum of hydration free energies", desc_hydration_energy);
    REG("AmphiphilicMoment", "Polar-hydrophobic center distance (1D)", desc_amphiphilic);
    REG("CarbonScaleEntropy", "Mass-weighted entropy proxy (C atoms)", desc_carbon_entropy);

    /* Group D: Hansen Solubility Parameters */
    REG("HansenDispersion", "Dispersion solubility parameter", desc_hansen_d);
    REG("HansenPolar", "Polar solubility parameter", desc_hansen_p);
    REG("HansenHBond", "H-bonding solubility parameter", desc_hansen_h);
    REG("CohesiveEnergyDensity", "Total cohesive energy (d^2+p^2+h^2)", desc_ced);
    REG("HLB", "Hydrophilic-lipophilic balance (20*polar/total)", desc_hlb);

    /* Group E: Linear Moments */
    REG("LinearCenterOfMass", "Mass-weighted center position", desc_linear_com);
    REG("LinearCenterOfPolarity", "Chi-weighted center position", desc_linear_cop);
    REG("LinearMomentOfInertia", "1D moment of inertia", desc_linear_moi);
    REG("PolarDistributionSpread", "Std dev of polar atom positions", desc_polar_spread);
    REG("MassDistributionAsymmetry", "Skewness of mass distribution", desc_mass_asym);

    /* Group F: Rotational & Vibrational Entropy */
    REG("RotationalFreedomCost", "Entropy cost of rotatable bonds", desc_rot_freedom);
    REG("RigidityPenalty", "(Rings + double bonds) / atoms", desc_rigidity);
    REG("KierFlexibility", "Kier kappa1*kappa2 flexibility index", desc_kier_flex);
    REG("TerminalMassFraction", "Mass of terminal atoms / total", desc_terminal_mass);
    REG("FoldingPotential", "Middle 50% / outer 25% atom ratio", desc_folding);
}
