/**
 * @file counts.c
 * @brief Count-based molecular descriptors
 *
 * This module provides optimized count-based descriptors computed from
 * parsed molecule structures. All counts are computed in O(n) time
 * where n is the number of atoms or bonds.
 *
 * Descriptors:
 * - Element counts: CarbonCount, HydrogenCount, OxygenCount, NitrogenCount,
 *                   PhosphorusCount, SulfurCount, ChlorineCount, BromineCount,
 *                   IodineCount, FluorineCount
 * - AtomCount: Total heavy atoms (non-hydrogen)
 * - HeavyAtomCount: Same as AtomCount
 * - TotalBondCount: Total number of bonds
 * - SingleBondCount: Number of single bonds
 * - DoubleBondCount: Number of double bonds
 * - TripleBondCount: Number of triple bonds
 * - AromaticBondCount: Number of aromatic bonds
 * - RingCount: Number of rings (SSSR)
 * - AromaticRingCount: Number of aromatic rings
 */

#include <string.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Batch Computation - All counts in single pass (thread-safe)
 * ============================================================================ */

#define NUM_COUNT_DESCRIPTORS 83

/* Batch compute ALL count descriptors in minimal passes */
int descriptors_compute_counts_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    const int n_atoms = mol->num_atoms;
    const int n_bonds = mol->num_bonds;
    const atom_t* atoms = mol->atoms;
    const bond_t* bonds = mol->bonds;

    /* Initialize all to zero */
    for (int i = 0; i < NUM_COUNT_DESCRIPTORS; i++) values[i].i = 0;

    /* === SINGLE PASS THROUGH ATOMS === */
    int64_t c_C=0, c_H=0, c_O=0, c_N=0, c_P=0, c_S=0, c_Cl=0, c_Br=0, c_I=0, c_F=0;
    int64_t c_heavy=0, c_hetero=0, c_halogen=0, c_aromatic_atom=0;
    int64_t c_sp3_c=0, c_sp2_c=0, c_sp_c=0;
    int64_t c_ring_atom=0, c_chain_atom=0, c_terminal=0, c_branch=0, c_quat_c=0;
    int64_t c_bridgehead=0, c_spiro=0;
    int64_t c_hbd=0, c_hba=0;
    int64_t c_charged=0, c_pos=0, c_neg=0, net_charge=0;
    int64_t c_chiral=0;
    int64_t c_prim_amine=0, c_sec_amine=0, c_tert_amine=0, c_arom_n=0;
    int64_t c_hydroxyl=0, c_thiol=0;
    int64_t c_methyl=0, c_methylene=0, c_methine=0;
    int64_t c_aliphatic_c=0, c_aromatic_c=0;
    int64_t c_degree2=0, c_degree3=0, c_degree4=0;
    int64_t c_polar=0, c_acidic_h=0, c_explicit_h=0, c_conjugated=0;

    for (int i = 0; i < n_atoms; i++) {
        const atom_t* a = &atoms[i];
        element_t e = a->element;

        /* Element counts */
        switch (e) {
            case ELEM_C: c_C++; break;
            case ELEM_O: c_O++; break;
            case ELEM_N: c_N++; break;
            case ELEM_P: c_P++; break;
            case ELEM_S: c_S++; break;
            case ELEM_Cl: c_Cl++; break;
            case ELEM_Br: c_Br++; break;
            case ELEM_I: c_I++; break;
            case ELEM_F: c_F++; break;
            case ELEM_H: c_explicit_h++; break;
            default: break;
        }

        /* Hydrogen count (explicit + implicit) */
        if (e == ELEM_H) c_H++;
        else c_H += a->implicit_h_count;

        /* Heavy atom */
        if (e != ELEM_H) {
            c_heavy++;

            /* Heteroatom */
            if (e != ELEM_C) c_hetero++;

            /* Halogen */
            if (e == ELEM_F || e == ELEM_Cl || e == ELEM_Br || e == ELEM_I) c_halogen++;

            /* Polar */
            if (e == ELEM_O || e == ELEM_N || e == ELEM_S) c_polar++;

            /* Aromatic */
            if (a->aromatic) c_aromatic_atom++;

            /* Ring vs chain */
            if (a->ring_count > 0) c_ring_atom++;
            else c_chain_atom++;

            /* Degree counting */
            int heavy_deg = 0;
            for (int j = 0; j < a->num_neighbors; j++) {
                if (atoms[a->neighbors[j]].element != ELEM_H) heavy_deg++;
            }

            if (heavy_deg == 1) c_terminal++;
            if (heavy_deg >= 3) c_branch++;
            if (heavy_deg == 2) c_degree2++;
            if (heavy_deg == 3) c_degree3++;
            if (heavy_deg == 4) c_degree4++;

            /* Carbon specifics */
            if (e == ELEM_C) {
                if (a->aromatic) c_aromatic_c++;
                else c_aliphatic_c++;

                int total_deg = heavy_deg + a->implicit_h_count;
                if (total_deg == 4 && !a->aromatic) c_sp3_c++;
                else if (total_deg == 3 || a->aromatic) c_sp2_c++;
                else if (total_deg == 2) c_sp_c++;

                if (heavy_deg == 4) c_quat_c++;

                /* Methyl/methylene/methine */
                if (heavy_deg == 1 && a->implicit_h_count == 3) c_methyl++;
                if (heavy_deg == 2 && a->implicit_h_count == 2 && !a->aromatic) c_methylene++;
                if (heavy_deg == 3 && a->implicit_h_count == 1 && !a->aromatic) c_methine++;
            }

            /* Nitrogen specifics */
            if (e == ELEM_N) {
                if (a->aromatic) c_arom_n++;
                int h_count = a->implicit_h_count;
                for (int j = 0; j < a->num_neighbors; j++) {
                    if (atoms[a->neighbors[j]].element == ELEM_H) h_count++;
                }
                if (!a->aromatic) {
                    if (h_count == 2 && heavy_deg == 1) c_prim_amine++;
                    else if (h_count == 1 && heavy_deg == 2) c_sec_amine++;
                    else if (h_count == 0 && heavy_deg == 3) c_tert_amine++;
                }
            }

            /* H-bond donors/acceptors */
            if (e == ELEM_O || e == ELEM_N) {
                c_hba++;
                int h_count = a->implicit_h_count;
                for (int j = 0; j < a->num_neighbors; j++) {
                    if (atoms[a->neighbors[j]].element == ELEM_H) h_count++;
                }
                if (h_count > 0) c_hbd++;
                c_acidic_h += h_count;
            }

            /* Hydroxyl/thiol */
            if (e == ELEM_O && heavy_deg == 1 && a->implicit_h_count > 0) c_hydroxyl++;
            if (e == ELEM_S && heavy_deg == 1 && a->implicit_h_count > 0) c_thiol++;

            /* Bridgehead/spiro */
            if (a->ring_count >= 2 && heavy_deg >= 3) c_bridgehead++;
            if (a->ring_count >= 2 && heavy_deg == 4) c_spiro++;

            /* Conjugated (simplified: aromatic or has double bond) */
            if (a->aromatic) c_conjugated++;
        }

        /* Charge */
        if (a->charge != 0) {
            c_charged++;
            if (a->charge > 0) c_pos++;
            else c_neg++;
            net_charge += a->charge;
        }

        /* Chirality */
        if (a->chirality != CHIRALITY_NONE) c_chiral++;
    }

    /* === SINGLE PASS THROUGH BONDS === */
    int64_t c_single=0, c_double=0, c_triple=0, c_arom_bond=0, c_ring_bond=0, c_rot=0;
    int64_t c_stereo_db=0;
    int64_t c_carbonyl=0, c_carboxyl=0, c_ether=0, c_ester=0, c_amide=0;
    int64_t c_nitrile=0, c_nitro=0, c_sulfone=0;

    for (int i = 0; i < n_bonds; i++) {
        const bond_t* b = &bonds[i];
        const atom_t* a1 = &atoms[b->atom1];
        const atom_t* a2 = &atoms[b->atom2];

        switch (b->type) {
            case BOND_SINGLE:
            case BOND_UP:
            case BOND_DOWN:
            case BOND_RING_SINGLE:
                c_single++;
                break;
            case BOND_DOUBLE:
                c_double++;
                /* Check carbonyl */
                if ((a1->element == ELEM_C && a2->element == ELEM_O) ||
                    (a1->element == ELEM_O && a2->element == ELEM_C)) {
                    c_carbonyl++;
                }
                break;
            case BOND_TRIPLE:
                c_triple++;
                /* Check nitrile */
                if ((a1->element == ELEM_C && a2->element == ELEM_N) ||
                    (a1->element == ELEM_N && a2->element == ELEM_C)) {
                    c_nitrile++;
                }
                break;
            case BOND_AROMATIC:
                c_arom_bond++;
                break;
            default:
                break;
        }

        /* Ring bond */
        if (b->in_ring) c_ring_bond++;

        /* Rotatable (simplified: single, not in ring, heavy-heavy) */
        if ((b->type == BOND_SINGLE || b->type == BOND_UP || b->type == BOND_DOWN) &&
            !b->in_ring &&
            a1->element != ELEM_H && a2->element != ELEM_H &&
            a1->num_neighbors > 1 && a2->num_neighbors > 1) {
            c_rot++;
        }

        /* Stereo double bond */
        if (b->type == BOND_DOUBLE && b->stereo != STEREO_NONE) c_stereo_db++;

        /* Ether: C-O-C single bond */
        if (b->type == BOND_SINGLE &&
            ((a1->element == ELEM_O && a2->element == ELEM_C) ||
             (a1->element == ELEM_C && a2->element == ELEM_O))) {
            const atom_t* o_atom = (a1->element == ELEM_O) ? a1 : a2;
            int c_neighbors = 0;
            for (int j = 0; j < o_atom->num_neighbors; j++) {
                if (atoms[o_atom->neighbors[j]].element == ELEM_C) c_neighbors++;
            }
            if (c_neighbors == 2 && o_atom->implicit_h_count == 0) c_ether++;
        }
    }

    /* === RING COUNTS (from molecule structure) === */
    int64_t c_rings = mol->num_rings;
    int64_t c_arom_rings=0, c_aliph_rings=0, c_hetero_rings=0, c_arom_hetero=0;
    int64_t c_ring3=0, c_ring4=0, c_ring5=0, c_ring6=0, c_large_ring=0;
    int64_t c_arom5=0, c_arom6=0, c_sat_rings=0, c_carbocycle=0;

    for (int r = 0; r < mol->num_rings; r++) {
        const ring_t* ring = &mol->rings[r];
        bool is_arom = ring->aromatic;
        bool has_hetero = false;
        bool all_saturated = true;
        bool all_carbon = true;

        for (int j = 0; j < ring->size; j++) {
            const atom_t* ra = &atoms[ring->atoms[j]];
            if (ra->element != ELEM_C) {
                has_hetero = true;
                all_carbon = false;
            }
            if (ra->aromatic) all_saturated = false;
        }

        if (is_arom) c_arom_rings++;
        else c_aliph_rings++;

        if (has_hetero) {
            c_hetero_rings++;
            if (is_arom) c_arom_hetero++;
        }

        if (all_carbon) c_carbocycle++;
        if (all_saturated) c_sat_rings++;

        switch (ring->size) {
            case 3: c_ring3++; break;
            case 4: c_ring4++; break;
            case 5: c_ring5++; if (is_arom) c_arom5++; break;
            case 6: c_ring6++; if (is_arom) c_arom6++; break;
            default: if (ring->size > 6) c_large_ring++; break;
        }
    }

    /* === STORE ALL VALUES === */
    int idx = 0;
    /* Element counts (10) */
    values[idx++].i = c_C;
    values[idx++].i = c_H;
    values[idx++].i = c_O;
    values[idx++].i = c_N;
    values[idx++].i = c_P;
    values[idx++].i = c_S;
    values[idx++].i = c_Cl;
    values[idx++].i = c_Br;
    values[idx++].i = c_I;
    values[idx++].i = c_F;

    /* Heteroatom counts (2) */
    values[idx++].i = c_hetero;
    values[idx++].i = c_halogen;

    /* Atom counts (3) */
    values[idx++].i = c_heavy;
    values[idx++].i = c_heavy;  /* HeavyAtomCount same as AtomCount */
    values[idx++].i = c_aromatic_atom;

    /* Carbon hybridization (3) */
    values[idx++].i = c_sp3_c;
    values[idx++].i = c_sp2_c;
    values[idx++].i = c_sp_c;

    /* Bond counts (7) */
    values[idx++].i = n_bonds;
    values[idx++].i = c_single;
    values[idx++].i = c_double;
    values[idx++].i = c_triple;
    values[idx++].i = c_arom_bond;
    values[idx++].i = c_ring_bond;
    values[idx++].i = c_rot;

    /* Ring counts (5) */
    values[idx++].i = c_rings;
    values[idx++].i = c_arom_rings;
    values[idx++].i = c_aliph_rings;
    values[idx++].i = c_hetero_rings;
    values[idx++].i = c_arom_hetero;

    /* Ring size counts (5) */
    values[idx++].i = c_ring3;
    values[idx++].i = c_ring4;
    values[idx++].i = c_ring5;
    values[idx++].i = c_ring6;
    values[idx++].i = c_large_ring;

    /* Connectivity (7) */
    values[idx++].i = c_ring_atom;
    values[idx++].i = c_chain_atom;
    values[idx++].i = c_terminal;
    values[idx++].i = c_branch;
    values[idx++].i = c_quat_c;
    values[idx++].i = c_bridgehead;
    values[idx++].i = c_spiro;

    /* H-bonding (2) */
    values[idx++].i = c_hbd;
    values[idx++].i = c_hba;

    /* Charge (4) */
    values[idx++].i = c_charged;
    values[idx++].i = c_pos;
    values[idx++].i = c_neg;
    values[idx++].i = net_charge;

    /* Stereochemistry (2) */
    values[idx++].i = c_chiral;
    values[idx++].i = c_stereo_db;

    /* Fragments and unsaturation (2) */
    values[idx++].i = mol->num_fragments > 0 ? mol->num_fragments : 1;
    values[idx++].i = c_C - c_H/2 + c_N/2 + 1;  /* Simplified DBE */

    /* Amine classification (4) */
    values[idx++].i = c_prim_amine;
    values[idx++].i = c_sec_amine;
    values[idx++].i = c_tert_amine;
    values[idx++].i = c_arom_n;

    /* Functional groups (10) */
    values[idx++].i = c_hydroxyl;
    values[idx++].i = c_carbonyl;
    values[idx++].i = c_carboxyl;
    values[idx++].i = c_ether;
    values[idx++].i = c_ester;
    values[idx++].i = c_amide;
    values[idx++].i = c_nitrile;
    values[idx++].i = c_nitro;
    values[idx++].i = c_thiol;
    values[idx++].i = c_sulfone;

    /* Carbon types (5) */
    values[idx++].i = c_methyl;
    values[idx++].i = c_methylene;
    values[idx++].i = c_methine;
    values[idx++].i = c_aliphatic_c;
    values[idx++].i = c_aromatic_c;

    /* Ring classification (4) */
    values[idx++].i = c_arom5;
    values[idx++].i = c_arom6;
    values[idx++].i = c_sat_rings;
    values[idx++].i = c_carbocycle;

    /* Degree distribution (3) */
    values[idx++].i = c_degree2;
    values[idx++].i = c_degree3;
    values[idx++].i = c_degree4;

    /* Polarity and acidity (4) */
    values[idx++].i = c_polar;
    values[idx++].i = c_acidic_h;
    values[idx++].i = c_explicit_h;
    values[idx++].i = c_conjugated;

    return NUM_COUNT_DESCRIPTORS;
}

/* ============================================================================
 * Element Count Descriptors
 * Ultra-optimized: single pass through atoms array
 * ============================================================================ */

/* Generic element count function */
static cchem_status_t count_element(const molecule_t* mol, element_t element, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_atoms;
    const atom_t* atoms = mol->atoms;

    /* Tight loop for cache efficiency */
    for (int i = 0; i < n; i++) {
        count += (atoms[i].element == element);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Carbon count (includes implicit hydrogens for total H) */
static cchem_status_t desc_carbon_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_C, value);
}

/* Hydrogen count (explicit + implicit) */
static cchem_status_t desc_hydrogen_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_atoms;
    const atom_t* atoms = mol->atoms;

    for (int i = 0; i < n; i++) {
        if (atoms[i].element == ELEM_H) {
            count++;
        } else {
            /* Add implicit hydrogens */
            count += atoms[i].implicit_h_count;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

static cchem_status_t desc_oxygen_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_O, value);
}

static cchem_status_t desc_nitrogen_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_N, value);
}

static cchem_status_t desc_phosphorus_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_P, value);
}

static cchem_status_t desc_sulfur_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_S, value);
}

static cchem_status_t desc_chlorine_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cl, value);
}

static cchem_status_t desc_bromine_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Br, value);
}

static cchem_status_t desc_iodine_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_I, value);
}

static cchem_status_t desc_fluorine_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_F, value);
}

/* Additional element counts - all elements in periodic table */
static cchem_status_t desc_helium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_He, value);
}

static cchem_status_t desc_lithium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Li, value);
}

static cchem_status_t desc_beryllium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Be, value);
}

static cchem_status_t desc_boron_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_B, value);
}

static cchem_status_t desc_neon_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ne, value);
}

static cchem_status_t desc_sodium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Na, value);
}

static cchem_status_t desc_magnesium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Mg, value);
}

static cchem_status_t desc_aluminum_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Al, value);
}

static cchem_status_t desc_silicon_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Si, value);
}

static cchem_status_t desc_argon_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ar, value);
}

static cchem_status_t desc_potassium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_K, value);
}

static cchem_status_t desc_calcium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ca, value);
}

static cchem_status_t desc_scandium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Sc, value);
}

static cchem_status_t desc_titanium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ti, value);
}

static cchem_status_t desc_vanadium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_V, value);
}

static cchem_status_t desc_chromium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cr, value);
}

static cchem_status_t desc_manganese_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Mn, value);
}

static cchem_status_t desc_iron_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Fe, value);
}

static cchem_status_t desc_cobalt_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Co, value);
}

static cchem_status_t desc_nickel_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ni, value);
}

static cchem_status_t desc_copper_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cu, value);
}

static cchem_status_t desc_zinc_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Zn, value);
}

static cchem_status_t desc_gallium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ga, value);
}

static cchem_status_t desc_germanium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ge, value);
}

static cchem_status_t desc_arsenic_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_As, value);
}

static cchem_status_t desc_selenium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Se, value);
}

static cchem_status_t desc_krypton_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Kr, value);
}

static cchem_status_t desc_rubidium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Rb, value);
}

static cchem_status_t desc_strontium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Sr, value);
}

static cchem_status_t desc_yttrium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Y, value);
}

static cchem_status_t desc_zirconium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Zr, value);
}

static cchem_status_t desc_niobium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Nb, value);
}

static cchem_status_t desc_molybdenum_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Mo, value);
}

static cchem_status_t desc_technetium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Tc, value);
}

static cchem_status_t desc_ruthenium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ru, value);
}

static cchem_status_t desc_rhodium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Rh, value);
}

static cchem_status_t desc_palladium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Pd, value);
}

static cchem_status_t desc_silver_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ag, value);
}

static cchem_status_t desc_cadmium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cd, value);
}

static cchem_status_t desc_indium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_In, value);
}

static cchem_status_t desc_tin_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Sn, value);
}

static cchem_status_t desc_antimony_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Sb, value);
}

static cchem_status_t desc_tellurium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Te, value);
}

static cchem_status_t desc_xenon_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Xe, value);
}

static cchem_status_t desc_cesium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cs, value);
}

static cchem_status_t desc_barium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ba, value);
}

static cchem_status_t desc_lanthanum_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_La, value);
}

static cchem_status_t desc_cerium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ce, value);
}

static cchem_status_t desc_praseodymium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Pr, value);
}

static cchem_status_t desc_neodymium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Nd, value);
}

static cchem_status_t desc_promethium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Pm, value);
}

static cchem_status_t desc_samarium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Sm, value);
}

static cchem_status_t desc_europium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Eu, value);
}

static cchem_status_t desc_gadolinium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Gd, value);
}

static cchem_status_t desc_terbium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Tb, value);
}

static cchem_status_t desc_dysprosium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Dy, value);
}

static cchem_status_t desc_holmium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ho, value);
}

static cchem_status_t desc_erbium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Er, value);
}

static cchem_status_t desc_thulium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Tm, value);
}

static cchem_status_t desc_ytterbium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Yb, value);
}

static cchem_status_t desc_lutetium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Lu, value);
}

static cchem_status_t desc_hafnium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Hf, value);
}

static cchem_status_t desc_tantalum_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ta, value);
}

static cchem_status_t desc_tungsten_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_W, value);
}

static cchem_status_t desc_rhenium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Re, value);
}

static cchem_status_t desc_osmium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Os, value);
}

static cchem_status_t desc_iridium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ir, value);
}

static cchem_status_t desc_platinum_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Pt, value);
}

static cchem_status_t desc_gold_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Au, value);
}

static cchem_status_t desc_mercury_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Hg, value);
}

static cchem_status_t desc_thallium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Tl, value);
}

static cchem_status_t desc_lead_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Pb, value);
}

static cchem_status_t desc_bismuth_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Bi, value);
}

static cchem_status_t desc_polonium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Po, value);
}

static cchem_status_t desc_astatine_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_At, value);
}

static cchem_status_t desc_radon_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Rn, value);
}

static cchem_status_t desc_francium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Fr, value);
}

static cchem_status_t desc_radium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ra, value);
}

static cchem_status_t desc_actinium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ac, value);
}

static cchem_status_t desc_thorium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Th, value);
}

static cchem_status_t desc_protactinium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Pa, value);
}

static cchem_status_t desc_uranium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_U, value);
}

static cchem_status_t desc_neptunium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Np, value);
}

static cchem_status_t desc_plutonium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Pu, value);
}

static cchem_status_t desc_americium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Am, value);
}

static cchem_status_t desc_curium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cm, value);
}

static cchem_status_t desc_berkelium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Bk, value);
}

static cchem_status_t desc_californium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cf, value);
}

static cchem_status_t desc_einsteinium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Es, value);
}

static cchem_status_t desc_fermium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Fm, value);
}

static cchem_status_t desc_mendelevium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Md, value);
}

static cchem_status_t desc_nobelium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_No, value);
}

static cchem_status_t desc_lawrencium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Lr, value);
}

static cchem_status_t desc_rutherfordium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Rf, value);
}

static cchem_status_t desc_dubnium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Db, value);
}

static cchem_status_t desc_seaborgium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Sg, value);
}

static cchem_status_t desc_bohrium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Bh, value);
}

static cchem_status_t desc_hassium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Hs, value);
}

static cchem_status_t desc_meitnerium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Mt, value);
}

static cchem_status_t desc_darmstadtium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ds, value);
}

static cchem_status_t desc_roentgenium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Rg, value);
}

static cchem_status_t desc_copernicium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Cn, value);
}

static cchem_status_t desc_nihonium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Nh, value);
}

static cchem_status_t desc_flerovium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Fl, value);
}

static cchem_status_t desc_moscovium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Mc, value);
}

static cchem_status_t desc_livermorium_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Lv, value);
}

static cchem_status_t desc_tennessine_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Ts, value);
}

static cchem_status_t desc_oganesson_count(const molecule_t* mol, descriptor_value_t* value) {
    return count_element(mol, ELEM_Og, value);
}

/* ============================================================================
 * Atom Count Descriptors
 * ============================================================================ */

/* Total atom count (heavy atoms only, excludes explicit H) */
static cchem_status_t desc_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_atoms;
    const atom_t* atoms = mol->atoms;

    for (int i = 0; i < n; i++) {
        count += (atoms[i].element != ELEM_H);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Heavy atom count (same as atom count) */
static cchem_status_t desc_heavy_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    return desc_atom_count(mol, value);
}

/* ============================================================================
 * Bond Count Descriptors
 * Ultra-optimized: single pass through bonds array
 * ============================================================================ */

/* Total bond count */
static cchem_status_t desc_total_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->i = mol->num_bonds;
    return CCHEM_OK;
}

/* Single bond count */
static cchem_status_t desc_single_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_bonds;
    const bond_t* bonds = mol->bonds;

    for (int i = 0; i < n; i++) {
        bond_type_t t = bonds[i].type;
        count += (t == BOND_SINGLE || t == BOND_UP || t == BOND_DOWN || t == BOND_RING_SINGLE);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Double bond count */
static cchem_status_t desc_double_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_bonds;
    const bond_t* bonds = mol->bonds;

    for (int i = 0; i < n; i++) {
        count += (bonds[i].type == BOND_DOUBLE);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Triple bond count */
static cchem_status_t desc_triple_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_bonds;
    const bond_t* bonds = mol->bonds;

    for (int i = 0; i < n; i++) {
        count += (bonds[i].type == BOND_TRIPLE);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aromatic bond count */
static cchem_status_t desc_aromatic_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_bonds;
    const bond_t* bonds = mol->bonds;

    for (int i = 0; i < n; i++) {
        count += bonds[i].aromatic;
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Ring Count Descriptors
 * ============================================================================ */

/* Ring count (SSSR) */
static cchem_status_t desc_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    /* Rings should already be computed during parsing */
    value->i = mol->num_rings;
    return CCHEM_OK;
}

/* Aromatic ring count */
static cchem_status_t desc_aromatic_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_rings;
    const ring_t* rings = mol->rings;

    for (int i = 0; i < n; i++) {
        count += rings[i].aromatic;
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aliphatic (non-aromatic) ring count */
static cchem_status_t desc_aliphatic_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_rings;
    const ring_t* rings = mol->rings;

    for (int i = 0; i < n; i++) {
        count += !rings[i].aromatic;
    }

    value->i = count;
    return CCHEM_OK;
}

/* Ring size specific counts */
static cchem_status_t desc_ring3_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        count += (mol->rings[i].size == 3);
    }
    value->i = count;
    return CCHEM_OK;
}

static cchem_status_t desc_ring4_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        count += (mol->rings[i].size == 4);
    }
    value->i = count;
    return CCHEM_OK;
}

static cchem_status_t desc_ring5_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        count += (mol->rings[i].size == 5);
    }
    value->i = count;
    return CCHEM_OK;
}

static cchem_status_t desc_ring6_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        count += (mol->rings[i].size == 6);
    }
    value->i = count;
    return CCHEM_OK;
}

static cchem_status_t desc_large_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        count += (mol->rings[i].size >= 7);
    }
    value->i = count;
    return CCHEM_OK;
}

/* Heterocycle count - rings containing non-C atoms */
static cchem_status_t desc_heterocycle_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        bool has_hetero = false;
        for (int j = 0; j < ring->size && !has_hetero; j++) {
            element_t elem = mol->atoms[ring->atoms[j]].element;
            has_hetero = (elem != ELEM_C && elem != ELEM_H);
        }
        count += has_hetero;
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aromatic heterocycle count */
static cchem_status_t desc_aromatic_heterocycle_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (!ring->aromatic) continue;

        bool has_hetero = false;
        for (int j = 0; j < ring->size && !has_hetero; j++) {
            element_t elem = mol->atoms[ring->atoms[j]].element;
            has_hetero = (elem != ELEM_C && elem != ELEM_H);
        }
        count += has_hetero;
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Heteroatom Descriptors
 * ============================================================================ */

/* Heteroatom count (non-C, non-H) */
static cchem_status_t desc_heteroatom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        count += (e != ELEM_C && e != ELEM_H);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Halogen count (F, Cl, Br, I) */
static cchem_status_t desc_halogen_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        count += (e == ELEM_F || e == ELEM_Cl || e == ELEM_Br || e == ELEM_I);
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Carbon Hybridization Counts
 * ============================================================================ */

/* sp3 Carbon count - carbons with 4 sigma bonds (tetrahedral) */
static cchem_status_t desc_sp3_carbon_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->aromatic) continue;

        /* Check bond orders - sp3 has only single bonds */
        bool has_multiple = false;
        for (int j = 0; j < atom->num_neighbors && !has_multiple; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            has_multiple = (bt == BOND_DOUBLE || bt == BOND_TRIPLE);
        }
        count += !has_multiple;
    }

    value->i = count;
    return CCHEM_OK;
}

/* sp2 Carbon count - carbons with double bond or aromatic */
static cchem_status_t desc_sp2_carbon_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        if (atom->aromatic) {
            count++;
            continue;
        }

        /* Check for double bond (but not triple) */
        bool has_double = false;
        bool has_triple = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt == BOND_DOUBLE) has_double = true;
            if (bt == BOND_TRIPLE) has_triple = true;
        }
        count += (has_double && !has_triple);
    }

    value->i = count;
    return CCHEM_OK;
}

/* sp Carbon count - carbons with triple bond */
static cchem_status_t desc_sp_carbon_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            if (mol->bonds[bond_idx].type == BOND_TRIPLE) {
                count++;
                break;
            }
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aromatic atom count */
static cchem_status_t desc_aromatic_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        count += mol->atoms[i].aromatic;
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Connectivity Descriptors
 * ============================================================================ */

/* Rotatable bond count - single non-ring bonds between heavy atoms */
static cchem_status_t desc_rotatable_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];

        /* Must be single bond, not in ring */
        if (bond->in_ring) continue;
        bond_type_t bt = bond->type;
        if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

        /* Both atoms must be heavy atoms with degree > 1 */
        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        if (a1->element == ELEM_H || a2->element == ELEM_H) continue;
        if (a1->num_neighbors <= 1 || a2->num_neighbors <= 1) continue;

        /* Exclude terminal groups like -CH3, -NH2, -OH */
        /* A rotatable bond should have substituents on both sides */
        int heavy1 = 0, heavy2 = 0;
        for (int j = 0; j < a1->num_neighbors; j++) {
            if (mol->atoms[a1->neighbors[j]].element != ELEM_H) heavy1++;
        }
        for (int j = 0; j < a2->num_neighbors; j++) {
            if (mol->atoms[a2->neighbors[j]].element != ELEM_H) heavy2++;
        }

        count += (heavy1 > 1 && heavy2 > 1);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Terminal atom count - heavy atoms with degree 1 */
static cchem_status_t desc_terminal_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;
        count += (atom->num_neighbors == 1);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Branch point count - atoms with degree >= 3 */
static cchem_status_t desc_branch_point_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;
        count += (atom->num_neighbors >= 3);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Quaternary carbon count - C with 4 non-H neighbors */
static cchem_status_t desc_quaternary_carbon_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        int heavy_neighbors = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) {
                heavy_neighbors++;
            }
        }
        count += (heavy_neighbors == 4);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Ring atom count */
static cchem_status_t desc_ring_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        count += (mol->atoms[i].ring_count > 0);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Non-ring atom count (chain atoms) */
static cchem_status_t desc_chain_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        count += (mol->atoms[i].ring_count == 0);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Ring bond count */
static cchem_status_t desc_ring_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        count += mol->bonds[i].in_ring;
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * H-Bond Descriptors
 * ============================================================================ */

/* H-bond donor count - OH, NH groups */
static cchem_status_t desc_hbond_donor_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        element_t e = atom->element;

        /* O and N with attached H are donors */
        if (e == ELEM_O || e == ELEM_N) {
            /* Count attached hydrogens (implicit + explicit) */
            int h_count = atom->implicit_h_count;
            for (int j = 0; j < atom->num_neighbors; j++) {
                if (mol->atoms[atom->neighbors[j]].element == ELEM_H) {
                    h_count++;
                }
            }
            count += (h_count > 0);
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* H-bond acceptor count - O and N atoms */
static cchem_status_t desc_hbond_acceptor_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        count += (e == ELEM_O || e == ELEM_N);
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Charge Descriptors
 * ============================================================================ */

/* Formal charge count - atoms with non-zero charge */
static cchem_status_t desc_formal_charge_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        count += (mol->atoms[i].charge != 0);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Positive charge count */
static cchem_status_t desc_positive_charge_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        count += (mol->atoms[i].charge > 0);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Negative charge count */
static cchem_status_t desc_negative_charge_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        count += (mol->atoms[i].charge < 0);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Net charge (sum of formal charges) */
static cchem_status_t desc_net_charge(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t charge = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        charge += mol->atoms[i].charge;
    }

    value->i = charge;
    return CCHEM_OK;
}

/* ============================================================================
 * Stereo Descriptors
 * ============================================================================ */

/* Chiral center count */
static cchem_status_t desc_chiral_center_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        chirality_t c = mol->atoms[i].chirality;
        count += (c != CHIRALITY_NONE);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Double bond stereo count (E/Z) */
static cchem_status_t desc_stereo_double_bond_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        count += (mol->bonds[i].stereo != STEREO_NONE);
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Fragment/Unsaturation Descriptors
 * ============================================================================ */

/* Fragment count (disconnected components) */
static cchem_status_t desc_fragment_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    /* If fragments computed, use that; otherwise count as 1 if connected */
    if (mol->num_fragments > 0) {
        value->i = mol->num_fragments;
    } else {
        value->i = (mol->num_atoms > 0) ? 1 : 0;
    }
    return CCHEM_OK;
}

/* Degree of unsaturation (DBE - double bond equivalents) */
/* DBE = (2C + 2 + N - H - X) / 2 where X = halogens */
static cchem_status_t desc_unsaturation_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t c_count = 0, h_count = 0, n_count = 0, x_count = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        element_t e = atom->element;

        switch (e) {
            case ELEM_C:
                c_count++;
                h_count += atom->implicit_h_count;
                break;
            case ELEM_H:
                h_count++;
                break;
            case ELEM_N:
                n_count++;
                h_count += atom->implicit_h_count;
                break;
            case ELEM_F:
            case ELEM_Cl:
            case ELEM_Br:
            case ELEM_I:
                x_count++;
                break;
            default:
                h_count += atom->implicit_h_count;
                break;
        }
    }

    /* DBE formula: (2C + 2 + N - H - X) / 2 */
    value->i = (2 * c_count + 2 + n_count - h_count - x_count) / 2;
    return CCHEM_OK;
}

/* Bridgehead atom count - atoms in multiple rings with degree >= 3 */
static cchem_status_t desc_bridgehead_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;
        /* Bridgehead: in multiple rings and has degree >= 3 */
        count += (atom->ring_count >= 2 && atom->num_neighbors >= 3);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Spiro atom count - atoms shared by exactly 2 rings with degree 4 */
static cchem_status_t desc_spiro_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;
        /* Spiro: in exactly 2 rings, degree 4, sp3 carbon */
        if (atom->ring_count == 2 && atom->num_neighbors == 4 && atom->element == ELEM_C) {
            /* Check that all 4 neighbors are in different positions of the two rings */
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Functional Group Counts
 * ============================================================================ */

/* Helper: count heavy (non-H) neighbors */
static int count_heavy_neighbors(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            count++;
        }
    }
    return count;
}

/* Helper: get total H count (implicit + explicit) */
static int get_total_h(const molecule_t* mol, const atom_t* atom) {
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) {
            h++;
        }
    }
    return h;
}

/* Helper: check if atom has double bond to element */
static bool has_double_bond_to(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (mol->bonds[bond_idx].type == BOND_DOUBLE) {
            if (mol->atoms[atom->neighbors[i]].element == elem) {
                return true;
            }
        }
    }
    return false;
}

/* Helper: check if atom has single bond to element with H */
static bool has_single_bond_to_with_h(const molecule_t* mol, const atom_t* atom, element_t elem) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        bond_type_t bt = mol->bonds[bond_idx].type;
        if (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) {
            const atom_t* neighbor = &mol->atoms[atom->neighbors[i]];
            if (neighbor->element == elem && get_total_h(mol, neighbor) > 0) {
                return true;
            }
        }
    }
    return false;
}

/* Primary amine count - NH2 groups */
static cchem_status_t desc_primary_amine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->aromatic) continue;

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Primary amine: 1 heavy neighbor, 2+ H */
        if (heavy == 1 && h >= 2) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Secondary amine count - R-NH-R groups */
static cchem_status_t desc_secondary_amine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->aromatic) continue;

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Secondary amine: 2 heavy neighbors, 1 H */
        if (heavy == 2 && h == 1) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Tertiary amine count - R3N groups */
static cchem_status_t desc_tertiary_amine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->aromatic) continue;

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Tertiary amine: 3 heavy neighbors, 0 H */
        if (heavy == 3 && h == 0) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aromatic nitrogen count */
static cchem_status_t desc_aromatic_nitrogen_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        count += (atom->element == ELEM_N && atom->aromatic);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Hydroxyl count - OH groups */
static cchem_status_t desc_hydroxyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Hydroxyl: 1 heavy neighbor, has H */
        if (heavy == 1 && h >= 1) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Carbonyl count - C=O groups */
static cchem_status_t desc_carbonyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_DOUBLE) continue;

        element_t e1 = mol->atoms[bond->atom1].element;
        element_t e2 = mol->atoms[bond->atom2].element;

        /* C=O double bond */
        if ((e1 == ELEM_C && e2 == ELEM_O) || (e1 == ELEM_O && e2 == ELEM_C)) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Carboxyl count - COOH groups */
static cchem_status_t desc_carboxyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        /* Check for C with =O and -OH */
        if (has_double_bond_to(mol, atom, ELEM_O) &&
            has_single_bond_to_with_h(mol, atom, ELEM_O)) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Ether count - C-O-C linkages */
static cchem_status_t desc_ether_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;
        if (atom->aromatic) continue;  /* Skip furan-type O */

        int h = get_total_h(mol, atom);
        if (h > 0) continue;  /* Has H, not ether */

        /* Count C neighbors via single bonds */
        int c_neighbors = 0;
        bool has_double = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt == BOND_DOUBLE) {
                has_double = true;
                break;
            }
            if (mol->atoms[atom->neighbors[j]].element == ELEM_C) {
                c_neighbors++;
            }
        }

        /* Ether: 2 C neighbors via single bonds, no double bond */
        if (!has_double && c_neighbors == 2) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Ester count - C(=O)O-C groups */
static cchem_status_t desc_ester_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        /* Check for C with =O and -O-C (no H on that O) */
        bool has_carbonyl_o = has_double_bond_to(mol, atom, ELEM_O);
        if (!has_carbonyl_o) continue;

        /* Check for single bond to O that connects to C (not H) */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt != BOND_SINGLE && bt != BOND_UP && bt != BOND_DOWN) continue;

            const atom_t* neighbor = &mol->atoms[atom->neighbors[j]];
            if (neighbor->element != ELEM_O) continue;

            int o_h = get_total_h(mol, neighbor);
            int o_heavy = count_heavy_neighbors(mol, neighbor);

            /* O with 2 heavy neighbors (C=O carbon + other C), no H */
            if (o_h == 0 && o_heavy == 2) {
                count++;
                break;
            }
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Amide count - C(=O)N groups */
static cchem_status_t desc_amide_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        /* Check for C with =O and -N */
        bool has_carbonyl_o = has_double_bond_to(mol, atom, ELEM_O);
        if (!has_carbonyl_o) continue;

        /* Check for single bond to N */
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt == BOND_DOUBLE || bt == BOND_TRIPLE) continue;

            if (mol->atoms[atom->neighbors[j]].element == ELEM_N) {
                count++;
                break;
            }
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Nitrile count - C#N groups */
static cchem_status_t desc_nitrile_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_TRIPLE) continue;

        element_t e1 = mol->atoms[bond->atom1].element;
        element_t e2 = mol->atoms[bond->atom2].element;

        /* C#N triple bond */
        if ((e1 == ELEM_C && e2 == ELEM_N) || (e1 == ELEM_N && e2 == ELEM_C)) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Nitro count - NO2 groups */
static cchem_status_t desc_nitro_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;

        /* Nitro: N with positive charge, 2 O neighbors */
        if (atom->charge > 0) {
            int o_count = 0;
            for (int j = 0; j < atom->num_neighbors; j++) {
                if (mol->atoms[atom->neighbors[j]].element == ELEM_O) {
                    o_count++;
                }
            }
            if (o_count == 2) {
                count++;
            }
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Thiol count - SH groups */
static cchem_status_t desc_thiol_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Thiol: 1 heavy neighbor, has H */
        if (heavy == 1 && h >= 1) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Sulfone count - SO2 groups */
static cchem_status_t desc_sulfone_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        /* Count double bonds to O */
        int double_o = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            if (mol->bonds[bond_idx].type == BOND_DOUBLE &&
                mol->atoms[atom->neighbors[j]].element == ELEM_O) {
                double_o++;
            }
        }

        /* Sulfone: 2 double bonds to O */
        if (double_o == 2) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Methyl count - CH3 groups */
static cchem_status_t desc_methyl_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Methyl: 1 heavy neighbor, 3 H */
        if (heavy == 1 && h == 3) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Methylene count - CH2 groups in chains */
static cchem_status_t desc_methylene_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->ring_count > 0) continue;  /* Not in ring */

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Methylene: 2 heavy neighbors, 2 H, not in ring */
        if (heavy == 2 && h == 2) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Methine count - CH groups */
static cchem_status_t desc_methine_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (atom->aromatic) continue;

        int heavy = count_heavy_neighbors(mol, atom);
        int h = get_total_h(mol, atom);

        /* Methine: 3 heavy neighbors, 1 H */
        if (heavy == 3 && h == 1) {
            count++;
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aliphatic carbon count - non-aromatic carbons */
static cchem_status_t desc_aliphatic_carbon_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        count += (atom->element == ELEM_C && !atom->aromatic);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aromatic carbon count */
static cchem_status_t desc_aromatic_carbon_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        count += (atom->element == ELEM_C && atom->aromatic);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aromatic 5-ring count */
static cchem_status_t desc_aromatic5_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        count += (ring->aromatic && ring->size == 5);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Aromatic 6-ring count */
static cchem_status_t desc_aromatic6_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        count += (ring->aromatic && ring->size == 6);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Saturated ring count - rings with no double bonds */
static cchem_status_t desc_saturated_ring_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        if (ring->aromatic) continue;

        /* Check if any bond in ring is double */
        bool has_double = false;
        for (int j = 0; j < ring->size && !has_double; j++) {
            if (mol->bonds[ring->bonds[j]].type == BOND_DOUBLE) {
                has_double = true;
            }
        }
        count += !has_double;
    }

    value->i = count;
    return CCHEM_OK;
}

/* Carbocycle count - rings containing only C atoms */
static cchem_status_t desc_carbocycle_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        bool all_carbon = true;
        for (int j = 0; j < ring->size && all_carbon; j++) {
            if (mol->atoms[ring->atoms[j]].element != ELEM_C) {
                all_carbon = false;
            }
        }
        count += all_carbon;
    }

    value->i = count;
    return CCHEM_OK;
}

/* Degree 2 atom count - atoms with exactly 2 heavy neighbors */
static cchem_status_t desc_degree2_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;
        int heavy = count_heavy_neighbors(mol, atom);
        count += (heavy == 2);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Degree 3 atom count - atoms with exactly 3 heavy neighbors */
static cchem_status_t desc_degree3_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;
        int heavy = count_heavy_neighbors(mol, atom);
        count += (heavy == 3);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Degree 4 atom count - atoms with exactly 4 heavy neighbors */
static cchem_status_t desc_degree4_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;
        int heavy = count_heavy_neighbors(mol, atom);
        count += (heavy == 4);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Polar atom count - O, N, S atoms */
static cchem_status_t desc_polar_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        element_t e = mol->atoms[i].element;
        count += (e == ELEM_O || e == ELEM_N || e == ELEM_S);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Acidic proton count - H on O or N */
static cchem_status_t desc_acidic_proton_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        element_t e = atom->element;

        if (e == ELEM_O || e == ELEM_N) {
            count += get_total_h(mol, atom);
        }
    }

    value->i = count;
    return CCHEM_OK;
}

/* Explicit H count - only explicit H atoms */
static cchem_status_t desc_explicit_h_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        count += (mol->atoms[i].element == ELEM_H);
    }

    value->i = count;
    return CCHEM_OK;
}

/* Conjugated atom count - atoms in alternating single-double pattern */
static cchem_status_t desc_conjugated_atom_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    /* Atoms are conjugated if aromatic OR have adjacent single and double bonds */
    int64_t count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        if (atom->aromatic) {
            count++;
            continue;
        }

        /* Check for both single and double bonds */
        bool has_single = false, has_double = false;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bond_idx = atom->neighbor_bonds[j];
            bond_type_t bt = mol->bonds[bond_idx].type;
            if (bt == BOND_SINGLE || bt == BOND_UP || bt == BOND_DOWN) has_single = true;
            if (bt == BOND_DOUBLE) has_double = true;
        }
        count += (has_single && has_double);
    }

    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

/* Helper macro for registering count descriptors */
#define REGISTER_COUNT_DESC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_COUNTS; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_counts(void) {
    /* Element counts */
    REGISTER_COUNT_DESC("CarbonCount", "Number of carbon atoms", desc_carbon_count);
    REGISTER_COUNT_DESC("HydrogenCount", "Number of hydrogen atoms (explicit + implicit)", desc_hydrogen_count);
    REGISTER_COUNT_DESC("OxygenCount", "Number of oxygen atoms", desc_oxygen_count);
    REGISTER_COUNT_DESC("NitrogenCount", "Number of nitrogen atoms", desc_nitrogen_count);
    REGISTER_COUNT_DESC("PhosphorusCount", "Number of phosphorus atoms", desc_phosphorus_count);
    REGISTER_COUNT_DESC("SulfurCount", "Number of sulfur atoms", desc_sulfur_count);
    REGISTER_COUNT_DESC("ChlorineCount", "Number of chlorine atoms", desc_chlorine_count);
    REGISTER_COUNT_DESC("BromineCount", "Number of bromine atoms", desc_bromine_count);
    REGISTER_COUNT_DESC("IodineCount", "Number of iodine atoms", desc_iodine_count);
    REGISTER_COUNT_DESC("FluorineCount", "Number of fluorine atoms", desc_fluorine_count);

    /* Additional element counts - Period 1-2 */
    REGISTER_COUNT_DESC("HeliumCount", "Number of helium atoms", desc_helium_count);
    REGISTER_COUNT_DESC("LithiumCount", "Number of lithium atoms", desc_lithium_count);
    REGISTER_COUNT_DESC("BerylliumCount", "Number of beryllium atoms", desc_beryllium_count);
    REGISTER_COUNT_DESC("BoronCount", "Number of boron atoms", desc_boron_count);
    REGISTER_COUNT_DESC("NeonCount", "Number of neon atoms", desc_neon_count);

    /* Period 3 */
    REGISTER_COUNT_DESC("SodiumCount", "Number of sodium atoms", desc_sodium_count);
    REGISTER_COUNT_DESC("MagnesiumCount", "Number of magnesium atoms", desc_magnesium_count);
    REGISTER_COUNT_DESC("AluminumCount", "Number of aluminum atoms", desc_aluminum_count);
    REGISTER_COUNT_DESC("SiliconCount", "Number of silicon atoms", desc_silicon_count);
    REGISTER_COUNT_DESC("ArgonCount", "Number of argon atoms", desc_argon_count);

    /* Period 4 */
    REGISTER_COUNT_DESC("PotassiumCount", "Number of potassium atoms", desc_potassium_count);
    REGISTER_COUNT_DESC("CalciumCount", "Number of calcium atoms", desc_calcium_count);
    REGISTER_COUNT_DESC("ScandiumCount", "Number of scandium atoms", desc_scandium_count);
    REGISTER_COUNT_DESC("TitaniumCount", "Number of titanium atoms", desc_titanium_count);
    REGISTER_COUNT_DESC("VanadiumCount", "Number of vanadium atoms", desc_vanadium_count);
    REGISTER_COUNT_DESC("ChromiumCount", "Number of chromium atoms", desc_chromium_count);
    REGISTER_COUNT_DESC("ManganeseCount", "Number of manganese atoms", desc_manganese_count);
    REGISTER_COUNT_DESC("IronCount", "Number of iron atoms", desc_iron_count);
    REGISTER_COUNT_DESC("CobaltCount", "Number of cobalt atoms", desc_cobalt_count);
    REGISTER_COUNT_DESC("NickelCount", "Number of nickel atoms", desc_nickel_count);
    REGISTER_COUNT_DESC("CopperCount", "Number of copper atoms", desc_copper_count);
    REGISTER_COUNT_DESC("ZincCount", "Number of zinc atoms", desc_zinc_count);
    REGISTER_COUNT_DESC("GalliumCount", "Number of gallium atoms", desc_gallium_count);
    REGISTER_COUNT_DESC("GermaniumCount", "Number of germanium atoms", desc_germanium_count);
    REGISTER_COUNT_DESC("ArsenicCount", "Number of arsenic atoms", desc_arsenic_count);
    REGISTER_COUNT_DESC("SeleniumCount", "Number of selenium atoms", desc_selenium_count);
    REGISTER_COUNT_DESC("KryptonCount", "Number of krypton atoms", desc_krypton_count);

    /* Period 5 */
    REGISTER_COUNT_DESC("RubidiumCount", "Number of rubidium atoms", desc_rubidium_count);
    REGISTER_COUNT_DESC("StrontiumCount", "Number of strontium atoms", desc_strontium_count);
    REGISTER_COUNT_DESC("YttriumCount", "Number of yttrium atoms", desc_yttrium_count);
    REGISTER_COUNT_DESC("ZirconiumCount", "Number of zirconium atoms", desc_zirconium_count);
    REGISTER_COUNT_DESC("NiobiumCount", "Number of niobium atoms", desc_niobium_count);
    REGISTER_COUNT_DESC("MolybdenumCount", "Number of molybdenum atoms", desc_molybdenum_count);
    REGISTER_COUNT_DESC("TechnetiumCount", "Number of technetium atoms", desc_technetium_count);
    REGISTER_COUNT_DESC("RutheniumCount", "Number of ruthenium atoms", desc_ruthenium_count);
    REGISTER_COUNT_DESC("RhodiumCount", "Number of rhodium atoms", desc_rhodium_count);
    REGISTER_COUNT_DESC("PalladiumCount", "Number of palladium atoms", desc_palladium_count);
    REGISTER_COUNT_DESC("SilverCount", "Number of silver atoms", desc_silver_count);
    REGISTER_COUNT_DESC("CadmiumCount", "Number of cadmium atoms", desc_cadmium_count);
    REGISTER_COUNT_DESC("IndiumCount", "Number of indium atoms", desc_indium_count);
    REGISTER_COUNT_DESC("TinCount", "Number of tin atoms", desc_tin_count);
    REGISTER_COUNT_DESC("AntimonyCount", "Number of antimony atoms", desc_antimony_count);
    REGISTER_COUNT_DESC("TelluriumCount", "Number of tellurium atoms", desc_tellurium_count);
    REGISTER_COUNT_DESC("XenonCount", "Number of xenon atoms", desc_xenon_count);

    /* Period 6 */
    REGISTER_COUNT_DESC("CesiumCount", "Number of cesium atoms", desc_cesium_count);
    REGISTER_COUNT_DESC("BariumCount", "Number of barium atoms", desc_barium_count);
    REGISTER_COUNT_DESC("LanthanumCount", "Number of lanthanum atoms", desc_lanthanum_count);
    REGISTER_COUNT_DESC("CeriumCount", "Number of cerium atoms", desc_cerium_count);
    REGISTER_COUNT_DESC("PraseodymiumCount", "Number of praseodymium atoms", desc_praseodymium_count);
    REGISTER_COUNT_DESC("NeodymiumCount", "Number of neodymium atoms", desc_neodymium_count);
    REGISTER_COUNT_DESC("PromethiumCount", "Number of promethium atoms", desc_promethium_count);
    REGISTER_COUNT_DESC("SamariumCount", "Number of samarium atoms", desc_samarium_count);
    REGISTER_COUNT_DESC("EuropiumCount", "Number of europium atoms", desc_europium_count);
    REGISTER_COUNT_DESC("GadoliniumCount", "Number of gadolinium atoms", desc_gadolinium_count);
    REGISTER_COUNT_DESC("TerbiumCount", "Number of terbium atoms", desc_terbium_count);
    REGISTER_COUNT_DESC("DysprosiumCount", "Number of dysprosium atoms", desc_dysprosium_count);
    REGISTER_COUNT_DESC("HolmiumCount", "Number of holmium atoms", desc_holmium_count);
    REGISTER_COUNT_DESC("ErbiumCount", "Number of erbium atoms", desc_erbium_count);
    REGISTER_COUNT_DESC("ThuliumCount", "Number of thulium atoms", desc_thulium_count);
    REGISTER_COUNT_DESC("YtterbiumCount", "Number of ytterbium atoms", desc_ytterbium_count);
    REGISTER_COUNT_DESC("LutetiumCount", "Number of lutetium atoms", desc_lutetium_count);
    REGISTER_COUNT_DESC("HafniumCount", "Number of hafnium atoms", desc_hafnium_count);
    REGISTER_COUNT_DESC("TantalumCount", "Number of tantalum atoms", desc_tantalum_count);
    REGISTER_COUNT_DESC("TungstenCount", "Number of tungsten atoms", desc_tungsten_count);
    REGISTER_COUNT_DESC("RheniumCount", "Number of rhenium atoms", desc_rhenium_count);
    REGISTER_COUNT_DESC("OsmiumCount", "Number of osmium atoms", desc_osmium_count);
    REGISTER_COUNT_DESC("IridiumCount", "Number of iridium atoms", desc_iridium_count);
    REGISTER_COUNT_DESC("PlatinumCount", "Number of platinum atoms", desc_platinum_count);
    REGISTER_COUNT_DESC("GoldCount", "Number of gold atoms", desc_gold_count);
    REGISTER_COUNT_DESC("MercuryCount", "Number of mercury atoms", desc_mercury_count);
    REGISTER_COUNT_DESC("ThalliumCount", "Number of thallium atoms", desc_thallium_count);
    REGISTER_COUNT_DESC("LeadCount", "Number of lead atoms", desc_lead_count);
    REGISTER_COUNT_DESC("BismuthCount", "Number of bismuth atoms", desc_bismuth_count);
    REGISTER_COUNT_DESC("PoloniumCount", "Number of polonium atoms", desc_polonium_count);
    REGISTER_COUNT_DESC("AstatineCount", "Number of astatine atoms", desc_astatine_count);
    REGISTER_COUNT_DESC("RadonCount", "Number of radon atoms", desc_radon_count);

    /* Period 7 */
    REGISTER_COUNT_DESC("FranciumCount", "Number of francium atoms", desc_francium_count);
    REGISTER_COUNT_DESC("RadiumCount", "Number of radium atoms", desc_radium_count);
    REGISTER_COUNT_DESC("ActiniumCount", "Number of actinium atoms", desc_actinium_count);
    REGISTER_COUNT_DESC("ThoriumCount", "Number of thorium atoms", desc_thorium_count);
    REGISTER_COUNT_DESC("ProtactiniumCount", "Number of protactinium atoms", desc_protactinium_count);
    REGISTER_COUNT_DESC("UraniumCount", "Number of uranium atoms", desc_uranium_count);
    REGISTER_COUNT_DESC("NeptuniumCount", "Number of neptunium atoms", desc_neptunium_count);
    REGISTER_COUNT_DESC("PlutoniumCount", "Number of plutonium atoms", desc_plutonium_count);
    REGISTER_COUNT_DESC("AmericiumCount", "Number of americium atoms", desc_americium_count);
    REGISTER_COUNT_DESC("CuriumCount", "Number of curium atoms", desc_curium_count);
    REGISTER_COUNT_DESC("BerkeliumCount", "Number of berkelium atoms", desc_berkelium_count);
    REGISTER_COUNT_DESC("CaliforniumCount", "Number of californium atoms", desc_californium_count);
    REGISTER_COUNT_DESC("EinsteiniumCount", "Number of einsteinium atoms", desc_einsteinium_count);
    REGISTER_COUNT_DESC("FermiumCount", "Number of fermium atoms", desc_fermium_count);
    REGISTER_COUNT_DESC("MendeleviumCount", "Number of mendelevium atoms", desc_mendelevium_count);
    REGISTER_COUNT_DESC("NobeliumCount", "Number of nobelium atoms", desc_nobelium_count);
    REGISTER_COUNT_DESC("LawrenciumCount", "Number of lawrencium atoms", desc_lawrencium_count);
    REGISTER_COUNT_DESC("RutherfordiumCount", "Number of rutherfordium atoms", desc_rutherfordium_count);
    REGISTER_COUNT_DESC("DubniumCount", "Number of dubnium atoms", desc_dubnium_count);
    REGISTER_COUNT_DESC("SeaborgiumCount", "Number of seaborgium atoms", desc_seaborgium_count);
    REGISTER_COUNT_DESC("BohriumCount", "Number of bohrium atoms", desc_bohrium_count);
    REGISTER_COUNT_DESC("HassiumCount", "Number of hassium atoms", desc_hassium_count);
    REGISTER_COUNT_DESC("MeitneriumCount", "Number of meitnerium atoms", desc_meitnerium_count);
    REGISTER_COUNT_DESC("DarmstadtiumCount", "Number of darmstadtium atoms", desc_darmstadtium_count);
    REGISTER_COUNT_DESC("RoentgeniumCount", "Number of roentgenium atoms", desc_roentgenium_count);
    REGISTER_COUNT_DESC("CoperniciumCount", "Number of copernicium atoms", desc_copernicium_count);
    REGISTER_COUNT_DESC("NihoniumCount", "Number of nihonium atoms", desc_nihonium_count);
    REGISTER_COUNT_DESC("FleroviumCount", "Number of flerovium atoms", desc_flerovium_count);
    REGISTER_COUNT_DESC("MoscoviumCount", "Number of moscovium atoms", desc_moscovium_count);
    REGISTER_COUNT_DESC("LivermoriumCount", "Number of livermorium atoms", desc_livermorium_count);
    REGISTER_COUNT_DESC("TennessineCount", "Number of tennessine atoms", desc_tennessine_count);
    REGISTER_COUNT_DESC("OganessonCount", "Number of oganesson atoms", desc_oganesson_count);

    /* Heteroatom counts */
    REGISTER_COUNT_DESC("HeteroatomCount", "Non-carbon, non-hydrogen atoms", desc_heteroatom_count);
    REGISTER_COUNT_DESC("HalogenCount", "Number of halogen atoms (F, Cl, Br, I)", desc_halogen_count);

    /* Atom counts */
    REGISTER_COUNT_DESC("AtomCount", "Total heavy atoms (non-hydrogen)", desc_atom_count);
    REGISTER_COUNT_DESC("HeavyAtomCount", "Total heavy atoms (non-hydrogen)", desc_heavy_atom_count);
    REGISTER_COUNT_DESC("AromaticAtomCount", "Number of aromatic atoms", desc_aromatic_atom_count);

    /* Carbon hybridization */
    REGISTER_COUNT_DESC("SP3CarbonCount", "sp3 hybridized carbons (tetrahedral)", desc_sp3_carbon_count);
    REGISTER_COUNT_DESC("SP2CarbonCount", "sp2 hybridized carbons (trigonal planar)", desc_sp2_carbon_count);
    REGISTER_COUNT_DESC("SPCarbonCount", "sp hybridized carbons (linear)", desc_sp_carbon_count);

    /* Bond counts */
    REGISTER_COUNT_DESC("TotalBondCount", "Total number of bonds", desc_total_bond_count);
    REGISTER_COUNT_DESC("SingleBondCount", "Number of single bonds", desc_single_bond_count);
    REGISTER_COUNT_DESC("DoubleBondCount", "Number of double bonds", desc_double_bond_count);
    REGISTER_COUNT_DESC("TripleBondCount", "Number of triple bonds", desc_triple_bond_count);
    REGISTER_COUNT_DESC("AromaticBondCount", "Number of aromatic bonds", desc_aromatic_bond_count);
    REGISTER_COUNT_DESC("RingBondCount", "Number of bonds in rings", desc_ring_bond_count);
    REGISTER_COUNT_DESC("RotatableBondCount", "Number of rotatable bonds", desc_rotatable_bond_count);

    /* Ring counts */
    REGISTER_COUNT_DESC("RingCount", "Number of rings (SSSR)", desc_ring_count);
    REGISTER_COUNT_DESC("AromaticRingCount", "Number of aromatic rings", desc_aromatic_ring_count);
    REGISTER_COUNT_DESC("AliphaticRingCount", "Number of non-aromatic rings", desc_aliphatic_ring_count);
    REGISTER_COUNT_DESC("HeterocycleCount", "Rings containing heteroatoms", desc_heterocycle_count);
    REGISTER_COUNT_DESC("AromaticHeterocycleCount", "Aromatic rings with heteroatoms", desc_aromatic_heterocycle_count);

    /* Ring size counts */
    REGISTER_COUNT_DESC("Ring3Count", "Number of 3-membered rings", desc_ring3_count);
    REGISTER_COUNT_DESC("Ring4Count", "Number of 4-membered rings", desc_ring4_count);
    REGISTER_COUNT_DESC("Ring5Count", "Number of 5-membered rings", desc_ring5_count);
    REGISTER_COUNT_DESC("Ring6Count", "Number of 6-membered rings", desc_ring6_count);
    REGISTER_COUNT_DESC("LargeRingCount", "Number of 7+ membered rings", desc_large_ring_count);

    /* Connectivity */
    REGISTER_COUNT_DESC("RingAtomCount", "Heavy atoms in rings", desc_ring_atom_count);
    REGISTER_COUNT_DESC("ChainAtomCount", "Heavy atoms not in rings", desc_chain_atom_count);
    REGISTER_COUNT_DESC("TerminalAtomCount", "Heavy atoms with degree 1", desc_terminal_atom_count);
    REGISTER_COUNT_DESC("BranchPointCount", "Atoms with degree >= 3", desc_branch_point_count);
    REGISTER_COUNT_DESC("QuaternaryCarbonCount", "Carbons with 4 non-H neighbors", desc_quaternary_carbon_count);
    REGISTER_COUNT_DESC("BridgeheadAtomCount", "Atoms in 2+ rings with degree >= 3", desc_bridgehead_atom_count);
    REGISTER_COUNT_DESC("SpiroAtomCount", "Spiro junction atoms", desc_spiro_atom_count);

    /* H-bonding */
    REGISTER_COUNT_DESC("HBondDonorCount", "H-bond donors (OH, NH)", desc_hbond_donor_count);
    REGISTER_COUNT_DESC("HBondAcceptorCount", "H-bond acceptors (O, N)", desc_hbond_acceptor_count);

    /* Charge */
    REGISTER_COUNT_DESC("FormalChargeCount", "Atoms with formal charge", desc_formal_charge_count);
    REGISTER_COUNT_DESC("PositiveChargeCount", "Atoms with positive charge", desc_positive_charge_count);
    REGISTER_COUNT_DESC("NegativeChargeCount", "Atoms with negative charge", desc_negative_charge_count);
    REGISTER_COUNT_DESC("NetCharge", "Sum of formal charges", desc_net_charge);

    /* Stereochemistry */
    REGISTER_COUNT_DESC("ChiralCenterCount", "Number of chiral centers", desc_chiral_center_count);
    REGISTER_COUNT_DESC("StereoDoubleBondCount", "Double bonds with E/Z stereo", desc_stereo_double_bond_count);

    /* Fragments and unsaturation */
    REGISTER_COUNT_DESC("FragmentCount", "Number of disconnected fragments", desc_fragment_count);
    REGISTER_COUNT_DESC("UnsaturationCount", "Degree of unsaturation (DBE)", desc_unsaturation_count);

    /* Amine classification */
    REGISTER_COUNT_DESC("PrimaryAmineCount", "Primary amine groups (-NH2)", desc_primary_amine_count);
    REGISTER_COUNT_DESC("SecondaryAmineCount", "Secondary amine groups (-NH-)", desc_secondary_amine_count);
    REGISTER_COUNT_DESC("TertiaryAmineCount", "Tertiary amine groups (-N<)", desc_tertiary_amine_count);
    REGISTER_COUNT_DESC("AromaticNitrogenCount", "Nitrogen atoms in aromatic rings", desc_aromatic_nitrogen_count);

    /* Functional groups */
    REGISTER_COUNT_DESC("HydroxylCount", "Hydroxyl groups (-OH)", desc_hydroxyl_count);
    REGISTER_COUNT_DESC("CarbonylCount", "Carbonyl groups (C=O)", desc_carbonyl_count);
    REGISTER_COUNT_DESC("CarboxylCount", "Carboxyl groups (-COOH)", desc_carboxyl_count);
    REGISTER_COUNT_DESC("EtherCount", "Ether linkages (C-O-C)", desc_ether_count);
    REGISTER_COUNT_DESC("EsterCount", "Ester groups (-COO-)", desc_ester_count);
    REGISTER_COUNT_DESC("AmideCount", "Amide groups (-CONH-)", desc_amide_count);
    REGISTER_COUNT_DESC("NitrileCount", "Nitrile groups (-C#N)", desc_nitrile_count);
    REGISTER_COUNT_DESC("NitroCount", "Nitro groups (-NO2)", desc_nitro_count);
    REGISTER_COUNT_DESC("ThiolCount", "Thiol groups (-SH)", desc_thiol_count);
    REGISTER_COUNT_DESC("SulfoneCount", "Sulfone groups (-SO2-)", desc_sulfone_count);

    /* Carbon types */
    REGISTER_COUNT_DESC("MethylCount", "Methyl groups (-CH3)", desc_methyl_count);
    REGISTER_COUNT_DESC("MethyleneCount", "Methylene groups in chains (-CH2-)", desc_methylene_count);
    REGISTER_COUNT_DESC("MethineCount", "Methine groups (>CH-)", desc_methine_count);
    REGISTER_COUNT_DESC("AliphaticCarbonCount", "Non-aromatic carbon atoms", desc_aliphatic_carbon_count);
    REGISTER_COUNT_DESC("AromaticCarbonCount", "Aromatic carbon atoms", desc_aromatic_carbon_count);

    /* Ring classification */
    REGISTER_COUNT_DESC("Aromatic5RingCount", "5-membered aromatic rings", desc_aromatic5_ring_count);
    REGISTER_COUNT_DESC("Aromatic6RingCount", "6-membered aromatic rings", desc_aromatic6_ring_count);
    REGISTER_COUNT_DESC("SaturatedRingCount", "Fully saturated rings", desc_saturated_ring_count);
    REGISTER_COUNT_DESC("CarbocycleCount", "Carbon-only rings", desc_carbocycle_count);

    /* Degree distribution */
    REGISTER_COUNT_DESC("Degree2AtomCount", "Atoms with 2 heavy neighbors", desc_degree2_atom_count);
    REGISTER_COUNT_DESC("Degree3AtomCount", "Atoms with 3 heavy neighbors", desc_degree3_atom_count);
    REGISTER_COUNT_DESC("Degree4AtomCount", "Atoms with 4 heavy neighbors", desc_degree4_atom_count);

    /* Polarity and acidity */
    REGISTER_COUNT_DESC("PolarAtomCount", "Polar atoms (O, N, S)", desc_polar_atom_count);
    REGISTER_COUNT_DESC("AcidicProtonCount", "Protons on O or N", desc_acidic_proton_count);
    REGISTER_COUNT_DESC("ExplicitHCount", "Explicit hydrogen atoms", desc_explicit_h_count);
    REGISTER_COUNT_DESC("ConjugatedAtomCount", "Atoms in conjugated systems", desc_conjugated_atom_count);
}
