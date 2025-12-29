/**
 * @file estate.c
 * @brief E-State atom type count descriptors
 *
 * Counts atoms by electrotopological state type.
 * Types use Hall-Kier notation:
 * - First letter(s): hybridization (s=sp3, d=sp2, t=sp, a=aromatic)
 * - Following letters: bonded atom types
 *
 * Reference: Hall & Kier, J. Chem. Inf. Comput. Sci. 1995
 *
 * All O(n) complexity.
 */

#include <string.h>
#include <stdbool.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/* Check if atom is sp3 (no double/triple/aromatic bonds) */
static bool is_sp3(const molecule_t* mol, const atom_t* atom) {
    if (atom->aromatic) return false;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bi = atom->neighbor_bonds[i];
        if (mol->bonds[bi].type == BOND_DOUBLE ||
            mol->bonds[bi].type == BOND_TRIPLE ||
            mol->bonds[bi].type == BOND_AROMATIC) {
            return false;
        }
    }
    return true;
}

/* Check if atom is sp2 (has double bond, not aromatic) */
static bool is_sp2(const molecule_t* mol, const atom_t* atom) {
    if (atom->aromatic) return false;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bi = atom->neighbor_bonds[i];
        if (mol->bonds[bi].type == BOND_DOUBLE) return true;
    }
    return false;
}

/* Check if atom is sp (has triple bond) */
static bool is_sp(const molecule_t* mol, const atom_t* atom) {
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bi = atom->neighbor_bonds[i];
        if (mol->bonds[bi].type == BOND_TRIPLE) return true;
    }
    return false;
}

/* Count heavy neighbors */
static int heavy_count(const molecule_t* mol, const atom_t* atom) {
    int count = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) count++;
    }
    return count;
}

/* Total H count */
static int h_count(const molecule_t* mol, const atom_t* atom) {
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h++;
    }
    return h;
}

/* ============================================================================
 * Carbon E-State Types
 * ============================================================================ */

/* ssCH3: sp3 CH3 (methyl) */
static cchem_status_t estate_ssCH3(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 1 && h_count(mol, atom) == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ssCH2: sp3 CH2 (methylene) */
static cchem_status_t estate_ssCH2(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 2 && h_count(mol, atom) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* sssCH: sp3 CH (methine) */
static cchem_status_t estate_sssCH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 3 && h_count(mol, atom) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ssssC: sp3 quaternary carbon */
static cchem_status_t estate_ssssC(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 4 && h_count(mol, atom) == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dCH2: sp2 =CH2 (terminal alkene) */
static cchem_status_t estate_dCH2(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp2(mol, atom)) continue;
        if (heavy_count(mol, atom) == 1 && h_count(mol, atom) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dsCH: sp2 =CH- (internal alkene with H) */
static cchem_status_t estate_dsCH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp2(mol, atom)) continue;
        if (heavy_count(mol, atom) == 2 && h_count(mol, atom) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dssC: sp2 =C< (substituted alkene) */
static cchem_status_t estate_dssC(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp2(mol, atom)) continue;
        if (heavy_count(mol, atom) == 3 && h_count(mol, atom) == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* aaCH: aromatic CH */
static cchem_status_t estate_aaCH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!atom->aromatic) continue;
        if (h_count(mol, atom) >= 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* aasC: aromatic C with substituent */
static cchem_status_t estate_aasC(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!atom->aromatic) continue;
        if (heavy_count(mol, atom) >= 3 && h_count(mol, atom) == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* tCH: sp CH (alkyne) */
static cchem_status_t estate_tCH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp(mol, atom)) continue;
        if (h_count(mol, atom) >= 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* tsC: sp C (internal alkyne) */
static cchem_status_t estate_tsC(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_C) continue;
        if (!is_sp(mol, atom)) continue;
        if (h_count(mol, atom) == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Nitrogen E-State Types
 * ============================================================================ */

/* sNH3: sp3 -NH3+ (primary amine, charged) */
static cchem_status_t estate_sNH3(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (atom->charge > 0 && h_count(mol, atom) == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* sNH2: sp3 -NH2 (primary amine) */
static cchem_status_t estate_sNH2(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 1 && h_count(mol, atom) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ssNH: sp3 >NH (secondary amine) */
static cchem_status_t estate_ssNH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 2 && h_count(mol, atom) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* sssN: sp3 >N- (tertiary amine) */
static cchem_status_t estate_sssN(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 3 && h_count(mol, atom) == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dNH: sp2 =NH (imine with H) */
static cchem_status_t estate_dNH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (!is_sp2(mol, atom)) continue;
        if (h_count(mol, atom) >= 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dsN: sp2 =N- (imine) */
static cchem_status_t estate_dsN(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (!is_sp2(mol, atom)) continue;
        if (h_count(mol, atom) == 0 && heavy_count(mol, atom) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* aaNH: aromatic NH (pyrrole-type) */
static cchem_status_t estate_aaNH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (!atom->aromatic) continue;
        if (h_count(mol, atom) >= 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* aaN: aromatic N (pyridine-type) */
static cchem_status_t estate_aaN(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (!atom->aromatic) continue;
        if (h_count(mol, atom) == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* tN: sp N (nitrile) */
static cchem_status_t estate_tN(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_N) continue;
        if (is_sp(mol, atom)) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Oxygen E-State Types
 * ============================================================================ */

/* sOH: sp3 -OH (hydroxyl) */
static cchem_status_t estate_sOH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;
        if (heavy_count(mol, atom) == 1 && h_count(mol, atom) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ssO: sp3 -O- (ether) */
static cchem_status_t estate_ssO(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;
        if (!is_sp3(mol, atom)) continue;
        if (heavy_count(mol, atom) == 2 && h_count(mol, atom) == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dO: sp2 =O (carbonyl) */
static cchem_status_t estate_dO(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;
        if (is_sp2(mol, atom) && heavy_count(mol, atom) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* aaO: aromatic O (furan-type) */
static cchem_status_t estate_aaO(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_O) continue;
        if (atom->aromatic) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Sulfur E-State Types
 * ============================================================================ */

/* sSH: sp3 -SH (thiol) */
static cchem_status_t estate_sSH(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;
        if (heavy_count(mol, atom) == 1 && h_count(mol, atom) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ssS: sp3 -S- (thioether) */
static cchem_status_t estate_ssS(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;
        if (is_sp3(mol, atom) && heavy_count(mol, atom) == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dS: sp2 =S (thione) */
static cchem_status_t estate_dS(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;
        if (is_sp2(mol, atom) && heavy_count(mol, atom) == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* aaS: aromatic S (thiophene-type) */
static cchem_status_t estate_aaS(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;
        if (atom->aromatic) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* dssS: sulfoxide S */
static cchem_status_t estate_dssS(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        /* Count double bonds to O */
        int dbl_o = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bi = atom->neighbor_bonds[j];
            if (mol->bonds[bi].type == BOND_DOUBLE &&
                mol->atoms[atom->neighbors[j]].element == ELEM_O) {
                dbl_o++;
            }
        }
        if (dbl_o == 1 && heavy_count(mol, atom) == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ddssS: sulfone S */
static cchem_status_t estate_ddssS(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element != ELEM_S) continue;

        int dbl_o = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            int bi = atom->neighbor_bonds[j];
            if (mol->bonds[bi].type == BOND_DOUBLE &&
                mol->atoms[atom->neighbors[j]].element == ELEM_O) {
                dbl_o++;
            }
        }
        if (dbl_o == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Halogen E-State Types
 * ============================================================================ */

/* sF: fluorine */
static cchem_status_t estate_sF(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_F) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* sCl: chlorine */
static cchem_status_t estate_sCl(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_Cl) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* sBr: bromine */
static cchem_status_t estate_sBr(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_Br) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* sI: iodine */
static cchem_status_t estate_sI(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_I) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_ESTATE(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_COUNTS; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_estate(void) {
    /* Carbon types */
    REGISTER_ESTATE("EState_ssCH3", "sp3 CH3 (methyl)", estate_ssCH3);
    REGISTER_ESTATE("EState_ssCH2", "sp3 CH2 (methylene)", estate_ssCH2);
    REGISTER_ESTATE("EState_sssCH", "sp3 CH (methine)", estate_sssCH);
    REGISTER_ESTATE("EState_ssssC", "sp3 quaternary C", estate_ssssC);
    REGISTER_ESTATE("EState_dCH2", "sp2 =CH2 (terminal alkene)", estate_dCH2);
    REGISTER_ESTATE("EState_dsCH", "sp2 =CH- (internal alkene)", estate_dsCH);
    REGISTER_ESTATE("EState_dssC", "sp2 =C< (substituted)", estate_dssC);
    REGISTER_ESTATE("EState_aaCH", "aromatic CH", estate_aaCH);
    REGISTER_ESTATE("EState_aasC", "aromatic C-R", estate_aasC);
    REGISTER_ESTATE("EState_tCH", "sp CH (alkyne)", estate_tCH);
    REGISTER_ESTATE("EState_tsC", "sp C (internal alkyne)", estate_tsC);

    /* Nitrogen types */
    REGISTER_ESTATE("EState_sNH3", "sp3 NH3+ (charged)", estate_sNH3);
    REGISTER_ESTATE("EState_sNH2", "sp3 -NH2 (primary amine)", estate_sNH2);
    REGISTER_ESTATE("EState_ssNH", "sp3 >NH (secondary amine)", estate_ssNH);
    REGISTER_ESTATE("EState_sssN", "sp3 >N- (tertiary amine)", estate_sssN);
    REGISTER_ESTATE("EState_dNH", "sp2 =NH (imine)", estate_dNH);
    REGISTER_ESTATE("EState_dsN", "sp2 =N- (imine)", estate_dsN);
    REGISTER_ESTATE("EState_aaNH", "aromatic NH (pyrrole)", estate_aaNH);
    REGISTER_ESTATE("EState_aaN", "aromatic N (pyridine)", estate_aaN);
    REGISTER_ESTATE("EState_tN", "sp N (nitrile)", estate_tN);

    /* Oxygen types */
    REGISTER_ESTATE("EState_sOH", "-OH (hydroxyl)", estate_sOH);
    REGISTER_ESTATE("EState_ssO", "-O- (ether)", estate_ssO);
    REGISTER_ESTATE("EState_dO", "=O (carbonyl)", estate_dO);
    REGISTER_ESTATE("EState_aaO", "aromatic O (furan)", estate_aaO);

    /* Sulfur types */
    REGISTER_ESTATE("EState_sSH", "-SH (thiol)", estate_sSH);
    REGISTER_ESTATE("EState_ssS", "-S- (thioether)", estate_ssS);
    REGISTER_ESTATE("EState_dS", "=S (thione)", estate_dS);
    REGISTER_ESTATE("EState_aaS", "aromatic S (thiophene)", estate_aaS);
    REGISTER_ESTATE("EState_dssS", "sulfoxide S", estate_dssS);
    REGISTER_ESTATE("EState_ddssS", "sulfone S", estate_ddssS);

    /* Halogen types */
    REGISTER_ESTATE("EState_sF", "fluorine", estate_sF);
    REGISTER_ESTATE("EState_sCl", "chlorine", estate_sCl);
    REGISTER_ESTATE("EState_sBr", "bromine", estate_sBr);
    REGISTER_ESTATE("EState_sI", "iodine", estate_sI);
}
