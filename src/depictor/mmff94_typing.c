/**
 * @file mmff94_typing.c
 * @brief MMFF94 atom type assignment
 *
 * Assigns MMFF94 symbolic atom types based on element, hybridization,
 * aromaticity, ring membership, and chemical environment.
 */

#include "cchem/depictor/mmff94_types.h"
#include "cchem/depictor/mmff94_params.h"
#include "cchem/canonicalizer/ring_finder.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Hybridization Detection
 * ============================================================================ */

mmff94_hybridization_t mmff94_get_hybridization(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Count double and triple bonds */
    int num_double = 0;
    int num_triple = 0;
    int num_aromatic = 0;

    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        const bond_t* bond = &mol->bonds[bond_idx];

        switch (bond->type) {
            case BOND_DOUBLE:
                num_double++;
                break;
            case BOND_TRIPLE:
                num_triple++;
                break;
            case BOND_AROMATIC:
                num_aromatic++;
                break;
            default:
                break;
        }
    }

    /* sp hybridization: triple bond or two double bonds (allene) */
    if (num_triple > 0) {
        return MMFF94_HYBRID_SP;
    }

    /* Check for linear geometry with two double bonds */
    if (num_double >= 2 && atom->num_neighbors == 2) {
        return MMFF94_HYBRID_SP;
    }

    /* sp2 hybridization: double bond, aromatic, or trigonal planar */
    if (num_double > 0 || num_aromatic > 0 || atom->aromatic) {
        return MMFF94_HYBRID_SP2;
    }

    /* Default to sp3 */
    return MMFF94_HYBRID_SP3;
}

/* ============================================================================
 * Helper Functions for Atom Classification
 * ============================================================================ */

/* Check if atom has a double bond to a specific element */
static bool has_double_bond_to(const molecule_t* mol, int atom_idx, element_t elem) {
    const atom_t* atom = &mol->atoms[atom_idx];

    for (int i = 0; i < atom->num_neighbors; i++) {
        int nb_idx = atom->neighbors[i];
        int bond_idx = atom->neighbor_bonds[i];
        const bond_t* bond = &mol->bonds[bond_idx];

        if (bond->type == BOND_DOUBLE && mol->atoms[nb_idx].element == elem) {
            return true;
        }
    }
    return false;
}

/* Check if atom is bonded to element with specific bond type */
static bool is_bonded_to_element(const molecule_t* mol, int atom_idx, element_t elem, bond_type_t btype) {
    const atom_t* atom = &mol->atoms[atom_idx];

    for (int i = 0; i < atom->num_neighbors; i++) {
        int nb_idx = atom->neighbors[i];
        int bond_idx = atom->neighbor_bonds[i];
        const bond_t* bond = &mol->bonds[bond_idx];

        if (mol->atoms[nb_idx].element == elem) {
            if (btype == BOND_NONE || bond->type == btype) {
                return true;
            }
        }
    }
    return false;
}

/* Count neighbors of a specific element */
static int count_neighbors_element(const molecule_t* mol, int atom_idx, element_t elem) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int count = 0;

    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == elem) {
            count++;
        }
    }
    return count;
}

/* Check if carbon is a carbonyl carbon (C=O) */
static bool is_carbonyl_carbon(const molecule_t* mol, int atom_idx) {
    return has_double_bond_to(mol, atom_idx, ELEM_O);
}

/* Check if carbon is a thione carbon (C=S) */
static bool is_thione_carbon(const molecule_t* mol, int atom_idx) {
    return has_double_bond_to(mol, atom_idx, ELEM_S);
}

/* Check if nitrogen is in an amide (N-C=O) */
static bool is_amide_nitrogen(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    for (int i = 0; i < atom->num_neighbors; i++) {
        int nb_idx = atom->neighbors[i];
        if (mol->atoms[nb_idx].element == ELEM_C) {
            if (is_carbonyl_carbon(mol, nb_idx)) {
                return true;
            }
        }
    }
    return false;
}

/* Check if oxygen is in a nitro group */
static bool is_nitro_oxygen(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    for (int i = 0; i < atom->num_neighbors; i++) {
        int nb_idx = atom->neighbors[i];
        if (mol->atoms[nb_idx].element == ELEM_N) {
            /* Check if N has 3 oxygens (nitro) */
            if (count_neighbors_element(mol, nb_idx, ELEM_O) >= 2) {
                return true;
            }
        }
    }
    return false;
}

/* Check if atom is in 5-membered ring */
static bool in_5_ring(const molecule_t* mol, int atom_idx) {
    for (int r = 0; r < mol->num_rings; r++) {
        if (mol->rings[r].size == 5) {
            for (int i = 0; i < mol->rings[r].size; i++) {
                if (mol->rings[r].atoms[i] == atom_idx) {
                    return true;
                }
            }
        }
    }
    return false;
}

/* Check if atom is in 6-membered ring */
static bool in_6_ring(const molecule_t* mol, int atom_idx) {
    for (int r = 0; r < mol->num_rings; r++) {
        if (mol->rings[r].size == 6) {
            for (int i = 0; i < mol->rings[r].size; i++) {
                if (mol->rings[r].atoms[i] == atom_idx) {
                    return true;
                }
            }
        }
    }
    return false;
}

/* ============================================================================
 * Carbon Type Assignment
 * ============================================================================ */

static mmff94_atom_type_t assign_carbon_type(const molecule_t* mol, int atom_idx,
                                              mmff94_hybridization_t hybrid) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Aromatic carbon */
    if (atom->aromatic) {
        return MMFF94_TYPE_C_AR;
    }

    /* sp carbon */
    if (hybrid == MMFF94_HYBRID_SP) {
        return MMFF94_TYPE_CSP;
    }

    /* sp2 carbon */
    if (hybrid == MMFF94_HYBRID_SP2) {
        /* Carbonyl carbon */
        if (is_carbonyl_carbon(mol, atom_idx)) {
            return MMFF94_TYPE_C_EQ_O;
        }
        /* Thione carbon */
        if (is_thione_carbon(mol, atom_idx)) {
            return MMFF94_TYPE_C_EQ_O;  /* Use same type for C=S */
        }
        /* Imine carbon (C=N) */
        if (has_double_bond_to(mol, atom_idx, ELEM_N)) {
            return MMFF94_TYPE_C_EQ_C;  /* Generic sp2 */
        }
        /* Alkene carbon (C=C) */
        if (has_double_bond_to(mol, atom_idx, ELEM_C)) {
            return MMFF94_TYPE_C_EQ_C;
        }
        /* Default sp2 */
        return MMFF94_TYPE_C_EQ_C;
    }

    /* sp3 carbon (default) */
    return MMFF94_TYPE_CR;
}

/* ============================================================================
 * Nitrogen Type Assignment
 * ============================================================================ */

static mmff94_atom_type_t assign_nitrogen_type(const molecule_t* mol, int atom_idx,
                                                mmff94_hybridization_t hybrid) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Aromatic nitrogen */
    if (atom->aromatic) {
        /* Check if in 5-membered ring (pyrrole-like vs pyridine-like) */
        if (in_5_ring(mol, atom_idx)) {
            /* Pyrrole-type has 2 neighbors, pyridine-type in 5-ring has 2 also
             * but with different bonding pattern */
            if (atom->num_neighbors == 2 && !has_double_bond_to(mol, atom_idx, ELEM_C)) {
                return MMFF94_TYPE_NPL3;  /* Pyrrole-like */
            }
        }
        return MMFF94_TYPE_N_AR;  /* Pyridine-like aromatic N */
    }

    /* Positively charged nitrogen */
    if (atom->charge > 0) {
        if (hybrid == MMFF94_HYBRID_SP3) {
            return MMFF94_TYPE_NR_PLUS;  /* Quaternary ammonium */
        }
        if (atom->aromatic || in_6_ring(mol, atom_idx)) {
            return MMFF94_TYPE_NPD_PLUS;  /* Pyridinium */
        }
        return MMFF94_TYPE_NR_PLUS;
    }

    /* Negatively charged nitrogen */
    if (atom->charge < 0) {
        return MMFF94_TYPE_NM;
    }

    /* sp nitrogen (nitrile) */
    if (hybrid == MMFF94_HYBRID_SP) {
        return MMFF94_TYPE_NSP;
    }

    /* sp2 nitrogen */
    if (hybrid == MMFF94_HYBRID_SP2) {
        /* Nitro group */
        if (count_neighbors_element(mol, atom_idx, ELEM_O) >= 2) {
            bool has_double_O = has_double_bond_to(mol, atom_idx, ELEM_O);
            if (has_double_O) {
                return MMFF94_TYPE_N2OX;
            }
        }

        /* N-oxide */
        if (count_neighbors_element(mol, atom_idx, ELEM_O) >= 1 &&
            has_double_bond_to(mol, atom_idx, ELEM_O)) {
            return MMFF94_TYPE_NPOX;
        }

        /* Imine (N=C) */
        if (has_double_bond_to(mol, atom_idx, ELEM_C)) {
            return MMFF94_TYPE_N_EQ_C;
        }

        /* Amide nitrogen (N-C=O) - trigonal planar */
        if (is_amide_nitrogen(mol, atom_idx)) {
            return MMFF94_TYPE_NPL3;
        }

        /* General trigonal N */
        return MMFF94_TYPE_NPL3;
    }

    /* sp3 nitrogen (amine) */
    return MMFF94_TYPE_NR;
}

/* ============================================================================
 * Oxygen Type Assignment
 * ============================================================================ */

static mmff94_atom_type_t assign_oxygen_type(const molecule_t* mol, int atom_idx,
                                              mmff94_hybridization_t hybrid) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Aromatic oxygen (furan) */
    if (atom->aromatic) {
        return MMFF94_TYPE_O_AR;
    }

    /* Negatively charged oxygen */
    if (atom->charge < 0) {
        return MMFF94_TYPE_OM;
    }

    /* Positively charged oxygen (oxonium) */
    if (atom->charge > 0) {
        return MMFF94_TYPE_O_PLUS;
    }

    /* sp2 oxygen */
    if (hybrid == MMFF94_HYBRID_SP2) {
        /* Carbonyl oxygen (C=O) */
        if (has_double_bond_to(mol, atom_idx, ELEM_C)) {
            return MMFF94_TYPE_O_EQ_C;
        }
        /* N-oxide oxygen (N=O) */
        if (has_double_bond_to(mol, atom_idx, ELEM_N)) {
            return MMFF94_TYPE_O_EQ_N;
        }
        /* Sulfoxide/sulfone oxygen (S=O) */
        if (has_double_bond_to(mol, atom_idx, ELEM_S)) {
            return MMFF94_TYPE_O_EQ_S;
        }
        /* Phosphate oxygen (P=O) */
        if (has_double_bond_to(mol, atom_idx, ELEM_P)) {
            return MMFF94_TYPE_O_EQ_S;  /* Same treatment */
        }
        return MMFF94_TYPE_O_EQ_C;  /* Default sp2 O */
    }

    /* Nitro oxygen (may appear as single bond in some representations) */
    if (is_nitro_oxygen(mol, atom_idx)) {
        return MMFF94_TYPE_O_EQ_N;
    }

    /* sp3 oxygen (ether, alcohol) */
    return MMFF94_TYPE_OR;
}

/* ============================================================================
 * Sulfur Type Assignment
 * ============================================================================ */

static mmff94_atom_type_t assign_sulfur_type(const molecule_t* mol, int atom_idx,
                                              mmff94_hybridization_t hybrid __attribute__((unused))) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Aromatic sulfur (thiophene) */
    if (atom->aromatic) {
        return MMFF94_TYPE_S_AR;
    }

    /* Negatively charged sulfur (thiolate) */
    if (atom->charge < 0) {
        return MMFF94_TYPE_S_MINUS;
    }

    /* Count oxygen neighbors */
    int num_oxygen = count_neighbors_element(mol, atom_idx, ELEM_O);

    /* Sulfone (SO2) */
    if (num_oxygen >= 2) {
        return MMFF94_TYPE_SO2;
    }

    /* Sulfoxide (SO) */
    if (num_oxygen == 1) {
        return MMFF94_TYPE_SO;
    }

    /* Thione (C=S) */
    if (has_double_bond_to(mol, atom_idx, ELEM_C)) {
        return MMFF94_TYPE_S_EQ_C;
    }

    /* Default: thiol/thioether */
    return MMFF94_TYPE_S;
}

/* ============================================================================
 * Phosphorus Type Assignment
 * ============================================================================ */

static mmff94_atom_type_t assign_phosphorus_type(const molecule_t* mol, int atom_idx) {
    /* Count oxygen neighbors */
    int num_oxygen = count_neighbors_element(mol, atom_idx, ELEM_O);

    /* Phosphate/phosphonate (3+ oxygens) */
    if (num_oxygen >= 3) {
        return MMFF94_TYPE_PO4;
    }

    /* Phosphonate (2 oxygens) */
    if (num_oxygen >= 2) {
        return MMFF94_TYPE_PO2;
    }

    /* Phosphine oxide (1 oxygen) */
    if (num_oxygen >= 1) {
        return MMFF94_TYPE_PO;
    }

    /* Phosphine (no oxygen) */
    return MMFF94_TYPE_P;
}

/* ============================================================================
 * Hydrogen Type Assignment
 * ============================================================================ */

static mmff94_atom_type_t assign_hydrogen_type(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    if (atom->num_neighbors == 0) {
        return MMFF94_TYPE_HC;  /* Free hydrogen, treat as HC */
    }

    /* Get the atom hydrogen is bonded to */
    int heavy_idx = atom->neighbors[0];
    const atom_t* heavy = &mol->atoms[heavy_idx];
    element_t heavy_elem = heavy->element;

    switch (heavy_elem) {
        case ELEM_C:
            return MMFF94_TYPE_HC;

        case ELEM_N:
            /* Check if nitrogen is positively charged */
            if (heavy->charge > 0) {
                return MMFF94_TYPE_HN_PLUS;
            }
            /* Check if amide nitrogen */
            if (is_amide_nitrogen(mol, heavy_idx)) {
                return MMFF94_TYPE_HNCO;
            }
            /* Check if aniline/enamine nitrogen */
            if (heavy->aromatic || is_bonded_to_element(mol, heavy_idx, ELEM_C, BOND_DOUBLE)) {
                return MMFF94_TYPE_HNCC;
            }
            return MMFF94_TYPE_HN;

        case ELEM_O:
            /* Check if oxygen is positively charged */
            if (heavy->charge > 0) {
                return MMFF94_TYPE_HO_PLUS;
            }
            /* Check if carboxylic acid */
            for (int i = 0; i < heavy->num_neighbors; i++) {
                int c_idx = heavy->neighbors[i];
                if (mol->atoms[c_idx].element == ELEM_C && is_carbonyl_carbon(mol, c_idx)) {
                    return MMFF94_TYPE_HOCO;
                }
            }
            return MMFF94_TYPE_HO;

        case ELEM_S:
            return MMFF94_TYPE_HOS;

        case ELEM_P:
            return MMFF94_TYPE_HP;

        case ELEM_Si:
            return MMFF94_TYPE_HSI;

        default:
            return MMFF94_TYPE_HC;  /* Default to carbon H */
    }
}

/* ============================================================================
 * Halogen Type Assignment
 * ============================================================================ */

static mmff94_atom_type_t assign_halogen_type(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Check for ionic form (no bonds, negative charge) */
    if (atom->num_neighbors == 0 && atom->charge < 0) {
        switch (atom->element) {
            case ELEM_F:  return MMFF94_TYPE_F_MINUS;
            case ELEM_Cl: return MMFF94_TYPE_CL_MINUS;
            case ELEM_Br: return MMFF94_TYPE_BR_MINUS;
            case ELEM_I:  return MMFF94_TYPE_I_MINUS;
            default: break;
        }
    }

    /* Covalent halogen */
    switch (atom->element) {
        case ELEM_F:  return MMFF94_TYPE_F;
        case ELEM_Cl: return MMFF94_TYPE_CL;
        case ELEM_Br: return MMFF94_TYPE_BR;
        case ELEM_I:  return MMFF94_TYPE_I;
        default:      return MMFF94_TYPE_UNKNOWN;
    }
}

/* ============================================================================
 * Main Type Assignment Function
 * ============================================================================ */

mmff94_atom_type_t mmff94_get_atom_type(const molecule_t* mol, int atom_idx,
                                         const mmff94_context_t* ctx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    mmff94_hybridization_t hybrid = ctx->atom_data[atom_idx].hybridization;

    switch (atom->element) {
        case ELEM_C:
            return assign_carbon_type(mol, atom_idx, hybrid);

        case ELEM_N:
            return assign_nitrogen_type(mol, atom_idx, hybrid);

        case ELEM_O:
            return assign_oxygen_type(mol, atom_idx, hybrid);

        case ELEM_S:
            return assign_sulfur_type(mol, atom_idx, hybrid);

        case ELEM_P:
            return assign_phosphorus_type(mol, atom_idx);

        case ELEM_H:
            return assign_hydrogen_type(mol, atom_idx);

        case ELEM_F:
        case ELEM_Cl:
        case ELEM_Br:
        case ELEM_I:
            return assign_halogen_type(mol, atom_idx);

        case ELEM_Si:
            return MMFF94_TYPE_SI;

        case ELEM_Li:
            return (atom->charge > 0) ? MMFF94_TYPE_LI_PLUS : MMFF94_TYPE_UNKNOWN;
        case ELEM_Na:
            return (atom->charge > 0) ? MMFF94_TYPE_NA_PLUS : MMFF94_TYPE_UNKNOWN;
        case ELEM_K:
            return (atom->charge > 0) ? MMFF94_TYPE_K_PLUS : MMFF94_TYPE_UNKNOWN;
        case ELEM_Zn:
            return MMFF94_TYPE_ZN_2PLUS;
        case ELEM_Ca:
            return MMFF94_TYPE_CA_2PLUS;
        case ELEM_Mg:
            return MMFF94_TYPE_MG_2PLUS;
        case ELEM_Fe:
            return (atom->charge >= 3) ? MMFF94_TYPE_FE_3PLUS : MMFF94_TYPE_FE_2PLUS;
        case ELEM_Cu:
            return (atom->charge >= 2) ? MMFF94_TYPE_CU_2PLUS : MMFF94_TYPE_CU_PLUS;

        default:
            return MMFF94_TYPE_UNKNOWN;
    }
}

/* ============================================================================
 * Context Management
 * ============================================================================ */

mmff94_context_t* mmff94_context_create(const molecule_t* mol) {
    if (!mol || mol->num_atoms == 0) return NULL;

    mmff94_context_t* ctx = calloc(1, sizeof(mmff94_context_t));
    if (!ctx) return NULL;

    ctx->num_atoms = mol->num_atoms;
    ctx->atom_data = calloc(mol->num_atoms, sizeof(mmff94_atom_data_t));
    if (!ctx->atom_data) {
        free(ctx);
        return NULL;
    }

    ctx->num_bonds = mol->num_bonds;
    ctx->types_assigned = false;
    ctx->charges_computed = false;

    return ctx;
}

void mmff94_context_free(mmff94_context_t* ctx) {
    if (!ctx) return;

    free(ctx->atom_data);
    free(ctx->bond_type_i);
    free(ctx->bond_type_j);
    free(ctx->pairs_14);
    free(ctx);
}

/* ============================================================================
 * Public Type Assignment Function
 * ============================================================================ */

cchem_status_t mmff94_assign_types(const molecule_t* mol, mmff94_context_t* ctx) {
    if (!mol || !ctx) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* First pass: compute hybridization for all atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        ctx->atom_data[i].hybridization = mmff94_get_hybridization(mol, i);
        ctx->atom_data[i].is_aromatic = mol->atoms[i].aromatic;
        ctx->atom_data[i].is_cation = (mol->atoms[i].charge > 0);
        ctx->atom_data[i].is_anion = (mol->atoms[i].charge < 0);

        /* Get smallest ring size */
        ctx->atom_data[i].ring_size = 0;
        for (int r = 0; r < mol->num_rings; r++) {
            for (int j = 0; j < mol->rings[r].size; j++) {
                if (mol->rings[r].atoms[j] == i) {
                    if (ctx->atom_data[i].ring_size == 0 ||
                        mol->rings[r].size < ctx->atom_data[i].ring_size) {
                        ctx->atom_data[i].ring_size = mol->rings[r].size;
                    }
                    break;
                }
            }
        }
    }

    /* Second pass: assign atom types */
    for (int i = 0; i < mol->num_atoms; i++) {
        ctx->atom_data[i].type = mmff94_get_atom_type(mol, i, ctx);
    }

    ctx->types_assigned = true;
    return CCHEM_OK;
}

/* ============================================================================
 * Type Name Lookup
 * ============================================================================ */

const char* mmff94_type_name(mmff94_atom_type_t type) {
    static const char* names[] = {
        [MMFF94_TYPE_UNKNOWN] = "UNKNOWN",
        [MMFF94_TYPE_CR]      = "CR",
        [MMFF94_TYPE_C_EQ_C]  = "C=C",
        [MMFF94_TYPE_C_EQ_O]  = "C=O",
        [MMFF94_TYPE_CSP]     = "CSP",
        [MMFF94_TYPE_HC]      = "HC",
        [MMFF94_TYPE_OR]      = "OR",
        [MMFF94_TYPE_O_EQ_C]  = "O=C",
        [MMFF94_TYPE_NR]      = "NR",
        [MMFF94_TYPE_N_EQ_C]  = "N=C",
        [MMFF94_TYPE_N_EQ_O]  = "N=O",
        [MMFF94_TYPE_F]       = "F",
        [MMFF94_TYPE_CL]      = "CL",
        [MMFF94_TYPE_BR]      = "BR",
        [MMFF94_TYPE_I]       = "I",
        [MMFF94_TYPE_S]       = "S",
        [MMFF94_TYPE_S_EQ_C]  = "S=C",
        [MMFF94_TYPE_SO]      = "SO",
        [MMFF94_TYPE_SO2]     = "SO2",
        [MMFF94_TYPE_SI]      = "SI",
        [MMFF94_TYPE_HO]      = "HO",
        [MMFF94_TYPE_HN]      = "HN",
        [MMFF94_TYPE_HOCO]    = "HOCO",
        [MMFF94_TYPE_P]       = "P",
        [MMFF94_TYPE_PO]      = "PO",
        [MMFF94_TYPE_HNCN]    = "HNCN",
        [MMFF94_TYPE_HNCO]    = "HNCO",
        [MMFF94_TYPE_HNCC]    = "HNCC",
        [MMFF94_TYPE_O_EQ_N]  = "O=N",
        [MMFF94_TYPE_HOS]     = "HOS",
        [MMFF94_TYPE_NR_PLUS] = "NR+",
        [MMFF94_TYPE_OM]      = "OM",
        [MMFF94_TYPE_HN_PLUS] = "HN+",
        [MMFF94_TYPE_C_AR]    = "CB",
        [MMFF94_TYPE_N_AR]    = "NPYD",
        [MMFF94_TYPE_NPL3]    = "NPL3",
        [MMFF94_TYPE_NC_EQ_N] = "NC=N",
        [MMFF94_TYPE_O_AR]    = "OFUR",
        [MMFF94_TYPE_NSP]     = "NSP",
        [MMFF94_TYPE_S_AR]    = "STHI",
        [MMFF94_TYPE_N2OX]    = "N2OX",
        [MMFF94_TYPE_NAZT]    = "NAZT",
        [MMFF94_TYPE_O_PLUS]  = "O+",
        [MMFF94_TYPE_HO_PLUS] = "HO+",
        [MMFF94_TYPE_N_IM_PLUS] = "NIM+",
        [MMFF94_TYPE_NGD_PLUS] = "NGD+",
        [MMFF94_TYPE_NPD_PLUS] = "NPD+",
        [MMFF94_TYPE_NM]      = "NM",
        [MMFF94_TYPE_NPOX]    = "NPOX",
        [MMFF94_TYPE_N3OX]    = "N3OX",
        [MMFF94_TYPE_O3_EQ_N] = "O3N",
        [MMFF94_TYPE_HP]      = "HP",
        [MMFF94_TYPE_S_MINUS] = "S-",
        [MMFF94_TYPE_SO3]     = "SO3",
        [MMFF94_TYPE_NC_AR]   = "NC",
    };

    if (type >= 0 && type < (int)(sizeof(names)/sizeof(names[0])) && names[type]) {
        return names[type];
    }
    return "UNKNOWN";
}
