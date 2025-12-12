/**
 * @file mmff94_charges.c
 * @brief MMFF94 partial charge computation using Bond Charge Increments
 *
 * Computes partial atomic charges using the MMFF94 bond charge increment (BCI)
 * method. Each bond transfers a specific amount of charge between atoms based
 * on their types and the bond type.
 */

#include "cchem/depictor/mmff94_types.h"
#include "cchem/depictor/mmff94_params.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * BCI Parameter Lookup
 * ============================================================================ */

const mmff94_bci_param_t* mmff94_lookup_bci_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, int bond_type) {

    /* Linear search through BCI parameters */
    for (int p = 0; p < (int)MMFF94_NUM_BCI_PARAMS; p++) {
        const mmff94_bci_param_t* param = &MMFF94_BCI_PARAMS[p];

        /* Check forward match */
        if (param->type_i == type_i && param->type_j == type_j &&
            param->bond_type == bond_type) {
            return param;
        }

        /* Check reverse match (negated BCI) */
        if (param->type_i == type_j && param->type_j == type_i &&
            param->bond_type == bond_type) {
            /* We'll handle the sign flip in the calling code */
            return param;
        }
    }

    return NULL;  /* No parameter found */
}

/* ============================================================================
 * Formal Charge Handling
 * ============================================================================ */

/* Get formal charge contribution to partial charge */
static double get_formal_charge_contribution(mmff94_atom_type_t type) {
    for (int p = 0; p < (int)MMFF94_NUM_FC_PARAMS; p++) {
        if (MMFF94_FC_PARAMS[p].type == type) {
            return MMFF94_FC_PARAMS[p].fc_increment;
        }
    }
    return 0.0;
}

/* ============================================================================
 * Empirical BCI Estimation
 *
 * For missing BCI parameters, estimate from electronegativity differences.
 * MMFF94 uses: BCI = 0.1 * (chi_j - chi_i) where chi is electronegativity
 * ============================================================================ */

/* Pauling electronegativity by element (common organic elements) */
static double get_electronegativity(element_t elem) {
    switch (elem) {
        case ELEM_H:  return 2.20;
        case ELEM_C:  return 2.55;
        case ELEM_N:  return 3.04;
        case ELEM_O:  return 3.44;
        case ELEM_F:  return 3.98;
        case ELEM_P:  return 2.19;
        case ELEM_S:  return 2.58;
        case ELEM_Cl: return 3.16;
        case ELEM_Br: return 2.96;
        case ELEM_I:  return 2.66;
        case ELEM_Si: return 1.90;
        default:      return 2.50;  /* Default ~carbon */
    }
}

/* Get element from MMFF94 atom type */
static element_t get_element_from_type(mmff94_atom_type_t type) {
    /* Carbon types */
    if (type == MMFF94_TYPE_CR || type == MMFF94_TYPE_C_EQ_C ||
        type == MMFF94_TYPE_C_EQ_O || type == MMFF94_TYPE_CSP ||
        type == MMFF94_TYPE_C_AR) {
        return ELEM_C;
    }

    /* Nitrogen types */
    if (type == MMFF94_TYPE_NR || type == MMFF94_TYPE_N_EQ_C ||
        type == MMFF94_TYPE_N_AR || type == MMFF94_TYPE_NSP ||
        type == MMFF94_TYPE_NPL3 || type == MMFF94_TYPE_NR_PLUS ||
        type == MMFF94_TYPE_NM || type == MMFF94_TYPE_N2OX ||
        type == MMFF94_TYPE_NPOX || type == MMFF94_TYPE_N3OX ||
        type == MMFF94_TYPE_NPD_PLUS || type == MMFF94_TYPE_N_IM_PLUS) {
        return ELEM_N;
    }

    /* Oxygen types */
    if (type == MMFF94_TYPE_OR || type == MMFF94_TYPE_O_EQ_C ||
        type == MMFF94_TYPE_O_AR || type == MMFF94_TYPE_O_EQ_N ||
        type == MMFF94_TYPE_OM || type == MMFF94_TYPE_O_PLUS ||
        type == MMFF94_TYPE_O3_EQ_N) {
        return ELEM_O;
    }

    /* Sulfur types */
    if (type == MMFF94_TYPE_S || type == MMFF94_TYPE_S_EQ_C ||
        type == MMFF94_TYPE_S_AR || type == MMFF94_TYPE_SO ||
        type == MMFF94_TYPE_SO2 || type == MMFF94_TYPE_SO3 ||
        type == MMFF94_TYPE_S_MINUS) {
        return ELEM_S;
    }

    /* Phosphorus types */
    if (type == MMFF94_TYPE_P || type == MMFF94_TYPE_PO ||
        type == MMFF94_TYPE_PO2 || type == MMFF94_TYPE_PO4) {
        return ELEM_P;
    }

    /* Halogens */
    if (type == MMFF94_TYPE_F || type == MMFF94_TYPE_F_MINUS) return ELEM_F;
    if (type == MMFF94_TYPE_CL || type == MMFF94_TYPE_CL_MINUS) return ELEM_Cl;
    if (type == MMFF94_TYPE_BR || type == MMFF94_TYPE_BR_MINUS) return ELEM_Br;
    if (type == MMFF94_TYPE_I || type == MMFF94_TYPE_I_MINUS) return ELEM_I;

    /* Hydrogen types */
    if (type == MMFF94_TYPE_HC || type == MMFF94_TYPE_HO ||
        type == MMFF94_TYPE_HN || type == MMFF94_TYPE_HOCO ||
        type == MMFF94_TYPE_HOS || type == MMFF94_TYPE_HN_PLUS ||
        type == MMFF94_TYPE_HO_PLUS || type == MMFF94_TYPE_HNCO ||
        type == MMFF94_TYPE_HNCC || type == MMFF94_TYPE_HNCN ||
        type == MMFF94_TYPE_HP) {
        return ELEM_H;
    }

    /* Silicon */
    if (type == MMFF94_TYPE_SI) return ELEM_Si;

    return ELEM_C;  /* Default */
}

/* Estimate BCI from electronegativity */
static double estimate_bci(mmff94_atom_type_t type_i, mmff94_atom_type_t type_j) {
    element_t elem_i = get_element_from_type(type_i);
    element_t elem_j = get_element_from_type(type_j);

    double chi_i = get_electronegativity(elem_i);
    double chi_j = get_electronegativity(elem_j);

    /* BCI = k * (chi_j - chi_i), k ~ 0.1 */
    return 0.10 * (chi_j - chi_i);
}

/* ============================================================================
 * Main Charge Computation
 * ============================================================================ */

cchem_status_t mmff94_compute_charges(const molecule_t* mol, mmff94_context_t* ctx) {
    if (!mol || !ctx || !ctx->types_assigned) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    int n = mol->num_atoms;

    /* Initialize partial charges from formal charges */
    for (int i = 0; i < n; i++) {
        double fc = (double)mol->atoms[i].charge;
        double fc_contrib = get_formal_charge_contribution(ctx->atom_data[i].type);

        /* For ionic species, use the formal charge directly
         * For neutral species with partial ionic character, use fc_contrib */
        if (fc != 0.0 && fabs(fc_contrib) > 0.01) {
            ctx->atom_data[i].partial_charge = fc_contrib;
        } else {
            ctx->atom_data[i].partial_charge = fc;
        }
    }

    /* Apply bond charge increments */
    for (int b = 0; b < mol->num_bonds; b++) {
        int i = mol->bonds[b].atom1;
        int j = mol->bonds[b].atom2;

        mmff94_atom_type_t type_i = ctx->atom_data[i].type;
        mmff94_atom_type_t type_j = ctx->atom_data[j].type;

        /* Get bond type: 0=single, 1=double, 2=triple, 3=aromatic */
        int bond_type = 0;
        switch (mol->bonds[b].type) {
            case BOND_SINGLE:
            case BOND_UP:
            case BOND_DOWN:
                bond_type = 0;
                break;
            case BOND_DOUBLE:
                bond_type = 1;
                break;
            case BOND_TRIPLE:
                bond_type = 2;
                break;
            case BOND_AROMATIC:
                bond_type = 3;
                break;
            default:
                bond_type = 0;
                break;
        }

        /* Look up BCI parameter */
        const mmff94_bci_param_t* param = mmff94_lookup_bci_param(type_i, type_j, bond_type);

        double bci;
        if (param) {
            /* Check if we need to flip the sign (reverse direction) */
            if (param->type_i == type_i && param->type_j == type_j) {
                bci = param->bci;  /* Forward: charge flows j -> i */
            } else {
                bci = -param->bci;  /* Reverse: negate BCI */
            }
        } else {
            /* Estimate BCI from electronegativity */
            bci = estimate_bci(type_i, type_j);
        }

        /* Transfer charge: positive BCI means j -> i (i becomes more negative) */
        ctx->atom_data[i].partial_charge += bci;
        ctx->atom_data[j].partial_charge -= bci;
    }

    /* Handle implicit hydrogens */
    /* In MMFF94, implicit H contribute to charge through C-H BCI */
    for (int i = 0; i < n; i++) {
        int impl_h = mol->atoms[i].implicit_h_count;
        if (impl_h > 0) {
            mmff94_atom_type_t type_i = ctx->atom_data[i].type;

            /* Get BCI for X-H bond */
            const mmff94_bci_param_t* param = mmff94_lookup_bci_param(type_i, MMFF94_TYPE_HC, 0);

            double bci_h;
            if (param) {
                if (param->type_i == type_i) {
                    bci_h = param->bci;
                } else {
                    bci_h = -param->bci;
                }
            } else {
                /* Estimate: H is less electronegative than most atoms */
                element_t elem_i = get_element_from_type(type_i);
                double chi_i = get_electronegativity(elem_i);
                double chi_h = get_electronegativity(ELEM_H);
                bci_h = 0.10 * (chi_h - chi_i);
            }

            /* Apply for each implicit H */
            ctx->atom_data[i].partial_charge += impl_h * bci_h;
        }
    }

    /* Apply charge equalization corrections for special cases */

    /* Carboxylate: distribute charge equally between two oxygens */
    for (int i = 0; i < n; i++) {
        if (mol->atoms[i].element == ELEM_C && ctx->atom_data[i].type == MMFF94_TYPE_C_EQ_O) {
            /* Check if this is a carboxylate (C with 2 oxygens, one charged) */
            int o_count = 0;
            int charged_o = -1;
            int neutral_o = -1;

            for (int k = 0; k < mol->atoms[i].num_neighbors; k++) {
                int nb = mol->atoms[i].neighbors[k];
                if (mol->atoms[nb].element == ELEM_O) {
                    o_count++;
                    if (mol->atoms[nb].charge < 0) {
                        charged_o = nb;
                    } else if (ctx->atom_data[nb].type == MMFF94_TYPE_O_EQ_C) {
                        neutral_o = nb;
                    }
                }
            }

            if (o_count == 2 && charged_o >= 0 && neutral_o >= 0) {
                /* Average the charges on the two oxygens */
                double avg = (ctx->atom_data[charged_o].partial_charge +
                             ctx->atom_data[neutral_o].partial_charge) / 2.0;
                ctx->atom_data[charged_o].partial_charge = avg;
                ctx->atom_data[neutral_o].partial_charge = avg;
            }
        }
    }

    /* Nitro group: distribute charge equally between two oxygens */
    for (int i = 0; i < n; i++) {
        if (ctx->atom_data[i].type == MMFF94_TYPE_N2OX) {
            double o_charge_sum = 0.0;
            int o_count = 0;
            int o_indices[3] = {-1, -1, -1};

            for (int k = 0; k < mol->atoms[i].num_neighbors; k++) {
                int nb = mol->atoms[i].neighbors[k];
                if (mol->atoms[nb].element == ELEM_O && o_count < 3) {
                    o_charge_sum += ctx->atom_data[nb].partial_charge;
                    o_indices[o_count++] = nb;
                }
            }

            if (o_count >= 2) {
                double avg = o_charge_sum / o_count;
                for (int k = 0; k < o_count; k++) {
                    if (o_indices[k] >= 0) {
                        ctx->atom_data[o_indices[k]].partial_charge = avg;
                    }
                }
            }
        }
    }

    /* Guanidinium: distribute positive charge among nitrogens */
    for (int i = 0; i < n; i++) {
        if (mol->atoms[i].element == ELEM_C) {
            /* Check for guanidinium: C bonded to 3 N */
            int n_count = 0;
            int n_indices[4] = {-1, -1, -1, -1};

            for (int k = 0; k < mol->atoms[i].num_neighbors; k++) {
                int nb = mol->atoms[i].neighbors[k];
                if (mol->atoms[nb].element == ELEM_N && n_count < 4) {
                    n_indices[n_count++] = nb;
                }
            }

            if (n_count == 3) {
                /* Check if any N is charged (guanidinium) */
                bool has_charged_n = false;
                for (int k = 0; k < n_count; k++) {
                    if (mol->atoms[n_indices[k]].charge > 0) {
                        has_charged_n = true;
                        break;
                    }
                }

                if (has_charged_n) {
                    /* Distribute charge equally */
                    double total_charge = 0.0;
                    for (int k = 0; k < n_count; k++) {
                        total_charge += ctx->atom_data[n_indices[k]].partial_charge;
                    }
                    double avg = total_charge / n_count;
                    for (int k = 0; k < n_count; k++) {
                        ctx->atom_data[n_indices[k]].partial_charge = avg;
                    }
                }
            }
        }
    }

    ctx->charges_computed = true;
    return CCHEM_OK;
}

/* ============================================================================
 * Charge Validation and Statistics
 * ============================================================================ */

/* Calculate total charge (should equal formal charge) */
double mmff94_total_charge(const molecule_t* mol, const mmff94_context_t* ctx) {
    if (!mol || !ctx || !ctx->charges_computed) return 0.0;

    double total = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        total += ctx->atom_data[i].partial_charge;

        /* Add implicit hydrogen charges */
        int impl_h = mol->atoms[i].implicit_h_count;
        if (impl_h > 0) {
            /* Implicit H charge is typically small positive */
            /* The main atom already has the BCI contribution included */
            /* For total charge, we don't double count */
        }
    }

    return total;
}

/* Get dipole moment components */
void mmff94_dipole_moment(const molecule_t* mol, const mmff94_context_t* ctx,
                          const mol_coords_t* coords,
                          double* mu_x, double* mu_y, double* mu_z) {
    if (!mol || !ctx || !coords || !ctx->charges_computed) {
        if (mu_x) *mu_x = 0.0;
        if (mu_y) *mu_y = 0.0;
        if (mu_z) *mu_z = 0.0;
        return;
    }

    double dx = 0.0, dy = 0.0, dz = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        double q = ctx->atom_data[i].partial_charge;
        dx += q * coords->coords_3d[i].x;
        dy += q * coords->coords_3d[i].y;
        dz += q * coords->coords_3d[i].z;
    }

    /* Convert to Debye (1 Debye = 3.336e-30 C*m = 0.2082 e*Angstrom) */
    const double DEBYE_CONV = 4.8032;  /* e*A to Debye */

    if (mu_x) *mu_x = dx * DEBYE_CONV;
    if (mu_y) *mu_y = dy * DEBYE_CONV;
    if (mu_z) *mu_z = dz * DEBYE_CONV;
}
