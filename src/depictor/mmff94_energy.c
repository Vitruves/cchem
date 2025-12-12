/**
 * @file mmff94_energy.c
 * @brief MMFF94 force field energy and gradient calculations
 *
 * Implements all seven MMFF94 energy terms:
 * 1. Bond stretching (quartic)
 * 2. Angle bending (cubic)
 * 3. Stretch-bend coupling
 * 4. Out-of-plane bending
 * 5. Torsion (3-term Fourier)
 * 6. Van der Waals (buffered 14-7)
 * 7. Electrostatics (buffered Coulomb)
 */

#include "cchem/depictor/mmff94_types.h"
#include "cchem/depictor/mmff94_params.h"
#include "cchem/depictor/types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Parameter Lookup Functions
 * ============================================================================ */

const mmff94_bond_param_t* mmff94_lookup_bond_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, int bond_type) {

    /* Canonicalize order: smaller type first */
    if (type_i > type_j) {
        mmff94_atom_type_t tmp = type_i;
        type_i = type_j;
        type_j = tmp;
    }

    for (int p = 0; p < (int)MMFF94_NUM_BOND_PARAMS; p++) {
        const mmff94_bond_param_t* param = &MMFF94_BOND_PARAMS[p];
        if ((param->type_i == type_i && param->type_j == type_j) ||
            (param->type_i == type_j && param->type_j == type_i)) {
            if (param->bond_type == bond_type) {
                return param;
            }
        }
    }
    return NULL;
}

const mmff94_angle_param_t* mmff94_lookup_angle_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, mmff94_atom_type_t type_k) {

    /* Canonicalize: smaller terminal type first */
    if (type_i > type_k) {
        mmff94_atom_type_t tmp = type_i;
        type_i = type_k;
        type_k = tmp;
    }

    for (int p = 0; p < (int)MMFF94_NUM_ANGLE_PARAMS; p++) {
        const mmff94_angle_param_t* param = &MMFF94_ANGLE_PARAMS[p];
        if (param->type_j == type_j) {
            if ((param->type_i == type_i && param->type_k == type_k) ||
                (param->type_i == type_k && param->type_k == type_i)) {
                return param;
            }
        }
    }
    return NULL;
}

const mmff94_strbnd_param_t* mmff94_lookup_strbnd_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, mmff94_atom_type_t type_k) {

    for (int p = 0; p < (int)MMFF94_NUM_STRBND_PARAMS; p++) {
        const mmff94_strbnd_param_t* param = &MMFF94_STRBND_PARAMS[p];
        if (param->type_j == type_j) {
            if ((param->type_i == type_i && param->type_k == type_k) ||
                (param->type_i == type_k && param->type_k == type_i)) {
                return param;
            }
        }
    }
    return NULL;
}

const mmff94_torsion_param_t* mmff94_lookup_torsion_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j,
    mmff94_atom_type_t type_k, mmff94_atom_type_t type_l) {

    /* First try exact match */
    for (int p = 0; p < (int)MMFF94_NUM_TORSION_PARAMS; p++) {
        const mmff94_torsion_param_t* param = &MMFF94_TORSION_PARAMS[p];
        if ((param->type_i == type_i && param->type_j == type_j &&
             param->type_k == type_k && param->type_l == type_l) ||
            (param->type_i == type_l && param->type_j == type_k &&
             param->type_k == type_j && param->type_l == type_i)) {
            return param;
        }
    }

    /* Try wildcard match (type 0 = any) */
    for (int p = 0; p < (int)MMFF94_NUM_TORSION_PARAMS; p++) {
        const mmff94_torsion_param_t* param = &MMFF94_TORSION_PARAMS[p];
        if ((param->type_i == 0 || param->type_i == type_i) &&
            (param->type_j == 0 || param->type_j == type_j) &&
            (param->type_k == 0 || param->type_k == type_k) &&
            (param->type_l == 0 || param->type_l == type_l)) {
            return param;
        }
        /* Reverse direction */
        if ((param->type_i == 0 || param->type_i == type_l) &&
            (param->type_j == 0 || param->type_j == type_k) &&
            (param->type_k == 0 || param->type_k == type_j) &&
            (param->type_l == 0 || param->type_l == type_i)) {
            return param;
        }
    }

    return NULL;
}

const mmff94_oop_param_t* mmff94_lookup_oop_param(
    mmff94_atom_type_t type_center,
    mmff94_atom_type_t type_j, mmff94_atom_type_t type_k, mmff94_atom_type_t type_l) {

    for (int p = 0; p < (int)MMFF94_NUM_OOP_PARAMS; p++) {
        const mmff94_oop_param_t* param = &MMFF94_OOP_PARAMS[p];
        if (param->type_i == type_center) {
            /* Check all permutations of j, k, l */
            if ((param->type_j == type_j && param->type_k == type_k && param->type_l == type_l) ||
                (param->type_j == type_j && param->type_k == type_l && param->type_l == type_k) ||
                (param->type_j == type_k && param->type_k == type_j && param->type_l == type_l) ||
                (param->type_j == type_k && param->type_k == type_l && param->type_l == type_j) ||
                (param->type_j == type_l && param->type_k == type_j && param->type_l == type_k) ||
                (param->type_j == type_l && param->type_k == type_k && param->type_l == type_j)) {
                return param;
            }
        }
    }
    return NULL;
}

const mmff94_vdw_param_t* mmff94_lookup_vdw_param(mmff94_atom_type_t type) {
    if (type >= 0 && type < (int)MMFF94_NUM_VDW_PARAMS) {
        return &MMFF94_VDW_PARAMS[type];
    }
    return &MMFF94_VDW_PARAMS[0];  /* Unknown type */
}

/* ============================================================================
 * VDW Combining Rules
 * ============================================================================ */

double mmff94_vdw_radius(const mmff94_vdw_param_t* param_i,
                         const mmff94_vdw_param_t* param_j) {
    /* MMFF94 combining rule for R* */
    double A_i = param_i->A > 0 ? param_i->A : 1.0;
    double A_j = param_j->A > 0 ? param_j->A : 1.0;
    double alpha_i = param_i->alpha > 0 ? param_i->alpha : 1.0;
    double alpha_j = param_j->alpha > 0 ? param_j->alpha : 1.0;

    double gamma = (param_i->G * param_j->G > 0) ?
                   (alpha_i - alpha_j) / (alpha_i + alpha_j) : 0.0;
    gamma = gamma * gamma;

    double R_i = A_i * pow(alpha_i, 0.25);
    double R_j = A_j * pow(alpha_j, 0.25);

    return 0.5 * (R_i + R_j) * (1.0 + 0.2 * (1.0 - exp(-12.0 * gamma)));
}

double mmff94_vdw_epsilon(const mmff94_vdw_param_t* param_i,
                          const mmff94_vdw_param_t* param_j,
                          double R_star) {
    /* MMFF94 combining rule for epsilon */
    double alpha_i = param_i->alpha > 0 ? param_i->alpha : 1.0;
    double alpha_j = param_j->alpha > 0 ? param_j->alpha : 1.0;
    double N_i = param_i->N_eff > 0 ? param_i->N_eff : 1.0;
    double N_j = param_j->N_eff > 0 ? param_j->N_eff : 1.0;

    if (R_star < 0.1) R_star = 0.1;

    double R_star3 = R_star * R_star * R_star;
    double epsilon = 181.16 * param_i->G * param_j->G * alpha_i * alpha_j /
                    ((sqrt(alpha_i / N_i) + sqrt(alpha_j / N_j)) * R_star3 * R_star3);

    return epsilon;
}

/* ============================================================================
 * Empirical Parameter Estimation
 * ============================================================================ */

void mmff94_estimate_bond_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, int bond_type,
    double* kb, double* r0) {

    /* Default bond parameters based on atomic radii */
    const mmff94_vdw_param_t* vdw_i = mmff94_lookup_vdw_param(type_i);
    const mmff94_vdw_param_t* vdw_j = mmff94_lookup_vdw_param(type_j);

    /* Estimate r0 from VDW radii (covalent ~ 0.7 * VDW) */
    double r_i = vdw_i->A > 0 ? pow(vdw_i->alpha, 0.25) * 0.7 : 0.75;
    double r_j = vdw_j->A > 0 ? pow(vdw_j->alpha, 0.25) * 0.7 : 0.75;

    *r0 = r_i + r_j;

    /* Adjust for bond order */
    switch (bond_type) {
        case 1: *r0 *= 0.87; break;  /* Double bond */
        case 2: *r0 *= 0.78; break;  /* Triple bond */
        case 3: *r0 *= 0.92; break;  /* Aromatic */
        default: break;
    }

    /* Estimate kb (larger for shorter bonds) */
    *kb = 5.0 / (*r0 * *r0);
}

void mmff94_estimate_angle_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, mmff94_atom_type_t type_k,
    double* ka, double* theta0) {

    (void)type_i; (void)type_k;  /* Suppress unused warnings */

    /* Default angle based on center atom hybridization */
    const mmff94_vdw_param_t* vdw_j = mmff94_lookup_vdw_param(type_j);

    /* Estimate theta0 based on expected geometry */
    /* This is a rough approximation */
    if (type_j == MMFF94_TYPE_CSP || type_j == MMFF94_TYPE_NSP) {
        *theta0 = 180.0;  /* Linear */
    } else if (type_j == MMFF94_TYPE_C_EQ_C || type_j == MMFF94_TYPE_C_EQ_O ||
               type_j == MMFF94_TYPE_C_AR || type_j == MMFF94_TYPE_N_AR ||
               type_j == MMFF94_TYPE_NPL3) {
        *theta0 = 120.0;  /* Trigonal */
    } else {
        *theta0 = 109.47;  /* Tetrahedral */
    }

    /* Estimate ka */
    *ka = 0.7;  /* Default force constant */

    (void)vdw_j;  /* Suppress unused warning */
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/* Check if atoms are 1-3 connected (share a common neighbor) */
static bool is_13_pair(const molecule_t* mol, int i, int j) {
    const atom_t* atom_i = &mol->atoms[i];

    for (int k = 0; k < atom_i->num_neighbors; k++) {
        int mid = atom_i->neighbors[k];
        const atom_t* atom_mid = &mol->atoms[mid];

        for (int l = 0; l < atom_mid->num_neighbors; l++) {
            if (atom_mid->neighbors[l] == j) {
                return true;
            }
        }
    }
    return false;
}

/* Check if atoms are 1-4 connected */
static bool is_14_pair(const molecule_t* mol, int i, int j) {
    const atom_t* atom_i = &mol->atoms[i];

    for (int k1 = 0; k1 < atom_i->num_neighbors; k1++) {
        int a1 = atom_i->neighbors[k1];
        const atom_t* atom_a1 = &mol->atoms[a1];

        for (int k2 = 0; k2 < atom_a1->num_neighbors; k2++) {
            int a2 = atom_a1->neighbors[k2];
            if (a2 == i) continue;
            const atom_t* atom_a2 = &mol->atoms[a2];

            for (int k3 = 0; k3 < atom_a2->num_neighbors; k3++) {
                if (atom_a2->neighbors[k3] == j) {
                    return true;
                }
            }
        }
    }
    return false;
}

/* Calculate dihedral angle between four atoms */
static double calc_dihedral(const mol_coords_t* coords, int i, int j, int k, int l) {
    /* Vectors */
    double b1x = coords->coords_3d[j].x - coords->coords_3d[i].x;
    double b1y = coords->coords_3d[j].y - coords->coords_3d[i].y;
    double b1z = coords->coords_3d[j].z - coords->coords_3d[i].z;

    double b2x = coords->coords_3d[k].x - coords->coords_3d[j].x;
    double b2y = coords->coords_3d[k].y - coords->coords_3d[j].y;
    double b2z = coords->coords_3d[k].z - coords->coords_3d[j].z;

    double b3x = coords->coords_3d[l].x - coords->coords_3d[k].x;
    double b3y = coords->coords_3d[l].y - coords->coords_3d[k].y;
    double b3z = coords->coords_3d[l].z - coords->coords_3d[k].z;

    /* Normal vectors n1 = b1 x b2, n2 = b2 x b3 */
    double n1x = b1y * b2z - b1z * b2y;
    double n1y = b1z * b2x - b1x * b2z;
    double n1z = b1x * b2y - b1y * b2x;

    double n2x = b2y * b3z - b2z * b3y;
    double n2y = b2z * b3x - b2x * b3z;
    double n2z = b2x * b3y - b2y * b3x;

    /* m1 = n1 x b2_normalized */
    double b2_len = sqrt(b2x*b2x + b2y*b2y + b2z*b2z);
    if (b2_len < 1e-10) return 0.0;

    double b2nx = b2x / b2_len;
    double b2ny = b2y / b2_len;
    double b2nz = b2z / b2_len;

    double m1x = n1y * b2nz - n1z * b2ny;
    double m1y = n1z * b2nx - n1x * b2nz;
    double m1z = n1x * b2ny - n1y * b2nx;

    /* Dihedral angle */
    double x = n1x * n2x + n1y * n2y + n1z * n2z;
    double y = m1x * n2x + m1y * n2y + m1z * n2z;

    return atan2(y, x);
}

/* Calculate Wilson out-of-plane angle */
static double calc_oop_angle(const mol_coords_t* coords, int center, int j, int k, int l) {
    /* Vector from center to each ligand */
    double rjx = coords->coords_3d[j].x - coords->coords_3d[center].x;
    double rjy = coords->coords_3d[j].y - coords->coords_3d[center].y;
    double rjz = coords->coords_3d[j].z - coords->coords_3d[center].z;

    double rkx = coords->coords_3d[k].x - coords->coords_3d[center].x;
    double rky = coords->coords_3d[k].y - coords->coords_3d[center].y;
    double rkz = coords->coords_3d[k].z - coords->coords_3d[center].z;

    double rlx = coords->coords_3d[l].x - coords->coords_3d[center].x;
    double rly = coords->coords_3d[l].y - coords->coords_3d[center].y;
    double rlz = coords->coords_3d[l].z - coords->coords_3d[center].z;

    /* Normal to plane defined by k and l */
    double nx = rky * rlz - rkz * rly;
    double ny = rkz * rlx - rkx * rlz;
    double nz = rkx * rly - rky * rlx;

    double n_len = sqrt(nx*nx + ny*ny + nz*nz);
    double rj_len = sqrt(rjx*rjx + rjy*rjy + rjz*rjz);

    if (n_len < 1e-10 || rj_len < 1e-10) return 0.0;

    /* sin(chi) = (rj . n) / (|rj| * |n|) */
    double sin_chi = (rjx*nx + rjy*ny + rjz*nz) / (rj_len * n_len);

    /* Clamp to [-1, 1] */
    if (sin_chi > 1.0) sin_chi = 1.0;
    if (sin_chi < -1.0) sin_chi = -1.0;

    return asin(sin_chi);
}

/* ============================================================================
 * Bond Stretching Energy (Quartic)
 * E = 143.9325 * kb * dr^2 * (1 + cs*dr + 7/12*cs^2*dr^2)
 * ============================================================================ */

double mmff94_calc_bond_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                const mol_coords_t* coords,
                                double* gx, double* gy, double* gz) {
    double energy = 0.0;
    const double cs = MMFF94_CUBIC_STRETCH;
    const double cs2 = cs * cs;

    for (int b = 0; b < mol->num_bonds; b++) {
        int i = mol->bonds[b].atom1;
        int j = mol->bonds[b].atom2;

        /* Get bond type */
        int bond_type = 0;
        switch (mol->bonds[b].type) {
            case BOND_DOUBLE: bond_type = 1; break;
            case BOND_TRIPLE: bond_type = 2; break;
            case BOND_AROMATIC: bond_type = 3; break;
            default: bond_type = 0; break;
        }

        /* Lookup parameters */
        const mmff94_bond_param_t* param = mmff94_lookup_bond_param(
            ctx->atom_data[i].type, ctx->atom_data[j].type, bond_type);

        double kb, r0;
        if (param) {
            kb = param->kb;
            r0 = param->r0;
        } else {
            mmff94_estimate_bond_param(ctx->atom_data[i].type, ctx->atom_data[j].type,
                                       bond_type, &kb, &r0);
        }

        /* Calculate distance */
        double dx = coords->coords_3d[j].x - coords->coords_3d[i].x;
        double dy = coords->coords_3d[j].y - coords->coords_3d[i].y;
        double dz = coords->coords_3d[j].z - coords->coords_3d[i].z;
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 0.01) r = 0.01;

        double dr = r - r0;
        double dr2 = dr * dr;

        /* Quartic energy */
        double factor = 1.0 + cs*dr + (7.0/12.0)*cs2*dr2;
        energy += MMFF94_BOND_CONST * kb * dr2 * factor;

        /* Gradient */
        if (gx) {
            double dfactor = cs + (7.0/6.0)*cs2*dr;
            double dE_dr = MMFF94_BOND_CONST * kb * 2.0 * dr * (factor + dr * dfactor);
            double g = dE_dr / r;

            gx[i] -= g * dx;
            gy[i] -= g * dy;
            gz[i] -= g * dz;
            gx[j] += g * dx;
            gy[j] += g * dy;
            gz[j] += g * dz;
        }
    }

    return energy;
}

/* ============================================================================
 * Angle Bending Energy (Cubic)
 * E = 0.043844 * ka * dtheta^2 * (1 + cb*dtheta)
 * ============================================================================ */

double mmff94_calc_angle_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                 const mol_coords_t* coords,
                                 double* gx, double* gy, double* gz) {
    double energy = 0.0;
    const double cb = MMFF94_CUBIC_BEND;
    const double DEG_TO_RAD = M_PI / 180.0;

    for (int j = 0; j < mol->num_atoms; j++) {
        if (mol->atoms[j].num_neighbors < 2) continue;

        for (int ni = 0; ni < mol->atoms[j].num_neighbors; ni++) {
            for (int nk = ni + 1; nk < mol->atoms[j].num_neighbors; nk++) {
                int i = mol->atoms[j].neighbors[ni];
                int k = mol->atoms[j].neighbors[nk];

                /* Lookup parameters */
                const mmff94_angle_param_t* param = mmff94_lookup_angle_param(
                    ctx->atom_data[i].type, ctx->atom_data[j].type, ctx->atom_data[k].type);

                double ka, theta0;
                if (param) {
                    ka = param->ka;
                    theta0 = param->theta0 * DEG_TO_RAD;
                } else {
                    mmff94_estimate_angle_param(ctx->atom_data[i].type, ctx->atom_data[j].type,
                                                ctx->atom_data[k].type, &ka, &theta0);
                    theta0 *= DEG_TO_RAD;
                }

                /* Calculate angle */
                double v1x = coords->coords_3d[i].x - coords->coords_3d[j].x;
                double v1y = coords->coords_3d[i].y - coords->coords_3d[j].y;
                double v1z = coords->coords_3d[i].z - coords->coords_3d[j].z;
                double r1 = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);

                double v2x = coords->coords_3d[k].x - coords->coords_3d[j].x;
                double v2y = coords->coords_3d[k].y - coords->coords_3d[j].y;
                double v2z = coords->coords_3d[k].z - coords->coords_3d[j].z;
                double r2 = sqrt(v2x*v2x + v2y*v2y + v2z*v2z);

                if (r1 < 0.01 || r2 < 0.01) continue;

                double cos_theta = (v1x*v2x + v1y*v2y + v1z*v2z) / (r1 * r2);
                if (cos_theta > 1.0) cos_theta = 1.0;
                if (cos_theta < -1.0) cos_theta = -1.0;

                double theta = acos(cos_theta);
                double dtheta = theta - theta0;

                /* Convert dtheta to degrees for cb (cb is in 1/deg) */
                double dtheta_deg = dtheta * 180.0 / M_PI;

                /* Cubic energy */
                energy += MMFF94_ANGLE_CONST * ka * dtheta * dtheta * (1.0 + cb * dtheta_deg);

                /* Gradient */
                if (gx) {
                    double sin_theta = sin(theta);
                    if (fabs(sin_theta) < 1e-6) continue;

                    double dE_dtheta = MMFF94_ANGLE_CONST * ka * 2.0 * dtheta *
                                      (1.0 + cb * dtheta_deg) +
                                      MMFF94_ANGLE_CONST * ka * dtheta * dtheta * cb * (180.0 / M_PI);

                    double factor = dE_dtheta / (-sin_theta);

                    double di_x = (v2x / (r1 * r2) - cos_theta * v1x / (r1 * r1));
                    double di_y = (v2y / (r1 * r2) - cos_theta * v1y / (r1 * r1));
                    double di_z = (v2z / (r1 * r2) - cos_theta * v1z / (r1 * r1));

                    double dk_x = (v1x / (r1 * r2) - cos_theta * v2x / (r2 * r2));
                    double dk_y = (v1y / (r1 * r2) - cos_theta * v2y / (r2 * r2));
                    double dk_z = (v1z / (r1 * r2) - cos_theta * v2z / (r2 * r2));

                    gx[i] += factor * di_x;
                    gy[i] += factor * di_y;
                    gz[i] += factor * di_z;

                    gx[k] += factor * dk_x;
                    gy[k] += factor * dk_y;
                    gz[k] += factor * dk_z;

                    gx[j] -= factor * (di_x + dk_x);
                    gy[j] -= factor * (di_y + dk_y);
                    gz[j] -= factor * (di_z + dk_z);
                }
            }
        }
    }

    return energy;
}

/* ============================================================================
 * Stretch-Bend Energy
 * E = 2.51124 * (kba_ijk * dr_ij + kba_kji * dr_jk) * dtheta
 * ============================================================================ */

double mmff94_calc_strbnd_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                  const mol_coords_t* coords,
                                  double* gx, double* gy, double* gz) {
    double energy = 0.0;
    const double DEG_TO_RAD = M_PI / 180.0;

    for (int j = 0; j < mol->num_atoms; j++) {
        if (mol->atoms[j].num_neighbors < 2) continue;

        for (int ni = 0; ni < mol->atoms[j].num_neighbors; ni++) {
            for (int nk = ni + 1; nk < mol->atoms[j].num_neighbors; nk++) {
                int i = mol->atoms[j].neighbors[ni];
                int k = mol->atoms[j].neighbors[nk];

                /* Lookup stretch-bend parameters */
                const mmff94_strbnd_param_t* sb_param = mmff94_lookup_strbnd_param(
                    ctx->atom_data[i].type, ctx->atom_data[j].type, ctx->atom_data[k].type);

                if (!sb_param) continue;  /* Skip if no stretch-bend parameters */

                double kba_ij = sb_param->kba_ijk;
                double kba_jk = sb_param->kba_kji;

                /* Get bond equilibrium lengths */
                int bond_ij_type = 0, bond_jk_type = 0;
                const bond_t* b_ij = molecule_get_bond_between_const(mol, i, j);
                const bond_t* b_jk = molecule_get_bond_between_const(mol, j, k);

                if (b_ij) {
                    switch (b_ij->type) {
                        case BOND_DOUBLE: bond_ij_type = 1; break;
                        case BOND_TRIPLE: bond_ij_type = 2; break;
                        case BOND_AROMATIC: bond_ij_type = 3; break;
                        default: bond_ij_type = 0; break;
                    }
                }
                if (b_jk) {
                    switch (b_jk->type) {
                        case BOND_DOUBLE: bond_jk_type = 1; break;
                        case BOND_TRIPLE: bond_jk_type = 2; break;
                        case BOND_AROMATIC: bond_jk_type = 3; break;
                        default: bond_jk_type = 0; break;
                    }
                }

                const mmff94_bond_param_t* bp_ij = mmff94_lookup_bond_param(
                    ctx->atom_data[i].type, ctx->atom_data[j].type, bond_ij_type);
                const mmff94_bond_param_t* bp_jk = mmff94_lookup_bond_param(
                    ctx->atom_data[j].type, ctx->atom_data[k].type, bond_jk_type);

                double r0_ij = bp_ij ? bp_ij->r0 : 1.5;
                double r0_jk = bp_jk ? bp_jk->r0 : 1.5;

                /* Get angle equilibrium */
                const mmff94_angle_param_t* a_param = mmff94_lookup_angle_param(
                    ctx->atom_data[i].type, ctx->atom_data[j].type, ctx->atom_data[k].type);
                double theta0 = a_param ? a_param->theta0 * DEG_TO_RAD : 109.47 * DEG_TO_RAD;

                /* Calculate current geometry */
                double v_ij_x = coords->coords_3d[i].x - coords->coords_3d[j].x;
                double v_ij_y = coords->coords_3d[i].y - coords->coords_3d[j].y;
                double v_ij_z = coords->coords_3d[i].z - coords->coords_3d[j].z;
                double r_ij = sqrt(v_ij_x*v_ij_x + v_ij_y*v_ij_y + v_ij_z*v_ij_z);

                double v_jk_x = coords->coords_3d[k].x - coords->coords_3d[j].x;
                double v_jk_y = coords->coords_3d[k].y - coords->coords_3d[j].y;
                double v_jk_z = coords->coords_3d[k].z - coords->coords_3d[j].z;
                double r_jk = sqrt(v_jk_x*v_jk_x + v_jk_y*v_jk_y + v_jk_z*v_jk_z);

                if (r_ij < 0.01 || r_jk < 0.01) continue;

                double cos_theta = (v_ij_x*v_jk_x + v_ij_y*v_jk_y + v_ij_z*v_jk_z) / (r_ij * r_jk);
                if (cos_theta > 1.0) cos_theta = 1.0;
                if (cos_theta < -1.0) cos_theta = -1.0;
                double theta = acos(cos_theta);

                double dr_ij = r_ij - r0_ij;
                double dr_jk = r_jk - r0_jk;
                double dtheta = theta - theta0;

                /* Stretch-bend energy */
                energy += MMFF94_STRBND_CONST * (kba_ij * dr_ij + kba_jk * dr_jk) * dtheta;

                /* Gradient computation is complex - simplified here */
                if (gx) {
                    /* This is an approximation; full gradient requires chain rule */
                    /* dE/dr_ij = STRBND_CONST * kba_ij * dtheta */
                    /* dE/dr_jk = STRBND_CONST * kba_jk * dtheta */
                    /* dE/dtheta = STRBND_CONST * (kba_ij * dr_ij + kba_jk * dr_jk) */

                    double dE_dr_ij = MMFF94_STRBND_CONST * kba_ij * dtheta;
                    double dE_dr_jk = MMFF94_STRBND_CONST * kba_jk * dtheta;

                    /* Bond gradient contribution */
                    gx[i] += dE_dr_ij * v_ij_x / r_ij;
                    gy[i] += dE_dr_ij * v_ij_y / r_ij;
                    gz[i] += dE_dr_ij * v_ij_z / r_ij;
                    gx[j] -= dE_dr_ij * v_ij_x / r_ij;
                    gy[j] -= dE_dr_ij * v_ij_y / r_ij;
                    gz[j] -= dE_dr_ij * v_ij_z / r_ij;

                    gx[k] += dE_dr_jk * v_jk_x / r_jk;
                    gy[k] += dE_dr_jk * v_jk_y / r_jk;
                    gz[k] += dE_dr_jk * v_jk_z / r_jk;
                    gx[j] -= dE_dr_jk * v_jk_x / r_jk;
                    gy[j] -= dE_dr_jk * v_jk_y / r_jk;
                    gz[j] -= dE_dr_jk * v_jk_z / r_jk;
                }
            }
        }
    }

    return energy;
}

/* ============================================================================
 * Torsion Energy (3-term Fourier)
 * E = 0.5 * (V1*(1+cos(phi)) + V2*(1-cos(2*phi)) + V3*(1+cos(3*phi)))
 * ============================================================================ */

double mmff94_calc_torsion_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                   const mol_coords_t* coords,
                                   double* gx, double* gy __attribute__((unused)), double* gz __attribute__((unused))) {
    double energy = 0.0;

    /* Iterate over all bonds (central bond of torsion) */
    for (int b = 0; b < mol->num_bonds; b++) {
        int j = mol->bonds[b].atom1;
        int k = mol->bonds[b].atom2;

        /* For each neighbor i of j (i != k) */
        for (int ni = 0; ni < mol->atoms[j].num_neighbors; ni++) {
            int i = mol->atoms[j].neighbors[ni];
            if (i == k) continue;

            /* For each neighbor l of k (l != j) */
            for (int nl = 0; nl < mol->atoms[k].num_neighbors; nl++) {
                int l = mol->atoms[k].neighbors[nl];
                if (l == j || l == i) continue;

                /* Lookup torsion parameters */
                const mmff94_torsion_param_t* param = mmff94_lookup_torsion_param(
                    ctx->atom_data[i].type, ctx->atom_data[j].type,
                    ctx->atom_data[k].type, ctx->atom_data[l].type);

                if (!param) continue;  /* Skip if no parameters */

                /* Calculate dihedral angle */
                double phi = calc_dihedral(coords, i, j, k, l);

                /* 3-term Fourier energy */
                double E = MMFF94_TORSION_CONST * (
                    param->V1 * (1.0 + cos(phi)) +
                    param->V2 * (1.0 - cos(2.0 * phi)) +
                    param->V3 * (1.0 + cos(3.0 * phi)));
                energy += E;

                /* Gradient (simplified - full torsion gradient is complex) */
                if (gx) {
                    double dE_dphi = MMFF94_TORSION_CONST * (
                        -param->V1 * sin(phi) +
                        2.0 * param->V2 * sin(2.0 * phi) -
                        3.0 * param->V3 * sin(3.0 * phi));

                    /* Approximate gradient using central differences on phi */
                    /* This is a placeholder - full implementation requires
                       analytical dihedral derivatives */
                    (void)dE_dphi;
                }
            }
        }
    }

    return energy;
}

/* ============================================================================
 * Out-of-Plane Energy
 * E = 0.043844 * koop * chi^2
 * ============================================================================ */

double mmff94_calc_oop_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                               const mol_coords_t* coords,
                               double* gx, double* gy, double* gz) {
    double energy = 0.0;

    /* Find trigonal planar centers */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].num_neighbors != 3) continue;

        /* Check if sp2 (planar) */
        if (ctx->atom_data[i].hybridization != MMFF94_HYBRID_SP2) continue;

        int j = mol->atoms[i].neighbors[0];
        int k = mol->atoms[i].neighbors[1];
        int l = mol->atoms[i].neighbors[2];

        /* Lookup OOP parameters */
        const mmff94_oop_param_t* param = mmff94_lookup_oop_param(
            ctx->atom_data[i].type,
            ctx->atom_data[j].type, ctx->atom_data[k].type, ctx->atom_data[l].type);

        /* Default value depends on aromaticity - aromatic gets stronger planarity */
        double default_koop = mol->atoms[i].aromatic ? 0.350 : 0.150;
        double koop = param ? param->koop : default_koop;

        /* Calculate OOP angle for each ligand */
        /* Wilson angle: angle of ligand from plane of other two ligands */
        double chi_j = calc_oop_angle(coords, i, j, k, l);
        double chi_k = calc_oop_angle(coords, i, k, j, l);
        double chi_l = calc_oop_angle(coords, i, l, j, k);

        /* Average the three Wilson angles */
        double chi_avg = (chi_j*chi_j + chi_k*chi_k + chi_l*chi_l) / 3.0;
        energy += MMFF94_OOP_CONST * koop * chi_avg;

        /* Gradient - use numerical approximation for stability */
        if (gx) {
            double dE = MMFF94_OOP_CONST * koop;
            double h = 0.001;  /* Small displacement */

            /* For each atom involved, compute gradient numerically */
            int atoms[4] = {i, j, k, l};
            double* grad_arrays[3] = {gx, gy, gz};

            for (int a = 0; a < 4; a++) {
                int atom_idx = atoms[a];

                for (int dim = 0; dim < 3; dim++) {
                    /* Save original position */
                    double orig;
                    if (dim == 0) orig = coords->coords_3d[atom_idx].x;
                    else if (dim == 1) orig = coords->coords_3d[atom_idx].y;
                    else orig = coords->coords_3d[atom_idx].z;

                    /* Displace forward */
                    mol_coords_t* mcoords = (mol_coords_t*)coords;  /* Cast away const for gradient calc */
                    if (dim == 0) mcoords->coords_3d[atom_idx].x = orig + h;
                    else if (dim == 1) mcoords->coords_3d[atom_idx].y = orig + h;
                    else mcoords->coords_3d[atom_idx].z = orig + h;

                    double chi_j_p = calc_oop_angle(coords, i, j, k, l);
                    double chi_k_p = calc_oop_angle(coords, i, k, j, l);
                    double chi_l_p = calc_oop_angle(coords, i, l, j, k);
                    double chi_avg_p = (chi_j_p*chi_j_p + chi_k_p*chi_k_p + chi_l_p*chi_l_p) / 3.0;

                    /* Displace backward */
                    if (dim == 0) mcoords->coords_3d[atom_idx].x = orig - h;
                    else if (dim == 1) mcoords->coords_3d[atom_idx].y = orig - h;
                    else mcoords->coords_3d[atom_idx].z = orig - h;

                    double chi_j_m = calc_oop_angle(coords, i, j, k, l);
                    double chi_k_m = calc_oop_angle(coords, i, k, j, l);
                    double chi_l_m = calc_oop_angle(coords, i, l, j, k);
                    double chi_avg_m = (chi_j_m*chi_j_m + chi_k_m*chi_k_m + chi_l_m*chi_l_m) / 3.0;

                    /* Restore original position */
                    if (dim == 0) mcoords->coords_3d[atom_idx].x = orig;
                    else if (dim == 1) mcoords->coords_3d[atom_idx].y = orig;
                    else mcoords->coords_3d[atom_idx].z = orig;

                    /* Central difference gradient */
                    double gradient = dE * (chi_avg_p - chi_avg_m) / (2.0 * h);
                    grad_arrays[dim][atom_idx] += gradient;
                }
            }
        }
    }

    return energy;
}

/* ============================================================================
 * Van der Waals Energy (Buffered 14-7)
 * E = epsilon * (1.07*Rstar/(R+0.07*Rstar))^7 * ((1.12*Rstar^7/(R^7+0.12*Rstar^7)) - 2)
 * ============================================================================ */

double mmff94_calc_vdw_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                               const mol_coords_t* coords,
                               double* gx, double* gy, double* gz) {
    double energy = 0.0;
    int n = mol->num_atoms;

    for (int i = 0; i < n; i++) {
        const mmff94_vdw_param_t* vdw_i = mmff94_lookup_vdw_param(ctx->atom_data[i].type);

        for (int j = i + 1; j < n; j++) {
            /* Skip bonded pairs (1-2) */
            if (molecule_get_bond_between_const(mol, i, j)) continue;

            /* Skip 1-3 pairs */
            if (is_13_pair(mol, i, j)) continue;

            const mmff94_vdw_param_t* vdw_j = mmff94_lookup_vdw_param(ctx->atom_data[j].type);

            /* Calculate distance */
            double dx = coords->coords_3d[j].x - coords->coords_3d[i].x;
            double dy = coords->coords_3d[j].y - coords->coords_3d[i].y;
            double dz = coords->coords_3d[j].z - coords->coords_3d[i].z;
            double R = sqrt(dx*dx + dy*dy + dz*dz);
            if (R < 0.5) R = 0.5;  /* Prevent singularity */

            /* Combining rules */
            double R_star = mmff94_vdw_radius(vdw_i, vdw_j);
            double epsilon = mmff94_vdw_epsilon(vdw_i, vdw_j, R_star);

            if (R_star < 0.1) R_star = 0.1;
            if (epsilon < 0.0) epsilon = 0.0;

            /* Buffered 14-7 potential */
            double rho = R / R_star;
            double rho7 = rho * rho * rho * rho * rho * rho * rho;

            double term1 = 1.07 / (rho + MMFF94_VDW_DELTA);
            double term1_7 = term1 * term1 * term1 * term1 * term1 * term1 * term1;

            double term2 = 1.12 / (rho7 + MMFF94_VDW_GAMMA) - 2.0;

            double E = epsilon * term1_7 * term2;

            /* Apply 1-4 scaling */
            if (is_14_pair(mol, i, j)) {
                E *= MMFF94_VDW_14_SCALE;
            }

            energy += E;

            /* Gradient */
            if (gx) {
                /* d(term1_7)/drho = 7 * term1^6 * (-1.07/(rho + delta)^2) / R_star */
                double dterm1_drho = -1.07 / ((rho + MMFF94_VDW_DELTA) * (rho + MMFF94_VDW_DELTA));
                double dterm1_7_drho = 7.0 * term1_7 / term1 * dterm1_drho;

                /* d(term2)/drho = -1.12 * 7 * rho^6 / (rho^7 + gamma)^2 */
                double rho6 = rho7 / rho;
                double denom = rho7 + MMFF94_VDW_GAMMA;
                double dterm2_drho = -1.12 * 7.0 * rho6 / (denom * denom);

                double dE_drho = epsilon * (dterm1_7_drho * term2 + term1_7 * dterm2_drho);
                double dE_dR = dE_drho / R_star;

                if (is_14_pair(mol, i, j)) {
                    dE_dR *= MMFF94_VDW_14_SCALE;
                }

                double g = dE_dR / R;
                gx[i] -= g * dx;
                gy[i] -= g * dy;
                gz[i] -= g * dz;
                gx[j] += g * dx;
                gy[j] += g * dy;
                gz[j] += g * dz;
            }
        }
    }

    return energy;
}

/* ============================================================================
 * Electrostatic Energy (Buffered Coulomb)
 * E = 332.0716 * qi*qj / (D*(R+delta))
 * ============================================================================ */

double mmff94_calc_elec_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                const mol_coords_t* coords,
                                double* gx, double* gy, double* gz) {
    double energy = 0.0;
    int n = mol->num_atoms;

    for (int i = 0; i < n; i++) {
        double qi = ctx->atom_data[i].partial_charge;
        if (fabs(qi) < 1e-8) continue;

        for (int j = i + 1; j < n; j++) {
            double qj = ctx->atom_data[j].partial_charge;
            if (fabs(qj) < 1e-8) continue;

            /* Skip bonded pairs (1-2) */
            if (molecule_get_bond_between_const(mol, i, j)) continue;

            /* Skip 1-3 pairs */
            if (is_13_pair(mol, i, j)) continue;

            /* Calculate distance */
            double dx = coords->coords_3d[j].x - coords->coords_3d[i].x;
            double dy = coords->coords_3d[j].y - coords->coords_3d[i].y;
            double dz = coords->coords_3d[j].z - coords->coords_3d[i].z;
            double R = sqrt(dx*dx + dy*dy + dz*dz);
            if (R < 0.1) R = 0.1;

            /* Buffered Coulomb */
            double E = MMFF94_ELEC_CONST * qi * qj /
                      (MMFF94_DIELECTRIC * (R + MMFF94_ELEC_BUFFER));

            /* Apply 1-4 scaling */
            bool is_14 = is_14_pair(mol, i, j);
            if (is_14) {
                E *= MMFF94_ELEC_14_SCALE;
            }

            energy += E;

            /* Gradient */
            if (gx) {
                double dE_dR = -E / (R + MMFF94_ELEC_BUFFER);
                double g = dE_dR / R;

                gx[i] -= g * dx;
                gy[i] -= g * dy;
                gz[i] -= g * dz;
                gx[j] += g * dx;
                gy[j] += g * dy;
                gz[j] += g * dz;
            }
        }
    }

    return energy;
}

/* ============================================================================
 * Total MMFF94 Energy
 * ============================================================================ */

double mmff94_calc_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                          const mol_coords_t* coords,
                          double* grad_x, double* grad_y, double* grad_z) {
    if (!mol || !ctx || !coords || !ctx->types_assigned) {
        return 0.0;
    }

    int n = mol->num_atoms;

    /* Initialize gradients to zero */
    if (grad_x) {
        memset(grad_x, 0, n * sizeof(double));
        memset(grad_y, 0, n * sizeof(double));
        memset(grad_z, 0, n * sizeof(double));
    }

    double energy = 0.0;

    /* Sum all energy components */
    energy += mmff94_calc_bond_energy(mol, ctx, coords, grad_x, grad_y, grad_z);
    energy += mmff94_calc_angle_energy(mol, ctx, coords, grad_x, grad_y, grad_z);
    energy += mmff94_calc_strbnd_energy(mol, ctx, coords, grad_x, grad_y, grad_z);
    energy += mmff94_calc_oop_energy(mol, ctx, coords, grad_x, grad_y, grad_z);
    energy += mmff94_calc_torsion_energy(mol, ctx, coords, grad_x, grad_y, grad_z);
    energy += mmff94_calc_vdw_energy(mol, ctx, coords, grad_x, grad_y, grad_z);

    /* Only include electrostatics if charges are computed */
    if (ctx->charges_computed) {
        energy += mmff94_calc_elec_energy(mol, ctx, coords, grad_x, grad_y, grad_z);
    }

    return energy;
}
