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
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Pre-computed Element Properties (Single Lookup Optimization)
 * ============================================================================ */

/* All properties for a single element - avoids 10+ separate table lookups */
typedef struct {
    double mw;          /* atomic weight */
    double en;          /* electronegativity */
    double rcov;        /* covalent radius */
    double vdw;         /* VdW volume */
    double polz;        /* polarizability */
    double ie;          /* ionization energy */
    double ea;          /* electron affinity */
    int8_t ox;          /* max oxidation state */
    int8_t val;         /* valence electrons */
    int8_t group;       /* periodic group */
    uint8_t flags;      /* bit flags: 0x01=halogen, 0x02=hetero, 0x04=heavy_z, 0x08=metalloid */
} elem_props_t;

#define ELEM_FLAG_HALOGEN   0x01
#define ELEM_FLAG_HETERO    0x02
#define ELEM_FLAG_HEAVY_Z   0x04
#define ELEM_FLAG_METALLOID 0x08

#define MAX_ELEM_IDX 128

/* Flags shorthand for readability */
#define F_HAL   ELEM_FLAG_HALOGEN
#define F_HET   ELEM_FLAG_HETERO
#define F_HVY   ELEM_FLAG_HEAVY_Z
#define F_MET   ELEM_FLAG_METALLOID

/* Static element property table - all 118 elements
 * Data sources: IUPAC 2021 (MW), Pauling (EN), Cordero 2008 (rcov),
 * Bondi (VdW radii -> volume), CRC Handbook (polz), NIST (IE, EA)
 * VdW volume = 4/3 * pi * r^3 in Angstrom^3 */
static const elem_props_t g_elem_props_data[MAX_ELEM_IDX] = {
    /* 0: Unknown - use carbon defaults */
    [0] = {12.011, 2.55, 0.77, 20.58, 1.76, 11.26, 1.26, 4, 4, 14, 0},

    /* Period 1 */
    [ELEM_H]  = {1.008,   2.20, 0.31,  7.24, 0.667, 13.60,  0.75, 1, 1, 1,  0},
    [ELEM_He] = {4.003,   0.00, 0.28,  8.80, 0.205, 24.59,  0.00, 0, 2, 18, F_HET},

    /* Period 2 */
    [ELEM_Li] = {6.941,   0.98, 1.28, 22.28, 24.3,   5.39,  0.62, 1, 1, 1,  F_HET},
    [ELEM_Be] = {9.012,   1.57, 0.96, 14.14, 5.60,   9.32,  0.00, 2, 2, 2,  F_HET},
    [ELEM_B]  = {10.811,  2.04, 0.84, 17.87, 3.03,   8.30,  0.28, 3, 3, 13, F_HET | F_MET},
    [ELEM_C]  = {12.011,  2.55, 0.77, 20.58, 1.76,  11.26,  1.26, 4, 4, 14, 0},
    [ELEM_N]  = {14.007,  3.04, 0.71, 15.60, 1.10,  14.53, -0.07, 5, 5, 15, F_HET},
    [ELEM_O]  = {15.999,  3.44, 0.66, 14.71, 0.802, 13.62,  1.46, 2, 6, 16, F_HET},
    [ELEM_F]  = {18.998,  3.98, 0.57, 13.31, 0.557, 17.42,  3.40, 1, 7, 17, F_HAL | F_HET},
    [ELEM_Ne] = {20.180,  0.00, 0.58, 12.33, 0.396, 21.56,  0.00, 0, 8, 18, F_HET | F_HVY},

    /* Period 3 */
    [ELEM_Na] = {22.990,  0.93, 1.66, 38.79, 24.1,   5.14,  0.55, 1, 1, 1,  F_HET | F_HVY},
    [ELEM_Mg] = {24.305,  1.31, 1.41, 25.08, 10.6,   7.65,  0.00, 2, 2, 2,  F_HET | F_HVY},
    [ELEM_Al] = {26.982,  1.61, 1.21, 32.52, 6.80,   5.99,  0.43, 3, 3, 13, F_HET | F_HVY},
    [ELEM_Si] = {28.086,  1.90, 1.11, 38.79, 5.38,   8.15,  1.39, 4, 4, 14, F_HET | F_HVY | F_MET},
    [ELEM_P]  = {30.974,  2.19, 1.07, 24.43, 3.63,  10.49,  0.75, 5, 5, 15, F_HET | F_HVY},
    [ELEM_S]  = {32.065,  2.58, 1.05, 24.43, 2.90,  10.36,  2.08, 6, 6, 16, F_HET | F_HVY},
    [ELEM_Cl] = {35.453,  3.16, 1.02, 22.45, 2.18,  12.97,  3.61, 7, 7, 17, F_HAL | F_HET | F_HVY},
    [ELEM_Ar] = {39.948,  0.00, 1.06, 23.59, 1.64,  15.76,  0.00, 0, 8, 18, F_HET | F_HVY},

    /* Period 4 */
    [ELEM_K]  = {39.098,  0.82, 2.03, 76.73, 43.4,   4.34,  0.50, 1, 1, 1,  F_HET | F_HVY},
    [ELEM_Ca] = {40.078,  1.00, 1.76, 53.81, 22.8,   6.11,  0.02, 2, 2, 2,  F_HET | F_HVY},
    [ELEM_Sc] = {44.956,  1.36, 1.70, 45.70, 17.8,   6.56,  0.19, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Ti] = {47.867,  1.54, 1.60, 39.95, 14.6,   6.83,  0.08, 4, 2, 4,  F_HET | F_HVY},
    [ELEM_V]  = {50.942,  1.63, 1.53, 35.16, 12.4,   6.75,  0.53, 5, 2, 5,  F_HET | F_HVY},
    [ELEM_Cr] = {51.996,  1.66, 1.39, 28.73, 11.6,   6.77,  0.68, 6, 1, 6,  F_HET | F_HVY},
    [ELEM_Mn] = {54.938,  1.55, 1.39, 28.73, 9.4,    7.43,  0.00, 7, 2, 7,  F_HET | F_HVY},
    [ELEM_Fe] = {55.845,  1.83, 1.32, 25.72, 8.4,    7.90,  0.15, 6, 2, 8,  F_HET | F_HVY},
    [ELEM_Co] = {58.933,  1.88, 1.26, 23.02, 7.5,    7.88,  0.66, 5, 2, 9,  F_HET | F_HVY},
    [ELEM_Ni] = {58.693,  1.91, 1.24, 22.10, 6.8,    7.64,  1.16, 4, 2, 10, F_HET | F_HVY},
    [ELEM_Cu] = {63.546,  1.90, 1.32, 25.72, 6.2,    7.73,  1.24, 4, 1, 11, F_HET | F_HVY},
    [ELEM_Zn] = {65.380,  1.65, 1.22, 21.19, 5.75,   9.39,  0.00, 2, 2, 12, F_HET | F_HVY},
    [ELEM_Ga] = {69.723,  1.81, 1.22, 32.52, 8.12,   6.00,  0.43, 3, 3, 13, F_HET | F_HVY},
    [ELEM_Ge] = {72.640,  2.01, 1.20, 36.73, 6.07,   7.90,  1.23, 4, 4, 14, F_HET | F_HVY | F_MET},
    [ELEM_As] = {74.922,  2.18, 1.19, 26.52, 4.31,   9.79,  0.81, 5, 5, 15, F_HET | F_HVY | F_MET},
    [ELEM_Se] = {78.960,  2.55, 1.20, 28.73, 3.77,   9.75,  2.02, 6, 6, 16, F_HET | F_HVY | F_MET},
    [ELEM_Br] = {79.904,  2.96, 1.20, 26.52, 3.05,  11.81,  3.36, 7, 7, 17, F_HAL | F_HET | F_HVY},
    [ELEM_Kr] = {83.798,  3.00, 1.16, 25.72, 2.48,  14.00,  0.00, 2, 8, 18, F_HET | F_HVY},

    /* Period 5 */
    [ELEM_Rb] = {85.468,  0.82, 2.20, 90.46, 47.3,   4.18,  0.49, 1, 1, 1,  F_HET | F_HVY},
    [ELEM_Sr] = {87.620,  0.95, 1.95, 67.85, 27.6,   5.69,  0.05, 2, 2, 2,  F_HET | F_HVY},
    [ELEM_Y]  = {88.906,  1.22, 1.90, 59.19, 22.7,   6.22,  0.31, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Zr] = {91.224,  1.33, 1.75, 48.31, 17.9,   6.63,  0.43, 4, 2, 4,  F_HET | F_HVY},
    [ELEM_Nb] = {92.906,  1.60, 1.64, 40.74, 15.7,   6.76,  0.89, 5, 1, 5,  F_HET | F_HVY},
    [ELEM_Mo] = {95.960,  2.16, 1.54, 35.16, 12.8,   7.09,  0.75, 6, 1, 6,  F_HET | F_HVY},
    [ELEM_Tc] = {98.000,  1.90, 1.47, 31.07, 11.4,   7.28,  0.55, 7, 2, 7,  F_HET | F_HVY},
    [ELEM_Ru] = {101.07,  2.20, 1.46, 30.32, 9.6,    7.36,  1.05, 8, 1, 8,  F_HET | F_HVY},
    [ELEM_Rh] = {102.91,  2.28, 1.42, 27.48, 8.6,    7.46,  1.14, 6, 1, 9,  F_HET | F_HVY},
    [ELEM_Pd] = {106.42,  2.20, 1.39, 28.73, 4.8,    8.34,  0.56, 4, 0, 10, F_HET | F_HVY},
    [ELEM_Ag] = {107.87,  1.93, 1.45, 29.58, 7.2,    7.58,  1.30, 3, 1, 11, F_HET | F_HVY},
    [ELEM_Cd] = {112.41,  1.69, 1.44, 28.73, 7.36,   8.99,  0.00, 2, 2, 12, F_HET | F_HVY},
    [ELEM_In] = {114.82,  1.78, 1.42, 39.95, 10.2,   5.79,  0.40, 3, 3, 13, F_HET | F_HVY},
    [ELEM_Sn] = {118.71,  1.96, 1.39, 41.78, 7.7,    7.34,  1.11, 4, 4, 14, F_HET | F_HVY},
    [ELEM_Sb] = {121.76,  2.05, 1.39, 38.79, 6.6,    8.61,  1.05, 5, 5, 15, F_HET | F_HVY | F_MET},
    [ELEM_Te] = {127.60,  2.10, 1.38, 36.73, 5.5,    9.01,  1.97, 6, 6, 16, F_HET | F_HVY | F_MET},
    [ELEM_I]  = {126.90,  2.66, 1.39, 32.52, 4.70,  10.45,  3.06, 7, 7, 17, F_HAL | F_HET | F_HVY},
    [ELEM_Xe] = {131.29,  2.60, 1.40, 35.16, 4.04,  12.13,  0.00, 8, 8, 18, F_HET | F_HVY},

    /* Period 6 */
    [ELEM_Cs] = {132.91,  0.79, 2.44, 117.1, 59.6,   3.89,  0.47, 1, 1, 1,  F_HET | F_HVY},
    [ELEM_Ba] = {137.33,  0.89, 2.15, 81.70, 39.7,   5.21,  0.14, 2, 2, 2,  F_HET | F_HVY},
    [ELEM_La] = {138.91,  1.10, 2.07, 70.97, 31.1,   5.58,  0.47, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Ce] = {140.12,  1.12, 2.04, 67.85, 29.6,   5.54,  0.50, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Pr] = {140.91,  1.13, 2.03, 66.55, 28.2,   5.47,  0.50, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Nd] = {144.24,  1.14, 2.01, 64.00, 31.4,   5.53,  0.50, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Pm] = {145.00,  1.13, 1.99, 61.52, 30.0,   5.58,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Sm] = {150.36,  1.17, 1.98, 60.30, 28.8,   5.64,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Eu] = {151.96,  1.20, 1.98, 60.30, 27.7,   5.67,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Gd] = {157.25,  1.20, 1.96, 57.92, 23.5,   6.15,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Tb] = {158.93,  1.20, 1.94, 55.60, 25.5,   5.86,  0.50, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Dy] = {162.50,  1.22, 1.92, 53.34, 24.5,   5.94,  0.50, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Ho] = {164.93,  1.23, 1.92, 53.34, 23.6,   6.02,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Er] = {167.26,  1.24, 1.89, 50.12, 22.7,   6.11,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Tm] = {168.93,  1.25, 1.90, 51.13, 21.8,   6.18,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Yb] = {173.05,  1.10, 1.87, 47.03, 21.0,   6.25,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Lu] = {174.97,  1.27, 1.87, 47.03, 21.9,   5.43,  0.50, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Hf] = {178.49,  1.30, 1.75, 48.31, 16.2,   6.83,  0.00, 4, 2, 4,  F_HET | F_HVY},
    [ELEM_Ta] = {180.95,  1.50, 1.70, 45.70, 13.1,   7.55,  0.32, 5, 2, 5,  F_HET | F_HVY},
    [ELEM_W]  = {183.84,  2.36, 1.62, 39.38, 11.1,   7.86,  0.82, 6, 2, 6,  F_HET | F_HVY},
    [ELEM_Re] = {186.21,  1.90, 1.51, 33.51, 9.7,    7.83,  0.15, 7, 2, 7,  F_HET | F_HVY},
    [ELEM_Os] = {190.23,  2.20, 1.44, 28.73, 8.5,    8.44,  1.10, 8, 2, 8,  F_HET | F_HVY},
    [ELEM_Ir] = {192.22,  2.20, 1.41, 26.79, 7.6,    8.97,  1.56, 6, 2, 9,  F_HET | F_HVY},
    [ELEM_Pt] = {195.08,  2.28, 1.36, 24.61, 6.5,    8.96,  2.13, 6, 1, 10, F_HET | F_HVY},
    [ELEM_Au] = {196.97,  2.54, 1.36, 24.61, 5.8,    9.23,  2.31, 5, 1, 11, F_HET | F_HVY},
    [ELEM_Hg] = {200.59,  2.00, 1.32, 25.72, 5.02,  10.44,  0.00, 4, 2, 12, F_HET | F_HVY},
    [ELEM_Tl] = {204.38,  1.62, 1.45, 32.52, 7.6,    6.11,  0.38, 3, 3, 13, F_HET | F_HVY},
    [ELEM_Pb] = {207.20,  2.33, 1.46, 34.32, 6.8,    7.42,  0.36, 4, 4, 14, F_HET | F_HVY},
    [ELEM_Bi] = {208.98,  2.02, 1.48, 36.73, 7.4,    7.29,  0.94, 5, 5, 15, F_HET | F_HVY},
    [ELEM_Po] = {209.00,  2.00, 1.40, 31.82, 6.8,    8.41,  1.90, 6, 6, 16, F_HET | F_HVY | F_MET},
    [ELEM_At] = {210.00,  2.20, 1.50, 35.16, 6.0,    9.30,  2.80, 7, 7, 17, F_HAL | F_HET | F_HVY},
    [ELEM_Rn] = {222.00,  2.20, 1.50, 38.79, 5.3,   10.75,  0.00, 2, 8, 18, F_HET | F_HVY},

    /* Period 7 */
    [ELEM_Fr] = {223.00,  0.70, 2.60, 136.3, 48.6,   4.07,  0.47, 1, 1, 1,  F_HET | F_HVY},
    [ELEM_Ra] = {226.00,  0.90, 2.21, 87.11, 38.3,   5.28,  0.10, 2, 2, 2,  F_HET | F_HVY},
    [ELEM_Ac] = {227.00,  1.10, 2.15, 81.70, 32.1,   5.17,  0.35, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Th] = {232.04,  1.30, 2.06, 68.45, 32.1,   6.31,  0.00, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Pa] = {231.04,  1.50, 2.00, 62.75, 25.9,   5.89,  0.00, 5, 2, 3,  F_HET | F_HVY},
    [ELEM_U]  = {238.03,  1.38, 1.96, 57.92, 24.9,   6.19,  0.00, 6, 2, 3,  F_HET | F_HVY},
    [ELEM_Np] = {237.00,  1.36, 1.90, 59.19, 24.8,   6.27,  0.00, 7, 2, 3,  F_HET | F_HVY},
    [ELEM_Pu] = {244.00,  1.28, 1.87, 55.60, 24.5,   6.03,  0.00, 7, 2, 3,  F_HET | F_HVY},
    [ELEM_Am] = {243.00,  1.30, 1.80, 50.97, 23.3,   5.97,  0.00, 6, 2, 3,  F_HET | F_HVY},
    [ELEM_Cm] = {247.00,  1.30, 1.69, 43.89, 23.0,   5.99,  0.00, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Bk] = {247.00,  1.30, 1.68, 43.18, 22.7,   6.20,  0.00, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Cf] = {251.00,  1.30, 1.68, 43.18, 20.5,   6.28,  0.00, 4, 2, 3,  F_HET | F_HVY},
    [ELEM_Es] = {252.00,  1.30, 1.65, 40.74, 19.7,   6.42,  0.00, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Fm] = {257.00,  1.30, 1.67, 42.48, 23.8,   6.50,  0.00, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Md] = {258.00,  1.30, 1.73, 46.48, 18.2,   6.58,  0.00, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_No] = {259.00,  1.30, 1.76, 48.93, 17.5,   6.65,  0.00, 3, 2, 3,  F_HET | F_HVY},
    [ELEM_Lr] = {262.00,  1.30, 1.61, 38.82, 16.8,   4.90,  0.00, 3, 2, 3,  F_HET | F_HVY},

    /* Period 7 d-block and beyond (superheavy elements - estimated values) */
    [ELEM_Rf] = {267.00,  1.30, 1.57, 36.46, 16.0,   6.00,  0.00, 4, 2, 4,  F_HET | F_HVY},
    [ELEM_Db] = {268.00,  1.30, 1.49, 31.82, 15.0,   6.80,  0.00, 5, 2, 5,  F_HET | F_HVY},
    [ELEM_Sg] = {269.00,  1.30, 1.43, 28.08, 14.0,   7.80,  0.00, 6, 2, 6,  F_HET | F_HVY},
    [ELEM_Bh] = {270.00,  1.30, 1.41, 26.79, 13.0,   7.70,  0.00, 7, 2, 7,  F_HET | F_HVY},
    [ELEM_Hs] = {269.00,  1.30, 1.34, 23.98, 12.0,   7.60,  0.00, 8, 2, 8,  F_HET | F_HVY},
    [ELEM_Mt] = {278.00,  1.30, 1.29, 21.67, 11.0,   8.50,  0.00, 6, 2, 9,  F_HET | F_HVY},
    [ELEM_Ds] = {281.00,  1.30, 1.28, 21.19, 10.0,   9.50,  0.00, 6, 2, 10, F_HET | F_HVY},
    [ELEM_Rg] = {282.00,  1.30, 1.21, 18.58, 9.0,   10.60,  0.00, 5, 1, 11, F_HET | F_HVY},
    [ELEM_Cn] = {285.00,  1.30, 1.22, 18.95, 8.0,   11.90,  0.00, 4, 2, 12, F_HET | F_HVY},
    [ELEM_Nh] = {286.00,  1.30, 1.36, 24.61, 7.4,    7.30,  0.00, 3, 3, 13, F_HET | F_HVY},
    [ELEM_Fl] = {289.00,  1.30, 1.43, 28.08, 7.0,    8.50,  0.00, 4, 4, 14, F_HET | F_HVY},
    [ELEM_Mc] = {290.00,  1.30, 1.62, 39.38, 6.5,    5.60,  0.00, 3, 5, 15, F_HET | F_HVY},
    [ELEM_Lv] = {293.00,  1.30, 1.75, 48.31, 6.0,    6.60,  0.00, 4, 6, 16, F_HET | F_HVY},
    [ELEM_Ts] = {294.00,  1.30, 1.65, 40.74, 5.5,    7.70,  0.00, 5, 7, 17, F_HAL | F_HET | F_HVY},
    [ELEM_Og] = {294.00,  1.30, 1.57, 36.46, 5.0,    8.90,  0.00, 6, 8, 18, F_HET | F_HVY},
};

#undef F_HAL
#undef F_HET
#undef F_HVY
#undef F_MET

/* Mutable copy for runtime (allows init check pattern) */
static elem_props_t g_elem_props[MAX_ELEM_IDX];
static bool g_elem_props_init = false;

static void init_elem_props(void) {
    if (g_elem_props_init) return;
    memcpy(g_elem_props, g_elem_props_data, sizeof(g_elem_props));
    g_elem_props_init = true;
}

/* Fast property lookup - single bounds check, returns pointer to all props */
static inline const elem_props_t* get_elem_props(element_t elem) {
    if (elem <= 0 || elem >= MAX_ELEM_IDX) return &g_elem_props[ELEM_C];
    return &g_elem_props[elem];
}

/* Check if element is a heavy atom (Z > 9, i.e., beyond first period main group) */
static inline bool is_heavy_z(element_t elem) {
    if (!g_elem_props_init) init_elem_props();
    if (elem <= 0 || elem >= MAX_ELEM_IDX) return false;
    return (g_elem_props[elem].flags & ELEM_FLAG_HEAVY_Z) != 0;
}

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

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

    /* Ensure property table is initialized */
    if (!g_elem_props_init) init_elem_props();

    const int n_atoms = mol->num_atoms;
    const int n_bonds = mol->num_bonds;

    /* Pre-compute hybridization for all atoms once (avoids repeated bond-type lookups) */
    int8_t hyb_cache[512];  /* Stack-allocated, max 512 atoms */
    const int max_cached = (n_atoms < 512) ? n_atoms : 512;
    for (int i = 0; i < max_cached; i++) {
        hyb_cache[i] = (int8_t)get_hybridization(mol, i);
    }

    /* First pass: collect basic counts and sums */
    for (int i = 0; i < n_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        element_t e = a->element;

        if (e == ELEM_H) continue;  /* Skip explicit H for most calculations */

        s->n_atoms++;
        s->n_heavy++;

        /* Single property lookup instead of 10+ separate lookups */
        const elem_props_t* p = get_elem_props(e);
        double mw = p->mw;
        double en = p->en;
        double rcov = p->rcov;
        double vdw = p->vdw;
        double polz = p->polz;
        double ie = p->ie;
        double ea = p->ea;
        int ox = p->ox;
        int val = p->val;
        int group = p->group;
        uint8_t flags = p->flags;

        s->total_mw += mw;
        s->sum_en += en;
        s->sum_en_mw += en * mw;
        s->sum_vdw_mw += vdw * mw;
        s->sum_rcov_mw += rcov * mw;

        /* Use flags for fast classification (single bit test vs function call) */
        bool is_hetero = (flags & ELEM_FLAG_HETERO) != 0;
        bool is_heavy = (flags & ELEM_FLAG_HEAVY_Z) != 0;
        bool is_metal = (flags & ELEM_FLAG_METALLOID) != 0;
        (void)flags; /* Avoid unused warning - flags used via is_* bools */

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
            default:      if (is_hetero) { s->mw_hetero += mw; s->n_hetero++; } break;
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
        if (is_heavy && ie < 11.0) s->n_low_ie_heavy++;

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

        /* Hybridization (use cached value) */
        int hyb = (i < max_cached) ? hyb_cache[i] : get_hybridization(mol, i);
        if (hyb == 1) s->n_sp++;
        else if (hyb == 2) s->n_sp2++;
        else if (hyb == 3) s->n_sp3++;

        /* Metalloid */
        if (is_metal) s->n_metalloid++;

        /* Group 16 */
        if (group == 16) s->n_group16++;

        /* Mixed property counts (use cached flags) */
        if (is_hetero && en > 2.55) s->n_hetero_polar++;
        if (is_heavy && hyb == 3) s->n_sp3_heavy++;
        if (hyb == 2 && en > 2.55) s->n_sp2_en_above++;
        if (!is_hetero && ea > 1.0) s->n_nonhetero_high_ea++;
        if (is_hetero && polz < 1.5) s->n_hetero_low_polz++;
        if (hyb == 1 && en > 2.5) s->n_sp_high_en++;
        if (a->ring_count > 0 && ox > 2) s->n_ring_high_ox++;
        if (is_heavy && a->charge != 0) s->n_heavy_formal++;
        if (hyb == 3 && polz > 3.5) s->n_sp3_high_polz++;
        if (is_heavy && ox < 0) s->n_heavy_ox_neg++;
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

    /* Second pass: EN above/below average (use cached property lookup) */
    for (int i = 0; i < n_atoms; i++) {
        const atom_t* a = &mol->atoms[i];
        if (a->element == ELEM_H) continue;
        const elem_props_t* p = get_elem_props(a->element);
        if (p->en > s->avg_en) s->n_en_above_avg++;
        else s->n_en_below_avg++;
    }

    /* Bond pass */
    s->n_bonds = n_bonds;
    for (int i = 0; i < n_bonds; i++) {
        const bond_t* b = &mol->bonds[i];
        element_t e1 = mol->atoms[b->atom1].element;
        element_t e2 = mol->atoms[b->atom2].element;

        if (e1 == ELEM_H || e2 == ELEM_H) continue;

        /* Use cached property lookups for both atoms */
        const elem_props_t* p1 = get_elem_props(e1);
        const elem_props_t* p2 = get_elem_props(e2);
        double en1 = p1->en;
        double en2 = p2->en;
        double en_diff = fabs(en1 - en2);

        /* Polar vs unpolar bonds */
        if (en_diff < 0.4) s->n_unpolar_bonds++;
        else s->n_polar_bonds++;

        /* Both atoms EN > 3.0 */
        if (en1 > 3.0 && en2 > 3.0) s->n_en_bonded++;

        /* C-C bond types (use cached hybridization) */
        if (e1 == ELEM_C && e2 == ELEM_C) {
            s->n_cc_total++;
            int hyb1 = (b->atom1 < max_cached) ? hyb_cache[b->atom1] : get_hybridization(mol, b->atom1);
            int hyb2 = (b->atom2 < max_cached) ? hyb_cache[b->atom2] : get_hybridization(mol, b->atom2);
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

        /* Heavy atom polar bonds (use flags) */
        bool is_heavy = ((p1->flags | p2->flags) & ELEM_FLAG_HEAVY_Z) != 0;
        if (is_heavy && en_diff >= 0.4) {
            s->n_heavy_polar_bond++;
        }

        /* Halogen polar bonds (use flags) */
        bool is_halo = ((p1->flags | p2->flags) & ELEM_FLAG_HALOGEN) != 0;
        if (is_halo && en_diff >= 0.4) {
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
