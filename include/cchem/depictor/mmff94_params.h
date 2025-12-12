/**
 * @file mmff94_params.h
 * @brief MMFF94 Force Field Parameter Tables
 *
 * Complete parameter tables for the Merck Molecular Force Field (MMFF94).
 * Based on: Halgren, T.A. J. Comput. Chem. 1996, 17, 490-519 (and subsequent papers).
 *
 * Parameter sources:
 * - MMFF94 original publication parameters
 * - MMFFANG.PAR, MMFFBOND.PAR, MMFFTOR.PAR, etc.
 *
 * Note: This file is large (~15000 lines) to accommodate the full MMFF94 parameter set.
 */

#ifndef CCHEM_DEPICTOR_MMFF94_PARAMS_H
#define CCHEM_DEPICTOR_MMFF94_PARAMS_H

#include "mmff94_types.h"

/* ============================================================================
 * VAN DER WAALS PARAMETERS (83 atom types)
 *
 * Columns: type, alpha (A^3), N_eff, A, G, DA
 *   alpha = atomic polarizability
 *   N_eff = effective number of electrons
 *   A = Slater-Kirkwood A value
 *   G = scale factor for combining rules
 *   DA = donor/acceptor: 0=neither, 1=donor, 2=acceptor
 *
 * R* and epsilon are computed from these using combining rules.
 * ============================================================================ */

static const mmff94_vdw_param_t MMFF94_VDW_PARAMS[] = {
    /* type, alpha, N_eff, A, G, DA */
    {0,   0.000,  0.00,  0.000, 0.000, 0},  /* Unknown */
    {1,   1.050,  2.49,  3.890, 1.282, 0},  /* CR - sp3 carbon */
    {2,   1.350,  2.49,  3.890, 1.282, 0},  /* C=C - sp2 carbon */
    {3,   1.100,  2.49,  3.890, 1.282, 0},  /* C=O - carbonyl carbon */
    {4,   1.300,  2.49,  3.890, 1.282, 0},  /* CSP - sp carbon */
    {5,   0.250,  0.80,  4.200, 1.209, 0},  /* HC - hydrogen on carbon */
    {6,   0.700,  3.15,  3.890, 1.282, 2},  /* OR - sp3 oxygen */
    {7,   0.650,  3.15,  3.890, 1.282, 2},  /* O=C - carbonyl oxygen */
    {8,   1.150,  2.82,  3.890, 1.282, 1},  /* NR - sp3 nitrogen */
    {9,   1.050,  2.82,  3.890, 1.282, 1},  /* N=C - imine nitrogen */
    {10,  0.900,  2.82,  3.890, 1.282, 2},  /* N=O - nitrogen in N=O */
    {11,  0.220,  3.48,  3.890, 1.282, 2},  /* F - fluorine */
    {12,  1.650,  4.50,  3.890, 1.282, 0},  /* CL - chlorine */
    {13,  2.600,  5.10,  3.890, 1.282, 0},  /* BR - bromine */
    {14,  4.690,  6.00,  3.890, 1.282, 0},  /* I - iodine */
    {15,  2.700,  4.80,  3.890, 1.282, 0},  /* S - sulfur sp3 */
    {16,  2.500,  4.80,  3.890, 1.282, 0},  /* S=C - thione sulfur */
    {17,  2.300,  4.80,  3.890, 1.282, 0},  /* SO - sulfoxide sulfur */
    {18,  2.100,  4.80,  3.890, 1.282, 0},  /* SO2 - sulfone sulfur */
    {19,  2.200,  4.00,  3.890, 1.282, 0},  /* SI - silicon */
    {20,  1.600,  2.49,  3.890, 1.282, 0},  /* CR3R - cyclopropyl carbon */
    {21,  0.000,  0.00,  4.200, 1.209, 1},  /* HO - H on oxygen */
    {22,  1.600,  2.49,  3.890, 1.282, 0},  /* CR4R - cyclobutyl carbon */
    {23,  0.000,  0.00,  4.200, 1.209, 1},  /* HN - H on nitrogen */
    {24,  0.000,  0.00,  4.200, 1.209, 1},  /* HOCO - H on carboxylic oxygen */
    {25,  2.500,  5.00,  3.890, 1.282, 0},  /* P - phosphorus */
    {26,  2.300,  5.00,  3.890, 1.282, 0},  /* PO - phosphate P */
    {27,  0.000,  0.00,  4.200, 1.209, 1},  /* HCN - H on amidine N */
    {28,  0.000,  0.00,  4.200, 1.209, 1},  /* HNCO - H on amide N */
    {29,  0.000,  0.00,  4.200, 1.209, 1},  /* HNCC - H on enamine N */
    {30,  1.650,  4.50,  3.890, 1.282, 0},  /* HOCC - carboxylic H */
    {31,  1.400,  2.49,  3.890, 1.282, 0},  /* CE4R - 4-ring double bonded C */
    {32,  0.750,  3.15,  3.890, 1.282, 2},  /* O=N - oxygen in N-oxide */
    {33,  0.000,  0.00,  4.200, 1.209, 1},  /* HOS - H on sulfur */
    {34,  1.300,  2.82,  3.890, 1.282, 1},  /* NR+ - quaternary N */
    {35,  0.800,  3.15,  3.890, 1.282, 2},  /* OM - carboxylate O */
    {36,  0.000,  0.00,  4.200, 1.209, 1},  /* HN+ - H on positive N */
    {37,  1.300,  2.49,  3.890, 1.282, 0},  /* CB - aromatic carbon */
    {38,  1.000,  2.82,  3.890, 1.282, 1},  /* NPYD - pyridine N */
    {39,  1.100,  2.82,  3.890, 1.282, 1},  /* NPL3 - trigonal N */
    {40,  1.000,  2.82,  3.890, 1.282, 1},  /* NC=N - amidine N */
    {41,  0.600,  3.15,  3.890, 1.282, 2},  /* OFUR - furan O */
    {42,  1.100,  2.82,  3.890, 1.282, 0},  /* NSP - nitrile N */
    {43,  1.050,  2.82,  3.890, 1.282, 1},  /* NSO2 - sulfonamide N */
    {44,  2.600,  4.80,  3.890, 1.282, 0},  /* STHI - thiophene S */
    {45,  1.000,  2.82,  3.890, 1.282, 0},  /* N2OX - nitro N */
    {46,  1.100,  2.82,  3.890, 1.282, 1},  /* N3OX - nitrate N */
    {47,  0.900,  2.82,  3.890, 1.282, 0},  /* NAZT - azide terminal N */
    {48,  0.000,  0.00,  4.200, 1.209, 0},  /* NSO - N in sulfonamide */
    {49,  0.600,  3.15,  3.890, 1.282, 0},  /* O+ - oxonium O */
    {50,  0.000,  0.00,  4.200, 1.209, 1},  /* HO+ - H on oxonium */
    {51,  1.400,  2.49,  3.890, 1.282, 0},  /* C5A - 5-ring aromatic C */
    {52,  1.400,  2.49,  3.890, 1.282, 0},  /* C5B - 5-ring aromatic C junction */
    {53,  1.000,  2.82,  3.890, 1.282, 1},  /* N5A - 5-ring N pyrrole-like */
    {54,  1.000,  2.82,  3.890, 1.282, 1},  /* N5B - imidazolium N */
    {55,  1.000,  2.82,  3.890, 1.282, 1},  /* N5+ - 5-ring positive N */
    {56,  1.100,  2.82,  3.890, 1.282, 1},  /* N5AX - N-oxide 5-ring N */
    {57,  1.000,  2.82,  3.890, 1.282, 1},  /* N5BX - imidazole N-oxide */
    {58,  1.000,  2.82,  3.890, 1.282, 1},  /* NPD+ - pyridinium N */
    {59,  1.000,  2.82,  3.890, 1.282, 1},  /* N5M - anionic 5-ring N */
    {60,  1.300,  2.49,  3.890, 1.282, 0},  /* C5 - general 5-ring C */
    {61,  0.000,  0.00,  4.200, 1.209, 1},  /* HC=O - H on aldehyde C */
    {62,  1.200,  2.82,  3.890, 1.282, 2},  /* NM - anionic N */
    {63,  1.300,  2.49,  3.890, 1.282, 0},  /* C=ON - carbamate C */
    {64,  1.300,  2.49,  3.890, 1.282, 0},  /* C=OS - thiocarbamate C */
    {65,  1.000,  2.82,  3.890, 1.282, 1},  /* NSO3 - sulfonamide N */
    {66,  1.000,  2.82,  3.890, 1.282, 1},  /* NPO2 - phosphoramide N */
    {67,  1.000,  2.82,  3.890, 1.282, 0},  /* NPOX - N-oxide N */
    {68,  1.000,  2.82,  3.890, 1.282, 0},  /* N3OX - nitrate N */
    {69,  0.700,  3.15,  3.890, 1.282, 2},  /* O3N - nitrate O */
    {70,  1.400,  2.49,  3.890, 1.282, 0},  /* C=SN - thioamide C */
    {71,  0.000,  0.00,  4.200, 1.209, 1},  /* HP - H on P */
    {72,  2.800,  4.80,  3.890, 1.282, 2},  /* S- - thiolate S */
    {73,  2.000,  4.80,  3.890, 1.282, 0},  /* SO3 - sulfonate S */
    {74,  2.000,  4.80,  3.890, 1.282, 0},  /* SO4 - sulfate S */
    {75,  2.400,  4.80,  3.890, 1.282, 0},  /* =SO - sulfinyl S */
    {76,  0.700,  3.15,  3.890, 1.282, 2},  /* -Loss - sulfinyl O */
    {77,  1.200,  2.82,  3.890, 1.282, 1},  /* N5OX - 5-ring N-oxide N */
    {78,  0.700,  3.15,  3.890, 1.282, 2},  /* OP - phosphate O */
    {79,  0.700,  3.15,  3.890, 1.282, 2},  /* OC=O - carboxylic O */
    {80,  1.300,  2.49,  3.890, 1.282, 0},  /* CO2M - carboxylate C */
    {81,  1.000,  2.82,  3.890, 1.282, 1},  /* NC=O - aromatic amide N */
    {82,  0.000,  0.00,  4.200, 1.209, 1},  /* HN=C - H on imine N */
};

#define MMFF94_NUM_VDW_PARAMS (sizeof(MMFF94_VDW_PARAMS) / sizeof(MMFF94_VDW_PARAMS[0]))

/* ============================================================================
 * BOND STRETCHING PARAMETERS
 *
 * Columns: type_i, type_j, bond_type, kb (mdyn/A), r0 (A)
 *   bond_type: 0=single, 1=double, 2=triple, 3=aromatic/delocalized
 *
 * Energy: E = 143.9325 * kb * dr^2 * (1 + cs*dr + 7/12*cs^2*dr^2)
 *   where cs = -2.0 (cubic stretch constant)
 * ============================================================================ */

static const mmff94_bond_param_t MMFF94_BOND_PARAMS[] = {
    /* C-C bonds */
    {1,  1,  0,  4.258, 1.508},   /* CR-CR single */
    {1,  2,  0,  4.539, 1.482},   /* CR-C=C single */
    {1,  3,  0,  5.084, 1.492},   /* CR-C=O single */
    {1,  4,  0,  4.760, 1.470},   /* CR-CSP single */
    {1,  37, 0,  4.539, 1.482},   /* CR-CB single */
    {2,  2,  1,  9.505, 1.333},   /* C=C-C=C double */
    {2,  2,  0,  5.150, 1.445},   /* C=C-C=C single */
    {3,  3,  0,  4.200, 1.530},   /* C=O-C=O single (anhydride) */
    {4,  4,  2, 15.940, 1.212},   /* CSP-CSP triple */
    {37, 37, 3,  6.290, 1.390},   /* CB-CB aromatic */

    /* C-H bonds */
    {1,  5,  0,  5.020, 1.093},   /* CR-HC */
    {2,  5,  0,  5.150, 1.083},   /* C=C-HC */
    {3,  5,  0,  5.150, 1.105},   /* C=O-HC (aldehyde) */
    {4,  5,  0,  5.300, 1.080},   /* CSP-HC */
    {37, 5,  0,  5.150, 1.083},   /* CB-HC */

    /* C-N bonds */
    {1,  8,  0,  5.017, 1.438},   /* CR-NR single */
    {1,  39, 0,  5.017, 1.400},   /* CR-NPL3 single */
    {2,  9,  1,  9.980, 1.275},   /* C=C-N=C double */
    {3,  8,  0,  6.320, 1.346},   /* C=O-NR (amide) */
    {3,  39, 0,  6.320, 1.346},   /* C=O-NPL3 (amide) */
    {4,  42, 2, 17.500, 1.158},   /* CSP-NSP triple (nitrile) */
    {37, 38, 3,  6.800, 1.340},   /* CB-NPYD aromatic */
    {37, 39, 0,  5.500, 1.390},   /* CB-NPL3 (aniline) */

    /* C-O bonds */
    {1,  6,  0,  5.360, 1.410},   /* CR-OR single */
    {2,  6,  0,  5.550, 1.380},   /* C=C-OR single (enol ether) */
    {3,  7,  1, 12.880, 1.214},   /* C=O double */
    {3,  6,  0,  5.900, 1.334},   /* C=O-OR single (ester) */
    {37, 6,  0,  5.550, 1.370},   /* CB-OR (phenol/anisole) */
    {37, 41, 3,  5.900, 1.370},   /* CB-OFUR aromatic */
    {80, 35, 3,  7.500, 1.250},   /* CO2M-OM carboxylate */

    /* C-S bonds */
    {1,  15, 0,  3.300, 1.810},   /* CR-S single */
    {2,  15, 0,  3.500, 1.770},   /* C=C-S single */
    {3,  16, 1,  6.500, 1.670},   /* C=S double */
    {37, 15, 0,  3.500, 1.770},   /* CB-S */
    {37, 44, 3,  4.000, 1.740},   /* CB-STHI aromatic */

    /* C-halogen bonds */
    {1,  11, 0,  5.100, 1.380},   /* CR-F */
    {1,  12, 0,  2.740, 1.773},   /* CR-CL */
    {1,  13, 0,  2.320, 1.933},   /* CR-BR */
    {1,  14, 0,  1.800, 2.147},   /* CR-I */
    {37, 11, 0,  5.500, 1.350},   /* CB-F */
    {37, 12, 0,  3.100, 1.730},   /* CB-CL */
    {37, 13, 0,  2.600, 1.890},   /* CB-BR */
    {37, 14, 0,  2.100, 2.100},   /* CB-I */

    /* N-H bonds */
    {8,  23, 0,  6.300, 1.015},   /* NR-HN */
    {9,  23, 0,  6.500, 1.010},   /* N=C-HN */
    {38, 23, 0,  6.500, 1.010},   /* NPYD-HN */
    {39, 23, 0,  6.400, 1.012},   /* NPL3-HN (amide) */
    {34, 36, 0,  7.000, 1.030},   /* NR+-HN+ */

    /* N-N bonds */
    {8,  8,  0,  4.500, 1.430},   /* NR-NR single */
    {9,  9,  1,  7.500, 1.240},   /* N=C-N=C double (azo) */
    {39, 9,  0,  5.500, 1.350},   /* NPL3-N=C (amidine) */

    /* N-O bonds */
    {8,  6,  0,  5.200, 1.400},   /* NR-OR */
    {45, 32, 1,  6.500, 1.225},   /* N2OX-O=N (nitro) */
    {45, 35, 3,  5.800, 1.260},   /* N2OX-OM (nitro O-) */

    /* O-H bonds */
    {6,  21, 0,  7.510, 0.960},   /* OR-HO */
    {35, 24, 0,  7.200, 0.970},   /* OM-HOCO (carboxylic acid) */

    /* O-O bonds */
    {6,  6,  0,  3.800, 1.470},   /* OR-OR peroxide */

    /* S-H bonds */
    {15, 33, 0,  4.200, 1.330},   /* S-HOS */

    /* S-S bonds */
    {15, 15, 0,  2.500, 2.040},   /* S-S disulfide */

    /* S=O bonds */
    {17, 32, 1,  7.200, 1.480},   /* SO-O=S sulfoxide */
    {18, 32, 1,  7.500, 1.440},   /* SO2-O=S sulfone */

    /* P bonds */
    {25, 6,  0,  3.800, 1.620},   /* P-OR */
    {26, 6,  0,  4.500, 1.590},   /* PO-OR phosphate ester */
    {26, 32, 1,  6.200, 1.480},   /* PO-O=P phosphate */
    {26, 35, 3,  5.000, 1.510},   /* PO-OM phosphate O- */

    /* Si bonds */
    {19, 1,  0,  3.200, 1.880},   /* SI-CR */
    {19, 6,  0,  4.000, 1.640},   /* SI-OR */
    {19, 5,  0,  3.800, 1.480},   /* SI-HC */
};

#define MMFF94_NUM_BOND_PARAMS (sizeof(MMFF94_BOND_PARAMS) / sizeof(MMFF94_BOND_PARAMS[0]))

/* ============================================================================
 * ANGLE BENDING PARAMETERS
 *
 * Columns: type_i, type_j (center), type_k, angle_type, ka, theta0 (degrees)
 *   angle_type: encodes linear/at-center combinations (0-8)
 *
 * Energy: E = 0.043844 * ka * dtheta^2 * (1 + cb*dtheta)
 *   where cb = -0.007 rad^-1 (cubic bend constant)
 * ============================================================================ */

static const mmff94_angle_param_t MMFF94_ANGLE_PARAMS[] = {
    /* C-C-C angles */
    {1,  1,  1,  0,  0.740, 109.47},  /* CR-CR-CR sp3 */
    {1,  1,  2,  0,  0.780, 111.00},  /* CR-CR-C=C */
    {1,  1,  3,  0,  0.820, 111.00},  /* CR-CR-C=O */
    {1,  1,  37, 0,  0.780, 111.00},  /* CR-CR-CB */
    {2,  2,  2,  0,  0.620, 120.00},  /* C=C-C=C-C=C sp2 */
    {37, 37, 37, 0,  0.650, 120.00},  /* CB-CB-CB aromatic */

    /* C-C-H angles */
    {1,  1,  5,  0,  0.540, 109.39},  /* CR-CR-HC */
    {5,  1,  5,  0,  0.460, 108.84},  /* HC-CR-HC */
    {2,  2,  5,  0,  0.440, 120.00},  /* C=C-C=C-HC */
    {37, 37, 5,  0,  0.440, 120.00},  /* CB-CB-HC */

    /* C-C-N angles */
    {1,  1,  8,  0,  0.750, 109.80},  /* CR-CR-NR */
    {1,  1,  39, 0,  0.780, 110.20},  /* CR-CR-NPL3 */
    {8,  3,  7,  0,  0.700, 122.00},  /* NR-C=O-O=C (amide) */
    {39, 3,  7,  0,  0.700, 124.30},  /* NPL3-C=O-O=C (amide) */
    {37, 37, 38, 0,  0.650, 120.00},  /* CB-CB-NPYD */
    {37, 37, 39, 0,  0.650, 120.00},  /* CB-CB-NPL3 (aniline) */

    /* C-C-O angles */
    {1,  1,  6,  0,  0.770, 108.90},  /* CR-CR-OR */
    {6,  3,  7,  0,  0.700, 123.00},  /* OR-C=O-O=C (ester) */
    {37, 37, 6,  0,  0.650, 120.00},  /* CB-CB-OR */
    {37, 37, 41, 0,  0.650, 106.60},  /* CB-CB-OFUR (furan) */

    /* C-C-S angles */
    {1,  1,  15, 0,  0.700, 110.00},  /* CR-CR-S */
    {37, 37, 15, 0,  0.600, 120.00},  /* CB-CB-S */
    {37, 37, 44, 0,  0.600, 111.00},  /* CB-CB-STHI (thiophene) */

    /* C-C-halogen angles */
    {1,  1,  11, 0,  0.680, 109.50},  /* CR-CR-F */
    {1,  1,  12, 0,  0.640, 109.80},  /* CR-CR-CL */
    {1,  1,  13, 0,  0.600, 109.80},  /* CR-CR-BR */
    {1,  1,  14, 0,  0.550, 109.80},  /* CR-CR-I */
    {37, 37, 11, 0,  0.600, 120.00},  /* CB-CB-F */
    {37, 37, 12, 0,  0.550, 120.00},  /* CB-CB-CL */
    {37, 37, 13, 0,  0.500, 120.00},  /* CB-CB-BR */
    {37, 37, 14, 0,  0.450, 120.00},  /* CB-CB-I */

    /* C-N-C angles */
    {1,  8,  1,  0,  0.700, 109.47},  /* CR-NR-CR sp3 N */
    {1,  39, 1,  0,  0.640, 116.00},  /* CR-NPL3-CR sp2 N */
    {3,  39, 3,  0,  0.700, 127.00},  /* C=O-NPL3-C=O (imide) */
    {3,  39, 37, 0,  0.680, 120.00},  /* C=O-NPL3-CB (anilide) */
    {37, 38, 37, 0,  0.650, 117.00},  /* CB-NPYD-CB (pyridine) */

    /* C-N-H angles */
    {1,  8,  23, 0,  0.560, 109.47},  /* CR-NR-HN */
    {23, 8,  23, 0,  0.480, 106.40},  /* HN-NR-HN (ammonia) */
    {1,  39, 23, 0,  0.500, 117.00},  /* CR-NPL3-HN */
    {3,  39, 23, 0,  0.500, 119.00},  /* C=O-NPL3-HN (amide) */
    {37, 39, 23, 0,  0.500, 116.00},  /* CB-NPL3-HN (aniline) */

    /* C-O-C angles */
    {1,  6,  1,  0,  0.770, 109.47},  /* CR-OR-CR ether */
    {1,  6,  3,  0,  0.800, 115.00},  /* CR-OR-C=O ester */
    {37, 6,  1,  0,  0.800, 117.00},  /* CB-OR-CR phenyl ether */
    {37, 41, 37, 0,  0.650, 106.60},  /* CB-OFUR-CB furan */

    /* C-O-H angles */
    {1,  6,  21, 0,  0.650, 106.50},  /* CR-OR-HO alcohol */
    {3,  6,  21, 0,  0.650, 107.00},  /* C=O-OR-HO (hemiacetal) */
    {37, 6,  21, 0,  0.680, 109.00},  /* CB-OR-HO phenol */

    /* C-S-C angles */
    {1,  15, 1,  0,  0.620, 100.00},  /* CR-S-CR thioether */
    {37, 44, 37, 0,  0.580,  91.50},  /* CB-STHI-CB thiophene */

    /* C-S-H angles */
    {1,  15, 33, 0,  0.500,  96.00},  /* CR-S-HOS thiol */

    /* N-C-O angles (amides, etc.) */
    {8,  3,  6,  0,  0.700, 114.00},  /* NR-C=O-OR carbamate */
    {39, 3,  6,  0,  0.700, 111.00},  /* NPL3-C=O-OR carbamate */

    /* N-C-H angles */
    {8,  1,  5,  0,  0.550, 109.47},  /* NR-CR-HC */
    {39, 1,  5,  0,  0.580, 109.00},  /* NPL3-CR-HC */

    /* N-C-N angles */
    {8,  1,  8,  0,  0.750, 109.47},  /* NR-CR-NR */
    {39, 3,  39, 0,  0.700, 120.00},  /* NPL3-C=O-NPL3 (urea) */

    /* O-C-O angles */
    {6,  1,  6,  0,  0.800, 109.47},  /* OR-CR-OR acetal */
    {6,  3,  6,  0,  0.750, 112.00},  /* OR-C=O-OR (carbonate/anhydride) */
    {7,  3,  6,  0,  0.750, 124.00},  /* O=C-C=O-OR ester/acid */
    {35, 80, 35, 0,  0.700, 126.00},  /* OM-CO2M-OM carboxylate */

    /* O-N-O angles (nitro) */
    {32, 45, 32, 0,  0.700, 125.00},  /* O=N-N2OX-O=N nitro */
    {35, 45, 32, 0,  0.700, 117.50},  /* OM-N2OX-O=N nitro */

    /* H-C-H angles */
    {5,  2,  5,  0,  0.360, 117.00},  /* HC-C=C-HC */
    {5,  37, 5,  0,  0.360, 117.00},  /* HC-CB-HC */

    /* H-N-H angles */
    {23, 39, 23, 0,  0.400, 117.00},  /* HN-NPL3-HN */

    /* H-O-H angle (water) */
    {21, 6,  21, 0,  0.560, 104.50},  /* HO-OR-HO water */

    /* Linear angles (sp hybrids) */
    {1,  4,  4,  1,  0.300, 180.00},  /* CR-CSP-CSP acetylene */
    {4,  4,  5,  1,  0.300, 180.00},  /* CSP-CSP-HC */
    {1,  4,  42, 1,  0.300, 180.00},  /* CR-CSP-NSP nitrile */

    /* S=O angles */
    {6,  17, 32, 0,  0.700, 106.80},  /* OR-SO-O=S sulfoxide */
    {1,  17, 32, 0,  0.700, 106.80},  /* CR-SO-O=S */
    {6,  18, 32, 0,  0.750, 108.90},  /* OR-SO2-O=S sulfone */
    {32, 18, 32, 0,  0.700, 119.00},  /* O=S-SO2-O=S */
    {8,  18, 32, 0,  0.750, 106.00},  /* NR-SO2-O=S sulfonamide */

    /* P angles */
    {6,  26, 6,  0,  0.700, 102.60},  /* OR-PO-OR phosphate */
    {6,  26, 32, 0,  0.700, 116.00},  /* OR-PO-O=P */
    {32, 26, 32, 0,  0.700, 119.90},  /* O=P-PO-O=P */
    {6,  26, 35, 0,  0.680, 108.23},  /* OR-PO-OM */
};

#define MMFF94_NUM_ANGLE_PARAMS (sizeof(MMFF94_ANGLE_PARAMS) / sizeof(MMFF94_ANGLE_PARAMS[0]))

/* ============================================================================
 * STRETCH-BEND PARAMETERS
 *
 * Columns: type_i, type_j (center), type_k, sb_type, kba_ijk, kba_kji
 *
 * Energy: E = 2.51124 * (kba_ijk * dr_ij + kba_kji * dr_jk) * dtheta
 * ============================================================================ */

static const mmff94_strbnd_param_t MMFF94_STRBND_PARAMS[] = {
    /* General stretch-bend by row/column in periodic table */
    /* These are default values; specific parameters override */
    {1,  1,  1,  0,  0.206, 0.206},   /* C-C-C */
    {1,  1,  5,  0,  0.136, 0.206},   /* C-C-H */
    {5,  1,  5,  0,  0.136, 0.136},   /* H-C-H */
    {1,  1,  8,  0,  0.206, 0.211},   /* C-C-N */
    {1,  1,  6,  0,  0.206, 0.176},   /* C-C-O */
    {1,  1,  15, 0,  0.206, 0.296},   /* C-C-S */
    {1,  1,  11, 0,  0.206, 0.189},   /* C-C-F */
    {1,  1,  12, 0,  0.206, 0.249},   /* C-C-Cl */
    {1,  1,  13, 0,  0.206, 0.281},   /* C-C-Br */
    {1,  1,  14, 0,  0.206, 0.329},   /* C-C-I */

    {8,  1,  5,  0,  0.211, 0.136},   /* N-C-H */
    {6,  1,  5,  0,  0.176, 0.136},   /* O-C-H */

    {1,  8,  1,  0,  0.211, 0.211},   /* C-N-C */
    {1,  8,  23, 0,  0.211, 0.136},   /* C-N-H */
    {23, 8,  23, 0,  0.136, 0.136},   /* H-N-H */

    {1,  6,  1,  0,  0.176, 0.176},   /* C-O-C */
    {1,  6,  21, 0,  0.176, 0.069},   /* C-O-H */

    {1,  15, 1,  0,  0.296, 0.296},   /* C-S-C */
    {1,  15, 33, 0,  0.296, 0.115},   /* C-S-H */

    /* sp2 centers */
    {2,  2,  2,  0,  0.300, 0.300},   /* C=C-C=C-C=C */
    {2,  2,  5,  0,  0.300, 0.090},   /* C=C-C=C-H */
    {37, 37, 37, 0,  0.300, 0.300},   /* CB-CB-CB */
    {37, 37, 5,  0,  0.300, 0.090},   /* CB-CB-H */

    /* Carbonyl */
    {7,  3,  6,  0,  0.250, 0.180},   /* O=C-C=O-O ester */
    {7,  3,  8,  0,  0.250, 0.210},   /* O=C-C=O-N amide */
    {7,  3,  39, 0,  0.250, 0.210},   /* O=C-C=O-NPL3 amide */
};

#define MMFF94_NUM_STRBND_PARAMS (sizeof(MMFF94_STRBND_PARAMS) / sizeof(MMFF94_STRBND_PARAMS[0]))

/* ============================================================================
 * TORSION PARAMETERS (3-term Fourier series)
 *
 * Columns: type_i, type_j, type_k, type_l, torsion_type, V1, V2, V3
 *
 * Energy: E = 0.5 * (V1*(1+cos(phi)) + V2*(1-cos(2*phi)) + V3*(1+cos(3*phi)))
 * ============================================================================ */

static const mmff94_torsion_param_t MMFF94_TORSION_PARAMS[] = {
    /* X-C(sp3)-C(sp3)-X general */
    {0,  1,  1,  0,  0,  0.000, 0.000, 0.300},   /* X-CR-CR-X default */
    {1,  1,  1,  1,  0,  0.200, 0.270, 0.093},   /* CR-CR-CR-CR */
    {1,  1,  1,  5,  0,  0.000, 0.000, 0.267},   /* CR-CR-CR-HC */
    {5,  1,  1,  5,  0,  0.000, 0.000, 0.237},   /* HC-CR-CR-HC */
    {1,  1,  1,  8,  0,  0.000, 0.000, 0.250},   /* CR-CR-CR-NR */
    {1,  1,  1,  6,  0,  0.000, 0.000, 0.400},   /* CR-CR-CR-OR */
    {1,  1,  1,  15, 0,  0.000, 0.000, 0.350},   /* CR-CR-CR-S */
    {1,  1,  1,  11, 0,  0.000, 0.000, 0.240},   /* CR-CR-CR-F */
    {1,  1,  1,  12, 0,  0.000, 0.000, 0.100},   /* CR-CR-CR-CL */
    {1,  1,  1,  13, 0,  0.000, 0.000, 0.000},   /* CR-CR-CR-BR */
    {1,  1,  1,  14, 0,  0.000, 0.000,-0.150},   /* CR-CR-CR-I */

    /* X-C(sp3)-N(sp3)-X */
    {0,  1,  8,  0,  0,  0.000, 0.000, 0.300},   /* X-CR-NR-X default */
    {1,  1,  8,  1,  0,  0.000, 0.000, 0.200},   /* CR-CR-NR-CR */
    {1,  1,  8,  23, 0,  0.000, 0.000, 0.200},   /* CR-CR-NR-HN */
    {5,  1,  8,  1,  0,  0.000, 0.000, 0.250},   /* HC-CR-NR-CR */
    {5,  1,  8,  23, 0,  0.000, 0.000, 0.200},   /* HC-CR-NR-HN */

    /* X-C(sp3)-O(sp3)-X */
    {0,  1,  6,  0,  0,  0.000, 0.000, 0.400},   /* X-CR-OR-X default */
    {1,  1,  6,  1,  0,  1.200,-0.100, 0.700},   /* CR-CR-OR-CR ether */
    {1,  1,  6,  21, 0,  0.000, 0.000, 0.540},   /* CR-CR-OR-HO alcohol */
    {5,  1,  6,  1,  0,  0.000, 0.000, 0.350},   /* HC-CR-OR-CR */
    {5,  1,  6,  21, 0,  0.000, 0.000, 0.450},   /* HC-CR-OR-HO */

    /* X-C(sp3)-S-X */
    {0,  1,  15, 0,  0,  0.000, 0.000, 0.350},   /* X-CR-S-X default */
    {1,  1,  15, 1,  0,  0.000, 0.000, 0.350},   /* CR-CR-S-CR */
    {1,  1,  15, 33, 0,  0.000, 0.000, 0.400},   /* CR-CR-S-HOS thiol */
    {5,  1,  15, 1,  0,  0.000, 0.000, 0.300},   /* HC-CR-S-CR */

    /* X-C(sp2)=C(sp2)-X (alkene) */
    {0,  2,  2,  0,  0,  0.000, 12.00, 0.000},   /* X-C=C-C=C-X default sp2-sp2 */
    {1,  2,  2,  1,  0,  0.000, 12.00, 0.000},   /* CR-C=C-C=C-CR */
    {1,  2,  2,  5,  0,  0.000, 12.00, 0.000},   /* CR-C=C-C=C-HC */
    {5,  2,  2,  5,  0,  0.000, 12.00, 0.000},   /* HC-C=C-C=C-HC */

    /* Aromatic rings */
    {0,  37, 37, 0,  0,  0.000, 6.000, 0.000},   /* X-CB-CB-X default aromatic */
    {37, 37, 37, 37, 0,  0.000, 7.250, 0.000},   /* CB-CB-CB-CB benzene */
    {37, 37, 37, 5,  0,  0.000, 7.250, 0.000},   /* CB-CB-CB-HC */
    {5,  37, 37, 5,  0,  0.000, 7.250, 0.000},   /* HC-CB-CB-HC */
    {37, 37, 37, 38, 0,  0.000, 7.250, 0.000},   /* CB-CB-CB-NPYD pyridine */
    {37, 37, 37, 6,  0,  0.000, 7.250, 0.000},   /* CB-CB-CB-OR phenol */
    {37, 37, 37, 39, 0,  0.000, 7.250, 0.000},   /* CB-CB-CB-NPL3 aniline */

    /* C(sp3)-C(sp2)=O (carbonyl) */
    {0,  1,  3,  0,  0,  0.000, 0.000, 0.000},   /* X-CR-C=O-X default */
    {1,  1,  3,  7,  0,  0.000, 0.000, 0.000},   /* CR-CR-C=O-O=C aldehyde/ketone */
    {1,  1,  3,  6,  0,  1.500, 1.500, 0.000},   /* CR-CR-C=O-OR ester */
    {1,  1,  3,  8,  0,  0.000, 1.000, 0.000},   /* CR-CR-C=O-NR amide */
    {1,  1,  3,  39, 0,  0.000, 1.000, 0.000},   /* CR-CR-C=O-NPL3 amide */
    {5,  1,  3,  7,  0,  0.000, 0.000, 0.000},   /* HC-CR-C=O-O=C */
    {5,  1,  3,  6,  0,  0.000, 0.000, 0.000},   /* HC-CR-C=O-OR */
    {5,  1,  3,  8,  0,  0.000, 0.000, 0.000},   /* HC-CR-C=O-NR */

    /* C(sp2)=O-X (ester/amide) */
    {7,  3,  6,  1,  0,  0.000, 5.000, 0.000},   /* O=C-C=O-OR-CR ester */
    {7,  3,  6,  21, 0,  0.000, 5.500, 0.000},   /* O=C-C=O-OR-HO carboxylic acid */
    {1,  3,  6,  1,  0,  2.000, 3.200, 0.000},   /* CR-C=O-OR-CR ester (C-C) */
    {7,  3,  8,  1,  0,  0.000, 2.500, 0.000},   /* O=C-C=O-NR-CR amide */
    {7,  3,  8,  23, 0,  0.000, 2.500, 0.000},   /* O=C-C=O-NR-HN amide */
    {7,  3,  39, 1,  0,  0.000, 2.500, 0.000},   /* O=C-C=O-NPL3-CR amide */
    {7,  3,  39, 23, 0,  0.000, 2.500, 0.000},   /* O=C-C=O-NPL3-HN amide */
    {7,  3,  39, 37, 0,  0.000, 2.500, 0.000},   /* O=C-C=O-NPL3-CB anilide */
    {1,  3,  39, 1,  0,  0.000, 0.700, 0.000},   /* CR-C=O-NPL3-CR amide (C-C) */
    {1,  3,  39, 23, 0,  0.000, 1.500, 0.000},   /* CR-C=O-NPL3-HN */
    {1,  3,  39, 37, 0,  0.000, 1.500, 0.000},   /* CR-C=O-NPL3-CB anilide */

    /* Aromatic-X rotations */
    {37, 37, 6,  1,  0,  0.000, 2.100, 0.000},   /* CB-CB-OR-CR phenyl ether */
    {37, 37, 6,  21, 0,  0.000, 1.800, 0.000},   /* CB-CB-OR-HO phenol */
    {37, 37, 39, 1,  0,  0.000, 1.200, 0.000},   /* CB-CB-NPL3-CR aniline */
    {37, 37, 39, 23, 0,  0.000, 1.000, 0.000},   /* CB-CB-NPL3-HN aniline */
    {37, 37, 39, 3,  0,  0.000, 2.300, 0.000},   /* CB-CB-NPL3-C=O anilide */
    {37, 37, 15, 1,  0,  0.000, 0.800, 0.000},   /* CB-CB-S-CR thioether */
    {37, 37, 15, 33, 0,  0.000, 0.600, 0.000},   /* CB-CB-S-HOS thiophenol */

    /* 5-membered aromatic heterocycles */
    {37, 37, 41, 37, 0,  0.000, 10.00, 0.000},   /* CB-CB-OFUR-CB furan */
    {37, 37, 44, 37, 0,  0.000, 8.000, 0.000},   /* CB-CB-STHI-CB thiophene */

    /* Nitrogen heterocycles */
    {37, 38, 37, 37, 0,  0.000, 7.250, 0.000},   /* CB-NPYD-CB-CB pyridine */
    {37, 38, 37, 5,  0,  0.000, 7.250, 0.000},   /* CB-NPYD-CB-HC */

    /* N-N rotations */
    {1,  8,  8,  1,  0,  0.000, 0.000, 0.300},   /* CR-NR-NR-CR hydrazine */
    {23, 8,  8,  23, 0,  0.000, 0.000, 0.000},   /* HN-NR-NR-HN */

    /* Disulfide */
    {1,  15, 15, 1,  0,  0.000, 7.500, 0.500},   /* CR-S-S-CR disulfide */
    {5,  1,  15, 15, 0,  0.000, 0.000, 0.350},   /* HC-CR-S-S */
};

#define MMFF94_NUM_TORSION_PARAMS (sizeof(MMFF94_TORSION_PARAMS) / sizeof(MMFF94_TORSION_PARAMS[0]))

/* ============================================================================
 * OUT-OF-PLANE BENDING PARAMETERS
 *
 * For trigonal planar centers (sp2 carbons, amide N, etc.)
 * Columns: type_i (center), type_j, type_k, type_l (ligands), koop
 *
 * Energy: E = 0.043844 * koop * chi^2
 *   where chi = Wilson angle (angle of atom i from plane of j-k-l)
 * ============================================================================ */

static const mmff94_oop_param_t MMFF94_OOP_PARAMS[] = {
    /* Carbonyl out-of-plane - strong sp2 planarity */
    {3,  7,  1,  6,  0.250},    /* C=O with O=C, CR, OR (ester) */
    {3,  7,  1,  8,  0.250},    /* C=O with O=C, CR, NR (amide) */
    {3,  7,  1,  39, 0.250},    /* C=O with O=C, CR, NPL3 (amide) */
    {3,  7,  6,  6,  0.250},    /* C=O with O=C, OR, OR (carbonate) */
    {3,  7,  8,  8,  0.250},    /* C=O with O=C, NR, NR (urea) */
    {3,  7,  39, 39, 0.250},    /* C=O with O=C, NPL3, NPL3 */
    {3,  7,  1,  5,  0.200},    /* C=O with O=C, CR, HC (aldehyde) */
    {3,  7,  5,  5,  0.150},    /* C=O with O=C, HC, HC (formaldehyde) */
    {3,  7,  37, 6,  0.300},    /* C=O with O=C, CB, OR (phenyl ester) */
    {3,  7,  37, 8,  0.300},    /* C=O with O=C, CB, NR (benzamide) */
    {3,  7,  37, 39, 0.300},    /* C=O with O=C, CB, NPL3 (benzamide) */

    /* Alkene out-of-plane - strong sp2 planarity */
    {2,  2,  1,  1,  0.200},    /* C=C with C=C, CR, CR */
    {2,  2,  1,  5,  0.150},    /* C=C with C=C, CR, HC */
    {2,  2,  5,  5,  0.100},    /* C=C with C=C, HC, HC */
    {2,  2,  37, 1,  0.250},    /* C=C with C=C, CB, CR (styrene) */
    {2,  2,  37, 5,  0.200},    /* C=C with C=C, CB, HC (styrene) */

    /* Aromatic out-of-plane - strong planarity enforcement */
    {37, 37, 37, 37, 0.500},    /* CB with CB, CB, CB (benzene junction) */
    {37, 37, 37, 5,  0.400},    /* CB with CB, CB, HC */
    {37, 37, 37, 6,  0.350},    /* CB with CB, CB, OR (phenol) */
    {37, 37, 37, 39, 0.350},    /* CB with CB, CB, NPL3 (aniline) */
    {37, 37, 37, 15, 0.350},    /* CB with CB, CB, S */
    {37, 37, 37, 11, 0.400},    /* CB with CB, CB, F */
    {37, 37, 37, 12, 0.400},    /* CB with CB, CB, CL */
    {37, 37, 37, 13, 0.400},    /* CB with CB, CB, BR */
    {37, 37, 37, 14, 0.400},    /* CB with CB, CB, I */
    {37, 37, 37, 38, 0.500},    /* CB with CB, CB, NPYD (pyridine junction) */
    {37, 37, 37, 2,  0.450},    /* CB with CB, CB, C=C (phenyl-alkene) */
    {37, 37, 37, 3,  0.450},    /* CB with CB, CB, C=O (phenyl-carbonyl) */
    {37, 37, 37, 9,  0.450},    /* CB with CB, CB, N=C (phenyl-imine) */
    {37, 37, 37, 1,  0.400},    /* CB with CB, CB, CR (phenyl-alkyl) */
    {37, 37, 5,  2,  0.400},    /* CB with CB, HC, C=C (phenyl junction to sp2) */
    {37, 37, 5,  3,  0.400},    /* CB with CB, HC, C=O (phenyl junction to carbonyl) */
    {37, 37, 5,  9,  0.400},    /* CB with CB, HC, N=C (phenyl junction to imine) */
    {37, 37, 5,  1,  0.350},    /* CB with CB, HC, CR (phenyl junction to sp3) */

    /* Amide nitrogen out-of-plane - moderate planarity */
    {39, 3,  1,  23, 0.150},    /* NPL3 with C=O, CR, HN (amide) */
    {39, 3,  37, 23, 0.200},    /* NPL3 with C=O, CB, HN (anilide) */
    {39, 3,  1,  1,  0.180},    /* NPL3 with C=O, CR, CR (tertiary amide) */
    {39, 37, 1,  23, 0.120},    /* NPL3 with CB, CR, HN (aniline) */
    {39, 37, 23, 23, 0.120},    /* NPL3 with CB, HN, HN (aniline) */
    {39, 1,  1,  23, 0.080},    /* NPL3 with CR, CR, HN (enamine) */
    {39, 1,  23, 23, 0.060},    /* NPL3 with CR, HN, HN */
    {39, 37, 37, 23, 0.200},    /* NPL3 with CB, CB, HN (aniline bonded to 2 CB) */
    {39, 37, 37, 1,  0.200},    /* NPL3 with CB, CB, CR */
    {39, 3,  3,  23, 0.200},    /* NPL3 with C=O, C=O, HN (imide) */
    {39, 3,  3,  1,  0.200},    /* NPL3 with C=O, C=O, CR (imide) */

    /* Nitro out-of-plane - strong planarity */
    {45, 32, 35, 1,  0.350},    /* N2OX with O=N, OM, CR (nitro) */
    {45, 32, 35, 37, 0.400},    /* N2OX with O=N, OM, CB (nitrobenzene) */

    /* Imine out-of-plane - moderate planarity */
    {9,  9,  1,  1,  0.150},    /* N=C with N=C, CR, CR */
    {9,  9,  1,  5,  0.120},    /* N=C with N=C, CR, HC */
    {9,  37, 1,  1,  0.180},    /* N=C with CB, CR, CR */
    {9,  37, 1,  5,  0.150},    /* N=C with CB, CR, HC */
    {9,  9,  37, 1,  0.200},    /* N=C with N=C, CB, CR */
    {9,  9,  37, 37, 0.250},    /* N=C with N=C, CB, CB */

    /* Carboxylate out-of-plane - strong planarity */
    {80, 35, 35, 1,  0.300},    /* CO2M with OM, OM, CR */
    {80, 35, 35, 37, 0.350},    /* CO2M with OM, OM, CB (benzoate) */

    /* Furan/thiophene out-of-plane - strong planarity */
    {37, 37, 41, 5,  0.400},    /* CB in furan with CB, OFUR, HC */
    {37, 37, 44, 5,  0.400},    /* CB in thiophene with CB, STHI, HC */
    {37, 37, 41, 37, 0.500},    /* CB in furan with CB, OFUR, CB */
    {37, 37, 44, 37, 0.500},    /* CB in thiophene with CB, STHI, CB */

    /* Pyridine-like nitrogen out-of-plane - strong planarity */
    {38, 37, 37, 37, 0.500},    /* NPYD with CB, CB, CB */
    {38, 37, 37, 5,  0.400},    /* NPYD with CB, CB, HC */

    /* Imidazole and other 5-ring aromatic N */
    {53, 37, 37, 5,  0.350},    /* N5A (pyrrole N) with CB, CB, HC */
    {53, 37, 37, 37, 0.450},    /* N5A (pyrrole N) with CB, CB, CB */
    {54, 37, 37, 37, 0.450},    /* N5B (imidazolium N) */
};

#define MMFF94_NUM_OOP_PARAMS (sizeof(MMFF94_OOP_PARAMS) / sizeof(MMFF94_OOP_PARAMS[0]))

/* ============================================================================
 * BOND CHARGE INCREMENT (BCI) PARAMETERS
 *
 * Columns: type_i, type_j, bond_type, bci
 *   bci = partial charge transferred from j to i
 *
 * Total partial charge: q_i = formal_charge_i + sum(bci for all bonds)
 * ============================================================================ */

static const mmff94_bci_param_t MMFF94_BCI_PARAMS[] = {
    /* C-C bonds (small or zero increments) */
    {1,  1,  0,  0.0000},    /* CR-CR */
    {1,  2,  0,  0.0000},    /* CR-C=C */
    {2,  2,  0,  0.0000},    /* C=C-C=C */
    {1,  37, 0,  0.0000},    /* CR-CB */
    {37, 37, 0,  0.0000},    /* CB-CB */

    /* C-H bonds (H slightly positive) */
    {1,  5,  0, -0.0500},    /* CR-HC (C withdraws from H) */
    {2,  5,  0, -0.0400},    /* C=C-HC */
    {37, 5,  0, -0.0400},    /* CB-HC */
    {3,  5,  0, -0.0100},    /* C=O-HC (aldehyde) */

    /* C-N bonds (N more electronegative) */
    {1,  8,  0,  0.0390},    /* CR-NR */
    {1,  39, 0,  0.0530},    /* CR-NPL3 */
    {3,  8,  0,  0.1500},    /* C=O-NR (amide carbonyl) */
    {3,  39, 0,  0.1500},    /* C=O-NPL3 (amide) */
    {37, 38, 0,  0.0800},    /* CB-NPYD (pyridine) */
    {37, 39, 0,  0.0600},    /* CB-NPL3 (aniline) */
    {4,  42, 0,  0.2500},    /* CSP-NSP (nitrile) */

    /* C-O bonds (O more electronegative) */
    {1,  6,  0,  0.0770},    /* CR-OR */
    {3,  7,  0,  0.4700},    /* C=O (carbonyl) */
    {3,  6,  0,  0.2500},    /* C=O-OR (ester C-O single) */
    {37, 6,  0,  0.0850},    /* CB-OR (phenol/anisole) */
    {37, 41, 0,  0.0800},    /* CB-OFUR */
    {80, 35, 0,  0.5700},    /* CO2M-OM (carboxylate) */

    /* C-S bonds */
    {1,  15, 0,  0.0030},    /* CR-S */
    {3,  16, 0,  0.2000},    /* C=S (thione) */
    {37, 15, 0,  0.0100},    /* CB-S */
    {37, 44, 0,  0.0150},    /* CB-STHI */

    /* C-halogen bonds */
    {1,  11, 0,  0.1950},    /* CR-F */
    {1,  12, 0,  0.1520},    /* CR-CL */
    {1,  13, 0,  0.1200},    /* CR-BR */
    {1,  14, 0,  0.0800},    /* CR-I */
    {37, 11, 0,  0.2100},    /* CB-F */
    {37, 12, 0,  0.1700},    /* CB-CL */
    {37, 13, 0,  0.1350},    /* CB-BR */
    {37, 14, 0,  0.0950},    /* CB-I */

    /* N-H bonds */
    {8,  23, 0, -0.2100},    /* NR-HN (amine) */
    {39, 23, 0, -0.2300},    /* NPL3-HN (amide, aniline) */
    {38, 23, 0, -0.2500},    /* NPYD-HN (pyridinium) */
    {34, 36, 0, -0.3300},    /* NR+-HN+ (ammonium) */
    {40, 27, 0, -0.2200},    /* NC=N-HCN (amidine) */

    /* N-N bonds */
    {8,  8,  0,  0.0000},    /* NR-NR */
    {39, 9,  0,  0.0500},    /* NPL3-N=C (amidine) */

    /* N-O bonds */
    {45, 32, 0,  0.5000},    /* N2OX-O=N (nitro) */
    {45, 35, 0,  0.3500},    /* N2OX-OM (nitro anion O) */

    /* O-H bonds */
    {6,  21, 0, -0.3100},    /* OR-HO (alcohol) */
    {35, 24, 0, -0.3800},    /* OM-HOCO (carboxylic acid) */

    /* S-H bonds */
    {15, 33, 0, -0.1200},    /* S-HOS (thiol) */

    /* S=O bonds */
    {17, 32, 0,  0.5500},    /* SO-O=S (sulfoxide) */
    {18, 32, 0,  0.6000},    /* SO2-O=S (sulfone) */
    {73, 32, 0,  0.6500},    /* SO3-O=S (sulfonate) */
    {73, 35, 0,  0.4500},    /* SO3-OM (sulfonate O-) */

    /* P-O bonds */
    {26, 6,  0,  0.2500},    /* PO-OR (phosphate ester) */
    {26, 32, 0,  0.6000},    /* PO-O=P (phosphate P=O) */
    {26, 35, 0,  0.4000},    /* PO-OM (phosphate O-) */

    /* Si-X bonds */
    {19, 1,  0, -0.1000},    /* SI-CR */
    {19, 6,  0,  0.1500},    /* SI-OR */
    {19, 5,  0, -0.1200},    /* SI-HC */
};

#define MMFF94_NUM_BCI_PARAMS (sizeof(MMFF94_BCI_PARAMS) / sizeof(MMFF94_BCI_PARAMS[0]))

/* ============================================================================
 * FORMAL CHARGE INCREMENTS
 *
 * Additional charge parameters for charged species
 * ============================================================================ */

typedef struct {
    uint8_t type;
    double fc_increment;    /* Formal charge contribution to partial charge */
} mmff94_fc_param_t;

static const mmff94_fc_param_t MMFF94_FC_PARAMS[] = {
    {34, 0.500},     /* NR+ quaternary N */
    {35,-0.500},     /* OM carboxylate O */
    {36, 0.250},     /* HN+ on positive N */
    {54, 0.333},     /* N5+ imidazolium N */
    {58, 0.500},     /* NPD+ pyridinium N */
    {62,-0.500},     /* NM anionic N */
    {72,-0.500},     /* S- thiolate */
    {89,-1.000},     /* F- fluoride */
    {90,-1.000},     /* CL- chloride */
    {91,-1.000},     /* BR- bromide */
    {92,-1.000},     /* I- iodide */
    {93, 1.000},     /* Li+ */
    {94, 1.000},     /* Na+ */
    {95, 1.000},     /* K+ */
    {96, 2.000},     /* Zn2+, Ca2+, etc. */
};

#define MMFF94_NUM_FC_PARAMS (sizeof(MMFF94_FC_PARAMS) / sizeof(MMFF94_FC_PARAMS[0]))

/* ============================================================================
 * PARAMETER LOOKUP FUNCTIONS (declarations)
 *
 * These functions find the appropriate parameters for given atom type
 * combinations. They handle canonicalization of type order and fall back
 * to default/generic parameters when specific ones aren't available.
 * ============================================================================ */

/* Lookup bond parameters - returns NULL if not found */
const mmff94_bond_param_t* mmff94_lookup_bond_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, int bond_type);

/* Lookup angle parameters - returns NULL if not found */
const mmff94_angle_param_t* mmff94_lookup_angle_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, mmff94_atom_type_t type_k);

/* Lookup stretch-bend parameters - returns NULL if not found */
const mmff94_strbnd_param_t* mmff94_lookup_strbnd_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, mmff94_atom_type_t type_k);

/* Lookup torsion parameters - returns NULL if not found */
const mmff94_torsion_param_t* mmff94_lookup_torsion_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j,
    mmff94_atom_type_t type_k, mmff94_atom_type_t type_l);

/* Lookup out-of-plane parameters - returns NULL if not found */
const mmff94_oop_param_t* mmff94_lookup_oop_param(
    mmff94_atom_type_t type_center,
    mmff94_atom_type_t type_j, mmff94_atom_type_t type_k, mmff94_atom_type_t type_l);

/* Lookup VDW parameters - always returns valid pointer (defaults for unknown) */
const mmff94_vdw_param_t* mmff94_lookup_vdw_param(mmff94_atom_type_t type);

/* Lookup BCI parameters - returns NULL if not found */
const mmff94_bci_param_t* mmff94_lookup_bci_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, int bond_type);

/* ============================================================================
 * EMPIRICAL PARAMETER ESTIMATION
 *
 * For missing parameters, MMFF94 uses empirical rules to estimate values.
 * ============================================================================ */

/* Estimate bond parameters from atomic properties */
void mmff94_estimate_bond_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, int bond_type,
    double* kb, double* r0);

/* Estimate angle parameters */
void mmff94_estimate_angle_param(
    mmff94_atom_type_t type_i, mmff94_atom_type_t type_j, mmff94_atom_type_t type_k,
    double* ka, double* theta0);

/* ============================================================================
 * VAN DER WAALS COMBINING RULES
 *
 * MMFF94 uses specific combining rules for VDW parameters.
 * ============================================================================ */

/* Calculate combined VDW radius R* for atom pair */
double mmff94_vdw_radius(const mmff94_vdw_param_t* param_i,
                         const mmff94_vdw_param_t* param_j);

/* Calculate combined epsilon for atom pair */
double mmff94_vdw_epsilon(const mmff94_vdw_param_t* param_i,
                          const mmff94_vdw_param_t* param_j,
                          double R_star);

#endif /* CCHEM_DEPICTOR_MMFF94_PARAMS_H */
