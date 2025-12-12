/**
 * @file mmff94_types.h
 * @brief MMFF94 force field type definitions
 *
 * Merck Molecular Force Field atom types and data structures.
 * Based on: Halgren, T.A. J. Comput. Chem. 1996, 17, 490-519.
 */

#ifndef CCHEM_DEPICTOR_MMFF94_TYPES_H
#define CCHEM_DEPICTOR_MMFF94_TYPES_H

#include <stdint.h>
#include <stdbool.h>
#include "cchem/canonicalizer/molecule.h"
#include "cchem/depictor/types.h"

/* ============================================================================
 * MMFF94 Symbolic Atom Types (83 types)
 * ============================================================================
 * Naming convention:
 *   Element + hybridization/environment descriptor
 *   AR = aromatic, R = sp3, SP = sp, EQ = double bond (=), IM = imidazole, etc.
 * ============================================================================ */

typedef enum {
    MMFF94_TYPE_UNKNOWN = 0,

    /* Carbon types (1-10) */
    MMFF94_TYPE_CR      = 1,    /* Alkyl carbon, sp3 */
    MMFF94_TYPE_C_EQ_C  = 2,    /* Vinylic carbon, sp2 (C=C) */
    MMFF94_TYPE_C_EQ_O  = 3,    /* Generic carbonyl carbon, sp2 (C=O) */
    MMFF94_TYPE_CSP     = 4,    /* Acetylenic carbon, sp (C#C, C#N) */
    MMFF94_TYPE_C_AR    = 37,   /* Aromatic carbon */

    /* Nitrogen types (11-30) */
    MMFF94_TYPE_NR      = 8,    /* Amine nitrogen, sp3 */
    MMFF94_TYPE_N_EQ_C  = 9,    /* Imine nitrogen, sp2 (N=C) */
    MMFF94_TYPE_N_EQ_O  = 10,   /* Nitrogen in N=O */
    MMFF94_TYPE_NC_EQ_O = 10,   /* Nitrogen in -N=C=O (isocyanate) */
    MMFF94_TYPE_NSP     = 42,   /* Nitrogen in nitrile, sp (N#C) */
    MMFF94_TYPE_N_AR    = 38,   /* Aromatic nitrogen (pyridine-like) */
    MMFF94_TYPE_NAZT    = 47,   /* Terminal nitrogen in azide */
    MMFF94_TYPE_NGD_PLUS = 56,  /* Guanidinium nitrogen */
    MMFF94_TYPE_NPD_PLUS = 58,  /* Pyridinium nitrogen */
    MMFF94_TYPE_NR_PLUS  = 34,  /* Quaternary nitrogen, sp3 */
    MMFF94_TYPE_N_IM_PLUS = 54, /* Imidazolium nitrogen */
    MMFF94_TYPE_NM      = 62,   /* Anionic nitrogen (deprotonated amide) */
    MMFF94_TYPE_N2OX    = 45,   /* Nitrogen in nitro group sp2 */
    MMFF94_TYPE_NPOX    = 67,   /* Nitrogen in N-oxide */
    MMFF94_TYPE_N3OX    = 68,   /* Nitrogen in nitrate */
    MMFF94_TYPE_NPL3    = 39,   /* Trigonal nitrogen, sp2 (amide, enamine) */
    MMFF94_TYPE_NC_EQ_N = 40,   /* Nitrogen in N=C-N (amidine) */
    MMFF94_TYPE_NC_AR   = 81,   /* Nitrogen in aromatic C-N */

    /* Oxygen types (31-40) */
    MMFF94_TYPE_OR      = 6,    /* Ether/alcohol oxygen, sp3 */
    MMFF94_TYPE_O_EQ_C  = 7,    /* Carbonyl oxygen (C=O) */
    MMFF94_TYPE_O_AR    = 41,   /* Aromatic oxygen (furan) */
    MMFF94_TYPE_O_EQ_N  = 32,   /* Oxygen in N-oxide */
    MMFF94_TYPE_O_EQ_S  = 32,   /* Oxygen in sulfoxide/sulfone */
    MMFF94_TYPE_OM      = 35,   /* Carboxylate oxygen */
    MMFF94_TYPE_OM2     = 35,   /* Anionic oxygen */
    MMFF94_TYPE_O_PLUS  = 49,   /* Oxonium oxygen */
    MMFF94_TYPE_O3_EQ_N = 69,   /* Oxygen in nitrate */
    MMFF94_TYPE_OX_EQ_C = 7,    /* Carbonyl oxygen in ester/acid */

    /* Sulfur types (41-50) */
    MMFF94_TYPE_S       = 15,   /* Thiol/thioether sulfur, sp3 */
    MMFF94_TYPE_S_EQ_C  = 16,   /* Thione sulfur (C=S) */
    MMFF94_TYPE_S_AR    = 44,   /* Aromatic sulfur (thiophene) */
    MMFF94_TYPE_SO      = 17,   /* Sulfoxide sulfur */
    MMFF94_TYPE_SO2     = 18,   /* Sulfone sulfur */
    MMFF94_TYPE_SO2N    = 18,   /* Sulfonamide sulfur */
    MMFF94_TYPE_SO3     = 73,   /* Sulfonate sulfur */
    MMFF94_TYPE_SO4     = 73,   /* Sulfate sulfur */
    MMFF94_TYPE_SNO     = 18,   /* Sulfur in -S(=O)N< */
    MMFF94_TYPE_S_MINUS = 72,   /* Anionic sulfur */
    MMFF94_TYPE_STHI    = 72,   /* Thiolate sulfur */

    /* Phosphorus types (51-55) */
    MMFF94_TYPE_P       = 25,   /* Phosphine phosphorus, sp3 */
    MMFF94_TYPE_PO      = 26,   /* Phosphate phosphorus */
    MMFF94_TYPE_PO2     = 26,   /* Phosphonate phosphorus */
    MMFF94_TYPE_PO3     = 26,   /* Phosphate phosphorus */
    MMFF94_TYPE_PO4     = 26,   /* Phosphate phosphorus */

    /* Halogen types (56-60) */
    MMFF94_TYPE_F       = 11,   /* Fluorine */
    MMFF94_TYPE_CL      = 12,   /* Chlorine */
    MMFF94_TYPE_BR      = 13,   /* Bromine */
    MMFF94_TYPE_I       = 14,   /* Iodine */
    MMFF94_TYPE_F_MINUS = 89,   /* Fluoride ion */
    MMFF94_TYPE_CL_MINUS = 90,  /* Chloride ion */
    MMFF94_TYPE_BR_MINUS = 91,  /* Bromide ion */
    MMFF94_TYPE_I_MINUS  = 92,  /* Iodide ion */

    /* Hydrogen types (61-70) */
    MMFF94_TYPE_HC      = 5,    /* Hydrogen on carbon */
    MMFF94_TYPE_HO      = 21,   /* Hydrogen on oxygen (alcohol, acid) */
    MMFF94_TYPE_HN      = 23,   /* Hydrogen on nitrogen */
    MMFF94_TYPE_HOCO    = 24,   /* Hydrogen on carboxylic acid oxygen */
    MMFF94_TYPE_HOS     = 33,   /* Hydrogen on sulfur (thiol) */
    MMFF94_TYPE_HN_PLUS = 36,   /* Hydrogen on positive nitrogen */
    MMFF94_TYPE_HO_PLUS = 50,   /* Hydrogen on positive oxygen */
    MMFF94_TYPE_HN_EQ_C = 28,   /* Hydrogen on imine nitrogen */
    MMFF94_TYPE_HN_EQ_N = 28,   /* Hydrogen on N=N */
    MMFF94_TYPE_HNCO    = 28,   /* Hydrogen on amide nitrogen */
    MMFF94_TYPE_HNCS    = 28,   /* Hydrogen on thioamide nitrogen */
    MMFF94_TYPE_HNCC    = 29,   /* Hydrogen on enamine/aniline nitrogen */
    MMFF94_TYPE_HNCN    = 27,   /* Hydrogen on amidine/guanidine nitrogen */
    MMFF94_TYPE_HNN_PLUS = 36,  /* Hydrogen on positively charged nitrogen */
    MMFF94_TYPE_HP      = 71,   /* Hydrogen on phosphorus */
    MMFF94_TYPE_HSI     = 5,    /* Hydrogen on silicon */

    /* Silicon types (71-72) */
    MMFF94_TYPE_SI      = 19,   /* Silicon, sp3 */
    MMFF94_TYPE_SIOX    = 19,   /* Silicon in siloxane */

    /* Metals/other (73-83) */
    MMFF94_TYPE_LI_PLUS  = 92,  /* Lithium cation */
    MMFF94_TYPE_NA_PLUS  = 93,  /* Sodium cation */
    MMFF94_TYPE_K_PLUS   = 94,  /* Potassium cation */
    MMFF94_TYPE_ZN_2PLUS = 95,  /* Zinc dication */
    MMFF94_TYPE_CA_2PLUS = 96,  /* Calcium dication */
    MMFF94_TYPE_CU_PLUS  = 96,  /* Copper cation */
    MMFF94_TYPE_CU_2PLUS = 96,  /* Copper dication */
    MMFF94_TYPE_FE_2PLUS = 87,  /* Iron(II) */
    MMFF94_TYPE_FE_3PLUS = 88,  /* Iron(III) */
    MMFF94_TYPE_MG_2PLUS = 99,  /* Magnesium dication */

    MMFF94_TYPE_MAX = 100
} mmff94_atom_type_t;

/* ============================================================================
 * Hybridization States
 * ============================================================================ */

typedef enum {
    MMFF94_HYBRID_SP3 = 3,      /* Tetrahedral */
    MMFF94_HYBRID_SP2 = 2,      /* Trigonal planar */
    MMFF94_HYBRID_SP  = 1,      /* Linear */
    MMFF94_HYBRID_UNKNOWN = 0
} mmff94_hybridization_t;

/* ============================================================================
 * Per-Atom MMFF94 Data
 * ============================================================================ */

typedef struct {
    mmff94_atom_type_t type;            /* Assigned MMFF94 atom type */
    double partial_charge;              /* Computed partial charge */
    mmff94_hybridization_t hybridization; /* sp, sp2, sp3 */
    bool is_cation;                     /* Formal positive charge */
    bool is_anion;                      /* Formal negative charge */
    bool is_aromatic;                   /* In aromatic ring */
    int ring_size;                      /* Smallest ring size (0 if not in ring) */
    int num_pi_electrons;               /* Pi electron contribution */
} mmff94_atom_data_t;

/* ============================================================================
 * MMFF94 Context for Molecule
 * ============================================================================ */

typedef struct {
    mmff94_atom_data_t* atom_data;      /* Per-atom typing and charges */
    int num_atoms;                      /* Number of atoms */
    bool types_assigned;                /* Have atom types been assigned */
    bool charges_computed;              /* Have partial charges been computed */

    /* Cached interaction lists for performance */
    int* bond_type_i;                   /* Bond atom i types */
    int* bond_type_j;                   /* Bond atom j types */
    int num_bonds;

    /* Cached 1-4 pairs for VDW/electrostatic scaling */
    int* pairs_14;                      /* Flattened 1-4 pair list */
    int num_pairs_14;
} mmff94_context_t;

/* ============================================================================
 * Parameter Table Entry Structures
 * ============================================================================ */

/* Bond stretching parameters */
typedef struct {
    uint8_t type_i;                     /* Atom type i */
    uint8_t type_j;                     /* Atom type j */
    uint8_t bond_type;                  /* 0=single, 1=double, 2=triple, 3=aromatic */
    double kb;                          /* Force constant (mdyn/A) */
    double r0;                          /* Equilibrium length (A) */
} mmff94_bond_param_t;

/* Angle bending parameters */
typedef struct {
    uint8_t type_i;                     /* Terminal atom type i */
    uint8_t type_j;                     /* Central atom type j */
    uint8_t type_k;                     /* Terminal atom type k */
    uint8_t angle_type;                 /* Angle type (0-8) */
    double ka;                          /* Force constant (mdyn*A/rad^2) */
    double theta0;                      /* Equilibrium angle (degrees) */
} mmff94_angle_param_t;

/* Stretch-bend parameters */
typedef struct {
    uint8_t type_i;
    uint8_t type_j;                     /* Central atom */
    uint8_t type_k;
    uint8_t sb_type;                    /* Stretch-bend type */
    double kba_ijk;                     /* Force constant for i-j stretch */
    double kba_kji;                     /* Force constant for k-j stretch */
} mmff94_strbnd_param_t;

/* Torsion parameters (3-term Fourier) */
typedef struct {
    uint8_t type_i;
    uint8_t type_j;                     /* Central bond atom 1 */
    uint8_t type_k;                     /* Central bond atom 2 */
    uint8_t type_l;
    uint8_t torsion_type;               /* Torsion type (0-5) */
    double V1;                          /* cos(phi) coefficient */
    double V2;                          /* cos(2*phi) coefficient */
    double V3;                          /* cos(3*phi) coefficient */
} mmff94_torsion_param_t;

/* Out-of-plane bending parameters */
typedef struct {
    uint8_t type_i;                     /* Central atom */
    uint8_t type_j;                     /* Ligand 1 */
    uint8_t type_k;                     /* Ligand 2 */
    uint8_t type_l;                     /* Ligand 3 */
    double koop;                        /* Force constant */
} mmff94_oop_param_t;

/* Van der Waals parameters */
typedef struct {
    uint8_t type;                       /* Atom type */
    double alpha;                       /* Atomic polarizability (A^3) */
    double N_eff;                       /* Effective electrons for VDW */
    double A;                           /* Slater-Kirkwood A parameter */
    double G;                           /* Scale factor */
    double DA;                          /* Donor/acceptor designation */
} mmff94_vdw_param_t;

/* Bond charge increment parameters */
typedef struct {
    uint8_t type_i;
    uint8_t type_j;
    uint8_t bond_type;                  /* Bond type (0-3) */
    double bci;                         /* Partial charge increment i->j */
} mmff94_bci_param_t;

/* ============================================================================
 * Constants
 * ============================================================================ */

/* Energy conversion factors */
#define MMFF94_BOND_CONST       143.9325    /* kcal*A^2/(mol*mdyn) */
#define MMFF94_ANGLE_CONST      0.043844    /* kcal/(mol*rad^2) */
#define MMFF94_STRBND_CONST     2.51124     /* kcal/(mol*A*rad) */
#define MMFF94_TORSION_CONST    0.5         /* Energy scaling */
#define MMFF94_OOP_CONST        0.043844    /* Same as angle */
#define MMFF94_ELEC_CONST       332.0716    /* e^2/A to kcal/mol */

/* Force field parameters */
#define MMFF94_CUBIC_STRETCH    (-2.0)      /* Cubic stretch constant cs */
#define MMFF94_CUBIC_BEND       (-0.007)    /* Cubic bend constant cb (1/deg) */
#define MMFF94_ELEC_BUFFER      0.05        /* Electrostatic buffering (A) */
#define MMFF94_DIELECTRIC       1.0         /* Vacuum dielectric constant */

/* VDW buffered 14-7 constants */
#define MMFF94_VDW_GAMMA        0.12        /* R^7 buffer term */
#define MMFF94_VDW_DELTA        0.07        /* Distance buffer term */

/* 1-4 scaling factors */
#define MMFF94_VDW_14_SCALE     0.5         /* VDW 1-4 interaction scale */
#define MMFF94_ELEC_14_SCALE    0.75        /* Electrostatic 1-4 scale */

/* ============================================================================
 * Function Declarations
 * ============================================================================ */

/* Context management */
mmff94_context_t* mmff94_context_create(const molecule_t* mol);
void mmff94_context_free(mmff94_context_t* ctx);

/* Atom type assignment */
cchem_status_t mmff94_assign_types(const molecule_t* mol, mmff94_context_t* ctx);
mmff94_atom_type_t mmff94_get_atom_type(const molecule_t* mol, int atom_idx,
                                         const mmff94_context_t* ctx);

/* Partial charge computation */
cchem_status_t mmff94_compute_charges(const molecule_t* mol, mmff94_context_t* ctx);

/* Energy and gradient calculation */
double mmff94_calc_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                          const mol_coords_t* coords,
                          double* grad_x, double* grad_y, double* grad_z);

/* Individual energy terms (for debugging/testing) */
double mmff94_calc_bond_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                const mol_coords_t* coords,
                                double* gx, double* gy, double* gz);
double mmff94_calc_angle_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                 const mol_coords_t* coords,
                                 double* gx, double* gy, double* gz);
double mmff94_calc_strbnd_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                  const mol_coords_t* coords,
                                  double* gx, double* gy, double* gz);
double mmff94_calc_torsion_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                   const mol_coords_t* coords,
                                   double* gx, double* gy, double* gz);
double mmff94_calc_oop_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                               const mol_coords_t* coords,
                               double* gx, double* gy, double* gz);
double mmff94_calc_vdw_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                               const mol_coords_t* coords,
                               double* gx, double* gy, double* gz);
double mmff94_calc_elec_energy(const molecule_t* mol, const mmff94_context_t* ctx,
                                const mol_coords_t* coords,
                                double* gx, double* gy, double* gz);

/* Utility functions */
const char* mmff94_type_name(mmff94_atom_type_t type);
mmff94_hybridization_t mmff94_get_hybridization(const molecule_t* mol, int atom_idx);

#endif /* CCHEM_DEPICTOR_MMFF94_TYPES_H */
