/**
 * @file types.h
 * @brief Core type definitions for SMILES canonicalizer
 */

#ifndef CCHEM_CANONICALIZER_TYPES_H
#define CCHEM_CANONICALIZER_TYPES_H

/* Platform compatibility - must be first */
#include "cchem/compat.h"

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

/* Status codes */
typedef enum {
    CCHEM_OK = 0,
    CCHEM_ERROR_INVALID_INPUT = -1,
    CCHEM_ERROR_MEMORY = -2,
    CCHEM_ERROR_PARSE = -3,
    CCHEM_ERROR_INVALID_SMILES = -4,
    CCHEM_ERROR_RING_CLOSURE = -5,
    CCHEM_ERROR_STEREO = -6,
    CCHEM_ERROR_FILE_IO = -7,
    CCHEM_ERROR_THREAD = -8,
    CCHEM_ERROR_NOT_IMPLEMENTED = -9
} cchem_status_t;

/* Bond types */
typedef enum {
    BOND_NONE = 0,
    BOND_SINGLE = 1,
    BOND_DOUBLE = 2,
    BOND_TRIPLE = 3,
    BOND_AROMATIC = 4,
    BOND_UP = 5,      /* Stereochemistry: / */
    BOND_DOWN = 6,    /* Stereochemistry: \ */
    BOND_RING_SINGLE = 7,
    BOND_RING_AROMATIC = 8
} bond_type_t;

/* Atom chirality */
typedef enum {
    CHIRALITY_NONE = 0,
    CHIRALITY_CW = 1,     /* @ - clockwise */
    CHIRALITY_CCW = 2,    /* @@ - counter-clockwise */
    CHIRALITY_TH1 = 3,    /* @TH1 */
    CHIRALITY_TH2 = 4,    /* @TH2 */
    CHIRALITY_AL1 = 5,    /* @AL1 */
    CHIRALITY_AL2 = 6,    /* @AL2 */
    CHIRALITY_SP1 = 7,    /* @SP1 */
    CHIRALITY_SP2 = 8,    /* @SP2 */
    CHIRALITY_SP3 = 9,    /* @SP3 */
    CHIRALITY_OH1 = 10,   /* @OH1-30 */
    CHIRALITY_TB1 = 11    /* @TB1-20 */
} chirality_t;

/* E/Z stereochemistry for double bonds */
typedef enum {
    STEREO_NONE = 0,
    STEREO_CIS = 1,       /* Z */
    STEREO_TRANS = 2      /* E */
} stereo_ez_t;

/* Maximum limits */
#define MAX_ATOMS 10000
#define MAX_BONDS 15000
#define MAX_RINGS 1000
#define MAX_RING_SIZE 100
#define MAX_NEIGHBORS 10
#define MAX_SMILES_LENGTH 100000
#define MAX_ELEMENT_LENGTH 3

/* Element atomic numbers */
typedef enum {
    ELEM_UNKNOWN = 0,
    ELEM_H = 1,
    ELEM_He = 2,
    ELEM_Li = 3,
    ELEM_Be = 4,
    ELEM_B = 5,
    ELEM_C = 6,
    ELEM_N = 7,
    ELEM_O = 8,
    ELEM_F = 9,
    ELEM_Ne = 10,
    ELEM_Na = 11,
    ELEM_Mg = 12,
    ELEM_Al = 13,
    ELEM_Si = 14,
    ELEM_P = 15,
    ELEM_S = 16,
    ELEM_Cl = 17,
    ELEM_Ar = 18,
    ELEM_K = 19,
    ELEM_Ca = 20,
    ELEM_Sc = 21,
    ELEM_Ti = 22,
    ELEM_V = 23,
    ELEM_Cr = 24,
    ELEM_Mn = 25,
    ELEM_Fe = 26,
    ELEM_Co = 27,
    ELEM_Ni = 28,
    ELEM_Cu = 29,
    ELEM_Zn = 30,
    ELEM_Ga = 31,
    ELEM_Ge = 32,
    ELEM_As = 33,
    ELEM_Se = 34,
    ELEM_Br = 35,
    ELEM_Kr = 36,
    ELEM_Rb = 37,
    ELEM_Sr = 38,
    ELEM_Y = 39,
    ELEM_Zr = 40,
    ELEM_Nb = 41,
    ELEM_Mo = 42,
    ELEM_Tc = 43,
    ELEM_Ru = 44,
    ELEM_Rh = 45,
    ELEM_Pd = 46,
    ELEM_Ag = 47,
    ELEM_Cd = 48,
    ELEM_In = 49,
    ELEM_Sn = 50,
    ELEM_Sb = 51,
    ELEM_Te = 52,
    ELEM_I = 53,
    ELEM_Xe = 54,
    ELEM_Cs = 55,
    ELEM_Ba = 56,
    ELEM_La = 57,
    ELEM_Ce = 58,
    ELEM_Pr = 59,
    ELEM_Nd = 60,
    ELEM_Pm = 61,
    ELEM_Sm = 62,
    ELEM_Eu = 63,
    ELEM_Gd = 64,
    ELEM_Tb = 65,
    ELEM_Dy = 66,
    ELEM_Ho = 67,
    ELEM_Er = 68,
    ELEM_Tm = 69,
    ELEM_Yb = 70,
    ELEM_Lu = 71,
    ELEM_Hf = 72,
    ELEM_Ta = 73,
    ELEM_W = 74,
    ELEM_Re = 75,
    ELEM_Os = 76,
    ELEM_Ir = 77,
    ELEM_Pt = 78,
    ELEM_Au = 79,
    ELEM_Hg = 80,
    ELEM_Tl = 81,
    ELEM_Pb = 82,
    ELEM_Bi = 83,
    ELEM_Po = 84,
    ELEM_At = 85,
    ELEM_Rn = 86,
    ELEM_Fr = 87,
    ELEM_Ra = 88,
    ELEM_Ac = 89,
    ELEM_Th = 90,
    ELEM_Pa = 91,
    ELEM_U = 92,
    ELEM_Np = 93,
    ELEM_Pu = 94,
    ELEM_Am = 95,
    ELEM_Cm = 96,
    ELEM_Bk = 97,
    ELEM_Cf = 98,
    ELEM_Es = 99,
    ELEM_Fm = 100,
    ELEM_Md = 101,
    ELEM_No = 102,
    ELEM_Lr = 103,
    ELEM_Rf = 104,
    ELEM_Db = 105,
    ELEM_Sg = 106,
    ELEM_Bh = 107,
    ELEM_Hs = 108,
    ELEM_Mt = 109,
    ELEM_Ds = 110,
    ELEM_Rg = 111,
    ELEM_Cn = 112,
    ELEM_Nh = 113,
    ELEM_Fl = 114,
    ELEM_Mc = 115,
    ELEM_Lv = 116,
    ELEM_Ts = 117,
    ELEM_Og = 118
} element_t;

#endif /* CCHEM_CANONICALIZER_TYPES_H */
