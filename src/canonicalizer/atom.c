/**
 * @file atom.c
 * @brief Atom data structure and operations implementation
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/bond.h"

/* Element data table */
static const struct {
    element_t elem;
    const char* symbol;
    int default_valences[5];
    double mass;
    bool organic_subset;
} ELEMENT_DATA[] = {
    {ELEM_UNKNOWN, "*",  {0, -1}, 0.0, false},
    {ELEM_H,  "H",  {1, -1}, 1.008, false},
    {ELEM_He, "He", {0, -1}, 4.003, false},
    {ELEM_Li, "Li", {1, -1}, 6.941, false},
    {ELEM_Be, "Be", {2, -1}, 9.012, false},
    {ELEM_B,  "B",  {3, -1}, 10.81, true},
    {ELEM_C,  "C",  {4, -1}, 12.01, true},
    {ELEM_N,  "N",  {3, 5, -1}, 14.01, true},
    {ELEM_O,  "O",  {2, -1}, 16.00, true},
    {ELEM_F,  "F",  {1, -1}, 19.00, true},
    {ELEM_Ne, "Ne", {0, -1}, 20.18, false},
    {ELEM_Na, "Na", {1, -1}, 22.99, false},
    {ELEM_Mg, "Mg", {2, -1}, 24.31, false},
    {ELEM_Al, "Al", {3, -1}, 26.98, false},
    {ELEM_Si, "Si", {4, -1}, 28.09, false},
    {ELEM_P,  "P",  {3, 5, -1}, 30.97, true},
    {ELEM_S,  "S",  {2, 4, 6, -1}, 32.07, true},
    {ELEM_Cl, "Cl", {1, -1}, 35.45, true},
    {ELEM_Ar, "Ar", {0, -1}, 39.95, false},
    {ELEM_K,  "K",  {1, -1}, 39.10, false},
    {ELEM_Ca, "Ca", {2, -1}, 40.08, false},
    {ELEM_Sc, "Sc", {3, -1}, 44.96, false},
    {ELEM_Ti, "Ti", {4, -1}, 47.87, false},
    {ELEM_V,  "V",  {5, -1}, 50.94, false},
    {ELEM_Cr, "Cr", {3, 6, -1}, 52.00, false},
    {ELEM_Mn, "Mn", {2, 4, 7, -1}, 54.94, false},
    {ELEM_Fe, "Fe", {2, 3, -1}, 55.85, false},
    {ELEM_Co, "Co", {2, 3, -1}, 58.93, false},
    {ELEM_Ni, "Ni", {2, -1}, 58.69, false},
    {ELEM_Cu, "Cu", {1, 2, -1}, 63.55, false},
    {ELEM_Zn, "Zn", {2, -1}, 65.38, false},
    {ELEM_Ga, "Ga", {3, -1}, 69.72, false},
    {ELEM_Ge, "Ge", {4, -1}, 72.63, false},
    {ELEM_As, "As", {3, 5, -1}, 74.92, false},
    {ELEM_Se, "Se", {2, 4, 6, -1}, 78.97, false},
    {ELEM_Br, "Br", {1, -1}, 79.90, true},
    {ELEM_Kr, "Kr", {0, -1}, 83.80, false},
    {ELEM_Rb, "Rb", {1, -1}, 85.47, false},
    {ELEM_Sr, "Sr", {2, -1}, 87.62, false},
    {ELEM_Y,  "Y",  {3, -1}, 88.91, false},
    {ELEM_Zr, "Zr", {4, -1}, 91.22, false},
    {ELEM_Nb, "Nb", {5, -1}, 92.91, false},
    {ELEM_Mo, "Mo", {6, -1}, 95.95, false},
    {ELEM_Tc, "Tc", {7, -1}, 98.00, false},
    {ELEM_Ru, "Ru", {4, -1}, 101.1, false},
    {ELEM_Rh, "Rh", {3, -1}, 102.9, false},
    {ELEM_Pd, "Pd", {2, 4, -1}, 106.4, false},
    {ELEM_Ag, "Ag", {1, -1}, 107.9, false},
    {ELEM_Cd, "Cd", {2, -1}, 112.4, false},
    {ELEM_In, "In", {3, -1}, 114.8, false},
    {ELEM_Sn, "Sn", {2, 4, -1}, 118.7, false},
    {ELEM_Sb, "Sb", {3, 5, -1}, 121.8, false},
    {ELEM_Te, "Te", {2, 4, 6, -1}, 127.6, false},
    {ELEM_I,  "I",  {1, -1}, 126.9, true},
    {ELEM_Xe, "Xe", {0, 2, -1}, 131.3, false},
    {ELEM_Cs, "Cs", {1, -1}, 132.9, false},
    {ELEM_Ba, "Ba", {2, -1}, 137.3, false},
    {ELEM_La, "La", {3, -1}, 138.9, false},
    {ELEM_Ce, "Ce", {3, 4, -1}, 140.1, false},
    {ELEM_Pr, "Pr", {3, -1}, 140.9, false},
    {ELEM_Nd, "Nd", {3, -1}, 144.2, false},
    {ELEM_Pm, "Pm", {3, -1}, 145.0, false},
    {ELEM_Sm, "Sm", {3, -1}, 150.4, false},
    {ELEM_Eu, "Eu", {3, -1}, 152.0, false},
    {ELEM_Gd, "Gd", {3, -1}, 157.3, false},
    {ELEM_Tb, "Tb", {3, -1}, 158.9, false},
    {ELEM_Dy, "Dy", {3, -1}, 162.5, false},
    {ELEM_Ho, "Ho", {3, -1}, 164.9, false},
    {ELEM_Er, "Er", {3, -1}, 167.3, false},
    {ELEM_Tm, "Tm", {3, -1}, 168.9, false},
    {ELEM_Yb, "Yb", {3, -1}, 173.0, false},
    {ELEM_Lu, "Lu", {3, -1}, 175.0, false},
    {ELEM_Hf, "Hf", {4, -1}, 178.5, false},
    {ELEM_Ta, "Ta", {5, -1}, 180.9, false},
    {ELEM_W,  "W",  {6, -1}, 183.8, false},
    {ELEM_Re, "Re", {7, -1}, 186.2, false},
    {ELEM_Os, "Os", {4, -1}, 190.2, false},
    {ELEM_Ir, "Ir", {4, -1}, 192.2, false},
    {ELEM_Pt, "Pt", {2, 4, -1}, 195.1, false},
    {ELEM_Au, "Au", {1, 3, -1}, 197.0, false},
    {ELEM_Hg, "Hg", {1, 2, -1}, 200.6, false},
    {ELEM_Tl, "Tl", {1, 3, -1}, 204.4, false},
    {ELEM_Pb, "Pb", {2, 4, -1}, 207.2, false},
    {ELEM_Bi, "Bi", {3, 5, -1}, 209.0, false},
    {ELEM_Po, "Po", {2, 4, -1}, 209.0, false},
    {ELEM_At, "At", {1, -1}, 210.0, false},
    {ELEM_Rn, "Rn", {0, -1}, 222.0, false},
    {ELEM_Fr, "Fr", {1, -1}, 223.0, false},
    {ELEM_Ra, "Ra", {2, -1}, 226.0, false},
    {ELEM_Ac, "Ac", {3, -1}, 227.0, false},
    {ELEM_Th, "Th", {4, -1}, 232.0, false},
    {ELEM_Pa, "Pa", {5, -1}, 231.0, false},
    {ELEM_U,  "U",  {6, -1}, 238.0, false},
    {ELEM_Np, "Np", {5, -1}, 237.0, false},
    {ELEM_Pu, "Pu", {4, -1}, 244.0, false},
    {ELEM_Am, "Am", {3, -1}, 243.0, false},
    {ELEM_Cm, "Cm", {3, -1}, 247.0, false},
    {ELEM_Bk, "Bk", {3, -1}, 247.0, false},
    {ELEM_Cf, "Cf", {3, -1}, 251.0, false},
    {ELEM_Es, "Es", {3, -1}, 252.0, false},
    {ELEM_Fm, "Fm", {3, -1}, 257.0, false},
    {ELEM_Md, "Md", {3, -1}, 258.0, false},
    {ELEM_No, "No", {3, -1}, 259.0, false},
    {ELEM_Lr, "Lr", {3, -1}, 262.0, false},
    {ELEM_Rf, "Rf", {4, -1}, 267.0, false},
    {ELEM_Db, "Db", {5, -1}, 270.0, false},
    {ELEM_Sg, "Sg", {6, -1}, 271.0, false},
    {ELEM_Bh, "Bh", {7, -1}, 270.0, false},
    {ELEM_Hs, "Hs", {8, -1}, 277.0, false},
    {ELEM_Mt, "Mt", {0, -1}, 276.0, false},
    {ELEM_Ds, "Ds", {0, -1}, 281.0, false},
    {ELEM_Rg, "Rg", {0, -1}, 282.0, false},
    {ELEM_Cn, "Cn", {0, -1}, 285.0, false},
    {ELEM_Nh, "Nh", {0, -1}, 286.0, false},
    {ELEM_Fl, "Fl", {0, -1}, 289.0, false},
    {ELEM_Mc, "Mc", {0, -1}, 290.0, false},
    {ELEM_Lv, "Lv", {0, -1}, 293.0, false},
    {ELEM_Ts, "Ts", {0, -1}, 294.0, false},
    {ELEM_Og, "Og", {0, -1}, 294.0, false}
};

static const int NUM_ELEMENTS = sizeof(ELEMENT_DATA) / sizeof(ELEMENT_DATA[0]);

/* Element functions */
const char* element_to_symbol(element_t elem) {
    if (elem >= 0 && elem < NUM_ELEMENTS) {
        return ELEMENT_DATA[elem].symbol;
    }
    return "*";
}

/* Fast element lookup using first character as hash */
element_t element_from_symbol(const char* symbol) {
    if (!symbol || !symbol[0]) return ELEM_UNKNOWN;

    char c1 = symbol[0];
    char c2 = symbol[1];

    /* Fast path for common single-letter elements */
    if (c2 == '\0') {
        switch (c1) {
            case 'H': return ELEM_H;
            case 'B': return ELEM_B;
            case 'C': return ELEM_C;
            case 'N': return ELEM_N;
            case 'O': return ELEM_O;
            case 'F': return ELEM_F;
            case 'P': return ELEM_P;
            case 'S': return ELEM_S;
            case 'K': return ELEM_K;
            case 'V': return ELEM_V;
            case 'Y': return ELEM_Y;
            case 'I': return ELEM_I;
            case 'W': return ELEM_W;
            case 'U': return ELEM_U;
            case '*': return ELEM_UNKNOWN;
        }
    }
    /* Fast path for common two-letter elements */
    else if (symbol[2] == '\0') {
        switch (c1) {
            case 'A': if (c2 == 'l') return ELEM_Al;
                      if (c2 == 'r') return ELEM_Ar;
                      if (c2 == 's') return ELEM_As;
                      if (c2 == 'g') return ELEM_Ag;
                      if (c2 == 'u') return ELEM_Au;
                      if (c2 == 't') return ELEM_At;
                      if (c2 == 'c') return ELEM_Ac;
                      if (c2 == 'm') return ELEM_Am;
                      break;
            case 'B': if (c2 == 'e') return ELEM_Be;
                      if (c2 == 'r') return ELEM_Br;
                      if (c2 == 'a') return ELEM_Ba;
                      if (c2 == 'i') return ELEM_Bi;
                      if (c2 == 'k') return ELEM_Bk;
                      if (c2 == 'h') return ELEM_Bh;
                      break;
            case 'C': if (c2 == 'l') return ELEM_Cl;
                      if (c2 == 'a') return ELEM_Ca;
                      if (c2 == 'r') return ELEM_Cr;
                      if (c2 == 'o') return ELEM_Co;
                      if (c2 == 'u') return ELEM_Cu;
                      if (c2 == 'd') return ELEM_Cd;
                      if (c2 == 's') return ELEM_Cs;
                      if (c2 == 'e') return ELEM_Ce;
                      if (c2 == 'f') return ELEM_Cf;
                      if (c2 == 'm') return ELEM_Cm;
                      if (c2 == 'n') return ELEM_Cn;
                      break;
            case 'D': if (c2 == 'y') return ELEM_Dy;
                      if (c2 == 'b') return ELEM_Db;
                      if (c2 == 's') return ELEM_Ds;
                      break;
            case 'E': if (c2 == 'u') return ELEM_Eu;
                      if (c2 == 'r') return ELEM_Er;
                      if (c2 == 's') return ELEM_Es;
                      break;
            case 'F': if (c2 == 'e') return ELEM_Fe;
                      if (c2 == 'r') return ELEM_Fr;
                      if (c2 == 'l') return ELEM_Fl;
                      if (c2 == 'm') return ELEM_Fm;
                      break;
            case 'G': if (c2 == 'a') return ELEM_Ga;
                      if (c2 == 'e') return ELEM_Ge;
                      if (c2 == 'd') return ELEM_Gd;
                      break;
            case 'H': if (c2 == 'e') return ELEM_He;
                      if (c2 == 'g') return ELEM_Hg;
                      if (c2 == 'f') return ELEM_Hf;
                      if (c2 == 'o') return ELEM_Ho;
                      if (c2 == 's') return ELEM_Hs;
                      break;
            case 'I': if (c2 == 'n') return ELEM_In;
                      if (c2 == 'r') return ELEM_Ir;
                      break;
            case 'K': if (c2 == 'r') return ELEM_Kr;
                      break;
            case 'L': if (c2 == 'i') return ELEM_Li;
                      if (c2 == 'a') return ELEM_La;
                      if (c2 == 'u') return ELEM_Lu;
                      if (c2 == 'r') return ELEM_Lr;
                      if (c2 == 'v') return ELEM_Lv;
                      break;
            case 'M': if (c2 == 'g') return ELEM_Mg;
                      if (c2 == 'n') return ELEM_Mn;
                      if (c2 == 'o') return ELEM_Mo;
                      if (c2 == 't') return ELEM_Mt;
                      if (c2 == 'd') return ELEM_Md;
                      if (c2 == 'c') return ELEM_Mc;
                      break;
            case 'N': if (c2 == 'e') return ELEM_Ne;
                      if (c2 == 'a') return ELEM_Na;
                      if (c2 == 'i') return ELEM_Ni;
                      if (c2 == 'b') return ELEM_Nb;
                      if (c2 == 'd') return ELEM_Nd;
                      if (c2 == 'p') return ELEM_Np;
                      if (c2 == 'o') return ELEM_No;
                      if (c2 == 'h') return ELEM_Nh;
                      break;
            case 'O': if (c2 == 's') return ELEM_Os;
                      if (c2 == 'g') return ELEM_Og;
                      break;
            case 'P': if (c2 == 'd') return ELEM_Pd;
                      if (c2 == 't') return ELEM_Pt;
                      if (c2 == 'b') return ELEM_Pb;
                      if (c2 == 'o') return ELEM_Po;
                      if (c2 == 'u') return ELEM_Pu;
                      if (c2 == 'r') return ELEM_Pr;
                      if (c2 == 'm') return ELEM_Pm;
                      if (c2 == 'a') return ELEM_Pa;
                      break;
            case 'R': if (c2 == 'b') return ELEM_Rb;
                      if (c2 == 'u') return ELEM_Ru;
                      if (c2 == 'h') return ELEM_Rh;
                      if (c2 == 'e') return ELEM_Re;
                      if (c2 == 'a') return ELEM_Ra;
                      if (c2 == 'n') return ELEM_Rn;
                      if (c2 == 'f') return ELEM_Rf;
                      if (c2 == 'g') return ELEM_Rg;
                      break;
            case 'S': if (c2 == 'i') return ELEM_Si;
                      if (c2 == 'e') return ELEM_Se;
                      if (c2 == 'r') return ELEM_Sr;
                      if (c2 == 'n') return ELEM_Sn;
                      if (c2 == 'b') return ELEM_Sb;
                      if (c2 == 'c') return ELEM_Sc;
                      if (c2 == 'm') return ELEM_Sm;
                      if (c2 == 'g') return ELEM_Sg;
                      break;
            case 'T': if (c2 == 'i') return ELEM_Ti;
                      if (c2 == 'e') return ELEM_Te;
                      if (c2 == 'c') return ELEM_Tc;
                      if (c2 == 'l') return ELEM_Tl;
                      if (c2 == 'a') return ELEM_Ta;
                      if (c2 == 'h') return ELEM_Th;
                      if (c2 == 'b') return ELEM_Tb;
                      if (c2 == 'm') return ELEM_Tm;
                      if (c2 == 's') return ELEM_Ts;
                      break;
            case 'X': if (c2 == 'e') return ELEM_Xe;
                      break;
            case 'Y': if (c2 == 'b') return ELEM_Yb;
                      break;
            case 'Z': if (c2 == 'n') return ELEM_Zn;
                      if (c2 == 'r') return ELEM_Zr;
                      break;
        }
    }

    /* Fallback: linear search for rare elements */
    for (int i = 0; i < NUM_ELEMENTS; i++) {
        if (strcmp(ELEMENT_DATA[i].symbol, symbol) == 0) {
            return ELEMENT_DATA[i].elem;
        }
    }
    return ELEM_UNKNOWN;
}

bool element_is_organic_subset(element_t elem) {
    if (elem >= 0 && elem < NUM_ELEMENTS) {
        return ELEMENT_DATA[elem].organic_subset;
    }
    return false;
}

int element_default_valence(element_t elem, int charge) {
    if (elem < 0 || elem >= NUM_ELEMENTS) return 0;

    int base_valence = ELEMENT_DATA[elem].default_valences[0];
    /* Positive charge increases valence (e.g., N+ can form 4 bonds)
     * Negative charge decreases valence (e.g., O- forms only 1 bond) */
    int adjusted = base_valence + charge;
    return (adjusted > 0) ? adjusted : 0;
}

double element_atomic_mass(element_t elem) {
    if (elem >= 0 && elem < NUM_ELEMENTS) {
        return ELEMENT_DATA[elem].mass;
    }
    return 0.0;
}

/* Atom functions */
atom_t* atom_create(element_t element) {
    atom_t* atom = (atom_t*)calloc(1, sizeof(atom_t));
    if (!atom) return NULL;

    atom_init(atom, element);
    return atom;
}

void atom_free(atom_t* atom) {
    if (atom) {
        free(atom);
    }
}

void atom_init(atom_t* atom, element_t element) {
    if (!atom) return;

    memset(atom, 0, sizeof(atom_t));
    atom->element = element;
    atom->isotope = 0;
    atom->charge = 0;
    atom->h_count = -1;  /* Unspecified */
    atom->aromatic = false;
    atom->chirality = CHIRALITY_NONE;
    atom->chirality_class = 0;
    atom->index = -1;
    atom->num_neighbors = 0;
    atom->ring_count = 0;
    atom->implicit_h_count = 0;
    atom->total_bond_order = 0;
    atom->invariant = 0;
    atom->canon_rank = -1;
    atom->visited = false;
    atom->atom_class = 0;
    atom->input_order = -1;
    atom->num_stereo_neighbors = 0;
}

void atom_reset(atom_t* atom) {
    if (atom) {
        element_t elem = atom->element;
        atom_init(atom, elem);
    }
}

cchem_status_t atom_add_neighbor(atom_t* atom, int neighbor_idx, int bond_idx) {
    if (!atom) return CCHEM_ERROR_INVALID_INPUT;
    if (atom->num_neighbors >= MAX_NEIGHBORS) return CCHEM_ERROR_INVALID_INPUT;

    /* Check for duplicate */
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (atom->neighbors[i] == neighbor_idx) {
            return CCHEM_OK;  /* Already connected */
        }
    }

    atom->neighbors[atom->num_neighbors] = neighbor_idx;
    atom->neighbor_bonds[atom->num_neighbors] = bond_idx;
    atom->num_neighbors++;

    return CCHEM_OK;
}

cchem_status_t atom_remove_neighbor(atom_t* atom, int neighbor_idx) {
    if (!atom) return CCHEM_ERROR_INVALID_INPUT;

    for (int i = 0; i < atom->num_neighbors; i++) {
        if (atom->neighbors[i] == neighbor_idx) {
            /* Shift remaining neighbors */
            for (int j = i; j < atom->num_neighbors - 1; j++) {
                atom->neighbors[j] = atom->neighbors[j + 1];
                atom->neighbor_bonds[j] = atom->neighbor_bonds[j + 1];
            }
            atom->num_neighbors--;
            return CCHEM_OK;
        }
    }

    return CCHEM_ERROR_INVALID_INPUT;  /* Not found */
}

int atom_get_bond_to(const atom_t* atom, int neighbor_idx) {
    if (!atom) return -1;

    for (int i = 0; i < atom->num_neighbors; i++) {
        if (atom->neighbors[i] == neighbor_idx) {
            return atom->neighbor_bonds[i];
        }
    }

    return -1;
}

bool atom_is_bonded_to(const atom_t* atom, int other_idx) {
    if (!atom) return false;

    for (int i = 0; i < atom->num_neighbors; i++) {
        if (atom->neighbors[i] == other_idx) {
            return true;
        }
    }

    return false;
}

int atom_calc_implicit_h(const atom_t* atom, const struct molecule* mol) {
    if (!atom || !mol) return 0;

    /* If explicit h_count specified, use it */
    if (atom->h_count >= 0) {
        return atom->h_count;
    }

    /* Calculate from valence */
    int valence = element_default_valence(atom->element, atom->charge);

    /* Sum of bond orders */
    int bond_sum = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        int bond_idx = atom->neighbor_bonds[i];
        if (bond_idx >= 0 && bond_idx < mol->num_bonds) {
            bond_sum += bond_get_int_order(&mol->bonds[bond_idx]);
        }
    }

    /* Aromatic atoms contribute 1 to bond sum for ring bond */
    if (atom->aromatic) {
        bond_sum += 1;  /* Aromatic contribution */
    }

    int implicit_h = valence - bond_sum;
    return (implicit_h > 0) ? implicit_h : 0;
}

int atom_get_degree(const atom_t* atom) {
    if (!atom) return 0;
    return atom->num_neighbors + atom->implicit_h_count;
}

int atom_get_heavy_degree(const atom_t* atom) {
    if (!atom) return 0;
    return atom->num_neighbors;
}

bool atom_in_ring(const atom_t* atom) {
    if (!atom) return false;
    return atom->ring_count > 0;
}

void atom_copy(atom_t* dest, const atom_t* src) {
    if (!dest || !src) return;
    memcpy(dest, src, sizeof(atom_t));
}
