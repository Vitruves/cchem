/**
 * @file colors.c
 * @brief CPK color scheme for molecular depiction
 */

#include "cchem/depictor/colors.h"

static const atom_display_info_t ATOM_DATA[] = {
    {ELEM_H,  {235, 235, 235}, 1.20, 0.31, "H"},
    {ELEM_He, {200, 230, 230}, 1.40, 0.28, "He"},
    {ELEM_Li, {180, 110, 220}, 1.82, 1.28, "Li"},
    {ELEM_Be, {170, 220, 50},  1.53, 0.96, "Be"},
    {ELEM_B,  {225, 160, 160}, 1.92, 0.84, "B"},
    {ELEM_C,  {45, 45, 45},    1.70, 0.76, "C"},
    {ELEM_N,  {42, 82, 190},   1.55, 0.71, "N"},
    {ELEM_O,  {220, 20, 60},   1.52, 0.66, "O"},
    {ELEM_F,  {120, 200, 70},  1.47, 0.57, "F"},
    {ELEM_Ne, {160, 210, 230}, 1.54, 0.58, "Ne"},
    {ELEM_Na, {150, 80, 210},  2.27, 1.66, "Na"},
    {ELEM_Mg, {120, 220, 50},  1.73, 1.41, "Mg"},
    {ELEM_Al, {175, 150, 150}, 1.84, 1.21, "Al"},
    {ELEM_Si, {210, 175, 140}, 2.10, 1.11, "Si"},
    {ELEM_P,  {204, 102, 0},   1.80, 1.07, "P"},
    {ELEM_S,  {204, 136, 0},   1.80, 1.05, "S"},
    {ELEM_Cl, {34, 139, 34},   1.75, 1.02, "Cl"},
    {ELEM_Ar, {110, 185, 200}, 1.88, 1.06, "Ar"},
    {ELEM_K,  {143, 64, 212},  2.75, 2.03, "K"},
    {ELEM_Ca, {61, 255, 0},    2.31, 1.76, "Ca"},
    {ELEM_Fe, {200, 90, 45},   1.26, 1.32, "Fe"},
    {ELEM_Cu, {184, 115, 51},  1.28, 1.32, "Cu"},
    {ELEM_Zn, {115, 118, 160}, 1.34, 1.22, "Zn"},
    {ELEM_Br, {139, 35, 35},   1.85, 1.20, "Br"},
    {ELEM_I,  {102, 51, 153},  1.98, 1.39, "I"},
    {ELEM_UNKNOWN, {180, 180, 180}, 1.70, 1.00, "?"}
};

#define NUM_ELEMENTS (sizeof(ATOM_DATA) / sizeof(ATOM_DATA[0]))

const atom_display_info_t* atom_get_display_info(element_t elem) {
    for (size_t i = 0; i < NUM_ELEMENTS - 1; i++) {
        if (ATOM_DATA[i].element == elem) {
            return &ATOM_DATA[i];
        }
    }
    return &ATOM_DATA[NUM_ELEMENTS - 1];
}

rgb_color_t atom_get_color(element_t elem) {
    return atom_get_display_info(elem)->color;
}

double atom_get_vdw_radius(element_t elem) {
    return atom_get_display_info(elem)->vdw_radius;
}

double atom_get_covalent_radius(element_t elem) {
    return atom_get_display_info(elem)->covalent_radius;
}
