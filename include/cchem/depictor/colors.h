/**
 * @file colors.h
 * @brief CPK color scheme for molecular depiction
 */

#ifndef CCHEM_DEPICTOR_COLORS_H
#define CCHEM_DEPICTOR_COLORS_H

#include "cchem/depictor/types.h"
#include "cchem/canonicalizer/types.h"

typedef struct {
    element_t element;
    rgb_color_t color;
    double vdw_radius;
    double covalent_radius;
    const char* name;
} atom_display_info_t;

const atom_display_info_t* atom_get_display_info(element_t elem);
rgb_color_t atom_get_color(element_t elem);
double atom_get_vdw_radius(element_t elem);
double atom_get_covalent_radius(element_t elem);

#endif /* CCHEM_DEPICTOR_COLORS_H */
