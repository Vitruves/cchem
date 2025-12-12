/**
 * @file element.h
 * @brief Element data and utilities
 */

#ifndef CCHEM_CANONICALIZER_ELEMENT_H
#define CCHEM_CANONICALIZER_ELEMENT_H

#include "types.h"

/* Element data structure */
typedef struct {
    element_t atomic_num;
    const char* symbol;
    const char* name;
    int default_valence[4];  /* Common valences, -1 terminated */
    int num_valences;
    double atomic_mass;
    int period;
    int group;
    bool is_metal;
} element_data_t;

/* Get element data by atomic number */
const element_data_t* element_get_data(element_t elem);

/* Get element from symbol string */
element_t element_from_symbol(const char* symbol);

/* Get element symbol string */
const char* element_to_symbol(element_t elem);

/* Check if element is in organic subset (can be written without brackets) */
bool element_is_organic_subset(element_t elem);

/* Get default valence for element */
int element_default_valence(element_t elem, int charge);

/* Check if valence is valid for element */
bool element_valence_is_valid(element_t elem, int valence, int charge);

/* Get atomic mass */
double element_atomic_mass(element_t elem);

#endif /* CCHEM_CANONICALIZER_ELEMENT_H */
