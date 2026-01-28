/**
 * @file cchem_api.c
 * @brief Library API functions for cchem
 *
 * These functions are the main entry points for the cchem library.
 * Separated from main.c for use in Python bindings and other embeddings.
 */

#include "cchem/cchem.h"

/* Version string is defined in cchem/cchem.h */

const char* cchem_version(void) {
    return CCHEM_VERSION_STRING;
}

cchem_status_t cchem_init(void) {
    return CCHEM_OK;
}

void cchem_cleanup(void) {
    /* Currently no cleanup needed */
}
