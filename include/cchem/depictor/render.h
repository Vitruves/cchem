/**
 * @file render.h
 * @brief Molecular rendering with Cairo
 */

#ifndef CCHEM_DEPICTOR_RENDER_H
#define CCHEM_DEPICTOR_RENDER_H

#include "cchem/depictor/types.h"
#include "cchem/canonicalizer/molecule.h"

typedef struct render_context render_context_t;

render_context_t* render_context_create(int width, int height, rgb_color_t background);
void render_context_free(render_context_t* ctx);
cchem_status_t render_molecule(render_context_t* ctx, const molecule_t* mol,
                               const mol_coords_t* coords, const depictor_options_t* opts);
cchem_status_t render_save_png(render_context_t* ctx, const char* filename);
cchem_status_t render_save_jpeg(render_context_t* ctx, const char* filename, int quality);

#endif /* CCHEM_DEPICTOR_RENDER_H */
