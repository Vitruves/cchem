/**
 * @file render.h
 * @brief Molecular rendering with Cairo
 */

#ifndef CCHEM_DEPICTOR_RENDER_H
#define CCHEM_DEPICTOR_RENDER_H

#include "cchem/depictor/types.h"
#include "cchem/canonicalizer/molecule.h"

typedef struct render_context render_context_t;

/* Create render context - for SVG output, pass IMG_FORMAT_SVG to render_context_create_ex */
render_context_t* render_context_create(int width, int height, rgb_color_t background);
render_context_t* render_context_create_ex(int width, int height, rgb_color_t background,
                                            image_format_t format, const char* svg_filename);
void render_context_free(render_context_t* ctx);

/* Render molecule to context */
cchem_status_t render_molecule(render_context_t* ctx, const molecule_t* mol,
                               const mol_coords_t* coords, const depictor_options_t* opts);

/* Save output */
cchem_status_t render_save_png(render_context_t* ctx, const char* filename);
cchem_status_t render_save_jpeg(render_context_t* ctx, const char* filename, int quality);
cchem_status_t render_save_svg(render_context_t* ctx, const char* filename);

/* Helper: calculate gap endpoints for bonds terminating at heteroatoms */
void calculate_bond_gap(point2d_t p1, point2d_t p2,
                        bool gap_at_p1, bool gap_at_p2,
                        double gap_factor, double font_size,
                        point2d_t* out_p1, point2d_t* out_p2);

#endif /* CCHEM_DEPICTOR_RENDER_H */
