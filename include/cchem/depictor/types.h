/**
 * @file types.h
 * @brief Core types for molecular depiction
 */

#ifndef CCHEM_DEPICTOR_TYPES_H
#define CCHEM_DEPICTOR_TYPES_H

#include <stdint.h>
#include <stdbool.h>
#include "cchem/canonicalizer/types.h"

/* ============== Geometric Types ============== */

typedef struct {
    double x;
    double y;
} point2d_t;

typedef struct {
    double x;
    double y;
    double z;
} point3d_t;

/* RGB color (0-255) */
typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
} rgb_color_t;

/* ============== Molecular Coordinates ============== */

typedef struct {
    point2d_t* coords_2d;
    point3d_t* coords_3d;
    int num_atoms;
    bool has_2d;
    bool has_3d;
} mol_coords_t;

/* ============== Depiction Options ============== */

typedef enum {
    IMG_FORMAT_PNG = 0,
    IMG_FORMAT_JPEG = 1,
    IMG_FORMAT_SVG = 2
} image_format_t;

typedef enum {
    LINE_CAP_ROUND = 0,
    LINE_CAP_BUTT = 1,
    LINE_CAP_SQUARE = 2
} line_cap_t;

typedef enum {
    DEPICT_STYLE_DEFAULT = 0,     /* Default cchem style with atom circles */
    DEPICT_STYLE_MODERN = 1       /* Modern publication style: no circles, bond gaps at heteroatoms */
} depict_style_preset_t;

typedef enum {
    DEPICT_MODE_2D = 0,
    DEPICT_MODE_3D = 1
} depict_mode_t;

typedef enum {
    RENDER_STYLE_WIREFRAME = 0,      /* Lines only, no atom spheres */
    RENDER_STYLE_STICKS = 1,         /* Thick bonds, small atom caps */
    RENDER_STYLE_BALLS_AND_STICKS = 2, /* Spheres (proportional) + bonds */
    RENDER_STYLE_SPACEFILL = 3,      /* VDW spheres, no explicit bonds */
    RENDER_STYLE_SURFACE = 4         /* Molecular surface representation */
} render_style_t;

typedef enum {
    SURFACE_COLOR_UNIFORM = 0,       /* Single uniform color (blue) */
    SURFACE_COLOR_ATOM = 1,          /* Color by nearest atom element */
    SURFACE_COLOR_POLARITY = 2       /* Color by partial charge (red=negative, blue=positive) */
} surface_color_mode_t;

typedef struct {
    depict_mode_t mode;
    render_style_t render_style;
    surface_color_mode_t surface_color;  /* Coloring mode for surface representation */
    depict_style_preset_t style_preset;  /* Overall style preset (default or chemdraw) */
    int width;
    int height;
    int margin;
    double bond_length;
    double atom_radius_scale;
    double bond_width;
    rgb_color_t background;
    bool show_carbons;
    bool show_hydrogens;
    bool draw_aromatic_circles;
    bool atom_filling;
    bool proportional_atoms;    /* Scale atoms by VDW radius (for 3D styles) */
    double font_size;
    int max_iterations;
    image_format_t format;
    int jpeg_quality;
    /* Publication-quality rendering options */
    double heteroatom_gap;      /* Gap factor for bonds at heteroatoms (0.0 = no gap, 1.0 = full) */
    line_cap_t line_cap;        /* Line cap style for bonds */
    double scale_factor;        /* Resolution scale factor (1.0 = normal, 2.0 = 2x resolution) */
    double double_bond_offset;  /* Offset for double bonds from center line */
    double triple_bond_offset;  /* Offset for triple bonds from center line */
    bool terminal_carbon_labels; /* Show labels on terminal carbons (CH3) */
} depictor_options_t;

/* Default: Publication-quality rendering */
#define DEPICTOR_OPTIONS_DEFAULT { \
    .mode = DEPICT_MODE_2D, \
    .render_style = RENDER_STYLE_WIREFRAME, \
    .surface_color = SURFACE_COLOR_UNIFORM, \
    .style_preset = DEPICT_STYLE_MODERN, \
    .width = 800, \
    .height = 800, \
    .margin = 50, \
    .bond_length = 35.0, \
    .atom_radius_scale = 0.60, \
    .bond_width = 2.5, \
    .background = {255, 255, 255}, \
    .show_carbons = false, \
    .show_hydrogens = false, \
    .draw_aromatic_circles = false, \
    .atom_filling = false, \
    .proportional_atoms = true, \
    .font_size = 3.5, \
    .max_iterations = 500, \
    .format = IMG_FORMAT_PNG, \
    .jpeg_quality = 95, \
    .heteroatom_gap = 1.0, \
    .line_cap = LINE_CAP_BUTT, \
    .scale_factor = 1.0, \
    .double_bond_offset = 0.12, \
    .triple_bond_offset = 0.18, \
    .terminal_carbon_labels = false \
}

/* ============== Coordinate Management ============== */

mol_coords_t* mol_coords_create(int num_atoms);
void mol_coords_free(mol_coords_t* coords);
mol_coords_t* mol_coords_copy(const mol_coords_t* coords);

/* ============== 2D Point Operations ============== */

point2d_t point2d_add(point2d_t a, point2d_t b);
point2d_t point2d_sub(point2d_t a, point2d_t b);
point2d_t point2d_scale(point2d_t p, double s);
double point2d_length(point2d_t p);
double point2d_distance(point2d_t a, point2d_t b);
point2d_t point2d_normalize(point2d_t p);
double point2d_dot(point2d_t a, point2d_t b);
double point2d_cross(point2d_t a, point2d_t b);
point2d_t point2d_rotate(point2d_t p, double angle);
point2d_t point2d_perp(point2d_t p);

/* ============== 3D Point Operations ============== */

point3d_t point3d_add(point3d_t a, point3d_t b);
point3d_t point3d_sub(point3d_t a, point3d_t b);
point3d_t point3d_scale(point3d_t p, double s);
double point3d_length(point3d_t p);
double point3d_distance(point3d_t a, point3d_t b);
point3d_t point3d_normalize(point3d_t p);
double point3d_dot(point3d_t a, point3d_t b);
point3d_t point3d_cross(point3d_t a, point3d_t b);

#endif /* CCHEM_DEPICTOR_TYPES_H */
