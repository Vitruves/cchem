/**
 * @file render.c
 * @brief Molecular rendering with Cairo
 */

#include "cchem/depictor/render.h"
#include "cchem/depictor/colors.h"
#include "cchem/canonicalizer/element.h"
#include "cchem/canonicalizer/bond.h"
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Base font size constant (internal) */
#define BASE_FONT_SIZE 3.5

struct render_context {
    cairo_surface_t* surface;
    cairo_t* cr;
    int width;
    int height;
    rgb_color_t background;
    image_format_t format;
    char* svg_filename;
};

/* ============== Context Management ============== */

render_context_t* render_context_create(int width, int height, rgb_color_t background) {
    return render_context_create_ex(width, height, background, IMG_FORMAT_PNG, NULL);
}

render_context_t* render_context_create_ex(int width, int height, rgb_color_t background,
                                            image_format_t format, const char* svg_filename) {
    render_context_t* ctx = calloc(1, sizeof(render_context_t));
    if (!ctx) return NULL;

    ctx->format = format;
    ctx->width = width;
    ctx->height = height;
    ctx->background = background;
    ctx->svg_filename = svg_filename ? strdup(svg_filename) : NULL;

    if (format == IMG_FORMAT_SVG && svg_filename) {
        ctx->surface = cairo_svg_surface_create(svg_filename, width, height);
    } else {
        ctx->surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    }

    if (cairo_surface_status(ctx->surface) != CAIRO_STATUS_SUCCESS) {
        free(ctx->svg_filename);
        free(ctx);
        return NULL;
    }

    ctx->cr = cairo_create(ctx->surface);
    if (cairo_status(ctx->cr) != CAIRO_STATUS_SUCCESS) {
        cairo_surface_destroy(ctx->surface);
        free(ctx->svg_filename);
        free(ctx);
        return NULL;
    }

    /* Enable best quality anti-aliasing */
    cairo_set_antialias(ctx->cr, CAIRO_ANTIALIAS_BEST);

    /* Fill background */
    cairo_set_source_rgb(ctx->cr,
                        background.r / 255.0,
                        background.g / 255.0,
                        background.b / 255.0);
    cairo_paint(ctx->cr);

    return ctx;
}

void render_context_free(render_context_t* ctx) {
    if (!ctx) return;
    if (ctx->cr) cairo_destroy(ctx->cr);
    if (ctx->surface) cairo_surface_destroy(ctx->surface);
    free(ctx->svg_filename);
    free(ctx);
}

/* ============== Drawing Primitives ============== */

static void set_color(cairo_t* cr, rgb_color_t color) {
    cairo_set_source_rgb(cr, color.r / 255.0, color.g / 255.0, color.b / 255.0);
}

static cairo_line_cap_t get_cairo_line_cap(line_cap_t cap) {
    switch (cap) {
        case LINE_CAP_BUTT:   return CAIRO_LINE_CAP_BUTT;
        case LINE_CAP_SQUARE: return CAIRO_LINE_CAP_SQUARE;
        case LINE_CAP_ROUND:
        default:              return CAIRO_LINE_CAP_ROUND;
    }
}

static void draw_line_ex(cairo_t* cr, point2d_t p1, point2d_t p2, double width, line_cap_t cap) {
    cairo_set_line_width(cr, width);
    cairo_set_line_cap(cr, get_cairo_line_cap(cap));
    cairo_move_to(cr, p1.x, p1.y);
    cairo_line_to(cr, p2.x, p2.y);
    cairo_stroke(cr);
}

/* Calculate bond gap endpoints for modern style (bond shortening at heteroatoms) */
void calculate_bond_gap(point2d_t p1, point2d_t p2,
                        bool gap_at_p1, bool gap_at_p2,
                        double gap_factor, double base_scale, double font_scale,
                        point2d_t* out_p1, point2d_t* out_p2) {
    double len = point2d_distance(p1, p2);
    if (len < 0.001) {
        *out_p1 = p1;
        *out_p2 = p2;
        return;
    }

    point2d_t dir = point2d_sub(p2, p1);
    dir = point2d_scale(dir, 1.0 / len);

    /* Gap size must be large enough to clear the atom label
     * Use the same scaling as render_atom_label_modern_ex */
    double fs = (font_scale > 0.0) ? font_scale : 1.0;
    double actual_font = BASE_FONT_SIZE * base_scale * 0.24 * fs;
    double gap_size = actual_font * 0.6 * gap_factor;

    double t1 = gap_at_p1 ? (gap_size / len) : 0.0;
    double t2 = gap_at_p2 ? (gap_size / len) : 0.0;

    /* Ensure we don't create too-short bonds */
    if (t1 + t2 >= 0.8) {
        double total = t1 + t2;
        t1 = t1 / total * 0.4;
        t2 = t2 / total * 0.4;
    }

    *out_p1 = point2d_add(p1, point2d_scale(dir, len * t1));
    *out_p2 = point2d_sub(p2, point2d_scale(dir, len * t2));
}

static void draw_text(cairo_t* cr, const char* text, point2d_t pos, double font_size) {
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, font_size);

    cairo_text_extents_t extents;
    cairo_text_extents(cr, text, &extents);

    cairo_move_to(cr,
                 pos.x - extents.width / 2 - extents.x_bearing,
                 pos.y - extents.height / 2 - extents.y_bearing);
    cairo_show_text(cr, text);
}

/* ============== Bond Rendering ============== */

static void render_wedge_bond(cairo_t* cr, point2d_t p1, point2d_t p2, double width) {
    rgb_color_t black = {0, 0, 0};
    set_color(cr, black);

    point2d_t dir = point2d_sub(p2, p1);
    point2d_t perp = point2d_normalize(point2d_perp(dir));
    double wedge_width = width * 5.0;

    point2d_t offset = point2d_scale(perp, wedge_width / 2.0);
    point2d_t p2a = point2d_add(p2, offset);
    point2d_t p2b = point2d_sub(p2, offset);

    cairo_move_to(cr, p1.x, p1.y);
    cairo_line_to(cr, p2a.x, p2a.y);
    cairo_line_to(cr, p2b.x, p2b.y);
    cairo_close_path(cr);
    cairo_fill(cr);
}

static void render_dash_bond(cairo_t* cr, point2d_t p1, point2d_t p2, double width) {
    rgb_color_t black = {0, 0, 0};
    set_color(cr, black);

    point2d_t dir = point2d_sub(p2, p1);
    double length = point2d_length(dir);
    if (length < 0.01) return;

    dir = point2d_scale(dir, 1.0 / length);
    point2d_t perp = point2d_perp(dir);

    int num_hashes = 8;
    for (int i = 0; i <= num_hashes; i++) {
        double t = (double)i / num_hashes;
        point2d_t center = point2d_add(p1, point2d_scale(dir, t * length));
        double hash_width = width * 2.0 * t + width;
        point2d_t offset = point2d_scale(perp, hash_width / 2.0);

        cairo_set_line_width(cr, width * 0.8);
        cairo_move_to(cr, center.x + offset.x, center.y + offset.y);
        cairo_line_to(cr, center.x - offset.x, center.y - offset.y);
        cairo_stroke(cr);
    }
}

/* ============== Atom Rendering ============== */

/* Get atom radius for rendering based on style and options */
static double get_atom_render_radius(const atom_t* atom, const depictor_options_t* opts,
                                      double base_scale) {
    double vdw = atom_get_vdw_radius(atom->element);
    double cov = atom_get_covalent_radius(atom->element);

    switch (opts->render_style) {
        case RENDER_STYLE_SPACEFILL:
            /* Full VDW radius */
            return vdw * opts->atom_radius_scale * base_scale;

        case RENDER_STYLE_BALLS_AND_STICKS:
            /* Proportional to VDW but smaller */
            if (opts->proportional_atoms) {
                return vdw * opts->atom_radius_scale * base_scale * 0.5;
            }
            return cov * opts->atom_radius_scale * base_scale * 0.6;

        case RENDER_STYLE_STICKS:
            /* Small caps, slightly proportional */
            if (opts->proportional_atoms && atom->element != ELEM_C) {
                return cov * opts->atom_radius_scale * base_scale * 0.4;
            }
            return opts->bond_width * 1.5;

        case RENDER_STYLE_SURFACE:
            /* VDW for surface calculation */
            return vdw * opts->atom_radius_scale * base_scale;

        case RENDER_STYLE_WIREFRAME:
        default:
            return 0;
    }
}

/* Check if atom is a terminal carbon (connected to only one heavy atom) */
static bool is_terminal_carbon(const molecule_t* mol, int atom_idx) {
    if (mol->atoms[atom_idx].element != ELEM_C) return false;

    int heavy_neighbors = 0;
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int other = -1;
        if (bond->atom1 == atom_idx) other = bond->atom2;
        else if (bond->atom2 == atom_idx) other = bond->atom1;

        if (other >= 0 && mol->atoms[other].element != ELEM_H) {
            heavy_neighbors++;
        }
    }
    return heavy_neighbors <= 1;
}

static bool should_show_label_ex(const atom_t* atom, const depictor_options_t* opts,
                                  const molecule_t* mol, int atom_idx) {
    /* In spacefill mode, labels are optional overlays */
    if (opts->render_style == RENDER_STYLE_SPACEFILL) {
        return opts->show_carbons || atom->element != ELEM_C;
    }
    if (opts->show_carbons) return true;
    if (atom->element != ELEM_C) return true;
    if (opts->show_hydrogens && atom->implicit_h_count > 0) return true;
    /* For modern style, optionally show terminal carbons as CH3 labels */
    if (opts->terminal_carbon_labels && mol && atom->element == ELEM_C) {
        if (is_terminal_carbon(mol, atom_idx) && atom->implicit_h_count >= 3) {
            return true;
        }
    }
    return false;
}

static bool should_show_label(const atom_t* atom, const depictor_options_t* opts) {
    return should_show_label_ex(atom, opts, NULL, -1);
}

/* Render atom as filled sphere with 3D shading effect */
static void render_atom_sphere(cairo_t* cr, const atom_t* atom, point2d_t pos,
                                double radius, const depictor_options_t* opts __attribute__((unused))) {
    rgb_color_t base_color = atom_get_color(atom->element);

    /* Create radial gradient for 3D effect */
    cairo_pattern_t* gradient = cairo_pattern_create_radial(
        pos.x - radius * 0.3, pos.y - radius * 0.3, radius * 0.1,
        pos.x, pos.y, radius
    );

    /* Highlight color (lighter) */
    double hr = fmin(1.0, base_color.r / 255.0 + 0.4);
    double hg = fmin(1.0, base_color.g / 255.0 + 0.4);
    double hb = fmin(1.0, base_color.b / 255.0 + 0.4);

    /* Shadow color (darker) */
    double sr = base_color.r / 255.0 * 0.6;
    double sg = base_color.g / 255.0 * 0.6;
    double sb = base_color.b / 255.0 * 0.6;

    cairo_pattern_add_color_stop_rgb(gradient, 0.0, hr, hg, hb);
    cairo_pattern_add_color_stop_rgb(gradient, 0.5, base_color.r / 255.0,
                                      base_color.g / 255.0, base_color.b / 255.0);
    cairo_pattern_add_color_stop_rgb(gradient, 1.0, sr, sg, sb);

    cairo_arc(cr, pos.x, pos.y, radius, 0, 2 * M_PI);
    cairo_set_source(cr, gradient);
    cairo_fill(cr);
    cairo_pattern_destroy(gradient);

    /* Subtle outline */
    cairo_arc(cr, pos.x, pos.y, radius, 0, 2 * M_PI);
    cairo_set_source_rgba(cr, 0, 0, 0, 0.3);
    cairo_set_line_width(cr, 0.5);
    cairo_stroke(cr);
}

/* Calculate average bond direction for an atom (for H label positioning) */
static point2d_t calculate_avg_bond_direction(const molecule_t* mol, int atom_idx,
                                               const mol_coords_t* coords) {
    point2d_t avg_dir = {0, 0};
    int bond_count = 0;

    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int other = -1;
        if (bond->atom1 == atom_idx) other = bond->atom2;
        else if (bond->atom2 == atom_idx) other = bond->atom1;

        if (other >= 0 && mol->atoms[other].element != ELEM_H) {
            point2d_t dir = point2d_sub(coords->coords_2d[other], coords->coords_2d[atom_idx]);
            avg_dir = point2d_add(avg_dir, dir);
            bond_count++;
        }
    }

    if (bond_count > 0) {
        avg_dir = point2d_scale(avg_dir, 1.0 / bond_count);
    }
    return avg_dir;
}

/* Format charge string for display */
static void format_charge_string(int charge, char* out, size_t out_size) {
    out[0] = '\0';
    if (charge == 1) {
        snprintf(out, out_size, "+");
    } else if (charge == -1) {
        snprintf(out, out_size, "-");
    } else if (charge > 1) {
        snprintf(out, out_size, "%d+", charge);
    } else if (charge < -1) {
        snprintf(out, out_size, "%d-", -charge);
    }
}

/* Render atom label in modern style (no circle background, just colored text)
 * Extended version with bond direction for H positioning */
static void render_atom_label_modern_ex(cairo_t* cr, const atom_t* atom, point2d_t pos,
                                         const depictor_options_t* opts,
                                         rgb_color_t bg __attribute__((unused)),
                                         const molecule_t* mol, int atom_idx,
                                         const mol_coords_t* coords, double base_scale) {
    const char* symbol = element_to_symbol(atom->element);
    char label[64];
    char charge_str[8];

    /* Determine if H should go on left or right based on bond direction */
    bool h_on_left = false;
    if (mol && coords && atom->implicit_h_count > 0) {
        point2d_t avg_dir = calculate_avg_bond_direction(mol, atom_idx, coords);
        /* If average bond direction points to the right (positive x), put H on left */
        h_on_left = (avg_dir.x > 0.1);
    }

    /* Format charge string */
    format_charge_string(atom->charge, charge_str, sizeof(charge_str));

    /* Build label with hydrogens if needed
     * Heteroatoms (non-C, non-H) should always show their H */
    bool is_heteroatom = (atom->element != ELEM_C && atom->element != ELEM_H);
    bool show_h = (is_heteroatom && atom->implicit_h_count > 0) ||
                  (opts->show_hydrogens && atom->implicit_h_count > 0) ||
                  (atom->element == ELEM_C && atom->implicit_h_count >= 3);

    if (show_h) {
        char h_part[16];
        if (atom->implicit_h_count == 1)
            snprintf(h_part, sizeof(h_part), "H");
        else
            snprintf(h_part, sizeof(h_part), "H%d", atom->implicit_h_count);

        if (h_on_left) {
            /* H before symbol: HO, HN, H2N, etc. */
            snprintf(label, sizeof(label), "%s%s%s", h_part, symbol, charge_str);
        } else {
            /* H after symbol: OH, NH, NH2, etc. */
            snprintf(label, sizeof(label), "%s%s%s", symbol, h_part, charge_str);
        }
    } else {
        snprintf(label, sizeof(label), "%s%s", symbol, charge_str);
    }

    /* Scale font size based on average bond length (base_scale)
     * This ensures consistent letter-to-bond proportions across molecule sizes */
    double fs = (opts->font_scale > 0.0) ? opts->font_scale : 1.0;
    double font_size = BASE_FONT_SIZE * base_scale * 0.24 * fs;

    cairo_text_extents_t extents;
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, font_size);
    cairo_text_extents(cr, label, &extents);

    /* No background - the bond gaps provide the clearance */

    /* Draw text - black by default, colored if colored_atoms is enabled */
    if (opts->colored_atoms) {
        rgb_color_t color = atom_get_color(atom->element);
        set_color(cr, color);
    } else {
        rgb_color_t black = {0, 0, 0};
        set_color(cr, black);
    }

    /* Center the text properly */
    cairo_move_to(cr,
                 pos.x - extents.width / 2.0 - extents.x_bearing,
                 pos.y - extents.height / 2.0 - extents.y_bearing);
    cairo_show_text(cr, label);
}

/* Render atom label in modern style (no circle background, just colored text) */
static void render_atom_label_modern(cairo_t* cr, const atom_t* atom, point2d_t pos,
                                      const depictor_options_t* opts,
                                      rgb_color_t bg, double base_scale) {
    /* Fallback without bond direction info - H goes on right */
    render_atom_label_modern_ex(cr, atom, pos, opts, bg, NULL, -1, NULL, base_scale);
}

static void render_atom_label(cairo_t* cr, const atom_t* atom, point2d_t pos,
                              const depictor_options_t* opts, rgb_color_t bg,
                              double base_scale) {
    bool show_sphere = (opts->render_style == RENDER_STYLE_BALLS_AND_STICKS ||
                        opts->render_style == RENDER_STYLE_SPACEFILL ||
                        opts->atom_filling);

    double radius = get_atom_render_radius(atom, opts, base_scale);

    if (show_sphere && radius > 0) {
        render_atom_sphere(cr, atom, pos, radius, opts);

        /* Show label on top for non-carbon or if requested */
        if (should_show_label(atom, opts) && opts->render_style != RENDER_STYLE_SPACEFILL) {
            double font_size = BASE_FONT_SIZE * opts->font_scale * 8.0;
            const char* symbol = element_to_symbol(atom->element);
            char label[64];
            char charge_str[8];

            format_charge_string(atom->charge, charge_str, sizeof(charge_str));

            /* Heteroatoms should always show their H */
            bool is_heteroatom = (atom->element != ELEM_C && atom->element != ELEM_H);
            bool show_h = (is_heteroatom && atom->implicit_h_count > 0) ||
                          (opts->show_hydrogens && atom->implicit_h_count > 0);

            if (show_h) {
                if (atom->implicit_h_count == 1)
                    snprintf(label, sizeof(label), "%sH%s", symbol, charge_str);
                else
                    snprintf(label, sizeof(label), "%sH%d%s", symbol, atom->implicit_h_count, charge_str);
            } else {
                snprintf(label, sizeof(label), "%s%s", symbol, charge_str);
            }

            /* White text on colored sphere */
            cairo_set_source_rgb(cr, 1, 1, 1);
            draw_text(cr, label, pos, font_size);
        }
        return;
    }

    /* Check if we should show this label */
    if (!should_show_label(atom, opts)) return;

    /* Modern style: no circle, just text with background */
    if (opts->style_preset == DEPICT_STYLE_MODERN) {
        render_atom_label_modern(cr, atom, pos, opts, bg, base_scale);
        return;
    }

    /* Default style: circle background */
    const char* symbol = element_to_symbol(atom->element);
    char label[64];
    char charge_str[8];

    format_charge_string(atom->charge, charge_str, sizeof(charge_str));

    /* Heteroatoms should always show their H */
    bool is_heteroatom = (atom->element != ELEM_C && atom->element != ELEM_H);
    bool show_h = (is_heteroatom && atom->implicit_h_count > 0) ||
                  (opts->show_hydrogens && atom->implicit_h_count > 0);

    if (show_h) {
        if (atom->implicit_h_count == 1)
            snprintf(label, sizeof(label), "%sH%s", symbol, charge_str);
        else
            snprintf(label, sizeof(label), "%sH%d%s", symbol, atom->implicit_h_count, charge_str);
    } else {
        snprintf(label, sizeof(label), "%s%s", symbol, charge_str);
    }

    double font_size = BASE_FONT_SIZE * opts->font_scale * 10.0;

    /* Background circle */
    cairo_text_extents_t extents;
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, font_size);
    cairo_text_extents(cr, label, &extents);

    double label_radius = (extents.width > extents.height ? extents.width : extents.height) / 2.0 + 4;

    set_color(cr, bg);
    cairo_arc(cr, pos.x, pos.y, label_radius, 0, 2 * M_PI);
    cairo_fill(cr);

    cairo_set_line_width(cr, 0.8);
    cairo_arc(cr, pos.x, pos.y, label_radius, 0, 2 * M_PI);
    rgb_color_t outline = {0, 0, 0};
    set_color(cr, outline);
    cairo_stroke(cr);

    /* Draw text - black by default, colored if colored_atoms is enabled */
    if (opts->colored_atoms) {
        rgb_color_t color = atom_get_color(atom->element);
        set_color(cr, color);
    } else {
        rgb_color_t black = {0, 0, 0};
        set_color(cr, black);
    }
    draw_text(cr, label, pos, font_size);
}

/* Extended version that handles terminal carbons for modern style */
static void render_atom_label_ex(cairo_t* cr, const molecule_t* mol, int atom_idx,
                                  point2d_t pos, const depictor_options_t* opts,
                                  rgb_color_t bg, double base_scale,
                                  const mol_coords_t* coords) {
    const atom_t* atom = &mol->atoms[atom_idx];

    bool show_sphere = (opts->render_style == RENDER_STYLE_BALLS_AND_STICKS ||
                        opts->render_style == RENDER_STYLE_SPACEFILL ||
                        opts->atom_filling);

    double radius = get_atom_render_radius(atom, opts, base_scale);

    if (show_sphere && radius > 0) {
        render_atom_sphere(cr, atom, pos, radius, opts);
        if (should_show_label_ex(atom, opts, mol, atom_idx) &&
            opts->render_style != RENDER_STYLE_SPACEFILL) {
            double font_size = BASE_FONT_SIZE * opts->font_scale * 8.0;
            const char* symbol = element_to_symbol(atom->element);
            char label[64];
            char charge_str[8];

            /* Format charge string */
            format_charge_string(atom->charge, charge_str, sizeof(charge_str));

            /* Heteroatoms should always show their H */
            bool is_heteroatom = (atom->element != ELEM_C && atom->element != ELEM_H);
            bool show_h = (is_heteroatom && atom->implicit_h_count > 0) ||
                          (opts->show_hydrogens && atom->implicit_h_count > 0);

            if (show_h) {
                if (atom->implicit_h_count == 1)
                    snprintf(label, sizeof(label), "%sH%s", symbol, charge_str);
                else
                    snprintf(label, sizeof(label), "%sH%d%s", symbol, atom->implicit_h_count, charge_str);
            } else {
                snprintf(label, sizeof(label), "%s%s", symbol, charge_str);
            }

            cairo_set_source_rgb(cr, 1, 1, 1);
            draw_text(cr, label, pos, font_size);
        }
        return;
    }

    if (!should_show_label_ex(atom, opts, mol, atom_idx)) return;

    if (opts->style_preset == DEPICT_STYLE_MODERN) {
        render_atom_label_modern_ex(cr, atom, pos, opts, bg, mol, atom_idx, coords, base_scale);
        return;
    }

    /* Fallback to standard label */
    render_atom_label(cr, atom, pos, opts, bg, base_scale);
}

/* ============== Stick Bond Rendering (colored by element) ============== */

static void render_stick_bond_ex(cairo_t* cr, const molecule_t* mol, const bond_t* bond,
                                  point2d_t p1, point2d_t p2, double width, line_cap_t cap) {
    /* Color each half of the bond by the atom's element color */
    point2d_t mid = {(p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0};

    rgb_color_t c1 = atom_get_color(mol->atoms[bond->atom1].element);
    rgb_color_t c2 = atom_get_color(mol->atoms[bond->atom2].element);

    cairo_set_line_width(cr, width);
    cairo_set_line_cap(cr, get_cairo_line_cap(cap));

    /* First half */
    set_color(cr, c1);
    cairo_move_to(cr, p1.x, p1.y);
    cairo_line_to(cr, mid.x, mid.y);
    cairo_stroke(cr);

    /* Second half */
    set_color(cr, c2);
    cairo_move_to(cr, mid.x, mid.y);
    cairo_line_to(cr, p2.x, p2.y);
    cairo_stroke(cr);
}

/* ============== Surface Rendering (Smooth molecular surface) ============== */

/*
 * Smooth blobby implicit surface using Gaussian metaballs.
 * The field function uses a smooth exponential falloff to create
 * a continuous surface that blends atoms together seamlessly.
 */

/* Gaussian-based field contribution for smooth blending */
static inline double gaussian_field(double d2, double r2) {
    /* Gaussian falloff: exp(-k * d^2 / r^2)
     * Lower k = more blending, atoms merge together
     * Higher k = more distinct atoms */
    double k = 0.8;  /* Low value for smooth blobby surface that fills rings */
    return exp(-k * d2 / r2);
}

/* Get weighted-average color from nearby atoms for atom coloring mode */
static void get_atom_weighted_color(double px, double py, const molecule_t* mol,
                                     const mol_coords_t* coords, const double* radii,
                                     double* r, double* g, double* b) {
    double total_weight = 0.0;
    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        double dx = px - coords->coords_2d[i].x;
        double dy = py - coords->coords_2d[i].y;
        double d2 = dx * dx + dy * dy;
        double r2 = radii[i] * radii[i];
        double weight = gaussian_field(d2, r2);

        if (weight > 0.001) {
            rgb_color_t c = atom_get_color(mol->atoms[i].element);
            sum_r += c.r * weight;
            sum_g += c.g * weight;
            sum_b += c.b * weight;
            total_weight += weight;
        }
    }

    if (total_weight > 0) {
        *r = sum_r / total_weight / 255.0;
        *g = sum_g / total_weight / 255.0;
        *b = sum_b / total_weight / 255.0;
    } else {
        /* Fallback to gray */
        *r = *g = *b = 0.5;
    }
}

/* Get Gasteiger partial charge for an atom (simplified estimation) */
static double estimate_partial_charge(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    double charge = 0.0;

    /* Base electronegativity differences */
    switch (atom->element) {
        case ELEM_O: charge = -0.4; break;  /* Electronegative */
        case ELEM_N: charge = -0.3; break;
        case ELEM_S: charge = -0.2; break;
        case ELEM_F: charge = -0.5; break;
        case ELEM_Cl: charge = -0.3; break;
        case ELEM_Br: charge = -0.2; break;
        case ELEM_C: charge = 0.0; break;
        case ELEM_H: charge = 0.1; break;
        default: charge = 0.0; break;
    }

    /* Adjust based on neighbors */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int other = -1;
        if (bond->atom1 == atom_idx) other = bond->atom2;
        else if (bond->atom2 == atom_idx) other = bond->atom1;

        if (other >= 0) {
            element_t neighbor = mol->atoms[other].element;
            /* Pull electron density toward electronegative atoms */
            if (neighbor == ELEM_O || neighbor == ELEM_N || neighbor == ELEM_F) {
                charge += 0.15;  /* More positive due to electron withdrawal */
            }
        }
    }

    return fmax(-1.0, fmin(1.0, charge));
}

/* Get polarity color: red for negative, white for neutral, blue for positive */
static void get_polarity_color(double px, double py, const molecule_t* mol,
                                const mol_coords_t* coords, const double* radii,
                                double* r, double* g, double* b) {
    double total_weight = 0.0;
    double weighted_charge = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        double dx = px - coords->coords_2d[i].x;
        double dy = py - coords->coords_2d[i].y;
        double d2 = dx * dx + dy * dy;
        double r2 = radii[i] * radii[i];
        double weight = gaussian_field(d2, r2);

        if (weight > 0.001) {
            double charge = estimate_partial_charge(mol, i);
            weighted_charge += charge * weight;
            total_weight += weight;
        }
    }

    double charge = (total_weight > 0) ? weighted_charge / total_weight : 0.0;

    /* Map charge to color: negative=red, neutral=white, positive=blue */
    if (charge < 0) {
        /* Red to white gradient for negative charge */
        double t = -charge;  /* 0 to 1 */
        *r = 1.0;
        *g = 1.0 - t * 0.7;
        *b = 1.0 - t * 0.7;
    } else {
        /* White to blue gradient for positive charge */
        double t = charge;  /* 0 to 1 */
        *r = 1.0 - t * 0.7;
        *g = 1.0 - t * 0.4;
        *b = 1.0;
    }
}

/* Compute ambient occlusion approximation (darker in cavities) */
static double compute_occlusion(double px, double py, const point2d_t* centers,
                                 const double* radii, int n, double probe_r) {
    /* Count how "buried" this point is by nearby atoms */
    double occlusion = 0.0;
    for (int i = 0; i < n; i++) {
        double dx = px - centers[i].x;
        double dy = py - centers[i].y;
        double d = sqrt(dx * dx + dy * dy);
        double r = radii[i] + probe_r;
        if (d < r) {
            /* Smooth falloff for occlusion contribution */
            double t = (r - d) / r;
            occlusion += t * t;  /* Quadratic falloff for softer shadows */
        }
    }
    /* Normalize and clamp - reduced effect to avoid dark spots */
    return fmin(0.6, occlusion * 0.15);
}

static void render_molecular_surface(cairo_t* cr, const molecule_t* mol,
                                      const mol_coords_t* coords,
                                      const depictor_options_t* opts,
                                      double base_scale) {
    if (mol->num_atoms == 0) return;

    /* Precompute atom radii and bounding box */
    double* radii = malloc(mol->num_atoms * sizeof(double));
    if (!radii) return;

    double min_x = 1e9, min_y = 1e9, max_x = -1e9, max_y = -1e9;
    double max_radius = 0;

    for (int i = 0; i < mol->num_atoms; i++) {
        double vdw = atom_get_vdw_radius(mol->atoms[i].element);
        radii[i] = vdw * opts->atom_radius_scale * base_scale;
        if (radii[i] > max_radius) max_radius = radii[i];

        double x = coords->coords_2d[i].x;
        double y = coords->coords_2d[i].y;
        double r = radii[i] * 2.0;
        if (x - r < min_x) min_x = x - r;
        if (y - r < min_y) min_y = y - r;
        if (x + r > max_x) max_x = x + r;
        if (y + r > max_y) max_y = y + r;
    }

    /* Calculate molecule center and extent */
    double cx = 0, cy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        cx += coords->coords_2d[i].x;
        cy += coords->coords_2d[i].y;
    }
    cx /= mol->num_atoms;
    cy /= mol->num_atoms;

    double mol_extent = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        double dx = coords->coords_2d[i].x - cx;
        double dy = coords->coords_2d[i].y - cy;
        double r = sqrt(dx * dx + dy * dy) + radii[i];
        if (r > mol_extent) mol_extent = r;
    }

    /* Grid sampling parameters */
    double grid_step = 1.0;
    int grid_w = (int)((max_x - min_x) / grid_step) + 4;
    int grid_h = (int)((max_y - min_y) / grid_step) + 4;

    /* Limit grid size */
    if (grid_w > 1000) { grid_step *= grid_w / 1000.0; grid_w = 1000; }
    if (grid_h > 1000) { grid_step *= grid_h / 1000.0; grid_h = 1000; }

    /* Compute field values and gradients */
    double* field = calloc(grid_w * grid_h, sizeof(double));
    if (!field) { free(radii); return; }

    /* Field computation with Gaussian blending */
    for (int gy = 0; gy < grid_h; gy++) {
        for (int gx = 0; gx < grid_w; gx++) {
            double px = min_x + gx * grid_step;
            double py = min_y + gy * grid_step;

            double sum = 0.0;
            for (int i = 0; i < mol->num_atoms; i++) {
                double dx = px - coords->coords_2d[i].x;
                double dy = py - coords->coords_2d[i].y;
                double d2 = dx * dx + dy * dy;
                double r2 = radii[i] * radii[i];
                sum += gaussian_field(d2, r2);
            }
            field[gy * grid_w + gx] = sum;
        }
    }

    /* Surface threshold - determines where the surface is drawn
     * Lower value = larger surface that fills ring centers
     * Higher value = tighter surface around atoms */
    double threshold = 0.25;

    /* Draw the surface with smooth shading */
    for (int gy = 0; gy < grid_h; gy++) {
        for (int gx = 0; gx < grid_w; gx++) {
            double v = field[gy * grid_w + gx];
            if (v < threshold) continue;

            double px = min_x + gx * grid_step;
            double py = min_y + gy * grid_step;

            /* Get base color based on surface coloring mode */
            double base_r, base_g, base_b;
            switch (opts->surface_color) {
                case SURFACE_COLOR_ATOM:
                    get_atom_weighted_color(px, py, mol, coords, radii, &base_r, &base_g, &base_b);
                    break;
                case SURFACE_COLOR_POLARITY:
                    get_polarity_color(px, py, mol, coords, radii, &base_r, &base_g, &base_b);
                    break;
                case SURFACE_COLOR_UNIFORM:
                default:
                    base_r = 0.15; base_g = 0.50; base_b = 0.85;
                    break;
            }

            /* Compute gradient for normal estimation */
            double gx_grad = 0, gy_grad = 0;
            if (gx > 0 && gx < grid_w - 1) {
                gx_grad = field[gy * grid_w + gx + 1] - field[gy * grid_w + gx - 1];
            }
            if (gy > 0 && gy < grid_h - 1) {
                gy_grad = field[(gy + 1) * grid_w + gx] - field[(gy - 1) * grid_w + gx];
            }

            /* Compute normal from gradient - always use gradient for 3D effect
             * Scale factor controls how much the surface "bulges" */
            double bulge_scale = 3.0;  /* Higher = more pronounced 3D */
            double nx = -gx_grad * bulge_scale;
            double ny = -gy_grad * bulge_scale;
            double nz = 1.0;

            /* Normalize the normal vector */
            double n_len = sqrt(nx * nx + ny * ny + nz * nz);
            nx /= n_len; ny /= n_len; nz /= n_len;

            /* Light direction (from upper-left, towards viewer) */
            double lx = -0.6, ly = -0.6, lz = 0.5;
            double l_len = sqrt(lx * lx + ly * ly + lz * lz);
            lx /= l_len; ly /= l_len; lz /= l_len;

            /* Diffuse lighting with wrap-around for softer shadows */
            double ndotl = nx * lx + ny * ly + nz * lz;
            double diffuse = fmax(0.0, ndotl);
            double wrap_diffuse = (ndotl + 0.5) / 1.5;  /* Soft wrap lighting */
            wrap_diffuse = fmax(0.0, wrap_diffuse);

            /* Specular highlight (Blinn-Phong) */
            double hx = lx, hy = ly, hz = lz + 1.0;  /* View is (0,0,1) */
            double h_len = sqrt(hx * hx + hy * hy + hz * hz);
            hx /= h_len; hy /= h_len; hz /= h_len;
            double spec_dot = fmax(0.0, nx * hx + ny * hy + nz * hz);
            double specular = pow(spec_dot, 30.0);

            /* Ambient occlusion - darken in cavities */
            double occlusion = compute_occlusion(px, py, coords->coords_2d,
                                                  radii, mol->num_atoms, max_radius * 0.3);

            /* Edge darkening based on proximity to boundary */
            double edge_factor = 1.0;
            if (v < threshold * 1.5) {
                double t = (v - threshold) / (threshold * 0.5);
                t = fmax(0.0, fmin(1.0, t));
                edge_factor = 0.6 + 0.4 * t;  /* Darken edges */
            }

            /* Combine lighting components with more contrast */
            double ambient = 0.25;
            double light = ambient + 0.45 * wrap_diffuse + 0.3 * diffuse;
            light *= (1.0 - occlusion * 0.4);
            light *= edge_factor;
            light = fmin(1.0, fmax(0.15, light));  /* Keep minimum brightness */

            /* Final color with lighting applied */
            double r = base_r * light + specular * 0.6;
            double g = base_g * light + specular * 0.6;
            double b = base_b * light + specular * 0.5;

            /* Fresnel effect - rim lighting on edges */
            double fresnel = 1.0 - nz;
            fresnel = pow(fresnel, 2.0) * 0.2;
            r += fresnel;
            g += fresnel;
            b += fresnel;

            r = fmin(1.0, fmax(0.0, r));
            g = fmin(1.0, fmax(0.0, g));
            b = fmin(1.0, fmax(0.0, b));

            cairo_set_source_rgb(cr, r, g, b);
            cairo_rectangle(cr, px, py, grid_step + 0.5, grid_step + 0.5);
            cairo_fill(cr);
        }
    }

    /* No explicit contour lines - the surface shading provides edge definition
     * through the edge_factor darkening near the threshold boundary */

    free(field);
    free(radii);
}

/* ============== Main Rendering ============== */

/* Helper: check if an atom should have a bond gap (is labeled) */
static bool needs_bond_gap(const molecule_t* mol, int atom_idx, const depictor_options_t* opts) {
    if (opts->heteroatom_gap <= 0.0) return false;

    const atom_t* atom = &mol->atoms[atom_idx];

    /* Heteroatoms always get gaps */
    if (atom->element != ELEM_C) return true;

    /* Terminal carbons with labels get gaps */
    if (opts->terminal_carbon_labels && is_terminal_carbon(mol, atom_idx) &&
        atom->implicit_h_count >= 3) {
        return true;
    }

    /* Show carbons option */
    if (opts->show_carbons) return true;

    return false;
}

cchem_status_t render_molecule(render_context_t* ctx, const molecule_t* mol,
                               const mol_coords_t* coords, const depictor_options_t* opts) {
    if (!ctx || !mol || !coords || !opts) return CCHEM_ERROR_INVALID_INPUT;
    if (!coords->has_2d) return CCHEM_ERROR_INVALID_INPUT;

    /* Calculate base scale from average bond length */
    double avg_bond_len = 0.0;
    if (mol->num_bonds > 0) {
        for (int b = 0; b < mol->num_bonds; b++) {
            avg_bond_len += point2d_distance(coords->coords_2d[mol->bonds[b].atom1],
                                             coords->coords_2d[mol->bonds[b].atom2]);
        }
        avg_bond_len /= mol->num_bonds;
    } else {
        avg_bond_len = opts->bond_length;
    }
    double base_scale = avg_bond_len / 1.5;  /* Scale relative to ~1.5 Angstrom avg bond */

    /* Surface mode renders differently */
    if (opts->render_style == RENDER_STYLE_SURFACE) {
        render_molecular_surface(ctx->cr, mol, coords, opts, base_scale);
        return CCHEM_OK;
    }

    /* Determine bond width based on style */
    double bond_width = opts->bond_width;
    if (opts->render_style == RENDER_STYLE_BALLS_AND_STICKS) {
        bond_width = opts->bond_width * 3.0;  /* Thick bonds for balls-and-sticks */
    } else if (opts->render_style == RENDER_STYLE_STICKS) {
        bond_width = opts->bond_width * 2.5;  /* Medium-thick for sticks */
    }

    /* Always apply bond gaps at heteroatoms when gap > 0 */
    bool apply_gaps = (opts->heteroatom_gap > 0.0);

    /* Draw bonds first (unless spacefill mode) */
    if (opts->render_style != RENDER_STYLE_SPACEFILL) {
        for (int b = 0; b < mol->num_bonds; b++) {
            const bond_t* bond = &mol->bonds[b];
            int i = bond->atom1;
            int j = bond->atom2;

            point2d_t p1 = coords->coords_2d[i];
            point2d_t p2 = coords->coords_2d[j];

            /* Calculate gaps at heteroatoms */
            bool gap1 = apply_gaps && needs_bond_gap(mol, i, opts);
            bool gap2 = apply_gaps && needs_bond_gap(mol, j, opts);
            point2d_t gp1, gp2;
            calculate_bond_gap(p1, p2, gap1, gap2, opts->heteroatom_gap, base_scale, opts->font_scale, &gp1, &gp2);

            int order = bond_get_int_order(bond);
            bool use_colored_sticks = (opts->render_style == RENDER_STYLE_STICKS ||
                                       opts->render_style == RENDER_STYLE_BALLS_AND_STICKS);

            double dbl_offset = opts->double_bond_offset * opts->bond_length;
            double trpl_offset = opts->triple_bond_offset * opts->bond_length;

            /* Check stereochemistry - only draw wedge/dash for tetrahedral chirality,
             * not for E/Z double bond indicators (which also use BOND_UP/DOWN) */
            bool is_chiral_bond = (mol->atoms[i].chirality != CHIRALITY_NONE ||
                                   mol->atoms[j].chirality != CHIRALITY_NONE);
            if (bond->stereo_type == BOND_UP && is_chiral_bond) {
                point2d_t stereo_gp1 = gp1, stereo_gp2 = gp2;
                if (bond->stereo_atom == j) {
                    stereo_gp1 = gp2; stereo_gp2 = gp1;
                }
                render_wedge_bond(ctx->cr, stereo_gp1, stereo_gp2, bond_width);
            }
            else if (bond->stereo_type == BOND_DOWN && is_chiral_bond) {
                point2d_t stereo_gp1 = gp1, stereo_gp2 = gp2;
                if (bond->stereo_atom == j) {
                    stereo_gp1 = gp2; stereo_gp2 = gp1;
                }
                render_dash_bond(ctx->cr, stereo_gp1, stereo_gp2, bond_width);
            }
            else if (use_colored_sticks && order == 1) {
                render_stick_bond_ex(ctx->cr, mol, bond, gp1, gp2, bond_width, opts->line_cap);
            }
            else if (order == 1) {
                rgb_color_t black = {0, 0, 0};
                set_color(ctx->cr, black);
                draw_line_ex(ctx->cr, gp1, gp2, bond_width, opts->line_cap);
            }
            else if (order == 2) {
                point2d_t dir = point2d_sub(p2, p1);
                point2d_t perp = point2d_normalize(point2d_perp(dir));

                /* For ring bonds, offset inner line toward ring center only */
                if (bond->in_ring && bond->num_rings > 0 && mol->rings_computed) {
                    /* Find ring center */
                    int ring_idx = bond->ring_ids[0];
                    point2d_t ring_center = {0, 0};
                    for (int ri = 0; ri < mol->rings[ring_idx].size; ri++) {
                        ring_center = point2d_add(ring_center, coords->coords_2d[mol->rings[ring_idx].atoms[ri]]);
                    }
                    ring_center = point2d_scale(ring_center, 1.0 / mol->rings[ring_idx].size);

                    /* Bond midpoint */
                    point2d_t bond_mid = point2d_scale(point2d_add(p1, p2), 0.5);

                    /* Determine which perpendicular direction points toward ring center */
                    point2d_t to_center = point2d_sub(ring_center, bond_mid);
                    double dot = perp.x * to_center.x + perp.y * to_center.y;
                    point2d_t inward = (dot > 0) ? perp : point2d_scale(perp, -1.0);

                    /* Outer line on bond edge, inner line offset toward center */
                    point2d_t inner_offset = point2d_scale(inward, dbl_offset * 2.4);
                    point2d_t gp1_inner = point2d_add(gp1, inner_offset);
                    point2d_t gp2_inner = point2d_add(gp2, inner_offset);

                    /* Shorten inner line slightly for better appearance */
                    point2d_t shorten = point2d_scale(dir, 0.12);
                    gp1_inner = point2d_add(gp1_inner, shorten);
                    gp2_inner = point2d_sub(gp2_inner, shorten);

                    if (use_colored_sticks) {
                        render_stick_bond_ex(ctx->cr, mol, bond, gp1, gp2, bond_width, opts->line_cap);
                        render_stick_bond_ex(ctx->cr, mol, bond, gp1_inner, gp2_inner, bond_width * 0.7, opts->line_cap);
                    } else {
                        rgb_color_t black = {0, 0, 0};
                        set_color(ctx->cr, black);
                        draw_line_ex(ctx->cr, gp1, gp2, bond_width, opts->line_cap);
                        draw_line_ex(ctx->cr, gp1_inner, gp2_inner, bond_width, opts->line_cap);
                    }
                } else {
                    /* Non-ring double bond: symmetric offset */
                    point2d_t offset = point2d_scale(perp, dbl_offset);

                    point2d_t gp1a = point2d_add(gp1, offset);
                    point2d_t gp2a = point2d_add(gp2, offset);
                    point2d_t gp1b = point2d_sub(gp1, offset);
                    point2d_t gp2b = point2d_sub(gp2, offset);

                    if (use_colored_sticks) {
                        render_stick_bond_ex(ctx->cr, mol, bond, gp1a, gp2a, bond_width * 0.7, opts->line_cap);
                        render_stick_bond_ex(ctx->cr, mol, bond, gp1b, gp2b, bond_width * 0.7, opts->line_cap);
                    } else {
                        rgb_color_t black = {0, 0, 0};
                        set_color(ctx->cr, black);
                        draw_line_ex(ctx->cr, gp1a, gp2a, bond_width, opts->line_cap);
                        draw_line_ex(ctx->cr, gp1b, gp2b, bond_width, opts->line_cap);
                    }
                }
            }
            else if (order == 3) {
                point2d_t dir = point2d_sub(p2, p1);
                point2d_t perp = point2d_normalize(point2d_perp(dir));
                point2d_t offset = point2d_scale(perp, trpl_offset);

                point2d_t gp1a = point2d_add(gp1, offset);
                point2d_t gp2a = point2d_add(gp2, offset);
                point2d_t gp1b = point2d_sub(gp1, offset);
                point2d_t gp2b = point2d_sub(gp2, offset);

                if (use_colored_sticks) {
                    render_stick_bond_ex(ctx->cr, mol, bond, gp1, gp2, bond_width * 0.6, opts->line_cap);
                    render_stick_bond_ex(ctx->cr, mol, bond, gp1a, gp2a, bond_width * 0.6, opts->line_cap);
                    render_stick_bond_ex(ctx->cr, mol, bond, gp1b, gp2b, bond_width * 0.6, opts->line_cap);
                } else {
                    rgb_color_t black = {0, 0, 0};
                    set_color(ctx->cr, black);
                    draw_line_ex(ctx->cr, gp1, gp2, bond_width, opts->line_cap);
                    draw_line_ex(ctx->cr, gp1a, gp2a, bond_width, opts->line_cap);
                    draw_line_ex(ctx->cr, gp1b, gp2b, bond_width, opts->line_cap);
                }
            }
            else {
                if (use_colored_sticks) {
                    render_stick_bond_ex(ctx->cr, mol, bond, gp1, gp2, bond_width, opts->line_cap);
                } else {
                    rgb_color_t black = {0, 0, 0};
                    set_color(ctx->cr, black);
                    draw_line_ex(ctx->cr, gp1, gp2, bond_width, opts->line_cap);
                }
            }
        }
    }

    /* Aromatic circles */
    if (opts->draw_aromatic_circles && mol->rings_computed) {
        for (int r = 0; r < mol->num_rings; r++) {
            if (!mol->rings[r].aromatic) continue;

            point2d_t center = {0, 0};
            for (int ri = 0; ri < mol->rings[r].size; ri++) {
                center = point2d_add(center, coords->coords_2d[mol->rings[r].atoms[ri]]);
            }
            center = point2d_scale(center, 1.0 / mol->rings[r].size);

            double radius = point2d_distance(center, coords->coords_2d[mol->rings[r].atoms[0]]) * 0.6;

            cairo_set_line_width(ctx->cr, bond_width);
            cairo_set_line_cap(ctx->cr, get_cairo_line_cap(opts->line_cap));
            cairo_arc(ctx->cr, center.x, center.y, radius, 0, 2 * M_PI);
            rgb_color_t black = {0, 0, 0};
            set_color(ctx->cr, black);
            cairo_stroke(ctx->cr);
        }
    }

    /* Atom labels/spheres - use extended version for modern style */
    for (int i = 0; i < mol->num_atoms; i++) {
        render_atom_label_ex(ctx->cr, mol, i, coords->coords_2d[i],
                            opts, ctx->background, base_scale, coords);
    }

    return CCHEM_OK;
}

/* ============== File Output ============== */

cchem_status_t render_save_png(render_context_t* ctx, const char* filename) {
    if (!ctx || !filename) return CCHEM_ERROR_INVALID_INPUT;

    cairo_surface_flush(ctx->surface);
    cairo_status_t status = cairo_surface_write_to_png(ctx->surface, filename);
    return (status == CAIRO_STATUS_SUCCESS) ? CCHEM_OK : CCHEM_ERROR_FILE_IO;
}

cchem_status_t render_save_jpeg(render_context_t* ctx, const char* filename, int quality) {
    if (!ctx || !filename) return CCHEM_ERROR_INVALID_INPUT;

    cairo_surface_flush(ctx->surface);

    char temp_png[1024];
    snprintf(temp_png, sizeof(temp_png), "%s.temp.png", filename);

    cairo_status_t status = cairo_surface_write_to_png(ctx->surface, temp_png);
    if (status != CAIRO_STATUS_SUCCESS) {
        return CCHEM_ERROR_FILE_IO;
    }

    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "magick \"%s\" -quality %d \"%s\" && rm \"%s\"",
            temp_png, quality, filename, temp_png);
    int ret = system(cmd);

    return (ret == 0) ? CCHEM_OK : CCHEM_ERROR_FILE_IO;
}

cchem_status_t render_save_svg(render_context_t* ctx, const char* filename) {
    if (!ctx || !filename) return CCHEM_ERROR_INVALID_INPUT;

    /* For SVG surfaces, the file is already written during rendering */
    if (ctx->format == IMG_FORMAT_SVG) {
        cairo_surface_flush(ctx->surface);
        cairo_surface_finish(ctx->surface);
        return CCHEM_OK;
    }

    /* For non-SVG contexts, we need to create a new SVG surface and re-render */
    /* This is a fallback - ideally SVG should be created with create_ex */
    return CCHEM_ERROR_INVALID_INPUT;
}
