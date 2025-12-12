/**
 * @file types.c
 * @brief Implementation of core depiction types
 */

#include "cchem/depictor/types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============== Coordinate Allocation ============== */

mol_coords_t* mol_coords_create(int num_atoms) {
    if (num_atoms <= 0) return NULL;

    mol_coords_t* coords = calloc(1, sizeof(mol_coords_t));
    if (!coords) return NULL;

    coords->num_atoms = num_atoms;
    coords->coords_2d = calloc(num_atoms, sizeof(point2d_t));
    coords->coords_3d = calloc(num_atoms, sizeof(point3d_t));

    if (!coords->coords_2d || !coords->coords_3d) {
        mol_coords_free(coords);
        return NULL;
    }

    return coords;
}

void mol_coords_free(mol_coords_t* coords) {
    if (!coords) return;
    free(coords->coords_2d);
    free(coords->coords_3d);
    free(coords);
}

mol_coords_t* mol_coords_copy(const mol_coords_t* coords) {
    if (!coords) return NULL;

    mol_coords_t* copy = mol_coords_create(coords->num_atoms);
    if (!copy) return NULL;

    if (coords->has_2d && coords->coords_2d) {
        memcpy(copy->coords_2d, coords->coords_2d,
               coords->num_atoms * sizeof(point2d_t));
        copy->has_2d = true;
    }

    if (coords->has_3d && coords->coords_3d) {
        memcpy(copy->coords_3d, coords->coords_3d,
               coords->num_atoms * sizeof(point3d_t));
        copy->has_3d = true;
    }

    return copy;
}

/* ============== 2D Point Operations ============== */

point2d_t point2d_add(point2d_t a, point2d_t b) {
    return (point2d_t){a.x + b.x, a.y + b.y};
}

point2d_t point2d_sub(point2d_t a, point2d_t b) {
    return (point2d_t){a.x - b.x, a.y - b.y};
}

point2d_t point2d_scale(point2d_t p, double s) {
    return (point2d_t){p.x * s, p.y * s};
}

double point2d_length(point2d_t p) {
    return sqrt(p.x * p.x + p.y * p.y);
}

double point2d_distance(point2d_t a, point2d_t b) {
    return point2d_length(point2d_sub(a, b));
}

point2d_t point2d_normalize(point2d_t p) {
    double len = point2d_length(p);
    if (len < 1e-10) return (point2d_t){1, 0};
    return point2d_scale(p, 1.0 / len);
}

double point2d_dot(point2d_t a, point2d_t b) {
    return a.x * b.x + a.y * b.y;
}

double point2d_cross(point2d_t a, point2d_t b) {
    return a.x * b.y - a.y * b.x;
}

point2d_t point2d_rotate(point2d_t p, double angle) {
    double c = cos(angle), s = sin(angle);
    return (point2d_t){p.x * c - p.y * s, p.x * s + p.y * c};
}

point2d_t point2d_perp(point2d_t p) {
    return (point2d_t){-p.y, p.x};
}

/* ============== 3D Point Operations ============== */

point3d_t point3d_add(point3d_t a, point3d_t b) {
    return (point3d_t){a.x + b.x, a.y + b.y, a.z + b.z};
}

point3d_t point3d_sub(point3d_t a, point3d_t b) {
    return (point3d_t){a.x - b.x, a.y - b.y, a.z - b.z};
}

point3d_t point3d_scale(point3d_t p, double s) {
    return (point3d_t){p.x * s, p.y * s, p.z * s};
}

double point3d_length(point3d_t p) {
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

double point3d_distance(point3d_t a, point3d_t b) {
    return point3d_length(point3d_sub(a, b));
}

point3d_t point3d_normalize(point3d_t p) {
    double len = point3d_length(p);
    if (len < 1e-10) return (point3d_t){1, 0, 0};
    return point3d_scale(p, 1.0 / len);
}

double point3d_dot(point3d_t a, point3d_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

point3d_t point3d_cross(point3d_t a, point3d_t b) {
    return (point3d_t){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}
