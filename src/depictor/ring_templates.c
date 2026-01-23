/**
 * @file ring_templates.c
 * @brief Ring template matching for 2D molecular layout
 *
 * Provides pre-defined templates for common ring systems and
 * matching/transformation algorithms.
 */

#include "cchem/depictor/templates.h"
#include "cchem/depictor/layout2d.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============== Pre-computed Ring Coordinates ============== */

/*
 * Pre-computed coordinates for regular polygons centered at origin with circumradius 1.0
 * Starting angle is PI/2 + PI/n (places first vertex at top)
 *
 * For n-gon: vertex i at angle = PI/2 + PI/n - i * 2*PI/n
 *
 * 3-gon (triangle):   cos/sin values for angles 90+60, 90+60-120, 90+60-240 degrees
 * 4-gon (square):     cos/sin values for angles 135, 45, -45, -135 degrees
 * 5-gon (pentagon):   cos/sin values for angles 126, 54, -18, -90, -162 degrees
 * 6-gon (hexagon):    cos/sin values for angles 120, 60, 0, -60, -120, -180 degrees
 * 7-gon (heptagon):   cos/sin values for angles ~115.7, ~64.3, ~12.9, -38.6, -90, -141.4, -192.9
 */

/* Pre-computed trigonometric values for common ring sizes */

/* 3-membered ring (triangle) */
#define TRI_R 1.0
#define TRI_A0_X (-0.5)          /* cos(150 deg) */
#define TRI_A0_Y (0.8660254038)  /* sin(150 deg) */
#define TRI_A1_X (0.5)           /* cos(30 deg) */
#define TRI_A1_Y (0.8660254038)  /* sin(30 deg) */
#define TRI_A2_X (0.0)           /* cos(-90 deg) */
#define TRI_A2_Y (-1.0)          /* sin(-90 deg) */

/* 4-membered ring (square) */
#define SQ_R 1.0
#define SQ_A0_X (-0.7071067812) /* cos(135 deg) */
#define SQ_A0_Y (0.7071067812)  /* sin(135 deg) */
#define SQ_A1_X (0.7071067812)  /* cos(45 deg) */
#define SQ_A1_Y (0.7071067812)  /* sin(45 deg) */
#define SQ_A2_X (0.7071067812)  /* cos(-45 deg) */
#define SQ_A2_Y (-0.7071067812) /* sin(-45 deg) */
#define SQ_A3_X (-0.7071067812) /* cos(-135 deg) */
#define SQ_A3_Y (-0.7071067812) /* sin(-135 deg) */

/* 5-membered ring (pentagon) */
#define PENT_R 1.0
#define PENT_A0_X (-0.5877852523) /* cos(126 deg) */
#define PENT_A0_Y (0.8090169944)  /* sin(126 deg) */
#define PENT_A1_X (0.5877852523)  /* cos(54 deg) */
#define PENT_A1_Y (0.8090169944)  /* sin(54 deg) */
#define PENT_A2_X (0.9510565163)  /* cos(-18 deg) */
#define PENT_A2_Y (-0.3090169944) /* sin(-18 deg) */
#define PENT_A3_X (0.0)           /* cos(-90 deg) */
#define PENT_A3_Y (-1.0)          /* sin(-90 deg) */
#define PENT_A4_X (-0.9510565163) /* cos(-162 deg) */
#define PENT_A4_Y (-0.3090169944) /* sin(-162 deg) */

/* 6-membered ring (hexagon) - flat-top orientation */
#define HEX_R 1.0
#define HEX_A0_X (-0.5)          /* cos(120 deg) */
#define HEX_A0_Y (0.8660254038)  /* sin(120 deg) */
#define HEX_A1_X (0.5)           /* cos(60 deg) */
#define HEX_A1_Y (0.8660254038)  /* sin(60 deg) */
#define HEX_A2_X (1.0)           /* cos(0 deg) */
#define HEX_A2_Y (0.0)           /* sin(0 deg) */
#define HEX_A3_X (0.5)           /* cos(-60 deg) */
#define HEX_A3_Y (-0.8660254038) /* sin(-60 deg) */
#define HEX_A4_X (-0.5)          /* cos(-120 deg) */
#define HEX_A4_Y (-0.8660254038) /* sin(-120 deg) */
#define HEX_A5_X (-1.0)          /* cos(180 deg) */
#define HEX_A5_Y (0.0)           /* sin(180 deg) */

/* 7-membered ring (heptagon) */
#define HEPT_R 1.0
#define HEPT_A0_X (-0.4338837391) /* cos(115.71 deg) */
#define HEPT_A0_Y (0.9009688679)  /* sin(115.71 deg) */
#define HEPT_A1_X (0.4338837391)  /* cos(64.29 deg) */
#define HEPT_A1_Y (0.9009688679)  /* sin(64.29 deg) */
#define HEPT_A2_X (0.9749279122)  /* cos(12.86 deg) */
#define HEPT_A2_Y (0.2225209340)  /* sin(12.86 deg) */
#define HEPT_A3_X (0.7818314825)  /* cos(-38.57 deg) */
#define HEPT_A3_Y (-0.6234898019) /* sin(-38.57 deg) */
#define HEPT_A4_X (0.0)           /* cos(-90 deg) */
#define HEPT_A4_Y (-1.0)          /* sin(-90 deg) */
#define HEPT_A5_X (-0.7818314825) /* cos(-141.43 deg) */
#define HEPT_A5_Y (-0.6234898019) /* sin(-141.43 deg) */
#define HEPT_A6_X (-0.9749279122) /* cos(-192.86 deg) = cos(167.14 deg) */
#define HEPT_A6_Y (0.2225209340)  /* sin(-192.86 deg) */

/* ============== Template Definitions ============== */

static const ring_template_t BUILTIN_TEMPLATES[] = {
    /* Cyclopropane (3-membered ring) */
    {
        .name = "cyclopropane",
        .num_atoms = 3,
        .ring_sizes = {3, 0, 0, 0, 0, 0, 0, 0},
        .num_rings = 1,
        .coords = {
            {TRI_A0_X, TRI_A0_Y},
            {TRI_A1_X, TRI_A1_Y},
            {TRI_A2_X, TRI_A2_Y},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        },
        .adjacency = {
            {1, 2, -1, -1}, {0, 2, -1, -1}, {1, 0, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}
        }
    },

    /* Cyclobutane (4-membered ring) */
    {
        .name = "cyclobutane",
        .num_atoms = 4,
        .ring_sizes = {4, 0, 0, 0, 0, 0, 0, 0},
        .num_rings = 1,
        .coords = {
            {SQ_A0_X, SQ_A0_Y},
            {SQ_A1_X, SQ_A1_Y},
            {SQ_A2_X, SQ_A2_Y},
            {SQ_A3_X, SQ_A3_Y},
            {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        },
        .adjacency = {
            {1, 3, -1, -1}, {0, 2, -1, -1}, {1, 3, -1, -1}, {2, 0, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}
        }
    },

    /* Cyclopentane (5-membered ring) */
    {
        .name = "cyclopentane",
        .num_atoms = 5,
        .ring_sizes = {5, 0, 0, 0, 0, 0, 0, 0},
        .num_rings = 1,
        .coords = {
            {PENT_A0_X, PENT_A0_Y},
            {PENT_A1_X, PENT_A1_Y},
            {PENT_A2_X, PENT_A2_Y},
            {PENT_A3_X, PENT_A3_Y},
            {PENT_A4_X, PENT_A4_Y},
            {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        },
        .adjacency = {
            {1, 4, -1, -1}, {0, 2, -1, -1}, {1, 3, -1, -1}, {2, 4, -1, -1}, {3, 0, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}
        }
    },

    /* Cyclohexane / Benzene (6-membered ring) */
    {
        .name = "cyclohexane",
        .num_atoms = 6,
        .ring_sizes = {6, 0, 0, 0, 0, 0, 0, 0},
        .num_rings = 1,
        .coords = {
            {HEX_A0_X, HEX_A0_Y},
            {HEX_A1_X, HEX_A1_Y},
            {HEX_A2_X, HEX_A2_Y},
            {HEX_A3_X, HEX_A3_Y},
            {HEX_A4_X, HEX_A4_Y},
            {HEX_A5_X, HEX_A5_Y},
            {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        },
        .adjacency = {
            {1, 5, -1, -1}, {0, 2, -1, -1}, {1, 3, -1, -1}, {2, 4, -1, -1}, {3, 5, -1, -1}, {4, 0, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}
        }
    },

    /* Cycloheptane (7-membered ring) */
    {
        .name = "cycloheptane",
        .num_atoms = 7,
        .ring_sizes = {7, 0, 0, 0, 0, 0, 0, 0},
        .num_rings = 1,
        .coords = {
            {HEPT_A0_X, HEPT_A0_Y},
            {HEPT_A1_X, HEPT_A1_Y},
            {HEPT_A2_X, HEPT_A2_Y},
            {HEPT_A3_X, HEPT_A3_Y},
            {HEPT_A4_X, HEPT_A4_Y},
            {HEPT_A5_X, HEPT_A5_Y},
            {HEPT_A6_X, HEPT_A6_Y},
            {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        },
        .adjacency = {
            {1, 6, -1, -1}, {0, 2, -1, -1}, {1, 3, -1, -1}, {2, 4, -1, -1},
            {3, 5, -1, -1}, {4, 6, -1, -1}, {5, 0, -1, -1},
            {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}
        }
    },

    /* Naphthalene (two fused 6-membered rings) */
    {
        .name = "naphthalene",
        .num_atoms = 10,
        .ring_sizes = {6, 6, 0, 0, 0, 0, 0, 0},
        .num_rings = 2,
        .coords = {
            /* First ring (left hexagon) */
            {-0.8660254038, 0.5},    /* 0 - top-left */
            {0.0, 1.0},              /* 1 - top-center */
            {0.8660254038, 0.5},     /* 2 - shared top */
            {0.8660254038, -0.5},    /* 3 - shared bottom */
            {0.0, -1.0},             /* 4 - bottom-center */
            {-0.8660254038, -0.5},   /* 5 - bottom-left */
            /* Second ring (right hexagon, unique atoms) */
            {1.7320508076, 1.0},     /* 6 - top-right */
            {2.5980762114, 0.5},     /* 7 - far-right top */
            {2.5980762114, -0.5},    /* 8 - far-right bottom */
            {1.7320508076, -1.0},    /* 9 - bottom-right */
            /* Padding */
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        },
        .adjacency = {
            {1, 5, -1, -1}, {0, 2, -1, -1}, {1, 3, 6, -1}, {2, 4, 9, -1},
            {3, 5, -1, -1}, {4, 0, -1, -1}, {2, 7, -1, -1}, {6, 8, -1, -1},
            {7, 9, -1, -1}, {8, 3, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}
        }
    },

    /* Indole / Benzimidazole (6+5 fused) */
    {
        .name = "indole",
        .num_atoms = 9,
        .ring_sizes = {6, 5, 0, 0, 0, 0, 0, 0},
        .num_rings = 2,
        .coords = {
            /* 6-membered ring */
            {-0.8660254038, 0.5},
            {0.0, 1.0},
            {0.8660254038, 0.5},     /* shared */
            {0.8660254038, -0.5},    /* shared */
            {0.0, -1.0},
            {-0.8660254038, -0.5},
            /* 5-membered ring unique atoms */
            {1.5388417686, 0.9048270525},
            {1.9021130326, 0.0},
            {1.5388417686, -0.9048270525},
            /* Padding */
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
            {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
        },
        .adjacency = {
            {1, 5, -1, -1}, {0, 2, -1, -1}, {1, 3, 6, -1}, {2, 4, 8, -1},
            {3, 5, -1, -1}, {4, 0, -1, -1}, {2, 7, -1, -1}, {6, 8, -1, -1}, {7, 3, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1},
            {-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}
        }
    }
};

static const int NUM_BUILTIN_TEMPLATES = sizeof(BUILTIN_TEMPLATES) / sizeof(BUILTIN_TEMPLATES[0]);

/* ============== Template Access ============== */

const ring_template_t* template_get_builtin(int* count) {
    if (count) *count = NUM_BUILTIN_TEMPLATES;
    return BUILTIN_TEMPLATES;
}

/* ============== Template Matching ============== */

/* Sort ring sizes for comparison */
static void sort_ring_sizes(int* sizes, int count) {
    for (int i = 0; i < count - 1; i++) {
        for (int j = i + 1; j < count; j++) {
            if (sizes[i] > sizes[j]) {
                int tmp = sizes[i];
                sizes[i] = sizes[j];
                sizes[j] = tmp;
            }
        }
    }
}

bool template_ring_signature_matches(const ring_system_t* system,
                                     const ring_template_t* templ,
                                     const molecule_t* mol) {
    if (system->num_rings != templ->num_rings) return false;

    /* Get ring sizes from system */
    int sys_sizes[MAX_TEMPLATE_RINGS];
    for (int i = 0; i < system->num_rings; i++) {
        sys_sizes[i] = mol->rings[system->ring_indices[i]].size;
    }
    sort_ring_sizes(sys_sizes, system->num_rings);

    /* Get ring sizes from template */
    int tmpl_sizes[MAX_TEMPLATE_RINGS];
    int tmpl_count = 0;
    for (int i = 0; i < MAX_TEMPLATE_RINGS && templ->ring_sizes[i] > 0; i++) {
        tmpl_sizes[tmpl_count++] = templ->ring_sizes[i];
    }
    sort_ring_sizes(tmpl_sizes, tmpl_count);

    if (system->num_rings != tmpl_count) return false;

    /* Compare sorted sizes */
    for (int i = 0; i < system->num_rings; i++) {
        if (sys_sizes[i] != tmpl_sizes[i]) return false;
    }

    return true;
}

const ring_template_t* template_find_match(const layout_context_t* ctx,
                                           const ring_system_t* system) {
    if (!ctx || !system) return NULL;

    const molecule_t* mol = ctx->mol;

    /* Try each built-in template */
    for (int i = 0; i < NUM_BUILTIN_TEMPLATES; i++) {
        if (template_ring_signature_matches(system, &BUILTIN_TEMPLATES[i], mol)) {
            /* Additional validation: atom count should match */
            if (system->num_atoms == BUILTIN_TEMPLATES[i].num_atoms) {
                return &BUILTIN_TEMPLATES[i];
            }
        }
    }

    return NULL;
}

/* ============== Subgraph Isomorphism for Atom Mapping ============== */

/* Build adjacency list for ring system atoms based on molecular connectivity */
static void build_system_adjacency(const layout_context_t* ctx,
                                   const ring_system_t* system,
                                   int adjacency[][4]) {
    const molecule_t* mol = ctx->mol;

    /* Initialize adjacency to -1 */
    for (int i = 0; i < system->num_atoms; i++) {
        for (int j = 0; j < 4; j++) {
            adjacency[i][j] = -1;
        }
    }

    /* Build mapping from molecule atom index to system atom index */
    int mol_to_sys[MAX_TEMPLATE_ATOMS];
    for (int i = 0; i < MAX_TEMPLATE_ATOMS; i++) mol_to_sys[i] = -1;
    for (int i = 0; i < system->num_atoms; i++) {
        if (system->all_atoms[i] < MAX_TEMPLATE_ATOMS) {
            mol_to_sys[system->all_atoms[i]] = i;
        }
    }

    /* For each atom in system, find its neighbors that are also in system */
    for (int i = 0; i < system->num_atoms; i++) {
        int mol_atom = system->all_atoms[i];
        const atom_t* atom = &mol->atoms[mol_atom];
        int adj_idx = 0;

        for (int n = 0; n < atom->num_neighbors && adj_idx < 4; n++) {
            int nb = atom->neighbors[n];
            if (nb < MAX_TEMPLATE_ATOMS && mol_to_sys[nb] >= 0) {
                adjacency[i][adj_idx++] = mol_to_sys[nb];
            }
        }
    }
}

/* Count neighbors in adjacency list */
static int count_adjacency(const int adj[4]) {
    int count = 0;
    for (int i = 0; i < 4 && adj[i] >= 0; i++) count++;
    return count;
}

/* Check if sys_atom's connectivity matches tmpl_atom's connectivity given current mapping */
static bool check_connectivity(int sys_atom, int tmpl_atom,
                               const int sys_adj[][4], const ring_template_t* templ,
                               const int* mapping, const int* reverse_mapping,
                               int num_sys_atoms) {
    /* Count neighbors */
    int sys_degree = count_adjacency(sys_adj[sys_atom]);
    int tmpl_degree = count_adjacency(templ->adjacency[tmpl_atom]);

    if (sys_degree != tmpl_degree) return false;

    /* Check that mapped neighbors match */
    for (int i = 0; i < 4 && sys_adj[sys_atom][i] >= 0; i++) {
        int sys_nb = sys_adj[sys_atom][i];
        if (sys_nb >= num_sys_atoms) continue;

        int mapped_tmpl_nb = mapping[sys_nb];
        if (mapped_tmpl_nb < 0) continue; /* Neighbor not yet mapped */

        /* Check if tmpl_atom has mapped_tmpl_nb as neighbor */
        bool found = false;
        for (int j = 0; j < 4 && templ->adjacency[tmpl_atom][j] >= 0; j++) {
            if (templ->adjacency[tmpl_atom][j] == mapped_tmpl_nb) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }

    /* Check reverse: template neighbors that are mapped should map to system neighbors */
    for (int i = 0; i < 4 && templ->adjacency[tmpl_atom][i] >= 0; i++) {
        int tmpl_nb = templ->adjacency[tmpl_atom][i];
        int mapped_sys_nb = reverse_mapping[tmpl_nb];
        if (mapped_sys_nb < 0) continue; /* Template neighbor not yet mapped */

        /* Check if sys_atom has mapped_sys_nb as neighbor */
        bool found = false;
        for (int j = 0; j < 4 && sys_adj[sys_atom][j] >= 0; j++) {
            if (sys_adj[sys_atom][j] == mapped_sys_nb) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }

    return true;
}

/* Recursive backtracking for subgraph isomorphism */
static bool isomorphism_backtrack(int sys_idx, int num_sys_atoms, int num_tmpl_atoms,
                                  const int sys_adj[][4], const ring_template_t* templ,
                                  int* mapping, int* reverse_mapping, bool* tmpl_used) {
    if (sys_idx >= num_sys_atoms) {
        return true; /* All system atoms mapped successfully */
    }

    int sys_degree = count_adjacency(sys_adj[sys_idx]);

    /* Try each unmapped template atom */
    for (int t = 0; t < num_tmpl_atoms; t++) {
        if (tmpl_used[t]) continue;

        /* Quick degree check */
        int tmpl_degree = count_adjacency(templ->adjacency[t]);
        if (sys_degree != tmpl_degree) continue;

        /* Check connectivity constraints */
        if (!check_connectivity(sys_idx, t, sys_adj, templ, mapping, reverse_mapping, num_sys_atoms)) {
            continue;
        }

        /* Try this mapping */
        mapping[sys_idx] = t;
        reverse_mapping[t] = sys_idx;
        tmpl_used[t] = true;

        if (isomorphism_backtrack(sys_idx + 1, num_sys_atoms, num_tmpl_atoms,
                                  sys_adj, templ, mapping, reverse_mapping, tmpl_used)) {
            return true;
        }

        /* Backtrack */
        mapping[sys_idx] = -1;
        reverse_mapping[t] = -1;
        tmpl_used[t] = false;
    }

    return false;
}

/* ============== Template Application ============== */

bool template_build_atom_mapping(const layout_context_t* ctx,
                                 const ring_system_t* system,
                                 const ring_template_t* templ,
                                 int* mapping) {
    /* Initialize mapping to unmapped */
    for (int i = 0; i < system->num_atoms; i++) {
        mapping[i] = -1;
    }

    /* For single rings, use ring order matching */
    if (system->num_rings == 1) {
        /* Get ring atoms in bond-connected order */
        const molecule_t* mol = ctx->mol;
        int ring_idx = system->ring_indices[0];
        const ring_t* ring = &mol->rings[ring_idx];

        int ordered[MAX_TEMPLATE_ATOMS];
        int num_ordered = 0;
        layout_get_ring_order(ring, mol, ring->atoms[0], ordered, &num_ordered);

        if (num_ordered == ring->size && num_ordered == templ->num_atoms) {
            /* Map ordered ring atoms to template atoms sequentially */
            for (int i = 0; i < num_ordered; i++) {
                /* Find system index for this ring atom */
                for (int j = 0; j < system->num_atoms; j++) {
                    if (system->all_atoms[j] == ordered[i]) {
                        mapping[j] = i;
                        break;
                    }
                }
            }
            return true;
        }

        /* Fallback: simple sequential */
        for (int i = 0; i < system->num_atoms && i < templ->num_atoms; i++) {
            mapping[i] = i;
        }
        return true;
    }

    /* For fused systems, use subgraph isomorphism */
    int sys_adj[MAX_TEMPLATE_ATOMS][4];
    build_system_adjacency(ctx, system, sys_adj);

    int reverse_mapping[MAX_TEMPLATE_ATOMS];
    bool tmpl_used[MAX_TEMPLATE_ATOMS];
    for (int i = 0; i < MAX_TEMPLATE_ATOMS; i++) {
        reverse_mapping[i] = -1;
        tmpl_used[i] = false;
    }

    if (isomorphism_backtrack(0, system->num_atoms, templ->num_atoms,
                              sys_adj, templ, mapping, reverse_mapping, tmpl_used)) {
        return true;
    }

    /* Fallback: simple sequential mapping if isomorphism fails */
    for (int i = 0; i < system->num_atoms && i < templ->num_atoms; i++) {
        mapping[i] = i;
    }

    return true;
}

/* Calculate optimal rotation angle to align source points with target points */
static double calculate_optimal_rotation(const point2d_t* source, const point2d_t* target,
                                         int num_points, point2d_t src_centroid,
                                         point2d_t tgt_centroid) {
    /*
     * For 2D rotation, we minimize sum of squared distances:
     * min_theta sum_i |R(theta) * (s_i - s_c) - (t_i - t_c)|^2
     *
     * The optimal angle is: theta = atan2(sum(cross), sum(dot))
     * where cross and dot are computed from centered points
     */
    double sum_cross = 0.0;
    double sum_dot = 0.0;

    for (int i = 0; i < num_points; i++) {
        /* Center the points */
        double sx = source[i].x - src_centroid.x;
        double sy = source[i].y - src_centroid.y;
        double tx = target[i].x - tgt_centroid.x;
        double ty = target[i].y - tgt_centroid.y;

        /* Accumulate cross and dot products */
        sum_cross += sx * ty - sy * tx;
        sum_dot += sx * tx + sy * ty;
    }

    if (fabs(sum_dot) < 1e-10 && fabs(sum_cross) < 1e-10) {
        return 0.0;
    }

    return atan2(sum_cross, sum_dot);
}

bool template_find_transformation(const layout_context_t* ctx,
                                  const ring_system_t* system,
                                  const ring_template_t* templ,
                                  point2d_t* out_center, double* out_rotation) {
    if (!ctx || !system || !templ) return false;

    /* Default: center at origin, no rotation */
    out_center->x = 0.0;
    out_center->y = 0.0;
    *out_rotation = 0.0;

    /* Build atom mapping first */
    int mapping[MAX_TEMPLATE_ATOMS];
    if (!template_build_atom_mapping(ctx, system, templ, mapping)) {
        return false;
    }

    /* Collect all placed atoms in this system */
    point2d_t placed_mol_pos[MAX_TEMPLATE_ATOMS];
    point2d_t placed_tmpl_pos[MAX_TEMPLATE_ATOMS];
    int num_placed = 0;

    for (int i = 0; i < system->num_atoms; i++) {
        int mol_atom = system->all_atoms[i];
        int tmpl_atom = mapping[i];

        if (ctx->placed[mol_atom] && tmpl_atom >= 0 && tmpl_atom < templ->num_atoms) {
            placed_mol_pos[num_placed] = ctx->coords->coords_2d[mol_atom];
            placed_tmpl_pos[num_placed].x = templ->coords[tmpl_atom].x * ctx->bond_length;
            placed_tmpl_pos[num_placed].y = templ->coords[tmpl_atom].y * ctx->bond_length;
            num_placed++;
        }
    }

    if (num_placed == 0) {
        /* No placed atoms - use default transformation */
        return true;
    }

    if (num_placed == 1) {
        /* Single placed atom - align template to that position, no rotation */
        out_center->x = placed_mol_pos[0].x - placed_tmpl_pos[0].x;
        out_center->y = placed_mol_pos[0].y - placed_tmpl_pos[0].y;
        *out_rotation = 0.0;
        return true;
    }

    /* Multiple placed atoms - calculate optimal rotation */

    /* Calculate centroids */
    point2d_t mol_centroid = {0, 0};
    point2d_t tmpl_centroid = {0, 0};
    for (int i = 0; i < num_placed; i++) {
        mol_centroid.x += placed_mol_pos[i].x;
        mol_centroid.y += placed_mol_pos[i].y;
        tmpl_centroid.x += placed_tmpl_pos[i].x;
        tmpl_centroid.y += placed_tmpl_pos[i].y;
    }
    mol_centroid.x /= num_placed;
    mol_centroid.y /= num_placed;
    tmpl_centroid.x /= num_placed;
    tmpl_centroid.y /= num_placed;

    /* Calculate optimal rotation to align template to molecule */
    *out_rotation = calculate_optimal_rotation(placed_tmpl_pos, placed_mol_pos,
                                               num_placed, tmpl_centroid, mol_centroid);

    /* Calculate center: after rotation, the rotated template centroid should align with mol centroid */
    double cos_r = cos(*out_rotation);
    double sin_r = sin(*out_rotation);

    /* Rotated template centroid */
    double rotated_tmpl_cx = tmpl_centroid.x * cos_r - tmpl_centroid.y * sin_r;
    double rotated_tmpl_cy = tmpl_centroid.x * sin_r + tmpl_centroid.y * cos_r;

    out_center->x = mol_centroid.x - rotated_tmpl_cx;
    out_center->y = mol_centroid.y - rotated_tmpl_cy;

    return true;
}

bool template_apply(layout_context_t* ctx, const ring_system_t* system,
                    const ring_template_t* templ, point2d_t center,
                    double scale, double rotation) {
    if (!ctx || !system || !templ) return false;

    /* Build atom mapping */
    int mapping[MAX_TEMPLATE_ATOMS];
    if (!template_build_atom_mapping(ctx, system, templ, mapping)) {
        return false;
    }

    /* Apply transformation to each atom */
    double cos_r = cos(rotation);
    double sin_r = sin(rotation);

    for (int i = 0; i < system->num_atoms; i++) {
        int mol_atom = system->all_atoms[i];
        int tmpl_atom = mapping[i];

        if (tmpl_atom < 0 || tmpl_atom >= templ->num_atoms) continue;

        /* Get template coordinate */
        double tx = templ->coords[tmpl_atom].x * scale;
        double ty = templ->coords[tmpl_atom].y * scale;

        /* Rotate */
        double rx = tx * cos_r - ty * sin_r;
        double ry = tx * sin_r + ty * cos_r;

        /* Translate */
        ctx->coords->coords_2d[mol_atom].x = center.x + rx;
        ctx->coords->coords_2d[mol_atom].y = center.y + ry;
        ctx->placed[mol_atom] = true;
    }

    return true;
}
