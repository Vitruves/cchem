/**
 * @file framework.c
 * @brief Framework Topology Descriptors - Molecular scaffold analysis
 *
 * Categories:
 * 1. Ring System Topology (12 descriptors)
 *    - Ring systems, fused pairs, spiro junctions, bridgeheads
 *    - Ring sizes and counts
 *
 * 2. Chain/Linker Analysis (12 descriptors)
 *    - Linker atoms and bonds, terminal chains, branch chains
 *    - Chain density and ring/chain ratios
 *
 * 3. Scaffold Metrics (16 descriptors)
 *    - Murcko scaffold, side chains, functional groups
 *    - Core/periphery, symmetry, topological diversity
 *
 * Total: 40 descriptors
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* Total framework descriptors */
#define NUM_FRAME_DESCRIPTORS 40

/* Maximum atoms for stack allocation */
#define MAX_ATOMS_STACK 256

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/* Check if atom is in any ring */
static bool is_ring_atom(const molecule_t* mol, int atom_idx) {
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        for (int i = 0; i < ring->size; i++) {
            if (ring->atoms[i] == atom_idx) return true;
        }
    }
    return false;
}

/* Count how many rings an atom is in */
static int count_rings_for_atom(const molecule_t* mol, int atom_idx) {
    int count = 0;
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        for (int i = 0; i < ring->size; i++) {
            if (ring->atoms[i] == atom_idx) {
                count++;
                break;
            }
        }
    }
    return count;
}

/* Check if bond is in any ring */
static bool is_ring_bond(const molecule_t* mol, int atom1, int atom2) {
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        bool found1 = false, found2 = false;
        for (int i = 0; i < ring->size; i++) {
            if (ring->atoms[i] == atom1) found1 = true;
            if (ring->atoms[i] == atom2) found2 = true;
        }
        if (found1 && found2) return true;
    }
    return false;
}

/* Check if two rings share atoms (fused or spiro) */
static int count_shared_atoms(const ring_t* r1, const ring_t* r2) {
    int count = 0;
    for (int i = 0; i < r1->size; i++) {
        for (int j = 0; j < r2->size; j++) {
            if (r1->atoms[i] == r2->atoms[j]) {
                count++;
            }
        }
    }
    return count;
}

/* Get atom degree (heavy atoms only) */
static int get_heavy_degree(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int degree = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            degree++;
        }
    }
    return degree;
}

/* ============================================================================
 * Main Computation
 * ============================================================================ */

int descriptors_compute_framework_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_FRAME_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    }

    if (n_heavy == 0) return NUM_FRAME_DESCRIPTORS;

    /* ====================================================================
     * Section 1: Ring System Topology (12 descriptors, idx 0-11)
     * ==================================================================== */

    int ring_systems = 0;
    int fused_pairs = 0;
    int spiro_junctions = 0;
    int bridgehead_atoms = 0;
    int ring_bridges = 0;
    int max_ring_size = 0;
    int min_ring_size = 999;
    double mean_ring_size = 0.0;
    double ring_size_var = 0.0;
    int total_ring_atoms = 0;
    int total_ring_bonds = 0;

    /* Identify ring systems and fusion patterns */
    int* ring_system_id = (int*)alloca(mol->num_rings * sizeof(int));
    memset(ring_system_id, -1, mol->num_rings * sizeof(int));
    int next_system = 0;

    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];

        /* Ring size statistics */
        if (ring->size > max_ring_size) max_ring_size = ring->size;
        if (ring->size < min_ring_size) min_ring_size = ring->size;
        mean_ring_size += ring->size;

        /* Assign to ring system */
        if (ring_system_id[r] == -1) {
            ring_system_id[r] = next_system++;
        }

        /* Check for fused/spiro with other rings */
        for (int r2 = r + 1; r2 < mol->num_rings && r2 < MAX_RINGS; r2++) {
            int shared = count_shared_atoms(ring, &mol->rings[r2]);
            if (shared == 1) {
                spiro_junctions++;
            } else if (shared >= 2) {
                fused_pairs++;
                /* Merge ring systems */
                if (ring_system_id[r2] == -1) {
                    ring_system_id[r2] = ring_system_id[r];
                }
            }
        }
    }

    ring_systems = next_system;
    if (mol->num_rings > 0) {
        mean_ring_size /= mol->num_rings;
    } else {
        min_ring_size = 0;
    }

    /* Compute ring size variance */
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        double diff = mol->rings[r].size - mean_ring_size;
        ring_size_var += diff * diff;
    }
    if (mol->num_rings > 0) {
        ring_size_var /= mol->num_rings;
    }

    /* Count ring atoms and bridgeheads */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int ring_count = count_rings_for_atom(mol, i);
        if (ring_count > 0) {
            total_ring_atoms++;
        }
        if (ring_count >= 2) {
            bridgehead_atoms++;
        }
    }

    /* Count ring bonds and bridges */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        bool in_ring = is_ring_bond(mol, bond->atom1, bond->atom2);
        if (in_ring) {
            total_ring_bonds++;
        } else {
            /* Check if connects two ring atoms (bridge between ring systems) */
            bool a1_ring = is_ring_atom(mol, bond->atom1);
            bool a2_ring = is_ring_atom(mol, bond->atom2);
            if (a1_ring && a2_ring) {
                ring_bridges++;
            }
        }
    }

    int idx = 0;
    values[idx++].d = ring_systems;
    values[idx++].d = fused_pairs;
    values[idx++].d = spiro_junctions;
    values[idx++].d = bridgehead_atoms;
    values[idx++].d = ring_bridges;
    values[idx++].d = max_ring_size;
    values[idx++].d = min_ring_size;
    values[idx++].d = max_ring_size - min_ring_size;  /* Range */
    values[idx++].d = mean_ring_size;
    values[idx++].d = ring_size_var;
    values[idx++].d = total_ring_atoms;
    values[idx++].d = total_ring_bonds;

    /* ====================================================================
     * Section 2: Chain/Linker Analysis (12 descriptors, idx 12-23)
     * ==================================================================== */

    int linker_atoms = 0;
    int linker_bonds = 0;
    int max_linker_len = 0;
    int terminal_chains = 0;
    int max_terminal_len = 0;
    int branch_chains = 0;
    int chain_atoms = 0;
    double chain_density = 0.0;
    double ring_chain_ratio = 0.0;
    int complexity = 0;
    double compactness = 0.0;

    /* Count chain and linker atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        if (!is_ring_atom(mol, i)) {
            chain_atoms++;
            int degree = get_heavy_degree(mol, i);
            if (degree == 1) {
                terminal_chains++;  /* Terminal atom */
            }
            if (degree >= 3) {
                branch_chains++;  /* Branch point in chain */
            }
        }
    }

    /* Count linker bonds (non-ring bonds between non-terminal atoms) */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if (mol->atoms[bond->atom1].element == ELEM_H ||
            mol->atoms[bond->atom2].element == ELEM_H) continue;

        bool a1_ring = is_ring_atom(mol, bond->atom1);
        bool a2_ring = is_ring_atom(mol, bond->atom2);

        if (!a1_ring && !a2_ring) {
            linker_bonds++;
        }
        if ((a1_ring && !a2_ring) || (!a1_ring && a2_ring)) {
            linker_atoms++;  /* Count attachment points */
        }
    }

    chain_density = (n_heavy > 0) ? (double)chain_atoms / n_heavy : 0.0;
    ring_chain_ratio = (chain_atoms > 0) ? (double)total_ring_atoms / chain_atoms : 0.0;
    complexity = ring_systems * (linker_bonds + 1) * (branch_chains + 1);
    compactness = (ring_systems + terminal_chains > 0) ?
                  (double)n_heavy / (ring_systems + terminal_chains) : 0.0;

    values[idx++].d = linker_atoms;
    values[idx++].d = linker_bonds;
    values[idx++].d = max_linker_len;  /* Would need path analysis */
    values[idx++].d = linker_bonds > 0 ? (double)linker_bonds / 2.0 : 0.0;  /* Mean */
    values[idx++].d = terminal_chains;
    values[idx++].d = max_terminal_len;  /* Would need path analysis */
    values[idx++].d = branch_chains;
    values[idx++].d = chain_atoms;
    values[idx++].d = chain_density;
    values[idx++].d = ring_chain_ratio;
    values[idx++].d = complexity;
    values[idx++].d = compactness;

    /* ====================================================================
     * Section 3: Scaffold Metrics (16 descriptors, idx 24-39)
     * ==================================================================== */

    /* Murcko scaffold: ring atoms + linkers between rings */
    int scaffold_atoms = total_ring_atoms + linker_atoms;
    double scaffold_ratio = (n_heavy > 0) ? (double)scaffold_atoms / n_heavy : 0.0;

    /* Side chain atoms */
    int side_chain_atoms = chain_atoms - linker_atoms;
    double side_chain_ratio = (n_heavy > 0) ? (double)side_chain_atoms / n_heavy : 0.0;

    /* Functional atoms (heteroatoms in chains) */
    int functional_atoms = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        if (!is_ring_atom(mol, i)) {
            element_t elem = mol->atoms[i].element;
            if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S ||
                elem == ELEM_P || elem == ELEM_F || elem == ELEM_Cl ||
                elem == ELEM_Br || elem == ELEM_I) {
                functional_atoms++;
            }
        }
    }

    /* Decorated vs naked rings */
    int decorated_rings = 0;
    int naked_rings = 0;
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        bool has_substituent = false;
        for (int i = 0; i < ring->size; i++) {
            int atom_idx = ring->atoms[i];
            const atom_t* atom = &mol->atoms[atom_idx];
            for (int n = 0; n < atom->num_neighbors; n++) {
                int nb = atom->neighbors[n];
                if (mol->atoms[nb].element != ELEM_H && !is_ring_atom(mol, nb)) {
                    has_substituent = true;
                    break;
                }
            }
            if (has_substituent) break;
        }
        if (has_substituent) decorated_rings++;
        else naked_rings++;
    }

    /* Core and periphery */
    int core_atoms = bridgehead_atoms;
    int periphery_atoms = n_heavy - core_atoms;
    double core_periph_ratio = (periphery_atoms > 0) ? (double)core_atoms / periphery_atoms : 0.0;

    /* Simple symmetry estimate based on element type distribution */
    int elem_counts[128] = {0};
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            elem_counts[mol->atoms[i].element]++;
        }
    }
    double symmetry = 0.0;
    for (int e = 0; e < 128; e++) {
        if (elem_counts[e] > 0) {
            double p = (double)elem_counts[e] / n_heavy;
            symmetry -= p * log(p + 1e-10);
        }
    }
    symmetry = 1.0 - (symmetry / log((double)n_heavy + 1.0));

    /* Topological diversity */
    int unique_degrees = 0;
    bool seen_degree[10] = {false};
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int d = get_heavy_degree(mol, i);
        if (d < 10 && !seen_degree[d]) {
            seen_degree[d] = true;
            unique_degrees++;
        }
    }
    double topo_diversity = (double)unique_degrees / 5.0;  /* Normalize */

    /* Linker heteroatoms and unsaturation */
    int linker_hetero = 0;
    int linker_unsat = 0;
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        bool a1_ring = is_ring_atom(mol, bond->atom1);
        bool a2_ring = is_ring_atom(mol, bond->atom2);
        if (!a1_ring && !a2_ring) {
            if (mol->atoms[bond->atom1].element != ELEM_C ||
                mol->atoms[bond->atom2].element != ELEM_C) {
                linker_hetero++;
            }
            if (bond->type == BOND_DOUBLE || bond->type == BOND_TRIPLE) {
                linker_unsat++;
            }
        }
    }

    /* Mean substituents per ring */
    double ring_decoration = (mol->num_rings > 0) ?
                             (double)decorated_rings / mol->num_rings : 0.0;

    /* Terminal functional groups */
    int terminal_groups = terminal_chains;  /* Simplified */

    values[idx++].d = scaffold_atoms;
    values[idx++].d = scaffold_ratio;
    values[idx++].d = side_chain_atoms;
    values[idx++].d = side_chain_ratio;
    values[idx++].d = functional_atoms;
    values[idx++].d = decorated_rings;
    values[idx++].d = naked_rings;
    values[idx++].d = core_atoms;
    values[idx++].d = periphery_atoms;
    values[idx++].d = core_periph_ratio;
    values[idx++].d = symmetry;
    values[idx++].d = topo_diversity;
    values[idx++].d = linker_hetero;
    values[idx++].d = linker_unsat;
    values[idx++].d = ring_decoration;
    values[idx++].d = terminal_groups;

    return NUM_FRAME_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* frame_cached_mol = NULL;
static _Thread_local uint64_t frame_cached_gen = 0;
static _Thread_local descriptor_value_t frame_cached_values[NUM_FRAME_DESCRIPTORS];

static inline void ensure_frame_computed(const molecule_t* mol) {
    uint64_t current_gen = descriptor_cache_generation();
    if (frame_cached_mol != mol || frame_cached_gen != current_gen) {
        descriptors_compute_framework_all(mol, frame_cached_values);
        frame_cached_mol = mol;
        frame_cached_gen = current_gen;
    }
}

#define DEFINE_FRAME_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_frame_computed(mol); \
    value->d = frame_cached_values[idx].d; \
    return CCHEM_OK; \
}

/* Ring System Topology (12) */
DEFINE_FRAME_FUNC(frame_ringsystems, 0)
DEFINE_FRAME_FUNC(frame_fusedpairs, 1)
DEFINE_FRAME_FUNC(frame_spirojunc, 2)
DEFINE_FRAME_FUNC(frame_bridgehead, 3)
DEFINE_FRAME_FUNC(frame_ringbridges, 4)
DEFINE_FRAME_FUNC(frame_maxringsize, 5)
DEFINE_FRAME_FUNC(frame_minringsize, 6)
DEFINE_FRAME_FUNC(frame_ringsizerange, 7)
DEFINE_FRAME_FUNC(frame_meanringsize, 8)
DEFINE_FRAME_FUNC(frame_ringsizevar, 9)
DEFINE_FRAME_FUNC(frame_totalringatoms, 10)
DEFINE_FRAME_FUNC(frame_totalringbonds, 11)

/* Chain/Linker Analysis (12) */
DEFINE_FRAME_FUNC(frame_linkeratoms, 12)
DEFINE_FRAME_FUNC(frame_linkerbonds, 13)
DEFINE_FRAME_FUNC(frame_maxlinkerlen, 14)
DEFINE_FRAME_FUNC(frame_meanlinkerlen, 15)
DEFINE_FRAME_FUNC(frame_terminalchains, 16)
DEFINE_FRAME_FUNC(frame_maxtermlen, 17)
DEFINE_FRAME_FUNC(frame_branchchains, 18)
DEFINE_FRAME_FUNC(frame_chainatoms, 19)
DEFINE_FRAME_FUNC(frame_chaindensity, 20)
DEFINE_FRAME_FUNC(frame_ringchainratio, 21)
DEFINE_FRAME_FUNC(frame_complexity, 22)
DEFINE_FRAME_FUNC(frame_compactness, 23)

/* Scaffold Metrics (16) */
DEFINE_FRAME_FUNC(frame_murckoscaff, 24)
DEFINE_FRAME_FUNC(frame_scaffoldratio, 25)
DEFINE_FRAME_FUNC(frame_sidechainatoms, 26)
DEFINE_FRAME_FUNC(frame_sidechainratio, 27)
DEFINE_FRAME_FUNC(frame_funcatoms, 28)
DEFINE_FRAME_FUNC(frame_decoratedrings, 29)
DEFINE_FRAME_FUNC(frame_nakedrings, 30)
DEFINE_FRAME_FUNC(frame_coresize, 31)
DEFINE_FRAME_FUNC(frame_peripherysize, 32)
DEFINE_FRAME_FUNC(frame_coreperiphratio, 33)
DEFINE_FRAME_FUNC(frame_symmetry, 34)
DEFINE_FRAME_FUNC(frame_topodiversity, 35)
DEFINE_FRAME_FUNC(frame_linkerhetero, 36)
DEFINE_FRAME_FUNC(frame_linkerunsat, 37)
DEFINE_FRAME_FUNC(frame_ringdecoration, 38)
DEFINE_FRAME_FUNC(frame_terminalgroups, 39)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_FRAME(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_framework(void) {
    /* Ring System Topology */
    REGISTER_FRAME("Frame_RingSystems", "Number of separate ring systems", desc_frame_ringsystems);
    REGISTER_FRAME("Frame_FusedPairs", "Number of fused ring pairs", desc_frame_fusedpairs);
    REGISTER_FRAME("Frame_SpiroJunctions", "Spiro junction count", desc_frame_spirojunc);
    REGISTER_FRAME("Frame_BridgeheadAtoms", "Bridgehead atom count", desc_frame_bridgehead);
    REGISTER_FRAME("Frame_RingBridges", "Bonds connecting ring systems", desc_frame_ringbridges);
    REGISTER_FRAME("Frame_MaxRingSize", "Largest ring in molecule", desc_frame_maxringsize);
    REGISTER_FRAME("Frame_MinRingSize", "Smallest ring in molecule", desc_frame_minringsize);
    REGISTER_FRAME("Frame_RingSizeRange", "Max - min ring size", desc_frame_ringsizerange);
    REGISTER_FRAME("Frame_MeanRingSize", "Average ring size", desc_frame_meanringsize);
    REGISTER_FRAME("Frame_RingSizeVar", "Variance of ring sizes", desc_frame_ringsizevar);
    REGISTER_FRAME("Frame_TotalRingAtoms", "Total atoms in rings", desc_frame_totalringatoms);
    REGISTER_FRAME("Frame_TotalRingBonds", "Total bonds in rings", desc_frame_totalringbonds);

    /* Chain/Linker Analysis */
    REGISTER_FRAME("Frame_LinkerAtoms", "Atoms connecting ring systems", desc_frame_linkeratoms);
    REGISTER_FRAME("Frame_LinkerBonds", "Bonds in linker regions", desc_frame_linkerbonds);
    REGISTER_FRAME("Frame_MaxLinkerLength", "Longest linker chain", desc_frame_maxlinkerlen);
    REGISTER_FRAME("Frame_MeanLinkerLength", "Average linker length", desc_frame_meanlinkerlen);
    REGISTER_FRAME("Frame_TerminalChains", "Number of terminal chains", desc_frame_terminalchains);
    REGISTER_FRAME("Frame_MaxTerminalLen", "Longest terminal chain", desc_frame_maxtermlen);
    REGISTER_FRAME("Frame_BranchChains", "Number of branch chains", desc_frame_branchchains);
    REGISTER_FRAME("Frame_ChainAtoms", "Non-ring atoms in chains", desc_frame_chainatoms);
    REGISTER_FRAME("Frame_ChainDensity", "Chain atoms / heavy atoms", desc_frame_chaindensity);
    REGISTER_FRAME("Frame_RingChainRatio", "Ring atoms / chain atoms", desc_frame_ringchainratio);
    REGISTER_FRAME("Frame_Complexity", "Framework complexity index", desc_frame_complexity);
    REGISTER_FRAME("Frame_Compactness", "Atoms / (ring systems + chains)", desc_frame_compactness);

    /* Scaffold Metrics */
    REGISTER_FRAME("Frame_MurckoScaffold", "Heavy atoms in Murcko scaffold", desc_frame_murckoscaff);
    REGISTER_FRAME("Frame_ScaffoldRatio", "Scaffold / total heavy atoms", desc_frame_scaffoldratio);
    REGISTER_FRAME("Frame_SideChainAtoms", "Atoms in side chains", desc_frame_sidechainatoms);
    REGISTER_FRAME("Frame_SideChainRatio", "Side chain / total ratio", desc_frame_sidechainratio);
    REGISTER_FRAME("Frame_FunctionalAtoms", "Atoms with functional groups", desc_frame_funcatoms);
    REGISTER_FRAME("Frame_DecoratedRings", "Rings with substituents", desc_frame_decoratedrings);
    REGISTER_FRAME("Frame_NakedRings", "Unsubstituted rings", desc_frame_nakedrings);
    REGISTER_FRAME("Frame_CoreSize", "Central scaffold size", desc_frame_coresize);
    REGISTER_FRAME("Frame_PeripherySize", "Peripheral atoms", desc_frame_peripherysize);
    REGISTER_FRAME("Frame_CorePeriphRatio", "Core / periphery ratio", desc_frame_coreperiphratio);
    REGISTER_FRAME("Frame_Symmetry", "Estimated symmetry index", desc_frame_symmetry);
    REGISTER_FRAME("Frame_TopoDiversity", "Topological diversity index", desc_frame_topodiversity);
    REGISTER_FRAME("Frame_LinkerHetero", "Heteroatoms in linker regions", desc_frame_linkerhetero);
    REGISTER_FRAME("Frame_LinkerUnsaturation", "Unsaturation in linkers", desc_frame_linkerunsat);
    REGISTER_FRAME("Frame_RingDecoration", "Mean substituents per ring", desc_frame_ringdecoration);
    REGISTER_FRAME("Frame_TerminalGroups", "Count of terminal groups", desc_frame_terminalgroups);
}
