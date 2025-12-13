/**
 * @file mqn.c
 * @brief Molecular Quantum Numbers (MQN) descriptors
 *
 * 42 integer descriptors for chemical space navigation.
 * Reference: Nguyen et al. ChemMedChem 2009, 4, 1803-1805
 *
 * Categories:
 * - Atom counts (12): C, F, Cl, Br, I, S, P, acyclic/cyclic N/O
 * - Bond counts (6): acyclic/cyclic single/double/triple
 * - Polarity (8): HBA, HBD, negative/positive charges, etc.
 * - Topology (16): rings, rotatable bonds, etc.
 *
 * All O(n) complexity where n is atoms or bonds.
 */

#include <string.h>
#include <stdbool.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * MQN Atom Count Descriptors (MQN1-12)
 * ============================================================================ */

/* MQN1: Carbon count */
static cchem_status_t mqn_c(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_C) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN2: Fluorine count */
static cchem_status_t mqn_f(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_F) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN3: Chlorine count */
static cchem_status_t mqn_cl(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_Cl) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN4: Bromine count */
static cchem_status_t mqn_br(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_Br) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN5: Iodine count */
static cchem_status_t mqn_i(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_I) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN6: Sulfur count */
static cchem_status_t mqn_s(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_S) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN7: Phosphorus count */
static cchem_status_t mqn_p(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_P) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN8: Acyclic nitrogen count */
static cchem_status_t mqn_an(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_N && atom->ring_count == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN9: Cyclic nitrogen count */
static cchem_status_t mqn_cn(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_N && atom->ring_count > 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN10: Acyclic oxygen count */
static cchem_status_t mqn_ao(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_O && atom->ring_count == 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN11: Cyclic oxygen count */
static cchem_status_t mqn_co(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_O && atom->ring_count > 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN12: Heavy atom count (non-H) */
static cchem_status_t mqn_hac(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * MQN Bond Count Descriptors (MQN13-18)
 * ============================================================================ */

/* MQN13: Acyclic single bonds */
static cchem_status_t mqn_asb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (!bond->in_ring && (bond->type == BOND_SINGLE ||
            bond->type == BOND_UP || bond->type == BOND_DOWN)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN14: Acyclic double bonds */
static cchem_status_t mqn_adb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (!bond->in_ring && bond->type == BOND_DOUBLE) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN15: Acyclic triple bonds */
static cchem_status_t mqn_atb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (!bond->in_ring && bond->type == BOND_TRIPLE) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN16: Cyclic single bonds */
static cchem_status_t mqn_csb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->in_ring && (bond->type == BOND_SINGLE ||
            bond->type == BOND_UP || bond->type == BOND_DOWN)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN17: Cyclic double bonds */
static cchem_status_t mqn_cdb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->in_ring && bond->type == BOND_DOUBLE) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN18: Cyclic triple bonds (rare) */
static cchem_status_t mqn_ctb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->in_ring && bond->type == BOND_TRIPLE) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * MQN Polarity Descriptors (MQN19-26)
 * ============================================================================ */

/* Helper: total H on atom */
static int total_h(const molecule_t* mol, const atom_t* atom) {
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h++;
    }
    return h;
}

/* MQN19: H-bond acceptors (N, O with lone pairs) */
static cchem_status_t mqn_hba(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_N || atom->element == ELEM_O) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN20: H-bond donors (NH, OH) */
static cchem_status_t mqn_hbd(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if ((atom->element == ELEM_N || atom->element == ELEM_O) &&
            total_h(mol, atom) > 0) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN21: Negative formal charges */
static cchem_status_t mqn_neg(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].charge < 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN22: Positive formal charges */
static cchem_status_t mqn_pos(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].charge > 0) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN23: Acyclic monovalent atoms (halogens, terminal H donors) */
static cchem_status_t mqn_hbam(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count > 0) continue;

        /* Count heavy neighbors */
        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 1 && (atom->element == ELEM_N || atom->element == ELEM_O ||
            atom->element == ELEM_F || atom->element == ELEM_Cl ||
            atom->element == ELEM_Br || atom->element == ELEM_I)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN24: Acyclic bivalent nodes (ether O, thioether S, etc.) */
static cchem_status_t mqn_hbab(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count > 0) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 2 && (atom->element == ELEM_N || atom->element == ELEM_O ||
            atom->element == ELEM_S)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN25: Acyclic trivalent nodes */
static cchem_status_t mqn_hbat(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count > 0) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 3 && (atom->element == ELEM_N || atom->element == ELEM_C)) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN26: Acyclic tetravalent nodes */
static cchem_status_t mqn_hbaq(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count > 0) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 4 && atom->element == ELEM_C) {
            count++;
        }
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * MQN Topology Descriptors (MQN27-42)
 * ============================================================================ */

/* MQN27: Cyclic monovalent nodes */
static cchem_status_t mqn_hbcm(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count == 0) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN28: Cyclic bivalent nodes */
static cchem_status_t mqn_hbcb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count == 0) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 2) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN29: Cyclic trivalent nodes */
static cchem_status_t mqn_hbct(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count == 0) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN30: Cyclic tetravalent nodes */
static cchem_status_t mqn_hbcq(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->ring_count == 0) continue;

        int heavy = 0;
        for (int j = 0; j < atom->num_neighbors; j++) {
            if (mol->atoms[atom->neighbors[j]].element != ELEM_H) heavy++;
        }

        if (heavy == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN31: Ring count */
static cchem_status_t mqn_rbc(const molecule_t* mol, descriptor_value_t* value) {
    value->i = mol->num_rings;
    return CCHEM_OK;
}

/* MQN32: 3-membered rings */
static cchem_status_t mqn_r3(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size == 3) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN33: 4-membered rings */
static cchem_status_t mqn_r4(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size == 4) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN34: 5-membered rings */
static cchem_status_t mqn_r5(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size == 5) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN35: 6-membered rings */
static cchem_status_t mqn_r6(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size == 6) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN36: 7-membered rings */
static cchem_status_t mqn_r7(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size == 7) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN37: 8-membered rings */
static cchem_status_t mqn_r8(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size == 8) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN38: 9+ membered rings */
static cchem_status_t mqn_r9(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_rings; i++) {
        if (mol->rings[i].size >= 9) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN39: Aromatic atoms */
static cchem_status_t mqn_ara(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].aromatic) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN40: Aromatic bonds */
static cchem_status_t mqn_arb(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        if (mol->bonds[i].type == BOND_AROMATIC) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN41: Rotatable bonds */
static cchem_status_t mqn_rob(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];
        if (bond->type != BOND_SINGLE && bond->type != BOND_UP &&
            bond->type != BOND_DOWN) continue;
        if (bond->in_ring) continue;

        const atom_t* a1 = &mol->atoms[bond->atom1];
        const atom_t* a2 = &mol->atoms[bond->atom2];

        /* Count heavy neighbors */
        int h1 = 0, h2 = 0;
        for (int j = 0; j < a1->num_neighbors; j++) {
            if (mol->atoms[a1->neighbors[j]].element != ELEM_H) h1++;
        }
        for (int j = 0; j < a2->num_neighbors; j++) {
            if (mol->atoms[a2->neighbors[j]].element != ELEM_H) h2++;
        }

        if (h1 > 1 && h2 > 1) count++;
    }
    value->i = count;
    return CCHEM_OK;
}

/* MQN42: Topological polar surface area proxy (N + O count weighted) */
static cchem_status_t mqn_tpsa(const molecule_t* mol, descriptor_value_t* value) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_N) count += 26;  /* ~26 A^2 per N */
        if (atom->element == ELEM_O) count += 20;  /* ~20 A^2 per O */
    }
    value->i = count;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_MQN(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_COUNTS; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_mqn(void) {
    /* Atom counts (MQN1-12) */
    REGISTER_MQN("MQN1", "Carbon atoms", mqn_c);
    REGISTER_MQN("MQN2", "Fluorine atoms", mqn_f);
    REGISTER_MQN("MQN3", "Chlorine atoms", mqn_cl);
    REGISTER_MQN("MQN4", "Bromine atoms", mqn_br);
    REGISTER_MQN("MQN5", "Iodine atoms", mqn_i);
    REGISTER_MQN("MQN6", "Sulfur atoms", mqn_s);
    REGISTER_MQN("MQN7", "Phosphorus atoms", mqn_p);
    REGISTER_MQN("MQN8", "Acyclic nitrogen atoms", mqn_an);
    REGISTER_MQN("MQN9", "Cyclic nitrogen atoms", mqn_cn);
    REGISTER_MQN("MQN10", "Acyclic oxygen atoms", mqn_ao);
    REGISTER_MQN("MQN11", "Cyclic oxygen atoms", mqn_co);
    REGISTER_MQN("MQN12", "Heavy atom count", mqn_hac);

    /* Bond counts (MQN13-18) */
    REGISTER_MQN("MQN13", "Acyclic single bonds", mqn_asb);
    REGISTER_MQN("MQN14", "Acyclic double bonds", mqn_adb);
    REGISTER_MQN("MQN15", "Acyclic triple bonds", mqn_atb);
    REGISTER_MQN("MQN16", "Cyclic single bonds", mqn_csb);
    REGISTER_MQN("MQN17", "Cyclic double bonds", mqn_cdb);
    REGISTER_MQN("MQN18", "Cyclic triple bonds", mqn_ctb);

    /* Polarity (MQN19-26) */
    REGISTER_MQN("MQN19", "H-bond acceptor atoms", mqn_hba);
    REGISTER_MQN("MQN20", "H-bond donor atoms", mqn_hbd);
    REGISTER_MQN("MQN21", "Negative charges", mqn_neg);
    REGISTER_MQN("MQN22", "Positive charges", mqn_pos);
    REGISTER_MQN("MQN23", "Acyclic monovalent nodes", mqn_hbam);
    REGISTER_MQN("MQN24", "Acyclic bivalent nodes", mqn_hbab);
    REGISTER_MQN("MQN25", "Acyclic trivalent nodes", mqn_hbat);
    REGISTER_MQN("MQN26", "Acyclic tetravalent nodes", mqn_hbaq);

    /* Topology (MQN27-42) */
    REGISTER_MQN("MQN27", "Cyclic monovalent nodes", mqn_hbcm);
    REGISTER_MQN("MQN28", "Cyclic bivalent nodes", mqn_hbcb);
    REGISTER_MQN("MQN29", "Cyclic trivalent nodes", mqn_hbct);
    REGISTER_MQN("MQN30", "Cyclic tetravalent nodes", mqn_hbcq);
    REGISTER_MQN("MQN31", "Ring count", mqn_rbc);
    REGISTER_MQN("MQN32", "3-membered rings", mqn_r3);
    REGISTER_MQN("MQN33", "4-membered rings", mqn_r4);
    REGISTER_MQN("MQN34", "5-membered rings", mqn_r5);
    REGISTER_MQN("MQN35", "6-membered rings", mqn_r6);
    REGISTER_MQN("MQN36", "7-membered rings", mqn_r7);
    REGISTER_MQN("MQN37", "8-membered rings", mqn_r8);
    REGISTER_MQN("MQN38", "9+ membered rings", mqn_r9);
    REGISTER_MQN("MQN39", "Aromatic atoms", mqn_ara);
    REGISTER_MQN("MQN40", "Aromatic bonds", mqn_arb);
    REGISTER_MQN("MQN41", "Rotatable bonds", mqn_rob);
    REGISTER_MQN("MQN42", "TPSA proxy (N+O weighted)", mqn_tpsa);
}
