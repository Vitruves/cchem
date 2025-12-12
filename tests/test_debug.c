#include <stdio.h>
#include <stdlib.h>
#include "cchem/cchem.h"
#include "cchem/canonicalizer/ring_finder.h"

int main(void) {
    char error_buf[256];
    /* Test chirality issue */
    const char* test_smiles = "CN1CCN(C2=CC=CC=C2)[C@H](C)[C@H]1C";

    printf("Testing %s parsing...\n", test_smiles);
    molecule_t* mol = smiles_to_molecule(test_smiles, error_buf, sizeof(error_buf));
    if (!mol) {
        printf("Failed to parse: %s\n", error_buf);
        return 1;
    }
    printf("Parsed: %d atoms, %d bonds\n", mol->num_atoms, mol->num_bonds);
    for (int i = 0; i < mol->num_atoms; i++) {
        printf("  Atom %d: elem=%d neighbors=%d [", i,
               mol->atoms[i].element, mol->atoms[i].num_neighbors);
        for (int j = 0; j < mol->atoms[i].num_neighbors; j++) {
            printf("%d ", mol->atoms[i].neighbors[j]);
        }
        printf("] bonds=[");
        for (int j = 0; j < mol->atoms[i].num_neighbors; j++) {
            int b = mol->atoms[i].neighbor_bonds[j];
            printf("%d:%d ", b, mol->bonds[b].type);
        }
        printf("]\n");
    }
    printf("Bonds:\n");
    for (int i = 0; i < mol->num_bonds; i++) {
        printf("  Bond %d: %d-%d type=%d in_ring=%d\n", i, mol->bonds[i].atom1, mol->bonds[i].atom2, mol->bonds[i].type, mol->bonds[i].in_ring);
    }

    printf("\nRing info: %d rings\n", mol->num_rings);
    for (int r = 0; r < mol->num_rings; r++) {
        printf("  Ring %d (size %d): atoms ", r, mol->rings[r].size);
        for (int i = 0; i < mol->rings[r].size; i++) {
            printf("%d ", mol->rings[r].atoms[i]);
        }
        printf("\n");
    }

    printf("\nAtom chirality and stereo neighbors:\n");
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].chirality != CHIRALITY_NONE || mol->atoms[i].num_stereo_neighbors > 0) {
            printf("  Atom %d: chirality=%d stereo_neighbors=[", i, mol->atoms[i].chirality);
            for (int j = 0; j < mol->atoms[i].num_stereo_neighbors; j++) {
                printf("%d ", mol->atoms[i].stereo_neighbors[j]);
            }
            printf("] neighbors=[");
            for (int j = 0; j < mol->atoms[i].num_neighbors; j++) {
                printf("%d ", mol->atoms[i].neighbors[j]);
            }
            printf("]\n");
        }
    }

    printf("\nCanonicalizing...\n");
    molecule_t* mol2 = molecule_clone(mol);
    cchem_status_t status = molecule_canonicalize(mol2, NULL);
    printf("Status: %d\n", status);

    printf("\nCanon order: ");
    if (mol2->canon_order) {
        for (int i = 0; i < mol2->num_atoms; i++) {
            printf("%d ", mol2->canon_order[i]);
        }
    }
    printf("\n");

    printf("\nAtom ranks:\n");
    for (int i = 0; i < mol2->num_atoms; i++) {
        printf("  Atom %d: rank=%d inv=0x%llx\n", i,
               mol2->atoms[i].canon_rank, (unsigned long long)mol2->atoms[i].invariant);
    }

    printf("\nGenerating SMILES...\n");
    smiles_output_options_t opts = SMILES_OUTPUT_CANONICAL;
    char* smiles = molecule_to_smiles(mol2, &opts);
    printf("Result: %s\n", smiles ? smiles : "NULL");
    if (smiles) free(smiles);

    molecule_free(mol2);
    molecule_free(mol);
    return 0;
}
