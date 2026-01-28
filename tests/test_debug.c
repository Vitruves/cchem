#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cchem/cchem.h"

void print_atom_details(molecule_t* mol) {
    printf("Chiral atoms:\n");
    for (int i = 0; i < mol->num_atoms; i++) {
        atom_t* atom = &mol->atoms[i];
        if (atom->chirality != CHIRALITY_NONE) {
            printf("  Atom %d: chiral=%d stereo_nbrs=[", i, atom->chirality);
            for (int j = 0; j < atom->num_stereo_neighbors; j++) {
                printf("%d", atom->stereo_neighbors[j]);
                if (j < atom->num_stereo_neighbors - 1) printf(",");
            }
            printf("] nbrs=[");
            for (int j = 0; j < atom->num_neighbors; j++) {
                printf("%d", atom->neighbors[j]);
                if (j < atom->num_neighbors - 1) printf(",");
            }
            printf("]\n");
        }
    }
}

int main(void) {
    char error_buf[256];

    /* Test case with stereo inversion issue */
    const char* test_cases[] = {
        /* Original Test 1 - ring opening at stereocenter */
        "C[C@@H]1C2=NN=C(C3=CSC(C4=CC=C(F)C=C4F)=N3)N2CCN1C(=O)C1=CC=C(C2=CC=CS2)C=C1",
        /* Original Test 2 - bicyclic with ring closures */
        "CC1(C)[C@H]2CC=C(CN3CCC(NC4=NC(C5=CC(C(F)(F)F)=CC(C(F)(F)F)=C5)=NS4)CC3)[C@@H]1C2",
        NULL
    };

    for (int i = 0; test_cases[i]; i++) {
        const char* smiles = test_cases[i];
        printf("\n=== Test %d ===\n", i + 1);
        printf("Input: %s\n", smiles);

        molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
        if (!mol) {
            printf("Parse failed: %s\n", error_buf);
            continue;
        }

        print_atom_details(mol);

        char* canonical = molecule_to_canonical_smiles(mol, NULL);
        if (canonical) {
            printf("\nCanonical: %s\n", canonical);
            free(canonical);
        } else {
            printf("\nCanonicalization failed\n");
        }

        molecule_free(mol);
    }

    return 0;
}
