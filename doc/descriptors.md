# Adding New Descriptors to cchem

This guide explains how to add new molecular descriptors to the cchem descriptor system.

## Architecture Overview

The descriptor system is modular and extensible:

```
include/cchem/descriptors.h    - Main header with types and registry API
src/descriptors/registry.c     - Descriptor registration and lookup
src/descriptors/counts.c       - Count-based descriptors (atoms, bonds, rings)
src/descriptors/ratios.c       - Ratio-based descriptors (to be implemented)
```

All descriptors are computed from **canonical SMILES** via parsed `molecule_t` structures.

## Adding a New Descriptor

### Step 1: Choose the Right File

- **counts.c** - Integer counts (atoms, bonds, rings, functional groups)
- **ratios.c** - Floating-point ratios and fractions
- **properties.c** - Molecular properties (MW, LogP, etc.)
- **New file** - Create a new file for a new category

### Step 2: Implement the Compute Function

Each descriptor needs a compute function with this signature:

```c
static cchem_status_t desc_my_descriptor(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    // Compute descriptor value
    int64_t count = 0;

    // Example: count something
    for (int i = 0; i < mol->num_atoms; i++) {
        // Your computation logic
    }

    value->i = count;  // For integers
    // OR
    value->d = 3.14;   // For doubles

    return CCHEM_OK;
}
```

### Step 3: Register the Descriptor

Add registration in the appropriate `descriptors_register_*` function:

```c
void descriptors_register_counts(void) {
    // ... existing registrations ...

    REGISTER_COUNT_DESC("MyDescriptor", "Description of what it computes", desc_my_descriptor);
}
```

The `REGISTER_COUNT_DESC` macro creates a descriptor with:
- **name**: CamelCase name (used in CLI and CSV headers)
- **description**: Human-readable description
- **function**: Your compute function

### Step 4: For New Categories

If creating a new descriptor category (e.g., `properties.c`):

1. **Create the source file** `src/descriptors/properties.c`:

```c
#include <string.h>
#include "cchem/descriptors.h"

// Helper macro for this category
#define REGISTER_PROPERTY_DESC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

// Implement descriptors
static cchem_status_t desc_molecular_weight(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;
    value->d = molecule_calc_weight(mol);
    return CCHEM_OK;
}

// Registration function
void descriptors_register_properties(void) {
    REGISTER_PROPERTY_DESC("MolecularWeight", "Molecular weight in Da", desc_molecular_weight);
}
```

2. **Declare the registration function** in `descriptors.h`:

```c
void descriptors_register_properties(void);
```

3. **Call it from `descriptors_init()`** in `registry.c`:

```c
void descriptors_init(void) {
    // ...
    descriptors_register_counts();
    descriptors_register_ratios();
    descriptors_register_properties();  // Add new call
}
```

4. **Add to CMakeLists.txt**:

```cmake
set(CCHEM_SOURCES
    # ...
    src/descriptors/properties.c
)
```

## Available Data Structures

The `molecule_t` structure provides:

```c
mol->atoms[i].element     // Element type (ELEM_C, ELEM_N, etc.)
mol->atoms[i].aromatic    // Is aromatic
mol->atoms[i].charge      // Formal charge
mol->atoms[i].h_count     // Explicit H count
mol->atoms[i].implicit_h_count  // Computed implicit H
mol->atoms[i].num_neighbors     // Number of bonded atoms
mol->atoms[i].ring_count        // Rings containing this atom

mol->bonds[i].type        // BOND_SINGLE, BOND_DOUBLE, BOND_TRIPLE, BOND_AROMATIC
mol->bonds[i].aromatic    // Is aromatic bond
mol->bonds[i].in_ring     // Is in a ring

mol->rings[i].size        // Ring size
mol->rings[i].aromatic    // Is aromatic ring

mol->num_atoms            // Total atoms
mol->num_bonds            // Total bonds
mol->num_rings            // Total rings (SSSR)
```

## Optimization Guidelines

For ultra-fast descriptor computation:

1. **Single-pass algorithms**: Iterate through atoms/bonds once
2. **Use lookup tables**: For element-specific computations
3. **Avoid allocations**: Use stack variables, preallocate if needed
4. **Cache-friendly loops**: Sequential access patterns

Example optimized element count:

```c
static cchem_status_t count_element(const molecule_t* mol, element_t element,
                                    descriptor_value_t* value) {
    int64_t count = 0;
    const int n = mol->num_atoms;
    const atom_t* atoms = mol->atoms;  // Local pointer for cache efficiency

    for (int i = 0; i < n; i++) {
        count += (atoms[i].element == element);  // Branchless counting
    }

    value->i = count;
    return CCHEM_OK;
}
```

## Descriptor Naming Convention

- Use **CamelCase** for descriptor names
- Names should be descriptive but concise
- Use consistent suffixes:
  - `*Count` for integer counts
  - `*Ratio` for ratios/fractions
  - `*Index` for indices

Examples:
- `CarbonCount`, `NitrogenCount` (element counts)
- `AromaticRingCount` (specific ring type)
- `HBondAcceptorCount` (functional group)
- `AromaticRatio` (fraction of aromatic atoms)

## Testing Your Descriptor

1. **Single molecule test**:
```bash
./cchem compute -S "CCO" -d MyDescriptor -v
```

2. **Batch test**:
```bash
./cchem compute -f tests/test_descriptors.csv -s SMILES -d MyDescriptor -o output.csv
```

3. **Add unit tests** in `tests/test_descriptors.c`:
```c
void test_my_descriptor(void) {
    // Test with known molecules
    molecule_t* mol = smiles_to_molecule("CCO", NULL, 0);
    descriptor_value_t value;

    assert(desc_my_descriptor(mol, &value) == CCHEM_OK);
    assert(value.i == expected_value);

    molecule_free(mol);
}
```

## Canonical SMILES Handling

All descriptors receive molecules parsed from **canonical SMILES** output. The canonicalization process:

1. Standardizes atom ordering
2. Perceives aromaticity
3. Assigns ring membership
4. Computes implicit hydrogens

This ensures consistent descriptor values regardless of input SMILES format.

## Value Types

Two value types are supported:

- `DESC_VALUE_INT` - For integer descriptors (counts)
- `DESC_VALUE_DOUBLE` - For floating-point descriptors (ratios, properties)

Set the appropriate `value_type` when registering:

```c
def.value_type = DESC_VALUE_INT;    // For counts
def.value_type = DESC_VALUE_DOUBLE; // For ratios/properties
```

## Complete Example: Adding HalogenCount

```c
// In counts.c

static cchem_status_t desc_halogen_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    int64_t count = 0;
    const int n = mol->num_atoms;
    const atom_t* atoms = mol->atoms;

    for (int i = 0; i < n; i++) {
        element_t e = atoms[i].element;
        count += (e == ELEM_F || e == ELEM_Cl || e == ELEM_Br || e == ELEM_I);
    }

    value->i = count;
    return CCHEM_OK;
}

// In descriptors_register_counts():
REGISTER_COUNT_DESC("HalogenCount", "Number of halogen atoms (F, Cl, Br, I)", desc_halogen_count);
```

Then rebuild:
```bash
cd build && make
./cchem compute -S "ClCCBr" -d HalogenCount
# Output: HalogenCount: 2
```
