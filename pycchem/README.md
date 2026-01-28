# pycchem

Python bindings for [cchem](https://github.com/Vitruves/cchem) - a high-performance cheminformatics library.

## Features

- **Fast SMILES canonicalization** - 2-3x faster than RDKit
- **1600+ molecular descriptors** - Comprehensive descriptor set
- **Molecular sanitization** - Salt removal, neutralization, aromatization
- **Pure C backend** - No heavy dependencies

## Installation

```bash
pip install pycchem
```

### From Source

```bash
git clone https://github.com/Vitruves/cchem.git
cd cchem/pycchem
pip install .
```

## Quick Start

### Canonicalization

```python
import pycchem

# Canonicalize a SMILES
canon = pycchem.canonicalize("C(C)O")
print(canon)  # CCO

# Check if two SMILES are equivalent
pycchem.are_equivalent("CCO", "OCC")  # True

# Validate SMILES
pycchem.validate("CCO")  # True
pycchem.validate("invalid")  # False
```

### Sanitization

```python
import pycchem

# Remove salts and neutralize
clean = pycchem.sanitize("[Na+].CC(=O)[O-]")
print(clean)  # CC(=O)O

# Specific operations
pycchem.sanitize("C1=CC=CC=C1", operations=["aromatize"])  # c1ccccc1
```

### Molecular Descriptors

```python
import pycchem

# Create a molecule
mol = pycchem.Molecule("CCO")

# Compute all descriptors (1600+)
descriptors = mol.descriptors()
print(descriptors["MolecularWeight"])  # 46.069
print(descriptors["CarbonCount"])  # 2

# Compute specific descriptors
selected = mol.descriptors(["MolecularWeight", "WCLogP", "TPSA"])

# List available descriptors
names = pycchem.list_descriptors()
print(len(names))  # ~1600
```

### Batch Processing

```python
import pycchem

smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]
results = pycchem.canonicalize_batch(smiles_list)
```

## API Reference

### Functions

| Function | Description |
|----------|-------------|
| `canonicalize(smiles)` | Canonicalize SMILES string |
| `canonicalize_batch(list)` | Batch canonicalization |
| `sanitize(smiles, operations)` | Sanitize molecule |
| `validate(smiles)` | Check SMILES validity |
| `are_equivalent(s1, s2)` | Check if SMILES are equivalent |
| `compute_descriptor(smiles, name)` | Compute single descriptor |
| `compute_descriptors(smiles, names)` | Compute multiple descriptors |
| `list_descriptors()` | List all descriptor names |
| `version()` | Get library version |

### Molecule Class

```python
mol = pycchem.Molecule(smiles)
mol.canonical_smiles  # Canonical form
mol.descriptors()     # All descriptors
mol.descriptor(name)  # Single descriptor
```

### Sanitization Operations

| Operation | Description |
|-----------|-------------|
| `complete` | Full sanitization (default) |
| `unsalt` | Remove salts/counter-ions |
| `aromatize` | Perceive aromaticity |
| `kekulize` | Convert to Kekule form |
| `neutralize` | Neutralize charges |
| `normalize` | Normalize functional groups |
| `remove_stereo` | Remove stereochemistry |
| `remove_isotopes` | Remove isotope labels |
| `remove_h` | Remove explicit hydrogens |

## Performance

Benchmarks on MacBook Air M3 (8-core, 16GB RAM):

| Library | Throughput | vs RDKit |
|---------|------------|----------|
| pycchem | 28,311/s | **2.52x** |
| RDKit | 11,229/s | 1.00x |

## Requirements

- Python 3.9+
- cffi >= 1.15.0

## License

Apache License 2.0
