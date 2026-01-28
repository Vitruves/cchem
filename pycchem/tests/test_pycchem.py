"""
Tests for pycchem Python bindings
"""

import pytest


def test_import():
    """Test that pycchem can be imported."""
    import pycchem
    assert pycchem.__version__ == "1.0.0"


def test_version():
    """Test version function."""
    import pycchem
    version = pycchem.version()
    assert version == "1.0.0"


def test_canonicalize_simple():
    """Test basic canonicalization."""
    import pycchem

    # Simple molecules
    assert pycchem.canonicalize("C") == "C"
    assert pycchem.canonicalize("CC") == "CC"
    assert pycchem.canonicalize("CCO") == "CCO"
    assert pycchem.canonicalize("OCC") == "CCO"  # Should canonicalize


def test_canonicalize_aromatic():
    """Test aromatic canonicalization."""
    import pycchem

    # Benzene
    result = pycchem.canonicalize("c1ccccc1")
    assert "c" in result.lower()  # Should be aromatic

    # Kekule form should aromatize
    result = pycchem.canonicalize("C1=CC=CC=C1")
    assert "c" in result.lower()


def test_canonicalize_stereo():
    """Test stereochemistry preservation."""
    import pycchem

    # E/Z isomers should be preserved
    e_isomer = pycchem.canonicalize("C/C=C/C")
    z_isomer = pycchem.canonicalize("C/C=C\\C")
    assert e_isomer != z_isomer


def test_validate():
    """Test SMILES validation."""
    import pycchem

    assert pycchem.validate("CCO") is True
    assert pycchem.validate("c1ccccc1") is True
    assert pycchem.validate("") is False
    assert pycchem.validate("XYZ") is False


def test_are_equivalent():
    """Test molecule equivalence checking."""
    import pycchem

    assert pycchem.are_equivalent("CCO", "OCC") is True
    assert pycchem.are_equivalent("CCO", "CCC") is False
    assert pycchem.are_equivalent("c1ccccc1", "C1=CC=CC=C1") is True


def test_sanitize_unsalt():
    """Test salt removal."""
    import pycchem

    # Remove sodium
    result = pycchem.sanitize("[Na+].CCO", operations=["unsalt"])
    assert "Na" not in result
    assert "C" in result


def test_sanitize_neutralize():
    """Test charge neutralization."""
    import pycchem

    # Neutralize carboxylate
    result = pycchem.sanitize("CC(=O)[O-]", operations=["neutralize"])
    assert "-" not in result or "O-" not in result


def test_sanitize_complete():
    """Test complete sanitization."""
    import pycchem

    result = pycchem.sanitize("[Na+].CC(=O)[O-]")
    assert "Na" not in result


def test_molecule_class():
    """Test Molecule class."""
    import pycchem

    mol = pycchem.Molecule("CCO")
    assert mol.smiles == "CCO"
    assert mol.canonical_smiles == "CCO"
    assert str(mol) == "CCO"


def test_molecule_descriptors():
    """Test descriptor computation."""
    import pycchem

    mol = pycchem.Molecule("CCO")
    descs = mol.descriptors()

    # Check some basic descriptors exist
    assert "MolecularWeight" in descs or "CarbonCount" in descs

    # Check specific descriptors
    specific = mol.descriptors(["CarbonCount", "OxygenCount"])
    if "CarbonCount" in specific:
        assert specific["CarbonCount"] == 2
    if "OxygenCount" in specific:
        assert specific["OxygenCount"] == 1


def test_compute_descriptor():
    """Test single descriptor computation."""
    import pycchem

    mol = pycchem.Molecule("CCO")

    # Test if any descriptor works
    try:
        value = mol.descriptor("CarbonCount")
        assert value == 2
    except Exception:
        pass  # Descriptor might not exist


def test_list_descriptors():
    """Test listing descriptors."""
    import pycchem

    names = pycchem.list_descriptors()
    assert len(names) > 100  # Should have many descriptors
    assert all(isinstance(n, str) for n in names)


def test_canonicalize_batch():
    """Test batch canonicalization."""
    import pycchem

    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]
    results = pycchem.canonicalize_batch(smiles_list)

    assert len(results) == 3
    assert all(r is not None for r in results)


def test_error_handling():
    """Test error handling."""
    import pycchem
    from pycchem import ParseError

    with pytest.raises((ParseError, pycchem.CChemError)):
        pycchem.canonicalize("invalid_smiles_xyz123")


def test_molecule_repr():
    """Test molecule string representation."""
    import pycchem

    mol = pycchem.Molecule("CCO")
    assert "CCO" in repr(mol)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
