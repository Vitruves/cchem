"""
pycchem - Python bindings for cchem cheminformatics library

High-performance SMILES canonicalization and molecular descriptor computation.

Example usage:
    >>> import pycchem
    >>> pycchem.canonicalize("c1ccccc1")
    'c1ccccc1'
    >>> mol = pycchem.Molecule("CCO")
    >>> mol.descriptors()["MolecularWeight"]
    46.069
"""

from __future__ import annotations

import atexit
from typing import Dict, List, Optional, Union

try:
    from ._cchem_ffi import ffi, lib
except ImportError as e:
    raise ImportError(
        "cchem C library not found. Make sure cchem is compiled and installed."
    ) from e


__version__ = "1.0.0"
__all__ = [
    "Molecule",
    "canonicalize",
    "canonicalize_batch",
    "sanitize",
    "validate",
    "are_equivalent",
    "compute_descriptor",
    "compute_descriptors",
    "list_descriptors",
    "version",
    "CChemError",
    "ParseError",
    "InvalidSmilesError",
]


# ============================================================================
# Exceptions
# ============================================================================

class CChemError(Exception):
    """Base exception for cchem errors."""
    pass


class ParseError(CChemError):
    """SMILES parsing error."""
    pass


class InvalidSmilesError(CChemError):
    """Invalid SMILES string."""
    pass


# Status code to exception mapping
_STATUS_EXCEPTIONS = {
    lib.CCHEM_ERROR_INVALID_INPUT: (CChemError, "Invalid input"),
    lib.CCHEM_ERROR_MEMORY: (CChemError, "Memory allocation failed"),
    lib.CCHEM_ERROR_PARSE: (ParseError, "Parse error"),
    lib.CCHEM_ERROR_INVALID_SMILES: (InvalidSmilesError, "Invalid SMILES"),
    lib.CCHEM_ERROR_RING_CLOSURE: (ParseError, "Ring closure error"),
    lib.CCHEM_ERROR_STEREO: (CChemError, "Stereochemistry error"),
    lib.CCHEM_ERROR_FILE_IO: (CChemError, "File I/O error"),
    lib.CCHEM_ERROR_THREAD: (CChemError, "Threading error"),
    lib.CCHEM_ERROR_NOT_IMPLEMENTED: (CChemError, "Not implemented"),
}


def _check_status(status: int, error_buf=None) -> None:
    """Check status code and raise appropriate exception."""
    if status == lib.CCHEM_OK:
        return
    exc_class, default_msg = _STATUS_EXCEPTIONS.get(status, (CChemError, "Unknown error"))
    if error_buf is not None:
        msg = ffi.string(error_buf).decode("utf-8")
        if msg:
            raise exc_class(msg)
    raise exc_class(default_msg)


def _encode_smiles(smiles: str) -> bytes:
    """Encode SMILES string to bytes."""
    if not isinstance(smiles, str):
        raise TypeError(f"Expected str, got {type(smiles).__name__}")
    return smiles.encode("utf-8")


# ============================================================================
# Library Initialization
# ============================================================================

# Initialize library and descriptors on import
lib.cchem_init()
lib.descriptors_init()


@atexit.register
def _cleanup():
    """Cleanup library resources on exit."""
    lib.descriptors_cleanup()
    lib.cchem_cleanup()


def version() -> str:
    """Get cchem library version."""
    return ffi.string(lib.cchem_version()).decode("utf-8")


# ============================================================================
# Molecule Class
# ============================================================================

class Molecule:
    """
    Molecular structure parsed from SMILES.

    Provides access to canonicalization, sanitization, and descriptor computation.

    Example:
        >>> mol = Molecule("CCO")
        >>> mol.canonical_smiles
        'CCO'
        >>> mol.num_atoms
        3
        >>> mol.descriptors()["MolecularWeight"]
        46.069
    """

    __slots__ = ("_mol", "_smiles", "_canonical_smiles", "_descriptors_cache")

    def __init__(self, smiles: str):
        """
        Create a Molecule from a SMILES string.

        Args:
            smiles: SMILES string to parse

        Raises:
            ParseError: If SMILES cannot be parsed
            InvalidSmilesError: If SMILES is invalid
        """
        self._smiles = smiles
        self._canonical_smiles: Optional[str] = None
        self._descriptors_cache: Optional[Dict[str, Union[int, float]]] = None

        error_buf = ffi.new("char[512]")
        self._mol = lib.smiles_to_molecule(_encode_smiles(smiles), error_buf, 512)

        if self._mol == ffi.NULL:
            msg = ffi.string(error_buf).decode("utf-8")
            raise ParseError(msg or f"Failed to parse SMILES: {smiles}")

    def __del__(self):
        """Free molecule memory."""
        if hasattr(self, "_mol") and self._mol != ffi.NULL:
            lib.molecule_free(self._mol)

    def __repr__(self) -> str:
        return f"Molecule({self._smiles!r})"

    def __str__(self) -> str:
        return self.canonical_smiles

    @property
    def smiles(self) -> str:
        """Original SMILES string."""
        return self._smiles

    @property
    def canonical_smiles(self) -> str:
        """Canonical SMILES representation."""
        if self._canonical_smiles is None:
            result = lib.molecule_to_canonical_smiles(self._mol, ffi.NULL)
            if result == ffi.NULL:
                raise CChemError("Failed to generate canonical SMILES")
            try:
                self._canonical_smiles = ffi.string(result).decode("utf-8")
            finally:
                lib.free(result)
        return self._canonical_smiles

    def descriptors(self, names: Optional[List[str]] = None) -> Dict[str, Union[int, float]]:
        """
        Compute molecular descriptors.

        Args:
            names: List of descriptor names to compute. If None, computes all descriptors.

        Returns:
            Dictionary mapping descriptor names to values.

        Example:
            >>> mol = Molecule("CCO")
            >>> mol.descriptors(["MolecularWeight", "CarbonCount"])
            {'MolecularWeight': 46.069, 'CarbonCount': 2}
        """
        if names is None:
            # Compute all descriptors (cached)
            if self._descriptors_cache is not None:
                return self._descriptors_cache.copy()

            num_desc = lib.descriptor_count()
            names_arr = ffi.new("const char*[]", num_desc)
            values_arr = ffi.new("descriptor_value_t[]", num_desc)

            count = lib.descriptors_compute_all(self._mol, names_arr, values_arr, num_desc)

            if count < 0:
                raise CChemError("Failed to compute descriptors")

            result = {}
            for i in range(count):
                name = ffi.string(names_arr[i]).decode("utf-8")
                desc_def = lib.descriptor_get(names_arr[i])
                if desc_def != ffi.NULL and desc_def.value_type == lib.DESC_VALUE_INT:
                    result[name] = values_arr[i].i
                else:
                    result[name] = values_arr[i].d

            self._descriptors_cache = result
            return result.copy()
        else:
            # Compute specific descriptors
            result = {}
            value = ffi.new("descriptor_value_t*")

            for name in names:
                status = lib.descriptor_compute(self._mol, _encode_smiles(name), value)
                if status == lib.CCHEM_OK:
                    desc_def = lib.descriptor_get(_encode_smiles(name))
                    if desc_def != ffi.NULL and desc_def.value_type == lib.DESC_VALUE_INT:
                        result[name] = value.i
                    else:
                        result[name] = value.d
                else:
                    result[name] = None

            return result

    def descriptor(self, name: str) -> Union[int, float]:
        """
        Compute a single descriptor.

        Args:
            name: Descriptor name

        Returns:
            Descriptor value

        Raises:
            CChemError: If descriptor computation fails
        """
        value = ffi.new("descriptor_value_t*")
        status = lib.descriptor_compute(self._mol, _encode_smiles(name), value)
        _check_status(status)

        desc_def = lib.descriptor_get(_encode_smiles(name))
        if desc_def != ffi.NULL and desc_def.value_type == lib.DESC_VALUE_INT:
            return value.i
        return value.d


# ============================================================================
# High-Level Functions
# ============================================================================

def canonicalize(
    smiles: str,
    *,
    isomeric: bool = True,
    kekulize: bool = False,
) -> str:
    """
    Canonicalize a SMILES string.

    Args:
        smiles: SMILES string to canonicalize
        isomeric: Include stereochemistry (default: True)
        kekulize: Convert aromatic to Kekule form (default: False)

    Returns:
        Canonical SMILES string

    Raises:
        ParseError: If SMILES cannot be parsed
        InvalidSmilesError: If SMILES is invalid

    Example:
        >>> canonicalize("C(C)O")
        'CCO'
        >>> canonicalize("c1ccccc1")
        'c1ccccc1'
    """
    error_buf = ffi.new("char[512]")

    # Use default options with modifications
    options = ffi.new("canon_options_t*")
    options.preserve_stereo = isomeric
    options.kekulize = kekulize
    options.remove_explicit_h = True
    options.remove_atom_maps = False
    options.canonical_stereo = True
    options.use_isomeric = isomeric
    options.use_isotopes = True
    options.use_charges = True

    result = lib.smiles_canonicalize(_encode_smiles(smiles), options, error_buf, 512)

    if result == ffi.NULL:
        msg = ffi.string(error_buf).decode("utf-8")
        raise ParseError(msg or f"Failed to canonicalize: {smiles}")

    try:
        return ffi.string(result).decode("utf-8")
    finally:
        lib.free(result)


def canonicalize_batch(
    smiles_list: List[str],
    *,
    isomeric: bool = True,
    n_threads: int = 1,
) -> List[Optional[str]]:
    """
    Canonicalize multiple SMILES strings.

    Args:
        smiles_list: List of SMILES strings
        isomeric: Include stereochemistry (default: True)
        n_threads: Number of threads (currently single-threaded in Python)

    Returns:
        List of canonical SMILES (None for failed entries)

    Example:
        >>> canonicalize_batch(["CCO", "c1ccccc1", "invalid"])
        ['CCO', 'c1ccccc1', None]
    """
    # TODO: Implement native C batch processing for better performance
    results = []
    for smi in smiles_list:
        try:
            results.append(canonicalize(smi, isomeric=isomeric))
        except CChemError:
            results.append(None)
    return results


# Sanitization flag mapping
_SANITIZE_FLAGS = {
    "unsalt": lib.SANITIZE_UNSALT,
    "aromatize": lib.SANITIZE_AROMATIZE,
    "kekulize": lib.SANITIZE_KEKULIZE,
    "neutralize": lib.SANITIZE_NEUTRALIZE,
    "remove_stereo": lib.SANITIZE_REMOVE_STEREO,
    "remove_isotopes": lib.SANITIZE_REMOVE_ISOTOPES,
    "remove_h": lib.SANITIZE_REMOVE_H,
    "normalize": lib.SANITIZE_NORMALIZE,
    "validate": lib.SANITIZE_VALIDATE,
    "complete": (
        lib.SANITIZE_UNSALT | lib.SANITIZE_AROMATIZE |
        lib.SANITIZE_NEUTRALIZE | lib.SANITIZE_NORMALIZE |
        lib.SANITIZE_VALIDATE
    ),
}


def sanitize(
    smiles: str,
    *,
    operations: Union[str, List[str]] = "complete",
) -> str:
    """
    Sanitize a SMILES string.

    Args:
        smiles: SMILES string to sanitize
        operations: Sanitization operations to apply. Can be:
            - "complete": Full sanitization (unsalt + aromatize + neutralize + normalize + validate)
            - List of operations: ["unsalt", "aromatize", "neutralize", "normalize",
                                   "kekulize", "remove_stereo", "remove_isotopes",
                                   "remove_h", "validate"]

    Returns:
        Sanitized SMILES string

    Example:
        >>> sanitize("[Na+].CC(=O)[O-]")
        'CC(=O)O'
        >>> sanitize("C1=CC=CC=C1", operations=["aromatize"])
        'c1ccccc1'
    """
    if isinstance(operations, str):
        if operations not in _SANITIZE_FLAGS:
            raise ValueError(f"Unknown sanitization operation: {operations}")
        flags = _SANITIZE_FLAGS[operations]
    else:
        flags = 0
        for op in operations:
            if op not in _SANITIZE_FLAGS:
                raise ValueError(f"Unknown sanitization operation: {op}")
            flags |= _SANITIZE_FLAGS[op]

    error_buf = ffi.new("char[512]")
    result = lib.smiles_sanitize(_encode_smiles(smiles), flags, error_buf, 512)

    if result == ffi.NULL:
        msg = ffi.string(error_buf).decode("utf-8")
        raise ParseError(msg or f"Failed to sanitize: {smiles}")

    try:
        return ffi.string(result).decode("utf-8")
    finally:
        lib.free(result)


def validate(smiles: str) -> bool:
    """
    Validate SMILES syntax.

    Args:
        smiles: SMILES string to validate

    Returns:
        True if valid, False otherwise

    Example:
        >>> validate("CCO")
        True
        >>> validate("invalid")
        False
    """
    error_buf = ffi.new("char[512]")
    status = lib.smiles_validate(_encode_smiles(smiles), error_buf, 512)
    return status == lib.CCHEM_OK


def are_equivalent(smiles1: str, smiles2: str) -> bool:
    """
    Check if two SMILES strings represent the same molecule.

    Args:
        smiles1: First SMILES string
        smiles2: Second SMILES string

    Returns:
        True if equivalent, False otherwise

    Example:
        >>> are_equivalent("CCO", "OCC")
        True
        >>> are_equivalent("CCO", "CCC")
        False
    """
    return lib.smiles_are_equivalent(_encode_smiles(smiles1), _encode_smiles(smiles2))


def compute_descriptor(smiles: str, name: str) -> Union[int, float]:
    """
    Compute a single descriptor for a SMILES string.

    Args:
        smiles: SMILES string
        name: Descriptor name

    Returns:
        Descriptor value

    Example:
        >>> compute_descriptor("CCO", "MolecularWeight")
        46.069
    """
    mol = Molecule(smiles)
    return mol.descriptor(name)


def compute_descriptors(
    smiles: str,
    names: Optional[List[str]] = None,
) -> Dict[str, Union[int, float]]:
    """
    Compute descriptors for a SMILES string.

    Args:
        smiles: SMILES string
        names: List of descriptor names (None for all)

    Returns:
        Dictionary of descriptor values

    Example:
        >>> compute_descriptors("CCO", ["MolecularWeight", "CarbonCount"])
        {'MolecularWeight': 46.069, 'CarbonCount': 2}
    """
    mol = Molecule(smiles)
    return mol.descriptors(names)


def list_descriptors() -> List[str]:
    """
    List all available descriptor names.

    Returns:
        List of descriptor names

    Example:
        >>> len(list_descriptors())
        1600
    """
    num_desc = lib.descriptor_count()
    descriptors = ffi.new("const descriptor_def_t*[]", num_desc)
    count = lib.descriptor_get_all(descriptors, num_desc)

    names = []
    for i in range(count):
        if descriptors[i] != ffi.NULL:
            names.append(ffi.string(descriptors[i].name).decode("utf-8"))

    return sorted(names)
