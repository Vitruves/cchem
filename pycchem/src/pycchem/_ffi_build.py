"""
cffi build script for pycchem bindings
"""

import os
from pathlib import Path
from cffi import FFI

# Find the cchem root directory (parent of pycchem)
_this_dir = Path(__file__).resolve().parent  # .../pycchem/src/pycchem
_pycchem_root = _this_dir.parent.parent  # .../pycchem/src/pycchem -> .../pycchem
_cchem_root = _pycchem_root.parent  # .../pycchem -> .../cchem

# Include and library directories
_include_dir = _cchem_root / "include"
_build_dir = _cchem_root / "build"

# Fallback to environment variables if set
if "CCHEM_INCLUDE_DIR" in os.environ:
    _include_dir = Path(os.environ["CCHEM_INCLUDE_DIR"])
if "CCHEM_LIBRARY_DIR" in os.environ:
    _build_dir = Path(os.environ["CCHEM_LIBRARY_DIR"])

ffi = FFI()

# C declarations - minimal interface for Python bindings
ffi.cdef("""
    /* Status codes */
    typedef enum {
        CCHEM_OK = 0,
        CCHEM_ERROR_INVALID_INPUT = -1,
        CCHEM_ERROR_MEMORY = -2,
        CCHEM_ERROR_PARSE = -3,
        CCHEM_ERROR_INVALID_SMILES = -4,
        CCHEM_ERROR_RING_CLOSURE = -5,
        CCHEM_ERROR_STEREO = -6,
        CCHEM_ERROR_FILE_IO = -7,
        CCHEM_ERROR_THREAD = -8,
        CCHEM_ERROR_NOT_IMPLEMENTED = -9
    } cchem_status_t;

    /* Opaque molecule handle */
    typedef struct molecule_t molecule_t;

    /* Descriptor value */
    typedef union {
        int64_t i;
        double d;
    } descriptor_value_t;

    typedef enum {
        DESC_VALUE_INT = 0,
        DESC_VALUE_DOUBLE = 1
    } descriptor_value_type_t;

    /* Descriptor definition (partial) */
    typedef struct {
        char name[64];
        char description[128];
        int category;
        descriptor_value_type_t value_type;
        ...;
    } descriptor_def_t;

    /* Sanitization flags */
    typedef enum {
        SANITIZE_NONE           = 0,
        SANITIZE_UNSALT         = 1,
        SANITIZE_AROMATIZE      = 2,
        SANITIZE_KEKULIZE       = 4,
        SANITIZE_NEUTRALIZE     = 8,
        SANITIZE_REMOVE_STEREO  = 16,
        SANITIZE_REMOVE_ISOTOPES = 32,
        SANITIZE_REMOVE_H       = 64,
        SANITIZE_NORMALIZE      = 128,
        SANITIZE_VALIDATE       = 256,
        SANITIZE_CLEANUP        = 512,
        SANITIZE_ALL            = 0xFFFF
    } sanitize_flags_t;

    /* Canonicalization options */
    typedef struct {
        bool preserve_stereo;
        bool kekulize;
        bool remove_explicit_h;
        bool remove_atom_maps;
        bool canonical_stereo;
        bool use_isomeric;
        bool use_isotopes;
        bool use_charges;
    } canon_options_t;

    extern const canon_options_t CANON_OPTIONS_DEFAULT;

    /* ========== Core Functions ========== */

    /* Library init/cleanup */
    const char* cchem_version(void);
    cchem_status_t cchem_init(void);
    void cchem_cleanup(void);

    /* Molecule parsing and freeing */
    molecule_t* smiles_to_molecule(const char* smiles, char* error_buf, size_t error_buf_size);
    void molecule_free(molecule_t* mol);

    /* SMILES validation */
    cchem_status_t smiles_validate(const char* smiles, char* error_buf, size_t error_buf_size);

    /* Canonicalization */
    char* smiles_canonicalize(const char* smiles, const canon_options_t* options,
                              char* error_buf, size_t error_buf_size);
    char* molecule_to_canonical_smiles(const molecule_t* mol, const canon_options_t* options);
    cchem_status_t molecule_canonicalize(molecule_t* mol, const canon_options_t* options);

    /* Equivalence check */
    bool smiles_are_equivalent(const char* smiles1, const char* smiles2);

    /* Sanitization */
    char* smiles_sanitize(const char* smiles, sanitize_flags_t flags,
                          char* error_buf, size_t error_buf_size);
    cchem_status_t sanitize_parse_flags(const char* str, sanitize_flags_t* flags);

    /* ========== Descriptors ========== */

    void descriptors_init(void);
    void descriptors_cleanup(void);
    int descriptor_count(void);
    const descriptor_def_t* descriptor_get(const char* name);

    cchem_status_t descriptor_compute(const molecule_t* mol,
                                      const char* name,
                                      descriptor_value_t* value);

    int descriptors_compute_all(const molecule_t* mol,
                                const char** out_names,
                                descriptor_value_t* out_values,
                                int max_descriptors);

    int descriptor_get_all(const descriptor_def_t** out_descriptors, int max_descriptors);

    /* ========== Memory ========== */
    void free(void* ptr);
""")

# Platform-specific library configuration
import sys
import glob as globmod

_libraries = []
_extra_link_args = []
_extra_compile_args = []
_sources = []

# Find all C source files from cchem
_src_dir = _cchem_root / "src"
if _src_dir.exists():
    # Compile cchem sources directly into the extension
    # Use relative paths from pycchem directory (required by setuptools)
    for pattern in ["canonicalizer/*.c", "descriptors/*.c", "utils/*.c", "splitter/*.c"]:
        full_pattern = str(_src_dir / pattern)
        for src_file in globmod.glob(full_pattern):
            # Convert absolute path to relative path from pycchem root
            rel_path = os.path.relpath(src_file, _pycchem_root)
            _sources.append(rel_path)
    # Add the API file
    _api_file = _src_dir / "cchem_api.c"
    if _api_file.exists():
        _sources.append(os.path.relpath(str(_api_file), _pycchem_root))
    # Exclude files that require external dependencies or aren't needed for Python bindings
    _exclude_patterns = [
        "/depictor/",      # Requires cairo
        "/commands/",      # CLI commands
        "parquet",         # Requires parquet lib
        "threading",       # Requires pthread
        "progress",        # Requires threading
        "csv",             # Requires threading for batch processing
        "splitter",        # Requires threading
    ]
    _sources = [s for s in _sources if not any(p in s for p in _exclude_patterns)]
    _extra_compile_args = ["-DCCHEM_BUILDING_EXTENSION"]

if sys.platform == "win32":
    pass
elif sys.platform == "darwin":
    _libraries.append("m")
else:
    _libraries.append("m")
    # Note: pthread not needed since threading files are excluded

ffi.set_source(
    "pycchem._cchem_ffi",
    """
    #include "cchem/cchem.h"
    #include <stdlib.h>
    """,
    sources=_sources,
    libraries=_libraries,
    include_dirs=[str(_include_dir)],
    library_dirs=[str(_build_dir)],
    extra_compile_args=_extra_compile_args,
    extra_link_args=_extra_link_args,
)

if __name__ == "__main__":
    ffi.compile(verbose=True)
