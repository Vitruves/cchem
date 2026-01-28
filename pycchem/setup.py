"""
Setup script for pycchem - Python bindings for cchem

This script handles:
1. Building the cchem C library if not already built
2. Building the cffi extension module
"""

import os
import subprocess
import sys
from pathlib import Path

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CChemBuildExt(build_ext):
    """Custom build_ext that builds cchem library first."""

    def run(self):
        # Get paths
        pycchem_dir = Path(__file__).parent.absolute()
        cchem_dir = pycchem_dir.parent  # cchem root directory
        build_dir = cchem_dir / "build"

        # Check if library exists
        lib_name = "libcchem_lib.a"
        if sys.platform == "win32":
            lib_name = "cchem_lib.lib"
        elif sys.platform == "darwin":
            lib_name = "libcchem_lib.a"

        lib_path = build_dir / lib_name

        if not lib_path.exists():
            print(f"Building cchem library...")
            build_dir.mkdir(exist_ok=True)

            # Configure
            cmake_args = [
                "cmake",
                str(cchem_dir),
                "-DCMAKE_BUILD_TYPE=Release",
                "-DENABLE_NATIVE_ARCH=OFF",
            ]

            if sys.platform == "win32":
                vcpkg_root = os.environ.get("VCPKG_INSTALLATION_ROOT", "")
                if vcpkg_root:
                    cmake_args.append(
                        f"-DCMAKE_TOOLCHAIN_FILE={vcpkg_root}/scripts/buildsystems/vcpkg.cmake"
                    )

            subprocess.check_call(cmake_args, cwd=build_dir)

            # Build
            build_args = ["cmake", "--build", ".", "--config", "Release"]
            if sys.platform != "win32":
                build_args.extend(["-j", str(os.cpu_count() or 1)])

            subprocess.check_call(build_args, cwd=build_dir)

        # Now build the cffi extension
        super().run()


def get_ext_modules():
    """Create cffi extension module."""
    from cffi import FFI

    ffi = FFI()

    # Read C declarations
    cdef = """
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

    typedef struct molecule_t molecule_t;

    typedef union {
        int64_t i;
        double d;
    } descriptor_value_t;

    typedef enum {
        DESC_VALUE_INT = 0,
        DESC_VALUE_DOUBLE = 1
    } descriptor_value_type_t;

    typedef struct {
        char name[64];
        char description[128];
        int category;
        descriptor_value_type_t value_type;
        ...;
    } descriptor_def_t;

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

    const char* cchem_version(void);
    cchem_status_t cchem_init(void);
    void cchem_cleanup(void);

    molecule_t* smiles_to_molecule(const char* smiles, char* error_buf, size_t error_buf_size);
    void molecule_free(molecule_t* mol);

    cchem_status_t smiles_validate(const char* smiles, char* error_buf, size_t error_buf_size);

    char* smiles_canonicalize(const char* smiles, const canon_options_t* options,
                              char* error_buf, size_t error_buf_size);
    char* molecule_to_canonical_smiles(const molecule_t* mol, const canon_options_t* options);
    cchem_status_t molecule_canonicalize(molecule_t* mol, const canon_options_t* options);

    bool smiles_are_equivalent(const char* smiles1, const char* smiles2);

    char* smiles_sanitize(const char* smiles, sanitize_flags_t flags,
                          char* error_buf, size_t error_buf_size);
    cchem_status_t sanitize_parse_flags(const char* str, sanitize_flags_t* flags);

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

    void free(void* ptr);
    """

    ffi.cdef(cdef)

    pycchem_dir = Path(__file__).parent.absolute()
    cchem_dir = pycchem_dir.parent
    build_dir = cchem_dir / "build"

    # Library and include paths
    include_dirs = [str(cchem_dir / "include")]
    library_dirs = [str(build_dir)]
    libraries = ["cchem_lib"]

    if sys.platform != "win32":
        libraries.append("m")

    # Platform-specific adjustments
    extra_link_args = []
    if sys.platform == "darwin":
        # Link against system frameworks
        extra_link_args.extend(["-framework", "CoreFoundation"])

    ffi.set_source(
        "pycchem._cchem_ffi",
        '#include "cchem/cchem.h"\n#include <stdlib.h>',
        libraries=libraries,
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        extra_link_args=extra_link_args,
    )

    return [ffi.distutils_extension()]


if __name__ == "__main__":
    setup(
        cffi_modules=["src/pycchem/_ffi_build.py:ffi"],
        cmdclass={"build_ext": CChemBuildExt},
    )
