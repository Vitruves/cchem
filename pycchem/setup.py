"""
Setup script for pycchem - Python bindings for cchem

The C sources are compiled directly into the cffi extension module,
so no pre-built library is required.
"""

from setuptools import setup

if __name__ == "__main__":
    setup(
        cffi_modules=["src/pycchem/_ffi_build.py:ffi"],
    )
