# Benchmarks CMakeLists.txt
#
# Build with:
#   cmake .. -DBUILD_BENCHMARKS=ON
#
# To enable RDKit/OpenBabel comparison:
#   cmake .. -DBUILD_BENCHMARKS=ON -DWITH_RDKIT=ON -DWITH_OPENBABEL=ON

cmake_minimum_required(VERSION 3.16)

# Options for external libraries
option(WITH_RDKIT "Build with RDKit support" OFF)
option(WITH_OPENBABEL "Build with OpenBabel support" OFF)

# C++ benchmark executable
add_executable(canonicalization_benchmark
    canonicalization_benchmark.cpp
)

# Link against cchem library
target_link_libraries(canonicalization_benchmark PRIVATE
    cchem_lib
    Threads::Threads
)

target_include_directories(canonicalization_benchmark PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# Try to find RDKit
if(WITH_RDKIT)
    # Manual search (more robust than CMake package which often has broken deps)
    set(RDKIT_FOUND FALSE)
        # Common RDKit installation paths
        find_path(RDKIT_INCLUDE_DIR
            NAMES GraphMol/GraphMol.h
            PATHS
                /usr/local/include/rdkit
                /usr/include/rdkit
                /opt/homebrew/include/rdkit
                $ENV{RDBASE}/Code
                $ENV{CONDA_PREFIX}/include/rdkit
        )

        find_library(RDKIT_GRAPHMOL_LIB
            NAMES RDKitGraphMol GraphMol
            PATHS
                /usr/local/lib
                /usr/lib
                /opt/homebrew/lib
                $ENV{RDBASE}/lib
                $ENV{CONDA_PREFIX}/lib
        )

        find_library(RDKIT_SMILESPARSE_LIB
            NAMES RDKitSmilesParse SmilesParse
            PATHS
                /usr/local/lib
                /usr/lib
                /opt/homebrew/lib
                $ENV{RDBASE}/lib
                $ENV{CONDA_PREFIX}/lib
        )

        find_library(RDKIT_RDGENERAL_LIB
            NAMES RDKitRDGeneral RDGeneral
            PATHS
                /usr/local/lib
                /usr/lib
                /opt/homebrew/lib
                $ENV{RDBASE}/lib
                $ENV{CONDA_PREFIX}/lib
        )

    if(RDKIT_INCLUDE_DIR AND RDKIT_GRAPHMOL_LIB AND RDKIT_SMILESPARSE_LIB)
        set(RDKIT_FOUND TRUE)
        set(RDKIT_INCLUDE_DIRS ${RDKIT_INCLUDE_DIR})
        set(RDKIT_LIBRARIES
            ${RDKIT_GRAPHMOL_LIB}
            ${RDKIT_SMILESPARSE_LIB}
            ${RDKIT_RDGENERAL_LIB}
        )
    endif()

    if(RDKIT_FOUND)
        message(STATUS "Found RDKit: ${RDKIT_INCLUDE_DIRS}")
        target_compile_definitions(canonicalization_benchmark PRIVATE HAVE_RDKIT)
        target_include_directories(canonicalization_benchmark PRIVATE ${RDKIT_INCLUDE_DIRS})
        target_link_libraries(canonicalization_benchmark PRIVATE ${RDKIT_LIBRARIES})
    else()
        message(WARNING "RDKit not found. Install with: brew install rdkit")
    endif()
endif()

# Try to find OpenBabel
if(WITH_OPENBABEL)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(OPENBABEL openbabel-3)
        if(NOT OPENBABEL_FOUND)
            pkg_check_modules(OPENBABEL openbabel)
        endif()
    endif()

    if(NOT OPENBABEL_FOUND)
        find_path(OPENBABEL_INCLUDE_DIR
            NAMES openbabel/mol.h
            PATHS
                /usr/local/include/openbabel3
                /usr/include/openbabel3
                /opt/homebrew/include/openbabel3
                /usr/local/include
                /usr/include
                /opt/homebrew/include
                $ENV{CONDA_PREFIX}/include
        )

        find_library(OPENBABEL_LIB
            NAMES openbabel openbabel-3
            PATHS
                /usr/local/lib
                /usr/lib
                /opt/homebrew/lib
                $ENV{CONDA_PREFIX}/lib
        )

        if(OPENBABEL_INCLUDE_DIR AND OPENBABEL_LIB)
            set(OPENBABEL_FOUND TRUE)
            set(OPENBABEL_INCLUDE_DIRS ${OPENBABEL_INCLUDE_DIR})
            set(OPENBABEL_LIBRARIES ${OPENBABEL_LIB})
        endif()
    endif()

    if(OPENBABEL_FOUND)
        message(STATUS "Found OpenBabel: ${OPENBABEL_INCLUDE_DIRS}")
        target_compile_definitions(canonicalization_benchmark PRIVATE HAVE_OPENBABEL)
        target_include_directories(canonicalization_benchmark PRIVATE ${OPENBABEL_INCLUDE_DIRS})
        target_link_libraries(canonicalization_benchmark PRIVATE ${OPENBABEL_LIBRARIES})
    else()
        message(WARNING "OpenBabel not found. Install with: brew install open-babel")
    endif()
endif()

# Set C++ standard for RDKit/OpenBabel compatibility
# RDKit requires C++20 for constexpr virtual functions
set_target_properties(canonicalization_benchmark PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)
