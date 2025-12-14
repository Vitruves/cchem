cmake_minimum_required(VERSION 3.16)
project(cchem VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Enable POSIX extensions on Linux (strdup, strcasecmp, madvise, usleep, etc.)
if(NOT APPLE)
    add_compile_definitions(_GNU_SOURCE)
endif()

# Enable Link Time Optimization (LTO) for Release builds
include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_OUTPUT)
if(LTO_SUPPORTED AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "LTO enabled")
else()
    message(STATUS "LTO not supported: ${LTO_OUTPUT}")
endif()

# Compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    # Base warnings
    add_compile_options(-Wall -Wextra -Wpedantic)

    # Option to enable sanitizers (causes malloc warnings on macOS)
    option(ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan (Debug only)" OFF)

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        # Debug build: optional sanitizers
        add_compile_options(-g -O0)
        if(ENABLE_SANITIZERS)
            message(STATUS "Sanitizers enabled (ASan + UBSan)")
            add_compile_options(-fsanitize=address,undefined)
            add_link_options(-fsanitize=address,undefined)
        endif()
    else()
        # Release/default build: aggressive optimizations
        # -O3: Aggressive optimization including vectorization
        # -march=native: Use all CPU instructions available (AVX/AVX2/NEON)
        # -ffast-math: Allow FP reordering (safe for descriptors)
        # -fno-strict-aliasing: Allow type punning patterns
        # -funroll-loops: Unroll small loops
        add_compile_options(-O3 -march=native -ffast-math -fno-strict-aliasing -funroll-loops)
    endif()

    # Suppress Homebrew LLVM libunwind reexport warnings on macOS
    if(APPLE)
        add_link_options(-Wl,-w)
    endif()
endif()

# Find required libraries
find_package(Threads REQUIRED)
find_package(JPEG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(IGRAPH REQUIRED igraph)
pkg_check_modules(CAIRO REQUIRED cairo)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${JPEG_INCLUDE_DIR})
include_directories(${IGRAPH_INCLUDE_DIRS})
include_directories(${CAIRO_INCLUDE_DIRS})

# Library sources
set(CCHEM_SOURCES
    src/canonicalizer/atom.c
    src/canonicalizer/bond.c
    src/canonicalizer/molecule.c
    src/canonicalizer/lexer.c
    src/canonicalizer/parser.c
    src/canonicalizer/invariant.c
    src/canonicalizer/canon.c
    src/canonicalizer/smiles_writer.c
    src/canonicalizer/ring_finder.c
    src/canonicalizer/stereo.c
    src/canonicalizer/sanitize.c
    src/arena.c
    src/mol_soa.c
    src/csv.c
    src/progress.c
    src/parallel.c
    src/splitter.c
    src/descriptors/registry.c
    src/descriptors/counts.c
    src/descriptors/logp.c
    src/descriptors/ratios.c
    src/descriptors/electronic.c
    src/descriptors/steric.c
    src/descriptors/energetics.c
    src/descriptors/bitstring.c
    src/descriptors/acidbase.c
    src/descriptors/solubility.c
    src/descriptors/topology.c
    src/descriptors/fractional.c
    src/descriptors/hash.c
    src/descriptors/graph.c
    src/descriptors/autocorrelations.c
    src/descriptors/logd74.c
    src/descriptors/functional.c
    src/descriptors/pharmacophore.c
    src/descriptors/mqn.c
    src/descriptors/estate.c
    src/descriptors/vsa.c
    src/descriptors/chi.c
    src/descriptors/atompairs.c
    src/descriptors/bcut.c
    src/descriptors/zagreb.c
    src/descriptors/infocontent.c
    src/descriptors/walkcounts.c
    src/descriptors/estate_sums.c
    src/descriptors/eta.c
    src/descriptors/ringcomplexity.c
    src/depictor/types.c
    src/depictor/colors.c
    src/depictor/coords2d.c
    src/depictor/coords3d.c
    src/depictor/render.c
    src/depictor/depictor.c
    src/depictor/mmff94_typing.c
    src/depictor/mmff94_charges.c
    src/depictor/mmff94_energy.c
)

# Create static library
add_library(cchem_lib STATIC ${CCHEM_SOURCES})
target_link_libraries(cchem_lib PRIVATE Threads::Threads m ${JPEG_LIBRARIES} ${IGRAPH_LIBRARIES} ${CAIRO_LIBRARIES})
target_link_directories(cchem_lib PRIVATE ${IGRAPH_LIBRARY_DIRS} ${CAIRO_LIBRARY_DIRS})

# Main executable
add_executable(cchem src/main.c)
target_link_libraries(cchem PRIVATE cchem_lib Threads::Threads m ${JPEG_LIBRARIES} ${IGRAPH_LIBRARIES} ${CAIRO_LIBRARIES})
target_link_directories(cchem PRIVATE ${IGRAPH_LIBRARY_DIRS})

# Tests
enable_testing()

add_executable(test_lexer tests/test_lexer.c)
target_link_libraries(test_lexer PRIVATE cchem_lib m)
add_test(NAME test_lexer COMMAND test_lexer)

add_executable(test_parser tests/test_parser.c)
target_link_libraries(test_parser PRIVATE cchem_lib m)
add_test(NAME test_parser COMMAND test_parser)

add_executable(test_molecule tests/test_molecule.c)
target_link_libraries(test_molecule PRIVATE cchem_lib m)
add_test(NAME test_molecule COMMAND test_molecule)

add_executable(test_canon tests/test_canon.c)
target_link_libraries(test_canon PRIVATE cchem_lib m)
add_test(NAME test_canon COMMAND test_canon)

add_executable(test_integration tests/test_integration.c)
target_link_libraries(test_integration PRIVATE cchem_lib m)
add_test(NAME test_integration COMMAND test_integration)

add_executable(test_stereo tests/test_stereo.c)
target_link_libraries(test_stereo PRIVATE cchem_lib m)
add_test(NAME test_stereo COMMAND test_stereo)

add_executable(test_debug tests/test_debug.c)
target_link_libraries(test_debug PRIVATE cchem_lib m)

add_executable(test_descriptors tests/test_descriptors.c)
target_link_libraries(test_descriptors PRIVATE cchem_lib m)
add_test(NAME test_descriptors COMMAND test_descriptors)

add_executable(test_mmff94 tests/test_mmff94.c)
target_link_libraries(test_mmff94 PRIVATE cchem_lib m)
add_test(NAME test_mmff94 COMMAND test_mmff94)

add_executable(test_sanitize tests/test_sanitize.c)
target_link_libraries(test_sanitize PRIVATE cchem_lib m)
add_test(NAME test_sanitize COMMAND test_sanitize)

# LogD 7.4 development test
#add_executable(logd74_test temp/logd74.c)
#target_link_libraries(logd74_test PRIVATE cchem_lib m)

# LogD 7.4 verification test
#add_executable(verify_logd74 temp/verify_logd74.c)
#target_link_libraries(verify_logd74 PRIVATE cchem_lib m)

# LogP neural network training
#add_executable(logp_nn temp/logp_nn.c)
#target_link_libraries(logp_nn PRIVATE cchem_lib m)

# LogP verification test
#add_executable(verify_logp temp/verify_logp.c)
#target_link_libraries(verify_logp PRIVATE cchem_lib m)

# Install
install(TARGETS cchem RUNTIME DESTINATION bin)
install(DIRECTORY include/cchem DESTINATION include)
