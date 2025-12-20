cmake_minimum_required(VERSION 3.16)
project(cchem VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Detect architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64|ARM64")
    set(CCHEM_ARCH_ARM64 TRUE)
    message(STATUS "Architecture: ARM64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
    set(CCHEM_ARCH_X64 TRUE)
    message(STATUS "Architecture: x86_64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "i[3-6]86|x86")
    set(CCHEM_ARCH_X86 TRUE)
    message(STATUS "Architecture: x86")
else()
    message(STATUS "Architecture: ${CMAKE_SYSTEM_PROCESSOR} (generic)")
endif()

# Platform-specific definitions
if(WIN32)
    # Windows: Enable POSIX-like functions and large file support
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
    if(MSVC)
        add_compile_definitions(_CRT_NONSTDC_NO_DEPRECATE)
    endif()
elseif(APPLE)
    # macOS: Darwin provides POSIX by default
    add_compile_definitions(_DARWIN_C_SOURCE)
else()
    # Linux/BSD: Enable GNU extensions for POSIX functions
    add_compile_definitions(_GNU_SOURCE)
endif()

# Enable Link Time Optimization (LTO) for Release builds
include(CheckIPOSupported)
check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_OUTPUT)
if(LTO_SUPPORTED AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "LTO enabled")
else()
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "LTO not supported: ${LTO_OUTPUT}")
    endif()
endif()

# Option to enable sanitizers (Debug only)
option(ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan (Debug only)" OFF)

# Option to build for native CPU (disable for portable binaries)
option(ENABLE_NATIVE_ARCH "Optimize for native CPU architecture" ON)

# Compiler-specific flags
if(MSVC)
    # Microsoft Visual C++
    message(STATUS "Compiler: MSVC ${MSVC_VERSION}")

    # Base warnings (W4 is high warning level)
    add_compile_options(/W4)

    # Disable specific noisy warnings
    add_compile_options(/wd4100)  # unreferenced formal parameter
    add_compile_options(/wd4244)  # conversion, possible loss of data
    add_compile_options(/wd4267)  # size_t to int conversion

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(/Od /Zi /RTC1)
    else()
        # Release optimizations
        add_compile_options(
            /O2           # Maximize speed
            /Oi           # Enable intrinsic functions
            /Ot           # Favor fast code
            /GL           # Whole program optimization
            /fp:fast      # Fast floating point
            /GS-          # Disable buffer security check for speed
        )
        add_link_options(/LTCG)  # Link-time code generation

        # Architecture-specific optimizations
        if(CCHEM_ARCH_X64 AND ENABLE_NATIVE_ARCH)
            add_compile_options(/arch:AVX2)
        elseif(CCHEM_ARCH_X86 AND ENABLE_NATIVE_ARCH)
            add_compile_options(/arch:SSE2)
        endif()
    endif()

elseif(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    # GCC, Clang, AppleClang
    message(STATUS "Compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")

    # Base warnings
    add_compile_options(-Wall -Wextra -Wpedantic)

    # Additional useful warnings
    add_compile_options(
        -Wformat=2
        -Wno-unused-parameter
        -Wno-sign-compare
    )

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        # Debug build
        add_compile_options(-g -O0)

        if(ENABLE_SANITIZERS)
            message(STATUS "Sanitizers enabled (ASan + UBSan)")
            add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
            add_link_options(-fsanitize=address,undefined)
        endif()
    else()
        # Release/default build: aggressive optimizations
        add_compile_options(
            -O3                    # Aggressive optimization
            -ffast-math            # Fast floating point (safe for descriptors)
            -fno-strict-aliasing   # Allow type punning patterns
            -funroll-loops         # Unroll small loops
            -ftree-vectorize       # Enable vectorization
        )

        # Architecture-specific optimizations
        if(ENABLE_NATIVE_ARCH)
            # Native architecture - best performance on build machine
            add_compile_options(-march=native -mtune=native)
        else()
            # Portable binaries with reasonable baseline
            if(CCHEM_ARCH_ARM64)
                if(APPLE)
                    # Apple Silicon baseline (M1+)
                    add_compile_options(-mcpu=apple-m1)
                else()
                    # Generic ARMv8.2-A with common extensions
                    add_compile_options(-march=armv8.2-a+fp16+dotprod)
                endif()
            elseif(CCHEM_ARCH_X64)
                # x86_64 baseline: SSE4.2 is universal on modern CPUs
                add_compile_options(-march=x86-64-v2)
            elseif(CCHEM_ARCH_X86)
                add_compile_options(-march=i686 -msse2)
            endif()
        endif()

        # Clang-specific optimizations
        if(CMAKE_C_COMPILER_ID MATCHES "Clang|AppleClang")
            add_compile_options(-fvectorize -fslp-vectorize)
        endif()

        # GCC-specific optimizations
        if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
            add_compile_options(-fipa-pta -fgraphite-identity)
        endif()
    endif()

    # Platform-specific linker flags
    if(APPLE)
        # Suppress Homebrew LLVM libunwind reexport warnings
        add_link_options(-Wl,-w)
    elseif(UNIX AND NOT APPLE)
        # Linux: Enable relro and now for security, disable for performance if needed
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
            add_link_options(-Wl,-O1 -Wl,--as-needed)
        endif()
    endif()

elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    # Intel compilers (icc/icx)
    message(STATUS "Compiler: Intel ${CMAKE_C_COMPILER_VERSION}")

    add_compile_options(-Wall)

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -O0)
    else()
        add_compile_options(
            -O3
            -xHost              # Optimize for host CPU
            -ipo                # Interprocedural optimization
            -fp-model=fast
            -qopt-zmm-usage=high
        )
    endif()
endif()

# Print build configuration summary
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Native arch optimization: ${ENABLE_NATIVE_ARCH}")

# Find required libraries
find_package(Threads REQUIRED)
find_package(JPEG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(CAIRO REQUIRED cairo)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${JPEG_INCLUDE_DIR})
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
    src/memory.c
    src/threading.c
    src/csv.c
    src/progress.c
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
    src/descriptors/cpsa.c
    src/descriptors/bcut_ext.c
    src/descriptors/moments.c
    src/descriptors/aromatic.c
    src/descriptors/atompairs_ext.c
    src/descriptors/framework.c
    src/descriptors/constitutional.c
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
target_link_libraries(cchem_lib PRIVATE Threads::Threads m ${JPEG_LIBRARIES} ${CAIRO_LIBRARIES})
target_link_directories(cchem_lib PRIVATE ${CAIRO_LIBRARY_DIRS})

# Main executable
add_executable(cchem src/main.c)
target_link_libraries(cchem PRIVATE cchem_lib Threads::Threads m ${JPEG_LIBRARIES} ${CAIRO_LIBRARIES})

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

add_executable(test_descriptors_extended tests/test_descriptors_extended.c)
target_link_libraries(test_descriptors_extended PRIVATE cchem_lib m)
add_test(NAME test_descriptors_extended COMMAND test_descriptors_extended)

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
