/**
 * @file compat.h
 * @brief Cross-platform compatibility layer for Windows/MSVC support
 */

#ifndef CCHEM_COMPAT_H
#define CCHEM_COMPAT_H

#ifdef _MSC_VER
/* ============================================================================
 * MSVC Compatibility
 * ============================================================================ */

#include <windows.h>

/* Thread-local storage */
#define __thread __declspec(thread)
#define _Thread_local __declspec(thread)

/* Function attributes - MSVC doesn't support GCC attributes */
#define __attribute__(x)

/* Case-insensitive string comparison */
#include <string.h>
#define strcasecmp _stricmp
#define strncasecmp _strnicmp

/* alloca - MSVC uses _alloca */
#include <malloc.h>
#define alloca _alloca

/* Builtin functions */
#include <intrin.h>
#include <stdlib.h>

static inline int __builtin_popcount(unsigned int x) {
    return (int)__popcnt(x);
}

static inline int __builtin_popcountll(unsigned long long x) {
    return (int)__popcnt64(x);
}

/* Use regular memset */
#define __builtin_memset memset

/* POSIX-like definitions */
#ifndef STDIN_FILENO
#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2
#endif

/* ssize_t for Windows */
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;

/* ============================================================================
 * C11 Atomics emulation for MSVC using Interlocked functions
 * ============================================================================ */

typedef volatile long atomic_size_t;
typedef volatile long atomic_int;

#define atomic_init(ptr, val) (*(ptr) = (long)(val))
#define atomic_store(ptr, val) InterlockedExchange((ptr), (long)(val))
#define atomic_load(ptr) InterlockedCompareExchange((ptr), 0, 0)
#define atomic_fetch_add(ptr, val) InterlockedExchangeAdd((ptr), (long)(val))
#define atomic_fetch_sub(ptr, val) InterlockedExchangeAdd((ptr), -(long)(val))

/* ============================================================================
 * Disable specific warnings
 * ============================================================================ */
#pragma warning(disable: 4100)  /* unreferenced formal parameter */
#pragma warning(disable: 4244)  /* conversion, possible loss of data */
#pragma warning(disable: 4267)  /* size_t to int conversion */
#pragma warning(disable: 4702)  /* unreachable code */

#else
/* ============================================================================
 * GCC/Clang - include standard POSIX headers
 * ============================================================================ */

#include <strings.h>
#include <alloca.h>
#include <stdatomic.h>

#endif /* _MSC_VER */

/* ============================================================================
 * Common includes
 * ============================================================================ */

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#endif /* CCHEM_COMPAT_H */
