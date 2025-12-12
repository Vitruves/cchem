/**
 * @file arena.h
 * @brief Linear memory arena for zero-allocation batch processing
 *
 * The arena allocator provides fast, bump-pointer allocation with
 * O(1) reset. Perfect for batch processing where all allocations
 * can be freed together at the end of each batch item.
 */

#ifndef CCHEM_ARENA_H
#define CCHEM_ARENA_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* Default arena size: 1MB per thread */
#define ARENA_DEFAULT_SIZE (1024 * 1024)

/* Arena structure */
typedef struct arena {
    uint8_t* data;       /* Base pointer */
    size_t size;         /* Total size */
    size_t offset;       /* Current allocation offset */
    size_t peak;         /* Peak usage (for monitoring) */
} arena_t;

/* Create arena with specified size */
arena_t* arena_create(size_t size);

/* Create arena with default size */
arena_t* arena_create_default(void);

/* Free arena and all memory */
void arena_free(arena_t* arena);

/* Reset arena to empty (O(1) - just resets offset) */
static inline void arena_reset(arena_t* arena) {
    if (arena) {
        if (arena->offset > arena->peak) {
            arena->peak = arena->offset;
        }
        arena->offset = 0;
    }
}

/* Allocate memory from arena (returns NULL if out of space) */
static inline void* arena_alloc(arena_t* arena, size_t size) {
    if (!arena) return NULL;

    /* Align to 8 bytes */
    size_t aligned_size = (size + 7) & ~(size_t)7;

    if (arena->offset + aligned_size > arena->size) {
        return NULL;  /* Out of arena space */
    }

    void* ptr = arena->data + arena->offset;
    arena->offset += aligned_size;
    return ptr;
}

/* Allocate and zero memory */
static inline void* arena_calloc(arena_t* arena, size_t count, size_t elem_size) {
    size_t total = count * elem_size;
    void* ptr = arena_alloc(arena, total);
    if (ptr) {
        __builtin_memset(ptr, 0, total);
    }
    return ptr;
}

/* Get remaining space in arena */
static inline size_t arena_remaining(const arena_t* arena) {
    if (!arena) return 0;
    return arena->size - arena->offset;
}

/* Get peak usage (for tuning arena size) */
static inline size_t arena_peak_usage(const arena_t* arena) {
    if (!arena) return 0;
    return arena->peak > arena->offset ? arena->peak : arena->offset;
}

/* Thread-local arena access */
arena_t* arena_get_thread_local(void);

#endif /* CCHEM_ARENA_H */
