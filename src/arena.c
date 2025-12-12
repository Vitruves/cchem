/**
 * @file arena.c
 * @brief Linear memory arena implementation
 */

#include <stdlib.h>
#include <string.h>
#include "cchem/arena.h"

arena_t* arena_create(size_t size) {
    arena_t* arena = (arena_t*)malloc(sizeof(arena_t));
    if (!arena) return NULL;

    arena->data = (uint8_t*)malloc(size);
    if (!arena->data) {
        free(arena);
        return NULL;
    }

    arena->size = size;
    arena->offset = 0;
    arena->peak = 0;

    return arena;
}

arena_t* arena_create_default(void) {
    return arena_create(ARENA_DEFAULT_SIZE);
}

void arena_free(arena_t* arena) {
    if (arena) {
        if (arena->data) free(arena->data);
        free(arena);
    }
}

/* Thread-local arena - one per thread */
static __thread arena_t* tl_arena = NULL;

arena_t* arena_get_thread_local(void) {
    if (!tl_arena) {
        tl_arena = arena_create_default();
    }
    return tl_arena;
}
