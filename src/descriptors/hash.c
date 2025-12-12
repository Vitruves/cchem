/**
 * @file hash.c
 * @brief Ultra-fast hash-based molecular descriptors using xxHash
 *
 * These descriptors extract structural information from SMILES strings
 * using xxHash for speed. Designed to capture patterns relevant for
 * LogD prediction: ionizable groups, aromaticity, polarity distribution.
 *
 * All descriptors are O(n) where n is SMILES length.
 */

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/cchem.h"

/* ============================================================================
 * xxHash32 Implementation (BSD-2-Clause, Yann Collet)
 * Embedded for portability - extremely fast non-cryptographic hash
 * ============================================================================ */

#define XXH_PRIME32_1  0x9E3779B1U
#define XXH_PRIME32_2  0x85EBCA77U
#define XXH_PRIME32_3  0xC2B2AE3DU
#define XXH_PRIME32_4  0x27D4EB2FU
#define XXH_PRIME32_5  0x165667B1U

static inline uint32_t xxh32_rotl(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

static inline uint32_t xxh32_round(uint32_t acc, uint32_t input) {
    acc += input * XXH_PRIME32_2;
    acc = xxh32_rotl(acc, 13);
    acc *= XXH_PRIME32_1;
    return acc;
}

static inline uint32_t xxh32_avalanche(uint32_t h32) {
    h32 ^= h32 >> 15;
    h32 *= XXH_PRIME32_2;
    h32 ^= h32 >> 13;
    h32 *= XXH_PRIME32_3;
    h32 ^= h32 >> 16;
    return h32;
}

static uint32_t xxhash32(const void* input, size_t len, uint32_t seed) {
    const uint8_t* p = (const uint8_t*)input;
    const uint8_t* end = p + len;
    uint32_t h32;

    if (len >= 16) {
        const uint8_t* limit = end - 16;
        uint32_t v1 = seed + XXH_PRIME32_1 + XXH_PRIME32_2;
        uint32_t v2 = seed + XXH_PRIME32_2;
        uint32_t v3 = seed + 0;
        uint32_t v4 = seed - XXH_PRIME32_1;

        do {
            uint32_t k1, k2, k3, k4;
            memcpy(&k1, p, 4); p += 4;
            memcpy(&k2, p, 4); p += 4;
            memcpy(&k3, p, 4); p += 4;
            memcpy(&k4, p, 4); p += 4;
            v1 = xxh32_round(v1, k1);
            v2 = xxh32_round(v2, k2);
            v3 = xxh32_round(v3, k3);
            v4 = xxh32_round(v4, k4);
        } while (p <= limit);

        h32 = xxh32_rotl(v1, 1) + xxh32_rotl(v2, 7) +
              xxh32_rotl(v3, 12) + xxh32_rotl(v4, 18);
    } else {
        h32 = seed + XXH_PRIME32_5;
    }

    h32 += (uint32_t)len;

    while (p + 4 <= end) {
        uint32_t k1;
        memcpy(&k1, p, 4);
        h32 += k1 * XXH_PRIME32_3;
        h32 = xxh32_rotl(h32, 17) * XXH_PRIME32_4;
        p += 4;
    }

    while (p < end) {
        h32 += (*p++) * XXH_PRIME32_5;
        h32 = xxh32_rotl(h32, 11) * XXH_PRIME32_1;
    }

    return xxh32_avalanche(h32);
}

/* Hash to normalized double [0, 1) */
static inline double hash_to_double(uint32_t h) {
    return (double)h / 4294967296.0;
}

/* ============================================================================
 * SMILES Generation Helper
 * ============================================================================ */

static char* get_canonical_smiles(const molecule_t* mol) {
    char* smiles = molecule_to_smiles(mol, NULL);
    return smiles;
}

/* ============================================================================
 * Hash Statistics Structure
 * ============================================================================ */

#define NUM_HASH_DESCRIPTORS 30
#define MINHASH_SIZE 8

typedef struct {
    /* Full SMILES hashes with different seeds */
    uint32_t full_hash_s0;
    uint32_t full_hash_s1;
    uint32_t full_hash_s2;

    /* N-gram hash sums */
    uint64_t bigram_sum;
    uint64_t trigram_sum;
    uint64_t quadgram_sum;
    int n_bigrams;
    int n_trigrams;
    int n_quadgrams;

    /* Pattern-specific hashes */
    uint32_t aromatic_hash;      /* lowercase chars only */
    uint32_t aliphatic_hash;     /* uppercase chars only */
    uint32_t hetero_hash;        /* N, O, S, P, F, Cl, Br, I patterns */
    uint32_t ring_hash;          /* digits (ring closures) */
    uint32_t branch_hash;        /* parentheses patterns */

    /* MinHash signatures for substructure similarity */
    uint32_t minhash[MINHASH_SIZE];

    /* Sliding window statistics */
    uint32_t window4_min;
    uint32_t window4_max;
    uint64_t window4_sum;
    int n_windows;

    /* Position-weighted hashes */
    uint64_t pos_weighted_sum;

    /* XOR chain for order sensitivity */
    uint32_t xor_chain;

    /* Character class transitions */
    uint32_t transition_hash;

    /* Ionizable pattern hashes */
    uint32_t acidic_hash;        /* COOH, SO3H patterns */
    uint32_t basic_hash;         /* amine patterns */

    size_t smiles_len;
} hash_stats_t;

static void collect_hash_stats(const char* smiles, hash_stats_t* s) {
    memset(s, 0, sizeof(hash_stats_t));

    if (!smiles || !smiles[0]) return;

    size_t len = strlen(smiles);
    s->smiles_len = len;

    /* Initialize minhash to max */
    for (int i = 0; i < MINHASH_SIZE; i++) {
        s->minhash[i] = UINT32_MAX;
    }
    s->window4_min = UINT32_MAX;
    s->window4_max = 0;

    /* Full SMILES hashes with different seeds */
    s->full_hash_s0 = xxhash32(smiles, len, 0);
    s->full_hash_s1 = xxhash32(smiles, len, 0x12345678);
    s->full_hash_s2 = xxhash32(smiles, len, 0xDEADBEEF);

    /* Build filtered strings for pattern hashes */
    char aromatic[512] = {0};
    char aliphatic[512] = {0};
    char hetero[512] = {0};
    char rings[64] = {0};
    char branches[128] = {0};
    char acidic[128] = {0};
    char basic[128] = {0};

    int ai = 0, ali = 0, hi = 0, ri = 0, bi = 0, aci = 0, basi = 0;

    char prev_class = 0;  /* 'a'=aromatic, 'A'=aliphatic, 'd'=digit, 'b'=branch, 'o'=other */

    for (size_t i = 0; i < len && i < 500; i++) {
        char c = smiles[i];
        char curr_class;

        /* Classify character */
        if (islower(c)) {
            if (ai < 510) aromatic[ai++] = c;
            curr_class = 'a';
        } else if (isupper(c)) {
            if (ali < 510) aliphatic[ali++] = c;
            curr_class = 'A';
        } else if (isdigit(c)) {
            if (ri < 62) rings[ri++] = c;
            curr_class = 'd';
        } else if (c == '(' || c == ')') {
            if (bi < 126) branches[bi++] = c;
            curr_class = 'b';
        } else {
            curr_class = 'o';
        }

        /* Track heteroatoms */
        if (c == 'N' || c == 'n' || c == 'O' || c == 'o' ||
            c == 'S' || c == 's' || c == 'P' || c == 'p' ||
            c == 'F' || c == 'I') {
            if (hi < 510) hetero[hi++] = c;
        }

        /* Acidic patterns: look for C(=O)O, S(=O)(=O)O */
        if (i >= 2 && c == 'O') {
            if (smiles[i-1] == ')' || (i >= 4 && smiles[i-2] == '=')) {
                if (aci < 126) acidic[aci++] = c;
            }
        }

        /* Basic patterns: N not in amide */
        if ((c == 'N' || c == 'n') && (i == 0 || smiles[i-1] != '=')) {
            if (basi < 126) basic[basi++] = c;
        }

        /* Track transitions */
        if (prev_class && prev_class != curr_class) {
            char trans[3] = {prev_class, curr_class, 0};
            s->transition_hash ^= xxhash32(trans, 2, (uint32_t)i);
        }
        prev_class = curr_class;

        /* XOR chain */
        s->xor_chain ^= xxhash32(&c, 1, (uint32_t)i * 31);

        /* Position-weighted hash */
        s->pos_weighted_sum += (uint64_t)xxhash32(&c, 1, 0) * (i + 1);
    }

    /* Hash pattern strings */
    if (ai > 0) s->aromatic_hash = xxhash32(aromatic, ai, 0xAA);
    if (ali > 0) s->aliphatic_hash = xxhash32(aliphatic, ali, 0xBB);
    if (hi > 0) s->hetero_hash = xxhash32(hetero, hi, 0xCC);
    if (ri > 0) s->ring_hash = xxhash32(rings, ri, 0xDD);
    if (bi > 0) s->branch_hash = xxhash32(branches, bi, 0xEE);
    if (aci > 0) s->acidic_hash = xxhash32(acidic, aci, 0xAC);
    if (basi > 0) s->basic_hash = xxhash32(basic, basi, 0xBA);

    /* N-gram hashes */
    for (size_t i = 0; i + 1 < len; i++) {
        uint32_t h = xxhash32(smiles + i, 2, 0x2222);
        s->bigram_sum += h;
        s->n_bigrams++;

        /* MinHash update */
        for (int m = 0; m < MINHASH_SIZE; m++) {
            uint32_t mh = xxhash32(smiles + i, 2, (uint32_t)m * 0x1111);
            if (mh < s->minhash[m]) s->minhash[m] = mh;
        }
    }

    for (size_t i = 0; i + 2 < len; i++) {
        s->trigram_sum += xxhash32(smiles + i, 3, 0x3333);
        s->n_trigrams++;
    }

    for (size_t i = 0; i + 3 < len; i++) {
        s->quadgram_sum += xxhash32(smiles + i, 4, 0x4444);
        s->n_quadgrams++;
    }

    /* Sliding window of 4 characters */
    for (size_t i = 0; i + 3 < len; i++) {
        uint32_t h = xxhash32(smiles + i, 4, 0xAAAA);
        s->window4_sum += h;
        s->n_windows++;
        if (h < s->window4_min) s->window4_min = h;
        if (h > s->window4_max) s->window4_max = h;
    }
}

/* ============================================================================
 * Batch Computation
 * ============================================================================ */

int descriptors_compute_hash_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    char* smiles = get_canonical_smiles(mol);
    if (!smiles) {
        for (int i = 0; i < NUM_HASH_DESCRIPTORS; i++) values[i].d = 0.0;
        return NUM_HASH_DESCRIPTORS;
    }

    hash_stats_t s;
    collect_hash_stats(smiles, &s);
    free(smiles);

    int idx = 0;

    /* 1-3: Full hash normalized (different seeds for diversity) */
    values[idx++].d = hash_to_double(s.full_hash_s0);      /* HashFull0 */
    values[idx++].d = hash_to_double(s.full_hash_s1);      /* HashFull1 */
    values[idx++].d = hash_to_double(s.full_hash_s2);      /* HashFull2 */

    /* 4-6: N-gram mean hashes */
    values[idx++].d = s.n_bigrams > 0 ?
        (double)(s.bigram_sum / s.n_bigrams) / 4294967296.0 : 0.0;   /* HashBigram */
    values[idx++].d = s.n_trigrams > 0 ?
        (double)(s.trigram_sum / s.n_trigrams) / 4294967296.0 : 0.0; /* HashTrigram */
    values[idx++].d = s.n_quadgrams > 0 ?
        (double)(s.quadgram_sum / s.n_quadgrams) / 4294967296.0 : 0.0; /* HashQuadgram */

    /* 7-11: Pattern-specific hashes normalized */
    values[idx++].d = hash_to_double(s.aromatic_hash);     /* HashAromatic */
    values[idx++].d = hash_to_double(s.aliphatic_hash);    /* HashAliphatic */
    values[idx++].d = hash_to_double(s.hetero_hash);       /* HashHetero */
    values[idx++].d = hash_to_double(s.ring_hash);         /* HashRing */
    values[idx++].d = hash_to_double(s.branch_hash);       /* HashBranch */

    /* 12-19: MinHash signatures (normalized) */
    for (int i = 0; i < MINHASH_SIZE; i++) {
        values[idx++].d = s.minhash[i] == UINT32_MAX ? 0.0 :
            hash_to_double(s.minhash[i]);                   /* MinHash0-7 */
    }

    /* 20-22: Sliding window statistics */
    values[idx++].d = s.window4_min == UINT32_MAX ? 0.0 :
        hash_to_double(s.window4_min);                      /* HashWindowMin */
    values[idx++].d = hash_to_double(s.window4_max);       /* HashWindowMax */
    values[idx++].d = s.n_windows > 0 ?
        (double)(s.window4_sum / s.n_windows) / 4294967296.0 : 0.0; /* HashWindowMean */

    /* 23: Position-weighted hash */
    values[idx++].d = s.smiles_len > 0 ?
        (double)(s.pos_weighted_sum / s.smiles_len) / 4294967296.0 : 0.0; /* HashPosWeight */

    /* 24: XOR chain */
    values[idx++].d = hash_to_double(s.xor_chain);         /* HashXorChain */

    /* 25: Transition hash */
    values[idx++].d = hash_to_double(s.transition_hash);   /* HashTransition */

    /* 26-27: Ionizable pattern hashes */
    values[idx++].d = hash_to_double(s.acidic_hash);       /* HashAcidic */
    values[idx++].d = hash_to_double(s.basic_hash);        /* HashBasic */

    /* 28: Hash entropy proxy (window range / mean) */
    double range = s.window4_max - (s.window4_min == UINT32_MAX ? 0 : s.window4_min);
    double mean = s.n_windows > 0 ? (double)s.window4_sum / s.n_windows : 1.0;
    values[idx++].d = mean > 0 ? range / mean : 0.0;       /* HashEntropy */

    /* 29: Hash collision indicator (minhash variance) */
    double mh_sum = 0, mh_sq_sum = 0;
    for (int i = 0; i < MINHASH_SIZE; i++) {
        double v = hash_to_double(s.minhash[i]);
        mh_sum += v;
        mh_sq_sum += v * v;
    }
    double mh_mean = mh_sum / MINHASH_SIZE;
    double mh_var = mh_sq_sum / MINHASH_SIZE - mh_mean * mh_mean;
    values[idx++].d = mh_var;                              /* HashMinVar */

    /* 30: Combined ionizable hash (acidic XOR basic weighted) */
    values[idx++].d = hash_to_double(s.acidic_hash ^ (s.basic_hash >> 1)); /* HashIonizable */

    return idx;
}

/* ============================================================================
 * Individual Descriptor Functions
 * ============================================================================ */

#define HASH_DESC_FN(name, stat_idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    descriptor_value_t vals[NUM_HASH_DESCRIPTORS]; \
    int n = descriptors_compute_hash_all(mol, vals); \
    if (n < 0 || stat_idx >= n) return CCHEM_ERROR_INVALID_INPUT; \
    *value = vals[stat_idx]; \
    return CCHEM_OK; \
}

HASH_DESC_FN(hash_full0, 0)
HASH_DESC_FN(hash_full1, 1)
HASH_DESC_FN(hash_full2, 2)
HASH_DESC_FN(hash_bigram, 3)
HASH_DESC_FN(hash_trigram, 4)
HASH_DESC_FN(hash_quadgram, 5)
HASH_DESC_FN(hash_aromatic, 6)
HASH_DESC_FN(hash_aliphatic, 7)
HASH_DESC_FN(hash_hetero, 8)
HASH_DESC_FN(hash_ring, 9)
HASH_DESC_FN(hash_branch, 10)
HASH_DESC_FN(minhash0, 11)
HASH_DESC_FN(minhash1, 12)
HASH_DESC_FN(minhash2, 13)
HASH_DESC_FN(minhash3, 14)
HASH_DESC_FN(minhash4, 15)
HASH_DESC_FN(minhash5, 16)
HASH_DESC_FN(minhash6, 17)
HASH_DESC_FN(minhash7, 18)
HASH_DESC_FN(hash_window_min, 19)
HASH_DESC_FN(hash_window_max, 20)
HASH_DESC_FN(hash_window_mean, 21)
HASH_DESC_FN(hash_pos_weight, 22)
HASH_DESC_FN(hash_xor_chain, 23)
HASH_DESC_FN(hash_transition, 24)
HASH_DESC_FN(hash_acidic, 25)
HASH_DESC_FN(hash_basic, 26)
HASH_DESC_FN(hash_entropy, 27)
HASH_DESC_FN(hash_min_var, 28)
HASH_DESC_FN(hash_ionizable, 29)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REG_HASH(dname, ddesc, fn) do { \
    memset(&def, 0, sizeof(def)); \
    strncpy(def.name, dname, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, ddesc, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_FINGERPRINT; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = fn; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_hash(void) {
    descriptor_def_t def;

    /* Full SMILES hashes */
    REG_HASH("HashFull0", "xxHash of full SMILES (seed 0)", desc_hash_full0);
    REG_HASH("HashFull1", "xxHash of full SMILES (seed 1)", desc_hash_full1);
    REG_HASH("HashFull2", "xxHash of full SMILES (seed 2)", desc_hash_full2);

    /* N-gram hashes */
    REG_HASH("HashBigram", "Mean hash of character bigrams", desc_hash_bigram);
    REG_HASH("HashTrigram", "Mean hash of character trigrams", desc_hash_trigram);
    REG_HASH("HashQuadgram", "Mean hash of character quadgrams", desc_hash_quadgram);

    /* Pattern hashes */
    REG_HASH("HashAromatic", "Hash of aromatic characters only", desc_hash_aromatic);
    REG_HASH("HashAliphatic", "Hash of aliphatic characters only", desc_hash_aliphatic);
    REG_HASH("HashHetero", "Hash of heteroatom characters", desc_hash_hetero);
    REG_HASH("HashRing", "Hash of ring closure digits", desc_hash_ring);
    REG_HASH("HashBranch", "Hash of branch parentheses", desc_hash_branch);

    /* MinHash signatures */
    REG_HASH("MinHash0", "MinHash signature 0", desc_minhash0);
    REG_HASH("MinHash1", "MinHash signature 1", desc_minhash1);
    REG_HASH("MinHash2", "MinHash signature 2", desc_minhash2);
    REG_HASH("MinHash3", "MinHash signature 3", desc_minhash3);
    REG_HASH("MinHash4", "MinHash signature 4", desc_minhash4);
    REG_HASH("MinHash5", "MinHash signature 5", desc_minhash5);
    REG_HASH("MinHash6", "MinHash signature 6", desc_minhash6);
    REG_HASH("MinHash7", "MinHash signature 7", desc_minhash7);

    /* Window statistics */
    REG_HASH("HashWindowMin", "Minimum sliding window hash", desc_hash_window_min);
    REG_HASH("HashWindowMax", "Maximum sliding window hash", desc_hash_window_max);
    REG_HASH("HashWindowMean", "Mean sliding window hash", desc_hash_window_mean);

    /* Order-sensitive hashes */
    REG_HASH("HashPosWeight", "Position-weighted hash mean", desc_hash_pos_weight);
    REG_HASH("HashXorChain", "XOR chain hash", desc_hash_xor_chain);
    REG_HASH("HashTransition", "Character class transition hash", desc_hash_transition);

    /* Ionizable pattern hashes */
    REG_HASH("HashAcidic", "Hash of acidic group patterns", desc_hash_acidic);
    REG_HASH("HashBasic", "Hash of basic group patterns", desc_hash_basic);

    /* Derived hash statistics */
    REG_HASH("HashEntropy", "Hash entropy proxy (range/mean)", desc_hash_entropy);
    REG_HASH("HashMinVar", "MinHash variance", desc_hash_min_var);
    REG_HASH("HashIonizable", "Combined ionizable hash", desc_hash_ionizable);
}
