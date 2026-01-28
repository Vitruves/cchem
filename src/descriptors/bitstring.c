/**
 * @file bitstring.c
 * @brief Bit-level SMILES string descriptors
 *
 * Ultra-fast descriptors computed directly from the ASCII representation
 * of canonical SMILES strings. These provide unique topological and
 * electronic fingerprints based on character encoding properties.
 *
 * Group A: Electronic Case Bits (aromaticity)
 * Group B: Vertical Bit Integration (columnar)
 * Group C: Horizontal Bit Dynamics (dipoles)
 * Group D: Atomic/Electronic Density (popcounts)
 * Group E: Structural Bit-masks
 * Group F: Information Entropy
 * Group G: Electronegativity Proxies
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/smiles_writer.h"
#include "cchem/utils/simd.h"

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

/* Population count for 8-bit byte - uses CPU POPCNT instruction */
static inline int popcount8(unsigned char b) {
    return __builtin_popcount(b);
}

/* Population count for all bits in a string - uses SIMD when available */
static int popcount_string(const char* s, int len) {
    return simd_popcount_bytes((const uint8_t*)s, len);
}

/* Hamming distance between two bytes */
static inline int hamming_dist(unsigned char a, unsigned char b) {
    return popcount8(a ^ b);
}

/* DJB2 hash function */
static unsigned long djb2_hash(const char* str) {
    unsigned long hash = 5381;
    int c;
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c;
    }
    return hash;
}

/* Check if character is lowercase (aromatic) */
static inline int is_aromatic_char(char c) {
    return (c >= 'a' && c <= 'z');
}

/* Check if character is uppercase (aliphatic) */
static inline int is_aliphatic_char(char c) {
    return (c >= 'A' && c <= 'Z');
}

/* Check if character is a digit */
static inline int is_digit_char(char c) {
    return (c >= '0' && c <= '9');
}

/* Check if character is a structural symbol */
static inline int is_structure_char(char c) {
    return (c == '(' || c == ')' || c == '[' || c == ']' ||
            c == '=' || c == '#' || c == '@' || c == '/' || c == '\\');
}

/* Get SMILES string from molecule */
static char* get_smiles(const molecule_t* mol) {
    return molecule_to_canonical_smiles_str(mol);
}

/* ============================================================================
 * Batch Computation Cache (thread-safe: stack-allocated per call)
 * ============================================================================ */

typedef struct {
    /* Pre-computed from single SMILES pass */
    int len;
    int lower_count;
    int upper_count;
    int digit_count;
    int structure_count;
    int case_alternations;
    int bit5_sum;
    int col_sums[8];
    int high_nibble_sum;
    int low_nibble_sum;
    int hamming_sum;
    int xor_result;
    int bit_transitions;
    int total_popcount;
    double popcount_sum;
    double popcount_sum_sq;
    int max_branch_depth;
    int char_counts[256];
    long char_sum;
    long hetero_char_sum;
    int carbon_count;
    int eneg_popcount;
    double weighted_elem_sum;
    int unique_chars;
    unsigned long djb2_hash;
    bool computed;
} bitstring_cache_t;

/* Compute all bitstring statistics in a single pass */
static void compute_bitstring_cache(const char* smiles, int len, bitstring_cache_t* cache) {
    memset(cache, 0, sizeof(bitstring_cache_t));
    if (!smiles || len == 0) {
        cache->computed = true;
        return;
    }

    cache->len = len;
    cache->djb2_hash = 5381;

    int prev_case = -1;
    int depth = 0;
    unsigned char prev_byte = 0;

    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)smiles[i];

        /* Character classification */
        if (c >= 'a' && c <= 'z') {
            cache->lower_count++;
            if (prev_case == 1) cache->case_alternations++;
            prev_case = 0;
        } else if (c >= 'A' && c <= 'Z') {
            cache->upper_count++;
            if (prev_case == 0) cache->case_alternations++;
            prev_case = 1;
        }

        if (c >= '0' && c <= '9') cache->digit_count++;
        if (c == '(' || c == ')' || c == '[' || c == ']' ||
            c == '=' || c == '#' || c == '@' || c == '/' || c == '\\') {
            cache->structure_count++;
        }

        /* Branch depth */
        if (c == '(') {
            depth++;
            if (depth > cache->max_branch_depth) cache->max_branch_depth = depth;
        } else if (c == ')') {
            depth--;
        }

        /* Bit operations */
        cache->bit5_sum += (c >> 5) & 1;
        for (int b = 0; b < 8; b++) {
            cache->col_sums[b] += (c >> b) & 1;
        }
        cache->high_nibble_sum += (c >> 4) & 0x0F;
        cache->low_nibble_sum += c & 0x0F;

        /* Popcount */
        int pc = popcount8(c);
        cache->total_popcount += pc;
        cache->popcount_sum += pc;
        cache->popcount_sum_sq += pc * pc;

        /* XOR chain */
        cache->xor_result ^= c;

        /* Hamming and transitions */
        if (i > 0) {
            cache->hamming_sum += popcount8(prev_byte ^ c);
            for (int b = 0; b < 8; b++) {
                int prev_bit = (prev_byte >> b) & 1;
                int curr_bit = (c >> b) & 1;
                if (prev_bit != curr_bit) cache->bit_transitions++;
            }
        }
        prev_byte = c;

        /* Character counts for entropy */
        if (cache->char_counts[c] == 0) cache->unique_chars++;
        cache->char_counts[c]++;

        /* Character sum */
        cache->char_sum += c;

        /* DJB2 hash */
        cache->djb2_hash = ((cache->djb2_hash << 5) + cache->djb2_hash) + c;

        /* Heteroatom and element analysis */
        char ch = smiles[i];
        if (ch == 'N' || ch == 'n' || ch == 'O' || ch == 'o' ||
            ch == 'S' || ch == 's' || ch == 'P' || ch == 'p' ||
            ch == 'F' || ch == 'I') {
            cache->hetero_char_sum += c;
        }
        if (i < len - 1) {
            if ((ch == 'C' && smiles[i+1] == 'l') ||
                (ch == 'B' && smiles[i+1] == 'r')) {
                cache->hetero_char_sum += c + (unsigned char)smiles[i+1];
            }
        }

        if (ch == 'C' || ch == 'c') cache->carbon_count++;

        /* Electronegative popcount */
        if (ch == 'N' || ch == 'n' || ch == 'O' || ch == 'o' ||
            ch == 'F' || ch == 'S' || ch == 's') {
            cache->eneg_popcount += pc;
        }

        /* Weighted element sum */
        if (ch == 'C' || ch == 'c') cache->weighted_elem_sum += 2.55;
        else if (ch == 'N' || ch == 'n') cache->weighted_elem_sum += 3.04;
        else if (ch == 'O' || ch == 'o') cache->weighted_elem_sum += 3.44;
        else if (ch == 'F') cache->weighted_elem_sum += 3.98;
        else if (ch == 'S' || ch == 's') cache->weighted_elem_sum += 2.58;
        else if (ch == 'P' || ch == 'p') cache->weighted_elem_sum += 2.19;
    }

    cache->computed = true;
}

/* Number of bitstring descriptors */
#define NUM_BITSTRING_DESCRIPTORS 30

/* Batch compute ALL bitstring descriptors - thread-safe (stack allocated cache) */
int descriptors_compute_bitstring_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    char* smiles = get_smiles(mol);
    if (!smiles) {
        for (int i = 0; i < NUM_BITSTRING_DESCRIPTORS; i++) values[i].d = 0.0;
        return NUM_BITSTRING_DESCRIPTORS;
    }

    int len = (int)strlen(smiles);
    bitstring_cache_t cache;
    compute_bitstring_cache(smiles, len, &cache);

    int idx = 0;

    /* Group A: Electronic Case Bits (4) */
    /* 1. LowercaseRatio */
    values[idx++].d = (len > 0) ? (double)cache.lower_count / len : 0.0;
    /* 2. CaseAlternation */
    values[idx++].d = (double)cache.case_alternations;
    /* 3. Bit5Sum */
    values[idx++].d = (double)cache.bit5_sum;
    /* 4. CaseEntropy */
    {
        int total = cache.lower_count + cache.upper_count;
        if (total > 0) {
            double p_low = (double)cache.lower_count / total;
            double p_up = (double)cache.upper_count / total;
            double H = 0.0;
            if (p_low > 0) H -= p_low * log2(p_low);
            if (p_up > 0) H -= p_up * log2(p_up);
            values[idx++].d = H;
        } else {
            values[idx++].d = 0.0;
        }
    }

    /* Group B: Vertical Bit Integration (5) */
    /* 5. BitColumnSum */
    {
        double total = 0;
        for (int b = 0; b < 8; b++) total += cache.col_sums[b];
        values[idx++].d = total;
    }
    /* 6. BitColumnVariance */
    {
        double mean = 0;
        for (int b = 0; b < 8; b++) mean += cache.col_sums[b];
        mean /= 8.0;
        double var = 0;
        for (int b = 0; b < 8; b++) {
            double diff = cache.col_sums[b] - mean;
            var += diff * diff;
        }
        values[idx++].d = var / 8.0;
    }
    /* 7. HighNibbleSum */
    values[idx++].d = (double)cache.high_nibble_sum;
    /* 8. LowNibbleSum */
    values[idx++].d = (double)cache.low_nibble_sum;
    /* 9. NibbleRatio */
    values[idx++].d = (cache.low_nibble_sum > 0) ? (double)cache.high_nibble_sum / cache.low_nibble_sum : 0.0;

    /* Group C: Horizontal Bit Dynamics (4) */
    /* 10. HammingDistanceSum */
    values[idx++].d = (double)cache.hamming_sum;
    /* 11. MeanHammingDistance */
    values[idx++].d = (len > 1) ? (double)cache.hamming_sum / (len - 1) : 0.0;
    /* 12. XORChainValue */
    values[idx++].d = (double)cache.xor_result;
    /* 13. BitTransitions */
    values[idx++].d = (double)cache.bit_transitions;

    /* Group D: Atomic/Electronic Density (4) */
    /* 14. TotalPopcount */
    values[idx++].d = (double)cache.total_popcount;
    /* 15. MeanPopcount */
    values[idx++].d = (len > 0) ? cache.popcount_sum / len : 0.0;
    /* 16. PopcountVariance */
    {
        if (len > 0) {
            double mean = cache.popcount_sum / len;
            values[idx++].d = (cache.popcount_sum_sq / len) - (mean * mean);
        } else {
            values[idx++].d = 0.0;
        }
    }
    /* 17. BitDensity */
    values[idx++].d = (len > 0) ? (double)cache.total_popcount / (8.0 * len) : 0.0;

    /* Group E: Structural Bit-masks (3) */
    /* 18. StructureCharCount */
    values[idx++].d = (double)cache.structure_count;
    /* 19. DigitCount */
    values[idx++].d = (double)cache.digit_count;
    /* 20. BranchDepth */
    values[idx++].d = (double)cache.max_branch_depth;

    /* Group F: Information Entropy (5) */
    /* 21. ShannonEntropy */
    {
        double H = 0.0;
        if (len > 0) {
            for (int i = 0; i < 256; i++) {
                if (cache.char_counts[i] > 0) {
                    double p = (double)cache.char_counts[i] / len;
                    H -= p * log2(p);
                }
            }
        }
        values[idx++].d = H;
    }
    /* 22. UniqueCharCount */
    values[idx++].d = (double)cache.unique_chars;
    /* 23. StringLength */
    values[idx++].d = (double)len;
    /* 24. DJB2Hash */
    values[idx++].d = (double)(cache.djb2_hash % 10000);
    /* 25. CharacterSum */
    values[idx++].d = (double)cache.char_sum;

    /* Group G: Electronegativity Proxies (5) */
    /* 26. HeteroCharSum */
    values[idx++].d = (double)cache.hetero_char_sum;
    /* 27. CarbonCharRatio */
    values[idx++].d = (len > 0) ? (double)cache.carbon_count / len : 0.0;
    /* 28. ENegPopcount */
    values[idx++].d = (double)cache.eneg_popcount;
    /* 29. WeightedElemSum */
    values[idx++].d = cache.weighted_elem_sum;
    /* 30. MeanCharValue */
    values[idx++].d = (len > 0) ? (double)cache.char_sum / len : 0.0;

    free(smiles);
    return NUM_BITSTRING_DESCRIPTORS;
}

/* ============================================================================
 * Group A: Electronic Case Bits (4 descriptors)
 * ============================================================================ */

/* 1. Lowercase Ratio (aromaticity indicator) */
static cchem_status_t desc_lowercase_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int lower = 0;
    for (int i = 0; i < len; i++) {
        if (is_aromatic_char(smiles[i])) lower++;
    }

    value->d = (len > 0) ? (double)lower / len : 0.0;
    free(smiles);
    return CCHEM_OK;
}

/* 2. Case Alternation Count */
static cchem_status_t desc_case_alternation(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int alt = 0;
    int prev_case = -1; /* -1=none, 0=lower, 1=upper */

    for (int i = 0; i < len; i++) {
        int curr_case = -1;
        if (is_aromatic_char(smiles[i])) curr_case = 0;
        else if (is_aliphatic_char(smiles[i])) curr_case = 1;

        if (curr_case >= 0 && prev_case >= 0 && curr_case != prev_case) {
            alt++;
        }
        if (curr_case >= 0) prev_case = curr_case;
    }

    value->d = (double)alt;
    free(smiles);
    return CCHEM_OK;
}

/* 3. Bit 5 Sum (case bit in ASCII) */
static cchem_status_t desc_bit5_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int sum = 0;
    for (int i = 0; i < len; i++) {
        sum += (smiles[i] >> 5) & 1;
    }

    value->d = (double)sum;
    free(smiles);
    return CCHEM_OK;
}

/* 4. Case Entropy */
static cchem_status_t desc_case_entropy(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int lower = 0, upper = 0;

    for (int i = 0; i < len; i++) {
        if (is_aromatic_char(smiles[i])) lower++;
        else if (is_aliphatic_char(smiles[i])) upper++;
    }

    int total = lower + upper;
    if (total > 0) {
        double p_low = (double)lower / total;
        double p_up = (double)upper / total;
        double H = 0.0;
        if (p_low > 0) H -= p_low * log2(p_low);
        if (p_up > 0) H -= p_up * log2(p_up);
        value->d = H;
    } else {
        value->d = 0.0;
    }

    free(smiles);
    return CCHEM_OK;
}

/* ============================================================================
 * Group B: Vertical Bit Integration (5 descriptors)
 * ============================================================================ */

/* 5. Bit Column Sum (sum of popcount for each bit position) */
static cchem_status_t desc_bit_column_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int col_sums[8] = {0};

    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)smiles[i];
        for (int b = 0; b < 8; b++) {
            col_sums[b] += (c >> b) & 1;
        }
    }

    double total = 0;
    for (int b = 0; b < 8; b++) {
        total += col_sums[b];
    }

    value->d = total;
    free(smiles);
    return CCHEM_OK;
}

/* 6. Bit Column Variance */
static cchem_status_t desc_bit_column_variance(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int col_sums[8] = {0};

    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)smiles[i];
        for (int b = 0; b < 8; b++) {
            col_sums[b] += (c >> b) & 1;
        }
    }

    double mean = 0;
    for (int b = 0; b < 8; b++) mean += col_sums[b];
    mean /= 8.0;

    double var = 0;
    for (int b = 0; b < 8; b++) {
        double diff = col_sums[b] - mean;
        var += diff * diff;
    }
    var /= 8.0;

    value->d = var;
    free(smiles);
    return CCHEM_OK;
}

/* 7. High Nibble Sum */
static cchem_status_t desc_high_nibble_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int sum = 0;
    for (int i = 0; i < len; i++) {
        sum += (smiles[i] >> 4) & 0x0F;
    }

    value->d = (double)sum;
    free(smiles);
    return CCHEM_OK;
}

/* 8. Low Nibble Sum */
static cchem_status_t desc_low_nibble_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int sum = 0;
    for (int i = 0; i < len; i++) {
        sum += smiles[i] & 0x0F;
    }

    value->d = (double)sum;
    free(smiles);
    return CCHEM_OK;
}

/* 9. Nibble Ratio */
static cchem_status_t desc_nibble_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int high = 0, low = 0;
    for (int i = 0; i < len; i++) {
        high += (smiles[i] >> 4) & 0x0F;
        low += smiles[i] & 0x0F;
    }

    value->d = (low > 0) ? (double)high / low : 0.0;
    free(smiles);
    return CCHEM_OK;
}

/* ============================================================================
 * Group C: Horizontal Bit Dynamics (4 descriptors)
 * ============================================================================ */

/* 10. Hamming Distance Sum (adjacent chars) */
static cchem_status_t desc_hamming_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int sum = 0;
    for (int i = 1; i < len; i++) {
        sum += hamming_dist((unsigned char)smiles[i-1], (unsigned char)smiles[i]);
    }

    value->d = (double)sum;
    free(smiles);
    return CCHEM_OK;
}

/* 11. Mean Hamming Distance */
static cchem_status_t desc_mean_hamming(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    if (len < 2) { value->d = 0.0; free(smiles); return CCHEM_OK; }

    int sum = 0;
    for (int i = 1; i < len; i++) {
        sum += hamming_dist((unsigned char)smiles[i-1], (unsigned char)smiles[i]);
    }

    value->d = (double)sum / (len - 1);
    free(smiles);
    return CCHEM_OK;
}

/* 12. XOR Chain Value */
static cchem_status_t desc_xor_chain(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    unsigned char result = 0;
    for (int i = 0; i < len; i++) {
        result ^= (unsigned char)smiles[i];
    }

    value->d = (double)result;
    free(smiles);
    return CCHEM_OK;
}

/* 13. Bit Transition Count */
static cchem_status_t desc_bit_transitions(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int trans = 0;
    unsigned char prev = 0;

    for (int i = 0; i < len; i++) {
        unsigned char curr = (unsigned char)smiles[i];
        for (int b = 0; b < 8; b++) {
            int prev_bit = (prev >> b) & 1;
            int curr_bit = (curr >> b) & 1;
            if (prev_bit != curr_bit) trans++;
        }
        prev = curr;
    }

    value->d = (double)trans;
    free(smiles);
    return CCHEM_OK;
}

/* ============================================================================
 * Group D: Atomic/Electronic Density (4 descriptors)
 * ============================================================================ */

/* 14. Total Popcount */
static cchem_status_t desc_total_popcount(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    value->d = (double)popcount_string(smiles, len);
    free(smiles);
    return CCHEM_OK;
}

/* 15. Mean Popcount per Character */
static cchem_status_t desc_mean_popcount(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int total = popcount_string(smiles, len);

    value->d = (len > 0) ? (double)total / len : 0.0;
    free(smiles);
    return CCHEM_OK;
}

/* 16. Popcount Variance */
static cchem_status_t desc_popcount_variance(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    if (len == 0) { value->d = 0.0; free(smiles); return CCHEM_OK; }

    double sum = 0, sum_sq = 0;
    for (int i = 0; i < len; i++) {
        int pc = popcount8((unsigned char)smiles[i]);
        sum += pc;
        sum_sq += pc * pc;
    }

    double mean = sum / len;
    value->d = (sum_sq / len) - (mean * mean);
    free(smiles);
    return CCHEM_OK;
}

/* 17. Bit Density (popcount / (8 * len)) */
static cchem_status_t desc_bit_density(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int total = popcount_string(smiles, len);

    value->d = (len > 0) ? (double)total / (8.0 * len) : 0.0;
    free(smiles);
    return CCHEM_OK;
}

/* ============================================================================
 * Group E: Structural Bit-masks (3 descriptors)
 * ============================================================================ */

/* 18. Structure Character Count */
static cchem_status_t desc_structure_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int count = 0;
    for (int i = 0; i < len; i++) {
        if (is_structure_char(smiles[i])) count++;
    }

    value->d = (double)count;
    free(smiles);
    return CCHEM_OK;
}

/* 19. Digit Count (ring closures) */
static cchem_status_t desc_digit_count(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int count = 0;
    for (int i = 0; i < len; i++) {
        if (is_digit_char(smiles[i])) count++;
    }

    value->d = (double)count;
    free(smiles);
    return CCHEM_OK;
}

/* 20. Branch Depth (max parenthesis nesting) */
static cchem_status_t desc_branch_depth(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int depth = 0, max_depth = 0;
    for (int i = 0; i < len; i++) {
        if (smiles[i] == '(') {
            depth++;
            if (depth > max_depth) max_depth = depth;
        } else if (smiles[i] == ')') {
            depth--;
        }
    }

    value->d = (double)max_depth;
    free(smiles);
    return CCHEM_OK;
}

/* ============================================================================
 * Group F: Information Entropy (5 descriptors)
 * ============================================================================ */

/* 21. Shannon Entropy (character distribution) */
static cchem_status_t desc_shannon_entropy(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    if (len == 0) { value->d = 0.0; free(smiles); return CCHEM_OK; }

    int counts[256] = {0};
    for (int i = 0; i < len; i++) {
        counts[(unsigned char)smiles[i]]++;
    }

    double H = 0.0;
    for (int i = 0; i < 256; i++) {
        if (counts[i] > 0) {
            double p = (double)counts[i] / len;
            H -= p * log2(p);
        }
    }

    value->d = H;
    free(smiles);
    return CCHEM_OK;
}

/* 22. Unique Character Count */
static cchem_status_t desc_unique_chars(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int seen[256] = {0};
    int unique = 0;

    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)smiles[i];
        if (!seen[c]) {
            seen[c] = 1;
            unique++;
        }
    }

    value->d = (double)unique;
    free(smiles);
    return CCHEM_OK;
}

/* 23. String Length */
static cchem_status_t desc_string_length(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    value->d = (double)strlen(smiles);
    free(smiles);
    return CCHEM_OK;
}

/* 24. DJB2 Hash (mod 10000) */
static cchem_status_t desc_djb2_hash(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    unsigned long hash = djb2_hash(smiles);
    value->d = (double)(hash % 10000);
    free(smiles);
    return CCHEM_OK;
}

/* 25. Character Sum */
static cchem_status_t desc_char_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    long sum = 0;
    for (int i = 0; i < len; i++) {
        sum += (unsigned char)smiles[i];
    }

    value->d = (double)sum;
    free(smiles);
    return CCHEM_OK;
}

/* ============================================================================
 * Group G: Electronegativity Proxies (5 descriptors)
 * ============================================================================ */

/* 26. Heteroatom Character Sum (N, O, S, P, F, Cl, Br, I chars) */
static cchem_status_t desc_hetero_char_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    long sum = 0;
    for (int i = 0; i < len; i++) {
        char c = smiles[i];
        if (c == 'N' || c == 'n' || c == 'O' || c == 'o' ||
            c == 'S' || c == 's' || c == 'P' || c == 'p' ||
            c == 'F' || c == 'I') {
            sum += (unsigned char)c;
        }
        /* Handle Cl, Br */
        if (i < len - 1) {
            if ((c == 'C' && smiles[i+1] == 'l') ||
                (c == 'B' && smiles[i+1] == 'r')) {
                sum += (unsigned char)c + (unsigned char)smiles[i+1];
            }
        }
    }

    value->d = (double)sum;
    free(smiles);
    return CCHEM_OK;
}

/* 27. Carbon Character Ratio */
static cchem_status_t desc_carbon_char_ratio(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int c_count = 0;
    for (int i = 0; i < len; i++) {
        if (smiles[i] == 'C' || smiles[i] == 'c') c_count++;
    }

    value->d = (len > 0) ? (double)c_count / len : 0.0;
    free(smiles);
    return CCHEM_OK;
}

/* 28. Electronegative Character Popcount */
static cchem_status_t desc_eneg_popcount(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    int sum = 0;
    for (int i = 0; i < len; i++) {
        char c = smiles[i];
        if (c == 'N' || c == 'n' || c == 'O' || c == 'o' ||
            c == 'F' || c == 'S' || c == 's') {
            sum += popcount8((unsigned char)c);
        }
    }

    value->d = (double)sum;
    free(smiles);
    return CCHEM_OK;
}

/* 29. Weighted Element Sum */
static cchem_status_t desc_weighted_elem_sum(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    double sum = 0.0;

    /* Weight by approximate electronegativity */
    for (int i = 0; i < len; i++) {
        char c = smiles[i];
        if (c == 'C' || c == 'c') sum += 2.55;
        else if (c == 'N' || c == 'n') sum += 3.04;
        else if (c == 'O' || c == 'o') sum += 3.44;
        else if (c == 'F') sum += 3.98;
        else if (c == 'S' || c == 's') sum += 2.58;
        else if (c == 'P' || c == 'p') sum += 2.19;
    }

    value->d = sum;
    free(smiles);
    return CCHEM_OK;
}

/* 30. Mean Character Value */
static cchem_status_t desc_mean_char_value(const molecule_t* mol, descriptor_value_t* value) {
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT;

    char* smiles = get_smiles(mol);
    if (!smiles) { value->d = 0.0; return CCHEM_OK; }

    int len = (int)strlen(smiles);
    if (len == 0) { value->d = 0.0; free(smiles); return CCHEM_OK; }

    long sum = 0;
    for (int i = 0; i < len; i++) {
        sum += (unsigned char)smiles[i];
    }

    value->d = (double)sum / len;
    free(smiles);
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REG(n, d, f) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, n, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, d, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_BITSTRING; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = f; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_bitstring(void) {
    /* Group A: Electronic Case Bits */
    REG("LowercaseRatio", "Ratio of lowercase (aromatic) chars", desc_lowercase_ratio);
    REG("CaseAlternation", "Case transitions in SMILES", desc_case_alternation);
    REG("Bit5Sum", "Sum of bit 5 (case bit) values", desc_bit5_sum);
    REG("CaseEntropy", "Shannon entropy of case distribution", desc_case_entropy);

    /* Group B: Vertical Bit Integration */
    REG("BitColumnSum", "Sum of all bit positions", desc_bit_column_sum);
    REG("BitColumnVariance", "Variance across bit columns", desc_bit_column_variance);
    REG("HighNibbleSum", "Sum of high nibbles", desc_high_nibble_sum);
    REG("LowNibbleSum", "Sum of low nibbles", desc_low_nibble_sum);
    REG("NibbleRatio", "High nibble / low nibble ratio", desc_nibble_ratio);

    /* Group C: Horizontal Bit Dynamics */
    REG("HammingDistanceSum", "Sum of adjacent Hamming distances", desc_hamming_sum);
    REG("MeanHammingDistance", "Mean adjacent Hamming distance", desc_mean_hamming);
    REG("XORChainValue", "XOR of all characters", desc_xor_chain);
    REG("BitTransitions", "Total bit transitions", desc_bit_transitions);

    /* Group D: Atomic/Electronic Density */
    REG("TotalPopcount", "Total popcount of SMILES", desc_total_popcount);
    REG("MeanPopcount", "Mean popcount per character", desc_mean_popcount);
    REG("PopcountVariance", "Variance of popcounts", desc_popcount_variance);
    REG("BitDensity", "Bit density (1s / total bits)", desc_bit_density);

    /* Group E: Structural Bit-masks */
    REG("StructureCharCount", "Count of structural symbols", desc_structure_count);
    REG("DigitCount", "Count of digits (ring closures)", desc_digit_count);
    REG("BranchDepth", "Maximum parenthesis nesting", desc_branch_depth);

    /* Group F: Information Entropy */
    REG("ShannonEntropy", "Character distribution entropy", desc_shannon_entropy);
    REG("UniqueCharCount", "Number of unique characters", desc_unique_chars);
    REG("StringLength", "SMILES string length", desc_string_length);
    REG("DJB2Hash", "DJB2 hash (mod 10000)", desc_djb2_hash);
    REG("CharacterSum", "Sum of character ASCII values", desc_char_sum);

    /* Group G: Electronegativity Proxies */
    REG("HeteroCharSum", "Sum of heteroatom char values", desc_hetero_char_sum);
    REG("CarbonCharRatio", "Carbon character ratio", desc_carbon_char_ratio);
    REG("ENegPopcount", "Electronegative char popcount", desc_eneg_popcount);
    REG("WeightedElemSum", "Chi-weighted element sum", desc_weighted_elem_sum);
    REG("MeanCharValue", "Mean character ASCII value", desc_mean_char_value);
}
