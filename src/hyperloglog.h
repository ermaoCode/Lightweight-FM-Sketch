#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdint.h>

#include "murmurhash.h"

typedef struct hll_cnt_ctx_s {
    int err;
    uint8_t log2m;
    uint32_t m;
    double alphaMM;
    uint8_t hf;
    uint8_t M[1];
}hll_cnt_ctx_t;
/**
 * Predefined error codes
 * */
enum {
    CCARD_OK = 0,                   /**< No error */
    CCARD_ERR_INVALID_CTX = -1,     /**< Invalid algorihm context */
    CCARD_ERR_MERGE_FAILED = -2,    /**< Merge failed */
    CCARD_ERR_INVALID_ARGUMENT = -3,    /**< Invalid argument */
    CCARD_ERR_PLACEHOLDER
};

/**
 * Algorithms
 * */
enum {
    CCARD_ALGO_ADAPTIVE = 1,
    CCARD_ALGO_HYPERLOGLOG = 2,
    CCARD_ALGO_LINEAR = 3,
    CCARD_ALGO_HYPERLOGLOGPLUS = 4,
    CCARD_ALGO_PLACEHOLDER
};

/**
 * Hash functions
 * */
enum {
    CCARD_HASH_MURMUR = 1,
    CCARD_HASH_LOOKUP3 = 2,
    CCARD_HASH_MURMUR64 = 3,
	CCARD_HASH_MURMUR3 = 4,
    CCARD_HASH_PLACEHOLDER
};

class Hyperloglog
{
public:
    Hyperloglog();
    int64_t hll_cnt_card();
    int hll_cnt_offer(const void *buf, uint32_t len);
    int hll_cnt_fini();
private:
    uint8_t num_of_trail_zeros(uint64_t i);
	
	//initialization
    hll_cnt_ctx_t *hll_cnt_raw_init(const void *obuf, uint32_t len_or_k, uint8_t hf);
    hll_cnt_ctx_t *hll_cnt_init(const void *obuf, uint32_t len_or_k, uint8_t hf);

    int hll_cnt_get_raw_bytes( void *buf, uint32_t *len);
    int hll_cnt_get_bytes( void *buf, uint32_t *len);
    int hll_cnt_reset();
    int hll_cnt_errnum(hll_cnt_ctx_t *ctx);
    const char *hll_cnt_errstr(int err);

    hll_cnt_ctx_t *ctx;



/*
    static const double POW_2_32 = 4294967296.0;
    static const double NEGATIVE_POW_2_32 = -4294967296.0;*/

	double POW_2_32;
    double NEGATIVE_POW_2_32;

    uint32_t murmurhash(void *buf, uint32_t len, uint32_t seed);
    uint32_t murmurhash_long(uint64_t data);
    uint64_t murmurhash64(void *buf, uint32_t len, uint32_t seed);
    uint64_t murmurhash64_no_seed(void *buf, uint32_t len);

    uint32_t lookup3(const uint32_t *k, uint32_t offset, uint32_t length, uint32_t initval);
    uint32_t lookup3ycs(const uint32_t *k, uint32_t offset, uint32_t length, uint32_t initval);
    uint32_t lookup3ycs_2(const char *s, uint32_t start, uint32_t end, uint32_t initval);
    uint64_t lookup3ycs64(const char *s, uint32_t start, uint32_t end, uint64_t initval);
    uint64_t lookup3ycs64_2(const char *s);

};

#endif // HYPERLOGLOG_H
