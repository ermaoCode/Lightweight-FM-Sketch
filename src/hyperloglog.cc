#include "hyperloglog.h"
#include <math.h>

Hyperloglog::Hyperloglog()
{
    ctx=this->hll_cnt_init (NULL, 1024, CCARD_HASH_MURMUR3);
	POW_2_32 = 4294967296.0;
    NEGATIVE_POW_2_32 = -4294967296.0;
}

// caculate the number of 0s starting from the right
uint8_t Hyperloglog::num_of_trail_zeros(uint64_t i)
{
    uint64_t y;
    uint8_t n = 63;

    if (i == 0)
        return 64;

    y = i << 32;    if (y != 0) { n -= 32; i = y; }
    y = i << 16;    if (y != 0) { n -= 16; i = y; }
    y = i << 8;     if (y != 0) { n -= 8; i = y; }
    y = i << 4;     if (y != 0) { n -= 4; i = y; }
    y = i << 2;     if (y != 0) { n -= 2; i = y; }

    return n - (uint8_t)((i << 1) >> 63);
}

// obuf: M[]
// len_or_k: length of M[], if buf is null, 1 << len_or_k  will be the length
hll_cnt_ctx_t * Hyperloglog::hll_cnt_raw_init(const void *obuf, uint32_t len_or_k, uint8_t hf)
{
    /*
     +--------------+---------+------------------------------+-----------+
     | algorithm[1] | hash[1] | bitmap length(base-2 log)[1] | bitmap[n] |
     +--------------+---------+------------------------------+-----------+
     */
//    struct hll_cnt_ctx_s {
//        int err;
//        uint8_t log2m;
//        uint32_t m;
//        double alphaMM;
//        uint8_t hf;
//        uint8_t M[1];
//    };
    hll_cnt_ctx_t *ctx;
    uint8_t *buf = (uint8_t *)obuf;
    uint8_t log2m = buf ? num_of_trail_zeros(len_or_k) : len_or_k;//obuf is NULL?
    //uint32_t m = pow(2, (int)log2m);
	uint32_t m = (uint32_t)1 << log2m;

    if (len_or_k == 0) {
        // invalid buffer length or k
        return NULL;
    }

    if (buf) {
        // initial bitmap was given
        if (len_or_k != (uint32_t)(1 << log2m)) {
            // invalid buffer size, its length must be a power of 2
            return NULL;
        }
		// m = pow((float)2, (int)log2m) = 1 << log2m = len_or_k
		// m: length of M[] in ctx
        ctx = (hll_cnt_ctx_t *)malloc(sizeof(hll_cnt_ctx_t) + m - 1);
        memcpy(ctx->M, buf, m);
    } else {
        // k was given
        ctx = (hll_cnt_ctx_t *)malloc(sizeof(hll_cnt_ctx_t) + m - 1);
        memset(ctx->M, 0, m);
    }
    ctx->err = CCARD_OK;
    ctx->log2m = log2m;
    ctx->m = m;
    ctx->hf = hf;

    /*
     * Description of the following magical numbers:
     *
     * In the HyperLogLog paper page 12-13, alphaMM is a_m*m^2, where:
     *
     *  a_m := 1/(m*J_0(m))
     *
     * Here J_s(m) is not the first-kind Bessel function, but defined as the
     * value of a special integrals:
     *
     *  J_s(m) := integral(u^s*f(u)^m, u=0..inf)
     *
     * where f(u) := log_2((2+u)/(1+u))
     *
     * After some deductions, we know that J_0(m) can be estimated by:
     *
     *  J_0(m) ~= 2*ln(2)/m*(1+(3*ln(2)-1)/m)
     *
     * As 1/(2*ln(2)) ~= 0.72135, 3*ln(2)-1 ~= 1.0794, thus:
     *
     *  a_m ~= 0.72135/(1+1.0794/m)
     *
     * When log_2(m)={4,5,6}, the corresponding a_m will be:
     *
     *  a_16 ~= 0.72135/(1+1.0794/16) = 0.67576
     *  a_32 ~= 0.72135/(1+1.0794/32) = 0.69781
     *  a_64 ~= 0.72135/(1+1.0794/64) = 0.70939
     *
     * There're small errors between calculated and actually used values,
     * because stream-lib copied those values from the pseudo code in page 14
     * directly. We had to keep compatibility with stream-lib so can't correct
     * these values.
     **/
    switch (log2m) {
        case 4:
            ctx->alphaMM = 0.673 * m * m;
            break;
        case 5:
            ctx->alphaMM = 0.697 * m * m;
            break;
        case 6:
            ctx->alphaMM = 0.709 * m * m;
            break;
        default:
            ctx->alphaMM = (0.7213 / (1 + 1.079 / m)) * m * m;
    }

    return ctx;
}

hll_cnt_ctx_t *Hyperloglog::hll_cnt_init(const void *obuf, uint32_t len_or_k, uint8_t hf)
{
    uint8_t *buf = (uint8_t *)obuf;

    if (buf) {
        // initial bitmap was given
		// log2m = num_of_trail_zeros(len_or_k);
        uint8_t log2m = buf ? num_of_trail_zeros(len_or_k) : len_or_k;

        if (buf[0] != CCARD_ALGO_HYPERLOGLOG ||
            buf[1] != hf ||
            buf[2] != log2m) {

            // counting algorithm, hash function or length not match
            return NULL;
        }

        return hll_cnt_raw_init(buf + 3, len_or_k, hf);
    }

    return hll_cnt_raw_init(NULL, len_or_k, hf);
}

int64_t Hyperloglog::hll_cnt_card()
{
    double sum = 0, estimate, zeros = 0;
    uint32_t j, z;

    if (!ctx) {
        return -1;
    }
    ctx->err = CCARD_OK;

    for (j = 0; j < ctx->m; j++) {
        sum += pow((float)2, (-1 * ctx->M[j]));
    }

    estimate = ctx->alphaMM * (1 / sum);

    if (estimate <= (5.0 / 2.0) * ctx->m) {
        /*
         * Small range correction:
         * Empty buckets may be too many, using linear counting estimator
         * instead. 
         * */
        for (z = 0; z < ctx->m; z++) {
            if (ctx->M[z] == 0) {
                zeros++;
            }
        }
        return (int64_t)floor(ctx->m * log(ctx->m / zeros) +0.5);
    } else if (estimate <= (1.0 / 30.0) * POW_2_32) {
        /* Intermediate range - no correction */
        return (int64_t)floor(estimate +0.5);
    } else {
        /* Large range correction */
        return (int64_t)floor((NEGATIVE_POW_2_32 * log(1.0 - (estimate / POW_2_32))) +0.5);
    }
}

// initialize the register array M[j] with one hash value,
// return a bool value that represent whether the M[j] has
// been changed
int Hyperloglog::hll_cnt_offer(const void *buf, uint32_t len)
{
    int modified = 0;
    uint64_t x, j;
    uint8_t r, hl;

    if (!ctx) {
        return -1;
    }
    switch (ctx->hf) {
		case CCARD_HASH_MURMUR3:
			x = (uint64_t)murmurhash3((const char*)buf, (uint32_t) strlen((const char *)buf), -1); // 0xb6d99cf8
			hl = 32;
			break;
        case CCARD_HASH_MURMUR:
            x = (uint64_t)murmurhash((void *)buf, len, -1);
            hl = 32;
            break;
        case CCARD_HASH_LOOKUP3:
            x = lookup3ycs64_2((const char *)buf);
            hl = 64;
            break;
        case CCARD_HASH_MURMUR64:
            x = (uint64_t)murmurhash64_no_seed((void *)buf, len);
            hl = 64;
            break;
        default:
            /* default to use murmurhash function */
            x = (uint64_t)murmurhash((void *)buf, len, -1);
            hl = 32;
    }

    j = x >> (hl - ctx->log2m);
    r = (uint8_t)(num_of_trail_zeros(x << (ctx->log2m + 64 - hl)) - (ctx->log2m + 64 - hl) + 1);
    if (ctx->M[j] < r) {
        ctx->M[j] = r;

        modified = 1;
    }

    ctx->err = CCARD_OK;
    return modified;
}

// ?
int Hyperloglog::hll_cnt_get_raw_bytes(void *buf, uint32_t *len)
{
    uint8_t *out = (uint8_t *)buf;

    if (!ctx || !len || (buf && *len < ctx->m)) {
        return -1;
    }

    if (out) {
        memcpy(out, ctx->M, ctx->m);
    }
    *len = ctx->m;

    return 0;
}

int Hyperloglog::hll_cnt_get_bytes(void *buf, uint32_t *len)
{
    /*
     +--------------+---------+------------------------------+-----------+
     | algorithm[1] | hash[1] | bitmap length(base-2 log)[1] | bitmap[n] |
     +--------------+---------+------------------------------+-----------+
     */
    uint8_t algo = CCARD_ALGO_HYPERLOGLOG;
    uint8_t *out = (uint8_t *)buf;

    if (!ctx || !len || (buf && *len < ctx->m + 3)) {
        return -1;
    }

    if (buf) {
        out[0] = algo;
        out[1] = ctx->hf;
        out[2] = ctx->log2m;
        memcpy(&out[3], ctx->M, ctx->m);
    }
    *len = ctx->m + 3;

    return 0;
}

int Hyperloglog::hll_cnt_reset()
{
    if (!ctx) {
        return -1;
    }

    ctx->err = CCARD_OK;
    memset(ctx->M, 0, ctx->m);

    return 0;
}

int Hyperloglog::hll_cnt_fini()
{
    if (ctx) {
        free(ctx);
        return 0;
    }

    return -1;
}

int Hyperloglog::hll_cnt_errnum(hll_cnt_ctx_t *ctx)
{
    if (ctx) {
        return ctx->err;
    }

    return CCARD_ERR_INVALID_CTX;
}

const char *Hyperloglog::hll_cnt_errstr(int err)
{
    static const char *msg[] = {
        "No error",
        "Invalid algorithm context",
        "Merge bitmap failed",
        NULL
    };

    if (-err >= 0 && -err < (int)(sizeof(msg) / sizeof(msg[0]) - 1)) {
        return msg[-err];
    }

    return "Invalid error number";
}

// vi:ft=c ts=4 sw=4 fdm=marker et

uint32_t Hyperloglog::murmurhash(void *buf, uint32_t len, uint32_t seed)
{
    uint8_t *data = (uint8_t *)buf;
    uint32_t m = 0x5bd1e995;
    uint32_t r = 24;
    uint32_t h = seed ^ len;
    uint32_t len_4 = len >> 2;
    uint32_t i;
    uint32_t len_m;
    uint32_t left;

    for(i = 0; i < len_4; i++) {
        uint32_t i_4 = i << 2;
        uint32_t k = data[i_4 + 3];
        k <<= 8;
        k |= data[i_4 + 2] & 0xff;
        k <<= 8;
        k |= data[i_4 + 1] & 0xff;
        k <<= 8;
        k |= data[i_4 + 0] & 0xff;
        k *= m;
        k ^= k >> r;
        k *= m;
        h *= m;
        h ^= k;
    }

    // avoid calculating modulo
    len_m = len_4 << 2;
    left = len - len_m;

    if (left != 0) {
        if (left >= 3) {
            h ^= data[len - 3] << 16;
        }
        if (left >= 2) {
            h ^= data[len - 2] << 8;
        }
        if (left >= 1) {
            h ^= data[len - 1];
        }
        h *= m;
    }

    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;

    return h;
}

uint32_t Hyperloglog::murmurhash_long(uint64_t data)
{
    uint32_t m = 0x5bd1e995;
    uint32_t r = 24;
    uint32_t h = 0;
    uint32_t k = (uint32_t)(data * m);

    k ^= k >> r;
    h ^= k * m;

    k = (data >> 32) * m;
    k ^= k >> r;
    h *= m;
    h ^= k * m;

    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;

    return h;
}

uint64_t Hyperloglog::murmurhash64(void *buf, uint32_t len, uint32_t seed)
{
    uint8_t *data = (uint8_t *)buf;
    uint64_t m = 0xc6a4a7935bd1e995L;
    uint32_t r = 47;
    uint64_t h = (seed & 0xffffffffl) ^ (len * m);
    uint32_t len8 = len / 8;
    uint32_t i;

    for (i = 0; i < len8; i++) {
        uint32_t i8 = i * 8;
        uint64_t k = ((uint64_t) data[i8 + 0] & 0xff) + (((uint64_t) data[i8 + 1] & 0xff) << 8)
                     + (((uint64_t) data[i8 + 2] & 0xff) << 16) + (((uint64_t) data[i8 + 3] & 0xff) << 24)
                     + (((uint64_t) data[i8 + 4] & 0xff) << 32) + (((uint64_t) data[i8 + 5] & 0xff) << 40)
                     + (((uint64_t) data[i8 + 6] & 0xff) << 48) + (((uint64_t) data[i8 + 7] & 0xff) << 56);

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    switch (len % 8) {
        case 7:
            h ^= (uint64_t) (data[(len & ~7) + 6] & 0xff) << 48;
        case 6:
            h ^= (uint64_t) (data[(len & ~7) + 5] & 0xff) << 40;
        case 5:
            h ^= (uint64_t) (data[(len & ~7) + 4] & 0xff) << 32;
        case 4:
            h ^= (uint64_t) (data[(len & ~7) + 3] & 0xff) << 24;
        case 3:
            h ^= (uint64_t) (data[(len & ~7) + 2] & 0xff) << 16;
        case 2:
            h ^= (uint64_t) (data[(len & ~7) + 1] & 0xff) << 8;
        case 1:
            h ^= (uint64_t) (data[len & ~7] & 0xff);
            h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

uint64_t Hyperloglog::murmurhash64_no_seed(void *buf, uint32_t len)
{
    return murmurhash64(buf, len, 0xe17a1465);
}

uint32_t Hyperloglog::lookup3(const uint32_t *k, uint32_t offset, uint32_t length, uint32_t initval)
{
    uint32_t a, b, c;
    uint32_t i = offset;
    a = b = c = 0xdeadbeef + (length << 2) + initval;

    while (length > 3) {
        a += k[i];
        b += k[i + 1];
        c += k[i + 2];

        a -= c;  a ^= (c << 4)  | (c >> (0x1f & -4));   c += b;
        b -= a;  b ^= (a << 6)  | (a >> (0x1f & -6));   a += c;
        c -= b;  c ^= (b << 8)  | (b >> (0x1f & -8));   b += a;
        a -= c;  a ^= (c << 16) | (c >> (0x1f & -16));  c += b;
        b -= a;  b ^= (a << 19) | (a >> (0x1f & -19));  a += c;
        c -= b;  c ^= (b << 4)  | (b >> (0x1f & -4));   b += a;

        length -= 3;
        i += 3;
    }

    switch(length) {
        case 3 : c += k[i + 2];  // fall through
        case 2 : b += k[i + 1];  // fall through
        case 1 : a += k[i + 0];  // fall through
            c ^= b; c -= (b << 14) | (b >> (0x1f & -14));
            a ^= c; a -= (c << 11) | (c >> (0x1f & -11));
            b ^= a; b -= (a << 25) | (a >> (0x1f & -25));
            c ^= b; c -= (b << 16) | (b >> (0x1f & -16));
            a ^= c; a -= (c << 4)  | (c >> (0x1f & -4));
            b ^= a; b -= (a << 14) | (a >> (0x1f & -14));
            c ^= b; c -= (b << 24) | (b >> (0x1f & -24));
        case 0:
            break;
    }
    return c;
}

uint32_t Hyperloglog::lookup3ycs(const uint32_t *k, uint32_t offset, uint32_t length, uint32_t initval)
{
    return lookup3(k, offset, length, initval - (length << 2));
}

uint32_t Hyperloglog::lookup3ycs_2(const char *s, uint32_t start, uint32_t end, uint32_t initval)
{
    uint32_t a, b, c;
    uint32_t i = start;
    a = b = c = 0xdeadbeef + initval;
    uint8_t mixed = 1;  // have the 3 state variables been adequately mixed?

    for(;;) {
        if (i >= end) break;
        mixed = 0;
        char ch;
        ch = s[i++];
        a += ch;
        if (i >= end) break;
        ch = s[i++];
        b += ch;
        if (i >= end) break;
        ch = s[i++];
        c += ch;
        if (i >= end) break;

        a -= c;  a ^= (c << 4)  | (c >> (0x1f & -4));   c += b;
        b -= a;  b ^= (a << 6)  | (a >> (0x1f & -6));   a += c;
        c -= b;  c ^= (b << 8)  | (b >> (0x1f & -8));   b += a;
        a -= c;  a ^= (c << 16) | (c >> (0x1f & -16));  c += b;
        b -= a;  b ^= (a << 19) | (a >> (0x1f & -19));  a += c;
        c -= b;  c ^= (b << 4)  | (b >> (0x1f & -4));   b += a;

        mixed = 1;
    }

    if (mixed == 0) {
        c ^= b; c -= (b << 14) | (b >> (0x1f & -14));
        a ^= c; a -= (c << 11) | (c >> (0x1f & -11));
        b ^= a; b -= (a << 25) | (a >> (0x1f & -25));
        c ^= b; c -= (b << 16) | (b >> (0x1f & -16));
        a ^= c; a -= (c << 4)  | (c >> (0x1f & -4));
        b ^= a; b -= (a << 14) | (a >> (0x1f & -14));
        c ^= b; c -= (b << 24) | (b >> (0x1f & -24));
    }

    return c;
}

uint64_t Hyperloglog::lookup3ycs64(const char *s, uint32_t start, uint32_t end, uint64_t initval)
{
    uint32_t a, b, c;
    uint32_t i = start;
    a = b = c = 0xdeadbeef + (uint32_t)initval;
    c += (uint32_t)(initval >> 32);

    uint8_t mixed = 1;  // have the 3 state variables been adequately mixed?
    for(;;) {
        if (i >= end) break;
        mixed = 0;
        char ch;
        ch = s[i++];
        a += ch;
        if (i >= end) break;
        ch = s[i++];
        b += ch;
        if (i >= end) break;
        ch = s[i++];
        c += ch;
        if (i >= end) break;

        a -= c;  a ^= ( c << 4)  | (c >> (0x1f & -4));   c += b;
        b -= a;  b ^= ( a << 6)  | (a >> (0x1f & -6));   a += c;
        c -= b;  c ^= ( b << 8)  | (b >> (0x1f & -8));   b += a;
        a -= c;  a ^= ( c << 16) | (c >> (0x1f & -16));  c += b;
        b -= a;  b ^= ( a << 19) | (a >> (0x1f & -19));  a += c;
        c -= b;  c ^= ( b << 4)  | (b >> (0x1f & -4));   b += a;

        mixed = 1;
    }

    if (mixed == 0) {
        c ^= b; c -= (b << 14) | (b >> (0x1f & -14));
        a ^= c; a -= (c << 11) | (c >> (0x1f & -11));
        b ^= a; b -= (a << 25) | (a >> (0x1f & -25));
        c ^= b; c -= (b << 16) | (b >> (0x1f & -16));
        a ^= c; a -= (c << 4)  | (c >> (0x1f & -4));
        b ^= a; b -= (a << 14) | (a >> (0x1f & -14));
        c ^= b; c -= (b << 24) | (b >> (0x1f & -24));
    }

    return c + ((uint64_t)b << 32);
}

uint64_t Hyperloglog::lookup3ycs64_2(const char *s)
{
    return lookup3ycs64(s, 0, strlen(s), -1);
}
