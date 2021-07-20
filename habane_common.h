
#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef unsigned int uint;

struct h_stereo_sample16 {
    int16_t l,r;
} typedef h_stereo_sample16;

struct h_stereo_sample32 {
    int32_t l,r;
} typedef h_stereo_sample32;

inline int16_t add16sat(int16_t a,int16_t b) {
    int32_t tmp = (int32_t)a + (int32_t)b;
    if (tmp < INT16_MIN) return INT16_MIN;
    if (tmp > INT16_MAX) return INT16_MAX;
    return tmp;
}

inline int h_sign(int x) {
    return (x > 0) - (x < 0);
}

inline int h_clamp(int val,int min,int max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

inline h_stereo_sample16 habane_mix16(h_stereo_sample16 a,h_stereo_sample16 b) {
    h_stereo_sample16 rv = {
        add16sat(a.l,b.l),add16sat(a.r,b.r)
    };
    return rv;
}

inline h_stereo_sample16 habane_clampto16(h_stereo_sample32 s) {
    h_stereo_sample16 rv = {
        h_clamp(s.l,INT16_MIN,INT16_MAX),h_clamp(s.r,INT16_MIN,INT16_MAX)
    };
    return rv;
}

inline h_stereo_sample32 habane_expand32(h_stereo_sample16 s) {
    h_stereo_sample32 rv = {
        s.l,s.r
    };
    return rv;
}

enum {
    HABANE_RATE = 48000,
    HABANE_LFSR_TAPS = 0x80000057u,
    HABANE_LFSR_INIT = 0xABCDEF00u,
    HABANE_KERNEL_SIZE = 64,
    HABANE_UNITLENGTH = 16,
    HABANE_FFTSIZE = 128,
    HABANE_UNITS_PER_BLOCK = 27,
    HABANE_BLOCKLENGTH_CODED = HABANE_UNITLENGTH*HABANE_UNITS_PER_BLOCK,
    HABANE_BLOCKLENGTH = HABANE_BLOCKLENGTH_CODED + 2,
};

struct habane_predictor {
    uint h1_sar,h2_sar;
};

inline uint8_t getnib1(uint8_t b) {
    return b >> 4;
}
inline uint8_t getnib0(uint8_t b) {
    return b & 0xf;
}
inline void setnib1(uint8_t* b,uint nib) {
    *b = (*b&0x0f)+((nib&0x0f)<<4);
}
inline void setnib0(uint8_t* b,uint nib) {
    *b = (*b&0xf0)+((nib&0x0f));
}

// note: odd-numbered predictors invert h2
extern const struct habane_predictor habane_predictors[16];

inline int32_t habane_predictor_value(uint predictor, int32_t h1, int32_t h2) {
    h1 <<= 2;
    h2 <<= 2;
    if (predictor&1) h2 = -h2;
    return (h1 >> habane_predictors[predictor].h1_sar) + (h2 >> habane_predictors[predictor].h2_sar);
}

struct habane_unit {
    uint8_t coding;
    uint8_t noiselev;
    uint8_t adpcm[HABANE_UNITLENGTH];
} __attribute__((packed));

struct habane_block {
    int16_t raw1l,raw1r;
    int16_t raw2l,raw2r;
    uint32_t branch_condition;
    int32_t  branch_offset;
    uint8_t padding[512-(16+sizeof(struct habane_unit)*HABANE_UNITS_PER_BLOCK)];
    struct habane_unit units[HABANE_UNITS_PER_BLOCK];
} __attribute__((packed));

