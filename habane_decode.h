
#pragma once
#include "habane_common.h"

struct habane_noise_state {
    uint32_t lfsr;
    int32_t filter1;
    int32_t filter2;
    int32_t filter3_l,filter3_r;
};

h_stereo_sample16 habane_noise_iterate(struct habane_noise_state *state,int level_l, int level_r);

void habane_noise_init(struct habane_noise_state *state);

void habane_noise_test(struct habane_noise_state *noise_state,int16_t *buffer,int samples);

void habane_decode_block(struct habane_noise_state *noise_state, struct habane_block *block, h_stereo_sample16 *buffer, bool no_noise, bool no_adpcm);

