
#include "habane_decode.h"

const struct habane_predictor habane_predictors[16] = {
    {15,15},
    {0,0},
    {0,15},
    {1,1},

    {0,5},
    {0,5},
    {0,4},
    {0,4}, // Apparently the best one?
    {0,3},
    {0,3},

    {0,0},
    {0,6},  
    {0,1},
    {0,7}, 
    {2,0},
};


h_stereo_sample16 habane_noise_iterate(struct habane_noise_state *state,int level_l, int level_r) {

    state->lfsr = (state->lfsr << 1) + (__builtin_parity(state->lfsr&HABANE_LFSR_TAPS) ? 1 : 0);

    int32_t lfsr_sample = ((int32_t)state->lfsr)>>19;
    int32_t diff = lfsr_sample - state->filter1;

    state->filter1 += diff;

    int32_t diff2 = diff - state->filter2;
    state->filter2 += diff2;

    int32_t diff3l = (diff2 >> level_l) - state->filter3_l;
    int32_t diff3r = (diff2 >> level_r) - state->filter3_r;

    state->filter3_l += diff3l;
    state->filter3_r += diff3r;

    h_stereo_sample16 rv = {(int16_t)diff3l,(int16_t)diff3r};
    //struct h_stereo_sample rv = {(int16_t)lfsr_sample,(int16_t)lfsr_sample};
    return rv;
}

void habane_noise_init(struct habane_noise_state *state) {
    state->lfsr = HABANE_LFSR_INIT;
    state->filter1 = 0;
    state->filter2 = 0;
    state->filter3_l = 0;
    state->filter3_r = 0;
}

void habane_noise_test(struct habane_noise_state *noise_state,int16_t *buffer,int samples) {
    for (int i=0;i<samples;i++) {
        struct h_stereo_sample16 rv = habane_noise_iterate(noise_state,0,0);
        buffer[i] = rv.l;
    }
}

void habane_decode_block(struct habane_noise_state *noise_state, struct habane_block *block, h_stereo_sample16 *buffer, bool no_noise, bool no_adpcm) {

    h_stereo_sample32 history[2];
    h_stereo_sample16 raw1 = {block->raw1l,block->raw1r}, raw2 = {block->raw2l,block->raw2r};
    history[1] = habane_expand32(raw1);
    history[0] = habane_expand32(raw2);

    uint initnoise_l = block->units[0].noiselev&15;
    uint initnoise_r = (block->units[0].noiselev>>4)&15;
    if (no_adpcm) {
        buffer[0] = habane_noise_iterate(noise_state,initnoise_l,initnoise_r);
        buffer[1] = habane_noise_iterate(noise_state,initnoise_l,initnoise_r);
        buffer[2] = habane_noise_iterate(noise_state,initnoise_l,initnoise_r);
        buffer[3] = habane_noise_iterate(noise_state,initnoise_l,initnoise_r);
    } else if (no_noise) {
        
        buffer[0] = raw1;
        buffer[1] = raw1;
        buffer[2] = raw2;
        buffer[3] = raw2;
    } else {
        buffer[0] = habane_mix16(raw1,habane_noise_iterate(noise_state,initnoise_l,initnoise_r));
        buffer[1] = habane_mix16(raw1,habane_noise_iterate(noise_state,initnoise_l,initnoise_r));
        buffer[2] = habane_mix16(raw2,habane_noise_iterate(noise_state,initnoise_l,initnoise_r));
        buffer[3] = habane_mix16(raw2,habane_noise_iterate(noise_state,initnoise_l,initnoise_r));
    }

    for (uint ui=0;ui<HABANE_UNITS_PER_BLOCK;ui++) {
        struct habane_unit *unit = &block->units[ui];
        uint8_t mode = unit->coding;
        uint scale = mode&15, predictor = (mode >> 4)&15;
        for (uint n=0;n<HABANE_UNITLENGTH;n++) {
            h_stereo_sample32 prediction = {
                habane_predictor_value(predictor,history[0].l,history[1].l),
                habane_predictor_value(predictor,history[0].r,history[1].r),
            };
            h_stereo_sample32 new_sample;
            if (scale>=13) { // Joint mode
                int32_t quantized_diff = (int8_t)unit->adpcm[n];
                new_sample.l = prediction.l + (quantized_diff<<((scale-13)<<2));
                new_sample.r = prediction.r + (quantized_diff<<((scale-13)<<2));
            } else { // Stereo mode
                uint8_t data = unit->adpcm[n];
                int32_t quantized_l = ((int8_t)(getnib0(data)<<4))>>4;
                int32_t quantized_r = ((int8_t)(getnib1(data)<<4))>>4;
                new_sample.l = prediction.l+(quantized_l<<scale);
                new_sample.r = prediction.r+(quantized_r<<scale);
            }
            history[1] = history[0];
            history[0] = new_sample;
            uint sample_idx = 2+ui*HABANE_UNITLENGTH+n;

            h_stereo_sample16 clamped_sample = {};
            if (!no_adpcm) clamped_sample = habane_clampto16(new_sample);
            if (no_noise) {
                buffer[sample_idx*2]   = clamped_sample;
                buffer[sample_idx*2+1] = clamped_sample;
            } else {
                buffer[sample_idx*2]   = habane_mix16(clamped_sample,habane_noise_iterate(noise_state,unit->noiselev&15,unit->noiselev>>4));
                buffer[sample_idx*2+1] = habane_mix16(clamped_sample,habane_noise_iterate(noise_state,unit->noiselev&15,unit->noiselev>>4));
            }
        }
    }
}

