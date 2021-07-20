
#pragma once
#include "habane_common.h"

void habane_downsample(
    h_stereo_sample16 *inbuffer,h_stereo_sample16 *inbuffer_prev,h_stereo_sample16 *inbuffer_next,
    uint blocklen, h_stereo_sample32 *lowbuffer,float *sinckernel);

void habane_encode_block(struct habane_block *block, 
    h_stereo_sample16 *inbuffer,h_stereo_sample16 *inbuffer_prev,h_stereo_sample16 *inbuffer_next,
    float *sinckernel, float *fftwin, double *noiseweight);

float* habane_makesinc();
float* habane_makefftwin();
double *habane_make_noiseweight();

