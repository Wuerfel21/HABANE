
#include "habane_encode.h"
#include <math.h>
#include <stdio.h>
#include "fft-complex.h"

inline int32_t round_sar(int32_t val,uint sar) {
    return val/(1<<sar);
    //return val>>sar;
}

inline h_stereo_sample16 get_sample(int n, int blocksize,h_stereo_sample16 *inbuffer,h_stereo_sample16 *inbuffer_prev,h_stereo_sample16 *inbuffer_next) {
    if (n>=0 && n<blocksize) return inbuffer[n];
    if (n>=(-blocksize)&& n < 0 &&  inbuffer_prev) return inbuffer_prev[n+blocksize];
    if (n >= blocksize && n<(2*blocksize) && inbuffer_next) return inbuffer_next[n-blocksize];

    h_stereo_sample16 nullsample = {};
    return nullsample;
}

static int double_compare(const void * a, const void * b) {
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;  
}


void habane_encode_block(struct habane_block *block, 
    h_stereo_sample16 *inbuffer,h_stereo_sample16 *inbuffer_prev,h_stereo_sample16 *inbuffer_next,
    float *sinckernel, float *fftwin, double *noiseweight) {

    // Downsample low band
    h_stereo_sample32 lowbuffer[HABANE_BLOCKLENGTH];
    for (int i=0;i<HABANE_BLOCKLENGTH;i++) {
        float accu_l = 0, accu_r = 0;
        for (int n=0;n<HABANE_KERNEL_SIZE;n++) {
            h_stereo_sample16 s = get_sample(i*2+n-(HABANE_KERNEL_SIZE/2),HABANE_BLOCKLENGTH*2,inbuffer,inbuffer_prev,inbuffer_next);
            float sinc = sinckernel[n];//n==32?1.f:0.f;//sinckernel[n];
            accu_l += sinc*s.l;
            accu_r += sinc*s.r;
        }
        h_stereo_sample32 os = {roundf(accu_l),roundf(accu_r)};
        lowbuffer[i] = os;
    }

    h_stereo_sample32 history[2] = {lowbuffer[1], lowbuffer[0]};
    block->raw1l = habane_clampto16(lowbuffer[0]).l;
    block->raw1r = habane_clampto16(lowbuffer[0]).r;
    block->raw2l = habane_clampto16(lowbuffer[1]).l;
    block->raw2r = habane_clampto16(lowbuffer[1]).r;

    for (uint ui=0;ui<HABANE_UNITS_PER_BLOCK;ui++) {
        uint unitbase = 2 + ui*HABANE_UNITLENGTH;

        h_stereo_sample32 new_history[256][2];
        uint8_t unitdata[256][HABANE_UNITLENGTH];
        uint64_t encoder_error[256];

        // Try encoding in all 256 possible modes
        for (uint mode=0;mode<256;mode++) {
            uint scale = mode&15, predictor = (mode >> 4)&15;
            new_history[mode][0] = history[0];
            new_history[mode][1] = history[1];
            encoder_error[mode] = 0;
            for (uint n=0;n<HABANE_UNITLENGTH;n++) {
                h_stereo_sample32 prediction = {
                    habane_predictor_value(predictor,new_history[mode][0].l,new_history[mode][1].l),
                    habane_predictor_value(predictor,new_history[mode][0].r,new_history[mode][1].r),
                };
                h_stereo_sample32 insample = lowbuffer[unitbase+n];
                int32_t diff_l = insample.l - prediction.l;
                int32_t diff_r = insample.r - prediction.r;
                h_stereo_sample32 new_sample;
                if (scale>=13) { // Joint mode
                    int32_t diff_joint = (diff_l+diff_r)>>1;
                    int32_t quantized_diff = h_clamp(round_sar(diff_joint,(scale-13)<<2),-128,127);
                    new_sample.l = prediction.l + (quantized_diff<<((scale-13)<<2));
                    new_sample.r = prediction.r + (quantized_diff<<((scale-13)<<2));
                    unitdata[mode][n] = (uint8_t)quantized_diff;

                } else { // Stereo mode
                    int32_t quantized_l = h_clamp(round_sar(diff_l,scale),-8,7);
                    int32_t quantized_r = h_clamp(round_sar(diff_r,scale),-8,7);
                    new_sample.l = prediction.l + (quantized_l<<scale);
                    new_sample.r = prediction.r + (quantized_r<<scale);
                    unitdata[mode][n] = ((uint8_t)quantized_l & 0xF) | ((uint8_t)quantized_r << 4);
                }
                new_history[mode][1] = new_history[mode][0];
                new_history[mode][0] = new_sample;

                int64_t error_l = (int64_t)new_sample.l - (int64_t)insample.l;
                int64_t error_r = (int64_t)new_sample.r - (int64_t)insample.r;
                encoder_error[mode] += (error_l*error_l);
                encoder_error[mode] += (error_r*error_r);
                //if (predictor!=2) encoder_error[mode] += 1000000000000U;
            }
        }

        // Find the best one
        uint best_coding = 0;
        uint64_t lowest_error = UINT64_MAX;
        for (uint mode=0;mode<256;mode++) {
            if (encoder_error[mode] < lowest_error) {
                best_coding = mode;
                lowest_error = encoder_error[mode];
            }
        }
        //best_coding = 0x2F; // DEBUG

        // put it in
        block->units[ui].coding = best_coding;
        memcpy(&block->units[ui].adpcm,&unitdata[best_coding],sizeof(uint8_t)*HABANE_UNITLENGTH);

        // prepare for next one
        history[0] = new_history[best_coding][0];
        history[1] = new_history[best_coding][1];

        // Gather samples for noise detector
        double complex fftbuf_l[HABANE_FFTSIZE],fftbuf_r[HABANE_FFTSIZE];
        for (int i=0;i<HABANE_FFTSIZE;i++) {
            int idx = unitbase + i - ((HABANE_FFTSIZE-HABANE_UNITLENGTH)/2);
            h_stereo_sample16 s = get_sample(idx,HABANE_BLOCKLENGTH*2,inbuffer,inbuffer_prev,inbuffer_next);
            float win = fftwin[i];
            fftbuf_l[i] = win*s.l;
            fftbuf_r[i] = win*s.r;
        }
        Fft_transform(fftbuf_l,HABANE_FFTSIZE,false);
        Fft_transform(fftbuf_r,HABANE_FFTSIZE,false);

        double binbuffer_l[HABANE_FFTSIZE/2],binbuffer_r[HABANE_FFTSIZE/2];
        //double noiseacc_l=0,noiseacc_r=0;
        for (int i=0;i<HABANE_FFTSIZE/2;i++) {
            binbuffer_l[i] = fabs(creal(fftbuf_l[i+(HABANE_FFTSIZE/2)]) * noiseweight[i+(HABANE_FFTSIZE/2)]);
            binbuffer_r[i] = fabs(creal(fftbuf_r[i+(HABANE_FFTSIZE/2)]) * noiseweight[i+(HABANE_FFTSIZE/2)]);
        }
        qsort(binbuffer_l,HABANE_FFTSIZE/2,sizeof(double),&double_compare);
        qsort(binbuffer_r,HABANE_FFTSIZE/2,sizeof(double),&double_compare);

        uint noiselev_l = h_clamp(round(16-log2f(0.00001f+powf(binbuffer_l[HABANE_FFTSIZE/4],1.0f))),0,15);
        uint noiselev_r = h_clamp(round(16-log2f(0.00001f+powf(binbuffer_r[HABANE_FFTSIZE/4],1.0f))),0,15);

        //uint noiselev_l = h_clamp(round(16-log2f(0.00001f+powf(noiseacc_l/HABANE_FFTSIZE,1.4f))),0,15);
        //uint noiselev_r = h_clamp(round(16-log2f(0.00001f+powf(noiseacc_r/HABANE_FFTSIZE,1.4f))),0,15);

        //printf("Unit %2d: %2d %2d\n",ui,noiselev_l,noiselev_r);

        block->units[ui].noiselev = noiselev_l + (noiselev_r<<4);
        

    }

}

float* habane_makefftwin() {
    float *buffer = calloc(HABANE_FFTSIZE,sizeof(float));
    for (int i=0;i<HABANE_FFTSIZE;i++) {
        float window = 0.35875f - 0.48829f*cosf((2*M_PI*i)/HABANE_KERNEL_SIZE) + 0.14128f*cosf((4*M_PI*i)/HABANE_KERNEL_SIZE) - 0.01168f*cosf((6*M_PI*i)/HABANE_KERNEL_SIZE);
        buffer[i] = window;
    }
    return buffer;
}

static float windowfun(int i,int size) {
    // Blackmann-Harris window function, I guess
    return 0.35875f - 0.48829f*cosf((2*M_PI*i)/size) + 0.14128f*cosf((4*M_PI*i)/size) - 0.01168f*cosf((6*M_PI*i)/size);
}

float* habane_makesinc() {
    float *buffer = calloc(HABANE_KERNEL_SIZE,sizeof(float));
    for (int i=0;i<HABANE_KERNEL_SIZE;i++) {
        float t = (i-(HABANE_KERNEL_SIZE/2))/2.f;
        float sinc = t == 0.0 ? 1.f : sinf(t*M_PI)/(t*M_PI);
        buffer[i] = windowfun(i,HABANE_KERNEL_SIZE) * sinc * 0.5f;
    }
    return buffer;
}

static double noisesens(double f) {
    // Psychoacoustic noise frequency response, according to ITU-R 468
    // Might try a more complex algorithm, but we're only interested in the high band, so eh
    double h1 = -4.737338981378384e-24*pow(f,6) + 2.04382833606125e-15*pow(f,4) - 1.363894795463638e-7*pow(f,2)+1;
    double h2 = 1.306612257412824e-19*pow(f,5) - 2.118150887518656e-15*pow(f,3) + 5.559488023498642e-4*f;
    double fsens = (1.246332637532143e-4*f)/sqrt(h1*h1+h2*h2);

    double scale = 5500;

    // Low bins shouldn't respond very much
    double wall = fmax(0,1-pow(2,-0.001*(f-23700)));

    return fsens*wall*scale;

}

double* habane_make_noiseweight(double f) {
    double *buffer = calloc(HABANE_FFTSIZE,sizeof(double));
    for (int i=0;i<HABANE_FFTSIZE;i++) {
        double f = (i+0.5) * ((float)HABANE_RATE/HABANE_FFTSIZE);
        buffer[i] = noisesens(f);
    }
    return buffer;
}


