
#include "habane_decode.h"
#include "habane_encode.h"
#include <stdio.h>


static int find_best_predictor(int32_t target,int32_t h1,int32_t h2) {
    uint best_predictor = 0;
    int32_t best_diff = INT32_MAX;
    for (uint p=0;p<512;p++) {
        uint h1sar = p&15;
        uint h2sar = (p>>4)&15;
        bool negh2 = p&256;

        int32_t prediction = (h1>>h1sar) + ((negh2?-h2:h2)>>h2sar);
        int32_t diff = abs(target-prediction);
        if (diff < best_diff) {
            best_diff = diff;
            best_predictor = p;
        }
    }
    return best_predictor;
}

int main(int argc,char** argv) {

    if (argc < 2) {
        printf("Error: no arguments\n");
        return -1;
    }

    if (!strcmp(argv[1],"noisetest")) {
        FILE* f = fopen("noisetest.raw","wb");
        struct habane_noise_state noise_state;
        int16_t *buffer = calloc(48000*120,sizeof(int16_t));
        habane_noise_init(&noise_state);
        habane_noise_test(&noise_state,buffer,48000*120);
        fwrite(buffer,sizeof(int16_t),48000*120,f);
        free(buffer);
        fclose(f);
    } else if (!strcmp(argv[1],"sinctest")) {
        float *sinckernel = habane_makesinc();
        for (int i=0;i<HABANE_KERNEL_SIZE;i++) {
            printf("%2d: %f\n",i,sinckernel[i]);
        }
        free(sinckernel);
    } else if (!strcmp(argv[1],"encode")) {
        if (argc < 4) {
            printf("Error: not enough arguments\n");
            return -1;
        }
        FILE* in = fopen(argv[2],"rb");
        FILE* out = fopen(argv[3],"wb");
        float *sinckernel = habane_makesinc();
        float *fftwin = habane_makefftwin();
        double *noiseweight = habane_make_noiseweight();
        h_stereo_sample16 buffer[3][HABANE_BLOCKLENGTH*2] = {};
        for (int i=0;;i++) {
            int read = fread(buffer[i%3],sizeof(h_stereo_sample16),HABANE_BLOCKLENGTH*2,in);
            struct habane_block block = {};
            if (i>=1) {
                habane_encode_block(&block,buffer[(i-1)%3],i>=2?buffer[(i-2)%3]:NULL,read?buffer[(i-0)%3]:NULL,sinckernel,fftwin,noiseweight);
                fwrite(&block,sizeof(struct habane_block),1,out);
            }
            if (!read) break;
        }
        free(sinckernel);
        free(fftwin);
        free(noiseweight);
        fclose(out);
        fclose(in);
    } else if (!strcmp(argv[1],"decode")||!strcmp(argv[1],"decode-noise")||!strcmp(argv[1],"decode-adpcm")) {
        if (argc < 4) {
            printf("Error: not enough arguments\n");
            return -1;
        }
        bool no_noise = !strcmp(argv[1],"decode-adpcm");
        bool no_adpcm = !strcmp(argv[1],"decode-noise");
        FILE* in = fopen(argv[2],"rb");
        FILE* out = fopen(argv[3],"wb");
        struct habane_noise_state noise_state;
        habane_noise_init(&noise_state);
        for (;;) {
            struct habane_block block = {};
            int read = fread(&block,sizeof(struct habane_block),1,in);
            if (read == 0) break;
            h_stereo_sample16 buffer[HABANE_BLOCKLENGTH*2];
            habane_decode_block(&noise_state,&block,buffer,no_noise,no_adpcm);
            fwrite(buffer,sizeof(h_stereo_sample16),HABANE_BLOCKLENGTH*2,out);
        }
        fclose(out);
        fclose(in);
    } else if (!strcmp(argv[1],"prediction_analyze")) {
        if (argc != 3) {
            printf("Error: not enough arguments\n");
            return -1;
        }
        FILE* in = fopen(argv[2],"rb");
        h_stereo_sample16 history[2] = {{},{}};
        uint32_t bestcount[512] = {};
        fread(&history[1],sizeof(h_stereo_sample16),1,in);
        fread(&history[0],sizeof(h_stereo_sample16),1,in);
        for (;;) {
            h_stereo_sample16 s;
            if (!fread(&s,sizeof(h_stereo_sample16),1,in)) break;
            bestcount[find_best_predictor(s.l,history[0].l,history[1].l)]++;
            bestcount[find_best_predictor(s.r,history[0].r,history[1].r)]++;
            history[1] = history[0];
            history[0] = s;
        }
        for (uint p=0;p<512;p++) {
            printf("%03X : %d\n",p,bestcount[p]);
        }
    } else {
        printf("Unknown command\n");
    }
    return 0;
}


