//
// Created by wentao on 2024/1/25.
//
#include <stdio.h>

#include "arch.h"
#include "vector.h"
#include "aec_kalman.h"
#include "dsp_common_func.h"
#include "io_utils.h"


int main() {
    int fs = 16000;
    int stepsize = 160;
    flt_t preemph = FLT(0.);
    flt_t notch_radius = FLT(0.982);
    flt_t notch_mem[2] = {0, 0};

    int filter_tap_length = 8000;
    int n_blocks = (filter_tap_length + stepsize - 1) / stepsize;
    printf("n_blocks=%d\n", n_blocks);
    CommonStates* cs = CommonStates_init(n_blocks, stepsize,
                                         stepsize+1, stepsize*2,
                                         fs, preemph);
    KalmanFilter* kf = KalmanFilter_init(n_blocks, stepsize, 0.998, cs);

    const char* mic_file = "../test_files/test_mic.pcm";
    const char* far_file = "../test_files/test_far.pcm";
    const char* out_file = "../test_files/results/aec_out_kalman.pcm";
    FILE* mic_fd = fopen(mic_file, "rb");
    check(mic_fd, "open mic file error.");
    FILE* far_fd = fopen(far_file, "rb");
    check(far_fd, "open far file error.");
    FILE* out_fd = fopen(out_file, "wb");
    check(out_fd, "open out file error.");

    short mic[stepsize];
    short far[stepsize];
    vector mic_notch, far_notch, out;
    auto_vector_construct(mic_notch, stepsize, Float);
    auto_vector_construct(far_notch, stepsize, Float);
    auto_vector_construct(out, stepsize, Float);

    int tot_samples = file_size(mic_fd) / 2;
    int tot_frames = (tot_samples+stepsize-1)/stepsize;

    int frame_num = 0;
    int n;
    while (1) {
        frame_num++;
        printf("Processing frame: %d/%d\r", frame_num, tot_frames);

        fread(mic, sizeof(short), stepsize, mic_fd);
        if (feof(mic_fd)) break;
        fread(far, sizeof(short), stepsize, far_fd);
        if (feof(far_fd)) break;
        for (int i = 0; i < stepsize; i++) vec_elem_set_FLT(&mic_notch, i, mic[i]);
        for (int i = 0; i < stepsize; i++) vec_elem_set_FLT(&far_notch, i, far[i]);
        dc_notch_vec(&notch_mem[0], &mic_notch, &mic_notch, notch_radius);

        KalmanFilter_process_frame(kf, &out, &mic_notch, &far_notch);

        vec_clip(&out, -32768, 32767);
        for (int i=0; i<stepsize; i++) mic[i] = (short)vec_elem_get_FLT(&out, i);

        n = fwrite(mic, sizeof(short), stepsize, out_fd);
        if (n < stepsize) fprintf(stderr, "Problems writing output.\n");
    }

    fclose(mic_fd);
    fclose(far_fd);
    fclose(out_fd);

    KalmanFilter_destroy(kf);
    CommonStates_destroy(cs);

    return 0;
}
