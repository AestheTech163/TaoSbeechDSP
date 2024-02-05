//
// Created by wentao on 2024/1/28.
//
#include "stdio.h"
#include "math.h"

#include "dsp_common_func.h"
#include "f0tracker.h"
#include "mmse_noise.h"
#include "arch.h"
#include "io_utils.h"

// TODO 使用rnnoise中的pitch增益公式对子带进行滤波，目前整个帧作为整体进行滤波

void get_comb_hann_window(vector* window) {
    check(window->size%2==1, "comb window must be odd length.");
    int half_len = window->size / 2;
    int len = half_len * 2 + 1;
    flt_t sum_ = 0.;
    flt_t* win_data = window->data;
    for (int i=1; i<half_len*2+2; i++) {
        win_data[i-1] = 0.5 - 0.5*cos(2.0*M_PI*i/(half_len*2.+2));
        sum_ += win_data[i-1];
    }
    for (int i=0; i<len; i++)
        win_data[i] /= sum_;
}

void get_comb_half_window(vector* window) {
    vector tmp_window;
    auto_vector_construct(tmp_window, window->size*2+1, Float);
    flt_t sum_ = 0;
    flt_t* window_data = window->data, *tmp_window_data = tmp_window.data;
    get_comb_hann_window(&tmp_window);
    for (int i=0; i<window->size; i++) {
        window_data[i] = tmp_window_data[window->size+i];
        sum_ += window_data[i];
    }
    for (int i=0; i<window->size; i++)
        window_data[i] /= sum_;
}

void comb_fitering(int n0, const vector* comb_window, int comb_order,
                   const vector* buf, int begin, vector* y, int N) {
    int half_comb_o = comb_order / 2;
    vector_view buf_vw;
    flt_t w;

    vec_fill(y, &flt_zero);
    for (int o=-half_comb_o; o<half_comb_o+1; o++) {
        vec_view_construct(buf, begin + n0 * o, N, &buf_vw);
        vec_elem_get(comb_window, half_comb_o + o, &w);
        vec_mul_scalar_accum(&buf_vw, &w, y);
    }
}

void comb_fitering_halfwindow(int n0, const vector* comb_window, int comb_order,
                              const vector* buf, int begin, vector* y, int N) {
    vector_view buf_vw;
    flt_t w;
    vec_fill(y, &flt_zero);
    for (int o=0; o<comb_order; o++) {
        vec_view_construct(buf, begin + n0 * o, N, &buf_vw);
        vec_elem_get(comb_window, o, &w);
        vec_mul_scalar_accum(&buf_vw, &w, y);
    }
}

int main() {
    FILE *in_fd, *out_fd, *f0_fd, *comb_fd;
//    const char* in_file = "/Users/wentao/PycharmProjects/audio11/main_white_snr10.pcm";
//    const char* in_file = "/Users/wentao/PycharmProjects/aec-scratch/testset/testf0_1.wav";
//    const char* in_file = "/Users/wentao/PycharmProjects/audio11/percepnet-noisy-16k.wav";
//    const char* in_file = "/Users/wentao/PycharmProjects/AEC-scratch2/log-mmse-ns/log-mmse-ns_man_white_snr10_7.wav";
    const char* in_file = "/Users/wentao/PycharmProjects/audio11/percepnet-clip-16k.wav";
    const char* out_file = "./percepnet-clip-16k_mmsenoise.pcm";
    const char* f0_file = "./percepnet-clip-16k_f0.bin";
    const char* comb_file = "./percepnet-clip-16k_mmse_comb_.pcm";
    in_fd = fopen(in_file, "rb");
    check(in_fd, "open input file error.");
    out_fd = fopen(out_file, "wb");
    check(out_fd, "open output file error.");
    f0_fd = fopen(f0_file, "wb");
    check(f0_fd, "open f0 file error.");
    comb_fd = fopen(comb_file, "wb");
    check(comb_fd, "open comb file error.");

    int fs = 16000;
    int stepsize = 160;
    flt_t notch_mem[2] = {0, 0};
    int comb_window_len = 5;
    vector comb_window;
    auto_vector_construct(comb_window, comb_window_len, Float);
    get_comb_hann_window(&comb_window);
//    get_comb_half_window(comb_window, comb_window_len);
    printf("comb_window: "), vec_print(&comb_window);

    int i, n0;
    vector notch_buf;
    auto_vector_construct(notch_buf, stepsize, Float);
    flt_t corr, pitch;
    short tmp_buf[stepsize];
    int comb_buf_len = stepsize * 8;
    int comb_overlap = comb_buf_len - stepsize;
    vector comb_buf, comb_y; //[comb_buf_len];
    auto_vector_construct(comb_buf, comb_buf_len, Float);
    auto_vector_construct(comb_y, stepsize, Float);
    vec_fill(&comb_buf, &flt_zero);
    vec_fill(&comb_y, &flt_zero);
    size_t n;

    MMSENoise* nst = init_MMSENoise(stepsize, fs, 24);
    PitchState* pst = init_pitch(fs, stepsize, stepsize*2,100, 1000);

    while (1) {
        fread(tmp_buf, sizeof(short), stepsize, in_fd);
        if (feof(in_fd)) break;
        for (i = 0; i < stepsize; i++) vec_elem_set_FLT(&notch_buf, i, tmp_buf[i]);
        dc_notch_vec(notch_mem, &notch_buf, &notch_buf, FLT(0.972));

        mmsenoise_process_frame(nst, &notch_buf, NULL);
        for (i=0; i<stepsize; i++) tmp_buf[i] = (short)vec_elem_get_FLT(&notch_buf, i);
        fwrite(tmp_buf, sizeof(short), stepsize, out_fd);

        vec_left_moving(&comb_buf, &notch_buf, comb_overlap, stepsize, 0);
//        f0tracker_process_frame(pst, comb_buf + comb_overlap - 2 * stepsize, &n0, &corr);
        vector_view comb_buf_vw;
        vec_view_construct(&comb_buf, comb_overlap-3*stepsize, stepsize, &comb_buf_vw);
        f0tracker_process_frame(pst, &comb_buf_vw, &n0, &corr);

        if (n0 > 0) pitch = (flt_t)fs / (flt_t)n0;
        else pitch = FLT(0.);
        n = fwrite(&pitch, sizeof(float), 1, f0_fd);
        if (n < 1) fprintf(stderr, "Problems writing pitch in FilterStream\n");

        if ((n0 > 31) && (n0 < 201)) {
            comb_fitering(n0, &comb_window, comb_window_len,
                          &comb_buf, comb_overlap-3*stepsize, &comb_y, stepsize);
        } else {
            vec_copy(&comb_buf_vw, &comb_y);
        }
//        for (int i=0; i<stepsize;i++) printf("%f,", comb_y[i]);
//        mmsenoise_process_frame(nst2, comb_y, NULL);
        vec_clip(&comb_y, -32767, 32767);
        for (i=0; i<stepsize; i++)
            tmp_buf[i] = (short)vec_elem_get_FLT(&comb_y, i);
        n = fwrite(tmp_buf, sizeof(short), stepsize, comb_fd);
        if (n < stepsize) fprintf(stderr, "Problems writing pitch in FilterStream\n");
    }

    free_pitch(pst);
    free_MMSENoise(nst);
    fclose(in_fd);
    fclose(out_fd);
    fclose(f0_fd);
    fclose(comb_fd);

    return 0;
}
