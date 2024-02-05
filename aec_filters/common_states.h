//
// Created by wentao on 2024/1/19.
//

#ifndef TAOSBEECHDSP_COMMON_STATE_H
#define TAOSBEECHDSP_COMMON_STATE_H

#include "arch.h"
#include "vector.h"

typedef struct CommonStates_ {
    flt_t preemph;
    flt_t pre_mem_x;
    flt_t pre_mem_y;
    flt_t depre_mem_o;
    flt_t notch_radius;
    flt_t notch_mem[2];
    int saturated;

    flt_t Sxx;
    vector* x;
    vector* y;
    vector* last_y;
    vector* X_fifo;
    vector* X2_fifo;
    vector* X;
    vector* X2;
    vector* Y;
//    vector* R2_mean;
//    vector* Y2_mean;

    flt_t min_leak;
    vector* Rf;
    vector* Eh;
    vector* Yf;
    vector* Yh;
    vector* Xf;
    flt_t Pey;
    flt_t Pyy;
    flt_t beta0;
    flt_t beta_max;
    flt_t leak_estimate;

} CommonStates;

CommonStates* CommonStates_init(int n_blocks, int n_framelen,
                                int n_freq, int n_fft, int sample_rate,
                                flt_t preemph);
void CommonStates_destroy(CommonStates* cs);
void compute_output(CommonStates* cs, vector* error, vector* out);

#endif //TAOSBEECHDSP_COMMON_STATE_H
