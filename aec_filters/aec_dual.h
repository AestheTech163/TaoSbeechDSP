//
// Created by wentao on 2024/1/14.
//

#ifndef TAOSBEECHDSP_AEC_DUAL_FILTER_H
#define TAOSBEECHDSP_AEC_DUAL_FILTER_H

#include "smallft.h"
#include "common_states.h"
#include "vector.h"
#include "mdf.h"
#include "aec_kalman.h"
#include "aec_pnlms.h"

typedef struct DualFilter_ {
    int fg_blocks;
    int bg_blocks;
//    int n_blocks;
//    int n_framelen;
//    int n_fft;
//    int n_freq;

    KalmanFilter* fg;
    pNLMSFilter* bg;
    SmallFFT* fft_handler;
    CommonStates* cs;
    get_residual_func get_residual;

    int saturated;
    int bg_better_cnt;
    int bg_update_cnt;
    int update_fg;
    int adapted;

    vector* window;

    vector_view left_wind;
    vector_view right_wind;

    flt_t Sxx;
    flt_t Davg1;
    flt_t Dvar1;
    flt_t ss;
    flt_t sum_adapt;
} DualFilter;


DualFilter* DualFilter_init(int fg_blocks, int bg_blocks, int n_framelen, flt_t A, CommonStates* cs);
void DualFilter_destroy(DualFilter* df);
void DualFilter_update_filter(void* filter, CommonStates* cs);
void DualFilter_process_frame(DualFilter* df, vector* out, const vector* y_frame, const vector* x_frame);

#endif //TAOSBEECHDSP_AEC_DUAL_FILTER_H
