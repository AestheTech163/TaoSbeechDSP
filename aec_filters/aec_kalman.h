//
// Created by wentao on 2024/1/12.
//

#ifndef TAOSBEECHDSP_AEC_KALMAN_H
#define TAOSBEECHDSP_AEC_KALMAN_H

#include "arch.h"
#include "mdf.h"

typedef struct KalmanFilter_ {
    MDFFilter* mdf;  // parent class instance
    CommonStates* cs;

    compute_error_func compute_error;
    compute_echo_func compute_echo;
    update_filter_func update_filter;

    flt_t A;  // 转移系数（isotropic转移矩阵）
    flt_t A2;
    flt_t A2_1;

    vector* E2_fifo;
    vector* P;
    vector* H2;
    vector* PHIss;
    vector* PHIee;
    vector* mu;
    vector* H_delta;

} KalmanFilter;

KalmanFilter* KalmanFilter_init(int n_blocks, int n_framelen, flt_t A, CommonStates* cs);
void KalmanFilter_destroy(KalmanFilter* kf);
void KalmanFilter_update_filter(void* filter, CommonStates* cs);
void KalmanFilter_process_frame(void* filter, vector* out, const vector* y_frame, const vector* x_frame);

#endif //TAOSBEECHDSP_AEC_KALMAN_H
