//
// Created by wentao on 2024/1/12.
//
#include "aec_kalman.h"
#include "common_states.h"
#include "numeric.h"
#include "arch.h"
#include "smallft.h"
#include "dsp_common_func.h"
#include "vector.h"
#include "reduce.h"
#include "mdf.h"

KalmanFilter* KalmanFilter_init(int n_blocks, int n_framelen, flt_t A, CommonStates* cs) {
    check((0 < A) && (A <= 1.0), "Invalid transition coefficient.");
    KalmanFilter* kf = (KalmanFilter*) malloc(sizeof(KalmanFilter));
    kf->mdf = mdf_init(n_blocks, n_framelen);  // construct parent instance
    kf->cs = cs;
    int n_freq = GET_nfreq(kf);

    kf->A = A;
    kf->A2 = A*A;
    kf->A2_1 = (flt_t)1. - kf->A2;

    int total_bins = n_blocks * n_freq;
//    kf->E2_fifo = zeros(total_bins, Float);
    kf->P = ones(total_bins, Float);
    kf->H2 = zeros(total_bins, Float);
    kf->PHIss = zeros(n_freq, Float);
    kf->PHIee = zeros(total_bins, Float);
    kf->mu = zeros(total_bins, Float);
    kf->H_delta = zeros(total_bins, Complex);

    SET_error_func(kf, kf->mdf->compute_error);
    SET_echo_func(kf, kf->mdf->compute_echo);
    SET_update_func(kf, KalmanFilter_update_filter);

    return kf;
}

void KalmanFilter_destroy(KalmanFilter* kf) {
    if (!kf) return;
    mdf_destroy(kf->mdf);
    BATCH_vec_free(/*kf->E2_fifo, */kf->P, kf->H2, kf->PHIss,
                                    kf->PHIee, kf->mu, kf->H_delta);
    free(kf);
}

//void update_P(vector* P, vector* A2, vector* mu, vector* X2, vector* H2) {
//    for (int i=0; i<P->size; i++) {
//
//    }
//}

void KalmanFilter_update_filter(void* filter, CommonStates* cs) {
    // 这里使用向量运算，不考虑内存和运算的优化
    KalmanFilter* kf = (KalmanFilter*)filter;
    const vector* X = cs->X_fifo;
    const vector* X2 = cs->X2_fifo;
    int n_freq = GET_nfreq(kf);
    const flt_t two = FLT(2.), neg_half = FLT(-0.5);
    vector* H_hat = GET_H_hat(kf);  // convenient reference of vectors
    vector* E = GET_Error(kf);
    vector* P = kf->P;
    vector* H2 = kf->H2;
    vector* PHIss = kf->PHIss;
    vector* PHIee = kf->PHIee;
    vector* mu = kf->mu;
    vector* H_delta = kf->H_delta;
    flt_t A2 = kf->A2;
    flt_t A2_1 = kf->A2_1;
    vector_view X_vw, X2_vw;
    vector X_conj, sum_X2, tmp_P;
    auto_vector_construct(sum_X2, n_freq, Float);
    auto_vector_construct(X_conj, VEC_SIZE(*H2), Complex);
    auto_vector_construct(tmp_P, VEC_SIZE(*P), Float);

    vec_abs2(H_hat, H2);
    vec_abs2(GET_Error(kf), PHIss);

    // np.sum(X2[:M, :], axis=0)
    vec_view_construct(X2, 0, VEC_SIZE(*H2), &X2_vw);  // 取前M个blocks，X2也许包含更多的数据
    sum2d(&X2_vw, n_freq, &sum_X2, AlongRow);
    // PHIee = st.P[:, :] * np.sum(st.X2_fifo[:-1, :], axis=0)
    vec_mul(P, &sum_X2, PHIee);
    // P[:, :] / (PHIee + 2 * PHIss)
    vec_mul_scalar_accum(PHIss, &two, PHIee);
    vec_div(P, PHIee, mu);

    // H_hat[:, :] = H_hat[:, :] + mu * np.conj(X[:M, :]) * self.E
    vec_view_construct(X, 0, VEC_SIZE(*H_hat), &X_vw);
    vec_conj(&X_vw, &X_conj);
    vec_mul(&X_conj, E, H_delta);
    vec_mul(H_delta, mu, H_delta);
    vec_add(H_delta, H_hat, H_hat);

    // P[:, :] = A2 * P[:, :] * (1.0 - 0.5 * mu * X2[:M, :]) + A2_1 * W2[:, :]
    vec_mul(mu, &X2_vw, &tmp_P);
    vec_mul_scalar(&tmp_P, &neg_half, &tmp_P);
    vec_add_scalar(&tmp_P, &flt_one, &tmp_P);  // (1.0 - 0.5 * mu * X2[:M, :])
    vec_mul(P, &tmp_P, P);
    vec_mul_scalar(P, &A2, P);  // A2 * P[:, :] * (1.0 - 0.5 * mu * X2[:M, :])
    vec_mul_scalar_accum(H2, &A2_1, P);

    matrix_view H_hat_mvw;
    mat_view_from_vec(H_hat, n_freq, &H_hat_mvw);
    linearize_convolution_2d(GET_fft_handler(kf), &H_hat_mvw);

//    float PHIee_bi, mu;
//    for (int b=0; b<M; b++) {
//        PHIee_bi = P[b*n_freq] * sum_X2[0];
//        mu = P[b*n_freq] / (PHIee_bi + 2 * PHIss[0]);
//        H_hat[b*N] += mu * X[b*N] * GET_Error(kf)[b*n_freq];
//        P[b*n_freq] = A2 * P[b*n_freq] * (1.f - 0.5*mu*X2[b*n_freq]) + A2_1 * H2[b*n_freq];
//
//        PHIee_bi = MTX_GET(kf->P, b, n_freq-1, n_freq) * sum_X2[n_freq-1];
//        mu = kf->P[b*n_freq] / (PHIee_bi + 2 * PHIss[n_freq-1]);
//        H_hat[(b+1)*N-1] += mu * X[(b+1)*N-1] * GET_Error(kf)[(b+1)*n_freq-1];
//        P[b*n_freq] = A2 * P[b*n_freq] * (1.f - 0.5*mu*X2[b*n_freq]) + A2_1 * H2[b*n_freq];
//
//        for (int i=1; i<n_freq-1; i++) {
//            PHIee_bi = MTX_GET(P, b, i, n_freq) * sum_X2[i];
//            mu = P[b*n_freq+i] / (PHIee_bi + 2 * PHIss[i]);
//            H_hat[b*N+2*i] += mu * conj(X[:M, :]) * GET_Error(kf)[b*n_freq+i];
//        }
//    }
}

void KalmanFilter_process_frame(void* filter,
                               vector* out,
                               const vector* y_frame,
                               const vector* x_frame) {
    KalmanFilter* kf = (KalmanFilter*)filter;
    CommonStates *cs = kf->cs;
    int M = GET_nblocks(kf);
    int F = GET_nframelen(kf);
    int n_freq = GET_nfreq(kf);
    vector_view X_head, X2_head;

    cs->saturated = signal_saturated(y_frame, -32000, 32000);
    vec_left_moving(cs->x, x_frame, F, F, 0);

    vec_right_moving(cs->X_fifo, NULL, (M - 1) * n_freq, n_freq, 0);
    vec_view_construct(cs->X_fifo, 0, n_freq, &X_head);
    small_drft_forward(GET_fft_handler(kf), cs->x, &X_head);

    vec_right_moving(cs->X2_fifo, NULL, (M - 1) * n_freq, n_freq, 0);
    vec_view_construct(cs->X2_fifo, 0, n_freq, &X2_head);
    vec_abs2(&X_head, &X2_head);

    kf->compute_echo(kf->mdf, cs->X_fifo);
    kf->compute_error(kf->mdf, y_frame);
    kf->update_filter(kf, cs);

    compute_output(cs, GET_error(kf), out);
//    compute_leak_coef(df, x_frame);

//    left_moving_vec(cs->last_y, NULL, F, F, 0);
//    vec_view_construct(cs->last_y, F, F, &vw);
//    vec_sub(y_frame, out, &vw);
}