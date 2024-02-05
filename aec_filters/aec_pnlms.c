//
// Created by wentao on 2024/1/13.
//
#include "stdlib.h"
#include "math.h"

#include "aec_pnlms.h"
#include "common_states.h"
#include "vector.h"
#include "arch.h"
#include "numeric.h"
#include "reduce.h"
#include "dsp_common_func.h"


pNLMSFilter* pNLMSFilter_init(int n_blocks, int n_framelen, CommonStates* common_states) {
    pNLMSFilter* pnlms = (pNLMSFilter*)malloc(sizeof(pNLMSFilter));
    pnlms->mdf = mdf_init(n_blocks, n_framelen);
    MDFFilter* mdf = pnlms->mdf;
    pnlms->cs = common_states;

    SET_error_func(pnlms, mdf->compute_error);
    SET_echo_func(pnlms, mdf->compute_echo);
    SET_update_func(pnlms, pNLMSFilter_update_filter);

    pnlms->X2_smth = zeros(mdf->n_freq, Float);
    pnlms->inv_X2_smth = ones(mdf->n_freq, Float);
    pnlms->prop = zeros(mdf->n_blocks, Float);
    pnlms->ss = FLT(0.35) / n_blocks;

    // Ratio of ~10 between adaptation rate of first and last block */
    flt_t decay = taodsp_exp(FLT(-2.4) / n_blocks);
    vec_elem_set_FLT(pnlms->prop, 0, FLT(0.7));
    flt_t sum_ = FLT(0.7);
    flt_t last = FLT(0.7);
    for (int n=1; n<n_blocks; n++) {
        last *= decay;
        vec_elem_set(pnlms->prop, n, &last);
        sum_ += last;
    }
    // self.prop[:] = 0.8 * self.prop / sum_
    flt_t normalizer = FLT(0.8)/sum_;
    vec_mul_scalar(pnlms->prop, &normalizer, pnlms->prop);

    return pnlms;
}

void pNLMSFilter_destroy(pNLMSFilter* pnlms) {
    if (!pnlms) return;
    mdf_destroy(pnlms->mdf);  // free parent
    BATCH_vec_free(pnlms->X2_smth, pnlms->inv_X2_smth, pnlms->prop);
    free(pnlms);
}

static void adjust_proportion(pNLMSFilter* nlms) {
    /* tmp = 1 + np.sum(self.H_hat * np.conj(self.H_hat), axis=1).real
    *  self.prop[:] = np.sqrt(tmp) */
    const int n_freq = GET_nfreq(nlms);
    vector* prop = nlms->prop;
    vector* H_hat = GET_H_hat(nlms);
    vector H2;
    auto_vector_construct(H2, VEC_SIZE(*H_hat), Float);
    vec_abs2(H_hat, &H2);
    sum2d(&H2, n_freq, prop, AlongCol);

    vec_add_scalar(prop, &flt_one, prop);
    vec_unary_FLT(prop, prop, taodsp_sqrt);
    /* add small portion of max_sum to every entry
     * max_norm = np.max(self.prop) + 1.0 */
    // flt_t tmp_norm = vec_reduce(prop, taodsp_max);
    flt_t tmp_norm;
    vec_max(prop, &tmp_norm);
    tmp_norm += FLT(1.);
    /* prop[:] += 0.1 * max_norm */
    tmp_norm *= FLT(0.1);
    vec_add_scalar(prop, &tmp_norm, prop);

    /* prop_sum = np.sum(self.prop) + 1.0
     * self.prop[:] = (0.99 / prop_sum) * self.prop */
    flt_t normalizer;
    sum(prop, &normalizer);
    normalizer += FLT(1.);
    normalizer = FLT(0.99) / normalizer;
    vec_mul_scalar(prop, &normalizer, prop);
}

void pNLMSFilter_update_filter(void* filter, CommonStates* cs) {
    pNLMSFilter* nlms = (pNLMSFilter*)filter;
    const vector* X = cs->X_fifo;
    const vector* X2 = cs->X2;
    // compute frequency domain recursive average energy of x
    const int N = GET_nfft(nlms);
    const int n_freq = GET_nfreq(nlms);
    flt_t ss = nlms->ss;
    flt_t ss_1 = flt_one - ss;
    flt_t tmp;
    flt_t lr = flt_zero;
    vector* prop = nlms->prop;

    /* X2_smth[:] = (1.0 - self.ss) * self.X2_smth + self.ss * X2 */
    vec_mul_scalar(nlms->X2_smth, &ss_1, nlms->X2_smth);
    vec_mul_scalar_accum(X2, &ss, nlms->X2_smth);
    // compute learning rate
    flt_t See = GET_See(nlms);
    See = MAX(See, 100 * N);
    vector* H_hat = GET_H_hat(nlms);

    // if far-end signal nearly silent, don't adapt
    if (cs->Sxx > FLT(1000.) * N) {
        flt_t tmp2 = MIN(See, cs->Sxx);
        tmp = FLT(0.25) * tmp2;
        lr = tmp / See;
    }

    /* inv_X2_smth[:] = lr / (10 + self.X2_smth) */
    flt_t ten = FLT(10.);
    vec_add_scalar(nlms->X2_smth, &ten, nlms->inv_X2_smth);
    scalar_div_vec(&lr, nlms->inv_X2_smth, nlms->inv_X2_smth);
    adjust_proportion(nlms);

    /*  self.H_hat[:, :] += self.prop[..., np.newaxis] * (self.inv_X2_smth * (np.conj(X[:M, :]) * self.E))
     *  wtmp = np.fft.irfft(self.H_hat, axis=1).real
     *  wtmp[:, F:] = 0.
     *  self.H_hat[:, :] = np.fft.rfft(wtmp, axis=1)
     */
    vector H_delta;
    auto_vector_construct(H_delta, VEC_SIZE(*H_hat), Complex);
    vector_view X_vw;
    vec_view_construct(X, 0, VEC_SIZE(*H_hat), &X_vw);
    vec_conj(&X_vw, &H_delta);
    vec_mul(&H_delta, GET_Error(nlms), &H_delta);
    vec_mul(&H_delta, nlms->inv_X2_smth, &H_delta);

    matrix_view H_delta_mvw, prop_mvw;
    mat_view_from_vec(&H_delta, n_freq, &H_delta_mvw);  // M x n_freq
    mat_view_from_vec(prop, 1, &prop_mvw);  // M x 1
    mat_elemwise_mul(&prop_mvw, &H_delta_mvw, &H_delta_mvw);
    vec_add(&H_delta, H_hat, H_hat);

    matrix_view H_hat_mvw;
    mat_view_from_vec(H_hat, n_freq, &H_hat_mvw);
    linearize_convolution_2d(GET_fft_handler(nlms), &H_hat_mvw);
}
void pNLMSFilter_process_frame(void* filter,
                               vector* out,
                               const vector* y_frame,
                               const vector* x_frame) {
    pNLMSFilter* nlms = (pNLMSFilter*)filter;
    CommonStates* cs = nlms->cs;
    int M = GET_nblocks(nlms);
    int F = GET_nframelen(nlms);
    int n_freq = GET_nfreq(nlms);

    cs->saturated = signal_saturated(y_frame, -32000, 32000);
    vec_dot(x_frame, x_frame, &(cs->Sxx));

    // moving update x, X queue and X2 queue
    vec_left_moving(cs->x, x_frame, F, F, 0);
    small_drft_forward(GET_fft_handler(nlms), cs->x, cs->X);
    vec_right_moving(cs->X_fifo, cs->X, (M-1)*n_freq, n_freq, 0);
    vec_abs2(cs->X, cs->X2);

    nlms->compute_echo(nlms->mdf, cs->X_fifo);
    nlms->compute_error(nlms->mdf, y_frame);
    if (cs->saturated == 0)
        nlms->update_filter(nlms, cs);
    else
        cs->saturated -= 1;

    compute_output(cs, GET_error(nlms), out);
}