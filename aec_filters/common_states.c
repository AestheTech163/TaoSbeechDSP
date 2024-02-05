//
// Created by wentao on 2024/1/19.
//

#include "common_states.h"

CommonStates* CommonStates_init(int n_blocks, int n_framelen,
                                int n_freq, int n_fft, int sample_rate,
                                flt_t preemph) {

    CommonStates* cs = (CommonStates*)malloc(sizeof(CommonStates));

    cs->saturated = 0;
    cs->preemph = preemph;

    cs->x = zeros(n_fft, Float);
    cs->y = zeros(n_fft, Float);
    cs->X = zeros(n_freq, Complex);
    cs->Y = zeros(n_freq, Complex);

    cs->X_fifo = zeros((n_blocks+0)*n_freq, Complex);
    cs->X2_fifo = zeros((n_blocks+0)*n_freq, Float);
    cs->X2 = zeros(n_freq, Float);

    cs->last_y = zeros(n_fft, Float);
    cs->min_leak = FLT(0.005);
    cs->Yf = zeros(n_freq, Float);
    cs->Rf = zeros(n_freq, Float);
    cs->Xf = zeros(n_freq, Float);
    cs->Eh = zeros(n_freq, Float);
    cs->Yh = zeros(n_freq, Float);

    cs->Sxx = flt_zero;
    cs->Pey = flt_one;
    cs->Pyy = flt_one;
    cs->beta0 = (FLT(2.) * n_framelen) / sample_rate;
    cs->beta_max = (FLT(.5) * n_framelen) / sample_rate;
    cs->leak_estimate = flt_zero;

    return cs;
}

void CommonStates_destroy(CommonStates* cs) {
    if (!cs) return;
    BATCH_vec_free(cs->x, cs->y, cs->last_y,
                   cs->Y, cs->X, cs->X2,
                   cs->X_fifo, cs->X2_fifo,
                   cs->Yf, cs->Rf, cs->Xf, cs->Eh, cs->Yh);
    free(cs);
}

void compute_output(CommonStates* cs, vector* error, vector* out) {
    int F = out->size;
    if (cs->preemph > 0) {
        flt_t tmp, preemph = cs->preemph;
        for (int n=0; n<F; n++) {
            tmp = vec_elem_get_FLT(error, F + n);
            tmp += preemph * cs->depre_mem_o;
            vec_elem_set_FLT(out, n, tmp);
            cs->depre_mem_o = tmp;
        }
    } else {
        vec_copy_slice(error, F, F, out, 0);
    }
}
