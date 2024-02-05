/* Dual filter for acoustic echo cancellation
 * filter 1 : multi-delay proportional-NLMS
 *            H[n+1] = H[n] + prop * mu * E[n-1] * conj(X[n-1]) / ema_X2
 *            with pre-emphsize & notch dc & proportional mu
 * filter 2 : multi-delay diagonal-Kalman-filter
 */
#include "stdlib.h"
#include "math.h"

#include "arch.h"
#include "common_states.h"
#include "aec_dual.h"
#include "smallft.h"
#include "dsp_common_func.h"
#include "vector.h"
#include "aec_pnlms.h"
#include "aec_kalman.h"


void DualFilter_get_residual(void* aec_obj, vector* resid);

DualFilter* DualFilter_init(int fg_blocks, int bg_blocks, int n_framelen, flt_t A, CommonStates* cs) {
    DualFilter* df = (DualFilter*)malloc(sizeof(DualFilter));

    df->fg_blocks = fg_blocks;
    df->bg_blocks = bg_blocks;
    int n_blocks = MAX(df->fg_blocks, df->bg_blocks);

    df->cs = cs;
    df->fg = KalmanFilter_init(df->fg_blocks, n_framelen, A, cs);
    df->bg = pNLMSFilter_init(df->bg_blocks, n_framelen, cs);
    int n_freq = GET_nfreq(df->fg);
    int n_fft = GET_nfft(df->fg);

    df->window = zeros(n_fft, Float);
    df->saturated = 0;
    df->bg_better_cnt = 0;
    df->bg_update_cnt = 0;
    df->update_fg = 0;

    for (int n=0; n<n_fft; n++)
        vec_elem_set_FLT(df->window, n, FLT(0.5)-FLT(0.5)*cos(2*M_PI*n/n_fft));
    vec_view_construct(df->window, n_framelen, n_framelen, &(df->left_wind));
    vec_view_construct(df->window, 0, n_framelen, &(df->right_wind));

    df->Davg1 = FLT(0.);
    df->Dvar1 = FLT(0.);
    df->ss = FLT(0.35) / n_blocks;
    df->sum_adapt = FLT(0.);
    df->adapted = 0;

    df->get_residual = DualFilter_get_residual;

    return df;
}

void DualFilter_destroy(DualFilter* df) {
    if (!df) return;
    KalmanFilter_destroy(df->fg);
    pNLMSFilter_destroy(df->bg);
    vec_free(df->window);
    free(df);
}

void compare_filter(DualFilter* df) {
    float Dbf;
    vector_view vw1, vw2;
    int F = GET_nframelen(df->fg);
    int n_freq = GET_nfreq(df->fg);
    int n_fft = GET_nfft(df->fg);
    df->update_fg = 0;

    if (GET_See(df->fg) * 0.5 > GET_See(df->bg)) df->bg_better_cnt += 1;
    else df->bg_better_cnt = 0;
    df->bg_update_cnt += 1;

    if ((df->bg_better_cnt >= 5) && (df->Sxx > n_fft * 1000)) {
        printf("update fg");
        df->update_fg = 1;
        df->bg_better_cnt = 0;
    }
//# elif (self.bg_update_cnt >= 20) and (self.fg.See * 1.5 < self.bg.See):
//#     self.bg_update_cnt = 0
//#     print("update bg")
//#     self.bg.H_hat[:, :] = self.fg.H_hat[:self.bg.n_blocks, :]
//#     self.bg.E[:] = self.fg.E

    vector e_diff;
    vector* fg_e = GET_error(df->fg);
    vector* bg_e = GET_error(df->bg);
    auto_vector_construct(e_diff, VEC_SIZE(*fg_e), Float);
    vec_sub(fg_e, bg_e, &e_diff);
    vec_dot(&e_diff, &e_diff, &Dbf);
    Dbf += 10;
    df->Davg1 = FLT(0.6) * df->Davg1 + FLT(0.4) * (GET_See(df->fg) - GET_See(df->bg));
    df->Dvar1 = FLT(0.16) * df->Dvar1 + FLT(0.36) * GET_See(df->fg) * Dbf;

    flt_t See_diff = GET_See(df->fg) - GET_See(df->bg);
    flt_t cond1 = See_diff * taodsp_abs(See_diff);
    flt_t thres1 = GET_See(df->fg) * Dbf;
    if (cond1 > thres1 * 1.0) {
        df->update_fg = 1;
    }

    if (df->Davg1 * taodsp_abs(df->Davg1) > FLT(1.) * df->Dvar1) {
        df->update_fg = 1;
    }

    if (df->update_fg) {
        int fg_H_size = VEC_SIZE(*GET_H_hat(df->fg));
        int bg_H_size = VEC_SIZE(*GET_H_hat(df->bg));
        int copy_len = MIN(fg_H_size, bg_H_size);
        int remain_len = fg_H_size - bg_H_size;
        vec_copy_slice(GET_H_hat(df->bg), 0, copy_len, GET_H_hat(df->fg), 0);
        if (remain_len > 0) {
            vec_view_construct(GET_H_hat(df->fg), copy_len, remain_len, &vw1);
            vec_fill(&vw1, &flt_zero);
        }
        // self.fg.e[F:] = self.window[F:] * self.fg.e[F:] + self.window[:F] * self.bg.e[F:]
        vec_view_construct(fg_e, F, F, &vw1);
        vec_mul(&(df->left_wind), &vw1, &vw1);
        vec_view_construct(bg_e, F, F, &vw2);
        vec_mul_accum(&(df->right_wind), &vw2, &vw1);
    }
}
//
//void compute_output(DualFilter* df, vector* out, int depre) {
//    int F = GET_nframelen(df->fg);
//    vector* error = GET_error(df->fg);
//    if (depre) {
//        T preemph = df->cs->preemph;
//        float tmp;
//        for (int n=0; n<F; n++) {
//            tmp = vec_get_real(error, F + n);
//            tmp += preemph * df->cs->depre_mem_o;
//            vec_set_real(out, n, tmp);
//            df->cs->depre_mem_o = tmp;
//        }
//    } else {
//        vec_copy_slice(error, F, F, out);
//    }
//}

void compute_leak_coef(DualFilter* df, const vector* x_data) {
    // from speex
    CommonStates* cs = df->cs;
    vector_view vw1;
    flt_t ss_1 = FLT(1.) - df->ss;
    flt_t Syy, tmp_Pey, tmp_Pyy;
    int F = GET_nframelen(df->fg);

    vec_abs2(GET_Error(df->fg), cs->Rf);
    vec_abs2(GET_Y_hat(df->fg), cs->Yf);
    // Sey = np.dot(self.fg.e[F:], self.fg.y_hat[F:])
    vec_view_construct(GET_y_hat(df->fg), F, F, &vw1);
    vec_dot(&vw1, &vw1, &Syy);
    // Sxx = np.dot(x_data, x_data)

    // Compute filtered spectra and (cross-)correlations
    vector tmp_Eh, tmp_Yh;
    auto_vector_construct(tmp_Eh, VEC_SIZE(*cs->Rf), Float);
    auto_vector_construct(tmp_Yh, VEC_SIZE(*cs->Rf), Float);
    vec_sub(cs->Rf, cs->Eh, &tmp_Eh);
    vec_sub(cs->Yf, cs->Yh, &tmp_Yh);
    vec_dot(&tmp_Eh, &tmp_Yh, &tmp_Pey);
    vec_dot(&tmp_Yh, &tmp_Yh, &tmp_Pyy);
    tmp_Pey += 1;
    tmp_Pyy += 1;

    // Eh = (1.-ss) * Eh + ss * Rf
//    vec_mul_scalar_add(df->Eh, 1.f - df->ss, df->Rf, df->ss, df->Eh);
    vec_mul_scalar(cs->Eh, &ss_1, cs->Eh);
    vec_mul_scalar_accum(cs->Rf, &df->ss, cs->Eh);
    // Yh = (1.-ss) * Yh + ss * Yf
//    vec_mul_scalar_add(df->Yh, 1.f - df->ss, df->Yf, df->ss, df->Yh);
    vec_mul_scalar(cs->Yh, &ss_1, cs->Yh);
    vec_mul_scalar_accum(cs->Yf, &df->ss, cs->Yh);

    tmp_Pyy = taodsp_sqrt(tmp_Pyy);
    tmp_Pey = tmp_Pey / tmp_Pyy;  // 使用偏相关系数

    /* Compute correlation update rate */
    flt_t tmp = cs->beta0 * Syy;  // beta0是个很小的值：128*2/8000=0.032
    tmp = MIN(tmp, cs->beta_max * GET_See(df->fg));
    flt_t alpha = tmp / GET_See(df->fg);
    flt_t alpha_1 = FLT(1.) - alpha;

    /* Update correlations (recursive average) */
    cs->Pey = alpha_1 * cs->Pey + alpha * tmp_Pey;
    cs->Pyy = alpha_1 * cs->Pyy + alpha * tmp_Pyy;
    cs->Pyy = MAX(flt_one, cs->Pyy);
    /* We don't really hope to get better than 33 dB (MIN_LEAK-3dB) attenuation anyway*/
    if (cs->Pey < cs->min_leak * cs->Pyy)
        cs->Pey = cs->min_leak * cs->Pyy;
    if (cs->Pey > cs->Pyy)
        cs->Pey = cs->Pyy;
    /* leak_estimate is the linear regression result */
    cs->leak_estimate = cs->Pey / cs->Pyy;
}

void DualFilter_get_residual(void* aec_obj, vector* resid) {
    DualFilter* df = aec_obj;
    vec_mul(df->window, df->cs->last_y, df->cs->y);
    small_drft_forward(GET_fft_handler(df->fg), df->cs->y, df->cs->Y);
    vec_abs2(df->cs->Y, resid);
    flt_t leak2 = MIN(1., 2*df->cs->leak_estimate);
    vec_mul_scalar(resid, &leak2, resid);
}

void DualFilter_process_frame(DualFilter* df,
                              vector* out,
                              const vector* y_frame,
                              const vector* x_frame) {
    int M = GET_nblocks(df->fg);
    int F = GET_nframelen(df->fg);
    int n_freq = GET_nfreq(df->fg);
    KalmanFilter* fg = df->fg;
    pNLMSFilter* bg = df->bg;
    CommonStates* cs = df->cs;
    vector_view vw;

    df->saturated = signal_saturated(y_frame, -32000, 32000);
    vec_dot(x_frame, x_frame, &(cs->Sxx));
    // moving update x, X queue and X2 queue
    vec_left_moving(cs->x, x_frame, F, F, 0);
    vec_right_moving(cs->X_fifo, NULL, (M - 1) * n_freq, n_freq, 0);
    vector_view X_head, X2_head;
    vec_view_construct(cs->X_fifo, 0, n_freq, &X_head);
    small_drft_forward(GET_fft_handler(fg), cs->x, &X_head);

    vec_right_moving(cs->X2_fifo, NULL, (M - 1) * n_freq, n_freq, 0);
    vec_view_construct(cs->X2_fifo, 0, n_freq, &X2_head);
    vec_abs2(&X_head, &X2_head);

    fg->compute_echo(fg->mdf, cs->X_fifo);
    fg->compute_error(fg->mdf, y_frame);
    bg->compute_echo(bg->mdf, cs->X_fifo);
    bg->compute_error(bg->mdf, y_frame);
    if (df->saturated == 0) {
        fg->update_filter(fg, cs);
        bg->update_filter(bg, cs);
    } else
        df->saturated -= 1;

    compare_filter(df);

    compute_output(cs, GET_error(df->fg), out);

    compute_leak_coef(df, x_frame);

    vec_left_moving(cs->last_y, NULL, F, F, 0);
    vec_view_construct(cs->last_y, F, F, &vw);
    vec_sub(y_frame, out, &vw);
}