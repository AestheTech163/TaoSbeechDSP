//
// Created by wentao on 2023/12/23.
//
// code from speexdsp
#include "math.h"
#include "stdlib.h"

#include "arch.h"
#include "mmse_noise.h"
#include "dsp_common_func.h"
#include "smallft.h"
#include "aec_dual.h"

flt_t table[21] = {0.82157, 1.02017, 1.20461, 1.37534, 1.53363, 1.68092, 1.81865,
        1.94811, 2.07038, 2.18638, 2.29688, 2.40255, 2.50391, 2.60144,
        2.69551, 2.78647, 2.87458, 2.96015, 3.04333, 3.12431, 3.20326};

void get_conj_window(vector* window) {
    int len = window->size;
    flt_t x, tmp;
    char inv;
    for (int i=0; i<len; i++) {
        x = 4.0 * i / len;
        inv = 0;
        if (x < 1.) ; // do nothing
        else if (x < 2.) {
            x = 2.0 - x;
            inv = 1;
        } else if (x < 3.) {
            x = x - 2.0;
            inv = 1;
        } else
            x = 4 - x;
        x = 1.271903 * x;
        tmp = (0.5 - 0.5 * cos(0.5 * M_PI * x));
        tmp *= tmp;
        if (inv) tmp = 1.0 - tmp;
        vec_elem_set_FLT(window, i, taodsp_sqrt(tmp));
    }
}

MMSENoise* init_MMSENoise(int frame_size, int fs, int n_bands) {
    int n_fft = 2 * frame_size;
    int n_freq = n_fft / 2 + 1;

    MMSENoise* mmse_ns = (MMSENoise*) malloc(sizeof(MMSENoise));
    mmse_ns->frame_size = frame_size;
    mmse_ns->n_freq = n_freq;
    mmse_ns->n_bands = n_bands;
    mmse_ns->n_adapt = 0;
    mmse_ns->min_count = 0;

    mmse_ns->ps = zeros(n_freq, Float);
    mmse_ns->xbands = zeros(n_bands, Float);
    mmse_ns->old_xbands = ones(n_bands, Float);
    mmse_ns->noise = ones(n_freq, Float);
    mmse_ns->noise_bands = ones(n_bands, Float);
//    mmse_ns->speech_prob = ones(N, Float);

    mmse_ns->echo_noise = zeros(n_freq, Float);
    mmse_ns->echo_bands = zeros(n_freq, Float);
    mmse_ns->resid_echo = zeros(n_freq, Float);

    mmse_ns->bands_gain = zeros(n_bands, Float);
    mmse_ns->bands_gain2 = zeros(n_bands, Float);
    mmse_ns->gain_floor = zeros(n_bands, Float);
    mmse_ns->linear_gain = zeros(n_freq, Float);

    mmse_ns->post = ones(n_bands, Float);
    mmse_ns->prior = ones(n_bands, Float);
    mmse_ns->zeta = zeros(n_bands, Float);

    mmse_ns->window = zeros(n_fft, Float);
    mmse_ns->frame = zeros(n_fft, Float);
    mmse_ns->ft = zeros(n_freq, Complex);
    mmse_ns->outbuf = zeros(frame_size, Float);
    mmse_ns->inbuf = zeros(frame_size, Float);

    mmse_ns->S = zeros(n_freq, Float);
    mmse_ns->Stmp = zeros(n_freq, Float);
    mmse_ns->Smin = zeros(n_freq, Float);

    mmse_ns->fft = small_drft_init(n_fft);
    mmse_ns->fbank = init_FilterBank(n_bands, fs, frame_size);

    get_conj_window(mmse_ns->window);
    mmse_ns->noise_floor = taodsp_exp(FLT(0.2302585) * noise_suppress);  // 10**(noise_suppress/10)

    return mmse_ns;
}

void free_MMSENoise(MMSENoise* mmse_ns) {
    if (mmse_ns == NULL) return;
    small_drft_free(mmse_ns->fft);
    free_FilterBank(mmse_ns->fbank);

    BATCH_vec_free(mmse_ns->S, mmse_ns->Stmp, mmse_ns->Smin,
                   mmse_ns->window, mmse_ns->ft,
                   mmse_ns->outbuf, mmse_ns->inbuf, mmse_ns->frame);
    BATCH_vec_free(mmse_ns->bands_gain, mmse_ns->bands_gain2, mmse_ns->linear_gain,
                   mmse_ns->gain_floor, mmse_ns->post, mmse_ns->prior,
                   mmse_ns->zeta, mmse_ns->ps, mmse_ns->xbands, mmse_ns->old_xbands,
                   mmse_ns->noise, mmse_ns->noise_bands, mmse_ns->echo_noise,
                   mmse_ns->echo_bands, mmse_ns->resid_echo);
//    vec_free(mmse_ns->speech_prob);
    free(mmse_ns);
}

void analysis(MMSENoise* nst, vector* x) {
    int N = nst->frame_size;
    SmallFFT* fft = nst->fft;

//    left_moving_overlap(nst->frame, x, N, N);
//    vec_assign(nst->frame, nst->inbuf, N);
//    vec_assign(nst->frame+N, x, N);
    vec_copy_slice(nst->inbuf, 0, N, nst->frame, 0);
    vec_copy_slice(x, 0, N, nst->frame, N);

//    vec_assign(nst->inbuf, x, N);
    vec_copy(x, nst->inbuf);
//    vec_mul(nst->frame, nst->window, nst->frame, 2*N);
    vec_mul(nst->frame, nst->window, nst->frame);
    small_drft_forward(fft, nst->frame, nst->ft);
    vec_abs2(nst->ft, nst->ps);
    psd_to_bank(nst->fbank, nst->ps, nst->xbands);
}

static void update_noise(MMSENoise* nst, flt_t beta) {
    // MCRA 噪声功率估计
    int N = nst->n_freq;
    check((nst->S->stride==1)&&(nst->ps->stride==1), "Need contiguous data layout.");
    flt_t* S = nst->S->data;
    flt_t* ps = nst->ps->data;
    flt_t* Smin = nst->Smin->data;
    flt_t* noise = nst->noise->data;
    int min_range, i;

    // 对当前帧的功率谱进行 时间和频率 的平滑
    for (i=1; i<N-1; i++)
        S[i] = FLT(0.8) * S[i] + FLT(0.05) * ps[i - 1] + FLT(0.1) * ps[i] + FLT(0.05) * ps[i + 1];
    S[0] = FLT(0.8) * S[0] + FLT(0.2) * ps[0];
    S[N - 1] = FLT(0.8) * S[N - 1] + FLT(0.2) * ps[N - 1];
    // st.S[0]   = 0.8*st.S[0] + 0.15*st.ps[0] + 0.05*st.ps[1]
    // st.S[N-1] = 0.8*st.S[N-1] + 0.15*st.ps[N-1] + 0.05*st.ps[N-2]

    if (nst->n_adapt == 1) {
        vec_fill(nst->Smin, &flt_zero);
        vec_fill(nst->Stmp, &flt_zero);
    }

    if (nst->n_adapt < 100) min_range = 15;
    else if (nst->n_adapt < 1000) min_range = 50;
    else if (nst->n_adapt < 10000) min_range = 150;
    else min_range = 300;

    // 对每个频点单独统计噪声功率,该算法能确保 Smin 在 2*min_range 后一定得到更新，
    // 并且 Smin 至少是 min_range 长的窗口内的统计值
    if (nst->min_count > min_range) {
        nst->min_count = 0;
//        for (i = 0; i < N; i++) {
//            nst->Smin[i] = MIN(nst->Stmp[i], nst->S[i]);
//            nst->Stmp[i] = nst->S[i];
//        }
        vec_minimum(nst->Stmp, nst->S, nst->Smin);
        vec_copy(nst->S, nst->Stmp);
    } else {
//        for (i=0; i<N; i++) {
//            nst->Smin[i] = MIN(nst->Smin[i], nst->S[i]);
//            nst->Stmp[i] = MIN(nst->Stmp[i], nst->S[i]);
//        }
        vec_minimum(nst->Smin, nst->S, nst->Smin);
        vec_minimum(nst->Stmp, nst->S, nst->Stmp);
    }
    for (i=0; i<N; i++) {
        if ((0.4 * S[i] < Smin[i]) ||  // 指示没有语音存在
            (ps[i] < noise[i])) {
            flt_t avg = (1-beta) * noise[i] + beta * ps[i];
            noise[i] = MAX(0, avg);  // 更新噪声功率统计值
        }
    }
    psd_to_bank(nst->fbank, nst->noise, nst->noise_bands);
}

static flt_t hypergeom_gain(flt_t x) {
    int integer = (int)(2 * x);
    if (integer < 0) return FLT(1.);
    if (integer > 19) return FLT(1.) + 0.1296 / x;
    flt_t frac = 2 * x - (flt_t)integer;
    return ((FLT(1.) - frac) * table[integer] + frac * table[integer + 1]) / taodsp_sqrt(x + FLT(.0001));
}

static inline flt_t qcurve(float x) {
    return FLT(1.0) / (1.0 + 0.15 / x);
}

static void compute_gain_floor(MMSENoise* nst, flt_t effective_echo_suppress) {
    flt_t echo_floor = taodsp_exp(FLT(0.2302585) * effective_echo_suppress);
//    for (int i=0; i<nst->n_bands; i++) {
//        flt_t floor_en = nst->noise_floor * nst->noise_bands[i] + echo_floor * nst->echo_bands[i];
//        nst->gain_floor[i] = sqrt(floor_en) / sqrt(1 + nst->noise_bands[i] + nst->echo_bands[i]);
//    }
    vector floor_en, mix_en;
    auto_vector_construct(floor_en, nst->n_bands, Float);
    auto_vector_construct(mix_en, nst->n_bands, Float);

    vec_mul_scalar(nst->noise_bands, &nst->noise_floor, &floor_en);
    vec_mul_scalar_accum(nst->echo_bands, &echo_floor, &floor_en);
    vec_sqrt(&floor_en, &floor_en);
    vec_add(nst->noise_bands, nst->echo_bands, &mix_en);
    vec_add_scalar(&mix_en, &flt_one, &mix_en);
    vec_sqrt(&mix_en, &mix_en);
    vec_div(&floor_en, &mix_en, nst->gain_floor);
}

void mmsenoise_process_frame(MMSENoise* nst, vector* x, void* aec_obj) {
    int i;
    int N = nst->frame_size;
    int M = nst->n_bands;
    nst->n_adapt += 1;
    nst->n_adapt = MIN(nst->n_adapt, 20000);
    nst->min_count += 1;
    flt_t alpha = FLT(0.6), tot_noise, gamma, last_g, beta = MAX(0.03, 1.0 / nst->n_adapt);

    flt_t* noise_bands = nst->noise_bands->data;
    flt_t* echo_bands = nst->echo_bands->data;
    flt_t* xbands = nst->xbands->data;
    flt_t* post = nst->post->data;
    flt_t* old_xbands = nst->old_xbands->data;
    flt_t* prior = nst->prior->data;
    flt_t* bands_gain = nst->bands_gain->data;
    flt_t* zeta = nst->zeta->data;
    flt_t* bands_gain2 = nst->bands_gain2->data;
    flt_t* gain_floor = nst->gain_floor->data;

    if (aec_obj) {
        ((DualFilter*)aec_obj)->get_residual(aec_obj, nst->resid_echo);
        vec_mul_scalar(nst->echo_noise, &alpha, nst->echo_noise);
        vec_maximum(nst->echo_noise, nst->resid_echo, nst->echo_noise);
//        for (i=0; i<N; i++)
//            nst->echo_noise[i] = MAX(0.6 * nst->echo_noise[i], nst->resid_echo[i]);  // 残留回声功率 的 递归平滑
        psd_to_bank(nst->fbank, nst->echo_noise, nst->echo_bands);  //
    } else {
        vec_fill(nst->echo_noise, &flt_zero);
        vec_fill(nst->echo_bands, &flt_zero);
    }

    analysis(nst, x);
    update_noise(nst, beta);

    if (nst->n_adapt == 1)
        vec_copy(nst->xbands, nst->old_xbands);

    for (i=0; i<M; i++) {
        // Total noise estimate including residual echo and reverberation
        tot_noise = 1 + noise_bands[i] + echo_bands[i];  // + nst->reverb_estimate[N+i];
        // A posteriori SNR = ps/noise - 1
        post[i] = xbands[i] / tot_noise - FLT(1.0);
        post[i] = CLIP(post[i], FLT(0.), FLT(100.));
        // Computing update gamma = 0.1 + 0.9*(old/(old+noise))^2
        last_g = old_xbands[i] / (old_xbands[i] + tot_noise);
        gamma = FLT(0.1) + FLT(0.89) * SQUARE(last_g);
        // A priori SNR update = gamma*max(0,post) + (1-gamma)*old/noise
        prior[i] = gamma * MAX(FLT(0.), post[i]) + (FLT(1.0)-gamma) * old_xbands[i] / tot_noise;
        prior[i] = MIN(prior[i], 100.);
    }
    flt_t alpha_z = FLT(0.7), alpha_z_1 = FLT(0.3);
    vec_mul_scalar(nst->zeta, &alpha_z, nst->zeta);
    vec_mul_scalar_accum(nst->prior, &alpha_z_1, nst->zeta);
    // Speech probability of presence for the entire frame is based on the average filterbank a priori SNR
    flt_t Zframe;
    mean(nst->zeta, &Zframe);
    flt_t Pframe = FLT(0.1) + 0.899 * qcurve(Zframe);

    // 在 -40 和 -15 之间做插值，语音存在的概率 越大，抑制程度 越低
    // Compute the gain floor based on different floors for the background noise and residual echo */
    float effective_echo_suppress = (FLT(1.) - Pframe) * echo_suppress + Pframe * echo_suppress_active;
    compute_gain_floor(nst, effective_echo_suppress);

    // Compute Ephraim & Malah gain, speech probability of presence for each critical band (Bark scale)
    for (i=0; i<M; i++) {
        flt_t prior_ratio = prior[i] / (prior[i] + 1);
        flt_t theta = prior_ratio * (FLT(1.) + post[i]);
        // TODO post在这里没有取 max(0,post), check
        flt_t MM = hypergeom_gain(theta);
        /* Gain with bound */
        bands_gain[i] = MIN(FLT(1.), prior_ratio * MM);
        /* Save old Bark power spectrum */
        old_xbands[i] = FLT(0.2) * old_xbands[i] + FLT(0.8) * SQUARE(bands_gain[i]) * xbands[i];
        flt_t P1 = FLT(0.199) + FLT(0.8) * qcurve(zeta[i]);  // 这里的 zeta 是 zeta_local, valin 没有使用 zeta_global 选项
        flt_t q = FLT(1.) - Pframe * P1;  // Pframe、P1 分别 是 全局、局部 语音存在的先验概率，因此 q 就是 语音不存在的先验概率
        bands_gain2[i] = FLT(1.) / (1. + q / (1. - q) * (1 + prior[i]) * exp(-theta));  // 给定当前帧的第 i 个频点的观测值下，语音存在的条件概率
    }
    for (i=0; i<M; i++) {
        flt_t p = bands_gain2[i];
        bands_gain[i] = MAX(bands_gain[i], gain_floor[i]);
        flt_t tmp = p * taodsp_sqrt(bands_gain[i]) + (FLT(1.)-p) * taodsp_sqrt(gain_floor[i]);
        bands_gain2[i] = tmp * tmp;
    }
    bank_to_psd(nst->fbank, nst->bands_gain2, nst->linear_gain);

    // Apply computed gain
//    for (i=1; i<N; i++) {
//        nst->ft[2*i-1] = nst->linear_gain[i] * nst->ft[2*i-1];
//        nst->ft[2*i] = nst->linear_gain[i] * nst->ft[2*i];
//    }
//    nst->ft[0] = nst->linear_gain[0] * nst->ft[0];
//    nst->ft[2*N-1] = nst->linear_gain[N-1] * nst->ft[2*N-1];

    vec_mul(nst->ft, nst->linear_gain, nst->ft);

    small_drft_backward(nst->fft, nst->ft, nst->frame);
    //self.frame[:] = self.frame * self.window
    vec_mul(nst->frame, nst->window, nst->frame);
//    vec_mul_scalar(nst->ft, 1./(2*N), nst->ft, 2*N);

//    for (i=0; i<N; i++) {
//        flt_t tmp = nst->outbuf[i] + nst->frame[i];
//        x[i] = CLIP(tmp, -32768, 32767);
//    }
    vector_view frame_vw;
    vec_view_construct(nst->frame, 0, N, &frame_vw);
    vec_add(nst->outbuf, &frame_vw, x);
    vec_clip(x, -32768, 32767);

//    vec_assign(nst->outbuf, nst->frame+N, N);
    vec_copy_slice(nst->frame, N, N, nst->outbuf, 0);
}