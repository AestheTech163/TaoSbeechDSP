//
// Created by wentao on 2023/12/20.
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "vector.h"
#include "arch.h"
#include "lpc.h"
#include "f0tracker.h"
#include "smallft.h"

ShabbySwapQueue* init_queue(int num_cand_per_list) {
    ShabbySwapQueue* queue = (ShabbySwapQueue*) malloc(sizeof(ShabbySwapQueue));
    check(queue, "malloc ShabbySwapQueue failed.");
    queue->list_len = num_cand_per_list;
    queue->this_len = 0;
    queue->that_len = 0;
    queue->this = queue->this_;
    queue->that = queue->that_;
//    queue->that_frame_index = -100;
//    queue->this_frame_index = 0;
    return queue;
}

void free_queue(ShabbySwapQueue* queue) {
    if (!queue) return;
    free(queue);
}

void swap_queue(ShabbySwapQueue* queue) {
    PeriodCandidate* tmpcand = queue->this;
    queue->this = queue->that;
    queue->that = tmpcand;
    queue->that_len = queue->this_len;
    queue->this_len = 0;
}

int this_empty(ShabbySwapQueue* queue) {
    return queue->this_len == 0;
}

void enqueue(ShabbySwapQueue* queue, char voiced,
             flt_t ncc_value, int ncc_period, flt_t local_cost) {
    // 这里就不做边界检查了，可以确保数量不超过 list_len
    PeriodCandidate* cand = (PeriodCandidate*)(queue->this + queue->this_len);
    cand->voiced = voiced;
    cand->ncc_value = ncc_value;
    cand->ncc_period = ncc_period;
    cand->local_cost = local_cost;
    queue->this_len += 1;
}

PitchState* init_pitch(int fs, int steplen, int windowlen,
                       flt_t rms_min_freq, flt_t rms_max_freq) {
    PitchState* pst = (PitchState*)malloc(sizeof(PitchState));
    check(pst, "malloc PitchState failed.");

    pst->rms_min_val = FLT(20.);
    pst->rms_max_val = FLT(50.);

    pst->in_buflen = steplen * 4;
    pst->in_buf = zeros(pst->in_buflen, Float);

    pst->steplen = steplen;
    pst->lpc_buflen = steplen * 4;
    pst->lpc_resid = zeros(pst->lpc_buflen, Float);

    int n_norm_frames_per_frame = 5;  // norm needs finer granularity
    pst->norm_steplen = steplen / n_norm_frames_per_frame;
    check((steplen%n_norm_frames_per_frame)==0, "steplen % n_norm_frames_per_frame should be zero.");
    pst->norm_buflen = steplen * 3;
    pst->norm_resid = zeros(pst->norm_buflen, Float);

    pst->fft = small_drft_init(steplen);

    int lpc_ord = get_lpc_order(fs);
    pst->lpc = init_lpc(lpc_ord, steplen, windowlen,
                        lpc_preemp_, lpc_noise_floor_);

    int n_future = 1;
    int n_past = (int)(time_span_ / internal_frame_interval_ + 0.5) - n_future;
    int rms_queue_len = n_future + n_past + 1;
    pst->rms_queue = zeros(rms_queue_len, Float);
    pst->rms_first_bin = roundup(steplen * rms_min_freq / fs);
    pst->rms_last_bin = roundup(steplen * rms_max_freq / fs);
    pst->ft_window = zeros(steplen, Float);
    flt_t arg = M_PI * FLT(2.) / steplen;
    for (int i = 0; i < steplen; i++)
        vec_elem_set_FLT(pst->ft_window, i, FLT(0.5) - FLT(0.5) * cos((FLT(0.5) + i) * arg));

    pst->peak_ind = roundup(peak_delay_ * fs);
    pst->skew_ind = roundup(skew_delay_ * fs);

    pst->corr_windowlen = (int)(fs * corr_window_dur_);  // 120;
    pst->first_nccf_lag = (int)(fs / max_f0_search_ + 0.5);  // 32
    pst->max_lag = (int)(fs / min_f0_search_ + 0.5);  // 200
    pst->n_nccf_lags = pst->max_lag - pst->first_nccf_lag;

    pst->lattice = init_queue(max_peaks_num_);

    return pst;
}

void free_pitch(PitchState* pst) {
    if (!pst) return;
    BATCH_vec_free(pst->in_buf, pst->rms_queue, pst->lpc_resid, pst->norm_resid, pst->ft_window);
    small_drft_free(pst->fft);
    free_lpc(pst->lpc);
    free_queue(pst->lattice);
    free(pst);
}

flt_t get_rms_voicing(PitchState* pst, flt_t frame_rms) {
    flt_t range_, prob_voiced, voiced_weight;
    if (frame_rms < pst->rms_min_val)
        pst->rms_min_val = frame_rms;
    else if (frame_rms > pst->rms_max_val)
        pst->rms_max_val = frame_rms;
    range_ = pst->rms_max_val - pst->rms_min_val;
    if (range_ < FLT(1.)) range_ = FLT(1.);
    prob_voiced = MAX((frame_rms - pst->rms_min_val) / range_, FLT(0.));
    voiced_weight = MAX((frame_rms - 40) / 60, FLT(0.));   // TODO
    return prob_voiced; //* voiced_weight;
}

void get_voice_transition(flt_t rms_new, flt_t rms_old, flt_t* prob_onset, flt_t* prob_offset) {
//    if first_frame:
//        return 0., 0.
    flt_t delta_rms = (rms_new - rms_old) / FLT(30.);
    if (delta_rms < 0.0) {
        *prob_onset = FLT(0.);
        *prob_offset = MIN(-delta_rms, FLT(1.));
    } else {
        *prob_onset = MIN(delta_rms, FLT(1.));
        *prob_offset = FLT(0.);
    }
}

void normalize_by_rms(vector* in, vector* out, flt_t* old_inv_rms) {
    flt_t* out_data = (flt_t*)(out->data);
    const flt_t* in_data = (flt_t*)(in->data);
    flt_t ref_energy, inv_rms, delta_inv_rms;
    int N = in->size;
    vec_dot(in, in, &ref_energy);
    ref_energy += FLT(1.);  // to prevent divz
    inv_rms = taodsp_sqrt((flt_t)N / ref_energy);
    delta_inv_rms = (inv_rms - *old_inv_rms) / N;
    for (int i=0; i<N; i++) {
        out_data[i] = in_data[i] * (*old_inv_rms);
        *old_inv_rms += delta_inv_rms;
    }
}

//float get_peak_quality(PitchState* pst, float* lpc_resid, int from, int to) {
//    float val;
//    from -= pst->peak_ind;
//    to -= pst->peak_ind;
//    for (int i=from; i<to; i++) {
//        val = lpc_resid[i];
//        if val >
//    }
//}

void norm_cross_corr(flt_t* ncc, const flt_t* data, int start, int first_lag, int n_lags, int cc_size) {
    flt_t cc, energy, lag_energy;
    const flt_t* input = data + start;
    arr_dot_FLT(input, input, cc_size, 1, 1, &energy);
    if (energy == FLT(0.)) return;
    arr_dot_FLT(input+first_lag, input+first_lag, cc_size, 1, 1, &lag_energy);
    int last_lag = first_lag + n_lags;
    for (int lag=first_lag, j = 0; lag<last_lag; lag++, j++) {
        arr_dot_FLT(input, input+lag, cc_size, 1, 1, &cc);
        if (lag_energy < FLT(1.)) lag_energy = FLT(1.);
        ncc[j] = cc / taodsp_sqrt(energy * lag_energy + FLT(1.));
        lag_energy = lag_energy - SQUARE(input[lag]) + SQUARE(input[lag+cc_size]);
    }
}

void find_ncc_peaks(int* peak_ind, flt_t* peak_val, int* n_peaks, int max_peaks,
                    const flt_t* ncc, int N, flt_t thresh) {
    *n_peaks = 0;
    for (int i=1, o=0; i<N-1; i++) {
        flt_t val = ncc[i];
        if ((val > thresh) && (val > ncc[i-1]) && (val > ncc[i+1])) {
            *n_peaks += 1;
            peak_ind[o] = i;
            peak_val[o] = val;
            o++;
            if (*n_peaks >= max_peaks) return;
        }
    }
}

void construct_lattice(ShabbySwapQueue* lattice,
                       const int* peaks_ind, const flt_t* peaks_val, int num_peaks,
                       flt_t prob_voiced, int ncc_lags) {
    flt_t max_ncc = FLT(0.);
    int max_period = 0;
    flt_t lowest_cost = FLT(1e30);
    flt_t uv_local_cost = FLT(0.);
    flt_t level_cost, period_cost, local_cost;

    if (num_peaks <= 0) return;
    swap_queue(lattice);

    for (int i=0; i<num_peaks; i++) {
        flt_t cc_value = peaks_val[i];
        if (cc_value <= 0) continue;

        level_cost = FLT(1.8) * (1.0 - prob_voiced);
        period_cost = FLT(0.002) * peaks_ind[i];  //  周期越大，惩罚越大
        // peak_qual_cost = 1.9 / (peak_quality)  // 1.9 / (peak_quality)
        local_cost = FLT(2.) * (FLT(1.) - cc_value) + level_cost + period_cost + reward_;
        enqueue(lattice, 1, cc_value, peaks_ind[i]+ncc_lags, local_cost);

        if (local_cost < lowest_cost) {
            lowest_cost = local_cost;
            max_ncc = cc_value;
            uv_local_cost = FLT(4.9) * cc_value + level_cost + FLT(0.9) + reward_;
            max_period = peaks_ind[i] + ncc_lags;
        }
    }
    if (!this_empty(lattice))
        enqueue(lattice, 0, max_ncc, max_period, uv_local_cost);
}

void dynamic_programming(ShabbySwapQueue* lattice, flt_t prob_onset, flt_t prob_offset) {
    /* TODO 间隔一定距离，从头开始计算？？
    if len(peak_queue) == 1:
    for p in peak_queue[0]:
        p.accum_cost = p.local_cost
        p.best_prev = -1
    */
    flt_t min_cost, trans_cost, sum_cost;
    int min_index;

//    if (first_frame) {
//        for (p in peak_queue[0])
//        p.accum_cost = p.local_cost
//        p.best_prev = -1
//    }

    for (int i=0; i<lattice->this_len; i++) {
        min_cost = FLT(1e30);
        min_index = 0;
        PeriodCandidate *that_cand, *this_cand = &(lattice->this[i]);
        if (lattice->that_len <= 0) break;
        for (int j=0; j<lattice->that_len; j++) {
            that_cand = &(lattice->that[j]);
            if (this_cand->voiced && that_cand->voiced) {
                trans_cost = FLT(1.9) * taodsp_abs(taodsp_log((flt_t)this_cand->ncc_period/(flt_t)that_cand->ncc_period));
                sum_cost = trans_cost + that_cand->accum_cost;
            } else if (this_cand->voiced) {
                trans_cost = FLT(0.4) * (FLT(1.) - prob_onset);
                sum_cost = that_cand->accum_cost + trans_cost;
            } else if (that_cand->voiced) {
                trans_cost = FLT(0.4) * (FLT(1.) - prob_offset);
                sum_cost = that_cand->accum_cost + trans_cost;
            } else {
                sum_cost = that_cand->accum_cost;
            }
            if (sum_cost < min_cost) {
                min_cost = sum_cost;
                min_index = j;
            }
        }  // for j
        this_cand->accum_cost = this_cand->local_cost + min_cost;
        this_cand->best_prev = min_index;
    }  // for i
}

void backtrack(ShabbySwapQueue* lattice, int* n0, flt_t* corr) {
    *n0 = 0;
    *corr = FLT(0.);
    flt_t min_cost = FLT(1.0e30);
    int min_ind = 0, is_voiced = 0;
    PeriodCandidate* end = lattice->this;
    // 目前，直接使用最后一帧的结果
    for (int i=0; i<lattice->this_len; i++) {
        if (end[i].accum_cost < min_cost) {
            min_cost = end[i].accum_cost;
            min_ind = i;
        }
    }
    *n0 = end[min_ind].ncc_period;
    *corr = end[min_ind].ncc_value;
}

void f0tracker_process_frame(PitchState *pst, const vector* in, int* n0, flt_t* corr) {
    flt_t rms, prob_voiced, prob_onset, prob_offset;
    int steplen = pst->steplen, overlap, buflen;
    int lpc_overlap = pst->lpc_buflen - steplen, norm_overlap;
    vector mixture;
    flt_t ncc[pst->n_nccf_lags];
    *n0 = 0;
    *corr = FLT(0.);
    auto_vector_construct(mixture, steplen*2, Float);

    vec_left_moving(pst->in_buf, in, pst->in_buflen - steplen, steplen, 0);

    rms = get_band_rms(pst->fft, in, pst->rms_first_bin, pst->rms_last_bin, pst->ft_window);
    vec_right_moving(pst->rms_queue, NULL, pst->rms_queue->size-1, 1, 0);
    vec_elem_set(pst->rms_queue, 0, &rms);

    prob_voiced = get_rms_voicing(pst, vec_elem_get_FLT(pst->rms_queue, 1));  // use last frame rms

    get_voice_transition(vec_elem_get_FLT(pst->rms_queue, 1),
                         vec_elem_get_FLT(pst->rms_queue, 3),
                         &prob_onset, &prob_offset);  // use last frame rms

    vec_left_moving(pst->lpc_resid, NULL, pst->lpc_buflen - steplen, steplen, 0);
//    lpc_overlap = pst->lpc_buflen - steplen;
//    for (int i=0; i<lpc_overlap; i++)
//        pst->lpc_resid[i] = pst->lpc_resid[i+steplen];  // lpc_resid : 4 * frame_step
    vector_view lpc_vw;
    vec_view_construct(pst->lpc_resid, lpc_overlap, steplen, &lpc_vw);
    get_lpc_residual(pst->lpc, pst->in_buf, pst->in_buflen-steplen, &lpc_vw);  // use new frame

//    norm_overlap = pst->norm_buflen - steplen;
//    for (int i=0; i<norm_overlap; i++)
//        pst->norm_resid[i] = pst->norm_resid[norm_overlap+i];
//    int n_norm_per_frame = steplen / pst->norm_steplen;
//    for (int n=0; n < n_norm_per_frame; n++) {
//        float* from = pst->lpc_resid + lpc_overlap + n*pst->norm_steplen;
//        float* to = pst->norm_resid + norm_overlap + n*pst->norm_steplen;
//        normalize_by_rms(from, to, pst->norm_steplen, &(pst->old_inv_rms));
//    }

    vector_view vw1, vw2;
    vec_view_construct(pst->lpc_resid, steplen*2, steplen*2, &vw1);
    vec_view_construct(pst->in_buf, steplen*2, steplen*2, &vw2);
    vec_mul_scalar(&vw1, &resid_fract_, &mixture);
    vec_mul_scalar_accum(&vw2, &pcm_fract_, &mixture);

    norm_cross_corr(ncc, mixture.data, 0, pst->first_nccf_lag, pst->n_nccf_lags, pst->corr_windowlen);

    int peak_ind[max_peaks_num_-1];  // 留一个给 unvoiced candidate
    flt_t peak_val[max_peaks_num_-1];
    int n_peaks;
    find_ncc_peaks(peak_ind, peak_val,
                   &n_peaks, max_peaks_num_-1,
                   ncc, pst->n_nccf_lags, corr_thresh_);

    construct_lattice(pst->lattice, peak_ind, peak_val,
                      n_peaks, prob_voiced, pst->first_nccf_lag);

    dynamic_programming(pst->lattice, prob_onset, prob_offset);

    backtrack(pst->lattice, n0, corr);
//    printf("%f_%d, ", prob_voiced, *n0);

    if (prob_voiced < output_thresh_) { *n0 = 0; *corr = FLT(0.); }
}

