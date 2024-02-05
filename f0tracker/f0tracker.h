//
// Created by wentao on 2023/12/20.
//
#ifndef TAOSBEECHDSP_F0TRACKER_H
#define TAOSBEECHDSP_F0TRACKER_H

#include "lpc.h"
#include "dsp_utils/dsp_common_func.h"
#include "dsp_utils/smallft.h"


#define lpc_preemp_ FLT(0.98)
#define lpc_noise_floor_ FLT(70.)

#define peak_delay_  FLT(0.0004)  // for measuring prominence
#define skew_delay_  FLT(0.00015) // for measuring shape
#define peak_val_wt_  FLT(0.1)
#define peak_prominence_wt_  FLT(0.3)
#define peak_skew_wt_  FLT(0.1)
#define peak_quality_floor_  FLT(0.01)
#define min_peak  FLT(-2.0)
#define period_deviation_wt_  FLT(1.0)
#define reward_  FLT(-2.5)
static flt_t resid_fract_ = FLT(0.7);
static flt_t pcm_fract_ = FLT(0.3);
#define max_peaks_num_  20

#define corr_window_dur_  FLT(0.0075)
#define corr_thresh_  FLT(0.2)

#define max_f0_search_  FLT(500.)
#define min_f0_search_  FLT(80.)

#define time_span_  FLT(0.02)
#define internal_frame_interval_  FLT(0.002)
#define level_change_den_ = FLT(30.)

#define output_thresh_ FLT(0.3)

typedef struct PeriodCandidate_ {
    char voiced;
    flt_t ncc_value;
    int ncc_period;
    flt_t local_cost;
    flt_t accum_cost;
    int best_prev;
} PeriodCandidate;

typedef struct ShabbySwapQueue_ {
//    PeriodCandidate* last;
//    PeriodCandidate* second_last;
    PeriodCandidate this_[max_peaks_num_];
    PeriodCandidate that_[max_peaks_num_];
    PeriodCandidate* this;
    PeriodCandidate* that;
    int this_len;
    int that_len;
    int list_len;
//    int that_frame_index;
//    int this_frame_index;
} ShabbySwapQueue;

typedef struct _PitchState {
    int steplen;
    vector* in_buf;
    int in_buflen;

    vector* ft_window;
    vector* rms_queue;
    int rms_first_bin;
    int rms_last_bin;
    flt_t rms_min_val;
    flt_t rms_max_val;

    int lpc_buflen;
    vector* lpc_resid;

    int norm_buflen;
    vector* norm_resid;
    int norm_steplen;  // TODO 目前先不增加复杂度
    flt_t old_inv_rms;
    int peak_ind;  // TODO 目前先不增加复杂度
    int skew_ind;

    int corr_windowlen;
    int first_nccf_lag;
    int max_lag;
    int n_nccf_lags;

    ShabbySwapQueue* lattice;
    LpcState* lpc;
    SmallFFT* fft;
} PitchState;

PitchState* init_pitch(int fs, int steplen, int windowlen, flt_t rms_min_freq, flt_t rms_max_freq);
void free_pitch(PitchState* pst);
void f0tracker_process_frame(PitchState *pst, const vector* in, int* n0, flt_t* corr);

#endif //TAOSBEECHDSP_F0TRACKER_H
