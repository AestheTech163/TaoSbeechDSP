//
// Created by wentao on 2023/12/23.
//
// from speexdsp

#ifndef TAOSBEECHDSP_MMSE_NOISE_H
#define TAOSBEECHDSP_MMSE_NOISE_H

#include "filterbank.h"
#include "smallft.h"

#define speech_prob_start 0.35
#define speech_prob_continue 0.20
#define noise_suppress (-35)
#define echo_suppress (-40)
#define echo_suppress_active (-25)


typedef struct MMSENoise_ {
    int frame_size;
    int sampling_rate;
    int n_bands;
    int n_adapt;
    int n_freq;
    int min_count;

    vector* S;
    vector* Smin;
    vector* Stmp;
    vector* ft;
    vector* frame;
    vector* ps;
    vector* xbands;
    vector* old_xbands;
    vector* bands_gain;
    vector* bands_gain2;
    vector* linear_gain;
    vector* gain_floor;
    vector* noise;
    vector* noise_bands;
//    vector* speech_prob;
    vector* echo_noise;
    vector* echo_bands;
    vector* resid_echo;
    vector* post;
    vector* prior;
    vector* zeta;

    vector* inbuf;
    vector* outbuf;
    vector* window;
    flt_t noise_floor;

    SmallFFT* fft;
    FilterBank* fbank;
} MMSENoise;

MMSENoise* init_MMSENoise(int frame_size, int sampling_rate, int n_bands);

void free_MMSENoise(MMSENoise* mmse_ns);

void mmsenoise_process_frame(MMSENoise* nst, vector* x, void* aec_st);

#endif //TAOSBEECHDSP_MMSE_NOISE_H
