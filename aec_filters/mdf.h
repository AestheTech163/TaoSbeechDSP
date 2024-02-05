//
// Created by wentao on 2024/1/11.
//

#ifndef TAOSBEECHDSP_MDF_H
#define TAOSBEECHDSP_MDF_H

#include "arch.h"
#include "smallft.h"
#include "vector.h"
#include "common_states.h"

struct MDFFilter_;
typedef struct MDFFilter_ MDFFilter;

typedef void (*init_func)(int n_blocks, int n_frame);
typedef void (*compute_echo_func)(void* mdf, const vector* X);
typedef void (*compute_error_func)(void* mdf, const vector* y);
typedef void (*update_filter_func)(void* mdf, CommonStates* cs);
typedef void (*get_residual_func)(void* aec_obj, vector* resid);

struct MDFFilter_ {
    // multi delay frequency domain adaptive filter
    int n_blocks;
    int filter_length;
    int n_freq;
    int n_fft;
    int n_framelen;
    int n_overlap;

    flt_t Syy;
    flt_t See;
    vector* y_hat;
    vector* e;
    vector* H_hat;
    vector* Y_hat;
    vector* E;
    SmallFFT* fft_handler;

    init_func init;
    compute_echo_func compute_echo;
    compute_error_func compute_error;
    update_filter_func update_filter;
};

//#define compute_echo(filter, X) (filter->mdf->compute(filter->mdf, X))
//#define compute_error(filter, X) (filter->mdf->compute(filter->mdf, X))
//#define update_filter(filter, X) (filter->mdf->compute(filter->mdf, X))
#define GET_fft_handler(filter) ((filter)->mdf->fft_handler)

#define GET_nblocks(filter) ((filter)->mdf->n_blocks)
#define GET_nfreq(filter) ((filter)->mdf->n_freq)
#define GET_nfft(filter) ((filter)->mdf->n_fft)
#define GET_nframelen(filter) ((filter)->mdf->n_framelen)

#define GET_error(filter) ((filter)->mdf->e)
#define GET_Error(filter) ((filter)->mdf->E)
#define GET_See(filter) ((filter)->mdf->See)
#define GET_H_hat(filter) ((filter)->mdf->H_hat)
#define GET_Syy(filter) ((filter)->mdf->Syy)
#define GET_y_hat(filter) ((filter)->mdf->y_hat)
#define GET_Y_hat(filter) ((filter)->mdf->Y_hat)

#define SET_error_func(filter, func) ((filter)->compute_error = (func))
#define SET_echo_func(filter, func) ((filter)->compute_echo = (func))
#define SET_update_func(filter, func) ((filter)->update_filter = (func))

MDFFilter* mdf_init(int n_blocks, int n_framelen);
void mdf_destroy(MDFFilter* mdf);

void mdf_compute_echo(void* filter, const vector* X);
void mdf_compute_error(void* filter, const vector* y);
void mdf_update_filter(void* filter, CommonStates* cs);

#endif //TAOSBEECHDSP_MDF_H
