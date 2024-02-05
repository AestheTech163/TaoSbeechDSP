//
// Created by wentao on 2024/2/1.
//

#ifndef TAOSBEECHDSP_DSP_COMMON_FUNC_H
#define TAOSBEECHDSP_DSP_COMMON_FUNC_H

#include "arch.h"
#include "complex.h"
#include "smallft.h"
#include "vector.h"
#include "reduce.h"

void left_moving(flt_t* buf, const flt_t* frm, int ovrlen, int frmlen, int reverse_frm);
void right_moving(flt_t* buf, const flt_t* frm, int ovrlen, int frmlen, int reverse_frm);
void vec_right_moving(vector* buf, const vector* frm, int ovrlen, int frmlen, int reverse_frm);
void vec_left_moving(vector* buf, const vector* frm, int ovrlen, int frmlen, int reverse_frm);

void biquad(flt_t* y, flt_t mem[2], const flt_t* x, const flt_t* b, const flt_t* a, int N);
void dc_notch(flt_t* mem, const flt_t* in, flt_t* out, int N, flt_t radius);
void dc_notch_vec(flt_t mem[2], const vector* in, vector* out, flt_t radius);

void linearize_convolution(SmallFFT* fft_handler, vector* H);
void linearize_convolution_2d(SmallFFT* fft_handler, matrix* H);

int signal_saturated(const vector* y, flt_t inf, flt_t sup);
flt_t get_band_rms(SmallFFT* fft, const vector* x, int first_bin, int last_bin, const vector* window);

#endif //TAOSBEECHDSP_DSP_COMMON_FUNC_H
