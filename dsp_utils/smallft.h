/********************************************************************
 *                                                                  *
 * THIS FILE IS PART OF THE OggVorbis SOFTWARE CODEC SOURCE CODE.   *
 * USE, DISTRIBUTION AND REPRODUCTION OF THIS LIBRARY SOURCE IS     *
 * GOVERNED BY A BSD-STYLE SOURCE LICENSE INCLUDED WITH THIS SOURCE *
 * IN 'COPYING'. PLEASE READ THESE TERMS BEFORE DISTRIBUTING.       *
 *                                                                  *
 * THE OggVorbis SOURCE CODE IS (C) COPYRIGHT 1994-2001             *
 * by the XIPHOPHORUS Company http://www.xiph.org/                  *
 *                                                                  */
/**
   @file smallft.h
   @brief Discrete Rotational Fourier Transform (DRFT)
*/

#ifndef _V_SMFT_H_
#define _V_SMFT_H_

#include "vector.h"

/** Discrete Rotational Fourier Transform lookup */
typedef struct _SmallFFT {
    int n;
    flt_t* trigcache;
    int* splitcache;
    flt_t* tmpbuf;
} SmallFFT;

extern SmallFFT* small_drft_init(int n);
extern void small_drft_free(SmallFFT* l);

void small_drft_unpack(const flt_t* in, vector* out);  // to n_fft/2+1 complex vector
void small_drft_pack(const vector* in, flt_t* out);    // to n_fft drft packed vector

extern void small_drft_forward(SmallFFT* l, const vector* in, vector* out);
extern void small_drft_backward(SmallFFT* l, const vector* in, vector* out);

extern void small_drft_forward_batch(SmallFFT* l, const matrix* in, matrix* out, int axis);
extern void small_drft_backward_batch(SmallFFT* l, const matrix* in, matrix* out, int axis);

//extern float get_band_rms(SmallFFT* fft, const float *x, int first_bin, int last_bin, float* window);

#endif
