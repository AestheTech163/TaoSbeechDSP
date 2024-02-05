//
// Created by wentao on 2024/1/11.
//
#include "dsp_common_func.h"
#include "vector.h"
#include "smallft.h"
#include "reduce.h"

#define GEN__right_moving(T)\
void right_moving_##T(T* buf, const T* frm, int ovrlen, int frmlen, int reverse_frm) {\
    int i;\
    for (i=ovrlen+frmlen-1; i>=frmlen; i--) buf[i] = buf[i-frmlen];\
    if (frm) {\
        if (reverse_frm) for (i = frmlen - 1; i >= 0; i--) buf[i] = frm[frmlen - 1 - i];\
        else for (i = 0; i < frmlen; i++) buf[i] = frm[i];\
    }\
}
GEN__right_moving(char)
GEN__right_moving(short)
GEN__right_moving(int)
GEN__right_moving(long)
GEN__right_moving(float)
GEN__right_moving(double)
GEN__right_moving(complex64)
GEN__right_moving(complex128)
void vec_right_moving(vector* buf, const vector* frm, int ovrlen, int frmlen, int reverse_frm) {
    void* frm_data = NULL;
    if (frm) frm_data = frm->data;
    SIMPLE_SWITCH_NO_RETURN(right_moving, buf->dtype, buf->data, frm_data, ovrlen, frmlen, reverse_frm)
}

#define GEN__left_moving(T)\
void left_moving_##T(T* buf, const T* frm, int ovrlen, int frmlen, int reverse_frm) {\
    int i;\
    for (i=0; i<ovrlen; i++) buf[i] = buf[frmlen+i];\
    if (frm) {\
        if (reverse_frm) for (i=0; i<frmlen; i++) buf[ovrlen+i] = frm[frmlen-1-i];\
        else for (i=0; i<frmlen; i++) buf[ovrlen+i] = frm[i];\
    }\
}
GEN__left_moving(char)
GEN__left_moving(short)
GEN__left_moving(int)
GEN__left_moving(long)
GEN__left_moving(float)
GEN__left_moving(double)
GEN__left_moving(complex64)
GEN__left_moving(complex128)
void vec_left_moving(vector* buf, const vector* frm, int ovrlen, int frmlen, int reverse_frm) {
    void* frm_data = NULL;
    if (frm) frm_data = frm->data;
    SIMPLE_SWITCH_NO_RETURN(left_moving, buf->dtype, buf->data, frm_data, ovrlen, frmlen, reverse_frm)
}

void biquad(flt_t *y, flt_t mem[2], const flt_t *x, const flt_t *b, const flt_t *a, int N) {
    int i;
    for (i=0; i<N; i++) {
        flt_t xi, yi;
        xi = x[i];
        yi = x[i] + mem[0];
        mem[0] = mem[1] + (b[0]*(flt_t)xi - a[0]*(flt_t)yi);
        mem[1] = (b[1]*(flt_t)xi - a[1]*(flt_t)yi);
        y[i] = yi;
    }
}

void dc_notch(flt_t *mem, const flt_t *in, flt_t *out, int N, flt_t radius) {
    flt_t den2 = radius * radius + 0.7 * (1 - radius) * (1 - radius);
    flt_t vin, vout;
    for (int i=0; i<N; i++) {
        vin = in[i];
        vout = mem[0] + vin;
        mem[0] = mem[1] + 2 * (radius * vout - vin);
        mem[1] = vin - den2 * vout;
        out[i] = radius * vout;
    }
}

void dc_notch_vec(flt_t mem[2], const vector* in, vector* out, flt_t radius) {
    flt_t den2 = radius * radius + FLT(0.7) * (1 - radius) * (1 - radius);
    flt_t vin, vout;
    flt_t* in_data = in->data;
    flt_t* out_data = out->data;
    for (int i=0; i<in->size; i++) {
        vin = in_data[i];
        vout = mem[0] + vin;
        mem[0] = mem[1] + 2 * (radius * vout - vin);
        mem[1] = vin - den2 * vout;
        out_data[i] = radius * vout;
    }
}

void linearize_convolution(SmallFFT* fft_handler, vector* H) {
    // H shape: 1 x n_freq
    int F = fft_handler->n / 2;
    vector wtmp;
    auto_vector_construct(wtmp, fft_handler->n, Float);

    vector_view wtmp_half;
    vec_view_construct(&wtmp, F, F, &wtmp_half);

    small_drft_backward(fft_handler, H, &wtmp);
    vec_fill(&wtmp_half, &flt_zero);
    small_drft_forward(fft_handler, &wtmp, H);
}

void linearize_convolution_2d(SmallFFT* fft_handler, matrix* H) {
    /* cyclic convolution to linear convolution*/
    // H shape: n_blocks x n_freq
    /* wtmp = np.fft.irfft(self.H_hat, axis=1).real
        wtmp[:, F:] = 0.
        H_hat[:, :] = np.fft.rfft(wtmp, axis=1) */
    vector_view wtmp_half, H_vw;
    int F = fft_handler->n / 2;
    vector wtmp;
    auto_vector_construct(wtmp, fft_handler->n, Float);
    for (int b=0; b<H->row; b++) {
        mat_get_row(H, b, &H_vw);
        small_drft_backward(fft_handler, &H_vw, &wtmp);
        vec_view_construct(&wtmp, F, F, &wtmp_half);
        vec_fill(&wtmp_half, &flt_zero);
        small_drft_forward(fft_handler, &wtmp, &H_vw);
    }
}

int signal_saturated(const vector* y, flt_t inf, flt_t sup) {
    flt_t* d = y->data;
    for (int i=0; i<y->size; i++) {
        if ((d[i] > sup) || (d[i] < inf))
            return 1;
    }
    return 0;
}

flt_t get_band_rms(SmallFFT* fft, const vector* x, int first_bin, int last_bin, const vector* window) {
    flt_t rms;
    vector wind_x, X, X2;
    auto_vector_construct(wind_x, fft->n, Float);
    auto_vector_construct(X, fft->n/2+1, Complex);
    auto_vector_construct(X2, last_bin-first_bin, Float);

    if (window != NULL) vec_mul(x, window, &wind_x);
    else vec_copy(x, &wind_x);

    small_drft_forward(fft, &wind_x, &X);
    vector_view X_vw;
    vec_view_construct(&X, first_bin, last_bin-first_bin, &X_vw);
    vec_abs2(&X_vw, &X2);
    sum(&X2, &rms);
    rms = FLT(20.) * log10(FLT(1.) + taodsp_sqrt(rms / (last_bin - first_bin + 1)));
    return rms;
}