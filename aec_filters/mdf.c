//
// Created by wentao on 2024/1/11.
//

#include "arch.h"
#include "stdlib.h"
#include "smallft.h"
#include "vector.h"
#include "mdf.h"


MDFFilter* mdf_init(int n_blocks, int n_framelen) {
    MDFFilter* mdf = (MDFFilter*)malloc(sizeof(MDFFilter));

    mdf->n_blocks = n_blocks;
    mdf->filter_length = n_blocks * n_framelen;
    mdf->n_framelen = n_framelen;
    mdf->n_fft = n_framelen * 2;
//    mdf->n_overlap = mdf->n_fft - mdf->n_framelen;
    mdf->n_freq = mdf->n_fft / 2 + 1;

    mdf->e = zeros(mdf->n_fft, Float);
    mdf->y_hat = zeros(mdf->n_fft, Float);
    mdf->H_hat = zeros(mdf->n_blocks * mdf->n_freq, Complex);
    mdf->Y_hat = zeros(mdf->n_freq, Complex);
    mdf->E = zeros(mdf->n_freq, Complex);
    mdf->Syy = FLT(0.);
    mdf->See = FLT(0.);

    mdf->fft_handler = small_drft_init(mdf->n_fft);
    mdf->compute_echo = mdf_compute_echo;
    mdf->compute_error = mdf_compute_error;
    mdf->update_filter = mdf_update_filter;

    return mdf;
}

void mdf_destroy(MDFFilter* mdf) {
    if (!mdf) return;
    BATCH_vec_free(mdf->y_hat, mdf->e, mdf->H_hat, mdf->Y_hat, mdf->E);
    small_drft_free(mdf->fft_handler);
    free(mdf);
}

void mdf_compute_echo(void* filter, const vector* X) {
    static int frame_idx = 0;
    frame_idx++;
    MDFFilter* mdf = (MDFFilter*)filter;
    const int flen = mdf->n_framelen;
    const int n_freq = mdf->n_freq;
    vector_view vw1, vw2;
    /*self.Y_hat[:] = np.sum(X[:M, :] *self.H_hat[:, :], axis = 0)
    self.y_hat[:] = np.fft.irfft(self.Y_hat).real
    self.Syy = np.dot(self.y_hat[flen:], self.y_hat[flen:]) */
    vec_fill(mdf->Y_hat, &cpx_zero);
    for (int b=0; b<mdf->n_blocks; b++) {
        vec_view_construct(mdf->H_hat, n_freq*b, n_freq, &vw1);  // H_hat[n_freq*b:n_freq*(b+1)]
        vec_view_construct(X, n_freq*b, n_freq, &vw2);  // X[n_freq*b:n_freq*(b+1)]
        vec_mul_accum(&vw1, &vw2, mdf->Y_hat);
    }
//    vector Y_tmp;
//    auto_vector_construct(Y_tmp, VEC_SIZE(*mdf->H_hat), Complex);
//    vec_mul(mdf->H_hat, X, &Y_tmp);
//    sum2d(&Y_tmp, n_freq, mdf->Y_hat, AlongRow);

    small_drft_backward(mdf->fft_handler, mdf->Y_hat, mdf->y_hat);
    vec_view_construct(mdf->y_hat, flen, flen, &vw1);  // y_hat[flen:]
    vec_dot(&vw1, &vw1, &(mdf->Syy));
}

void mdf_compute_error(void* filter, const vector* y) {
    MDFFilter* mdf = (MDFFilter*)filter;
    const int flen = mdf->n_framelen;
    vector_view vw1, vw2;
    /*self.e[frame_len:] = y - self.y_hat[frame_len:]
    self.e[:frame_len] = 0.
    self.E[:] = np.fft.rfft(self.e)
    self.See = np.dot(self.e[frame_len:], self.e[frame_len:])*/
    vec_view_construct(mdf->y_hat, flen, flen, &vw1);  // y_hat[ovl:]
    vec_view_construct(mdf->e, flen, flen, &vw2);  // e[ovl:]
    vec_sub(y, &vw1, &vw2);

    vec_view_construct(mdf->e, 0, flen, &vw2);
    vec_fill(&vw2, &flt_zero);  // e[:flen] = 0
    small_drft_forward(mdf->fft_handler, mdf->e, mdf->E);

    vec_view_construct(mdf->e, flen, flen, &vw2);  // e[ovl:]
    vec_dot(&vw2, &vw2, &(mdf->See));
}

void mdf_update_filter(void* filter, CommonStates* cs) {
    // MDFFilter* mdf = (MDFFilter*)filter;
    // TODO
}
