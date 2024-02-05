//
// Created by wentao on 2023/12/20.
//
#include <math.h>
#include <stdlib.h>

#include "lpc.h"
#include "arch.h"
#include "dsp_common_func.h"


LpcState* init_lpc(int lpc_ord, int stepsize, int windowsize, flt_t preemp, flt_t noise_floor) {
    LpcState* lpc = (LpcState*) malloc(sizeof(LpcState));
    check(lpc, "malloc lpc failed.");
    lpc->lpc_ord = lpc_ord;
    lpc->stepsize = stepsize;
    lpc->windowsize = windowsize;
    lpc->preemp = preemp;
    lpc->pre_mem = FLT(0.);
    lpc->window = vec_alloc(windowsize, Float);  // (T*) malloc(sizeof(float)*windowsize);
    lpc->lpc_buf = vec_alloc(windowsize, Float);  // (float*) malloc(sizeof(float)*windowsize);
    lpc->wind_data = vec_alloc(windowsize, Float);  // (float*) malloc(sizeof(float)*windowsize);
    lpc->old_lpc = vec_alloc(lpc_ord+1, Float);  // (float*) malloc(sizeof(float)*(lpc_ord+1));
    lpc->delta_lpc = vec_alloc(lpc_ord+1, Float);  // (float*)malloc(sizeof(float)*(lpc_ord+1));
    flt_t arg = M_PI * FLT(2.) / windowsize;
    for (int i = 0; i < windowsize; i++)
        vec_elem_set_FLT(lpc->window, i, FLT(0.5) - FLT(0.5) * cos((FLT(0.5) + i) * arg));
    lpc->ffact = FLT(1.) / (FLT(1.) + taodsp_exp((-noise_floor / FLT(20.)) * taodsp_log(FLT(10.))));
    return lpc;
}

void free_lpc(LpcState* lpc) {
    if (!lpc) return;
    BATCH_vec_free(lpc->lpc_buf, lpc->window, lpc->wind_data, lpc->old_lpc, lpc->delta_lpc);
    free(lpc);
}

void autocorr(int windowsize, flt_t* s, int p, flt_t* r) {
    int i, j;
    flt_t e, *q = s, *t, sum, sum0 = 0.0;

    for (i=0; i<windowsize; i++) {
        sum = *q++;
        sum0 += sum*sum;
    }
    *r = 1.;
    if (sum0 == 0.0) {
        for (i = 1; i <= p; i++) r[i] = 0.;
        return;
    }
//    e = sqrt(sum0 / windowsize);
    sum0 = 1.0 / sum0;
    for (i = 1; i <= p; i++) {
        for (sum = 0.0, j = windowsize - i, q = s, t = s + i; j--;)
            sum += (*q++) * (*t++);
        *(++r) = sum * sum0;
    }
}

void durbin(flt_t* r, flt_t* a, int p) {
    flt_t k[p];
    flt_t bb[p];
    int i, j;
    flt_t e, s, *b = bb;

    e = *r;
    *k = -r[1] / e;
    *a = *k;
    e *= (1.0 - (*k) * (*k));

    for (i = 1; i < p; i++) {
        s = 0;
        for (j = 0; j < i; j++) s -= a[j] * r[i - j];
        k[i] = (s - r[i + 1]) / e;
        a[i] = k[i];
        for (j = 0; j <= i; j++) b[j] = a[j];
        for (j = 0; j < i; j++) a[j] += k[i] * b[i - j - 1];
        e *= (1.0 - (k[i] * k[i]));
    }
}

void left_moving_with_preem(vector* buf, const vector* frm, flt_t preemph, flt_t* premem) {
    int i;
    int frmlen = frm->size;
    int ovrlen = buf->size - frm->size;
    flt_t* buf_data = buf->data;
    flt_t* frm_data = frm->data;
    for (i=0; i<ovrlen; i++) buf_data[i] = buf_data[frmlen+i];
    for (i=0; i<frmlen; i++) {
        buf_data[ovrlen+i] = frm_data[i] - preemph * *premem;
        *premem = frm_data[i];
    }
}

void compute_lpc(LpcState* st, const vector* data, vector* lpca) {
    int lpc_ord = st->lpc_ord;
    int stepsize = st->stepsize;
    int overlap = st->windowsize - st->stepsize;
    int i;
    flt_t preemp = st->preemp;
//    T ac[st->lpc_ord + 1];
    vector ac;
    auto_vector_construct(ac, st->lpc_ord+1, Float);

//    for (i = 0; i < overlap; i++)
//        st->lpc_buf[i] = st->lpc_buf[stepsize + i];
//    for (i = 0; i < stepsize; i++) {
//        st->lpc_buf[overlap + i] = data[i] - preemp * st->pre_mem;
//        st->pre_mem = data[i];
//    }
    left_moving_with_preem(st->lpc_buf, data, preemp, &st->pre_mem);

    vec_mul(st->lpc_buf, st->window, st->wind_data);
    autocorr(st->windowsize, st->wind_data->data, st->lpc_ord, ac.data);
    vector_view vw;
    vec_view_construct(&ac, 1, lpc_ord, &vw);
    vec_mul_scalar(&vw, &st->ffact, &vw);

    durbin(ac.data, lpca->data, lpc_ord);

    vec_view_construct(lpca, 0, lpc_ord, &vw);
    vec_reverse(&vw);  // reverse lpc coef
    vec_elem_set(lpca, lpc_ord, &flt_one);
}

void get_lpc_residual(LpcState* st, vector* in_buf, int in_buf_overlap, vector* resid) {
    int cur_i;
    int stepsize = st->stepsize;
    int lpc_ord = st->lpc_ord;
    flt_t new_gain, delta_gain, sum_;
    vector* old_lpc = st->old_lpc;
    vector lpc;  //T lpc[lpc_ord + 1];
    auto_vector_construct(lpc, lpc_ord+1, Float);

    vector_view vw;
//    printf("in_buf %d\n", in_buf->size);
    vec_view_construct(in_buf, in_buf_overlap, stepsize, &vw);
    compute_lpc(st, &vw, &lpc);

    sum(&lpc, &new_gain);
    if (new_gain < 0.) new_gain = 1.;
    delta_gain = (new_gain - st->old_gain) / st->stepsize;

    vec_sub(&lpc, st->old_lpc, st->delta_lpc);
    flt_t norm = 1. / st->stepsize;
    vec_mul_scalar(st->delta_lpc, &norm, st->delta_lpc);

    for (int n = 0; n < st->stepsize; n++) {
        // cur_i = lpc_overlap // 2 - lpc_ord + n
        // cur_i = frame_step//2 - lpc_ord + n
        cur_i = in_buf_overlap - lpc_ord + n;  // TODO
//        T* y = (T*)(in_buf->data) + cur_i;
//            printf("in_buf %d\n", in_buf->size);
        vec_view_construct(in_buf, cur_i, lpc_ord+1, &vw);
        vec_dot(old_lpc, &vw, &sum_);
        vec_add(old_lpc, st->delta_lpc, old_lpc);
        flt_t tmp = sum_/st->old_gain;
        vec_elem_set(resid, n, &tmp);
        st->old_gain += delta_gain;
    }
}

int get_lpc_order(int fs) {
    return (int)(2.5f + ((float)fs / 1000.f));
}
