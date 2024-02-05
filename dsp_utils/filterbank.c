//
// Created by wentao on 2023/12/23.
//

#include "stdlib.h"
#include "memory.h"
#include "math.h"

#include "arch.h"
#include "filterbank.h"

static inline flt_t to_bank(flt_t hz) {
    return FLT(13.1) * taodsp_atan(FLT(0.00074) * hz) + FLT(2.24) * taodsp_atan(hz * hz * FLT(1.85e-8)) + FLT(1e-4) * hz;
}

FilterBank* init_FilterBank(int n_bank, int sampling_rate, int n_freq) {
    flt_t df = (flt_t)sampling_rate / (FLT(2.) * (flt_t)n_freq);
    flt_t max_mel = to_bank((flt_t)sampling_rate / 2);
    flt_t mel_interval = max_mel / (n_bank - 1);
    flt_t val;
    int id1, id2, i;

    FilterBank* fbank = (FilterBank*)malloc(sizeof(FilterBank));
    check(fbank, "malloc fbank failed.");

    fbank->n_freq = n_freq;
    fbank->n_bank = n_bank;
    fbank->scaling = zeros(n_bank, Float);
    fbank->bank_left = zeros(n_freq, Int32); // (int*)malloc(sizeof(int)*n_freq);
    fbank->bank_right = zeros(n_freq, Int32); // (int*)malloc(sizeof(int)*n_freq);
    fbank->filter_left = zeros(n_freq, Float); // (float*)malloc(sizeof(float)*n_freq);
    fbank->filter_right = zeros(n_freq, Float); // (float*)malloc(sizeof(float)*n_freq);

    for (i=0; i<n_freq; i++) {
        flt_t curr_freq = i * df;
        flt_t mel = to_bank(curr_freq);
        if (mel > max_mel) break;
        id1 = (int)(mel / mel_interval);
        if (id1 > n_bank - 2){
            id1 = n_bank - 2;
            val = FLT(1.0);
        } else
            val = (mel - id1 * mel_interval) / mel_interval;
        id2 = id1 + 1;
        vec_elem_set_int(fbank->bank_left, i, id1);
        vec_elem_set_FLT(fbank->filter_left, i, FLT(1.)-val);
        vec_elem_set_int(fbank->bank_right, i, id2);
        vec_elem_set_FLT(fbank->filter_right, i, val);
    }
    for (i=0; i<fbank->n_bank; i++) vec_elem_set_FLT(fbank->scaling, i, FLT(0.));
    for (i=0; i<fbank->n_freq; i++) {
        int id = vec_elem_get_int(fbank->bank_left, i);
        vec_elem_incr_FLT(fbank->scaling, id, vec_elem_get_FLT(fbank->filter_left, i));
        id = vec_elem_get_int(fbank->bank_right, i);
        vec_elem_incr_FLT(fbank->scaling, id, vec_elem_get_FLT(fbank->filter_right, i));
    }
    for (i=0; i<fbank->n_bank; i++) {
        flt_t tmp = vec_elem_get_FLT(fbank->scaling, i);
        vec_elem_set_FLT(fbank->scaling, i, FLT(1.) / (tmp));
    }
    return fbank;
}

void free_FilterBank(FilterBank* fbank) {
    if (!fbank) return;
    BATCH_vec_free(fbank->scaling, fbank->bank_left,
                   fbank->bank_right, fbank->filter_left,
                   fbank->filter_right);
    free(fbank);
}

void psd_to_bank(FilterBank* fbank, const vector* ps, vector* bank) {
    int* bank_left = fbank->bank_left->data;
    int* bank_right = fbank->bank_right->data;
    flt_t* filter_left = fbank->filter_left->data;
    flt_t* filter_right = fbank->filter_right->data;
    flt_t* ps_data = ps->data;
    flt_t* bank_data = bank->data;

    flt_t zero = 0;
    vec_fill(bank, &zero);
    for (int i=0; i<fbank->n_freq; i++) {
        int idx = bank_left[i];
        bank_data[idx] += filter_left[i] * ps_data[i];
        idx = bank_right[i];
        bank_data[idx] += filter_right[i] * ps_data[i];
    }
//    for (int i=0; i<fbank->n_bank; i++)
//        bank[i] *= fbank->scaling[i];
}

void bank_to_psd(FilterBank* fbank, const vector* bank, vector* ps) {
    int* bank_left = fbank->bank_left->data;
    int* bank_right = fbank->bank_right->data;
    flt_t* filter_left = fbank->filter_left->data;
    flt_t* filter_right = fbank->filter_right->data;
    flt_t* ps_data = ps->data;
    flt_t* bank_data = bank->data;

    for (int i=0; i<fbank->n_freq; i++) {
        int idx = bank_left[i];
        ps_data[i] = bank_data[idx] * filter_left[i];
        idx = bank_right[i];
        ps_data[i] += bank_data[idx] * filter_right[i];
    }
}
