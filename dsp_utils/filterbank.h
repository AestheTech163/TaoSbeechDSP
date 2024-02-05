//
// Created by wentao on 2023/12/23.
//
#ifndef TAOSBEECHDSP_FILTERBANK_H
#define TAOSBEECHDSP_FILTERBANK_H

#include "vector.h"

typedef struct FilterBank_ {
    int n_bank;
    int n_freq;
    vector* bank_left;
    vector* bank_right;
    vector* filter_left;
    vector* filter_right;
    vector* scaling;
} FilterBank;

FilterBank* init_FilterBank(int n_bank, int sampling_rate, int n_freq);
void free_FilterBank(FilterBank* fbank);
void psd_to_bank(FilterBank* fbank, const vector* ps, vector* bank);
void bank_to_psd(FilterBank* fbank, const vector* bank, vector* ps);

#endif //TAOSBEECHDSP_FILTERBANK_H
