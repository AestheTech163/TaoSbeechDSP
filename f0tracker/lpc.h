//
// Created by wentao on 2023/12/20.
//

#ifndef MMSENOISE_LPC_H
#define MMSENOISE_LPC_H

#include "arch.h"
#include "vector.h"

typedef struct _LpcState {
    int lpc_ord;
    int windowsize;
    int stepsize;
    vector* window;
    vector* wind_data;
    vector* lpc_buf;
    flt_t preemp;
    flt_t pre_mem;
    flt_t ffact;
    vector* old_lpc;
    vector* delta_lpc;
    flt_t old_gain;
} LpcState;

LpcState* init_lpc(int lpc_ord, int stepsize, int windowsize, flt_t preemp, flt_t noise_floor);

void free_lpc(LpcState* lpc);

void compute_lpc(LpcState* st, const vector* data, vector* lpca);

void get_lpc_residual(LpcState* st, vector* in_buf, int in_buf_overlap, vector* resid);

int get_lpc_order(int fs);

#endif //MMSENOISE_LPC_H
