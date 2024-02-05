//
// Created by wentao on 2024/1/13.
//

#ifndef TAOSBEECHDSP_AEC_PNLMS_H
#define TAOSBEECHDSP_AEC_PNLMS_H

#include "mdf.h"

typedef struct pNLMSFilter_ {
    MDFFilter* mdf;  // parent class instance
    CommonStates* cs;

    compute_echo_func compute_echo;
    compute_error_func compute_error;
    update_filter_func update_filter;

    int n_freq;
    int n_blocks;
    int n_fft;

    vector* X2_smth;
    vector* inv_X2_smth;
    vector* prop;

    flt_t ss;
} pNLMSFilter;

pNLMSFilter* pNLMSFilter_init(int n_blocks, int n_frame, CommonStates* common_states);
void pNLMSFilter_destroy(pNLMSFilter* pnlms);
void pNLMSFilter_update_filter(void* filter, CommonStates* cs);
void pNLMSFilter_process_frame(void* filter,
                               vector* out,
                               const vector* y_frame,
                               const vector* x_frame);

#endif //TAOSBEECHDSP_AEC_PNLMS_H
