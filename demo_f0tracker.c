//
// Created by wentao on 2024/1/26.
//
#include "math.h"
#include "stdio.h"
#include "vector.h"
#include "dsp_common_func.h"
#include "f0tracker.h"
#include "arch.h"


int main(int ac, char** av) {
    FILE* in_fd, *f0_fd;
    const char* in_file = "../test_files/p225_242_mic1_16k.pcm";
    const char* tmp_f0_filepath = "../test_files/results/p225_242_mic1_16k.f0bin";

    int fs = 16000;
    int stepsize = 160;
    int i, n0;
    vector notch_buf;
    auto_vector_construct(notch_buf, stepsize, Float);
    flt_t corr, pitch;
    flt_t notch_mem[2] = {0., 0.};
    short tmp_buf[stepsize];

    PitchState* pst = init_pitch(fs, stepsize, stepsize*2,100, 1000);

    in_fd = fopen(in_file, "rb");
    check(in_fd, "open input file error.");
    f0_fd = fopen(tmp_f0_filepath, "wb");
    check(f0_fd, "open fw_f0 file error.");

    while (1) {
        fread(tmp_buf, sizeof(short), stepsize, in_fd);
        if (feof(in_fd)) break;
        for (i = 0; i < stepsize; i++) vec_elem_set_FLT(&notch_buf, i, tmp_buf[i]);
        dc_notch_vec(notch_mem, &notch_buf, &notch_buf, 0.972);

        f0tracker_process_frame(pst, &notch_buf, &n0, &corr);

        if ((n0 >= 32) && (n0 <=200)) pitch = (flt_t)fs / (flt_t)n0;
        else pitch = FLT(0.);
        i = fwrite(&pitch, sizeof(flt_t), 1, f0_fd);
        check(i == 1, "Problems writing pitch into file stream.");
    }

    fclose(in_fd);
    fclose(f0_fd);
    free_pitch(pst);
    printf("Done.\n");

    return 0;
}
