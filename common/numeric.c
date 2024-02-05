//
// Created by wentao on 2024/2/2.
//
#include "numeric.h"

#define GEN__complex_sqrt(ComplexType, RealType)\
void ComplexType##_sqrt(ComplexType* b, const ComplexType* a) {\
    RealType abs_a, abs_b, phase_a, phase_b;\
    abs_a = ComplexType##_abs(*a);\
    atan2_##RealType(phase_a, (*a).im, (*a).re);\
    sqrt_##RealType(abs_b, abs_a);\
    phase_b = phase_a / (RealType)2;\
    cos_##RealType(phase_a, phase_b);\
    (*b).re = abs_b * phase_a;\
    sin_##RealType(phase_a, phase_b);\
    (*b).im = abs_b * phase_a;\
}
GEN__complex_sqrt(complex64, float)
GEN__complex_sqrt(complex128, double)

