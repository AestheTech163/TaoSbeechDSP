//
// Created by wentao on 2024/1/11.
//

#ifndef TAOSBEECHDSP_NUMERIC_H
#define TAOSBEECHDSP_NUMERIC_H

#include "arch.h"
#include "float.h"

static inline int roundupf(float val) { return (int)(val + 0.5f); }
static inline int roundupd(double val) { return (int)(val + 0.5); }

#ifndef M_PI
#define M_PI (3.14159265359)
#endif

typedef enum Dtype_ {
    Complex64 = 0,  // 复数
    Complex128,
    Float32,     // 实数
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
} Dtype;

#define IS_VALID_DTYPE(x) ((x)>=Complex64&&(x)<=Int64)
#define IS_REAL(x) ((x)>=Float32&&(x)<=Float64)
#define IS_COMPLEX(x) ((x)>=Complex64&&(x)<=Complex128)
#define IS_COMPARABLE(x) ((x)>Complex128&&(x)<=Int64)

#define AlmostEqualF(a, b) (fabsf(a - b) <= FLT_EPSILON)
#define AlmostEqualD(a, b) (fabs(a - b) <= DBL_EPSILON)

typedef struct complex64_ {
    float re, im;
} complex64;

typedef struct complex128_ {
    double re, im;
} complex128;

#define complex64_abs2(x) (SQUARE((x).re) + SQUARE((x).im))
#define complex128_abs2 complex64_abs2
#define complex64_abs(x) sqrtf(complex64_abs2(x))
#define complex128_abs(x) sqrt(complex128_abs2(x))
void complex64_sqrt(complex64* b, const complex64* a);
void complex128_sqrt(complex128* b, const complex128* a);

#define complex64_assign(b, a) ((b) = (a))
#define complex64_add(c, a, b) {(c).re=(a).re+(b).re; (c).im=(a).im+(b).im;}
#define complex64_add_assign(b, a) {(b).re += (a).re, (b).im += (a).im;}
#define complex64_sub(c, a, b) {(c).re=(a).re-(b).re; (c).im=(a).im-(b).im;}
#define complex64_sub_assign(b, a) {(b).re -= (a).re, (b).im -= (a).im;}
#define complex64_mul(c, a, b) {complex64 _a=(a),_b=(b);(c).re=_a.re*_b.re-_a.im*_b.im;(c).im=_a.re*_b.im+_a.im*_b.re;}
#define complex64_mul_assign(b, a) {complex64 _b=(b);(b).re=(a).re*_b.re-(a).im*_b.im;(b).im=(a).re*_b.im+(a).im*_b.re;}
#define complex64_div(c, a, b) {complex64 _a = a, _b = b;float den_recip = 1.f/complex64_abs2(_b); \
    (c).re = den_recip*(_a.re*_b.re+_a.im*_b.im);(c).im = den_recip*(_a.im*_b.re-_a.re*_b.im);}
#define complex64_div_assign(b, a) { \
    complex64 _b = b; float den_recip = 1.f/complex64_abs2(a); \
    (b).re = den_recip*(_a.re*_b.re+_a.im*_b.im);(b).im = den_recip*(_a.im*_b.re-_a.re*_b.im);}
#define complex64_conj(b, a) {(b).re = (a).re; (b).im = -(a).im;}

#define complex128_assign complex64_assign
#define complex128_add complex64_add
#define complex128_add_assign complex64_add_assign
#define complex128_sub complex64_sub
#define complex128_sub_assign complex64_sub_assign
#define complex128_mul(c, a, b) {complex128 _a=(a),_b=(b);(c).re=_a.re*_b.re-_a.im*_b.im;(c).im=_a.re*_b.im+_a.im*_b.re;}
#define complex128_mul_assign(b, a) {complex128 _b=(b);(b).re=(a).re*_b.re-(a).im*_b.im;(b).im=(a).re*_b.im+(a).im*_b.re;}
#define complex128_div(c, a, b) {complex128 _a = a, _b = b;float den_recip = 1.f/complex128_abs2(_b); \
    (c).re = den_recip*(_a.re*_b.re+_a.im*_b.im);(c).im = den_recip*(_a.im*_b.re-_a.re*_b.im);}
#define complex128_div_assign(b, a) { \
    complex128 _b = b; float den_recip = 1.f/complex128_abs2(a); \
    (b).re = den_recip*(_a.re*_b.re+_a.im*_b.im);(b).im = den_recip*(_a.im*_b.re-_a.re*_b.im);}
#define complex128_conj complex64_conj

// mixed type operations
#define complex64_add_float(c, a, b) {(c).re=(a).re+(b); (c).im=(a).im;}
#define float_add_complex64(c, a, b) {(c).re=(a)+(b).re; (c).im=(b).im;}
#define complex64_mul_float(c, a, b) {(c).re=(a).re*(b); (c).im=(a).im*(b);}
#define float_mul_complex64(c, a, b) {(c).re=(a)*(b).re; (c).im=(a)*(b).im;}
#define complex64_sub_float(c, a, b) {(c).re=(a).re-(b); (c).im=(a).im;}
#define float_sub_complex64(c, a, b) {(c).re=(a)-(b).re; (c).im=-(b).im;}
#define complex64_div_float(c, a, b) {(c).re=(a).re/(b); (c).im=(a).im/(b);}
#define float_div_complex64(c, a, b) {float den_recip=1.f/complex64_abs2(b);(c).re=(a)*(b).re*den_recip;(c).im=-(a)*(b).im*den_recip;}

#define complex128_add_double complex64_add_float
#define double_add_complex128 float_add_complex64
#define complex128_sub_double complex64_sub_float
#define double_sub_complex128 float_sub_complex64
#define complex128_mul_double complex64_mul_float
#define double_mul_complex128 float_mul_complex64
#define complex128_div_double complex64_div_float
#define double_div_complex128(c, a, b) {double den_recip=1./complex128_abs2(b); (c).re=(a)*(b).re*den_recip;(c).im=-(a)*(b).im*den_recip;}

// Complex exponential, using Euler's formula
#define exp_complex64(b, a) {float tmp=expf((a).re); (b).re=tmp*cosf((a).im); (b).im=tmp*sinf((a).im);}
#define exp_complex128(b, a) {double tmp=exp((a).re); (b).re=tmp*cos((a).im); (b).im=tmp*sin((a).im);}


// convenient macros
#define int_assign(b, a) ((b) = (a))
#define int_add(c, a, b) ((c) = (a) + (b))
#define int_add_assign(b, a) ((b) += (a))
#define int_sub(c, a, b) ((c) = (a) - (b))
#define int_mul(c, a, b) ((c) = (a) * (b))
#define int_div(c, a, b) ((c) = (a) / (b))

#define char_assign int_assign
#define char_add int_add
#define char_add_assign int_add_assign
#define char_sub int_sub
#define char_mul int_mul
#define char_div int_div

#define short_assign int_assign
#define short_add int_add
#define short_add_assign int_add_assign
#define short_sub int_sub
#define short_mul int_mul
#define short_div int_div

#define long_assign(b, a) ((b) = (a))
#define long_add_assign(b, a) ((b) += (a))
#define long_add(c, a, b) ((c) = (a) + (b))
#define long_sub(c, a, b) ((c) = (a) - (b))
#define long_mul(c, a, b) ((c) = (a) * (b))
#define long_div(c, a, b) ((c) = (a) / (b))

#define float_assign(b, a) ((b) = (a))
#define float_add_assign(b, a) ((b) += (a))
#define float_add(c, a, b) ((c) = (a) + (b))
#define float_sub(c, a, b) ((c) = (a) - (b))
#define float_mul(c, a, b) ((c) = (a) * (b))
#define float_div(c, a, b) ((c) = (a) / (b))

#define double_assign(b, a) ((b) = (a))
#define double_add_assign(b, a) ((b) += (a))
#define double_add(c, a, b) ((c) = (a) + (b))
#define double_sub(c, a, b) ((c) = (a) - (b))
#define double_mul(c, a, b) ((c) = (a) * (b))
#define double_div(c, a, b) ((c) = (a) / (b))

#define char_greater(a, b) ((a)>(b))
#define short_greater(a, b) ((a)>(b))
#define int_greater(a, b) ((a)>(b))
#define long_greater(a, b) ((a)>(b))
#define float_greater(a, b) ((a)>(b))
#define double_greater(a, b) ((a)>(b))

#define char_less(a, b) ((a)<(b))
#define short_less(a, b) ((a)<(b))
#define int_less(a, b) ((a)<(b))
#define long_less(a, b) ((a)<(b))
#define float_less(a, b) ((a)<(b))
#define double_less(a, b) ((a)<(b))

#define int_maximum(c, a, b) ((c)=(a)>(b)?(a):(b))
#define long_maximum(c, a, b) ((c)=(a)>(b)?(a):(b))
#define float_maximum(c, a, b) ((c)=(a)>(b)?(a):(b))
#define double_maximum(c, a, b) ((c)=(a)>(b)?(a):(b))
#define char_maximum int_maximum
#define short_maximum int_maximum

#define int_minimum(c, a, b) ((c)=(a)<(b)?(a):(b))
#define long_minimum(c, a, b) ((c)=(a)<(b)?(a):(b))
#define float_minimum(c, a, b) ((c)=(a)<(b)?(a):(b))
#define double_minimum(c, a, b) ((c)=(a)<(b)?(a):(b))
#define char_minimum int_minimum
#define short_minimum int_minimum


// define neutral number of add and mul
#define char_ZERO 0
#define short_ZERO 0
#define int_ZERO 0
#define long_ZERO 0
#define float_ZERO 0.f
#define double_ZERO 0.
#define complex64_ZERO {0.f, 0.f}
#define complex128_ZERO {0., 0.}
#define char_ONE 1
#define short_ONE 1
#define int_ONE 1
#define long_ONE 1
#define float_ONE 1.f
#define double_ONE 1.
#define complex64_ONE {1.f, 0.f}
#define complex128_ONE {1., 0.}

static const char char_zero = 0;
static const short short_zero = 0;
static const int int_zero = 0;
static const long long_zero = 0;
static const float float_zero = 0.f;
static const double double_zero = 0.;
static const complex64 complex64_zero = {0.f, 0.f};
static const complex128 complex128_zero = {0., 0.};
static const char char_one = 1;
static const short short_one = 1;
static const int int_one = 1;
static const long long_one = 1;
static const float float_one = 1.f;
static const double double_one = 1.;
static const complex64 complex64_one = {1.f, 0.f};
static const complex128 complex128_one = {1., 0.};

#endif //TAOSBEECHDSP_NUMERIC_H
