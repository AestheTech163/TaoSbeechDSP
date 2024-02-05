//
// Created by wentao on 2024/1/11.
//
#ifndef TAOSBEECHDSP_ARCH_H
#define TAOSBEECHDSP_ARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "numeric.h"

#define VecGen_int
#define VecGen_long
#define VecGen_float
#define VecGen_double
#define VecGen_complex64
#define VecGen_complex128

#define VEC_PRINT_LIMIT_REAL 16
#define VEC_PRINT_LIMIT_COMPLEX 8

typedef char bool;

// basic functions alias
#define sqrt_float(b, a) ((b) = sqrtf(a))
#define sqrt_double(b, a) ((b) = sqrt(a))
#define exp_float(b, a) ((b) = expf(a))
#define exp_double(b, a) ((b) = exp(a))
#define log_float(b, a) ((b) = logf(a))
#define log_double(b, a) ((b) = log(a))
#define cos_float(b, a) ((b) = cosf(a))
#define cos_double(b, a) ((b) = cos(a))
#define sin_float(b, a) ((b) = sinf(a))
#define sin_double(b, a) ((b) = sin(a))
#define atan_float(b, a) ((b) = atanf(a))
#define atan_double(b, a) ((b) = atan(a))
#define atan2_float(c, a, b) ((c) = atan2f(a, b))
#define atan2_double(c, a, b) ((c) = atan2(a, b))


#ifdef USE_FLOAT
    typedef float flt_t;  // type of floating number
//    typedef complex64 cpx_t;  // type of complex number
    #define cpx_t complex64
    #define Float Float32
    #define Complex Complex64
    #define FLT(x) (x##f)
    #define taodsp_sqrt sqrtf
    #define taodsp_exp expf
    #define taodsp_abs fabsf
    #define taodsp_log logf
    #define taodsp_atan atanf

    #define flt_zero float_zero
    #define flt_one float_one
    #define cpx_zero complex64_zero
    #define cpx_one complex64_one
    #define vec_elem_set_FLT vec_elem_set_float
    #define vec_elem_incr_FLT vec_elem_incr_float
    #define vec_elem_get_FLT vec_elem_get_float
    #define vec_elem_set_CPX vec_elem_set_complex64
    #define vec_elem_set_CPX_parts vec_elem_set_complex64_parts
    #define arr_dot_FLT arr_dot_float
    #define roundup roundupf
    #define vec_unary_FLT vec_unary_float

#else
    typedef double flt_t;  // type of float number
//    typedef complex128 cpx_t;  // type of complex number
    #define cpx_t complex64
    #define Float Float64
    #define Complex Complex128
    #define FLT(x) (x)
    #define taodsp_sqrt sqrt
    #define taodsp_exp exp
    #define taodsp_abs fabs
    #define taodsp_log log
    #define taodsp_atan atan

    #define flt_zero double_zero
    #define flt_one double_one
    #define cpx_zero complex128_zero
    #define cpx_one complex128_one
    #define vec_elem_set_FLT vec_elem_set_double
    #define vec_elem_incr_FLT vec_elem_incr_double
    #define vec_elem_get_FLT vec_elem_get_double
    #define vec_elem_set_CPX vec_elem_set_complex128
    #define vec_elem_set_CPX_parts vec_elem_set_complex128_parts
    #define arr_dot_FLT arr_dot_double
    #define vec_unary_FLT vec_unary_double
    #define roundup roundupd

#endif


#define check_do(cond, action) {if (cond) {action;}}
#define check_free(mem) {if (mem) {free(mem);}}

#define DEBUF_FILE stdout
#ifdef DEBUG
    #define check(cond, str) if (!(cond)) fatal(str)
    //#define assert(cond, message) {if (!(cond)) {fatal("assertion failed: " #cond "\n" message);}}
    #define fatal(str) _fatal(str, __FILE__, __LINE__)
#else
    #define check(cond, str)
    #define fatal(str)
#endif

static inline void _fatal(const char *str, const char *file, int line) {
    fprintf (DEBUF_FILE, "Fatal (internal) error in %s, line %d: %s\n", file, line, str);
    fflush(DEBUF_FILE);
    abort();
}

#define ONES(size, T) ({ T* obj = (T*)malloc(sizeof(T)*(size)); for(int i=0; i<size; i++) obj[i]=1; obj; })
#define C_ONES(size, T) ({ T* obj = (T*)malloc(sizeof(T)*(size)); for(int i=0; i<size; i++) {obj[i].re=1;obj[i].im=0;} obj; })
#define ZEROS(size, T) ({ T* obj = (T*)calloc(size, sizeof(T)); obj; })

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SQUARE(a) ((a)*(a))
#define CLIP(x, a, b) (((x)>(b))?(b):(((x)<(a))?(a):(x)))


#endif //TAOSBEECHDSP_ARCH_H
