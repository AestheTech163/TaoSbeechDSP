//
// Created by wentao on 2024/1/30.
//

#ifndef TAOSBEECHDSP_REDUCE_H
#define TAOSBEECHDSP_REDUCE_H

#include "arch.h"
#include "vector.h"
#include "numeric.h"

typedef enum ReduceAxis_ {
    AlongRow = 0,
    AlongCol = 1,
    FullDim = 2,
} ReduceAxis;

typedef void (*unary)(const void* a, Dtype dtype, void* b);
typedef void (*binary)(const void* a, const void* b, void* c);
typedef int (*unary_int)(int x);
typedef long (*unary_long)(long x);
typedef float (*unary_float)(float x);
typedef double (*unary_double)(double x);

//typedef T (*real_unary_with_mem)(T x, void* mem);
typedef void (*complex_unary)(const complex64 x, complex64* y);
//typedef T (*real_binary)(T x, T y);
typedef void (*complex_binary)(const complex64 x, const complex64 y, complex64* z);

//inline T ssdsp_square(T x) { return x*x; }
//inline T ssdsp_add(T x, T y) { return x+y; }
int taodsp_max_int(int x, int y);
long taodsp_max_long(long x, long y);
float taodsp_max_float(float x, float y);
double taodsp_max_double(double x, double y);

void vec_unary_int(const vector* x, vector* y, unary_int u);
void vec_unary_long(const vector* x, vector* y, unary_long u);
void vec_unary_float(const vector* x, vector* y, unary_float u);
void vec_unary_double(const vector* x, vector* y, unary_double u);
void vec_binary(const vector* x, const vector* y, vector* z, binary b);
void vec_reduce(const vector* x, binary bifunc, void* out);
void vec_unary_accum(const vector* x, unary unary_func, binary accum_func, void* out);

void sum2d(const vector* x, int C, vector* y, ReduceAxis axis);
void mat_sum(const matrix* m, vector* y, int axis);
void sum(const vector* x, void* a);
void mean(const vector* x, void* a);

void vec_max(const vector* x, void* a);

#endif //TAOSBEECHDSP_REDUCE_H
