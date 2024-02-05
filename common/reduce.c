//
// Created by wentao on 2024/1/17.
//

#include "reduce.h"
#include "memory.h"


#define SIMPLE_SWITCH_CODE(op, dtype, ...) \
switch (dtype) {\
    case Int32: op##_int(__VA_ARGS__); break;\
    case Int64: op##_long(__VA_ARGS__); break;\
    case Float32: op##_float(__VA_ARGS__); break;\
    case Float64: op##_double(__VA_ARGS__); break;\
    case Complex64:  op##_complex64(__VA_ARGS__); break;\
    case Complex128: op##_complex128(__VA_ARGS__); break;\
    default: ;\
}

#define COMPARABLE_SWITCH_CODE(op, dtype, ...) \
switch (dtype) {\
    case Int32: op##_int(__VA_ARGS__); break;\
    case Int64: op##_long(__VA_ARGS__); break;\
    case Float32: op##_float(__VA_ARGS__); break;\
    case Float64: op##_double(__VA_ARGS__); break;\
    default: ;\
}

void vec_unary(const vector* x, vector* y, unary u) {
    const char* x_data = (char*)(x->data);
    char* y_data = (char*)(y->data);
    int x_step = x->elemsize;
    int y_step = y->elemsize;
    for (int i=0; i<x->size; i++) {
        u(x_data, x->dtype, y_data);
        x_data += x_step;
        y_data += y_step;
    }
}

#define GEN__vec_unary(T)\
void vec_unary_##T(const vector* x, vector* y, unary_##T u) {\
    const T* x_data = (T*)(x->data);\
    T* y_data = (T*)(y->data);\
    int x_step = x->stride;\
    int y_step = y->stride;\
    for (int i=0; i<x->size; i++) {\
        *y_data = u(*x_data);\
        x_data += x_step;\
        y_data += y_step;\
    }\
}
GEN__vec_unary(int)
GEN__vec_unary(long)
GEN__vec_unary(float)
GEN__vec_unary(double)

void vec_binary(const vector* x, const vector* y, vector* z, binary b) {
    const char* x_data = (char*)(x->data);
    const char* y_data = (char*)(y->data);
    char* z_data = (char*)(z->data);
    int x_step = x->elemsize;
    int y_step = y->elemsize;
    int z_step = z->elemsize;
    for (int i=0; i<x->size; i++) {
        b(x_data, y_data, z_data);
        x_data += x_step;
        y_data += y_step;
        z_data += z_step;
    }
}

void vec_reduce(const vector* x, binary bifuc, void* out) {
    const char* x_data = (char*)(x->data);
    int elemsize = x->elemsize;
    if (x->size == 1) {
        memcpy(out, x->data, x->elemsize);
        return;
    }
    bifuc(x_data, x_data+elemsize, out);
    x_data += elemsize*2;
    for (int i=2; i<x->size; i++) {
        bifuc(x_data, out, out);
        x_data += elemsize;
    }
}

//void vec_unary_accum(const vector* x, unary unary_func, binary accum_func, void* out) {
//    char* x_data = (char*)(x->data);
//    char tmp[x->elemsize];
//    for (int i=0; i<x->size; i++) {
//        unary_func(x_data, tmp);
//        accum_func(tmp, out, out);
//    }
//}

//void vec_complex_unary(const vector* x, vector* y, complex_unary cb) {
//    complex* x_data = GET_COMPLEX_DATA(x);
//    for (int i=0; i<x->size; i++) {
//        complex tmp;
//        cb(x_data[i], &tmp);
//        vec_set_complex(y, i, tmp.re, tmp.im);
//    }
//}

#define GEN__sum(T) \
void sum_##T(const vector* x, void* a) {\
    const T* x_data = (T*)(x->data);\
    T a_ = T##_zero;\
    for (int i=0; i<x->size; i++) T##_add(a_, x_data[i], a_);\
    T* ap = (T*)a;\
    *ap = a_;\
}
GEN__sum(int)
GEN__sum(long)
GEN__sum(float)
GEN__sum(double)
GEN__sum(complex64)
GEN__sum(complex128)
// wrapper function
void sum(const vector* x, void* a) {
    SIMPLE_SWITCH_CODE(sum, x->dtype, x, a)
}

void mean(const vector* x, void* m) {
    SIMPLE_SWITCH_CODE(sum, x->dtype, x, m)
    switch (x->dtype) {
        case Int32:
        case Float32:
        case Complex64: {
            float* m_ = (float*)m;
            *m_ = *m_ / x->size;
            return;
        }
        case Int64:
        case Float64:
        case Complex128:{
            double* m_ = (double*)m;
            *m_ = *m_ / x->size;
            return;
        }
    }
}

#define GEN__sum2d(T) \
void sum2d_##T(const vector* x, int C, vector* y, ReduceAxis axis) {\
    int R = x->size / C;\
    const T* x_data = (T*)(x->data);\
    T* y_data = (T*)(y->data);\
    T zero = T##_zero;\
    vec_fill_##T(y, zero);\
    int x_stride = x->stride;\
    int y_stride = y->stride;\
    if (axis == AlongRow) {\
        for (int r=0; r<R; r++) {\
            y_data = (T*)(y->data);\
            for (int c=0; c<C; c++) {\
                T##_add(*y_data, *x_data, *y_data);\
                x_data += x_stride;\
                y_data += y_stride;\
            }\
        }\
    } else if (axis == AlongCol) {\
        for (int r=0; r<R; r++) {\
            for (int c=0; c<C; c++) {\
                T##_add(*y_data, *x_data, *y_data);\
                x_data += x_stride;\
            }\
            y_data += y_stride;\
        }\
    }\
}
GEN__sum2d(int)
GEN__sum2d(long)
GEN__sum2d(float)
GEN__sum2d(double)
GEN__sum2d(complex64)
GEN__sum2d(complex128)
// wrapper function
void sum2d(const vector* x, int C, vector* y, ReduceAxis axis) {
    check(x->dtype==y->dtype, "Inconsistent data type of arguments of sum2d.");
    SIMPLE_SWITCH_CODE(sum2d, x->dtype, x, C, y, axis)
}

void mat_sum(const matrix* m, vector* y, int axis) {
    check(y->size==((axis==AlongCol)?m->row:m->col), "Inconsistent data shape of arguments of mat_sum.");
    vector_view vw;
    mat_flatten(m, &vw);
    sum2d(&vw, m->col, y, axis);
}

#define GEN__vec_max(T) \
void vec_max_##T(const vector* x, void* a) {\
    T* x_data = (T*)(x->data);\
    T tmp = x_data[0];\
    for (int i=1; i<x->size; i++) if (tmp < x_data[i]) tmp = x_data[i];\
    T* ap = (T*)a;\
    *ap = tmp;\
}
GEN__vec_max(int)
GEN__vec_max(long)
GEN__vec_max(float)
GEN__vec_max(double)

void vec_max(const vector* x, void* a) {
    check(IS_VEC_COMPARABLE(*x), "Invalid data type for .");
    if (x->size == 0) return;
    COMPARABLE_SWITCH_CODE(vec_max, x->dtype, x, a);
}
//
//int taodsp_max_int(int x, int y) { if (x > y) return x; return y; }
//long taodsp_max_long(long x, long y) { if (x > y) return x; return y; }
//float taodsp_max_float(float x, float y) { if (x > y) return x; return y; }
//double taodsp_max_double(double x, double y) { if (x > y) return x; return y; }
