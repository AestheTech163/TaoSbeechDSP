//
// Created by wentao on 2024/1/28.
//
#ifndef TAOSBEECHDSP_VECTOR_H
#define TAOSBEECHDSP_VECTOR_H

#include "arch.h"
#include "numeric.h"

typedef struct vector_ {
    Dtype dtype;
    int size;
    int elemsize;
    int stride;
    int byte_stride;
//    char data[0];
    void* data;
} vector;

typedef struct matrix_ {
    Dtype dtype;
    int row;
    int col;
    int elemsize;
    int byte_stride;
    int stride;
    int col_stride;
    void* data;
} matrix;

typedef vector vector_view;
typedef matrix matrix_view;

#define IS_VEC_COMPARABLE(v) IS_COMPARABLE((v).dtype)
#define IS_VEC_COMPLEX(v) IS_COMPLEX((v).dtype)
#define IS_VEC_REAL(v) IS_REAL((v).dtype)
#define IS_VEC_VALID_DTYPE(v) IS_VALID_DTYPE((v).dtype)
#define VEC_SIZE(v) ((v).size)

vector* vec_alloc(int size, Dtype dtype);
void vec_free(vector* v);
void batch_vec_free(vector* v, ...);
#define BATCH_vec_free(v, ...) batch_vec_free(v, __VA_ARGS__, NULL)
vector* zeros(int size, Dtype dtype);
vector* ones(int size, Dtype dtype);
vector* zeros_like(const vector* x);
vector* ones_like(const vector* x);
vector* arange(int end, Dtype dtype);
vector* arange2(int start, int end, int step, Dtype dtype);

extern void vec_view_construct(const vector* v, int offset, int size, vector_view* vw);
extern void vec_view_construct2(const vector* v, int offset, int size, int stride, vector_view* vw);
extern void vec_view_destruct(vector_view* v);

extern void vec_fill_int(vector*, int);
extern void vec_fill_long(vector*, long);
extern void vec_fill_float(vector*, float);
extern void vec_fill_double(vector*, double);
extern void vec_fill_complex64(vector*, complex64);
extern void vec_fill_complex128(vector*, complex128);
extern void vec_fill_complex64_parts(vector*, float, float);
extern void vec_fill_complex128_parts(vector*, double, double);
extern void vec_fill(vector* x, const void* a);

#define vec_elem_set_char(x, i, a) (((char*)((x)->data))[i] = (a))
#define vec_elem_set_short(x, i, a) (((short*)((x)->data))[i] = (a))
#define vec_elem_set_int(x, i, a) (((int*)((x)->data))[i] = (a))
#define vec_elem_set_long(x, i, a) (((long*)((x)->data))[i] = (a))
#define vec_elem_set_float(x, i, a) (((float*)((x)->data))[i] = (a))
#define vec_elem_set_double(x, i, a) (((double*)((x)->data))[i] = (a))
#define vec_elem_set_complex64(x, i, a) (((complex64*)((x)->data))[i] = (a))
#define vec_elem_set_complex128(x, i, a) (((complex128*)((x)->data))[i] = (a))
#define vec_elem_set_complex64_parts(x, i, real, imag) {complex64 __tmp={real, imag};((complex64*)((x)->data))[i] = __tmp;}
#define vec_elem_set_complex128_parts(x, i, real, imag) {complex128 __tmp={real, imag};((complex128*)((x)->data))[i] = __tmp;}
void vec_elem_set(vector* x, int i, const void* a);

#define vec_elem_incr_short(x, i, a) (((short*)((x)->data))[i] += (a))
#define vec_elem_incr_int(x, i, a) (((int*)((x)->data))[i] += (a))
#define vec_elem_incr_long(x, i, a) (((long*)((x)->data))[i] += (a))
#define vec_elem_incr_float(x, i, a) (((float*)((x)->data))[i] += (a))
#define vec_elem_incr_double(x, i, a) (((double*)((x)->data))[i] += (a))
#define vec_elem_incr_complex64(x, i, a) {complex64* x_d=(complex64*)((x)->data);complex64_add_assign(x_d[i],(a));}
#define vec_elem_incr_complex128(x, i, a) {complex128* x_d = (complex128*)((x)->data); complex128_add_assign(x_d[i],(a));}
#define vec_elem_incr_complex64_parts(x, i, real, imag) {complex64* x_d=(complex64*)((x)->data);x_d[i].re+=(real);x_d[i].im+=(imag);}
#define vec_elem_incr_complex128_parts(x, i, real, imag) {complex128* x_d=(complex128*)((x)->data);x_d[i].re+=(real);x_d[i].im+=(imag);}
//void vec_elem_incr(vector* x, int i, const void* a);

#define vec_elem_get_short(x, i) (((short*)((x)->data))[i])
#define vec_elem_get_int(x, i) (((int*)((x)->data))[i])
#define vec_elem_get_long(x, i) (((long*)((x)->data))[i])
#define vec_elem_get_float(x, i) (((float*)((x)->data))[i])
#define vec_elem_get_double(x, i) (((double*)((x)->data))[i])
#define vec_elem_get_complex64(x, i) (((complex64*)((x)->data))[i])
#define vec_elem_get_complex128(x, i) (((complex128*)((x)->data))[i])
void vec_elem_get(const vector* x, int i, void* a);

/* 内积运算
 * 对于实数 Int32,Int64,Float32,Float64:
 *      dot(x,y) = x1*y1+x2*y2+...+xn*yn;
 * 对于复数 Complex64,Complex128:
 *      dot(x, y) = x1*conj(y1)+x2*conj(y2)+...+xn*conj(yn);
 */
void vec_dot(const vector* x, const vector* y, void* a);
void arr_dot_int(const int* x, const int* y, int size, int xstride, int ystride, void* xy);
void arr_dot_long(const long* x, const long* y, int size, int xstride, int ystride, void* xy);
void arr_dot_float(const float* x, const float* y, int size, int xstride, int ystride, void* xy);
void arr_dot_double(const double* x, const double* y, int size, int xstride, int ystride, void* xy);
void arr_dot_complex64(const complex64* x, const complex64* y, int size, int xstride, int ystride, void* xy);
void arr_dot_complex128(const complex128* x, const complex128* y, int size, int xstride, int ystride, void* xy);

void vec_add_scalar(const vector *x, const void *a, vector *y);
void vec_sub_scalar(const vector *x, const void *a, vector *y);
void vec_mul_scalar(const vector *x, const void *a, vector *y);
void vec_div_scalar(const vector *x, const void *a, vector *y);

void scalar_sub_vec(const void *a, const vector *x, vector *y);
void scalar_div_vec(const void *a, const vector *x, vector *y);

void vec_less(const vector* x, const vector* y, vector* z);
void vec_greater(const vector* x, const vector* y, vector* z);
bool vec_all(const vector* x);
bool vec_any(const vector* x);

void vec_maximum(const vector* x, const vector* y, vector* z);
void vec_minimum(const vector* x, const vector* y, vector* z);

void vec_add(const vector *x, const vector *y, vector *z);
void vec_sub(const vector *x, const vector *y, vector *z);
void vec_mul(const vector *x, const vector *y, vector *z);
void vec_div(const vector *x, const vector *y, vector *z);
void vec_mul_accum(const vector* x, const vector* y, vector* z);
void vec_mul_scalar_accum(const vector* x, const void* a, vector* y);

void vec_reverse(vector* x);
void vec_conj(const vector* x, vector* y);
void vec_real_part(const vector* x, vector_view* y);
void vec_imag_part(const vector* x, vector_view* y);
void vec_abs2(const vector* x, vector* y);
void vec_norm2(const vector* x, void* a);
void vec_norm(const vector* x, void* a);

// basic functions vectorization
void vec_abs(const vector* x, vector* y);
void vec_sqrt(const vector* x, vector* y);
void vec_exp(const vector* x, vector* y);
void vec_log(const vector* x, vector* y);

extern void vec_clip(vector* v, double inf, double sup);
void vec_copy_slice(const vector* x, int xoffset, int size, vector* y, int yoffset);
void vec_copy(const vector* x, vector* y);
void vec_print(const vector* v);

void mat_view_construct(const matrix* v, matrix_view* vw);
void mat_view_from_vec_slice(const vector* v, int offset, int size, int col, matrix_view* vw);
void mat_view_from_vec(const vector* v, int col, matrix_view* vw);
void mat_view_destruct(matrix_view* v);
void mat_flatten(const matrix* m, vector_view* v);
void mat_get_col(const matrix* m, int r, vector_view* vw);
void mat_get_row(const matrix* m, int r, vector_view* vw);
void mat_print(matrix* m);

void mat_elemwise_add(const matrix *x, const matrix *y, matrix *z);
void mat_elemwise_mul(const matrix *x, const matrix *y, matrix *z);
void mat_elemwise_sub(const matrix *x, const matrix *y, matrix *z);
void mat_elemwise_div(const matrix *x, const matrix *y, matrix *z);

#define auto_vector_construct(vname, vsize, vdtype) \
    int size_bytes_##vname, elemsize_##vname;\
    switch (vdtype) {\
        case Int32:   size_bytes_##vname = (vsize) * sizeof(int); elemsize_##vname = sizeof(int); break;\
        case Int64:   size_bytes_##vname = (vsize) * sizeof(long); elemsize_##vname = sizeof(long); break;\
        case Float32: size_bytes_##vname = (vsize) * sizeof(float); elemsize_##vname = sizeof(float); break;\
        case Float64: size_bytes_##vname = (vsize) * sizeof(double); elemsize_##vname = sizeof(double); break;\
        case Complex64:  size_bytes_##vname = (vsize) * sizeof(complex64); elemsize_##vname = sizeof(complex64); break;\
        case Complex128: size_bytes_##vname = (vsize) * sizeof(complex128); elemsize_##vname = sizeof(complex128); break;\
        default: check(0, "Not implemented data type for auto_vector_construct.");\
    }\
    char arr_##vname[size_bytes_##vname];\
    vname.data=&arr_##vname[0];\
    vname.size=(vsize);\
    vname.dtype=(vdtype);\
    vname.elemsize=elemsize_##vname;\
    vname.stride=1;\
    vname.byte_stride=elemsize_##vname;

#define auto_vector_like(x, y) auto_vector_construct(x, (y)->size, (y)->dtype)

#define auto_vector_construct_int(vname, vsize) \
    int arr_##vname[vsize];\
    vname.data=&arr_##vname[0];\
    vname.size=vsize;\
    vname.dtype=Int32;\
    vname.elemsize=sizeof(int);\
    vname.byte_stride=vname.elemsize;\
    vname.stride=1;
#define auto_vector_construct_long(vname, vsize) \
    long arr_##vname[vsize];\
    vname.data=&arr_##vname[0];\
    vname.size=vsize;\
    vname.dtype=Int64;\
    vname.elemsize=sizeof(long);\
    vname.byte_stride=vname.elemsize;\
    vname.stride=1;
#define auto_vector_construct_float(vname, vsize) \
    float arr_##vname[vsize];\
    vname.data=&arr_##vname[0];\
    vname.size=vsize;\
    vname.dtype=Float32;\
    vname.elemsize=sizeof(float);\
    vname.byte_stride=vname.elemsize;\
    vname.stride=1;
#define auto_vector_construct_double(vname, vsize) \
    double arr_##vname[vsize];\
    vname.data=&arr_##vname[0];\
    vname.size=vsize;\
    vname.dtype=Float64;\
    vname.elemsize=sizeof(double);\
    vname.byte_stride=vname.elemsize;\
    vname.stride=1;
#define auto_vector_construct_complex64(vname, vsize) \
    complex64 arr_##vname[vsize];\
    vname.data=&arr_##vname[0];\
    vname.size=vsize;\
    vname.dtype=Complex64;\
    vname.elemsize=sizeof(complex64);\
    vname.byte_stride=vname.elemsize;\
    vname.stride=1;
#define auto_vector_construct_complex128(vname, vsize) \
    complex128 arr_##vname[vsize];\
    vname.data=&arr_##vname[0];\
    vname.size=vsize;\
    vname.dtype=Complex128;\
    vname.elemsize=sizeof(complex128);\
    vname.byte_stride=vname.elemsize;\
    vname.stride=1;

#define SIMPLE_SWITCH_NO_RETURN(op, dtype, ...) \
switch (dtype) {\
    case Int8:  op##_char(__VA_ARGS__); break;\
    case Int16: op##_short(__VA_ARGS__); break;\
    case Int32: op##_int(__VA_ARGS__); break;\
    case Int64: op##_long(__VA_ARGS__); break;\
    case Float32: op##_float(__VA_ARGS__); break;\
    case Float64: op##_double(__VA_ARGS__); break;\
    case Complex64:  op##_complex64(__VA_ARGS__); break;\
    case Complex128: op##_complex128(__VA_ARGS__); break;\
    default: ;\
}

#define COMPARABLE_SWITCH_NO_RETURN(op, dtype, ...) \
switch (dtype) {\
    case Int8:  op##_char(__VA_ARGS__); break;\
    case Int16: op##_short(__VA_ARGS__); break;\
    case Int32: op##_int(__VA_ARGS__); break;\
    case Int64: op##_long(__VA_ARGS__); break;\
    case Float32: op##_float(__VA_ARGS__); break;\
    case Float64: op##_double(__VA_ARGS__); break;\
    default: ;\
}

#define SWITCH_CODE_AT_LEAST_4BYTES(op, dtype, ...) \
switch (dtype) {\
    case Int32: op##_int(__VA_ARGS__); break;\
    case Int64: op##_long(__VA_ARGS__); break;\
    case Float32: op##_float(__VA_ARGS__); break;\
    case Float64: op##_double(__VA_ARGS__); break;  \
    case Complex64:  op##_complex64(__VA_ARGS__); break;\
    case Complex128: op##_complex128(__VA_ARGS__); break;\
    default: ;\
}
#define SWITCH_DATA_OFFSET(v, offset, vw) \
    switch (v->dtype) {\
        case Int32:   vw->data = (int*)(v->data) + offset*v->stride; break;\
        case Int64:   vw->data = (long*)(v->data) + offset*v->stride; break;\
        case Float32: vw->data = (float*)(v->data) + offset*v->stride; break;\
        case Float64: vw->data = (double*)(v->data) + offset*v->stride; break;\
        case Complex64:  vw->data = (complex64*)(v->data) + offset*v->stride; break;\
        case Complex128: vw->data = (complex128*)(v->data) + offset*v->stride; break;\
    }

#endif //TAOSBEECHDSP_VECTOR_H
