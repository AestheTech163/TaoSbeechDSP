//
// Created by wentao on 2024/1/28.
//
#include "math.h"
#include "stdarg.h"
#include "string.h"

#include "vector.h"
#include "arch.h"
#include "numeric.h"
//#include "reduce.h"

#define VECTORIZE_PREDICATE(Predicate, T) \
void vec_##Predicate##_##T(const vector* x, const vector* y, vector* z) {\
    const T* x_data = x->data;\
    const T* y_data = y->data;\
    char* z_data = z->data;\
    int xstride = x->stride, ystride = y->stride, zstride = z->stride;\
    for (int i=0; i<x->size; i++) {\
        *z_data = T##_##Predicate(*x_data, *y_data);\
        x_data += xstride;\
        y_data += ystride;\
        z_data += zstride;\
    }\
}
VECTORIZE_PREDICATE(less, char)
VECTORIZE_PREDICATE(less, short)
VECTORIZE_PREDICATE(less, int)
VECTORIZE_PREDICATE(less, long)
VECTORIZE_PREDICATE(less, float)
VECTORIZE_PREDICATE(less, double)
VECTORIZE_PREDICATE(greater, char)
VECTORIZE_PREDICATE(greater, short)
VECTORIZE_PREDICATE(greater, int)
VECTORIZE_PREDICATE(greater, long)
VECTORIZE_PREDICATE(greater, float)
VECTORIZE_PREDICATE(greater, double)
// wrapper function
void vec_less(const vector* x, const vector* y, vector* z) {
    check(IS_VEC_COMPARABLE(*x)&&IS_VEC_COMPARABLE(*y), "Not implemented data type for vec_less.");
    check(x->size==y->size, "Size of x and y should be equal.");
    check(z->dtype==Int8, "Data type of z should be Int8.");
    COMPARABLE_SWITCH_NO_RETURN(vec_less, x->dtype, x, y, z)
}
void vec_greater(const vector* x, const vector* y, vector* z) {
    check(IS_VEC_COMPARABLE(*x)&&IS_VEC_COMPARABLE(*y), "Not implemented data type for vec_greater.");
    check(x->size==y->size, "Size of x and y should be equal.");
    check(z->dtype==Int8, "Data type of z should be Int8.");
    COMPARABLE_SWITCH_NO_RETURN(vec_greater, x->dtype, x, y, z)
}

bool vec_all(const vector* x) {
    switch (x->dtype) {
#define TMP_VEC_ALL(T) {T* x_d=x->data; for(int i=0;i<x->size;i++) if (!x_d[i]) return 0; return 1;}
        case Int8:  TMP_VEC_ALL(char); break;
        case Int16: TMP_VEC_ALL(short); break;
        case Int32: TMP_VEC_ALL(int); break;
        case Int64: TMP_VEC_ALL(long); break;
        case Float32: TMP_VEC_ALL(float); break;
        case Float64: TMP_VEC_ALL(double); break;
#undef TMP_VEC_ALL
        default: check(0, "Invalid data type for vec_all.");
    }
    return 0;  // never reached
}

bool vec_any(const vector* x) {
    switch (x->dtype) {
#define TMP_VEC_ALL(T) {T* x_d=x->data; for(int i=0;i<x->size;i++) if (x_d[i]) return 1; return 0;}
        case Int8:  TMP_VEC_ALL(char); break;
        case Int16: TMP_VEC_ALL(short); break;
        case Int32: TMP_VEC_ALL(int); break;
        case Int64: TMP_VEC_ALL(long); break;
        case Float32: TMP_VEC_ALL(float); break;
        case Float64: TMP_VEC_ALL(double); break;
#undef TMP_VEC_ALL
        default: check(0, "Invalid data type for vec_all.");
    }
    return 0;  // never reached
}

#define VECTORIZE_BINARY_FUNC(Binary, T) \
void vec_##Binary##_##T(const vector* x, const vector* y, vector* z) {\
    const T* x_data = x->data;\
    const T* y_data = y->data;\
    T* z_data = z->data;\
    int xstride = x->stride, ystride = y->stride, zstride = z->stride;\
    for (int i=0; i<x->size; i++) {\
        T##_##Binary(*z_data, *x_data, *y_data);\
        x_data += xstride;\
        y_data += ystride;\
        z_data += zstride;\
    }\
}
VECTORIZE_BINARY_FUNC(maximum, char)
VECTORIZE_BINARY_FUNC(maximum, short)
VECTORIZE_BINARY_FUNC(maximum, int)
VECTORIZE_BINARY_FUNC(maximum, long)
VECTORIZE_BINARY_FUNC(maximum, float)
VECTORIZE_BINARY_FUNC(maximum, double)
VECTORIZE_BINARY_FUNC(minimum, char)
VECTORIZE_BINARY_FUNC(minimum, short)
VECTORIZE_BINARY_FUNC(minimum, int)
VECTORIZE_BINARY_FUNC(minimum, long)
VECTORIZE_BINARY_FUNC(minimum, float)
VECTORIZE_BINARY_FUNC(minimum, double)

// wrapper function
void vec_maximum(const vector* x, const vector* y, vector* z) {
    check(IS_VEC_COMPARABLE(*x), "Not implemented data type for vec_maximum.");
    check(x->size==y->size, "Size of x and y should be equal.");
    COMPARABLE_SWITCH_NO_RETURN(vec_maximum, x->dtype, x, y, z)
}
void vec_minimum(const vector* x, const vector* y, vector* z) {
    check(IS_VEC_COMPARABLE(*x), "Not implemented data type for vec_minimum.");
    check(x->size==y->size, "Size of x and y should be equal.");
    COMPARABLE_SWITCH_NO_RETURN(vec_minimum, x->dtype, x, y, z)
}

#define GEN__vec_reverse(T)\
void vec_reverse_##T(vector* x) {\
    T tmp;\
    T* x_data = (T*)(x->data);\
    for (int i=0; i<x->size/2; i++) {\
        tmp = x_data[i];\
        x_data[i] = x_data[x->size-1-i];\
        x_data[x->size-1-i] = tmp;\
    }\
}
GEN__vec_reverse(char);
GEN__vec_reverse(short);
GEN__vec_reverse(int);
GEN__vec_reverse(long);
GEN__vec_reverse(float);
GEN__vec_reverse(double);
GEN__vec_reverse(complex64);
GEN__vec_reverse(complex128);
// wrapper function
void vec_reverse(vector* x) {
    SIMPLE_SWITCH_NO_RETURN(vec_reverse, x->dtype, x)
}

#define GEN__arr_dot_real(T)\
void arr_dot_##T(const T* x, const T* y, int size, int xstride, int ystride, void* xy) {\
    T s = T##_ZERO;\
    for (int i=0; i<size; i++) {\
        T tmp;\
        T##_mul(tmp, *x, *y);\
        T##_add_assign(s, tmp);\
        x += xstride;\
        y += ystride;\
    }\
    T##_assign(*(T*)xy, s);\
}
#define GEN__arr_dot_complex(T)\
void arr_dot_##T(const T* x, const T* y, int size, int xstride, int ystride, void* xy) { \
    T s = T##_ZERO;\
    for (int i=0; i<size; i++) {\
        T x_conj, tmp;\
        T##_conj(x_conj, *x);\
        T##_mul(tmp, x_conj, *y);\
        T##_add_assign(s, tmp);\
        x += xstride;\
        y += ystride;\
    }\
    *(T*)xy = s;\
}
GEN__arr_dot_real(char)
GEN__arr_dot_real(short)
GEN__arr_dot_real(int)
GEN__arr_dot_real(long)
GEN__arr_dot_real(float)
GEN__arr_dot_real(double)
GEN__arr_dot_complex(complex64)
GEN__arr_dot_complex(complex128)
// vec_dot wrapper function
void vec_dot(const vector* x, const vector* y, void* a) {
    SIMPLE_SWITCH_NO_RETURN(arr_dot, x->dtype, x->data, y->data, x->size, x->stride, y->stride, a)
}

void vec_conj(const vector* x, vector* y) {
    check(x->size==y->size, "<vec_conj> size of x and y should be equal.\n");
    switch (x->dtype) {
        case Complex64: {
            const complex64* x_data = (complex64*)(x->data);
            complex64* y_data = (complex64*)(y->data);
            for (int i=0; i<x->size; i++) {
                complex64_conj(*y_data, *x_data);
                x_data += x->stride;
                y_data += y->stride;
            }
            return;
        }
        case Complex128: {
            const complex128* x_data = (complex128*)(x->data);
            complex128* y_data = (complex128*)(y->data);
            for (int i=0; i<x->size; i++) {
                complex128_conj(*y_data, *x_data);
                x_data += x->stride;
                y_data += y->stride;
            }
            return;
        }
        case Int32:
        case Int64:
        case Float32:
        case Float64:
            if (x->data != y->data) {
                // memcpy(y->data, x->data, sizeof(T)*x->size);
                const char* x_data = (char*)(x->data);
                char* y_data = (char*)(y->data);
                for (int i=0; i<x->size; i++) {
                    memcpy(y_data, x_data, x->elemsize);
                    x_data += x->byte_stride;
                    y_data += y->byte_stride;
                }
            }
            return;
        default: check(0, "Not implemented data type for vec_conj.");
    }
}

#define GEN__vec_clip(T)\
void vec_clip_##T(vector* v, double inf, double sup) {\
    check(IS_VEC_COMPARABLE(*v), "Invalid data type for vec_clip.");\
    T* v_data = (T*)(v->data);\
    for (int i = 0; i < v->size; i++)\
        v_data[i] = CLIP(v_data[i], inf, sup);\
}
GEN__vec_clip(char)
GEN__vec_clip(short)
GEN__vec_clip(int)
GEN__vec_clip(float)
GEN__vec_clip(long)
GEN__vec_clip(double)
void vec_clip(vector* v, double inf, double sup) {
    COMPARABLE_SWITCH_NO_RETURN(vec_clip, v->dtype, v, inf, sup)
}

#define GEN__vec_op_scalar(op, T)\
void vec_##op##_scalar_##T(const vector* x, const void* a, vector* y) {\
    const T* x_data = (T*)(x->data);\
    T* y_data = (T*)(y->data);\
    const T a_ = *(T*)a;\
    for (int i=0; i<x->size; i++)\
        T##_##op(y_data[i], x_data[i], a_);\
}
GEN__vec_op_scalar(add, int)
GEN__vec_op_scalar(add, long)
GEN__vec_op_scalar(add, float)
GEN__vec_op_scalar(add, double)
GEN__vec_op_scalar(add, complex64)
GEN__vec_op_scalar(add, complex128)

GEN__vec_op_scalar(sub, int)
GEN__vec_op_scalar(sub, long)
GEN__vec_op_scalar(sub, float)
GEN__vec_op_scalar(sub, double)
GEN__vec_op_scalar(sub, complex64)
GEN__vec_op_scalar(sub, complex128)

GEN__vec_op_scalar(mul, int)
GEN__vec_op_scalar(mul, long)
GEN__vec_op_scalar(mul, float)
GEN__vec_op_scalar(mul, double)
GEN__vec_op_scalar(mul, complex64)
GEN__vec_op_scalar(mul, complex128)

GEN__vec_op_scalar(div, int)
GEN__vec_op_scalar(div, long)
GEN__vec_op_scalar(div, float)
GEN__vec_op_scalar(div, double)
GEN__vec_op_scalar(div, complex64)
GEN__vec_op_scalar(div, complex128)

#define GEN_WRAPPER__vec_op_scalar(op)\
void vec_##op##_scalar(const vector* x, const void* a, vector* y) {\
    check(x->size > 0 && y->size > 0, "Invalid vector size.");\
    check(((x->size % y->size)==0) || ((y->size % x->size)==0), "Invalid vector size.");\
    SWITCH_CODE_AT_LEAST_4BYTES(vec_##op##_scalar, x->dtype, x, a, y)\
}
GEN_WRAPPER__vec_op_scalar(add)
GEN_WRAPPER__vec_op_scalar(sub)
GEN_WRAPPER__vec_op_scalar(mul)
GEN_WRAPPER__vec_op_scalar(div)

#define GEN__scalar_op_vec(op, T)\
void scalar_##op##_vec_##T(const void* a, const vector* x, vector* y) {\
    const T* x_data = (T*)(x->data);\
    T* y_data = (T*)(y->data);\
    const T a_ = *(T*)a;\
    for (int i=0; i<x->size; i++)\
        T##_##op(y_data[i], a_, x_data[i]);\
}
// Only gen Non-Commutative operations
GEN__scalar_op_vec(sub, char)
GEN__scalar_op_vec(sub, short)
GEN__scalar_op_vec(sub, int)
GEN__scalar_op_vec(sub, long)
GEN__scalar_op_vec(sub, float)
GEN__scalar_op_vec(sub, double)
GEN__scalar_op_vec(sub, complex64)
GEN__scalar_op_vec(sub, complex128)

GEN__scalar_op_vec(div, char)
GEN__scalar_op_vec(div, short)
GEN__scalar_op_vec(div, int)
GEN__scalar_op_vec(div, long)
GEN__scalar_op_vec(div, float)
GEN__scalar_op_vec(div, double)
GEN__scalar_op_vec(div, complex64)
GEN__scalar_op_vec(div, complex128)

#define GEN_WRAPPER__scalar_op_vec(op)\
void scalar_##op##_vec(const void* a, const vector* x, vector* y) {\
    check(x->size > 0 && y->size > 0, "Invalid vector size.");\
    check(((x->size % y->size)==0) || ((y->size % x->size)==0), "Invalid vector size.");\
    SIMPLE_SWITCH_NO_RETURN(scalar_##op##_vec, x->dtype, a, x, y)\
}
GEN_WRAPPER__scalar_op_vec(sub)
GEN_WRAPPER__scalar_op_vec(div)

#define GEN__vec_elem_wise_op(op, T)\
void vec_##op##_##T(const vector* x, const vector* y, vector* z) {\
    const T* x_data = (T*)(x->data);\
    const T* y_data = (T*)(y->data);\
    T* z_data = (T*)(z->data);\
    int xstride = x->stride;\
    int ystride = y->stride;\
    int zstride = z->stride;\
    if (x->size >= y->size) {\
        for (int b=0; b<x->size/y->size; b++) {\
            for (int i = 0; i < y->size; i++) {\
                T##_##op(*z_data, *x_data, *y_data); \
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            };\
            y_data = (T*)(y->data);\
        }\
    } else {\
        for (int b=0; b<y->size/x->size; b++) {\
            for (int i=0; i<x->size; i++) {\
                T##_##op(*z_data, *x_data, *y_data); \
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            }\
            x_data = (T*)(x->data);\
        }\
    }\
}
#ifdef VecGen_int
GEN__vec_elem_wise_op(sub, int);
GEN__vec_elem_wise_op(add, int);
GEN__vec_elem_wise_op(mul, int);
GEN__vec_elem_wise_op(div, int);
#endif
#ifdef VecGen_long
GEN__vec_elem_wise_op(sub, long);
GEN__vec_elem_wise_op(add, long);
GEN__vec_elem_wise_op(mul, long);
GEN__vec_elem_wise_op(div, long);
#endif
#ifdef VecGen_float
GEN__vec_elem_wise_op(sub, float);
GEN__vec_elem_wise_op(add, float);
GEN__vec_elem_wise_op(mul, float);
GEN__vec_elem_wise_op(div, float);
#endif
#ifdef VecGen_double
GEN__vec_elem_wise_op(sub, double);
GEN__vec_elem_wise_op(add, double);
GEN__vec_elem_wise_op(mul, double);
GEN__vec_elem_wise_op(div, double);
#endif
#ifdef VecGen_complex64
GEN__vec_elem_wise_op(sub, complex64);
GEN__vec_elem_wise_op(add, complex64);
GEN__vec_elem_wise_op(mul, complex64);
GEN__vec_elem_wise_op(div, complex64);
#endif
#ifdef VecGen_complex128
GEN__vec_elem_wise_op(sub, complex128);
GEN__vec_elem_wise_op(add, complex128);
GEN__vec_elem_wise_op(mul, complex128);
GEN__vec_elem_wise_op(div, complex128);
#endif

#define GEN_MIX__vec_elem_wise_op(op, T1, T2, T3)\
void vec_##T1##_##op##_##T2(const vector* x, const vector* y, vector* z) {\
    const T1* x_data = (T1*)(x->data);\
    const T2* y_data = (T2*)(y->data);\
    T3* z_data = (T3*)(z->data);\
    int xstride = x->stride;\
    int ystride = y->stride;\
    int zstride = z->stride;\
    if (x->size >= y->size) {\
        for (int b=0; b<x->size/y->size; b++) {\
            y_data = (T2*)(y->data);\
            for (int i = 0; i < y->size; i++) {\
                T1##_##op##_##T2(*z_data, *x_data, *y_data);\
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            }\
        }\
    } else {\
        for (int b=0; b<y->size/x->size; b++) {\
            x_data = (T1*)(x->data);\
            for (int i = 0; i < x->size; i++) {\
                T1##_##op##_##T2(*z_data, *x_data, *y_data);\
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            }\
        }\
    }\
}
#if defined(VecGen_complex64)&&defined(VecGen_float)
// Commutative operations
GEN_MIX__vec_elem_wise_op(add, complex64, float, complex64);
GEN_MIX__vec_elem_wise_op(mul, complex64, float, complex64);
GEN_MIX__vec_elem_wise_op(add, float, complex64, complex64);
GEN_MIX__vec_elem_wise_op(mul, float, complex64, complex64);
// Non-Commutative operations
GEN_MIX__vec_elem_wise_op(sub, complex64, float, complex64);
GEN_MIX__vec_elem_wise_op(div, complex64, float, complex64);
GEN_MIX__vec_elem_wise_op(sub, float, complex64, complex64);
GEN_MIX__vec_elem_wise_op(div, float, complex64, complex64);
#endif
#if defined(VecGen_complex128)&&defined(VecGen_double)
// Commutative operations
GEN_MIX__vec_elem_wise_op(add, complex128, double, complex128);
GEN_MIX__vec_elem_wise_op(mul, complex128, double, complex128);
GEN_MIX__vec_elem_wise_op(add, double, complex128, complex128);
GEN_MIX__vec_elem_wise_op(mul, double, complex128, complex128);
// Non-Commutative operations
GEN_MIX__vec_elem_wise_op(sub, complex128, double, complex128);
GEN_MIX__vec_elem_wise_op(div, complex128, double, complex128);
GEN_MIX__vec_elem_wise_op(sub, double, complex128, complex128);
GEN_MIX__vec_elem_wise_op(div, double, complex128, complex128);
#endif

#define GEN_WRAPPER__vec_op(op)\
void vec_##op(const vector* x, const vector* y, vector* z) {\
    check(MAX(x->size,y->size)==z->size, "Size of z must equal maximum of x and y.");\
    check(x->size>0&&y->size>0&&(((x->size%y->size)==0)||((y->size%x->size)==0)), "Invalid vector size.");\
    switch (x->dtype) {\
        case Complex64: {\
            if (y->dtype==Float32) return vec_complex64_##op##_float(x, y, z);\
            else return vec_##op##_complex64(x, y, z);\
        }\
        case Complex128: {\
            if (y->dtype==Float64) return vec_complex128_##op##_double(x, y, z);\
            else return vec_##op##_complex128(x, y, z);\
        }\
        case Float32:{\
            if (y->dtype==Complex64) return vec_float_##op##_complex64(x, y, z);\
            else return vec_##op##_float(x, y, z);\
        }\
        case Float64: {\
            if (y->dtype==Complex128) return vec_double_##op##_complex128(x, y, z);\
            else return vec_##op##_double(x, y, z);\
        }\
        case Int32: return vec_##op##_int(x, y, z);\
        case Int64: return vec_##op##_long(x, y, z);\
        default: check(0, "Not implemented data type for vec_"#op".");\
    }\
}
GEN_WRAPPER__vec_op(add);
GEN_WRAPPER__vec_op(mul);
GEN_WRAPPER__vec_op(sub);
GEN_WRAPPER__vec_op(div);

#define GEN__vec_mul_accum(T)\
void vec_mul_accum_##T(const vector* x, const vector* y, vector* z) {\
    const T* x_data = (T*)(x->data);\
    const T* y_data = (T*)(y->data);\
    T* z_data = (T*)(z->data);\
    int xstride = x->stride;\
    int ystride = y->stride;\
    int zstride = z->stride;\
    T tmp;\
    if (x->size >= y->size) {\
        for (int b=0; b<x->size/y->size; b++) {\
            y_data = (T*)(y->data);\
            for (int i = 0; i < y->size; i++) {\
                T##_mul(tmp, *x_data, *y_data);\
                T##_add_assign(*z_data, tmp);\
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            };\
        }\
    } else {\
        for (int b=0; b<y->size/x->size; b++) {\
            x_data = (T*)(x->data);\
            for (int i = 0; i < y->size; i++) {\
                T##_mul(tmp, *x_data, *y_data);\
                T##_add_assign(*z_data, tmp);\
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            }\
        }\
    }\
}
GEN__vec_mul_accum(int)
GEN__vec_mul_accum(long)
GEN__vec_mul_accum(float)
GEN__vec_mul_accum(double)
GEN__vec_mul_accum(complex64)
GEN__vec_mul_accum(complex128)

#define GEN_MIX__vec_mul_accum(ComplexType, RealType)\
void vec_##ComplexType##_mul_##RealType##_accum(const vector* x, const vector* y, vector* z) {\
    const ComplexType* x_data = (ComplexType*)(x->data);\
    const RealType* y_data = (RealType*)(y->data);\
    ComplexType* z_data = (ComplexType*)(z->data);\
    int xstride = x->stride;\
    int ystride = y->stride;\
    int zstride = z->stride;\
    ComplexType tmp;\
    if (x->size >= y->size) {\
        for (int b=0; b<x->size/y->size; b++) {\
            y_data = (RealType*)(y->data);\
            for (int i = 0; i < y->size; i++) {\
                ComplexType##_mul_##RealType(tmp, *x_data, *y_data);\
                ComplexType##_add_assign(*z_data, tmp);\
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            }\
        }\
    } else {\
        for (int b=0; b<y->size/x->size; b++) {\
            x_data = (ComplexType*)(x->data);\
            for (int i = 0; i < x->size; i++) {\
                ComplexType##_mul_##RealType(tmp, *x_data, *y_data);\
                ComplexType##_add_assign(*z_data, tmp);\
                x_data += xstride;\
                y_data += ystride;\
                z_data += zstride;\
            }\
        }\
    }\
}
#if defined(VecGen_complex64)&&defined(VecGen_float)
GEN_MIX__vec_mul_accum(complex64 , float)
#endif
#if defined(VecGen_complex128)&&defined(VecGen_double)
GEN_MIX__vec_mul_accum(complex128 , double)
#endif

void vec_mul_accum(const vector* x, const vector* y, vector* z) {
//    SIMPLE_SWITCH_CODE(vec_mul_accum, x->dtype, x, y, z)
    check(IS_VEC_VALID_DTYPE(*x), "Not implemented data type for vec_mul_accum.");
    check(z->size==MAX(x->size,y->size), "Inconsistent data shape of arguments.");
    switch (x->dtype) {
        case Int32: return vec_mul_accum_int(x, y, z);
        case Int64: return vec_mul_accum_long(x, y, z);
        case Float32:
            if (y->dtype == Complex64) return vec_complex64_mul_float_accum(y, x, z);
            return vec_mul_accum_float(x, y, z);
        case Float64:
            if (y->dtype == Complex128) return vec_complex128_mul_double_accum(y, x, z);
            return vec_mul_accum_double(x, y, z);
        case Complex64:
            if (y->dtype == Float32) return vec_complex64_mul_float_accum(x, y, z);
            return vec_mul_accum_complex64(x, y, z);
        case Complex128:
            if (y->dtype == Float64) return vec_complex128_mul_double_accum(x, y, z);
            return vec_mul_accum_complex128(x, y, z);
        default:;
    }
}

#define GEN__vec_mul_scalar_accum(T)\
void vec_mul_scalar_accum_##T(const vector* x, const void* a, vector* y) {\
    check((y->size%x->size)==0, "Size of y should be multiple of x.");\
    const T* x_data = (T*)(x->data);\
    T* y_data = (T*)(y->data);\
    const T s = *(T*)a;\
    T tmp;\
    int xstride = x->stride;\
    int ystride = y->stride;\
    for (int b=0; b<y->size/x->size; b++) {\
        for (int i=0; i<x->size; i++) {\
            T##_mul(tmp, *x_data, s);\
            T##_add_assign(*y_data, tmp);\
            x_data += xstride;\
            y_data += ystride;\
        }\
        x_data = (T*)(x->data);\
    }\
}
GEN__vec_mul_scalar_accum(char)
GEN__vec_mul_scalar_accum(short)
GEN__vec_mul_scalar_accum(int)
GEN__vec_mul_scalar_accum(float)
GEN__vec_mul_scalar_accum(long)
GEN__vec_mul_scalar_accum(double)
GEN__vec_mul_scalar_accum(complex64)
GEN__vec_mul_scalar_accum(complex128)

#define GEN_MIX_vec_mul_scalar_accum(op, ComplexType, RealType)\
void vec_##ComplexType##_mul_scalar_##RealType##_accum(const vector* x, RealType a, vector* y) {\
    check((y->size%x->size)==0, "Invalid vector shape.");\
    ComplexType* y_data = (ComplexType*)(y->data);\
    const ComplexType* x_data;\
    int xstride = x->stride;\
    int ystride = y->stride;\
    ComplexType tmp;\
    for (int b=0; b<y->size/x->size; b++) {\
        x_data = (ComplexType*)(x->data);\
        for (int i=0; i<x->size; i++) {\
            ComplexType##_mul_##RealType(tmp, *x_data, a);\
            ComplexType##_add_assign(*y_data, tmp);\
            x_data += xstride;\
            y_data += ystride;\
        }\
    }\
}
#if defined(VecGen_complex64)&&defined(VecGen_float)
GEN_MIX_vec_mul_scalar_accum(add, complex64 , float);
#endif
#if defined(VecGen_complex128)&&defined(VecGen_double)
GEN_MIX_vec_mul_scalar_accum(add, complex128 , double);
#endif

void vec_mul_scalar_accum(const vector* x, const void* a, vector* y) {
    check(IS_VEC_VALID_DTYPE(*x), "Not implemented data type for vec_mul_scalar_accum.");

    COMPARABLE_SWITCH_NO_RETURN(vec_mul_scalar_accum, x->dtype, x, a, y)

    switch (x->dtype) {
        case Complex64: {
            complex64* a_ = (complex64*)a;
            if (AlmostEqualF(a_->im, 0.f)) return vec_complex64_mul_scalar_float_accum(x, a_->re, y);
            else return vec_mul_scalar_accum_complex64(x, a, y);
        }
        case Complex128: {
            complex64* a_ = (complex64*)a;
            if (AlmostEqualD(a_->im, 0.)) return vec_complex128_mul_scalar_double_accum(x, a_->re, y);
            else return vec_mul_scalar_accum_complex128(x, a, y);
        }
        default:;
    }
}

#define VECTORIZE_UNARY(func, T)\
void vec_##func##_##T(const vector* x, vector* y) {\
    T* x_data = (T*)(x->data);\
    T* y_data = (T*)(y->data);\
    int xstride = x->stride;\
    int ystride = y->stride;\
    for (int i=0; i<x->size; i++) {\
        func##_##T(*y_data, *x_data);\
        x_data += xstride;\
        y_data += ystride;\
    }\
}
VECTORIZE_UNARY(sqrt, float)
VECTORIZE_UNARY(sqrt, double)
VECTORIZE_UNARY(exp, float)
VECTORIZE_UNARY(exp, double)
VECTORIZE_UNARY(exp, complex64)
VECTORIZE_UNARY(exp, complex128)
VECTORIZE_UNARY(log, float)
VECTORIZE_UNARY(log, double)
// wrapper functions
void vec_sqrt(const vector* x, vector* y) {
    switch (x->dtype) {
        case Float32: return vec_exp_float(x, y);
        case Float64: return vec_exp_double(x, y);
        default: check(0, "Not implemented data type for vec_sqrt.");
    }
}

void vec_exp(const vector* x, vector* y) {
    switch (x->dtype) {
        case Float32: return vec_exp_float(x, y);
        case Float64: return vec_exp_double(x, y);
        case Complex64: return vec_exp_complex64(x, y);
        case Complex128: return vec_exp_complex128(x, y);
        default: check(0, "Not implemented data type for vec_exp.");
    }
}

void vec_log(const vector* x, vector* y) {
    switch (x->dtype) {
        case Float32: return vec_exp_float(x, y);
        case Float64: return vec_exp_double(x, y);
        default: check(0, "Not implemented data type for vec_log.");
    }
}


#define GEN__arr_print(T, fmt)\
void arr_print_##T(T* arr, int N, const char* prefix, const char* sufix) {\
    if (prefix) printf("%s", prefix);\
    printf("[");\
    int i=0;\
    if (N<=VEC_PRINT_LIMIT_REAL-1) for (; i<N-1; i++) printf(#fmt",", arr[i]); \
    else {\
        for (; i<VEC_PRINT_LIMIT_REAL/2; i++) printf(#fmt",", arr[i]);    \
        printf("...,");\
        for (i=N-VEC_PRINT_LIMIT_REAL/2; i<N-1; i++) printf(#fmt",", arr[i]);    \
    }\
    printf(#fmt"]", arr[i]);\
    if (sufix) printf("%s", sufix);\
}
GEN__arr_print(char, %d);
GEN__arr_print(short, %d);
GEN__arr_print(int, %d);
GEN__arr_print(long, %ld);
GEN__arr_print(float, %e);
GEN__arr_print(double, %e);

#define GEN__carr_print(ComplexType, fmt)\
void arr_print_##ComplexType(ComplexType* arr, int N, const char* prefix, const char* sufix) { \
    if (prefix) printf("%s", prefix);\
    printf("[");\
    int i=0;\
    if (N<=VEC_PRINT_LIMIT_COMPLEX-1) for (; i<N-1; i++) {\
        if (arr[i].im >= 0) printf(fmt"+j"fmt", ", arr[i].re, arr[i].im);\
        else printf(fmt"-j"fmt", ", arr[i].re, -arr[i].im);\
    }\
    else {\
        for (; i<VEC_PRINT_LIMIT_COMPLEX/2; i++) {\
            if (arr[i].im >= 0) printf(fmt"+j"fmt", ", arr[i].re, arr[i].im);\
            else printf(fmt"-j"fmt", ", arr[i].re, -arr[i].im);\
        }\
        printf("...,");\
        for (i=N-VEC_PRINT_LIMIT_COMPLEX/2; i<N-1; i++) {\
            if (arr[i].im >= 0) printf(fmt"+j"fmt", ", arr[i].re, arr[i].im);\
            else printf(fmt"-j"fmt", ", arr[i].re, -arr[i].im);\
        }\
    }\
    if (arr[N-1].im >= 0) printf(fmt"+j"fmt"]", arr[N-1].re, arr[N-1].im);\
    else printf(fmt"-j"fmt"]", arr[N-1].re, -arr[N-1].im);\
    if (sufix) printf("%s", sufix);\
}
GEN__carr_print(complex64, "%e");
GEN__carr_print(complex128, "%e");

void vec_print_xfix(const vector* v, const char* prefix, const char* sufix) {
    check(v->data, "Uninitialised vector data.");
    check(IS_VEC_VALID_DTYPE(*v), "Invalid data type for vec_print_xfix.");
    SIMPLE_SWITCH_NO_RETURN(arr_print, v->dtype, v->data, v->size, prefix, sufix)
}

void vec_print(const vector* v) { vec_print_xfix(v, "", "\n"); }

void mat_print(matrix* m) {
    printf("[");
    vector_view vw;
    int i;
    mat_get_row(m, 0, &vw);
    vec_print_xfix(&vw, NULL, "\n");
    for (i=1; i<m->row-1; i++) {
        mat_get_row(m, i, &vw);
        vec_print_xfix(&vw, " ", "\n");
    }
    mat_get_row(m, i, &vw);
    vec_print_xfix(&vw, " ", NULL);
    printf("]\n");
}

vector* vec_alloc(int size, Dtype dtype) {
    if (size <= 0) return NULL;
    vector* v = (vector*)malloc(sizeof(vector));
    check(v, "vector alloc failed.");
    int size_bytes = size;
    switch (dtype) {
        #define TMP_CAL_SIZE(T) size_bytes *= sizeof(T); v->elemsize = sizeof(T)
        case Int8:   TMP_CAL_SIZE(char); break;
        case Int16:   TMP_CAL_SIZE(short); break;
        case Int32:   TMP_CAL_SIZE(int); break;
        case Int64:   TMP_CAL_SIZE(long); break;
        case Float32: TMP_CAL_SIZE(float); break;
        case Float64: TMP_CAL_SIZE(double); break;
        case Complex64:  TMP_CAL_SIZE(complex64); break;
        case Complex128: TMP_CAL_SIZE(complex128); break;
        #undef TMP_CAL_SIZE
        default: check(0, "Invalid data type for vec_alloc.");
    }
    v->size = size;
    v->dtype = dtype;
    v->stride = 1;
    v->byte_stride = v->elemsize * v->stride;
    v->data = malloc(size_bytes);
    check(v->data, "vector data alloc failed.");
    return v;
}

void vec_free(vector* v) {
    if (v) {
        if (v->data) { free(v->data); v->data = NULL; }
        free(v);
    }
}

void batch_vec_free(vector* v, ...) {
    va_list ap;
    va_start(ap, v);
    vector* tmp;
    while ((tmp=va_arg(ap, vector*))!=NULL) {
        vec_free(tmp);
    }
    vec_free(v);
    va_end(ap);
}

#define GEN__vec_fill(T)\
void vec_fill_##T(vector* x, T a) {\
    T* x_data = (T*)(x->data);\
    for (int i = 0; i < x->size; i++) {T##_assign(*x_data, a); x_data += x->stride;}\
}
#define GEN__vec_fill_complex_parts(ComplexType, RealType)\
void vec_fill_##ComplexType##_parts(vector* x, RealType re, RealType im) {\
    ComplexType* x_data = (ComplexType*)(x->data);  \
    ComplexType tmp = {re, im};\
    for (int i = 0; i < x->size; i++) {ComplexType##_assign(*x_data, tmp); x_data += x->stride;}\
}
GEN__vec_fill(char)
GEN__vec_fill(short)
GEN__vec_fill(int)
GEN__vec_fill(long)
GEN__vec_fill(float)
GEN__vec_fill(double)
GEN__vec_fill(complex64)
GEN__vec_fill(complex128)
GEN__vec_fill_complex_parts(complex64, float)
GEN__vec_fill_complex_parts(complex128, double)

void vec_fill(vector* x, const void* a) {
    switch (x->dtype) {
#define TMP_FILL(T) vec_fill_##T(x, *(T*)a)
        case Int8:   TMP_FILL(char); break;
        case Int16:   TMP_FILL(short); break;
        case Int32:   TMP_FILL(int); break;
        case Int64:   TMP_FILL(long); break;
        case Float32: TMP_FILL(float); break;
        case Float64: TMP_FILL(double); break;
        case Complex64:  TMP_FILL(complex64); break;
        case Complex128: TMP_FILL(complex128); break;
#undef TMP_FILL
        default: check(0, "Invalid data type for vec_fill_real.");
    }
}

void vec_copy_slice(const vector* x, int xoffset, int size, vector* y, int yoffset) {
    check((x->dtype==y->dtype)&&(x->size>=xoffset+size)&&(y->size>=yoffset+size),
          "Inconsistent data-type or size.");
    int xoffset_bytes = xoffset, yoffset_bytes = yoffset, size_bytes = size;
    switch (x->dtype) {
#define TMP_CAL_SIZE(T) {xoffset_bytes*=sizeof(T);yoffset_bytes*=sizeof(T);size_bytes*=sizeof(T);}
        case Int32:   TMP_CAL_SIZE(int); break;
        case Int64:   TMP_CAL_SIZE(long); break;
        case Float32: TMP_CAL_SIZE(float); break;
        case Float64: TMP_CAL_SIZE(double); break;
        case Complex64:  TMP_CAL_SIZE(complex64); break;
        case Complex128: TMP_CAL_SIZE(complex128); break;
#undef TMP_CAL_SIZE
        default: check(0, "Invalid data type for vec_copy_slice.");
    }
    char* x_data = (char*)(x->data) + xoffset_bytes;
    char* y_data = (char*)(y->data) + yoffset_bytes;
    if (x->stride == 1) { // for contiguous memory
        memcpy(y_data, x_data, size_bytes);
    } else {  // for striding data
        for (int i=0; i<size; i++) {
            for (int j=0; j<x->elemsize; j++) *y_data = *x_data;
            y_data += y->byte_stride;
            x_data += y->byte_stride;
        }
    }
}

void vec_copy(const vector* src, vector* dst) {
    vec_copy_slice(src, 0, src->size, dst, 0);
}

vector* zeros(int size, Dtype dtype) {
    if (size <= 0) return NULL;
    vector* v = vec_alloc(size, dtype);
    switch (dtype) {
        case Int8:   {char a=char_zero; vec_fill(v, &a); break;}
        case Int16:   {short a=short_zero; vec_fill(v, &a); break;}
        case Int32:   {int a=int_zero; vec_fill(v, &a); break;}
        case Int64:   {long a=long_zero; vec_fill(v, &a); break;}
        case Float32: {float a=float_zero; vec_fill(v, &a); break;}
        case Float64: {double a=double_zero; vec_fill(v, &a); break;}
        case Complex64:  {complex64 a=complex64_zero; vec_fill(v, &a); break;}
        case Complex128: {complex128 a=complex128_zero; vec_fill(v, &a); break;}
        default:;
    }
    return v;
}

vector* zeros_like(const vector* x) { return zeros(x->size, x->dtype); }

vector* ones(int size, Dtype dtype) {
    if (size <= 0) return NULL;
    vector* v = vec_alloc(size, dtype);
    switch (dtype) {
        case Int8:   {char a=char_one; vec_fill(v, &a); break;}
        case Int16:   {short a=short_one; vec_fill(v, &a); break;}
        case Int32:   {int a=int_one; vec_fill(v, &a); break;}
        case Int64:   {long a=long_one; vec_fill(v, &a); break;}
        case Float32: {float a=float_one; vec_fill(v, &a); break;}
        case Float64: {double a=double_one; vec_fill(v, &a); break;}
        case Complex64:  {complex64 a=complex64_one; vec_fill(v, &a); break;}
        case Complex128: {complex128 a=complex128_one; vec_fill(v, &a); break;}
        default: check(0, "Invalid data type for ones.");
    }
    return v;
}

vector* ones_like(const vector* x) { return ones(x->size, x->dtype); }

vector* arange(int end, Dtype dtype) {
    vector* v = vec_alloc(end, dtype);
    switch (dtype) {
#define ARANGE(v, T) for (int i=0; i<(v)->size; i++) vec_elem_set_##T((v), i, (T)i)
        case Int16:   ARANGE(v, short); break;
        case Int32:   ARANGE(v, int); break;
        case Int64:   ARANGE(v, long); break;
        case Float32: ARANGE(v, float); break;
        case Float64: ARANGE(v, double); break;
#undef ARANGE
        default: check(0, "Invalid data type for arange.");
    }
    return v;
}

vector* arange2(int start, int end, int step, Dtype dtype) {
    check(0, "Not implemented.");
    vector* v = vec_alloc((end-start)/step, dtype);
    // TODO
    return v;
}

void vec_view_construct(const vector* v, int offset, int size, vector_view* vw) {
    check(v->size>=offset+size, "Invalid slice bounds for vec_view_construct.");
    vw->dtype = v->dtype;
    vw->size = size;
    vw->elemsize = v->elemsize;
    vw->stride = v->stride;
    vw->byte_stride = v->byte_stride;
    SWITCH_DATA_OFFSET(v, offset, vw)
}

void vec_view_destruct(vector_view* v) { v->data = NULL; /* view won't own any data. */ }

void vec_elem_set(vector* x, int i, const void* a) {
    check(i<x->size, "Invalid index for vec_set.");
    switch (x->dtype) {
#define TMP_ELEM_SET(T) vec_elem_set_##T(x, i, *(T*)a)
        case Int8:  TMP_ELEM_SET(char); break;
        case Int16: TMP_ELEM_SET(short); break;
        case Int32: TMP_ELEM_SET(int); break;
        case Int64: TMP_ELEM_SET(long); break;
        case Float32: TMP_ELEM_SET(float); break;
        case Float64: TMP_ELEM_SET(double); break;
        case Complex64:  TMP_ELEM_SET(complex64); break;
        case Complex128: TMP_ELEM_SET(complex128); break;
#undef TMP_ELEM_SET
        default:;
    }
}

void vec_elem_incr(vector* x, int i, const void* a) {
    check(i<x->size, "Invalid index for vec_set.");
    switch (x->dtype) {
#define TMP_ELEM_INCR(T) vec_elem_incr_##T(x, i, *(T*)a)
        case Int16:   TMP_ELEM_INCR(short); break;
        case Int32:   TMP_ELEM_INCR(int); break;
        case Int64:   TMP_ELEM_INCR(long); break;
        case Float32: TMP_ELEM_INCR(float); break;
        case Float64: TMP_ELEM_INCR(double); break;
        case Complex64:  TMP_ELEM_INCR(complex64); break;
        case Complex128: TMP_ELEM_INCR(complex128); break;
#undef TMP_ELEM_INCR
    }
}

void vec_elem_get(const vector* x, int i, void* a) {
    check(x->data && i<x->size && i>=0, "Invalid index or uninitialised data for vec_get.");
    switch (x->dtype) {
#define TMP_ELEM_GET(T) *(T*)a = vec_elem_get_##T(x, i)
        case Int16:   TMP_ELEM_GET(short); break;
        case Int32:   TMP_ELEM_GET(int); break;
        case Int64:   TMP_ELEM_GET(long); break;
        case Float32: TMP_ELEM_GET(float); break;
        case Float64: TMP_ELEM_GET(double); break;
        case Complex64:  TMP_ELEM_GET(complex64); break;
        case Complex128: TMP_ELEM_GET(complex128); break;
#undef TMP_ELEM_GET
        default:;
    }
}

void mat_view_construct(const matrix* v, matrix_view* vw) { /* TODO */ }

void mat_view_from_vec_slice(const vector* v, int offset, int size, int col, matrix_view* vw) {
    check((size%col)==0, "Invalid view shape.");
    vw->dtype = v->dtype;
    vw->col = col;
    vw->row = size/col;
    vw->stride = v->stride;
    vw->col_stride = v->stride * vw->col;
    vw->byte_stride = v->byte_stride;
    vw->elemsize = v->elemsize;
    SWITCH_DATA_OFFSET(v, offset, vw);
}

void mat_view_from_vec(const vector* v, int col, matrix_view* vw) {
    mat_view_from_vec_slice(v, 0, v->size, col, vw);
}

void mat_view_destruct(matrix_view* v) { v->data = NULL; /* view won't own any data. */ }

void mat_flatten(const matrix* m, vector_view* v) {
    v->dtype = m->dtype;
    v->size = m->row * m->col;
    v->data = m->data;
    v->byte_stride = m->elemsize * m->stride;
    v->stride = m->stride;
    v->elemsize = m->elemsize;
}

void vec_real_part(const vector* x, vector_view* y) {
 y->dtype = (x->dtype==Complex64) ? Float32 : Float64;
 y->size = x->size;
 y->elemsize = x->elemsize/2;
 y->byte_stride = x->byte_stride;
 y->stride = x->stride*2;
 y->data = x->data;
}

void vec_imag_part(const vector* x, vector_view* y) {
    y->dtype = (x->dtype==Complex64) ? Float32 : Float64;
    y->size = x->size;
    y->elemsize = x->elemsize/2;
    y->byte_stride = x->byte_stride;
    y->stride = x->stride*2;
    y->data = (char*)(x->data) + y->elemsize;
}

void vec_abs2(const vector* x, vector* y) {
    check(x->size==y->size, "Size of x and y should be equal.");
    switch (x->dtype) {
        case Int16:
        case Int32:
        case Int64:
        case Float32:
        case Float64:
            vec_mul(x, x, y);
            break;
        case Complex64:
        case Complex128: {
            vector_view real_part, image_part;
            vec_real_part(x, &real_part);
            vec_imag_part(x, &image_part);
            vec_mul(&real_part, &real_part, y);
            vec_mul_accum(&image_part, &image_part, y);
            break;
        }
        default: check(0, "Invalid data type for vec_abs2.");
    }
}

//void vec_abs(const vector* x, vector* y) {
//    vec_abs2(x, y);
//    switch (x->dtype) {
//        case Float32: return vec_unary_float(y, y, sqrt_float);
//        case Float64: return vec_unary_double(y, y, sqrt_double);
//        default: check(0, "Unimplemented data type for vec_abs.");
//    }
//}
void vec_abs(const vector* x, vector* y) {
    vec_abs2(x, y);
    vec_sqrt(y, y);
}

void vec_norm2(const vector* x, void* a) {
//    return vec_unary_accum_float(x, ssdsp_square, ssdsp_add);
    switch (x->dtype) {
        case Complex64: {
            // 虽然可以使用vec_dot实现norm2，但是为了节约一半的计算，在这里单独实现
            const complex64* x_data = (complex64*)(x->data);
            float norm2 = 0.f;
            for (int i=0; i<x->size; i++) norm2 += complex64_abs2(x_data[i]);
            *(float*)a = norm2;
            break;
        }
        case Complex128: {
            const complex128* x_data = (complex128*)(x->data);
            double norm2 = 0.;
            for (int i=0; i<x->size; i++) norm2 += complex128_abs2(x_data[i]);
            *(double*)a = norm2;
            break;
        }
        case Int32:
        case Int64:
        case Float32:
        case Float64:
            return vec_dot(x, x, a);
        default:;
    }
}

void vec_norm(const vector* x, void* a) {
    vec_norm2(x, a);
    switch (x->dtype) {
        case Int32:
        case Float32: {
            float norm = sqrtf(*(float *) a);
            *(float*)a = norm;
            break;
        }
        case Int64:
        case Float64: {
            double norm = sqrt(*(double *) a);
            *(double*)a = norm;
            break;
        }
        default: check(0, "Invalid data type for vec_norm.");
    }
}

void mat_get_col(const matrix* m, int c, vector_view* vw) {
    vw->dtype = m->dtype;
    vw->size = m->row;
    vw->stride = m->col_stride;
    vw->byte_stride = m->col * m->elemsize * m->stride;
    int offset = c * m->stride;
    SWITCH_DATA_OFFSET(m, offset, vw)
}

void mat_get_row(const matrix* m, int r, vector_view* vw) {
    vw->dtype = m->dtype;
    vw->size = m->col;
    vw->stride = m->stride;
    vw->byte_stride = m->elemsize * m->stride;
    vw->elemsize = m->elemsize;
    int offset = r * m->col_stride;
    SWITCH_DATA_OFFSET(m, offset, vw)
}

// matrix has row-major data layout
#define GEN__mat_elemwise_op(op)\
void mat_elemwise_##op(const matrix* x, const matrix* y, matrix* z) { \
    int shape_cond_1 = ((x->row==y->row)&&((x->col%y->col==0)||(y->col%x->col==0)));\
    int shape_cond_2 = ((x->col==y->col)&&((x->row==1)||(y->row==1)));\
    check(shape_cond_1 || shape_cond_2, "Inconsistent matrix shape for elemwise op.");\
    if (shape_cond_1) {        \
        for (int i=0; i<x->row; i++) {\
            vector_view vw_x_row, vw_y_row, vw_z_row;\
            mat_get_row(x, i, &vw_x_row);\
            mat_get_row(y, i, &vw_y_row);\
            mat_get_row(z, i, &vw_z_row);\
            vec_##op(&vw_x_row, &vw_y_row, &vw_z_row);\
        }\
    } else if (shape_cond_2) {\
        vector_view x_flat, y_flat, z_flat;\
        mat_flatten(x, &x_flat);\
        mat_flatten(y, &y_flat);\
        mat_flatten(z, &z_flat);\
        vec_##op(&x_flat, &y_flat, &z_flat);\
    }\
}
GEN__mat_elemwise_op(add)
GEN__mat_elemwise_op(sub)
GEN__mat_elemwise_op(mul)
GEN__mat_elemwise_op(div)
