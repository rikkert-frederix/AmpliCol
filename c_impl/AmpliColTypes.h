#pragma once

// AmpliColTypes.h
// Basic types for AmpliCol C implementation.
// Dual C/C++/CUDA header: include from any language.
//   C    : AC_CX = double complex        (C99 <complex.h>)
//   C++  : AC_CX = std::complex<double>  (same ABI, C++11 §26.4)
//   CUDA : AC_CX = cuda::std::complex<double> (same ABI; operators are
//          __host__ __device__, unlike std::complex which is __host__ only)

#ifdef __CUDACC__
#  include <cuda/std/complex>
#  include <cmath>
#  include <cstddef>
#  ifdef AC_FP_SINGLE
     typedef float                        AC_FP;
     typedef cuda::std::complex<float>    AC_CX;
#  else
     typedef double                       AC_FP;
     typedef cuda::std::complex<double>   AC_CX;
#  endif
     typedef double                       AC_D_FP;
     typedef cuda::std::complex<double>   AC_D_CX;
#  define AC_HOSTDEV __host__ __device__
#elif defined(__cplusplus)
#  include <complex>
#  include <cmath>
#  include <cstddef>
#  ifdef AC_FP_SINGLE
     typedef float                        AC_FP;
     typedef std::complex<float>          AC_CX;
#  else
     typedef double                       AC_FP;
     typedef std::complex<double>         AC_CX;
#  endif
     typedef double                       AC_D_FP;
     typedef std::complex<double>         AC_D_CX;
#  define AC_HOSTDEV
#else
#  include <complex.h>
#  include <math.h>
#  include <stddef.h>
#  ifdef AC_FP_SINGLE
     typedef float                        AC_FP;
     typedef float complex                AC_CX;
#  else
     typedef double                       AC_FP;
     typedef double complex               AC_CX;
#  endif
   typedef double                         AC_D_FP;
   typedef double complex                 AC_D_CX;
#  define AC_HOSTDEV
#endif

// ─── inner products ────────────────────────────────────────────────────────

// First 4 components only (wavefunctions are size-6 to accommodate tensor gluons)
static AC_HOSTDEV inline AC_CX inprod_wf(const AC_CX* a, const AC_CX* b) {
    AC_CX sum = 0;
    for (int i = 0; i < 4; ++i)
        sum += a[i] * b[i];
    return sum;
}

static AC_HOSTDEV inline AC_CX inprod_cx(const AC_CX* a, const AC_CX* b, int n) {
    AC_CX sum = 0;
    for (int i = 0; i < n; ++i)
        sum += a[i] * b[i];
    return sum;
}

static AC_HOSTDEV inline AC_FP inprod_fp(const AC_FP* a, const AC_FP* b, int n) {
    AC_FP sum = 0;
    for (int i = 0; i < n; ++i)
        sum += a[i] * b[i];
    return sum;
}

// ─── ordered_sum ──────────────────────────────────────────────────────────
//
// Sum rows of a 2D array selected by index[min..max].
// Two conventions preserved from the C++ templates:
//   _excl : range is [min, max)  (matches the generic template)
//   _incl : range is [min, max]  (matches the AC_CX/6 specialisation)

// AC_CX rows of width 6, exclusive upper bound
static AC_HOSTDEV inline void ordered_sum_cx6_excl(
    const AC_CX src[][6], const int* index, int min, int max, AC_CX dst[6])
{
    for (int i = 0; i < 6; ++i) dst[i] = 0;
    for (int j = min; j < max; ++j)
        for (int i = 0; i < 6; ++i)
            dst[i] += src[index[j]][i];
}

// AC_CX rows of width 6, inclusive upper bound
static AC_HOSTDEV inline void ordered_sum_cx6(
    const AC_CX src[][6], const int* index, int min, int max, AC_CX dst[6])
{
    for (int i = 0; i < 6; ++i) dst[i] = 0;
    for (int j = min; j <= max; ++j)
        for (int i = 0; i < 6; ++i)
            dst[i] += src[index[j]][i];
}

// AC_FP rows of width N, exclusive upper bound (generic via runtime width)
static AC_HOSTDEV inline void ordered_sum_fp_excl(
    const AC_FP* src, int row_width,
    const int* index, int min, int max, AC_FP* dst)
{
    for (int i = 0; i < row_width; ++i) dst[i] = 0;
    for (int j = min; j < max; ++j)
        for (int i = 0; i < row_width; ++i)
            dst[i] += src[index[j] * row_width + i];
}

struct AC_AMP{ void (*eval)(const AC_D_FP[][4], AC_D_CX*); };

// ─── reshape helpers ───────────────────────────────────────────────────────
//
// reshape_colmajor: copy flat[N1*N2] → dst[N1][N2]  (element order preserved)
// In C the pointer-cast equivalent is:  AC_CX (*mat)[N2] = (AC_CX(*)[N2])flat;
// Use that pattern directly for zero-copy aliasing.

static AC_HOSTDEV inline void reshape_colmajor_int(
    const int* src, size_t n1, size_t n2, int* dst)
{
    for (size_t i = 0; i < n1 * n2; ++i)
        dst[(i / n2) * n2 + (i % n2)] = src[i];
}

static AC_HOSTDEV inline void reshape_colmajor_cx(
    const AC_CX* src, size_t n1, size_t n2, AC_CX* dst)
{
    for (size_t i = 0; i < n1 * n2; ++i)
        dst[(i / n2) * n2 + (i % n2)] = src[i];
}

static AC_HOSTDEV inline void reshape_colmajor_fp(
    const AC_FP* src, size_t n1, size_t n2, AC_FP* dst)
{
    for (size_t i = 0; i < n1 * n2; ++i)
        dst[(i / n2) * n2 + (i % n2)] = src[i];
}
