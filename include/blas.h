#ifndef BLAS__H
#define BLAS__H
#include <stddef.h>
#include <string.h>

#define blas_axpy(n, a, x, y)                                                  \
  _Generic(*(x),                                                \
      float: blas_axpyf,                                                       \
      double: blas_axpy,                                                       \
      long double: blas_axpyl)(n, a, x, y)

float *blas_axpyf(size_t n, float alpha, float const *restrict x,
                  float *restrict y);
double *(blas_axpy)(size_t n, double alpha, double const *restrict x,
                    double *restrict y);
long double *blas_axpyl(size_t n, long double alpha,
                        long double const *restrict x, long double *restrict y);

#define blas_scal(n, a, x)                                                     \
  _Generic(*(x),                                                       \
      float: blas_scalf,                                                       \
      double: blas_scal,                                                       \
      long double: blas_scall)(n, a, x)

float *blas_scalf(size_t n, float alpha, float *x);
double *(blas_scal)(size_t n, double alpha, double *x);
long double *blas_scall(size_t n, long double alpha, long double *x);

#define blas_copy(n, x, y)                                                     \
  _Generic(*(x),                                                      \
      float: blas_copyf,                                                       \
      double: blas_copy,                                                       \
      long double: blas_copyl)(n, x, y)

inline void blas_copyf(size_t n, float const *restrict x, float *restrict y) {
  memcpy(y, x, n * sizeof *x);
}

inline void(blas_copy)(size_t n, double const *restrict x, double *restrict y) {
  memcpy(y, x, n * sizeof *x);
}

inline void blas_copyl(size_t n, long double const *restrict x,
                long double *restrict y) {
  memcpy(y, x, n * sizeof *x);
}

#define blas_swap(n, x, y)                                                     \
  _Generic(*(x),                                                        \
      float: blas_swapf,                                                       \
      double: blas_swap,                                                       \
      long double: blas_swapl)(n, x, y);

void blas_swapf(size_t n, float *restrict x, float *restrict y);
void(blas_swap)(size_t n, double *restrict x, double *restrict y);
void blas_swapl(size_t n, long double *restrict x, long double *restrict y);

#define blas_dot(n, x, y)                                                      \
  _Generic(*(x),                                                      \
      float: blas_dotf,                                                        \
      double: blas_dot,                                                        \
      long double: blas_dotl)(n, x, y)

float blas_dotf(size_t n, float const *x, float const *y);
double(blas_dot)(size_t n, double const *x, double const *y);
long double blas_dotl(size_t n, long double const *x, long double const *y);

#define blas_nrm2(n, x)                                                        \
  _Generic(*(x),                                                             \
      float: blas_nrm2f,                                                       \
      double: blas_nrm2,                                                       \
      long double: blas_nrm2l)(n, x)

float blas_nrm2f(size_t n, float const *x);
double(blas_nrm2)(size_t n, double const *x);
long double blas_nrm2l(size_t n, long double const *x);

#define blas_asum(n, x)                                                        \
  _Generic(*(x),                                                             \
      float: blas_asumf,                                                       \
      double: blas_asum,                                                       \
      long double: blas_asuml)(n, x)

float blas_asumf(size_t n, float const *x);
double(blas_asum)(size_t n, double const *x);
long double blas_asuml(size_t n, long double const *x);

#endif
