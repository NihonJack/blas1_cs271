#define _POSIX_C_SOURCE 200809L
#define _DEFAULT_SOURCE

#include <stdlib.h>

#include "blas_ref.h"

extern inline float *blas_ref_axpyf(size_t n, float alpha,
                                    float const *restrict x, float *restrict y);
extern inline double *(blas_ref_axpy)(size_t n, double alpha,
                                      double const *restrict x,
                                      double *restrict y);
extern inline long double *blas_ref_axpyl(size_t n, long double alpha,
                                          long double const *restrict x,
                                          long double *restrict y);

extern inline float *blas_ref_scalf(size_t n, float alpha, float *x);
extern inline double *(blas_ref_scal)(size_t n, double alpha, double *x);
extern inline long double *blas_ref_scall(size_t n, long double alpha,
                                          long double *x);

extern inline void blas_ref_copyf(size_t n, float const *restrict x,
                                  float *restrict y);
extern inline void(blas_ref_copy)(size_t n, double const *restrict x,
                                  double *restrict y);
extern inline void blas_ref_copyl(size_t n, long double const *restrict x,
                                  long double *restrict y);

extern inline void blas_ref_swapf(size_t n, float *restrict x,
                                  float *restrict y);
extern inline void(blas_ref_swap)(size_t n, double *restrict x,
                                  double *restrict y);
extern inline void blas_ref_swapl(size_t n, long double *restrict x,
                                  long double *restrict y);

extern inline float blas_ref_dotf(size_t n, float const *x, float const *y);
extern inline double(blas_ref_dot)(size_t n, double const *x, double const *y);
extern inline long double blas_ref_dotl(size_t n, long double const *x,
                                        long double const *y);

extern inline float blas_ref_nrm2f(size_t n, float const *x);
extern inline double(blas_ref_nrm2)(size_t n, double const *x);
extern inline long double blas_ref_nrm2l(size_t n, long double const *x);

extern inline float blas_ref_asumf(size_t n, float const *x);
extern inline double(blas_ref_asum)(size_t n, double const *x);
extern inline long double blas_ref_asuml(size_t n, long double const *x);

extern inline float *blas_ref_emptyf(size_t n, float **x);
extern inline double *(blas_ref_empty)(size_t n, double **x);
extern inline long double *blas_ref_emptyl(size_t n, long double **x);
extern inline float *blas_ref_zerosf(size_t n, float **x);
extern inline double *(blas_ref_zeros)(size_t n, double **x);
extern inline long double *blas_ref_zerosl(size_t n, long double **x);
extern inline float *blas_ref_onesf(size_t n, float **x);
extern inline double *(blas_ref_ones)(size_t n, double **x);
extern inline long double *blas_ref_onesl(size_t n, long double **x);
extern inline float *blas_ref_randf(size_t n, float **x);
extern inline double *(blas_ref_rand)(size_t n, double **x);
extern inline long double *blas_ref_randl(size_t n, long double **x);
extern inline void blas_ref_free(void *p);
extern inline float blas_ref_msef(size_t n, float const *x, float const *y);
extern inline double(blas_ref_mse)(size_t n, double const *x, double const *y);
extern inline long double blas_ref_msel(size_t n, long double const *x,
                                        long double const *y);
