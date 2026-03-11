#ifndef BLAS_REF__H
#define BLAS_REF__H
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#define blas_ref_axpy(n, a, x, y)                                              \
  _Generic(*(x),                                                               \
      float: blas_ref_axpyf,                                                   \
      double: blas_ref_axpy,                                                   \
      long double: blas_ref_axpyl)(n, a, x, y)

inline float *blas_ref_axpyf(size_t n, float alpha, float const *restrict x,
                             float *restrict y) {
  for (size_t i = 0; i < n; ++i) {
    y[i] += alpha * x[i];
  }
  return y;
}
inline double *(blas_ref_axpy)(size_t n, double alpha, double const *restrict x,
                               double *restrict y) {
  for (size_t i = 0; i < n; ++i) {
    y[i] += alpha * x[i];
  }
  return y;
}
inline long double *blas_ref_axpyl(size_t n, long double alpha,
                                   long double const *restrict x,
                                   long double *restrict y) {
  for (size_t i = 0; i < n; ++i) {
    y[i] += alpha * x[i];
  }
  return y;
}

#define blas_ref_scal(n, a, x)                                                 \
  _Generic(*(x),                                                               \
      float: blas_ref_scalf,                                                   \
      double: blas_ref_scal,                                                   \
      long double: blas_ref_scall)(n, a, x)

inline float *blas_ref_scalf(size_t n, float alpha, float *x) {
  for (size_t i = 0; i < n; ++i) {
    x[i] *= alpha;
  }
  return x;
}
inline double *(blas_ref_scal)(size_t n, double alpha, double *x) {
  for (size_t i = 0; i < n; ++i) {
    x[i] *= alpha;
  }
  return x;
}
inline long double *blas_ref_scall(size_t n, long double alpha,
                                   long double *x) {
  for (size_t i = 0; i < n; ++i) {
    x[i] *= alpha;
  }
  return x;
}

#define blas_ref_copy(n, x, y)                                                 \
  _Generic(*(x),                                                               \
      float: blas_ref_copyf,                                                   \
      double: blas_ref_copy,                                                   \
      long double: blas_ref_copyl)(n, x, y)

inline void blas_ref_copyf(size_t n, float const *restrict x,
                           float *restrict y) {
  memcpy(y, x, n * sizeof *x);
}

inline void(blas_ref_copy)(size_t n, double const *restrict x,
                           double *restrict y) {
  memcpy(y, x, n * sizeof *x);
}

inline void blas_ref_copyl(size_t n, long double const *restrict x,
                           long double *restrict y) {
  memcpy(y, x, n * sizeof *x);
}

#define blas_ref_swap(n, x, y)                                                 \
  _Generic(*(x),                                                               \
      float: blas_ref_swapf,                                                   \
      double: blas_ref_swap,                                                   \
      long double: blas_ref_swapl)(n, x, y) {

inline void blas_ref_swapf(size_t n, float *restrict x, float *restrict y) {
  float tmp;
  for (size_t i = 0; i < n; ++i) {
    tmp = x[i];
    x[i] = y[i];
    y[i] = tmp;
  }
}

inline void(blas_ref_swap)(size_t n, double *restrict x, double *restrict y) {
  double tmp;
  for (size_t i = 0; i < n; ++i) {
    tmp = x[i];
    x[i] = y[i];
    y[i] = tmp;
  }
}

inline void blas_ref_swapl(size_t n, long double *restrict x,
                           long double *restrict y) {
  long double tmp;
  for (size_t i = 0; i < n; ++i) {
    tmp = x[i];
    x[i] = y[i];
    y[i] = tmp;
  }
}

#define blas_ref_dot(n, x, y)                                                  \
  _Generic(*(x),                                                      \
      float: blas_ref_dotf,                                                        \
      double: blas_ref_dot,                                                        \
      long double: blas_ref_dotl)(n, x, y)

inline float blas_ref_dotf(size_t n, float const *x, float const *y) {
  float res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += x[i] * y[i];
  }
  return res;
}

inline double(blas_ref_dot)(size_t n, double const *x, double const *y) {
  double res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += x[i] * y[i];
  }
  return res;
}

inline long double blas_ref_dotl(size_t n, long double const *x,
                                 long double const *y) {
  long double res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += x[i] * y[i];
  }
  return res;
}

#define blas_ref_nrm2(n, x)                                                    \
  _Generic(*(x),                                                               \
      float: blas_ref_nrm2f,                                                   \
      double: blas_ref_nrm2,                                                   \
      long double: blas_ref_nrm2l)(n, x)

inline float blas_ref_nrm2f(size_t n, float const *x) {
  float res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += x[i] * x[i];
  }
  return sqrt(res);
}
inline double(blas_ref_nrm2)(size_t n, double const *x) {
  double res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += x[i] * x[i];
  }
  return sqrt(res);
}
inline long double blas_ref_nrm2l(size_t n, long double const *x) {
  long double res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += x[i] * x[i];
  }
  return sqrt(res);
}

#define blas_ref_asum(n, x)                                                    \
  _Generic(*(x),                                                               \
      float: blas_ref_asumf,                                                   \
      double: blas_ref_asum,                                                   \
      long double: blas_ref_asuml)(n, x)

inline float blas_ref_asumf(size_t n, float const *x) {
  float res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += fabs(x[i]);
  }
  return res;
}

inline double(blas_ref_asum)(size_t n, double const *x) {
  double res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += fabs(x[i]);
  }
  return res;
}

inline long double blas_ref_asuml(size_t n, long double const *x) {
  long double res = 0;
  for (size_t i = 0; i < n; ++i) {
    res += fabs(x[i]);
  }
  return res;
}

#define blas_ref_empty(n, x)                                                   \
  _Generic(*(x),                                                               \
      float *: blas_ref_emptyf,                                                \
      double *: blas_ref_empty,                                                \
      long double *: blas_ref_emptyl)(n, x)

inline float *blas_ref_emptyf(size_t n, float **x) {
  void *tmp = realloc(*x, n * sizeof(**x));
  if (!tmp)
    return 0;
  return *x = tmp;
}
inline double *(blas_ref_empty)(size_t n, double **x) {
  void *tmp = realloc(*x, n * sizeof(**x));
  if (!tmp)
    return 0;
  return *x = tmp;
}

inline long double *blas_ref_emptyl(size_t n, long double **x) {
  void *tmp = realloc(*x, n * sizeof(**x));
  if (!tmp)
    return 0;
  return *x = tmp;
}

#define blas_ref_zeros(n, x)                                                   \
  _Generic(*(x),                                                               \
      float *: blas_ref_zerosf,                                                \
      double *: blas_ref_zeros,                                                \
      long double *: blas_ref_zerosl)(n, x)

inline float *blas_ref_zerosf(size_t n, float **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = 0;
  }
  return *x;
}

inline double *(blas_ref_zeros)(size_t n, double **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = 0;
  }
  return *x;
}

inline long double *blas_ref_zerosl(size_t n, long double **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = 0;
  }
  return *x;
}

#define blas_ref_ones(n, x)                                                    \
  _Generic(*(x),                                                               \
      float *: blas_ref_onesf,                                                 \
      double *: blas_ref_ones,                                                 \
      long double *: blas_ref_onesl)(n, x)

inline float *blas_ref_onesf(size_t n, float **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = 1;
  }
  return *x;
}

inline double *(blas_ref_ones)(size_t n, double **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = 1;
  }
  return *x;
}

inline long double *blas_ref_onesl(size_t n, long double **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = 1;
  }
  return *x;
}

#define blas_ref_rand(n, x)                                                    \
  _Generic(*(x),                                                               \
      float *: blas_ref_randf,                                                   \
      double *: blas_ref_rand,                                                   \
      long double *: blas_ref_randl)(n, x)

inline float *blas_ref_randf(size_t n, float **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = drand48() * 2 - 1;
  }
  return *x;
}

inline double *(blas_ref_rand)(size_t n, double **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = drand48() * 2 - 1;
  }
  return *x;
}

inline long double *blas_ref_randl(size_t n, long double **x) {
  blas_ref_empty(n, x);
  for (size_t i = 0; i < n; ++i) {
    (*x)[i] = drand48() * 2 - 1;
  }
  return *x;
}

inline void blas_ref_free(void *p) { free(p); }

#define blas_ref_mse(n, x, y)                                                  \
  _Generic(*(x),                                                               \
      float: blas_ref_msef,                                                    \
      double: blas_ref_mse,                                                    \
      long double: blas_ref_msel)(n, x, y)

inline float blas_ref_msef(size_t n, float const *x, float const *y) {
  float err = 0;
  for (size_t i = 0; i < n; ++i) {
    err += pow(x[i] - y[i], 2);
  }
  return err / n;
}

inline double(blas_ref_mse)(size_t n, double const *x, double const *y) {
  double err = 0;
  for (size_t i = 0; i < n; ++i) {
    err += pow(x[i] - y[i], 2);
  }
  return err / n;
}

inline long double blas_ref_msel(size_t n, long double const *x,
                                 long double const *y) {
  float err = 0;
  for (size_t i = 0; i < n; ++i) {
    err += pow(x[i] - y[i], 2);
  }
  return err / n;
}

#endif
