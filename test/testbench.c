#define _POSIX_C_SOURCE 200809L
#define _DEFAULT_SOURCE
#include <err.h>
#include <getopt.h>
#include <libgen.h>
#include <setjmp.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "blas.h"
#include "blas_ref.h"
#include "test_utils.h"

volatile sig_atomic_t last_signal = 0;
static sigjmp_buf env;
static void handler(int sig);
static void install_handlers(void);

#ifndef ep
#define ep 1.0e-5
#endif

static size_t test_vec_sizes[] = {
    1,  2,  3,  4,  5,  6,  7,  8,   9,    10,   11,   12,   13,  14,
    15, 16, 17, 18, 19, 20, 64, 100, 1024, 1025, 1026, 1027, 4096};

static int verbosity = 0;

#define X(name, T)                                                             \
  static int name(void) {                                                      \
    log_header(#name " (" #T ")");                                             \
                                                                               \
    T *x = 0;                                                                  \
    T *y_in = 0;                                                               \
    T *y = 0;                                                                  \
    T *yref = 0;                                                               \
    T *alpha = 0;                                                              \
    bool pass = true;                                                          \
    for (size_t i = 0; i < sizeof test_vec_sizes / sizeof *test_vec_sizes;     \
         ++i) {                                                                \
      size_t n = test_vec_sizes[i];                                            \
      blas_ref_rand(n, &x);                                                    \
      blas_ref_rand(1, &alpha);                                                \
      blas_ref_rand(n, &y_in);                                                 \
      blas_ref_empty(n, &y);                                                   \
      blas_ref_empty(n, &yref);                                                \
      blas_ref_copy(n, y_in, y);                                               \
      blas_ref_copy(n, y_in, yref);                                            \
      blas_axpy(n, alpha[0], x, y);                                            \
      blas_ref_axpy(n, alpha[0], x, yref);                                     \
      T mse = blas_ref_mse(n, y, yref);                                        \
      if (!(mse < ep)) {                                                       \
        pass = false;                                                          \
        if (verbosity > 0) {                                                   \
          log_info("N=%zu, MSE (%Lf) < ep (%Lf)", n, (long double)mse,         \
                   (long double)ep);                                           \
        }                                                                      \
        if (verbosity > 1) {                                                   \
          fputs("x = [", stderr);                                              \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)x[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fprintf(stderr, "alpha = %Lf\n", (long double)alpha[0]);             \
          fputs("y_in      = [", stderr);                                      \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)y_in[i],                     \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fputs("y      = [", stderr);                                         \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)y[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fputs("yref = [", stderr);                                           \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)yref[i],                     \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    blas_ref_free(x);                                                          \
    blas_ref_free(alpha);                                                      \
    blas_ref_free(y_in);                                                       \
    blas_ref_free(y);                                                          \
    blas_ref_free(yref);                                                       \
    if (pass) {                                                                \
      log_pass(0);                                                             \
      return 0;                                                                \
    } else {                                                                   \
      log_fail(0);                                                             \
      return -1;                                                               \
    }                                                                          \
  }
X(axpyf_test, float);
X(axpy_test, double);
X(axpyl_test, long double);
#undef X

#define X(name, T)                                                             \
  static int name(void) {                                                      \
    log_header(#name " (" #T ")");                                             \
                                                                               \
    T *x = 0;                                                                  \
    bool pass = true;                                                          \
    for (size_t i = 0; i < sizeof test_vec_sizes / sizeof *test_vec_sizes;     \
         ++i) {                                                                \
      size_t n = test_vec_sizes[i];                                            \
      blas_ref_rand(n, &x);                                                    \
      T res = blas_asum(n, x);                                                 \
      T ref_res = blas_ref_asum(n, x);                                         \
      T mse = blas_ref_mse(1, &res, &ref_res);                                 \
      if (!(mse < ep*n)) {                                                       \
        pass = false;                                                          \
        if (verbosity > 0) {                                                   \
          log_info("N=%zu, MSE (%Lf) < ep (%Lf)", n, (long double)mse,         \
                   (long double)ep);                                           \
        }                                                                      \
        if (verbosity > 1) {                                                   \
          fputs("x = [", stderr);                                              \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)x[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fprintf(stderr, "result = %Lf\n", (long double)res);                 \
          fprintf(stderr, "reference result  = %Lf\n", (long double)ref_res);  \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    blas_ref_free(x);                                                          \
    if (pass) {                                                                \
      log_pass(0);                                                             \
      return 0;                                                                \
    } else {                                                                   \
      log_fail(0);                                                             \
      return -1;                                                               \
    }                                                                          \
  }

X(asumf_test, float);
X(asum_test, double);
X(asuml_test, long double);
#undef X

#define X(name, T)                                                             \
  static int name(void) {                                                      \
    log_header(#name " (" #T ")");                                             \
                                                                               \
    T *x = 0;                                                                  \
    T *y_in = 0;                                                               \
    T *y = 0;                                                                  \
    T *yref = 0;                                                               \
    bool pass = true;                                                          \
    for (size_t i = 0; i < sizeof test_vec_sizes / sizeof *test_vec_sizes;     \
         ++i) {                                                                \
      size_t n = test_vec_sizes[i];                                            \
      blas_ref_rand(n, &x);                                                    \
      blas_ref_rand(n, &y);                                                    \
      T res = blas_dot(n, x, y);                                               \
      T ref_res = blas_ref_dot(n, x, y);                                       \
      T mse = blas_ref_mse(1, &res, &ref_res);                                 \
      if (!(mse < ep*n)) {                                                       \
        pass = false;                                                          \
        if (verbosity > 0) {                                                   \
          log_info("N=%zu, MSE (%Lf) < ep (%Lf)", n, (long double)mse,         \
                   (long double)ep);                                           \
        }                                                                      \
        if (verbosity > 1) {                                                   \
          fputs("x = [", stderr);                                              \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)x[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fputs("y = [", stderr);                                              \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)y[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fprintf(stderr, "result = %Lf\n", (long double)res);                 \
          fprintf(stderr, "reference result  = %Lf\n", (long double)ref_res);  \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    blas_ref_free(x);                                                          \
    blas_ref_free(y_in);                                                       \
    blas_ref_free(y);                                                          \
    blas_ref_free(yref);                                                       \
    if (pass) {                                                                \
      log_pass(0);                                                             \
      return 0;                                                                \
    } else {                                                                   \
      log_fail(0);                                                             \
      return -1;                                                               \
    }                                                                          \
  }
X(dotf_test, float);
X(dot_test, double);
X(dotl_test, long double);
#undef X

#define X(name, T)                                                             \
  static int name(void) {                                                      \
    log_header(#name " (" #T ")");                                             \
                                                                               \
    T *x = 0;                                                                  \
    bool pass = true;                                                          \
    for (size_t i = 0; i < sizeof test_vec_sizes / sizeof *test_vec_sizes;     \
         ++i) {                                                                \
      size_t n = test_vec_sizes[i];                                            \
      blas_ref_rand(n, &x);                                                    \
      T res = blas_nrm2(n, x);                                                 \
      T ref_res = blas_ref_nrm2(n, x);                                         \
      T mse = blas_ref_mse(1, &res, &ref_res);                                 \
      if (!(mse < ep*n)) {                                                       \
        pass = false;                                                          \
        if (verbosity > 0) {                                                   \
          log_info("N=%zu, MSE (%Lf) < ep (%Lf)", n, (long double)mse,         \
                   (long double)ep);                                           \
        }                                                                      \
        if (verbosity > 1) {                                                   \
          fputs("x = [", stderr);                                              \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)x[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fprintf(stderr, "result = %Lf\n", (long double)res);                 \
          fprintf(stderr, "reference result  = %Lf\n", (long double)ref_res);  \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    blas_ref_free(x);                                                          \
    if (pass) {                                                                \
      log_pass(0);                                                             \
      return 0;                                                                \
    } else {                                                                   \
      log_fail(0);                                                             \
      return -1;                                                               \
    }                                                                          \
  }

X(nrm2f_test, float);
X(nrm2_test, double);
X(nrm2l_test, long double);
#undef X

#define X(name, T)                                                             \
  static int name(void) {                                                      \
    log_header(#name " (" #T ")");                                             \
                                                                               \
    T *x_in = 0;                                                               \
    T *x = 0;                                                                  \
    T *xref = 0;                                                               \
    T *alpha = 0;                                                              \
    bool pass = true;                                                          \
    for (size_t i = 0; i < sizeof test_vec_sizes / sizeof *test_vec_sizes;     \
         ++i) {                                                                \
      size_t n = test_vec_sizes[i];                                            \
      blas_ref_rand(1, &alpha);                                                \
      blas_ref_rand(n, &x_in);                                                 \
      blas_ref_empty(n, &x);                                                   \
      blas_ref_empty(n, &xref);                                                \
      blas_ref_copy(n, x_in, x);                                               \
      blas_ref_copy(n, x_in, xref);                                            \
      blas_scal(n, alpha[0], x);                                               \
      blas_ref_scal(n, alpha[0], xref);                                        \
      T mse = blas_ref_mse(n, x, xref);                                        \
      if (!(mse < ep)) {                                                       \
        pass = false;                                                          \
        if (verbosity > 0) {                                                   \
          log_info("N=%zu, MSE (%Lf) < ep (%Lf)", n, (long double)mse,         \
                   (long double)ep);                                           \
        }                                                                      \
        if (verbosity > 1) {                                                   \
          fputs("x_in = [", stderr);                                              \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)x_in[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fprintf(stderr, "alpha = %Lf\n", (long double)alpha[0]);             \
          fputs("x      = [", stderr);                                         \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)x[i],                        \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
          fputs("xref = [", stderr);                                           \
          for (size_t i = 0; i < n; ++i) {                                     \
            if (verbosity < 3 && i > 10 && n - 10 > i) {                       \
              fprintf(stderr, "..., ");                                        \
              i = n - 10;                                                      \
            }                                                                  \
            fprintf(stderr, "%Lf%s", (long double)xref[i],                     \
                    i + 1 < n ? ", " : "]\n");                                 \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    blas_ref_free(alpha);                                                      \
    blas_ref_free(x_in);                                                       \
    blas_ref_free(x);                                                          \
    blas_ref_free(xref);                                                       \
    if (pass) {                                                                \
      log_pass(0);                                                             \
      return 0;                                                                \
    } else {                                                                   \
      log_fail(0);                                                             \
      return -1;                                                               \
    }                                                                          \
  }
X(scalf_test, float);
X(scal_test, double);
X(scall_test, long double);
#undef X

/*
double blas_asum(size_t n, double const *x);
float *blas_axpy(size_t n, float alpha, float const *x, float *y);
double blas_dot(size_t n, double const *x, double const *y);
double blas_nrm2(size_t n, double const *x);
float *blas_scal(size_t n, float alpha, float *x);
*/

int main(int argc, char *argv[]) {
  srand48(time(0));
  int opt;
  while ((opt = getopt(argc, argv, "hv")) != -1) {
    switch (opt) {
    case 'h':
      fprintf(stderr, "usage: %s [-hv] [test_name...]\n", basename(argv[0]));
      exit(0);
    case 'v':
      verbosity += 1;
      break;
    default:
      fprintf(stderr, "usage: %s [-hv] [test_name...]\n", basename(argv[0]));
      errx(1, "unrecognized option: %c", opt);
    }
  }

#define checktest(name)                                                        \
  if (strcmp(argv[i], #name) == 0) {                                           \
    test_utils_test(name##_test);                                               \
  }

  int i = optind;
  if (optind < argc) {
    for (int i = 1; i < argc; ++i) {
      checktest(asumf);
      checktest(asum);
      checktest(asuml);

      checktest(axpyf);
      checktest(axpy);
      checktest(axpyl);

      checktest(dotf);
      checktest(dot);
      checktest(dotl);

      checktest(nrm2f);
      checktest(nrm2);
      checktest(nrm2l);

      checktest(scalf);
      checktest(scal);
      checktest(scall);
    }
  } else {
    test_utils_test(asumf_test);
    test_utils_test(asum_test);
    test_utils_test(asuml_test);

    test_utils_test(axpyf_test);
    test_utils_test(axpy_test);
    test_utils_test(axpyl_test);

    test_utils_test(dotf_test);
    test_utils_test(dot_test);
    test_utils_test(dotl_test);

    test_utils_test(nrm2f_test);
    test_utils_test(nrm2_test);
    test_utils_test(nrm2l_test);

    test_utils_test(scalf_test);
    test_utils_test(scal_test);
    test_utils_test(scall_test);
  }
  return 0;
}
