#define POSIX_C_SOURCE 200809L
#define _DEFAULT_SOURCE
#include <err.h>
#include <setjmp.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "color.h"
#include "test_utils.h"

static volatile sig_atomic_t signum = 0;
static sigjmp_buf sigenv;

static void log_headerv(char const *fmt, va_list ap);
static void log_infov(char const *fmt, va_list ap);
static void log_failv(char const *fmt, va_list ap);
static void log_passv(char const *fmt, va_list ap);
static void sighandler(int sig);

static void log_headerv(char const *fmt, va_list ap) {
  fputs(C_BOLD "=== ", stderr);
  vfprintf(stderr, fmt, ap);
  fputs(C_RESET " ===\n", stderr);
}

void log_header(char const *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  log_headerv(fmt, ap);
  va_end(ap);
}

static void log_infov(char const *fmt, va_list ap) {
  fputs(C_F_YELLOW "INFO: ", stderr);
  vfprintf(stderr, fmt, ap);
  fputs(C_RESET "\n", stderr);
}

void log_info(char const *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  log_infov(fmt, ap);
  va_end(ap);
}

static void log_failv(char const *fmt, va_list ap) {
  if (fmt) {
    fputs(C_F_RED "FAIL: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputs(C_RESET "\n", stderr);
  } else {
    fputs(C_F_RED "FAIL\n" C_RESET, stderr);
  }
}

void log_fail(char const *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  log_failv(fmt, ap);
  va_end(ap);
}

static void log_passv(char const *fmt, va_list ap) {
  if (fmt) {
    fputs(C_F_GREEN "PASS: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputs(C_RESET "\n", stderr);
  } else {
    fputs(C_F_GREEN "PASS\n" C_RESET, stderr);
  }
}

void log_pass(char const *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  log_passv(fmt, ap);
  va_end(ap);
}

static void sighandler(int sig) {
  signum = sig;
  siglongjmp(sigenv, sig);
}

int test_utils_test(int (*f)(void)) {
  struct sigaction old_sa[NSIG];
  struct sigaction sa;
  sa.sa_handler = sighandler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = 0;

  for (int i = 1; i < NSIG; i++) {
    if (i == SIGKILL || i == SIGSTOP)
      continue;
    sigaction(i, &sa, &old_sa[i]);
  }
  volatile int res = -1;
  switch (sigsetjmp(sigenv, 1)) {
  case 0: {
    res = f();
  } break;
  default:
    warnx("signal caught %s", strsignal(signum));
    break;
  }
  for (int i = 1; i < NSIG; i++) {
    if (i == SIGKILL || i == SIGSTOP)
      continue;

    sigaction(i, &old_sa[i], NULL);
  }
  return res;
}
