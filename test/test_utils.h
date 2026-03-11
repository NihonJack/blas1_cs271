#ifndef TESTBENCH__H
#define TESTBENCH__H

extern void log_header(char const *fmt, ...);
extern void log_info(char const *fmt, ...);
extern void log_fail(char const *fmt, ...);
extern void log_pass(char const *fmt, ...);
extern int test_utils_test(int (*f)(void));

#endif
