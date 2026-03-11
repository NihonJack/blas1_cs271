AC_INIT([blas], [0.1], [ryan.gambord@oregonstate.edu])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_AR

AX_CHECK_COMPILE_FLAG([-std=c11],
  [CFLAGS="$CFLAGS -std=c11"],
  [AC_MSG_ERROR([Your compiler does not support the required C11 standard. Please use a C11-compliant compiler.])],
  [-Werror]
)


AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile])

LT_INIT

AC_OUTPUT
