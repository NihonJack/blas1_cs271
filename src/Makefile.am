# Use libtool to build both static and shared
lib_LTLIBRARIES = libblas.la
libblas_la_SOURCES = asm/asum.s asm/axpy.s asm/dot.s asm/nrm2.s asm/scal.s asm/copy.s asm/swap.s
libblas_la_CFLAGS = -I$(top_srcdir)/include
libblas_la_LDFLAGS = -version-info 0:0:0
# -all-static will allow building static (.a) too
# -version-info is required for shared libraries

# Install headers
include_HEADERS = $(top_srcdir)/include/blas.h

CCAS = $(CC)
CCASFLAGS =
