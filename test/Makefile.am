# Include the headers from the project
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -lm

# Name of the test executable
bin_PROGRAMS = testbench

# Source files for the test program
testbench_SOURCES = testbench.c test_utils.c blas_ref.c

# Link against the library we built in src/
testbench_LDADD = ../src/libblas.la
