#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.56])
AC_INIT([robospect], [2.1], [watersc1@ifa.hawaii.edu])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_MAKE_SET

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lc':
AC_CHECK_LIB([c], [getgid])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [cos])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([Unable to compile without GSL])])
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [gsl_blas_dgemm], [], [AC_MSG_ERROR([Unable to compile without GSL])])
# FIXME: Replace `main' with a function in `-lplplotd':
AC_CHECK_LIB([plplotd], [c_plinit])

AC_CHECK_LIB([cfitsio], [ffdkopn])
# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])
AH_TEMPLATE([ODD_CFITSIO], [fitsio.h is located in odd location])
AC_CHECK_HEADERS(cfitsio/fitsio.h,AC_DEFINE([ODD_CFITSIO],[0],), )

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt strcasecmp strdup])

AC_CONFIG_FILES([Makefile src/Makefile])

AM_INIT_AUTOMAKE

AC_OUTPUT
