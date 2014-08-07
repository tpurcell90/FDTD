#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([FDTD], 1.0, [ ])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([fdtd_config.h])

AC_LANG([C])
AC_PREFIX_DEFAULT(["/usr/local/fdtd"])

m4_include([m4/acx_mpi.m4])

AC_PROG_CC([gcc])
AC_PROG_CXX([g++])

AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_sf_coupling_3j])

AX_CXX_COMPILE_STDCXX_11([noext])

AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH([include], [AS_HELP_STRING([--with-include],[include flags])], [with_include=$withval], [with_include=no])
if test "x${with_include}" != xno; then
    CFLAGS="${with_include} ${CFLAGS}"
    CXXFLAGS="${with_include} ${CXXFLAGS}"
    CPPFLAGS="${with_include} ${CPPFLAGS}"
fi

#lapack and blas
AC_ARG_ENABLE(mkl, [AS_HELP_STRING([--disable-mkl],[disable MKL extensions.])], [use_mkl=no], [use_mkl=yes])

AX_BOOST_BASE([1.40], [], [AC_MSG_ERROR("BOOST version >= 1.4.0 required.")])

AC_LANG_PUSH([C++])

if test x${use_mkl} = xyes; then
    AC_CHECK_HEADERS([mkl.h], [], [AC_MSG_ERROR([mkl.h not found or not working])], [])
    AC_CHECK_LIB([mkl_intel_lp64], main, LIBS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -liomp5 -lgsl -lpthread -lm $LIBS",
    AC_MSG_ERROR([Cannot find Intel Math Kernel Library]), [-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -liomp5 -lpthread -lm])
fi

AC_CONFIG_FILES([Makefile
                 src/Makefile])

AC_OUTPUT