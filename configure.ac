#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FDTD],[1.0],[ ])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([src/fdtd_config.h])

AC_LANG([C])
AC_PREFIX_DEFAULT(["/usr/local/fdtd"])

AM_INIT_AUTOMAKE([subdir-objects])

m4_include([m4/acx_mpi.m4])

AC_PROG_CC([gcc])
AC_PROG_CXX([g++])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_sf_coupling_3j])

AC_CHECK_LIB(boost_system, main, , [
    AC_CHECK_LIB(boost_system-mt, main, , [
        AC_MSG_ERROR("Linking against boost::system library failed.")
    ])
])
AC_CHECK_LIB(boost_filesystem, main, , [
    AC_CHECK_LIB(boost_filesystem-mt, main, , [
        AC_MSG_ERROR("Linking against boost::filesystem library failed.")
    ])
])
AC_CHECK_LIB(boost_serialization, main, , [
    AC_CHECK_LIB(boost_serialization-mt, main, , [
        AC_MSG_ERROR("Linking against boost::serialization library failed.")
    ])
])
AC_CHECK_LIB(boost_mpi, main, , [
    AC_CHECK_LIB(boost_mpi-mt, main, , [
        AC_MSG_ERROR("Linking against boost::mpi library failed.")
    ])
])

AX_CXX_COMPILE_STDCXX_11([noext])

AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH([include], [AS_HELP_STRING([--with-include],[include flags])], [with_include=$withval], [with_include=no])
if test "x${with_include}" != xno; then
    CFLAGS="${with_include} ${CFLAGS}"
    CXXFLAGS="${with_include} ${CXXFLAGS}"
    CPPFLAGS="${with_include} ${CPPFLAGS}"
fi
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([sqrt])
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
#LDFLAGS="$LDFLAGS -lboost_system -lboost_filesystem -lfftw3 -lm"
#lapack and blas
AC_ARG_ENABLE(mkl, [AS_HELP_STRING([--disable-mkl],[disable MKL extensions.])], [use_mkl=no], [use_mkl=yes])

AX_BOOST_BASE([1.40], [], [AC_MSG_ERROR("BOOST version >= 1.4.0 required.")])

AC_LANG_PUSH([C++])


AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],[Use MPI. Please make sure mpicxx and mpif90 are in $PATH])], [with_mpi=$withval], [with_mpi=no])
if test x${with_mpi} = xopenmpi ; then
    echo "..openmpi requested.."
elif test x${with_mpi} = xintel ; then
    echo "..Intel mpi requested.."
elif test x${with_mpi} = xmvapich ; then
    echo "..mvapich requested.."
elif test x${with_mpi} != xno ; then
    AC_MSG_ERROR("Please specify MPI implementation (openmpi, intel, mvapich)")
fi
if test x${with_mpi} != xno; then
    ACX_MPI
    AC_PROG_CC([$MPICC])
    AC_PROG_CXX([$MPICXX])
    LIBS="${MPILIBS} ${LIBS}"
    CXXFLAGS="${CXXFLAGS}"
fi

# LT_INIT([disable-static])
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH([include], [AS_HELP_STRING([--with-include],[include flags])], [with_include=$withval], [with_include=no])
if test "x${with_include}" != xno; then
    CFLAGS="${with_include} ${CFLAGS}"
    CXXFLAGS="${with_include} ${CXXFLAGS}"
    CPPFLAGS="${with_include} ${CPPFLAGS}"
fi

# AC_ARG_ENABLE(mkl, [AS_HELP_STRING([--enable-mkl],[enable MKL extensions.])], [use_mkl=yes], [use_mkl=no])
AC_ARG_ENABLE(scalapack, [AS_HELP_STRING([--enable-scalapack],[enable scalapack - used for parallel diagonalization])], [use_scalapack=yes], [use_scalapack=no])
# AC_ARG_ENABLE(acml, [AS_HELP_STRING([--enable-acml],[enable AMD ACML extensions.])], [use_acml=yes], [use_acml=no])

AC_ARG_WITH(blas, [AS_HELP_STRING([--with-blas],[Make sure the BLAS library is in$LDLIBRARY_PATH])], [blas=$withval], [blas=no])
if test x${blas} = xmkl ; then
    echo "..Intel's MKL requested.."
    AC_CHECK_HEADERS([mkl.h], [], [AC_MSG_ERROR([mkl.h not found or not working])], [])
    # AC_CHECK_LIB([mkl_core], main, LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm $LIBS",
    AC_CHECK_LIB([mkl_core], main, LIBS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -liomp5 -lgsl -lpthread -lm $LIBS",
    AC_MSG_ERROR([Cannot find Intel Math Kernel Library]), [-lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lpthread -lm])

   if test "x${use_scalapack}" = xyes; then
        if test x${with_mpi} = xopenmpi; then
            AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 $LIBS",
                [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
                [-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
            # AC_CHECK_LIB(mpiblacs,  main,  , [AC_MSG_ERROR("Linking against mpiblacs library failed.")])

        else
            AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 $LIBS",
                [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
                [-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
            # AC_CHECK_LIB(mpiblacs,  main,  , [AC_MSG_ERROR("Linking against mpiblacs library failed.")])

        fi
    fi
    CFLAGS="-DMKL ${CFLAGS}"
    CXXFLAGS="-DMKL ${CXXFLAGS}"
    CPPFLAGS="-DMKL ${CPPFLAGS}"
elif test x${blas} = xacml ; then
    echo "..AMD's ACML requested.."
    AH_TEMPLATE([HAVE_ACML], [the acml library will be linked.])
    AC_CHECK_LIB(acml, main,  [AC_DEFINE([HAVE_ACML]) LIBS="-lacml $LIBS"], [AC_MSG_ERROR("Linking against acml library failed.")])
    CFLAGS="-DACML ${CFLAGS}"
    CXXFLAGS="-DACML ${CXXFLAGS}"
    CPPFLAGS="-DACML ${CPPFLAGS}"
elif test x${blas} != xno ; then
    AC_MSG_ERROR("Please specify a BLAS Library (acml and mkl are known to work)")
fi

if test "x${use_scalapack}" = xyes; then
    if test x${with_mpi} != xno; then
        AH_TEMPLATE([HAVE_SCALAPACK], [a scalapack library will be linked.])
        if test x${blas} != xmkl; then
            AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])
            AC_CHECK_LIB(mpiblacs,  main,  , [AC_MSG_ERROR("Linking against mpiblacs library failed.")])
        fi
    else
        AC_MSG_ERROR("scalapack should be used with mpi")
    fi
fi

# if test x${use_mkl} = xyes; then
#     AC_CHECK_HEADERS([mkl.h], [], [AC_MSG_ERROR([mkl.h not found or not working])], [])
#     # AC_CHECK_LIB([mkl_core], main, LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm $LIBS",
#    AC_CHECK_LIB([mkl_core], main, LIBS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -liomp5 -lgsl -lpthread -lm $LIBS",
#     AC_MSG_ERROR([Cannot find Intel Math Kernel Library]), [-lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lpthread -lm])

#    if test "x${use_scalapack}" = xyes; then
#         if test x${with_mpi} = xopenmpi; then
#             AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 $LIBS",
#                 [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
#                 [-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
#             # AC_CHECK_LIB(mpiblacs,  main,  , [AC_MSG_ERROR("Linking against mpiblacs library failed.")])

#         else
#             AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 $LIBS",
#                 [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
#                 [-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
#             # AC_CHECK_LIB(mpiblacs,  main,  , [AC_MSG_ERROR("Linking against mpiblacs library failed.")])

#         fi
#     fi
# fi

if test "x${use_acml}" = xyes; then
    AH_TEMPLATE([HAVE_ACML], [the acml library will be linked.])
    AC_CHECK_LIB(acml, main,  [AC_DEFINE([HAVE_ACML]) LIBS="-lacml $LIBS"], [AC_MSG_ERROR("Linking against acml library failed.")])
fi

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 ])

AC_OUTPUT

