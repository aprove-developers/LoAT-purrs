# Autoconf source file for the Parma Recurrence Relation Solver.
# Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# Process this file with autoconf to produce a configure script.

# Every other copy of the package version number gets its value from here.
AC_INIT([the Parma Recurrence Relation Solver], 0.0.3, purrs-devel@cs.unipr.it, purrs)

# Minimum autoconf version required.
AC_PREREQ(2.61)

# Make sure the sources are there.
AC_CONFIG_SRCDIR(src/Recurrence.cc)

# Specify the location of additional local Autoconf macros.
AC_CONFIG_MACRO_DIR([m4])

# Determine the system type and set output variables to the names of
# the canonical system types.
AC_CANONICAL_HOST

# Use Automake.
AM_INIT_AUTOMAKE([foreign dist-bzip2 dist-zip 1.5])
AC_SUBST(VERSION)

# Generate a configuration header file.
AM_CONFIG_HEADER(config.h)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AH_TOP([
#ifndef _purrs_config_h
#define _purrs_config_h 1
])

AH_BOTTOM([

// Neutralize relops for GCC 2.96.
#ifndef __SGI_STL_INTERNAL_RELOPS
#define __SGI_STL_INTERNAL_RELOPS
#endif

#if __GNUC__ >= 3
// Flexible array members at the end of a struct
// are defined as in
//   Type array[];
#define PURRS_FLEXIBLE_ARRAY
#else
// Flexible array members at the end of a struct
// are defined as in
//   Type array[0];
#define PURRS_FLEXIBLE_ARRAY 0
#endif

#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#define BIG_ENDIAN 4321
#endif

#ifndef __BYTE_ORDER
#ifdef WORDS_BIGENDIAN
#define __BYTE_ORDER BIG_ENDIAN
#else
#define __BYTE_ORDER LITTLE_ENDIAN
#endif
#endif

#define BOUNDARY_TYPE_IS_DOUBLE 1

#define INLINE inline

#endif
])

# C compiler

AC_ARG_WITH(cc,
[  --with-cc=xxx           use XXX as the C compiler],
CC=$with_cc)

AC_LANG(C)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# C++ compiler

AC_ARG_WITH(cxx,
[  --with-cxx=xxx          use XXX as the C++ compiler],
CXX=$with_cxx)

AC_LANG(C++)

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_MAKE_SET
AC_PROG_INSTALL

# GMP library

AC_CHECK_GMP

if test x"$have_gmp" = xno
then
  AC_MSG_ERROR([Cannot find GMP version 4 or higher.
GMP is the GNU Multi-Precision library:
see http://www.swox.com/gmp/ for more information.
Do not forget to enable the C++ interface: add --enable-cxx
to the configuration options.])
fi

AC_CHECK_SIZEOF(mp_limb_t, 4, [#include <gmp.h>])

extra_includes="${extra_includes} ${gmp_includes_option}"
extra_libraries="${extra_libraries} ${gmp_library_option}"


# CLN library

AC_LIB_LINKFLAGS([cln])
extra_libraries="${extra_libraries} ${LTLIBCLN}"

# GiNaC library

AC_LIB_LINKFLAGS([ginac])
extra_libraries="${extra_libraries} ${LTLIBGINAC}"

# Compilation flags

COMP_FLAGS=""
OPT_FLAGS=""

enableval=yes
AC_MSG_CHECKING([whether to compile with debug info])
AC_ARG_ENABLE(debugging,
[  --enable-debugging      compile with debugging information])
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  OPT_FLAGS="-g"
  ;;
no)
  AC_MSG_RESULT(no)
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-debugging, needs yes or no])
  ;;
esac

enableval=no
AC_MSG_CHECKING([whether to compile for profiling])
AC_ARG_ENABLE(profiling,
[  --enable-profiling      compile for profiling])
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  OPT_FLAGS="-g"
  COMP_FLAGS="$COMP_FLAGS -pg"
  ;;
no)
  AC_MSG_RESULT(no)
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-profiling, needs yes or no])
  ;;
esac

enableval=no
AC_MSG_CHECKING([whether to compile for test coverage])
AC_ARG_ENABLE(coverage,
[  --enable-coverage       compile for test coverage])
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  OPT_FLAGS="-g"
  COMP_FLAGS="$COMP_FLAGS -fprofile-arcs -ftest-coverage"
  ;;
no)
  AC_MSG_RESULT(no)
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-coverage, needs yes or no])
  ;;
esac

enableval=yes
AC_MSG_CHECKING([whether to enable checking of run-time assertions])
AC_ARG_ENABLE(assertions,
[  --enable-assertions     check run-time assertions])
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  ;;
no)
  AC_MSG_RESULT(no)
  AC_DEFINE(NDEBUG, 1, [Assertions are disabled when this is defined])
  debug_flag="-DNDEBUG=1"
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} --enable-assertions, needs yes or no])
  ;;
esac

arch=no
enableval=no
AC_MSG_CHECKING([whether to enable optimizations])
AC_ARG_ENABLE(optimization,
[  --enable-optimization   enable compiler optimizations])
case "${enableval}" in
sspeed)
  AC_MSG_RESULT(sspeed)
  OPT_FLAGS="$OPT_FLAGS -O3 -fomit-frame-pointer"
  arch=yes
  ;;
speed)
  AC_MSG_RESULT(speed)
  OPT_FLAGS="$OPT_FLAGS -O3"
  arch=yes
  ;;
size)
  AC_MSG_RESULT(speed)
  OPT_FLAGS="$OPT_FLAGS -Os"
  arch=yes
  ;;
standard)
  AC_MSG_RESULT(standard)
  OPT_FLAGS="$OPT_FLAGS -O2"
  ;;
yes)
  AC_MSG_RESULT(standard)
  OPT_FLAGS="$OPT_FLAGS -O2"
  ;;
no)
  AC_MSG_RESULT(no)
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} --enable-optimization, needs sspeed, speed, size, standard, yes or no])
  ;;
esac

enableval=$arch
AC_MSG_CHECKING([for which architecture to optimize])
AC_ARG_ENABLE(arch,
[  --enable-arch[=arch]      optimize for architecture arch])
case "${enableval}" in
yes)
  m=`uname -m`
  case $m in
  i?86 | k6 | athlon)
    AC_MSG_RESULT($m)
    OPT_FLAGS="$OPT_FLAGS -march=$m"
    ;;
  *)
    AC_MSG_RESULT(default)
    ;;
  esac
  ;;
no)
  AC_MSG_RESULT(default)
  ;;
*)
  AC_MSG_RESULT($enableval)
  OPT_FLAGS="$OPT_FLAGS -march=$enableval"
  ;;
esac

CFLAGS="$COMP_FLAGS $OPT_FLAGS"
CXXFLAGS="$COMP_FLAGS $OPT_FLAGS"

# Allow additions to C and C++ compilation flags

AC_ARG_WITH(cflags,
[  --with-cflags=xxx       add XXX to the options for the C compiler],
CFLAGS="$CFLAGS $with_cflags")

AC_ARG_WITH(cxxflags,
[  --with-cxxflags=xxx     add XXX to the options for the C++ compiler],
CXXFLAGS="$CXXFLAGS $with_cxxflags")


# Checks for machine architecture features
AC_C_BIGENDIAN

# Libtool.
# Turn off shared libraries during beta-testing, since they
# make debugging harder and the build process takes too long.
AC_DISABLE_SHARED
AC_LIBTOOL_DLOPEN
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

# Checks for typedefs, structures, and compiler characteristics.
# Nothing for the time being.

# Checks whether the readline library should be used.
AC_CHECK_USE_READLINE

# Checks for library functions.
AC_HEADER_TIME
AC_CHECK_HEADERS(unistd.h sys/resource.h sys/time.h)
AC_CHECK_FUNCS(nice getrusage getrlimit setrlimit)

# If we are using GCC we want to compile with warnings enabled.
if test x"$GCC" = xyes
then
  CFLAGS="$CFLAGS -W -Wall"
fi
if test x"$GXX" = xyes
then
  CXXFLAGS="$CXXFLAGS -W -Wall"
fi

AC_SUBST(extra_includes)
AC_SUBST(extra_cpp_flags)
AC_SUBST(extra_libraries)
AC_SUBST(debug_flag)

AC_OUTPUT(Makefile
          src/Makefile src/Intervals/Makefile src/Intervals/sysdep/Makefile
          m4/Makefile
          doc/Makefile
          tests/Makefile
          demos/Makefile
          doc/user-browse.doxyconf-latex doc/devref-browse.doxyconf-latex
          doc/user-print.doxyconf-latex  doc/devref-print.doxyconf-latex
          doc/user.doxyconf-html         doc/devref.doxyconf-html)
#          purrs.lsm purrs.spec

if test x"${ac_cv_header_readline_readline_h}" = xno -o x"${ac_cv_header_readline_history_h}" = xno -o x"${ac_cv_lib_readline_readline}" = xno
then
  AC_MSG_WARN([+++ CANNOT FIND THE READLINE LIBRARY:])
  AC_MSG_WARN([+++ USING SOME TEST PROGRAMS WILL BE UNEASY.])
else
  AC_DEFINE(USE_READLINE, 1, [The readline library is used if this is defined])
fi
