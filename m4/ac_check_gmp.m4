dnl A function to check for the existence and usability of GMP.
dnl Copyright (C) 2001, 2002 Roberto Bagnara <bagnara@cs.unipr.it>
dnl  
dnl This file is part of the Parma Polyhedra Library (PPL).
dnl 
dnl The PPL is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by the
dnl Free Software Foundation; either version 2 of the License, or (at your
dnl option) any later version.
dnl 
dnl The PPL is distributed in the hope that it will be useful, but WITHOUT
dnl ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
dnl FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl for more details.
dnl 
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
dnl USA.
dnl 
dnl For the most up-to-date information see the Parma Polyhedra Library
dnl site: http://www.cs.unipr.it/ppl/ .
dnl
AC_DEFUN([AC_CHECK_GMP],
[
AC_ARG_WITH(gmp-includes,
            [  --with-gmp-includes=DIR GMP include files are in DIR],
            gmp_includes=${with_gmp_includes}
            gmp_includes_option="-I${gmp_includes}")

gmp_library_option="-lgmp -lgmpxx"
AC_ARG_WITH(gmp-dir,
            [  --with-gmp-dir=DIR      GMP library files are in DIR],
            gmp_dir=${with_gmp_dir}
            gmp_library_option="-L${gmp_dir} ${gmp_library_option}")

ac_save_LIBS="$LIBS"
LIBS="${gmp_library_option} $LIBS"
ac_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="${gmp_includes_option} $CPPFLAGS"
AC_LANG_SAVE
AC_LANG_CPLUSPLUS

AC_MSG_CHECKING([for the GMP library])
AC_TRY_RUN([
#include <gmpxx.h>

using namespace std;

int main() {
  mpz_class pie("3141592653589793238462643383279502884");
  exit(0);
}
],
  AC_MSG_RESULT(yes)
  ac_cv_have_gmp=yes,
  AC_MSG_RESULT(no)
  ac_cv_have_gmp=no,
  AC_MSG_RESULT(no)
  ac_cv_have_gmp=no)

have_gmp=${ac_cv_have_gmp}

if test x"$ac_cv_have_gmp" = xyes
then

AC_MSG_CHECKING([whether GMP has been compiled with support for exceptions])
AC_TRY_RUN([
#include <gmpxx.h>
#include <new>
#include <cstddef>
#include <cstdlib>

using namespace std;

static void*
x_malloc(size_t) {
  throw bad_alloc();
}

static void*
x_realloc(void*, size_t, size_t) {
  throw bad_alloc();
}

static void
x_free(void*, size_t) {
}

int main() {
  mp_set_memory_functions(x_malloc, x_realloc, x_free);
  try {
    mpz_class pie("3141592653589793238462643383279502884");
  }
  catch (bad_alloc) {
    exit(0);
  }
  exit(1);
}
],
  AC_MSG_RESULT(yes)
  ac_cv_gmp_supports_exceptions=yes,
  AC_MSG_RESULT(no)
  ac_cv_gmp_supports_exceptions=no,
  AC_MSG_RESULT(no)
  ac_cv_gmp_supports_exceptions=no)

gmp_supports_exceptions=${ac_cv_gmp_supports_exceptions}
if test x"$gmp_supports_exceptions" = xyes
then
  value=1
else
  value=0
fi
AC_DEFINE_UNQUOTED(GMP_SUPPORTS_EXCEPTIONS, $value,
  [Not zero if GMP has been compiled with support for exceptions.])

fi

AC_LANG_RESTORE
CPPFLAGS="$ac_save_CPPFLAGS"
LIBS="$ac_save_LIBS"
])
