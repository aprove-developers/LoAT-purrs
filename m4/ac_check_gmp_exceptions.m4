dnl Copyright (C) 2001 Roberto Bagnara <bagnara@cs.unipr.it>
dnl  
dnl This file is free software; as a special exception the author gives
dnl unlimited permission to copy and/or distribute it, with or without 
dnl modifications, as long as this notice is preserved.
dnl 
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
dnl implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
dnl
AC_DEFUN([AC_CHECK_GMP_EXCEPTIONS],
[AC_CACHE_CHECK([whether GMP has been compiled with support for exceptions],
ac_cv_gmp_supports_exceptions,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_RUN([
#include <gmp.h>
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
    mpz_t pie;
    mpz_init_set_str(pie, "3141592653589793238462643383279502884", 10);
  }
  catch (bad_alloc) {
    exit(0);
  }
  exit(1);
}
],
 ac_cv_gmp_supports_exceptions=yes,
 ac_cv_gmp_supports_exceptions=no,
 ac_cv_gmp_supports_exceptions=no)
 AC_LANG_RESTORE
])
gmp_supports_exceptions=${ac_cv_gmp_supports_exceptions}
if test x"$gmp_supports_exceptions" = xyes
then
  value=1
else
  value=0
fi
AC_DEFINE_UNQUOTED(GMP_SUPPORTS_EXCEPTIONS, $value,
  [Not zero if GMP has been compiled with support for exceptions.])
])
