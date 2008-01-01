dnl Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>
dnl  
dnl This file is free software; as a special exception the author gives
dnl unlimited permission to copy and/or distribute it, with or without 
dnl modifications, as long as this notice is preserved.
dnl 
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
dnl implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
dnl
AC_DEFUN([AC_CHECK_USE_READLINE],
[
enableval=yes
AC_ARG_ENABLE(readline,
[  --disable-readline      do not use the GNU readline library])
case "${enableval}" in
yes)
  ac_cv_enable_readline=yes
  ;;
no)
  ac_cv_enable_readline=no
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-readline, needs yes or no])
  ;;
esac

if test x"$ac_cv_enable_readline" = xyes
then
  # Check for the header files.
  AC_CHECK_HEADERS(readline/readline.h readline/history.h)
  # Try to find the termcap library functions (tgetent et al.).
  AC_MSG_CHECKING([for termcap library functions])
  AC_MSG_RESULT([in termcap or curses or ncurses...])
  AC_CHECK_LIB(termcap, tgetent, ac_cv_termcap_lib=-ltermcap,
    [AC_CHECK_LIB(curses,  tgetent, ac_cv_termcap_lib=-lcurses,
      [AC_CHECK_LIB(ncurses, tgetent, ac_cv_termcap_lib=-lncurses,
        ac_cv_termcap_lib='')
      ])
    ])
  # Now check for the readline library.
  AC_CHECK_LIB(readline, readline,
    ac_cv_have_readline=yes, ac_cv_have_readline=no, $ac_cv_termcap_lib)
fi

if test x"$ac_cv_enable_readline" = xyes -a x"$ac_cv_have_readline" = xyes
then
  termcap_library="$ac_cv_termcap_lib"
  readline_libraries="-lreadline $termcap_library"
  AC_DEFINE(USE_READLINE, 1, [GNU readline will be used when this is defined])
else
  termcap_library=""
  readline_libraries=""
fi

AC_SUBST(termcap_library)
AC_SUBST(readline_libraries)
])
