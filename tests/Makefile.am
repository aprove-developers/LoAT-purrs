# Copyright (C) 2001 Roberto Bagnara <bagnara@cs.unipr.it>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

srcdir = @srcdir@
VPATH = @srcdir@

@SET_MAKE@

SUBDIRS =

INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/CoStLy @extra_includes@

LDADD = $(top_builddir)/src/libpurrs.a ehandlers.o $(top_builddir)/CoStLy/libCoStLy.a -lstdc++ @extra_libraries@ -lntl @readline_libraries@

CPPFLAGS = @extra_cpp_flags@

TESTS =

noinst_HEADERS = ehandlers.hh

EXTRA_DIST = ehandlers.cc

BUGS =

noinst_PROGRAMS = eqs_driver rrs_driver simplify_driver gosper_driver

eqs_driver_SOURCES = eqs_driver.cc
rrs_driver_SOURCES = rrs_driver.cc
simplify_driver_SOURCES = simplify_driver.cc
gosper_driver_SOURCES = gosper_driver.cc

check_PROGRAMS = $(TESTS) $(BUGS)

check-local: eqs_driver polynomials rrs_driver recurrences
	./eqs_driver <$(srcdir)/polynomials
	./rrs_driver <$(srcdir)/recurrences
	./rrs_driver -r <$(srcdir)/heap
