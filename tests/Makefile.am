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

INCLUDES = \
-I$(top_builddir)/src \
-I$(top_srcdir)/src \
@extra_includes@

LDADD = $(top_builddir)/src/libpurrs.a ehandlers.o -lstdc++ -lntl @extra_libraries@ @readline_libraries@

CPPFLAGS = @extra_cpp_flags@

TESTS =

noinst_HEADERS = ehandlers.hh readlinebuf.hh

EXTRA_DIST = ehandlers.cc simpl_tests polynomials recurrences heap

BUGS =

noinst_PROGRAMS = simplify_driver eqs_driver rrs_driver gosper_driver rr_filter zeilberger_driver

simplify_driver_SOURCES = simplify_driver.cc
eqs_driver_SOURCES = eqs_driver.cc
rrs_driver_SOURCES = rrs_driver.cc
gosper_driver_SOURCES = gosper_driver.cc
rr_filter_SOURCES = rr_filter.cc
zeilberger_driver_SOURCES = zeilberger_driver.cc

check_PROGRAMS = $(TESTS) $(BUGS)

check-local: simplify_driver simpl_tests eqs_driver polynomials rrs_driver recurrences heap
	./simplify_driver -r <$(srcdir)/simpl_tests
	./eqs_driver <$(srcdir)/polynomials
	./rrs_driver <$(srcdir)/recurrences
	./rrs_driver -r <$(srcdir)/heap
