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

INCLUDES = -I$(top_builddir)/src @extra_includes@

LDADD = $(top_builddir)/src/libpurrs.a -lginac -lcln @gmp_library_option@

CPPFLAGS = @extra_cpp_flags@

TESTS =

noinst_PROGRAMS = eqs_driver

eqs_driver_SOURCES = eqs_driver.cc

check_PROGRAMS = $(TESTS) $(BUGS)

# noinst_HEADERS = ehandlers.hh print.hh

BUGS = bug1

bug1_SOURCES = bug1.cc
