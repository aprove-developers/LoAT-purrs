# Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>
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

LDADD = $(top_builddir)/src/libpurrs.a -lstdc++ -lntl @extra_libraries@ -lcgicc

CPPFLAGS = @extra_cpp_flags@

noinst_PROGRAMS = purrs_eqs.cgi purrs_rrs.cgi

purrs_eqs_cgi_SOURCES = purrs_eqs_cgi.cc
purrs_rrs_cgi_SOURCES = purrs_rrs_cgi.cc
