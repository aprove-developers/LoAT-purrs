# Copyright (C) 2002 Roberto Bagnara <bagnara@cs.unipr.it>
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

CPPFLAGS = @extra_cpp_flags@

noinst_LIBRARIES = libCoStLy.a

libCoStLy_a_SOURCES = \
cimath.cpp \
cimath.h \
cinterval.cpp \
cinterval.h \
error.h

EXTRA_DIST = COPYING README
