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

SUBDIRS = sysdep

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/Intervals/sysdep @extra_includes@

CPPFLAGS = @extra_cpp_flags@

noinst_LIBRARIES = libIntervals.a

libIntervals_a_SOURCES = \
IEEE754.h \
round.h \
Bool3.cc \
Bool3.types.hh \
Bool3.defs.hh \
Bool3.inlines.hh \
Boundary.cc \
Boundary.types.hh \
Boundary.defs.hh \
Boundary.inlines.hh \
Interval.cc \
Interval.types.hh \
Interval.defs.hh \
Interval.inlines.hh
