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

SUBDIRS = . Intervals

INCLUDES = -I$(top_builddir) -I$(top_srcdir)/src @extra_includes@

CPPFLAGS = -DPURRS_INTERNAL=public @extra_cpp_flags@

lib_LIBRARIES = libpurrs.a

# The ordering of the following list *matters*!
INCLUDE_FILES = \
globals.hh \
alg_eq_solver.hh \
rr_solver.hh

libpurrs_a_SOURCES = \
alg_eq_solver.cc \
rr_solver.cc \
globals.cc \
util.cc \
util.hh \
$(INCLUDE_FILES)

# The ordering of the following list *matters*!
DOC_FILES = \
alg_eq_solver.hh

libpurrs_a_LIBADD =

noinst_HEADERS = purrs.hh purrs_install.hh

BUILT_SOURCES = purrs.hh purrs_install.hh

MORE_INCLUDE_FILES = /dev/null

purrs.hh:	$(DOC_FILES)
	(cd $(srcdir) ; cat $(DOC_FILES)) \
		| sed -e "s/PURRS_INTERNAL/private/g" >purrs.hh

purrs_install.hh:	purrs_header.top $(top_builddir)/config.h purrs_header.middle $(INCLUDE_FILES) $(MORE_INCLUDE_FILES) purrs_header.bottom
	mkdir -p tmp
	(cd $(srcdir) ; cat $(MORE_INCLUDE_FILES) $(INCLUDE_FILES)) >tmp/concat
	grep -e "^#include <" <tmp/concat | sort -u >tmp/includes
	sed -e "/^#include/d" -e "s/PURRS_INTERNAL/public/g" <tmp/concat >tmp/noincludes.cc
	cat $(srcdir)/purrs_header.top $(top_builddir)/config.h \
		$(srcdir)/purrs_header.middle tmp/includes >purrs_install.hh
	$(CXXCPP) -P @debug_flag@ tmp/noincludes.cc >>purrs_install.hh
	cat $(srcdir)/purrs_header.bottom >>purrs_install.hh
	rm -rf tmp*

# purrs.hh and purrs_install.hh are not distributed.
dist-hook:
	rm -f $(distdir)/purrs.hh $(distdir)/purrs_install.hh

install-data-local: purrs_install.hh
	$(mkinstalldirs) $(DESTDIR)$(includedir)
	$(INSTALL_DATA) purrs_install.hh $(DESTDIR)$(includedir)/purrs.hh

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/purrs.hh

EXTRA_DIST = purrs_header.top purrs_header.middle purrs_header.bottom

MOSTLYCLEANFILES = $(BUILT_SOURCES)
