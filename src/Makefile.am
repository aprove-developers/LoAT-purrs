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

INCLUDES = -I$(top_builddir) -I$(top_srcdir)/src @extra_includes@

CPPFLAGS = @extra_cpp_flags@

lib_LIBRARIES = libpurrs.a

# The ordering of the following list *matters*!
INCLUDE_FILES = \
globals.hh \
Bool3.types.hh \
Expr.types.hh \
Constant.types.hh \
Symbol.types.hh \
Number.types.hh \
Expr_List.types.hh \
Init.types.hh \
Matrix.types.hh \
Finite_Order_Info.types.hh \
Functional_Equation_Info.types.hh \
Blackboard.types.hh \
Cached_Expr.types.hh \
Non_Linear_Info.types.hh \
Weighted_Average_Info.types.hh \
Recurrence.types.hh \
Symbol.defs.hh \
Symbol.inlines.hh \
Number.defs.hh \
Number.inlines.hh \
Constant.defs.hh \
Constant.inlines.hh \
Expr.defs.hh \
Expr_List.defs.hh \
Expr.inlines.hh \
Expr_List.inlines.hh \
Bool3.defs.hh \
Bool3.inlines.hh \
Blackboard.defs.hh \
Blackboard.inlines.hh \
Matrix.defs.hh \
Matrix.inlines.hh \
Init.defs.hh \
Init.inlines.hh \
initializer.hh \
size_norm.hh \
size_norm_impl.hh \
alg_eq_solver.hh \
poly_factor.hh \
zeilberger.hh \
gosper.hh \
simplify.hh \
factorize.hh \
factorize_giac.hh \
compute_prod.hh \
ep_decomp.hh \
fact_decomp.hh \
sum_poly.hh \
finite_order.hh \
functional_equation.hh \
util.hh \
Cached_Expr.defs.hh \
Cached_Expr.inlines.hh \
Recurrence.defs.hh \
Finite_Order_Info.defs.hh \
Finite_Order_Info.inlines.hh \
Functional_Equation_Info.defs.hh \
Functional_Equation_Info.inlines.hh \
Non_Linear_Info.defs.hh \
Non_Linear_Info.inlines.hh \
Weighted_Average_Info.defs.hh \
Weighted_Average_Info.inlines.hh \
Recurrence.inlines.hh

libpurrs_a_SOURCES = \
Finite_Order_Info.cc \
Functional_Equation_Info.cc \
Blackboard.cc \
Cached_Expr.cc \
Non_Linear_Info.cc \
Weighted_Average_Info.cc \
Bool3.cc \
Constant.cc \
Expr.cc \
Expr_List.cc \
Init.cc \
Matrix.cc \
Number.cc \
Recurrence.cc \
Symbol.cc \
alg_eq_solver.cc \
size_norm.cc \
globals.cc \
zeilberger.cc \
gosper.cc \
factorize.cc \
factorize_giac.cc \
poly_factor.cc \
finite_order.cc \
functional_equation.cc \
rr_classifier.cc \
simplify.cc \
compute_prod.cc \
ep_decomp.cc \
fact_decomp.cc \
sum_poly.cc \
compute_bounds.cc \
compute_exact_solution.cc \
verify.cc \
util.cc \
$(INCLUDE_FILES)

# The ordering of the following list *matters*!
DOC_FILES = \
globals.hh \
Symbol.types.hh \
Number.types.hh \
Constant.types.hh \
Expr.types.hh \
Recurrence.types.hh \
Symbol.defs.hh \
Number.defs.hh \
Constant.defs.hh \
Expr.defs.hh \
Blackboard.defs.hh \
Recurrence.defs.hh

libpurrs_a_LIBADD =

noinst_HEADERS = purrs.hh purrs_install.hh

BUILT_SOURCES = purrs.hh purrs_install.hh

MORE_INCLUDE_FILES = /dev/null

purrs.hh:	$(DOC_FILES) Makefile
	(cd $(srcdir) ; cat $(DOC_FILES)) \
		>purrs.hh

purrs_install.hh:	purrs_header.top $(top_builddir)/config.h purrs_header.middle $(INCLUDE_FILES) $(MORE_INCLUDE_FILES) purrs_header.bottom Makefile
	mkdir -p tmp
	(cd $(srcdir) ; cat $(MORE_INCLUDE_FILES) $(INCLUDE_FILES)) >tmp/concat
	grep -e "^#include <" <tmp/concat | sort -u >tmp/includes
	sed -e "/^#include/d" <tmp/concat >tmp/noincludes.cc
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
