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

INCLUDES = -I$(top_builddir) -I$(top_srcdir)/src -I$(top_srcdir)/CoStLy @extra_includes@

CPPFLAGS = @extra_cpp_flags@

lib_LIBRARIES = libpurrs.a

# The ordering of the following list *matters*!
INCLUDE_FILES = \
CInterval.types.hh \
complint.hh \
Complex_Interval.types.hh \
Bool3.types.hh \
Constant.types.hh \
Expr.types.hh \
Expr_List.types.hh \
Init.types.hh \
Matrix.types.hh \
Number.types.hh \
Finite_Order_Info.types.hh \
Functional_Equation_Info.types.hh \
Blackboard.types.hh \
Cached_Expr.types.hh \
Non_Linear_Info.types.hh \
Infinite_Order_Info.types.hh \
Recurrence.types.hh \
Symbol.types.hh \
Complex_Interval.defs.hh \
Complex_Interval.inlines.hh \
Bool3.defs.hh \
Bool3.inlines.hh \
Constant.defs.hh \
Constant.inlines.hh \
Number.defs.hh \
Number.inlines.hh \
Symbol.defs.hh \
Symbol.inlines.hh \
Expr.defs.hh \
Expr_List.defs.hh \
Expr.inlines.hh \
Expr_List.inlines.hh \
Matrix.defs.hh \
Matrix.inlines.hh \
Init.defs.hh \
Init.inlines.hh \
initializer.hh \
approximate.hh \
approximate_impl.hh \
size_norm.hh \
size_norm_impl.hh \
alg_eq_solver.hh \
poly_factor.hh \
zeilberger.hh \
gosper.hh \
simplify.hh \
numerator_denominator.hh \
factorize.hh \
factorize_giac.hh \
compute_prod.hh \
ep_decomp.hh \
sum_poly.hh \
finite_order.hh \
functional_equation.hh \
util.hh \
Finite_Order_Info.defs.hh \
Finite_Order_Info.inlines.hh \
Functional_Equation_Info.defs.hh \
Functional_Equation_Info.inlines.hh \
Blackboard.defs.hh \
Blackboard.inlines.hh \
Cached_Expr.defs.hh \
Cached_Expr.inlines.hh \
Non_Linear_Info.defs.hh \
Non_Linear_Info.inlines.hh \
Infinite_Order_Info.defs.hh \
Infinite_Order_Info.inlines.hh \
Recurrence.defs.hh \
Recurrence.inlines.hh

libpurrs_a_SOURCES = \
Finite_Order_Info.cc \
Functional_Equation_Info.cc \
Blackboard.cc \
Cached_Expr.cc \
Non_Linear_Info.cc \
Infinite_Order_Info.cc \
Complex_Interval.cc \
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
approximate.cc \
size_norm.cc \
complint.cc \
zeilberger.cc \
gosper.cc \
numerator_denominator.cc \
factorize.cc \
factorize_giac.cc \
poly_factor.cc \
finite_order.cc \
functional_equation.cc \
rr_classifier.cc \
simplify.cc \
compute_prod.cc \
ep_decomp.cc \
sum_poly.cc \
compute_bounds.cc \
compute_exact_solution.cc \
util.cc \
$(INCLUDE_FILES)

# The ordering of the following list *matters*!
DOC_FILES = \
Expr.types.hh \
Symbol.types.hh \
Number.types.hh \
Constant.types.hh \
Recurrence.types.hh \
Expr.defs.hh \
Symbol.defs.hh \
Number.defs.hh \
Constant.defs.hh \
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
