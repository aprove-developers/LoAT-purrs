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
srcpath = @srcpath@
VPATH = @srcdir@

@SET_MAKE@

LICENSES = \
gpl.pdf \
gpl.ps.gz \
fdl.pdf \
fdl.ps.gz

HTML_DIRS = \
purrs-user-@VERSION@-html \
purrs-devref-@VERSION@-html

ALL_TARGETS = \
$(HTML_DIRS) \
$(LICENSES) \
purrs-user-@VERSION@-browse.pdf \
purrs-user-@VERSION@-print.pdf \
purrs-user-@VERSION@-print.ps.gz \
purrs-user-@VERSION@-html.tar.gz \
purrs-devref-@VERSION@-browse.pdf \
purrs-devref-@VERSION@-print.pdf \
purrs-devref-@VERSION@-print.ps.gz \
purrs-devref-@VERSION@-html.tar.gz

# We distribute some preprocessed user's documentation,
# both for browsing and for printing.
DIST_DOCS = \
purrs-user-@VERSION@-browse.pdf \
purrs-user-@VERSION@-print.pdf \
purrs-user-@VERSION@-print.ps.gz

# Do nothing by default.
all:

dist-hook: $(DIST_DOCS)

user-browse: purrs-user-@VERSION@-browse.pdf

world: $(ALL_TARGETS)

PURRS_SOURCE_FILES = $(wildcard $(srcdir)/../src/*.hh $(srcdir)/../src/*.cc)

DOX_FILES = \
gpl.dox \
fdl.dox

USER_STUFF = $(DOX_FILES) user.tex $(top_builddir)/src/purrs.hh

DEVREF_STUFF = $(DOX_FILES) devref.tex $(PURRS_SOURCE_FILES)

LATEX_DIRS = \
user-browse.latex-dir \
user-print.latex-dir \
devref-browse.latex-dir \
devref-print.latex-dir

TEX_ENV = TEXINPUTS=$(srcpath):$(TEXINPUTS)

.SECONDARY: $(LATEX_DIRS)

user-%.latex-dir: user-%.doxyconf-latex $(USER_STUFF)
	rm -rf $@
	$(TEX_ENV) doxygen $<

devref-%.latex-dir: devref-%.doxyconf-latex $(DEVREF_STUFF)
	rm -rf $@
	$(TEX_ENV) doxygen $<

purrs-user-@VERSION@-%.pdf: user-%.latex-dir
	$(TEX_ENV) make -C $< pdf && mv -f $</refman.pdf $@

purrs-devref-@VERSION@-%.pdf: devref-%.latex-dir
	$(TEX_ENV) make -C $< pdf && mv -f $</refman.pdf $@

purrs-user-@VERSION@-%.ps: user-%.latex-dir
	$(TEX_ENV) make -C $< ps && mv -f $</refman.ps $@

purrs-devref-@VERSION@-%.ps: devref-%.latex-dir
	$(TEX_ENV) make -C $< ps && mv -f $</refman.ps $@

purrs-%-@VERSION@-html: %.doxyconf-html $(DEVREF_STUFF)
	rm -rf html
	$(TEX_ENV) doxygen $< && rm -rf $@ && mv -f html $@

%.tar: %
	tar cf $@ $<

%.gz: %
	gzip --best --force $<

%.pdf: %.tex
	$(TEX_ENV) pdflatex $<

%.dvi: %.tex
	$(TEX_ENV) latex $<

%.ps:	%.dvi
	dvips -o $@ $<

EXTRA_DIST= \
$(DIST_DOCS) \
$(DOX_FILES) \
user-browse.doxyconf-latex \
devref-browse.doxyconf-latex \
user-print.doxyconf-latex \
devref-print.doxyconf-latex \
user.doxyconf-html \
devref.doxyconf-html \
user.tex \
devref.tex \
gpl.tex \
fdl.tex \
purrs.sty

mostlyclean-local:
	rm -rf $(LATEX_DIRS) $(HTML_DIRS)

CLEANFILES = $(DIST_DOCS)
