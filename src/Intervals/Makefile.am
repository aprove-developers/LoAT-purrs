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

#strtod.cc \

#pow.cc \
#id_double.cc \
#Rel.cc \
#Bool3.cc

# EXTRA_DIST = \
# Action.cc \
# ActionSet.cc \
# BindingCell.cc \
# BindingCellBag.cc \
# CNet.cc \
# CNet1.cc \
# CNet2.cc \
# CNet_out.cc \
# CacheCell.cc \
# CacheCellBag.cc \
# Expr.cc \
# Expr_out.cc \
# ExprRefSet.cc \
# HistoryCell.cc \
# HistoryCellBag.cc \
# LinConstraintCell.cc \
# LinConstraintCellBag.cc \
# Numeric.cc \
# Qty_Flag.cc \
# QtyArithCell.cc \
# QtyArithCellBag.cc \
# QtyCell.cc \
# QtyCellBag.cc \
# QtyCellBag_Iterator.cc \
# QtyCell_Key.cc \
# QtyCell_Key_Map.cc \
# QtyCell_Key_Set.cc \
# RelArithCell.cc \
# RelArithCellBag.cc \
# RelCell.cc \
# RelCellBag.cc \
# RelCellBag_Iterator.cc \
# RelCell_Key.cc \
# RelCell_Key_Set.cc \
# TStamp.cc \
# Action.types.hh \
# ActionSet.types.hh \
# BindingCell.types.hh \
# BindingCellBag.types.hh \
# Bool3.types.hh \
# Boundary.types.hh \
# CNet.types.hh \
# CacheCell.types.hh \
# CacheCellBag.types.hh \
# Expr.types.hh \
# ExprRefSet.types.hh \
# HistoryCell.types.hh \
# HistoryCellBag.types.hh \
# Interval.types.hh \
# LinConstraintCell.types.hh \
# LinConstraintCellBag.types.hh \
# Qty.types.hh \
# Qty_Flag.types.hh \
# QtyArithCell.types.hh \
# QtyArithCellBag.types.hh \
# QtyCell.types.hh \
# QtyCellBag.types.hh \
# QtyCellBag_Iterator.types.hh \
# QtyCell_Key.types.hh \
# QtyCell_Key_Map.types.hh \
# QtyCell_Key_Set.types.hh \
# Rel.types.hh \
# RelArithCell.types.hh \
# RelArithCellBag.types.hh \
# RelCell.types.hh \
# RelCellBag.types.hh \
# RelCellBag_Iterator.types.hh \
# RelCell_Key.types.hh \
# RelCell_Key_Set.types.hh \
# TStamp.types.hh \
# Action.defs.hh \
# ActionSet.defs.hh \
# BindingCell.defs.hh \
# BindingCellBag.defs.hh \
# Bool3.defs.hh \
# Boundary.defs.hh \
# CNet.defs.hh \
# CacheCell.defs.hh \
# CacheCellBag.defs.hh \
# Expr.defs.hh \
# ExprRefSet.defs.hh \
# HistoryCell.defs.hh \
# HistoryCellBag.defs.hh \
# Interval.defs.hh \
# LinConstraintCell.defs.hh \
# LinConstraintCellBag.defs.hh \
# Numeric.defs.hh \
# Ptr.defs.hh \
# Qty.defs.hh \
# Qty_Flag.defs.hh \
# QtyArithCell.defs.hh \
# QtyArithCellBag.defs.hh \
# QtyCell.defs.hh \
# QtyCellBag.defs.hh \
# QtyCellBag_Iterator.defs.hh \
# QtyCell_Key.defs.hh \
# QtyCell_Key_Map.defs.hh \
# QtyCell_Key_Set.defs.hh \
# Rel.defs.hh \
# RelArithCell.defs.hh \
# RelArithCellBag.defs.hh \
# RelCell.defs.hh \
# RelCellBag.defs.hh \
# RelCellBag_Iterator.defs.hh \
# RelCell_Key.defs.hh \
# RelCell_Key_Set.defs.hh \
# TStamp.defs.hh \
# Action.inlines.hh \
# ActionSet.inlines.hh \
# BindingCell.inlines.hh \
# BindingCellBag.inlines.hh \
# Bool3.inlines.hh \
# Boundary.inlines.hh \
# CNet.inlines.hh \
# CacheCell.inlines.hh \
# CacheCellBag.inlines.hh \
# Expr.inlines.hh \
# ExprRefSet.inlines.hh \
# HistoryCell.inlines.hh \
# HistoryCellBag.inlines.hh \
# Interval.inlines.hh \
# LinConstraintCell.inlines.hh \
# LinConstraintCellBag.inlines.hh \
# Numeric.inlines.hh \
# Qty.inlines.hh \
# Qty_Flag.inlines.hh \
# QtyArithCell.inlines.hh \
# QtyArithCellBag.inlines.hh \
# QtyCell.inlines.hh \
# QtyCellBag.inlines.hh \
# QtyCellBag_Iterator.inlines.hh \
# QtyCell_Key.inlines.hh \
# QtyCell_Key_Map.inlines.hh \
# QtyCell_Key_Set.inlines.hh \
# Rel.inlines.hh \
# RelArithCell.inlines.hh \
# RelArithCellBag.inlines.hh \
# RelCell.inlines.hh \
# RelCellBag.inlines.hh \
# RelCellBag_Iterator.inlines.hh \
# RelCell_Key.inlines.hh \
# RelCell_Key_Set.inlines.hh \
# TStamp.inlines.hh
