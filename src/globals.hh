/* Declaration of global object.
   Copyright (C) 2002 Roberto Bagnara <bagnara@cs.unipr.it>

This file is part of the Parma University's Recurrence Relation
Solver (PURRS).

The PURRS is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The PURRS is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#ifndef PURRS_globals_hh
#define PURRS_globals_hh 1

#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Expr_List.defs.hh"
#include "Number.defs.hh"
#include "Matrix.defs.hh"
#include "Constant.defs.hh"

#include <ginac/ginac.h>
#include <cln/complex.h>

#if 0
typedef GiNaC::ex GExpr;
typedef GiNaC::symbol GSymbol;
typedef GiNaC::lst GList;
typedef GiNaC::numeric GNumber;
typedef GiNaC::matrix GMatrix;
#else
typedef Parma_Recurrence_Relation_Solver::Expr GExpr;
typedef Parma_Recurrence_Relation_Solver::Symbol GSymbol;
typedef Parma_Recurrence_Relation_Solver::Expr_List GList;
typedef Parma_Recurrence_Relation_Solver::Number GNumber;
typedef Parma_Recurrence_Relation_Solver::Matrix GMatrix;
typedef Parma_Recurrence_Relation_Solver::Constant constant;
#endif

namespace GiNaC {
ex x_eval(const ex& e);
ex x_evalf(const ex& e);
ex x_deriv(const ex&, unsigned int);

DECLARE_FUNCTION_1P(x);
}

#include <Interval.h>
#include <cinterval.h>
typedef filib::Interval Interval;
typedef cinterval CInterval;

#endif
