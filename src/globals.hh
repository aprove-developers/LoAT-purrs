
#ifndef PURRS_globals_hh
#define PURRS_globals_hh 1

#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Expr_List.defs.hh"
#include "Number.defs.hh"
#include "Matrix.defs.hh"

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
#endif

GExpr x_eval(const GExpr& e);
GExpr x_evalf(const GExpr& e);
GExpr x_deriv(const GExpr&, unsigned int);

namespace GiNaC {
  DECLARE_FUNCTION_1P(x);
}

#include <Interval.h>
#include <cinterval.h>
typedef filib::Interval Interval;
typedef cinterval CInterval;

#endif
