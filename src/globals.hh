
#ifndef _globals_hh
#define _globals_hh 1

#include <ginac/ginac.h>
#include <cln/complex.h>

typedef GiNaC::ex GExpr;
typedef GiNaC::symbol GSymbol;
typedef GiNaC::lst GList;
typedef GiNaC::numeric GNumber;
typedef GiNaC::matrix GMatrix;

GExpr x_eval(const GExpr& e);
GExpr x_evalf(const GExpr& e);
GExpr x_deriv(const GExpr&, unsigned int);

namespace GiNaC {
  DECLARE_FUNCTION_1P(x);
}

// #include <Interval.h>
// #include <complex>
// typedef filib::Interval Interval;
// typedef std::complex<Interval> CInterval;

#endif
