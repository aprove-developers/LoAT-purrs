
#include "globals.hh"

using namespace GiNaC;
using namespace Parma_Recurrence_Relation_Solver;

GExpr
x_eval(const GExpr& e) {
  return x(e).hold();
}

GExpr
x_evalf(const GExpr& e) {
  return x(e).hold();
}

GExpr
x_deriv(const GExpr&, unsigned int) {
  abort();
}

namespace GiNaC {
  REGISTER_FUNCTION(x,
                    eval_func(x_eval).
                    evalf_func(x_evalf).
                    derivative_func(x_deriv));
} 
