
#include "globals.hh"

namespace GiNaC {
ex
x_eval(const ex& e) {
  return x(e).hold();
}

ex
x_evalf(const ex& e) {
  return x(e).hold();
}

ex
x_deriv(const ex&, unsigned int) {
  abort();
}

REGISTER_FUNCTION(x,
		  eval_func(x_eval).
		  evalf_func(x_evalf).
		  derivative_func(x_deriv));
}
