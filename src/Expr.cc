/* Expr class implementation (non-inline functions).
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

#ifndef NOISY
#define NOISY 0
#endif

#include <config.h>

#include "util.hh"
#include "Blackboard.defs.hh"
#include "Expr.defs.hh"
#include "simplify.hh"
#include <cln/rational_io.h>
#include <cln/rational_ring.h>


namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace GiNaC {

/* x() function */

ex
x1_eval(const ex& e) {
  return x(e).hold();
}

ex
x1_evalf(const ex& e) {
  return x(e).hold();
}

ex
x1_deriv(const ex&, unsigned int) {
  abort();
}

ex
x2_eval(const ex& index, const ex& arg_list) {
  return x(index, arg_list).hold();
}

ex
x2_evalf(const ex& index, const ex& arg_list) {
  return x(index, arg_list).hold();
}

ex
x2_deriv(const ex&, const ex&, unsigned int) {
  abort();
}

  // We can't use the standard REGISTER_FUNCTION macro since we are
  // overloading x().

  unsigned x1_SERIAL::serial =
  function::register_new(function_options("x", 1).
			 eval_func(x1_eval).
			 evalf_func(x1_evalf).
			 derivative_func(x1_deriv).
			 overloaded(2));

  unsigned x2_SERIAL::serial =
  function::register_new(function_options("x", 2).
			 eval_func(x2_eval).
			 evalf_func(x2_evalf).
			 derivative_func(x2_deriv).
			 overloaded(2));



/* floor() function */

//! Evaluation of the <CODE>floor(x)</CODE>.
ex
floor_eval(const ex& x) {
  if (is_a<numeric>(x)) {
    numeric num_x = ex_to<numeric>(x);
    if (num_x.is_integer())
      return x;
    else if (num_x.is_rational()) {
      const cln::cl_I tmp = cln::floor1(cln::the<cln::cl_RA>(num_x.to_cl_N()));
      return numeric(tmp);
    }
    else
      return floor(x).hold();
  }
  else
    return floor(x).hold();
}

ex
floor_evalf(const ex& x) {
  return floor(x).hold();
}

ex
floor_deriv(const ex&, unsigned int) {
  abort();
}

/*!
  We define the symbolic function
  \f[
    floor(x).
  \f]
  It gives the largest integer less than or equal to x.
*/
REGISTER_FUNCTION(floor,
		  eval_func(floor_eval).
		  evalf_func(floor_evalf).
		  derivative_func(floor_deriv));


/* Sc() function */

ex
Sc_eval(const ex& x, const ex& y) {
  return Sc(x, y).hold();
}

ex
Sc_evalf(const ex& x, const ex& y) {
  return Sc(x, y).hold();
}

ex
Sc_deriv(const ex&, const ex&, unsigned int) {
  abort();
}

/*!
  We define the symbolic function
  \f[
    Sc(x, y) = \lfloor
                 \frac{x}{y^{\lfloor \frac{\log x}{\log y} \rfloor}}
               \rfloor.
  \f]
*/
REGISTER_FUNCTION(Sc,
		  eval_func(Sc_eval).
		  evalf_func(Sc_evalf).
		  derivative_func(Sc_deriv));


/* mod() function */

//! Evaluation of the <CODE>mod(x, y)</CODE>.
/*!
  mod_eval(const ex& x, const ex& y)

  \param x    The expression that we want divide.
  \param y    The expression for which \p x is divided.

  We apply the following properties:
  \f[
    \begin{cases}
      mod(x, y) = Number::mod(x, y),
        \quad \text{if x and y are both numerics}; \\
      mod(x - h y, y) = mod(x, y),
          \quad \text{if y and h are integers}.
    \end{cases}
  \f]
*/
ex
mod_eval(const ex& x, const ex& y) {
  if (is_a<numeric>(x) && is_a<numeric>(y))
    return mod(ex_to<numeric>(x), ex_to<numeric>(y));
  else
    if (is_a<numeric>(y)) {
      numeric k = ex_to<numeric>(y);
      if (k.is_integer() && is_a<add>(x) && x.nops() == 2) {
	if (is_a<symbol>(x.op(0)) && is_a<numeric>(x.op(1))) {
	  numeric h = ex_to<numeric>(x.op(1));
	  if (mod(h, k) == 0)
	    return mod(x.op(0), y).hold();
	}
	if (is_a<symbol>(x.op(1)) && is_a<numeric>(x.op(0))) {
	  numeric h = ex_to<numeric>(x.op(0));
	  if (mod(h, k) == 0)
	    return mod(x.op(1), y).hold();
	}
      }
    }
  return mod(x, y).hold();
}

ex
mod_evalf(const ex& x, const ex& y) {
  return mod(x, y).hold();
}

ex
mod_deriv(const ex&, const ex&, unsigned int) {
  abort();
}

/*!
  We define the general symbolic function
  \f[
    mod(n, k).
  \f]
*/
REGISTER_FUNCTION(mod,
		  eval_func(mod_eval).
		  evalf_func(mod_evalf).
		  derivative_func(mod_deriv));


/* binom() function */

//! ...
ex
binom_eval(const ex& m, const ex& k) {
  if (is_a<numeric>(m) && is_a<numeric>(k)) {
    numeric num_k = ex_to<numeric>(k);;
    if (num_k < 0 || !num_k.is_integer())
      throw std::range_error("We do not know how to evaluate\n"
			     "`binom(m, k)' with `k' not non-negative "
			     "integer.");
    else {
      // If `k == 0' then `prod(h, m - num_k + 1, m, h) / factorial(num_k)'
      // is equal to 1.
      symbol h;
      return prod(ex(h), m-num_k+1, m, ex(h)) / factorial(num_k);
    }
  }
  else
    return binom(m, k).hold();
}

ex
binom_evalf(const ex& m, const ex& k) {
  return binom(m, k).hold();
}

ex
binom_deriv(const ex&, const ex&, unsigned int) {
  abort();
}

/*!
  We define the symbolic function
  \f[
    \binom{m}{k}, \text{for } m \in \Cset, k \in \Zset, k \geq 0.
  \f]
*/
REGISTER_FUNCTION(binom,
		  eval_func(binom_eval).
		  evalf_func(binom_evalf).
		  derivative_func(binom_deriv));


/* sum() function */

/*!
  Let \f$ e(x) \f$ be the expression in \p x contained in \p e.
  This function computes two expressions \f$ e_1 \f$ and \f$ e_2 \f$
  such that \f$ e = e_1 \cdot e_2 \f$: \f$ e_1 \f$ contains all factors of
  \f$ e \f$ that do not depend from the symbol \p x; \f$ e_2 \f$
  contains all factors of \f$ e \f$ that depend from the symbol \p x.
*/
void
get_out_factors_from_argument(const ex& e, const ex& x, ex& in, ex& out) {
  if (is_a<mul>(e))
    for (unsigned int i = e.nops(); i-- > 0; ) {
      const ex& factor = e.op(i);
      if (factor.has(x))
	in *= factor;
      else
	out *= factor;
    }
  else
    if (e.has(x))
      in *= e;
    else
      out *= e;
}

//! Evaluation of the <CODE>sum(index, lower, upper, summand)</CODE>.
/*!
  sum_eval(const ex& index, const ex& lower, const ex& upper,
           const ex& summand)

  \param index    The symbol we are summing over.
  \param lower    The lower limit of the sum (we want integer lower limit).
  \param upper    The upper limit of the sum.
  \param summand  The summand.

  We explicitly allow parameters in the fourth argument, and therefore we
  have to specify the symbol we are summing over.

  We apply the following properties:
  \f[
    \begin{cases}
      \sum_{k = a}^b f(k) = 0, \quad \text{if } a > b; \\
      \sum_{k = a}^b f(k) = f(a), \quad \text{if } a = b; \\
      \sum_{k = a}^b f(k) = f(a) + f(a+1) + \cdots + f(b),
        \quad \text{if } a < b; \\
      \sum_{k = a}^b \alpha f(k) = \alpha \sum_{k = a}^b f(k),
        \quad \text{where } \alpha \text{does not depend from } k, \\
      \sum_{k = a}^b 1 = b - a + 1.
    \end{cases}
  \f]

  \exception std::invalid_argument thrown if \f$ lower \f$ in not an
                                   integer number.
  \exception std::invalid_argument thrown if \f$ upper \f$ is a number but
                                   not integer.
*/
ex
sum_eval(const ex& index, const ex& lower, const ex& upper,
	 const ex& summand) {
  if (!is_a<numeric>(lower))
    throw std::invalid_argument("The lower limit of a sum must be a number");
  numeric num_lower = ex_to<numeric>(lower);
  if (!num_lower.is_integer())
    throw std::invalid_argument("The lower limit of a sum must be an integer");
  ex s = 0;
  // `upper' is a number.
  if (is_a<numeric>(upper)) {
    numeric num_upper = ex_to<numeric>(upper);
    if (!num_upper.is_integer())
      throw std::invalid_argument("If the upper limit of a sum is a number, "
				  "it must be an integer");
    if (num_lower > num_upper)
      return 0;
    else if (num_lower == num_upper)
      return summand.subs(index == lower);
    else
      for (numeric j = num_lower; j <= num_upper; ++j)
	s += summand.subs(index == j);
  }
  // `upper' is not a number.
  else {
    // `summand' is equal to `1'.
    if (summand == 1)
      return upper - lower + 1;
    ex factors_in = 1;
    ex factors_out = 1;
    get_out_factors_from_argument(summand, index, factors_in, factors_out);
    if (factors_in == 1)
      return factors_out * (upper - lower + 1);
    else
      return factors_out * sum(index, lower, upper, factors_in).hold();
  }
  return s;
}

ex
sum_evalf(const ex& index, const ex& lower, const ex& upper,
	  const ex& summand) {
  return sum(index, lower, upper, summand).hold();
}

ex
sum_deriv(const ex&, const ex&, const ex&, const ex&, unsigned int) {
  abort();
}

/*!
  We define the general symbolic function summation
  \f[
    \sum_{k = a}^b f(k).
  \f]
*/
REGISTER_FUNCTION(sum,
		  eval_func(sum_eval).
		  evalf_func(sum_evalf).
		  derivative_func(sum_deriv));


/* prod() function */

//! Evaluation of the <CODE>prod(index, lower, upper, factor)</CODE>.
/*!
  prod_eval(const ex& index, const ex& lower, const ex& upper,
  const ex& factor)

  \param index    The symbol we are multiplying over.
  \param lower    The lower limit of the prod (we consider only rational
                  number like lower limit).
  \param upper    The upper limit of the prod.
  \param factor   The argument of the product.

  We explicitly allow parameters in the fourth argument, and therefore we
  have to specify the symbol we are multiplying over.

  We apply the following properties:
  \f[
    \begin{cases}
      \prod_{k = a}^b f(k) = 1, \quad \text{if } a > b; \\
      \prod_{k = a}^b f(k) = f(a), \quad \text{if } a = b; \\
      \prod_{k = a}^b f(k) = f(a) \cdot f(a+1) \cdots f(c),
        \quad \text{where } 0 \leq b - c \le 1, \text{if } a < b; \\
      \prod_{k = a}^b f(k)
        = \prod_{k = a}^n f(k) \cdot f(n)^(-1) \cdot f(n-1)^(-1)
	  \cdots f(n-j+1)^(-1),
        \quad \text{if } b = n + j \text{and j is a positive integer}; \\
      \prod_{k = a}^b f(k) = \prod_{k = a}^n f(k) \cdot f(n+1) \cdots f(n+j),
        \quad \text{if } b = n + j \text{and j is a negative integer}; \\
      \prod_{k = a}^b \alpha f(k) = \alpha^{b-a+1} \prod_{k = a}^b f(k),
        \quad \text{where } \alpha \text{does not depend from } k; \\
      \prod_{k = a}^b 1 = 1.
    \end{cases}
  \f]

  \exception std::invalid_argument thrown if \f$ lower \f$ in not a
                                   rational number.
  \exception std::invalid_argument thrown if \f$ upper \f$ is a number but
                                   not rational.
*/
ex
prod_eval(const ex& index, const ex& lower, const ex& upper,
	  const ex& factor) {
  if (!is_a<numeric>(lower))
    throw std::invalid_argument("The lower limit of a product "
				"must be a number");
  numeric num_lower = ex_to<numeric>(lower);
  if (!num_lower.is_rational())
    throw std::invalid_argument("The lower limit of a product "
				"must be a rational number");
  ex p = 1;
  // `upper' is a number.
  if (is_a<numeric>(upper)) {
    numeric num_upper = ex_to<numeric>(upper);
    if (!num_upper.is_rational())
      throw std::invalid_argument("If the upper limit of a product is a "
				  "number, it must be rational");
    if (num_lower > num_upper)
      return 1;
    else if (num_lower == num_upper)
      return factor.subs(index == lower);
    else
      for (numeric j = num_lower; j <= num_upper; ++j)
	p *= factor.subs(index == j);
  }
  // `upper' is not a number.
  else {
    // `factor' is equal to `1'.
    if (factor == 1)
      return 1;
    // `upper' is a sum of two addends.
    if (is_a<add>(upper) && upper.nops() == 2) {
      ex first_term = upper.op(0);
      ex second_term = upper.op(1);
      numeric numeric_term;
      symbol symbolic_term;
      if (is_a<numeric>(first_term) && is_a<symbol>(second_term)) {
 	numeric_term = ex_to<numeric>(first_term);
 	symbolic_term = ex_to<symbol>(second_term);
      }
      else if (is_a<symbol>(first_term) && is_a<numeric>(second_term)) {
	numeric_term = ex_to<numeric>(second_term);
 	symbolic_term = ex_to<symbol>(first_term);
      }
      else {
	ex factors_in = 1;
	ex factors_out = 1;
	get_out_factors_from_argument(factor, index, factors_in, factors_out);
	if (factors_in == 1)
	  return power(factors_out, upper-lower+1);
	else
	  return power(factors_out, upper-lower+1)
	    * prod(index, lower, upper, factors_in).hold();
      }
      if (numeric_term.is_integer()) {
	ex factors_in = 1;
	ex factors_out = 1;
	get_out_factors_from_argument(factor, index, factors_in, factors_out);
	if (factors_in == 1)
	  p *= power(factors_out, upper-lower+1);
	else
	  p *= power(factors_out, upper-lower+1)
	    * prod(index, lower, ex(symbolic_term), factors_in).hold();
	if (numeric_term.is_pos_integer())
	  for (numeric j = 1; j <= numeric_term; ++j)
	    p *= factor.subs(index == symbolic_term + j);
	else
	  for (numeric j = numeric_term + 1; j <= 0 ; ++j)
	    p *= 1 / factor.subs(index == symbolic_term + j);
      }
      else {
	ex factors_in = 1;
	ex factors_out = 1;
	get_out_factors_from_argument(factor, index, factors_in, factors_out);
	if (factors_in == 1)
	  return power(factors_out, upper-lower+1);
	else
	  return power(factors_out, upper-lower+1)
	    * prod(index, lower, upper, factors_in).hold();
      }
    }
    else {
      ex factors_in = 1;
      ex factors_out = 1;
      get_out_factors_from_argument(factor, index, factors_in, factors_out);
      if (factors_in == 1)
	return power(factors_out, upper-lower+1);
      else
	return power(factors_out, upper-lower+1)
	  * prod(index, lower, upper, factors_in).hold();
    }
  }
  return p;
}

ex
prod_evalf(const ex& index, const ex& lower, const ex& upper,
	   const ex& factor) {
  return prod(index, lower, upper, factor).hold();
}

ex
prod_deriv(const ex&, const ex&, const ex&, const ex&, unsigned int) {
  abort();
}

/*!
  We define the general symbolic function product
  \f[
    \prod_{k = a}^b f(k).
  \f]
*/
REGISTER_FUNCTION(prod,
		  eval_func(prod_eval).
		  evalf_func(prod_evalf).
		  derivative_func(prod_deriv));

/* max() function */

//! Evaluation of the <CODE>max(a, b)</CODE>.
/*!
  max_eval(const ex& first, const ex& second)

  \param first    The first quantity to compare.
  \param lower    The second quantity to compare.

  All symbols occurring in the two expressions are assumed
  to stand for nonnegative quantities.

*/
ex
max_eval(const ex& first, const ex& second) {
  PURRS::Expr diff(first - second);
  if (diff.is_a_number())
    if (compare(diff, 0) == 1)
      return first;
    else return second;
  else if (diff.preserves_nonnegativity())
    return first;
  else if ((-diff).preserves_nonnegativity())
    return second;
  else
    return max(first, second).hold();
}

ex
max_evalf(const ex& first, const ex& second) {
  return max(first, second).hold();
}

ex
max_deriv(const ex&, const ex&, unsigned int) {
  abort();
}

/*!
  We define the max function.
  \f[
    \max(a,b).
  \f]
*/
REGISTER_FUNCTION(max,
		  eval_func(max_eval).
		  evalf_func(max_evalf).
		  derivative_func(max_deriv));

} // namespace GiNaC

namespace {
using namespace PURRS;

//! \brief
//! Substitute every critical expression contained in \p e with
//! an arbitrary symbol.
/*!
  A <EM>critical expression</EM> is a polynomial expression that
  for GiNaC is not a polynomial. More rigorously: 
  -  every constant function with respect to the symbol contained
  in \p y, is a critical expression;
  -  every power with not numeric exponent or numeric not integer
     exponent, is a critical expression.
*/
Expr
substitute_critical_ex(const Expr& e, const Expr_List& y,
		       Blackboard& blackboard) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += substitute_critical_ex(e.op(i), y, blackboard);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten *= substitute_critical_ex(e.op(i), y, blackboard);
  }
  else if (e.is_a_power()) {
    const Expr& exponent = e.arg(1);
    Number num;
    if (!exponent.is_a_number()
	|| (exponent.is_a_number(num) && !num.is_integer()))
      return blackboard.insert_definition(e);
    else
      return pwr(substitute_critical_ex(e.arg(0), y, blackboard),
		 substitute_critical_ex(exponent, y, blackboard));
  }
  else if (e.is_a_function()) {
    for (unsigned int i = y.nops(); i-- > 0; )
      assert(e.is_a_constant_function(y.op(i).ex_to_symbol()));
    return blackboard.insert_definition(e);
  }
  else
    e_rewritten = e;
  return e_rewritten;
}

} // anonymous namespace

PURRS::Expr
PURRS::sqrfree(const Expr& x, const Expr_List& y) {
  for (unsigned int i = y.nops(); i-- > 0; )
    assert(x.is_polynomial(y.op(i).ex_to_symbol()));
  // FIXME: temporary!
  // Substitute every critical expression with an arbitrary symbol
  // and at the end resubstitute the symbol with the original value.
  // A critical expression is a polynomial expression that for GiNaC
  // is not a polynomial (i.e. log(2)).
  Blackboard blackboard;
  Expr x_subs = substitute_critical_ex(x, y, blackboard);
  x_subs = GiNaC::sqrfree(x_subs, y.l);
  return blackboard.rewrite(x_subs);
}

/*! \relates Parma_Recurrence_Realtion_Solver::Expr */
int
PURRS::compare(const Expr& x, const Expr& y) {
  // FIXME: to be improved without using `Expr::unsafe_fp_approximation()'.
  Expr diff = x - y;
  if (diff == 0)
    return 0;
  else {
    static const Number threshold = exact_pwr(10, -6); 
    diff = diff.unsafe_fp_approximation();
    assert(diff.is_a_number());
    Number numeric_diff = diff.ex_to_number();
    if (numeric_diff.is_positive() && numeric_diff > threshold)
      return 1;
    else if (numeric_diff.is_negative() && numeric_diff < -threshold)
      return -1;
    else
      return 2;
  }
}

bool
PURRS::Expr::is_a_constant_power(const Symbol& x) const {
  const Expr& e = *this;
  if (e.is_a_power() && !e.has(x))
    return true;
  else
    return false;
}

bool
PURRS::Expr::is_a_constant_function(const Symbol& x) const {
  const Expr& e = *this;
  if (e.is_a_function() && !e.has(x))
    return true;
  else
    return false;
}

PURRS::Expr
PURRS::Expr::substitute(const Expr& s, const Expr& r) const {
  const Expr& e = *this;
  Expr e_substituted;
  if (e == s)
    return r;
  else if (e.is_a_add()) {
    e_substituted = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_substituted += e.op(i).substitute(s, r);
  }
  else if (e.is_a_mul()) {
    e_substituted = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_substituted *= e.op(i).substitute(s, r);
  }
  else if (e.is_a_power())
    return pwr(e.arg(0).substitute(s, r), e.arg(1).substitute(s, r));
  else if (e.is_a_function())
    if (e.nops() == 1)
      return apply(e.functor(), e.arg(0).substitute(s, r));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int j = 0; j < num_argument; ++j)
	argument[j] = e.arg(j).substitute(s, r);
      return apply(e.functor(), argument);
    }
  else if (e.is_a_Expr_List()) {
    unsigned int num_arguments = e.nops();
    GiNaC::lst l;
    for (unsigned int i = 0; i < num_arguments; ++i)
      l.append(static_cast<const GiNaC::ex>(e.op(i).substitute(s, r)));
    e_substituted = static_cast<const GiNaC::ex>(l);
  }
  else
    return e;
  return e_substituted;
}

namespace {
using namespace PURRS;

Expr
distribute_mul_over_add_factor(const Expr& e) {
  assert(e.is_a_mul());
  Expr distributed_e = e.op(0);
  for (unsigned int i = e.nops(); i-- > 1; ) {
    Expr factor = e.op(i);
    Expr tmp = 0;
    if (factor.is_a_add())
      for (unsigned int j = factor.nops(); j-- > 0; )
	if (distributed_e.is_a_add())
	  for (unsigned int h = distributed_e.nops(); h-- > 0; )
	    tmp += factor.op(j) * distributed_e.op(h);
	else
	  tmp += factor.op(j) * distributed_e;
    else
      if (distributed_e.is_a_add())
	for (unsigned int h = distributed_e.nops(); h-- > 0; )
	  tmp += factor * distributed_e.op(h);
      else
	tmp += factor * distributed_e;
    distributed_e = tmp;
  }
  return distributed_e;
}

} // anonymous namespace

PURRS::Expr
PURRS::Expr::distribute_mul_over_add() const {
  const Expr& e = *this;
  Expr distributed_e;
  if (e.is_a_add()) {
    distributed_e = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      distributed_e += e.op(i).distribute_mul_over_add();
  }
  else if (e.is_a_mul())
    distributed_e = distribute_mul_over_add_factor(e);
  else
    distributed_e = e;
  return distributed_e;
}

bool
PURRS::Expr::is_scalar_representation(const Symbol& x) const {
  const Expr& e = *this;
  Number num;
  if (e.is_a_number(num))
    return true;
  else if (e.is_a_constant())
    return true;
  else if (e.is_a_symbol() && e != x)
    return true;
  else if (e.is_a_power())
    return e.arg(0).is_scalar_representation(x)
      && e.arg(1).is_scalar_representation(x);
  else if (e.is_a_function()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.arg(i).is_scalar_representation(x))
	return false;
    return true;
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_scalar_representation(x))
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_polynomial(const Symbol& x) const {
  const Expr& e = *this;
  if (e.is_scalar_representation(x))
    return true;
  else if (e == x)
    return true;
  else if (e.is_a_power()) {
    if (e.arg(0).is_polynomial(x)) {
      Number exponent;
      if (e.arg(1).is_a_number(exponent) && exponent.is_positive_integer()) 
	return true;
    }
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_polynomial(x))
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_integer_scalar_representation(const Symbol& x) const {
  const Expr& e = *this;
  Number num;
  if (e.is_a_number(num) && num.is_integer())
    return true;
  else if (e.is_a_symbol() && e != x)
    return true;
  else if (e.is_a_power()) {
    if (e.arg(0).is_integer_scalar_representation(x)) {
      Number exponent;
      if ((e.arg(1).is_a_number(exponent) && exponent.is_positive_integer())
	  || (e.is_a_symbol() && e != x)) 
	return true;
    }
  }
  else if (e.is_a_function()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!(e.arg(i).is_a_symbol() && e != x))
	return false;
    return true;
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_integer_scalar_representation(x))
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_integer_polynomial(const Symbol& x) const {
  const Expr& e = *this;
  if (e.is_integer_scalar_representation(x))
    return true;
  else if (e == x)
    return true;
  else if (e.is_a_power()) {
    if (e.arg(0).is_integer_polynomial(x)) {
      Number exponent;
      if (e.arg(1).is_a_number(exponent) && exponent.is_positive_integer()) 
	return true;
    }
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_integer_polynomial(x))
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_rational_scalar_representation(const Symbol& x) const {
  const Expr& e = *this;
  Number num;
  if (e.is_a_number(num) && num.is_rational())
    return true;
  else if (e.is_a_symbol() && e != x)
    return true;
  else if (e.is_a_power())
    if (e.arg(0).is_rational_scalar_representation(x)) {
      Number exponent;
      if (e.arg(1).is_a_number(exponent)) {
	if (exponent.is_positive_integer() || !e.arg(0).is_a_number())
	  return true;
      }
      else
	return true;
    }
  else if (e.is_a_function()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!(e.arg(i).is_a_symbol() && e != x))
	return false;
    return true;
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_rational_scalar_representation(x))
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_rational_polynomial(const Symbol& x) const {
  const Expr& e = *this;
  if (e.is_rational_scalar_representation(x))
    return true;
  else if (e == x)
    return true;
  else if (e.is_a_power()) {
    if (e.arg(0).is_rational_polynomial(x)) {
      Number exponent;
      if (e.arg(1).is_a_number(exponent) && exponent.is_positive_integer()) 
	return true;
    }
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_rational_polynomial(x))
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_scalar_representation() const {
  const Expr& e = *this;
  if (e.is_a_number())
    return true;
  else if (e.is_a_constant())
    return true;
  else if (e.is_a_power())
    return e.arg(0).is_scalar_representation()
      && e.arg(1).is_scalar_representation();
  else if (e.is_a_function()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.arg(i).is_scalar_representation())
	return false;
    return true;
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_scalar_representation())
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_multivariate_polynomial() const {
  const Expr& e = *this;
  if (e.is_scalar_representation())
    return true;
  else if (e.is_a_symbol())
    return true;
  else if (e.is_a_power()) {
    if (e.arg(0).is_multivariate_polynomial()) {
      Number exponent;
      if (e.arg(1).is_a_number(exponent) && exponent.is_positive_integer()) 
	return true;
    }
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!e.op(i).is_multivariate_polynomial())
	return false;
    return true;
  }
  return false;
}

bool
PURRS::Expr::is_rational_function(const Symbol& x) const {
  const Expr& e = *this;
  Expr numerator;
  Expr denominator;
  e.numerator_denominator(numerator, denominator);
  if (numerator.is_polynomial(x) && denominator.is_polynomial(x))
    return true;
  return false;
}

bool
PURRS::Expr::has_non_rational_numbers() const {
  const Expr& e = *this;
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (e.op(i).has_non_rational_numbers())
	return true;
  }
  else if (e.is_a_power()) {
    if (e.arg(0).has_non_rational_numbers()
	|| e.arg(1).has_non_rational_numbers())
      return true;
  }
  else if (e.is_a_function()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (e.arg(i).has_non_rational_numbers())
	return true;
  }
  else {
    Number num;
    if (e.is_a_number(num) && !num.is_rational())
      return true;
  }
  return false;
}

bool
PURRS::Expr::has_x_function() const {
  const Expr& e = *this;
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (e.op(i).has_x_function())
	return true;
  }
  else if (e.is_a_power()) {
    if (e.arg(0).has_x_function()
	|| e.arg(1).has_x_function())
      return true;
  }
  else if (e.is_a_function())
    if (e.is_the_x_function())
	return true;
    else
      for (unsigned int i = e.nops(); i-- > 0; )
	if (e.arg(i).has_x_function())
	  return true;
  return false;
}

bool
PURRS::Expr::has_x_function(const Expr& y) const {
  const Expr& e = *this;
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (e.op(i).has_x_function(y))
	return true;
  }
  else if (e.is_a_power()) {
    if (e.arg(0).has_x_function(y)
	|| e.arg(1).has_x_function(y))
      return true;
  }
  else if (e.is_a_function())
    if (e.is_the_x_function()) {
      if (e.arg(0).has(y))
	return true;
    }
    else
      for (unsigned int i = e.nops(); i-- > 0; )
	if (e.arg(i).has_x_function(y))
	  return true;
  return false;
}

bool
PURRS::Expr::has_sum_or_prod_function() const {
  const Expr& e = *this;
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (e.op(i).has_sum_or_prod_function())
	return true;
  }
  else if (e.is_a_power()) {
    if (e.arg(0).has_sum_or_prod_function()
	|| e.arg(1).has_sum_or_prod_function())
      return true;
  }
  else if (e.is_a_function())
    if (e.is_the_sum_function() || is_the_prod_function())
      return true;
    else
      for (unsigned int i = e.nops(); i-- > 0; )
	if (e.arg(i).has_sum_or_prod_function())
	  return true;
  return false;
}

bool
PURRS::Expr::preserves_nonnegativity() const {
  const Expr& e = *this;
  // All symbols are assumed to stand for nonnegative quantities.
  if (e.is_a_symbol())
    return true;
  else if (e.is_a_number())
    return e.ex_to_number().is_positive();
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; ) {
      const Expr& operand = e.op(i);
      if (!operand.preserves_nonnegativity()) {
	return false;
      }
    }
    return true;
  }
  else if (e.is_a_power())
    return e.arg(0).preserves_nonnegativity();
  else
    // If we were unable to prove nonnegativity, return false.
    return false;
}

Expr
PURRS::Expr::collect_term(const Expr& x, Expr& coeff_x) const {
  assert((*this).is_expanded());
  assert(!x.is_a_add() && !x.is_a_mul() && !x.is_a_power()
	 && !x.is_a_number());
  const Expr& e = simplify_ex_for_input(*this, true);
  Expr e_rewritten = 0;
  if (e.is_a_add()) {
    Expr coeffs_of_x = 0;
    for (unsigned int i = e.nops(); i-- > 0; ) {
      const Expr& addend = e.op(i);
      // The i-th addend is a product.
      if (addend.is_a_mul()) {
	bool found_x = false;
	Expr mul_for_x = 1;
	for (unsigned int i = addend.nops(); i-- > 0; ) {
	  const Expr& factor = addend.op(i);
	  if (factor == x)
	    found_x = true;
	  else if (factor.is_a_power() && factor.arg(0) == x) {
	    found_x = true;
	    mul_for_x *= pwr(x, factor.arg(1)-1);
	  }
	  else
	    mul_for_x *= addend.op(i);
	}
	if (found_x)
	  coeffs_of_x += mul_for_x;
	else
	  e_rewritten += addend;
      }
      // The i-th addend is equal to `x'.
      else if (addend == x)
	coeffs_of_x += 1;
      // The i-th addend is a power.
      else if (addend.is_a_power() && addend.arg(0) == x)
	coeffs_of_x += pwr(x, addend.arg(1)-1);
      // It is not possible to collect `x' from the i-th addend.
      else
	e_rewritten += addend;
    }
    e_rewritten += x * coeffs_of_x;
    coeff_x = coeffs_of_x;
    assert(e.expand() == e_rewritten.expand());
    return e_rewritten;
  }
  else {
    coeff_x = 1;
    return e;
  }
}

void
PURRS::Expr::collect_symbols(Symbol::SymbolSet& system_generated_symbols,
			     Symbol::SymbolSet& new_symbols) const {
  const Expr& e = *this;
  Symbol s;
  if (e.is_a_add() || e.is_a_mul())
    for (unsigned int i = e.nops(); i-- > 0; )
      e.op(i).collect_symbols(system_generated_symbols, new_symbols);
  else if (e.is_a_power()) {
    e.arg(0).collect_symbols(system_generated_symbols, new_symbols);
    e.arg(1).collect_symbols(system_generated_symbols, new_symbols);
  }
  else if (e.is_a_function())
    for (unsigned int i = e.nops(); i-- > 0; )
      e.arg(i).collect_symbols(system_generated_symbols, new_symbols);
  else if (e.is_a_symbol(s))
    if (s.is_system_generated())
      system_generated_symbols.insert(s);
    else
      new_symbols.insert(s);
}

void
PURRS::mathml_print(const GiNaC::ex & e, std::ostream& s)
{
  std::string open_tag = "";
  std::string middle_tag = "";
  std::string close_tag = "";
  bool fenced;
  if (GiNaC::is_a<GiNaC::function>(e)) {
    std::string function_name=GiNaC::ex_to<GiNaC::function>(e).get_name();
    s << "<mi>" << function_name << "</mi>";
  }
  else {
    std::string class_name=GiNaC::ex_to<GiNaC::basic>(e).class_name();
    if (class_name=="numeric") {
      if (e.info(GiNaC::info_flags::integer))
	fenced=false;
      else
	fenced=true;
      open_tag = (fenced?"<mn><mfenced>":"<mn>");
      middle_tag = "";
      close_tag = (fenced?"</mfenced></mn>":"</mn>");
      if (e.info(GiNaC::info_flags::rational) && !e.info(GiNaC::info_flags::integer)) {
	open_tag="";
	middle_tag="";
	close_tag="";
      }
    }
    else if (class_name=="symbol") {
      open_tag = "<mi>";
      middle_tag = "";
      close_tag = "</mi>";
    }
    else if (class_name=="add") {
      open_tag = "<mrow>";
      middle_tag = "<mo>+</mo>";
      close_tag = "</mrow>";
    }
    else if (class_name=="mul") {
      open_tag = "<mrow>";
      middle_tag = "<mo>*</mo>";
      close_tag = "</mrow>";
    }
    else if (class_name=="power") {
      open_tag = "<msup><mrow>";
      middle_tag = "</mrow><mrow>";
      close_tag = "</mrow></msup>";
    }
    else
      s << "UNKNOWN:" << class_name;
  }
  s << open_tag;
  size_t n = e.nops();
  if (n)
    for (size_t i=0; i<n; i++) {
      size_t child_nops=e.op(i).nops();
      if (child_nops > 1)
	s << "<mfenced>";
      mathml_print(e.op(i), s);
      if (child_nops > 1)
	s << "</mfenced>";
      if (i != n-1)
	//s << ",";
	s << middle_tag;
    }
  else if (e.info(GiNaC::info_flags::rational) && !e.info(GiNaC::info_flags::integer))
    s << "<mfrac>" << "<mn>" << e.numer() << "</mn><mn>" << e.denom() << "</mn></mfrac> ";
  else
    s << e;
  s << close_tag << std::endl;
}

void
PURRS::Expr::mathml_output(std::ostream& s) const {
  const Expr e = *this;
  s << "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">" << std::endl;
  mathml_print(e, s);
  s << "</math>" << std::endl;
}
