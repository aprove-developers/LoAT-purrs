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

#include <config.h>

#include "Expr.defs.hh"

#include "numerator_denominator.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace GiNaC {

//! ...
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


//! Evaluation of the <CODE>sum(index, lower, upper, summand)</CODE>.
/*!
  \fn GiNaC::ex
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
      \sum_{k = a}^b f(k) =
        \sum_{k = a}^n f(k) - f(n) - f(n-1) - \cdots - f(n-j+1),
        \quad \text{if } b = n + j \text{and j is a positive integer}; \\
      \sum_{k = a}^b f(k) = \sum_{k = a}^n f(k) + f(n+1) + \cdots + f(n+j),
        \quad \text{if } b = n + j \text{and j is a negative integer}.
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
  if (is_a<numeric>(upper)) {
    numeric num_upper = ex_to<numeric>(upper);
    if (!num_upper.is_integer())
      throw std::invalid_argument("If the upper limit of a sum is a number,"
				  "it must be an integer");
    if (num_lower > num_upper)
      return 0;
    else if (num_lower == num_upper)
      return summand.subs(index == lower);
    else
      for (numeric j = num_lower; j <= num_upper; ++j)
	s += summand.subs(index == j);
  }
  else {
    ex upper_limit = wild(0) + wild(1);
    lst substitution;
    if (upper.match(upper_limit, substitution)) {
      ex first_term = substitution.op(0).rhs();
      ex second_term = substitution.op(1).rhs();
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
      else
	return sum(index, lower, upper, summand).hold();
      if (numeric_term.is_integer()) {
	s += sum(index, lower, ex(symbolic_term), summand).hold();
	if (numeric_term.is_pos_integer())
	  for (numeric j = 1; j <= numeric_term; ++j)
	    s += summand.subs(index == symbolic_term + j);
	else
	  for (numeric j = numeric_term + 1; j <= 0 ; ++j)
	    s -= summand.subs(index == symbolic_term + j);
      }
      else
	return sum(index, lower, upper, summand).hold();
    }
    else
      return sum(index, lower, upper, summand).hold();
  }
  return s;
}

ex
sum_evalf(const ex&, const ex&, const ex&, const ex&) {
  abort();
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


//! Evaluation of the <CODE>prod(index, lower, upper, factor)</CODE>.
/*!
  \fn GiNaC::ex
  prod_eval(const ex& index, const ex& lower, const ex& upper,
  const ex& factor)

  \param index    The symbol we are multiplying over.
  \param lower    The lower limit of the prod (we want integer lower limit).
  \param upper    The upper limit of the prod.
  \param factor   The argument of the product.

  We explicitly allow parameters in the fourth argument, and therefore we
  have to specify the symbol we are multiplying over.

  We apply the following properties:
  \f[
    \begin{cases}
      \prod_{k = a}^b f(k) = 0, \quad \text{if } a > b; \\
      \prod_{k = a}^b f(k) = f(a), \quad \text{if } a = b; \\
      \prod_{k = a}^b f(k) = f(a) \cdot f(a+1) \cdots f(b),
        \quad \text{if } a < b; \\
      \prod_{k = a}^b f(k)
        = \prod_{k = a}^n f(k) \cdot f(n)^(-1) \cdot f(n-1)^(-1)
	  \cdots f(n-j+1)^(-1),
        \quad \text{if } b = n + j \text{and j is a positive integer}; \\
      \prod_{k = a}^b f(k) = \prod_{k = a}^n f(k) \cdot f(n+1) \cdots f(n+j),
        \quad \text{if } b = n + j \text{and j is a negative integer}.
    \end{cases}
  \f]

  \exception std::invalid_argument thrown if \f$ lower \f$ in not an
                                   integer number.
  \exception std::invalid_argument thrown if \f$ upper \f$ is a number but
                                   not integer.
*/
ex
prod_eval(const ex& index, const ex& lower, const ex& upper,
	  const ex& factor) {
  if (!is_a<numeric>(lower))
    throw std::invalid_argument("The lower limit of a product"
				"must be a number");
  numeric num_lower = ex_to<numeric>(lower);
  if (!num_lower.is_integer())
    throw std::invalid_argument("The lower limit of a product"
				"must be an integer");
  ex p = 1;
  if (is_a<numeric>(upper)) {
    numeric num_upper = ex_to<numeric>(upper);
    if (!num_upper.is_integer())
      throw std::invalid_argument("If the upper limit of a product is a"
				  "number, it must be an integer");
    if (num_lower > num_upper)
      return 1;
    else if (num_lower == num_upper)
      return factor.subs(index == lower);
    else
      for (numeric j = num_lower; j <= num_upper; ++j)
	p *= factor.subs(index == j);
  }
  else {
    ex upper_limit = wild(0) + wild(1);
    lst substitution;
    if (upper.match(upper_limit, substitution)) {
      ex first_term = substitution.op(0).rhs();
      ex second_term = substitution.op(1).rhs();
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
      else
	return prod(index, lower, upper, factor).hold();
      if (numeric_term.is_integer()) {
	p *= prod(index, lower, ex(symbolic_term), factor).hold();
	if (numeric_term.is_pos_integer())
	  for (numeric j = 1; j <= numeric_term; ++j)
	    p *= factor.subs(index == symbolic_term + j);
	else
	  for (numeric j = numeric_term + 1; j <= 0 ; ++j)
	    p *= 1 / factor.subs(index == symbolic_term + j);
      }
      else
	return prod(index, lower, upper, factor).hold();
    }
    else
      return prod(index, lower, upper, factor).hold();
  }
  return p;
}

ex
prod_evalf(const ex&, const ex&, const ex&, const ex&) {
  abort();
}

ex
prod_deriv(const ex&, const ex&, const ex&, const ex&, unsigned int) {
  abort();
}

/*!
  We define the general symbolic function summation
  \f[
    \prod_{k = a}^b f(k).
  \f]
*/
REGISTER_FUNCTION(prod,
		  eval_func(prod_eval).
		  evalf_func(prod_evalf).
		  derivative_func(prod_deriv));

} // namespace GiNaC

static PURRS::Expr
distribute_mul_over_add_factor(const PURRS::Expr& e) {
  assert(e.is_a_mul());
  PURRS::Expr distributed_e = e.op(0);
  for (unsigned i = e.nops(); i-- > 1; ) {
    PURRS::Expr factor = e.op(i);
    PURRS::Expr tmp = 0;
    if (factor.is_a_add())
      for (unsigned j = factor.nops(); j-- > 0; )
	if (distributed_e.is_a_add())
	  for (unsigned h = distributed_e.nops(); h-- > 0; )
	    tmp += factor.op(j) * distributed_e.op(h);
	else
	  tmp += factor.op(j) * distributed_e;
    else
      if (distributed_e.is_a_add())
	for (unsigned h = distributed_e.nops(); h-- > 0; )
	  tmp += factor * distributed_e.op(h);
      else
	tmp += factor * distributed_e;
    distributed_e = tmp;
  }
  return distributed_e;
}

PURRS::Expr
PURRS::Expr::distribute_mul_over_add() const {
  const Expr& e = *this;
  Expr distributed_e;
  if (e.is_a_add()) {
    distributed_e = 0;
    for (unsigned i = e.nops(); i-- > 0; )
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
  if (e.is_a_number())
    return true;
  else if (e.is_a_constant())
    return true;
  else if (e.is_a_symbol() && e != x)
    return true;
  else if (e.is_a_power())
    return e.arg(0).is_scalar_representation(x)
      && e.arg(1).is_scalar_representation(x);
  else if (e.is_a_function()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (!e.arg(i).is_scalar_representation(x))
	return false;
    return true;
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
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
    for (unsigned i = e.nops(); i-- > 0; )
      if (!e.op(i).is_polynomial(x))
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
  numerator_denominator_purrs(e, numerator, denominator);
  if (numerator.is_polynomial(x) && denominator.is_polynomial(x))
    return true;
  return false;
}

namespace {

struct Norm_1 {
  unsigned size_norm(const PURRS::Symbol&) const {
    return 1;
  }
};

}

unsigned
PURRS::size_norm(const Expr& e) {
  return generic_size_norm(e, Norm_1());
}
