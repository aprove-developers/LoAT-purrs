/* To be written.
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

#include "compute_prod.hh"

#include "util.hh"
#include "sum_poly.hh"
#include "ep_decomp.hh"
#include "numerator_denominator.hh"
#include "factorize.hh"
#include "Number.defs.hh"
#include "Expr.defs.hh"
#include "Recurrence.defs.hh"

// TEMPORARY
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

PURRS::Expr
PURRS::compute_product(const Expr& e, const Number& lower, const Expr& upper,
		       bool is_denominator);

namespace {
using namespace PURRS;

//! \brief
//! When possible, computes \f$ \prod_{k=lower}^upper e(k) \f$
//! if \f$ e \f$ is a sum of terms, otherwise returns the symbolic product.
Expr
compute_product_on_add(const Expr& e, const Number& lower, const Expr& upper,
		       bool is_denominator) {
  Expr e_prod;
  bool e_prod_computed = false;
  if (e.is_a_add() && e.nops() == 2) {
    // `e' is of the form a + b.
    const Expr& a = e.op(0);
    const Expr& b = e.op(1);
    Number num;
    // `e' is of the form n + p with p a number.
    if ((a == Recurrence::n && b.is_a_number(num))
	|| (b == Recurrence::n && a.is_a_number(num)) && num.is_integer())
      if (num.is_positive_integer()) {
	e_prod = factorial(e) / factorial(lower + num - 1);
	e_prod_computed = true;
      }
      else
	if (lower > -num) {
	  e_prod = factorial(e) / factorial(lower + num - 1);
	  e_prod_computed = true;
	}
	else
	  if (is_denominator)
	    throw std::domain_error("Cannot compute a product at the "
				    "denominator if one of the factor "
				    "is zero");
	  else {
	    e_prod = 0;
	    e_prod_computed = true;
	  }
  }
  else if (e == 2*Recurrence::n+1) {
    e_prod = factorial(2*Recurrence::n+1) * pwr(2, -Recurrence::n)
      / factorial(Recurrence::n);
    e_prod_computed = true;
  }
  else {
    // Allows to compute `\prod_{k=lower}^upper e(k)' for function as `a*n+a*b'
    // (`a' not rational).
    Expr a = e.content(Recurrence::n);
    if (a != 1) {
      e_prod = compute_product(e.primpart(Recurrence::n), lower, upper)
	* compute_product(a, lower, upper);
      e_prod_computed = true;
    }
    // To compute numerator and denominator is useful because allows
    // to solve cases as `a/b * n + c/d': infact consider separately
    // `a*n + c*d' (that we are able to solve if `a = 1 && c/d is
    // positive integer' or `a = 2 && c*d = 1) and `b*d'.
    Expr numerator;
    Expr denominator;
    numerator_denominator_purrs(e, numerator, denominator);
    if (denominator != 1) {
      e_prod = compute_product(numerator, lower, upper)
	* pwr(compute_product(denominator, lower, upper), -1);
      e_prod_computed = true;
    }
  }
  if (!e_prod_computed) {
    Symbol h;
    e_prod = PURRS::prod(h, lower, upper, e.substitute(Recurrence::n, h));
  }
  return e_prod;
}

//! \brief
//! When possible, computes \f$ \prod_{k=lower}^upper e(k) \f$
//! if \f$ e \f$ is a power, otherwise returns the symbolic product.
Expr
compute_product_on_power(const Expr& e,
			 const Number& lower, const Expr& upper) {
  assert(e.is_a_power());
  const Expr& base_e = e.arg(0);
  const Expr& exponent_e = e.arg(1);
  Expr e_prod;
  bool e_prod_computed = false;
  if (base_e.has(Recurrence::n)) {
    Number exponent;
    if (exponent_e.is_a_number(exponent)) {
      if (exponent.is_positive_integer())
	e_prod = pwr(compute_product(base_e, lower, upper), exponent_e);
      else
	e_prod
	  = pwr(compute_product(base_e, lower, upper, true), exponent_e);
      e_prod_computed = true;
    }
  }
  // In this case `\prod_{k=lower}^upper e(k) = k^{\sum_{h=lower}^upper f(h)}'.
  else {
    std::vector<Expr> base_of_exps;
    std::vector<Expr> exp_poly_coeff;
    std::vector<Expr> exp_no_poly_coeff;
    exp_poly_decomposition(exponent_e,
			   base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
    Expr new_exponent = 0;
    // `f(h)' is a polynomial or a product of a polynomial times an
    // exponential.
    if (vector_not_all_zero(exp_poly_coeff)) {
      Symbol k("k");
      for (unsigned i = base_of_exps.size(); i-- > 0; ) {
	Expr coeff_k = exp_poly_coeff[i].substitute(Recurrence::n, k);
	new_exponent += sum_poly_times_exponentials(coeff_k, k,
						    base_of_exps[i]);
	// `sum_poly_times_exponentials' computes the sum from 0, whereas
	// we want that the sum start from `1'.
	new_exponent -= coeff_k.substitute(k, 0);
      }
      e_prod = pwr(base_e, new_exponent);
      e_prod_computed = true;
    }
    // FIXME: aggiungere anche 
    // if (vector_not_all_zero(exp_no_poly_coeff)) {...}
    // per risolvere altre sommatorie risolvibili solo con gosper.
  }
  if (!e_prod_computed) {
    Symbol h;
    e_prod = PURRS::prod(h, lower, upper, e.substitute(Recurrence::n, h));
  }
  return e_prod;
}

} // anonymous namespace

//! \brief
//! Let \f$ e(n) \f$ be an expression in the variable \f$ n \f$.
//! This function computes \f$ e!(n) \f$ defined as follows:
//! \f[
//!   e!(0) \defeq 1,
//!   \qquad
//!   e!(n) \defeq \prod_{k=lower}^upper e(k).
//! \f]
/*!
  When possible to find the closed form for \f$ \prod_{k=lower}^upper e(k) \f$,
  we compute it; when it is not possible we returns the symbolic function
  for the product.
  We observe that if also \f$ upper \f$ is a number, in particular it must be
  an integer number, than the product is always computable: so the following
  definition is applied only when \f$ upper \f$ is not a number.
  We defined inductively \f$ \prod_{k=lower}^upper e(k) \f$ as follows:
  - if \f$ e \f$ is a constant, i.e. it not contains \f$ n \f$,
    then \f$ \prod_{k=lower}^upper e(k) = e^{upper - lower + 1} \f$;
  - if \f$ e = n \f$ then
      if \f$ lower > 0 \f$ then
        \f$ \prod_{k=lower}^upper e(k) = upper! / (lower - 1)! \f$;
      else \f$ \prod_{k=lower}^upper e(k) = 0 \f$;
  - if \f$ e = n + k \f$ where \f$ k \in \Zset \f$
      if \f$ lower > -k \f$
        \f$ e_prod = e! / (lower + k - 1)! \f$;
      else \f$ \prod_{k=lower}^upper e(k) = 0 \f$;
  - if \f$ e = 2*n+1 \f$,
    then \f$ \prod_{k=lower}^upper e(k) = \frac{(2*n + 1)!}{2^n * n} \f$;
  - if \f$ e \f$ is a power there are two cases.
    We consider \f$ a \f$ and \f$ b \f$ so that \f$ e = a^b \f$, 
    - if \f$ a \f$ contains \f$ n \f$ and \f$ b \f$ is a number,
      then \f$ \prod_{k=lower}^upper e(k) = (\prod_{k=lower}^upper a(k))^b;
    - if \f$ a \f$ not contains \f$ n, i.e. \f$ a \f$ is a constant,
      then \f$ \prod_{k=lower}^upper e(k) = k^{\sum_{h=lower}^upper f(h)} \f$;
  - if \f$ e = e_1 \cdots e_m \f$, where \f$ e_i \f$,
    for \f$ i = 1, \dots, m \f$, is one of the previous case,
    then \f$ \prod_{k=lower}^upper e(k) =  \prod_{k=lower}^upper e_1(k) \cdots
    \prod_{k=lower}^upper e_m(k) \f$.

  Note that \p e must be normalized.  
*/
PURRS::Expr
PURRS::compute_product(const Expr& e, const Number& lower, const Expr& upper,
		       bool is_denominator) {
  D_MSGVAR("*** ", e);
  assert(lower.is_integer());
  if (upper.is_a_number()) {
    Number num_upper = upper.ex_to_number();
    assert(num_upper.is_integer());
    if (lower > num_upper)
      return 1;
    else if (lower == num_upper)
      return e.substitute(Recurrence::n, lower);
    else {
      Expr tmp = 1;
      for (Number i = lower; i <= num_upper; ++i)
	tmp *= e.substitute(Recurrence::n, i);
      return tmp;
    }
  }
  Expr common_factor;
  Expr rem;
  factorize(e, Recurrence::n, common_factor, rem);
  D_VAR(common_factor);
  D_VAR(rem);
  Expr e_prod;
  // `e' has been factorized: `e = common_factor * rem'.
  if (common_factor != 1) {
    e_prod = compute_product(common_factor, lower, upper);
    e_prod *= compute_product(rem, lower, upper);
  }
  else {
    if (!rem.has(Recurrence::n))
      e_prod = pwr(rem, upper - lower + 1);
    else if (rem == Recurrence::n) {
      if (lower > 0)
	e_prod = factorial(upper) / factorial(lower - 1);
      else
	if (is_denominator)
	  throw std::domain_error("Cannot compute a product at the "
				  "denominator if one of the factor "
				  "is zero");
	else
	  e_prod = 0;
    }
    else if (rem.is_a_add())
      e_prod = compute_product_on_add(rem, lower, upper, is_denominator);
    else if (rem.is_a_power())
      e_prod = compute_product_on_power(rem, lower, upper);
    else if (rem.is_a_mul()) {
      e_prod = 1;
      for (unsigned i = rem.nops(); i-- > 0; )
	e_prod *= compute_product(e.op(i), lower, upper);
    }
    else {
      Symbol h;
      e_prod = PURRS::prod(h, lower, upper, rem.substitute(Recurrence::n, h));
    }
  }
  D_MSGVAR("*** ", e_prod);
  return e_prod;
}
