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

#include <config.h>

#ifndef NOISY
#define NOISY 0
#endif

#include "functional_equation.hh"

#include "util.hh"
#include "numerator_denominator.hh"
#include "ep_decomp.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Recurrence.defs.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

//! \brief;
//! Returns <CODE>true</CODE> if \p e, valuated in \f$ i \f$, is a
//! non-negative number; returns <CODE>false</CODE> otherwise.
bool
is_non_negative(const Expr& e, const Symbol& x, Number& i) {
  largest_positive_int_zero(numerator(e), i);
  largest_positive_int_zero(denominator(e), i);
  D_VAR(i);
  if (i == -1)
    ++i;
  Number num;
  if (numerator(e).substitute(x, i).is_a_number(num) && !num.is_negative()) {
    D_MSG("true");
    return true;
  }
  else {
    D_MSG("false");
    return false;
  }
}

//! \brief
//! Returns <CODE>true</CODE> if the polynomial function in \p x \p e
//! is non-decreasing; returns <CODE>false</CODE> otherwise.
/*!
  This function checks "heuristically" if the \f$ e(x) \f$ is a non-negative,
  non-decreasing function in \f$ x \f$, where \f$ x \f$ is any symbol.
  The function works in an inductive way as follows:
  - every number or constant non-decreasing function in \f$ x \f$;
  - \f$ x \f$ is a non-decreasing function in \f$ x \f$;
  - if \f$ a \f$ is a non-decreasing function in \f$ x \f$ and
    \f$ b \f$ is a positive integer, then \f$ a^b \f$ is a
    non-decreasing function in \f$ x \f$;
  - if \f$ a \f$ and \f$ b \f$ are non-decreasing functions
    in \f$ x \f$ then \f$ a + b \f$ and \f$ a * b \f$ are
    non-decreasing functions in \f$ x \f$.
*/
bool
is_non_decreasing_poly(const Expr& e, const Symbol& x) {
  if (e.is_a_number() || e.is_a_constant())
    return true;
  else if (e == x)
    return true;
  else if (e.is_a_power()) {
    if (is_non_decreasing_poly(e.arg(0), x)) {
      Number exponent;
      if (e.arg(1).is_a_number(exponent) && exponent.is_positive_integer()) 
	return true;
    }
  }
  else if (e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (!is_non_decreasing_poly(e.op(i), x))
	return false;
    return true;
  }
  else if (e.is_a_add()) {
    for (unsigned i = e.nops(); i-- > 0; ) {
      Number num;
      if ((e.op(i).is_a_number(num) && !num.is_positive())
	  || !is_non_decreasing_poly(e.op(i), x))
	return false;
    }
    return true;
  }
  return false;
}

//! \brief
//! Returns <CODE>true</CODE> if the non polynomial function in \p x \p e
//! is non-decreasing; returns <CODE>false</CODE> otherwise.
/*!
  We consider separately the mathematical functions:
  - a logarithm with a base grater or equal than \f$ 1 \f$
    is a non-decreasing function in \f$ x \f$.
*/
bool
is_non_decreasing_no_poly(const Expr& e, const Symbol& x) {
  if (e.is_polynomial(x)) {
    if (is_non_decreasing_poly(e, x))
      return true;
  }
  else if (e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (!is_non_decreasing_no_poly(e.op(i), x))
	return false;
    return true;
  }
  else if (e.is_a_add()) {
    for (unsigned i = e.nops(); i-- > 0; ) {
      Number num;
      if ((e.op(i).is_a_number(num) && !num.is_positive())
	  || !is_non_decreasing_no_poly(e.op(i), x))
	return false;
    }
    return true;;
  }
  else if (e.is_the_log_function())
    return true;
  return false;
}

} // anonymous namespace

//! \brief
//! Returns <CODE>true</CODE> if \p f is a non-negative, non-decreasing
//! function in \p x; returns <CODE>false</CODE> otherwise.
bool
PURRS::is_non_negative_non_decreasing(const Expr& f, const Symbol& x,
				      Number& condition) {
  D_VAR(f);
  // We search exponentials in `n' (for this the expression
  // `f' must be expanded).
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `f'.
  // In the `i'-th position of the vectors
  // `exp_poly_coeff' and `exp_no_poly_coeff' there are respectively
  // the polynomial part and possibly non polynomial part of the coefficient
  // of the exponential with the base in `i'-th position of `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of base_of_exps[i]^n.
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(f.expand(), base_of_exps,
			 exp_poly_coeff, exp_no_poly_coeff);
  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  for (unsigned i = base_of_exps.size(); i-- > 0; ) {
    // First step: checks that all exponentials'bases are greater or equal
    // than `1'.
    Number num_base;
    if (!base_of_exps[i].is_a_number(num_base)
	|| !num_base.is_positive())
      return false;
    // Second step: checks the polynomial part of the exponential's
    // coefficient.
    if (!exp_poly_coeff[i].is_zero())
      if (!is_non_negative(exp_poly_coeff[i], x, condition))
	if (!is_non_decreasing_poly(exp_poly_coeff[i], x)
	    && !is_non_negative_non_decreasing
	    (f.substitute(Recurrence::n, Recurrence::n+1)-f, x, condition))
	  return false;
    // Third step: checks the possibly non polynomial part of the
    // exponential's coefficient.
    D_MSG("STEP 3");
    if (!exp_no_poly_coeff[i].is_zero())
      if (!is_non_negative(exp_no_poly_coeff[i], x, condition)
	  || !is_non_decreasing_no_poly(exp_no_poly_coeff[i], x))
	return false;
  }
  return true;
}
