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
#include "simplify.hh"
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
  // Find the integer `i' starting from which the function `e' is
  // well defined and non-negative.
  if (!find_domain_in_N(e, x, i))
    return false;
  D_VAR(i);
  if (i == -1)
    ++i;
  return true;
}

//! \brief
//! Returns <CODE>true</CODE> if the polynomial function in \p x \p e
//! is non-decreasing; returns <CODE>false</CODE> otherwise.
/*!
  This function checks "heuristically" if the \f$ e(x) \f$ is a non-negative,
  non-decreasing function in \f$ x \f$, where \f$ x \f$ is any symbol.
  The function works in an inductive way as follows:
  - every number is a non-decreasing function in \f$ x \f$;
  - \f$ x \f$ is a non-decreasing function in \f$ x \f$;
  - if \f$ a \f$ is a non-decreasing function in \f$ x \f$ and
    \f$ b \f$ is a positive integer, then \f$ a^b \f$ is a
    non-decreasing function in \f$ x \f$;
  - every \f$ f(x) \f$ constant function, i.e., the argument does not
    depend from \p x, is a non-decreasing function in \f$ x \f$.
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
    if (e.is_a_constant_power(x))
      return true;
    if (is_non_decreasing_poly(e.arg(0), x)) {
      Number exponent;
      if (e.arg(1).is_a_number(exponent) && exponent.is_positive_integer()) 
	return true;
    }
  }
  else if (e.is_a_function() && e.nops() == 1)
    if (e.is_a_constant_function(x))
      return true;
    else
      return false;
  else if (e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!is_non_decreasing_poly(e.op(i), x))
	return false;
    return true;
  }
  else if (e.is_a_add()) {
    for (unsigned int i = e.nops(); i-- > 0; ) {
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
    is a non-decreasing function in \f$ x \f$;
  - the function \f$ p(x)^{q(x)} = e^{q(x) log(p(x))} \f$ is a
    non-decreasing function in \f$ x \f$ if \f$ p(x) \f$ and
    \f$ q(x) \f$ are non-decreasing functions in \f$ x \f$.
*/
bool
is_non_decreasing_no_poly(const Expr& e, const Symbol& x) {
  Number k;
  if (e.is_polynomial(x)) {
    if (is_non_decreasing_poly(e, x))
      return true;
  }
  else if (e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!is_non_decreasing_no_poly(e.op(i), x))
	return false;
    return true;
  }
  else if (e.is_a_add()) {
    for (unsigned int i = e.nops(); i-- > 0; ) {
      Number num;
      if ((e.op(i).is_a_number(num) && !num.is_positive())
	  || !is_non_decreasing_no_poly(e.op(i), x))
	return false;
    }
    return true;
  }
  else if (e.is_a_power() && e.arg(0) == Recurrence::n
	   && e.arg(1).is_a_number(k) && k.is_positive())
    return true;
  else if (e.is_the_log_function())
    return true;
  else if (e.is_a_power()) {
    const Expr& base = e.arg(0);
    const Expr& exp = e.arg(1);
    if (base.is_polynomial(x) && exp.is_polynomial(x)) {
      Number num;
      if (exp.is_a_number(num)) {
	if (num.is_positive_integer())
	  return true;
      }
      else
	if (is_non_decreasing_poly(base, x)
	    && is_non_decreasing_poly(exp, x))
	  return true;
    }
  }
  return false;
}

} // anonymous namespace

//! \brief
//! Returns <CODE>true</CODE> if \f$ f = base^n coefficient \f$ is a
//! non-negative, non-decreasing function in \p x;
//! returns <CODE>false</CODE> otherwise.
bool
PURRS::
is_non_negative_non_decreasing(const Number& base, const Expr& coefficient,
			       bool poly, const Symbol& x, bool first_time,
			       Number& condition) {
  D_VAR(base);
  D_VAR(coefficient);
  // First step: checks that exponential's base is greater or equal
  // than `1'.
  if (!base.is_positive())
    return false;
  // Second step: checks the polynomial part of the exponential's
  // coefficient.
  if (poly) {
    if (!is_non_negative(coefficient, x, condition))
      return false;
    if (!is_non_decreasing_poly(coefficient, x)) {
      Expr tmp = simplify_ex_for_input
	(coefficient.substitute(Recurrence::n, Recurrence::n+1) - coefficient,
	 true);
      if (first_time) {
	if (!is_non_negative_non_decreasing(base, tmp, poly, x, false,
					    condition))
	  return false;
      }
      else
	return false;
    }
  }
  else {
    // Third step: checks the possibly non polynomial part of the
    // exponential's coefficient.
    if (!is_non_negative(coefficient, x, condition)
	|| !is_non_decreasing_no_poly(coefficient, x))
      return false;
  }
  return true;
}
