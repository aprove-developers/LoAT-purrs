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


#include "approx_expr.hh"
#include "Number.defs.hh"
#include "Complex_Interval.defs.hh"
#include "cimath.h"
#include <cln/rational.h>
#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

Interval
approximate_integer(const Number& n) {
  // Kludge!!!
  return Interval(n.to_int());
}

Interval
approximate_rational(const Number& n) {
  if (n.is_integer())
    return approximate_integer(n);
  else if (n.is_rational())
    return
      approximate_integer(n.numerator())
      / approximate_integer(n.denominator());
  else
    abort();
}

CInterval
approximate(const Number& n) {
  if (n.is_real())
    return CInterval(approximate_rational(n.real()),
		     0);
  else
    return CInterval(approximate_rational(n.real()),
		     approximate_rational(n.imaginary()));
}

#if 0
CInterval
approximate(const Expr& e) {
  CInterval r;
  if (e.is_a_number())
    return approximate(e.ex_to_number());
  else if (e.is_a_add()) {
    r = CInterval(0, 0);
    for (unsigned i = 0, n = e.nops(); i < n; ++i)
      r += approximate(e.op(i));
  }
  else if (e.is_a_mul()) {
    r = CInterval(1, 0);
    for (unsigned i = 0, n = e.nops(); i < n; ++i)
      r *= approximate(e.op(i));
  }
  else if (e.is_a_power()) {
    static Expr one_half = Number(1)/2;
    const Expr& base = e.op(0);
    const Expr& exponent = e.op(1);
#if 0
    return pow(approximate(base), approximate(exponent));
#else
    std::cout << base << "^" << exponent << std::endl;
    CInterval abase = approximate(base);
    CInterval aexponent = approximate(exponent);
    std::cout << abase << "^" << aexponent << " = ";
    CInterval result = pow(abase, aexponent);
    std::cout << result << std::endl;
    return result;
#endif
  }
  else if (e.is_a_function()) {
    const Expr& arg = e.op(0);
    CInterval aarg = approximate(arg);
#if 0
    if (e.is_the_function(abs))
      return abs(aarg);
    else
#endif
    if (e.is_the_exp_function())
      return exp(aarg);
    else if (e.is_the_log_function())
      return ln(aarg);
    else if (e.is_the_sin_function())
      return sin(aarg);
    else if (e.is_the_cos_function())
      return cos(aarg);
    else if (e.is_the_tan_function())
      return tan(aarg);
    else if (e.is_the_acos_function())
      return acos(aarg);
    else
      abort();
      abort();
  }
  else if (e.is_a_constant()) {
    if (e == Constant::Pi)
      return CInterval(Interval::PI(), Interval::ZERO());
#if 0
    else if (e == Constant::Euler)
      return CInterval(Interval(E_lower_bound.d, E_upper_bound.d),
                       Interval::ZERO());
#endif
    else
      abort();
  }
  else
    abort();

  return r;
}
#endif

bool
approximate(const Expr& e, Expr& ae, CInterval& aci) {
  static Expr operand_ae;
  static CInterval operand_aci;
  bool interval_result = true;
  if (e.is_a_number())
    aci = approximate(e.ex_to_number());
#if 0
  else if (e.is_a_complex_interval())
    aci = e.ex_to_complex_interval();
#endif
  else if (e.is_a_add()) {
    Expr accumulated_ae = 0;
    CInterval accumulated_aci(Interval::ZERO(), Interval::ZERO());
    bool non_trivial_interval = false;
    for (unsigned i = e.nops(); i-- > 0; ) {
      if (approximate(e.op(i), operand_ae, operand_aci)) {
	accumulated_aci += operand_aci;
	non_trivial_interval = true;
      }
      else {
	accumulated_ae *= operand_ae;
	interval_result = false;
      }
    }
    if (interval_result)
      aci = accumulated_aci;
    else if (non_trivial_interval)
      ae = Complex_Interval(accumulated_aci) + accumulated_ae;
    else
      ae = accumulated_ae;
  }
  else if (e.is_a_mul()) {
    Expr accumulated_ae = 1;
    CInterval accumulated_aci(Interval::ONE(), Interval::ZERO());
    bool non_trivial_interval = false;
    for (unsigned i = e.nops(); i-- > 0; ) {
      if (approximate(e.op(i), operand_ae, operand_aci)) {
	accumulated_aci *= operand_aci;
	non_trivial_interval = true;
      }
      else {
	accumulated_ae *= operand_ae;
	interval_result = false;
      }
    }
    if (interval_result)
      aci = accumulated_aci;
    else if (non_trivial_interval)
      ae = Complex_Interval(accumulated_aci) * accumulated_ae;
    else
      ae = accumulated_ae;
  }
  else if (e.is_a_power()) {
    Expr base_ae;
    CInterval base_aci;
    bool base_interval_result = approximate(e.arg(0),
					    base_ae, base_aci);
    Expr exponent_ae;
    CInterval exponent_aci;
    bool exponent_interval_result = approximate(e.arg(1),
						exponent_ae, exponent_aci);
    if (base_interval_result)
      if (exponent_interval_result) 
	aci = pow(base_aci, exponent_aci);
      else
	ae = pwr(Complex_Interval(base_aci), exponent_ae);
    else
      if (exponent_interval_result) 
	ae = pwr(base_ae, Complex_Interval(exponent_aci));
      else
	ae = pwr(base_ae, exponent_ae);
    if (!base_interval_result || !exponent_interval_result)
      interval_result = false;
  }
  else if (e.is_a_function()) {
    switch (e.nops()) {
    case 1:
      {
	if (approximate(e.arg(0), operand_ae, operand_aci))
	  if (e.is_the_exp_function())
	    aci = exp(operand_aci);
	  else if (e.is_the_log_function())
	    aci = ln(operand_aci);
	  else if (e.is_the_sin_function())
	    aci = sin(operand_aci);
	  else if (e.is_the_cos_function())
	    aci = cos(operand_aci);
	  else if (e.is_the_tan_function())
	    aci = tan(operand_aci);
	  else if (e.is_the_acos_function())
	    aci = acos(operand_aci);
	  else
	    abort();
	else
	  ae = apply(e.functor(), operand_ae);
      }
      break;
    case 3:
      assert(e.is_the_sum_function() || e.is_the_prod_function());
      ae = e;
      interval_result = false;
      break;
    default:
      abort();
    }
  }
  else if (e.is_a_constant()) {
    if (e == Constant::Pi)
      aci = CInterval(Interval::PI(), Interval::ZERO());
    else if (e == Constant::Euler)
      abort();
    else
      abort();
  }
  else
    abort();

  return interval_result;
}

Expr
approximate(const Expr& e) {
  Expr ae;
  CInterval aci;
  if (approximate(e, ae, aci))
    return Complex_Interval(aci);
  else
    return ae;
}

} // namespace Parma_Recurrence_Relation_Solver
