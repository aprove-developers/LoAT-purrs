/* Inline functions implementing the approximation of expressions.
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

#ifndef PURRS_approximate_impl_hh
#define PURRS_approximate_impl_hh 1

#include "Expr.defs.hh"
#include "Number.defs.hh"
#include "complint.hh"
#include <cimath.h>
#include <cln/rational.h>
#include <cmath>

namespace Parma_Recurrence_Relation_Solver {

inline Interval
approximate_integer(const Number& n) {
  // Kludge!!!
  try {
    return Interval(n.to_int());
  }
  catch (...) {
    return Interval(n.to_double());
  }
}

inline Interval
approximate_rational(const Number& n) {
  if (n.is_integer())
    return approximate_integer(n);
  else if (n.is_rational())
    return
      approximate_integer(n.numerator())
      / approximate_integer(n.denominator());
  else
    return Interval(n.to_double());
}

inline CInterval
approximate(const Number& n) {
  if (n.is_real())
    return CInterval(approximate_rational(n.real()),
		     Interval::ZERO());
  else
    return CInterval(approximate_rational(n.real()),
		     approximate_rational(n.imaginary()));
}

// Kludge to get around what is likely to be a CoStLy bug: see
// http://www.cs.unipr.it/pipermail/purrs-devel/2002-November/000724.html
inline CInterval
mypow(const CInterval& base, const CInterval& exponent) {
  if (exponent.im() == Interval::ZERO()) {
    const Interval& exponent_re = exponent.re();
    if (exponent_re.sup() < 0.0)
      if (exponent_re == -Interval::ONE())
	return CInterval(Interval::ONE(), Interval::ZERO()) / base;
      else
	return CInterval(Interval::ONE(), Interval::ZERO())
	  / mypow(base, -exponent);
    else if (exponent_re.isPoint()) {
      double ep = exponent_re.sup();
      if (ep == 2.0)
	return sqr(base);
      else if (ep == 0.5)
	return sqrt(base);
      else if (ep == floor(ep) && base.re().sup() < 0.0) {
	CInterval r = pow(-base, exponent);
	if ((long(ep) % 2) != 0)
	  r = -r;
	return r;
      }
    }
  }
  return pow(base, exponent);
}

template <typename SymbolHandler>
bool
generic_approximate(const Expr& e, const SymbolHandler& sh,
		    Expr& ae, CInterval& aci) {
#if 0
  static unsigned indent = 0;
  for (unsigned int i = 0; i < indent; ++i)
    std::cout << ' ';
  std::cout << "approx " << e << std::endl;
  ++indent;
#endif
  static Expr operand_ae;
  static CInterval operand_aci;
  bool interval_result = true;
  if (e.is_a_number())
    aci = approximate(e.ex_to_number());
  else if (e.is_a_complex_interval())
    aci = e.ex_to_complex_interval().get_interval();
  else if (e.is_a_add()) {
    Expr accumulated_ae = 0;
    CInterval accumulated_aci(Interval::ZERO(), Interval::ZERO());
    bool non_trivial_interval = false;
    for (unsigned i = e.nops(); i-- > 0; ) {
      if (generic_approximate(e.op(i), sh, operand_ae, operand_aci)) {
	accumulated_aci += operand_aci;
	non_trivial_interval = true;
      }
      else {
	accumulated_ae += operand_ae;
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
      if (generic_approximate(e.op(i), sh, operand_ae, operand_aci)) {
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
    bool base_interval_result
      = generic_approximate(e.arg(0), sh, base_ae, base_aci);
    Expr exponent_ae;
    CInterval exponent_aci;
    bool exponent_interval_result
      = generic_approximate(e.arg(1), sh, exponent_ae, exponent_aci);
    if (base_interval_result)
      if (exponent_interval_result) {
#if 1
	aci = mypow(base_aci, exponent_aci);
#else
	std::cout << "Trying " << base_aci << " ^ " << exponent_aci << std::endl;
	aci = pow(base_aci, exponent_aci);
#endif
      }
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
	if (e.is_the_x_function()
	    || e.is_the_factorial_function()) {
	  ae = e;
	  interval_result = false;
	}
	else if (generic_approximate(e.arg(0), sh, operand_ae, operand_aci))
	  if (e.is_the_exp_function())
	    aci = exp(operand_aci);
	  else if (e.is_the_log_function())
	    aci = log(operand_aci);
	  else if (e.is_the_sin_function())
	    aci = sin(operand_aci);
	  else if (e.is_the_cos_function())
	    aci = cos(operand_aci);
	  else if (e.is_the_tan_function())
	    aci = tan(operand_aci);
	  else if (e.is_the_acos_function())
	    aci = acos(operand_aci);
	  else {
	    ae = apply(e.functor(), Complex_Interval(operand_aci));
	    interval_result = false;
	  }
	else {
	  ae = apply(e.functor(), operand_ae);
	  interval_result = false;
	}
      }
      break;
    case 2:
      assert(e.is_the_mod_function());
      ae = e;
      interval_result = false;
      break;
    case 4:
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
  else if (e.is_a_symbol()) {
    if (sh.approximate(e.ex_to_symbol(), operand_ae, operand_aci))
      aci = operand_aci;
    else {
      ae = operand_ae;
      interval_result = false;
    }
  }
  else
    abort();

#if 0
  --indent;
  for (unsigned int i = 0; i < indent; ++i)
    std::cout << ' ';
  if (interval_result)
    std::cout << "gives " << aci << std::endl;
  else
    std::cout << "gives " << ae << std::endl;
#endif

  return interval_result;
}

template <typename SymbolHandler>
Expr
generic_approximate(const Expr& e, const SymbolHandler& sh) {
  Expr ae;
  CInterval aci;
  if (generic_approximate(e, sh, ae, aci))
    return Complex_Interval(aci);
  else
    return ae;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_approximate_impl_hh)
