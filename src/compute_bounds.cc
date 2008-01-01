/* Recurrence class implementation: methods (and associated auxiliary
   functions) associated to the computation of lower and upper bounds.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#ifndef NOISY
#define NOISY 0
#endif

#include "util.hh"
#include "simplify.hh"
#include "sum_poly.hh"
#include "ep_decomp.hh"
#include "gosper.hh"
#include "functional_equation.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Cached_Expr.defs.hh"
#include "Recurrence.defs.hh"
#include "Recurrence.inlines.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

//! Kinds of bounds for the solution to compute.
enum Type_Bound {
  /*! An upper bound. */
  UPPER_BOUND,
    
  /*! A lower bound. */
  LOWER_BOUND
};

bool
compute_bounds_for_exp_function(Type_Bound type, const Number& coeff,
				const Number& divisor, const Number& base,
				Expr& bound) {
  
  // Lower bound.
  if (type == LOWER_BOUND)
    bound = pwr(base, Recurrence::n);
  // Upper bound.
  else {
    Expr constant;
    int cmp = compare(pwr(base, pwr(divisor, 2)), coeff * pwr(base, divisor));
    if (cmp == 1) {
      const Expr& tmp = pwr(base, exact_pwr(divisor, 2));
      constant = coeff * tmp * pwr(tmp - coeff * pwr(base, divisor), -1);
    }
    else if (cmp == 0 || cmp == -1)
      constant = coeff * exact_pwr(coeff - 1, -1) * pwr(base, divisor)
	* pwr(log(base), -1 - log(coeff) * pwr(log(divisor), -1))
	* gamma(1 + log(coeff) * pwr(log(divisor), -1));
    else
      // cmp == 2.
      return false;
    bound = pwr(base, Recurrence::n)
      + constant * pwr(base, Recurrence::n / divisor);
  }
  return true;
}

bool
compute_bounds_for_power_of_n(Type_Bound type,
			      const Number& coeff, const Number& divisor,
			      const Number& k, Expr& bound) {
  assert(!k.is_negative());
  const Expr& pow_div_k = pwr(divisor, k);
  int diff_coeff_pow_div_k = compare(coeff, pow_div_k);
  const Expr& pow_n_k = pwr(Recurrence::n, k);
  const Expr& frac_log = log(Recurrence::n) / log(divisor);
  // Lower bound.
  if (type == LOWER_BOUND) {
    // FIXME: for the moment use only the theorem in order to compute
    // lower bound in the case `g(n) = n^k'. 
    Expr mu_b;
    Expr c_b;
    if (divisor.is_positive_integer()) {
      mu_b = 1;
      c_b = 1;
    }
    else {
      mu_b = log(2 * divisor / (divisor - 1)) / log(divisor);
      c_b = divisor / (divisor - 1);
    }
    const Expr& pow_frac_log = pwr(Recurrence::n, log(coeff) / log(divisor));
    // k == 0.
    if (k == 0) {
      Expr lambda_b;
      if (divisor >= 2)
	lambda_b = 1;
      else
	// 1 < divisor < 2.
	lambda_b = log((2 * divisor - 1)/(divisor - 1)) / log(divisor);
      if (coeff == 1)
	bound = frac_log - lambda_b;
      else
	bound = (1 - pwr(Recurrence::n, log(coeff) / log(divisor))
		 * pwr(coeff, -lambda_b)) / (1 - coeff);
    }
    // k >= 1.
    else if (k >= 1) {
      const Expr& pow_div_k_minus = pwr(divisor, k - 1);
      const Expr& pow_n_k_minus = pwr(Recurrence::n, k - 1);
      if (compare(coeff, pow_div_k_minus) == -1)
	bound = pow_div_k / (pow_div_k - coeff)
	  * (pow_n_k - pow_frac_log * pwr(pow_div_k / coeff, mu_b));
      else if (compare(coeff, pow_div_k_minus) == 0)
	bound = divisor / (divisor - 1) * pow_n_k
	  - (pwr(divisor, 1 + mu_b) / (divisor - 1) + c_b * k * frac_log - 1)
	  * pow_n_k_minus;
      else if (diff_coeff_pow_div_k == -1)
	bound = pow_div_k / (pow_div_k - coeff) * pow_n_k
	  - (pow_div_k / (pow_div_k - coeff) * pwr(pow_div_k / coeff, mu_b)
	     + c_b * k * pow_div_k_minus / (coeff - pow_div_k_minus))
	  * pow_frac_log
	  + c_b * k * coeff * pow_n_k_minus / (coeff - pow_div_k_minus);
      else if (diff_coeff_pow_div_k == 0)
	bound = pow_n_k * (frac_log - mu_b - c_b * k / (divisor - 1))
	  + c_b * k * pow_n_k_minus * divisor / (divisor - 1);
      else if (diff_coeff_pow_div_k == 1) {
	// FIXME: theorem 3.8 (alpha > beta^k).
#if 0
	bound = pow_frac_log
	  * (pow_div_k / (coeff - pow_div_k) * pwr(pow_div_k / coeff, mu_b)
	     - c_b * k * pow_div_k_minus / (coeff - pow_div_k_minus))
	  - pow_div_k / (coeff - pow_div_k) * pow_n_k
	  + c_b * k * coeff * pow_n_k_minus
	  / (coeff - pow_div_k_minus);
#else
	return false;
#endif
      }
      else
	// diff_coeff_pow_div_k == 2.
	return false;
    }
    // 0 < k < 1.
    else {
      const Expr& pow_c_b_k = pwr(c_b, k);
      if (diff_coeff_pow_div_k == -1)
	bound = pow_n_k * pow_div_k / (pow_div_k - coeff)
	  - pow_div_k / (pow_div_k - coeff) * pwr(pow_div_k / coeff, 1 + mu_b)
	  * pow_frac_log
	  - pow_c_b_k * (log(Recurrence::n) / log(divisor) - 1);
      else if (diff_coeff_pow_div_k == 0)
	bound = (log(Recurrence::n) / log(divisor) - mu_b)
	  * (pow_n_k - pow_c_b_k) + pow_c_b_k;
      else if (diff_coeff_pow_div_k == 1)
	bound = pow_div_k / (coeff - pow_div_k)
	  * pwr(pow_div_k / coeff, 1 + mu_b)
	  * pow_frac_log
	  - pow_div_k / (coeff - pow_div_k) * pow_n_k
	  - pow_c_b_k * (log(Recurrence::n) / log(divisor) - 1);
      else
	// diff_coeff_pow_div_k == 2.
	return false;
    }
  }
  // Upper bound.
  else
    if (diff_coeff_pow_div_k == -1)
      bound = pow_n_k * pow_div_k / (pow_div_k - coeff);
    else if (diff_coeff_pow_div_k == 0)
      bound = pow_n_k * frac_log;
    else if (diff_coeff_pow_div_k == 1)
      bound = pow_div_k
	* (pwr(Recurrence::n, log(coeff) / log(divisor)) - pow_n_k)
	/ (coeff - pow_div_k);
    else
      // diff_coeff_pow_div_k == 2.
      return false;
  return true;
}

bool
compute_bounds_for_logarithm_function(Type_Bound type, const Number& coeff,
				      const Number& divisor, Expr& bound) {
  // Lower bound.
  if (type == LOWER_BOUND)
    // FIXME: for the moment use only the theorem in order to compute
    // lower bound in the case `g(n) = log(n)'.
#if 0
    if (divisor.is_positive_integer()) {
      const Expr& frac_log = log(Recurrence::n) / log(divisor);
      if (coeff < 1)
	bound = log(divisor) * exact_pwr(coeff - 1, -2)
	  * ((1 - coeff) * frac_log - 1
	     + pwr(Recurrence::n, log(coeff) / log(divisor)));
      else if (coeff == 1)
	bound = log(Recurrence::n) / (2 * coeff) * (frac_log - 1);
      else {
	assert(coeff > 1);
	bound = log(divisor) * exact_pwr(coeff - 1, -2)
	  * (pwr(Recurrence::n, log(coeff) / log(divisor))
	     - (coeff - 1) * frac_log - coeff);
      }
    }
    else
#endif
      return false;
  // Upper bound.
  else {
    if (coeff == 1)
      bound = Number(1, 2) * log(Recurrence::n)
	* (log(Recurrence::n) / log(divisor) + 1);
    else
      bound = log(divisor) / exact_pwr(coeff - 1, 2)
	* (coeff * pwr(Recurrence::n, log(coeff) / log(divisor))
	   - (coeff - 1) * log(Recurrence::n) / log(divisor) - coeff);
  }
  return true;
}

bool
compute_bounds_for_power_times_logarithm_function(Type_Bound type,
						  const Number& coeff,
						  const Number& divisor,
						  const Number& k,
						  Expr& bound) {
  const Expr& pow_div_k = pwr(divisor, k);
  int diff_coeff_pow_div_k = compare(coeff, pow_div_k);
  const Expr& pow_n_k = pwr(Recurrence::n, k);
  // Lower bound.
  if (type == LOWER_BOUND)
    // FIXME: for the moment use only the theorem in order to compute
    // lower bound in the case `g(n) = n^k log(n)'.
#if 0
    if (divisor.is_positive_integer()) {
      const Expr& frac_log = log(Recurrence::n) / log(divisor);
      if (diff_coeff_pow_div_k == -1)
	bound = log(divisor) * pwr(coeff - pow_div_k, -2)
	  * (((pow_div_k - coeff) * frac_log - pow_div_k) * pow_n_k
	     + pow_div_k * pwr(Recurrence::n, log(coeff) / log(divisor)));
      else if (diff_coeff_pow_div_k == 0)
	bound = pow_n_k * log(Recurrence::n) / (2 * coeff)
	  * (frac_log - 1);
      else if (diff_coeff_pow_div_k == 1)
	bound = pow_div_k * log(divisor) * pwr(coeff - pow_div_k, -2)
	  * (pwr(Recurrence::n, log(coeff) / log(divisor))
	     - ((coeff - pow_div_k) * frac_log + coeff) * pow_n_k);
      else
	// diff_coeff_pow_div_k == 2.
	return false;
    }
    else
#endif
      return false;
  // Upper bound.
  else {
    const Expr& pow_n_k_log = pow_n_k * log(Recurrence::n);
    if (diff_coeff_pow_div_k == 0)
      bound = Number(1, 2) * pow_n_k_log
	* (log(Recurrence::n) / log(divisor) + 1);
    else if (diff_coeff_pow_div_k == -1 || diff_coeff_pow_div_k == 1)
      bound = (pow_div_k * log(divisor) / pwr(coeff - pow_div_k, 2)) 
	* (coeff * pwr(Recurrence::n, log(coeff) / log(divisor))
	   - (coeff - pow_div_k) * pow_n_k_log / log(divisor) - coeff);
    else
      // diff_coeff_pow_div_k == 2.
      return false;
  }
  return true;
}

// Swap the type of the bound.
Type_Bound swap_type(Type_Bound type) {
  return type == LOWER_BOUND ? UPPER_BOUND : LOWER_BOUND;
}

/*!
  g(n) = a * n^k.
*/
bool
sharper_bounds_for_polynomial_function(Type_Bound type, const Expr& poly_coeff,
				       const Number& coeff,
				       const Number& divisor,
				       Expr& bound) {
  // Case `g(n) = a n^0', `a' number.
  if (poly_coeff.is_a_number()
      || ((poly_coeff.is_a_constant_function(Recurrence::n)
	   || poly_coeff.is_a_constant_power(Recurrence::n))
	  && !has_parameters(poly_coeff))) {
    int cmp = compare(poly_coeff, 0);
    Expr tmp_bound;
    if (cmp == 1) {
      // `a' positive number.
      if (!compute_bounds_for_power_of_n(type, coeff, divisor, 0, tmp_bound))
	return false;
    }
    else if (cmp == -1) {
      // `a' negative number -> swap lower with upper or upper with lower.
      if (!compute_bounds_for_power_of_n(swap_type(type), coeff, divisor, 0,
					 tmp_bound))
	return false;
    }
    else
      return false;
    bound += poly_coeff * tmp_bound;
    return true;
  }
  // Case `g(n) = n^1'.
  else if (poly_coeff == Recurrence::n) {
    Expr tmp_bound;
    if (!compute_bounds_for_power_of_n(type, coeff, divisor, 1, tmp_bound))
      return false;
    bound += tmp_bound;
    return true;
  }
  // Case `g(n) = -n^1' -> swap lower with upper or upper with lower.
  else if (poly_coeff == -Recurrence::n) {
    Expr tmp_bound;
    if (!compute_bounds_for_power_of_n(swap_type(type), coeff, divisor, 1,
				       tmp_bound))
      return false;
    bound += tmp_bound;
    return true;
  }
  // Case `g(n) = n^k' with `k > 1' integer.
  else if (poly_coeff.is_a_power()
	   && poly_coeff.arg(0) == Recurrence::n
	   && poly_coeff.arg(1).is_a_number()) {
    Number k = poly_coeff.arg(1).ex_to_number();
    if (k.is_positive()) {
      Expr tmp_bound;
      if (!compute_bounds_for_power_of_n(type, coeff, divisor, k, tmp_bound))
	return false;
      bound += tmp_bound;
      return true;
    }
  }
  // Case `g(n) = a n^k'.
  if (poly_coeff.is_a_mul()) {
    Number k = 0;
    Expr factor = 1;
    for (unsigned int i = poly_coeff.nops(); i-- > 0; ) {
      const Expr& poly_coeff_i = poly_coeff.op(i);
      if (poly_coeff_i == Recurrence::n)
	k = 1;
      else if (poly_coeff_i.is_a_power()
	       && poly_coeff_i.arg(0) == Recurrence::n
	       && poly_coeff_i.arg(1).is_a_number()) {
	assert(poly_coeff_i.arg(1).ex_to_number().is_positive_integer());
	k = poly_coeff_i.arg(1).ex_to_number();
      }
      else if (poly_coeff_i.is_a_number()
	       || ((poly_coeff_i.is_a_constant_function(Recurrence::n)
		    || poly_coeff_i.is_a_constant_power(Recurrence::n))
		   && !has_parameters(poly_coeff_i)))
	factor *= poly_coeff_i;
      else
	return false;
    }
    int cmp = compare(factor, 0);
    Expr tmp_bound;
    // `a' positive number.
    if (cmp == 1) {
      if (!compute_bounds_for_power_of_n(type, coeff, divisor, k, tmp_bound))
	return false;
    }
    // `a' negative number -> swap lower with upper or upper with lower.
    else if (cmp == -1) {
      if (!compute_bounds_for_power_of_n(swap_type(type), coeff, divisor, k,
					 tmp_bound))
	return false;
    }
    else
      return false;
    bound += factor * tmp_bound;
    return true;
  }
  return false;
}

/*!
  g(n) = a n^k or g(n) = a log(n) or g(n) = a * n^k log(n).
*/
bool
sharper_bounds_for_no_polynomial_function(Type_Bound type,
					  const Expr& no_poly_coeff,
					  const Number& coeff,
					  const Number& divisor,
					  Expr& bound) {
  Number k;
  // Case `g(n) = log(n)'.
  Expr tmp_bound;
  if (no_poly_coeff.is_the_log_function()
      && no_poly_coeff.arg(0) == Recurrence::n
      && compute_bounds_for_logarithm_function(type, coeff, divisor,
					       tmp_bound)) {
    bound += tmp_bound;
    return true;
  }
  // Case `g(n) = n^k' with `k' rational number.
  else if (no_poly_coeff.is_a_power()
	   && no_poly_coeff.arg(0) == Recurrence::n
	   && no_poly_coeff.arg(1).is_a_number(k)
	   && k.is_positive()) {
    Expr tmp_bound;
    if (!compute_bounds_for_power_of_n(type, coeff, divisor, k, tmp_bound))
      return false;
    bound += tmp_bound;
    return true;
  }
  else if (no_poly_coeff.is_a_mul()) {
    // `a n^k log(n)'
    Expr factor = 1;
    Number k = 0;
    Expr log_part = 0;
    for (unsigned int i = no_poly_coeff.nops(); i-- > 0; ) {
      const Expr& no_poly_coeff_i = no_poly_coeff.op(i);
      if (no_poly_coeff_i == log(Recurrence::n))
	log_part = no_poly_coeff_i;
      else if (no_poly_coeff_i == Recurrence::n)
	k = 1;
      else if (no_poly_coeff_i.is_a_power()
	       && no_poly_coeff_i.arg(0) == Recurrence::n
	       && no_poly_coeff_i.arg(1).is_a_number()) {
	assert(no_poly_coeff_i.arg(1).ex_to_number().is_positive());
	k = no_poly_coeff_i.arg(1).ex_to_number();
      }
      else if (no_poly_coeff_i.is_a_number()
	       || ((no_poly_coeff_i.is_a_constant_function(Recurrence::n)
		    || no_poly_coeff_i.is_a_constant_power(Recurrence::n))
		   && !has_parameters(no_poly_coeff_i)))
	factor *= no_poly_coeff_i;
      else
	return false;
    }
    int cmp = compare(factor, 0);
    // `a log(n)'.
    if (k == 0) {
      assert(!log_part.is_zero());
      // `a' positive number.
      Expr tmp_bound;
      if (cmp == 1) {
	if (!compute_bounds_for_logarithm_function(type, coeff, divisor,
						   tmp_bound))
	  return false;
      }
      // `a' negative number -> swap lower with upper or upper with lower.
      else if (cmp == -1) {
	if (!compute_bounds_for_logarithm_function(swap_type(type), coeff,
						   divisor, tmp_bound))
	  return false;
      }
      else
	return false;
      bound += factor * tmp_bound;
      return true;
    }
    else if (log_part == 0) {
      assert(k != 0);
      Expr tmp_bound;
      // `a' positive number.
      if (cmp == 1) {
	if (!compute_bounds_for_power_of_n(type, coeff, divisor, k,
					   tmp_bound))
	  return false;
      }
      // `a' negative number -> swap lower with upper or upper with lower.
      else if (cmp == -1) {
	if (!compute_bounds_for_power_of_n(swap_type(type), coeff, divisor, k,
					   tmp_bound))
	  return false;
      }
      else
	return false;
      bound += factor * tmp_bound;
      return true;
    }
    else {
      assert(k != 0 && !log_part.is_zero());
      Expr tmp_bound;
      // `a' positive number.
      if (cmp == 1) {
	if (!compute_bounds_for_power_times_logarithm_function(type, coeff,
							       divisor, k,
							       tmp_bound))
	  return false;
      }
      // `a' negative number -> swap lower with upper or upper with lower.
      else if (cmp == -1) {
	if (!compute_bounds_for_power_times_logarithm_function(swap_type(type),
							       coeff, divisor,
							       k, tmp_bound))
	  return false;
      }
      else
	return false;
      bound += factor * tmp_bound;
      return true;
    }
  }
  return false;
}

/*!  
  g(n) = a c^n, c > 1, a number.
*/
bool
sharper_bounds_for_exponential(Type_Bound type,
			       const Number& base, const Expr& poly_coeff,
			       const Number& coeff, const Number& divisor,
			       Expr& bound) {
  if (coeff >= 1) {
    // `a' number or constant function or constant power. 
    if (poly_coeff.is_a_number()
	|| ((poly_coeff.is_a_constant_function(Recurrence::n)
	     || poly_coeff.is_a_constant_power(Recurrence::n))
	    && !has_parameters(poly_coeff))) {
      int cmp = compare(poly_coeff, 0);
      Expr tmp_bound;
      if (cmp == 1) {
	// `a' positive number.
	if (!compute_bounds_for_exp_function(type, coeff, divisor, base,
					     tmp_bound))
	  return false;
      }
      else if (cmp == -1) {
	// `a' negative number -> swap lower with upper or upper with lower.
	if (!compute_bounds_for_exp_function(swap_type(type), coeff, divisor,
					     base, tmp_bound))
	  return false;
      }
      else
	return false;
      bound += poly_coeff * tmp_bound;
      return true;
    }
    // `a' is a product of numbers, constant functions and constant powers.
    if (poly_coeff.is_a_mul()) {
      Expr factor = 1;
      for (unsigned int i = poly_coeff.nops(); i-- > 0; ) {
	const Expr& poly_coeff_i = poly_coeff.op(i);
	if (poly_coeff_i.is_a_number()
	    || ((poly_coeff_i.is_a_constant_function(Recurrence::n)
		 || poly_coeff_i.is_a_constant_power(Recurrence::n))
		&& !has_parameters(poly_coeff_i)))
	  factor *= poly_coeff_i;
	else
	  return false;
      }
      int cmp = compare(factor, 0);
      // `a' positive number.
      Expr tmp_bound;
      if (cmp == 1) {
	if (!compute_bounds_for_exp_function(type, coeff, divisor, base,
					     tmp_bound))
	return false;
      }
      else if (cmp == -1) {
	// `a' negative number -> swap lower with upper or upper with lower.
	if (!compute_bounds_for_exp_function(swap_type(type), coeff, divisor,
					     base, tmp_bound))
	  return false;
      }
      else
	return false;
      bound += factor * tmp_bound;
      return true;
    }
  }
  return false;
}

/*!
  If \p type is <CODE>LOWER_BOUND</CODE> then this function computes
  the sum \f$ G(n) = sum_{k = 1}^n a^{n-k} p(b^k) \f$, otherwise, when
  \p type is <CODE>UPPER_BOUND</CODE> computes the sum
  \f$ H(n) = sum_{k = 2}^n a^{n-k} p(b^k-1) \f$,
  where \f$ a \f$ is stored in \p coefficient, \f$ b \f$ is stored
  in \p divisor and \f$ p(n) \f$ is stored in \p summand.
*/
bool
compute_sum(Type_Bound type, const Expr& summand,
	    const Number& coefficient, const Number& divisor,
	    Expr& sum) {
  std::vector<Expr> bases_of_exp;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  Expr tmp;
  if (type == LOWER_BOUND)
    tmp = summand.substitute(Recurrence::n,
			     pwr(divisor, Recurrence::n));
  else
    tmp = summand.substitute(Recurrence::n,
			     pwr(divisor, Recurrence::n) - 1);
  exp_poly_decomposition(tmp.expand(), Recurrence::n,
			 bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
  D_VAR(tmp);
  D_VEC(bases_of_exp, 0, bases_of_exp.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);
  // FIXME: temporary!
  if (vector_not_all_zero(exp_poly_coeff)
      && vector_not_all_zero(exp_no_poly_coeff))
    return false;
  else if (vector_not_all_zero(exp_poly_coeff))
    for (unsigned int i = bases_of_exp.size(); i-- > 0; ) {
      Symbol k;
      sum += sum_poly_times_exponentials(exp_poly_coeff[i]
					 .substitute(Recurrence::n, k), k,
					 Recurrence::n,
					 bases_of_exp[i] / coefficient);
      // `sum_poly_times_exponentials' computes the sum from 0, whereas
      // we want that the sum starts from `1' in the case of lower bound
      // (type = LOWER_BOUND) or from `2' in the case of upper bound
      // (type = UPPER_BOUND).
      sum -= exp_poly_coeff[i].substitute(Recurrence::n, 0);
      if (type == UPPER_BOUND)
	sum -= exp_poly_coeff[i].substitute(Recurrence::n, 1)
	  * bases_of_exp[i] / coefficient;
    }
  else {
    Number index = type == LOWER_BOUND ? 1 : 2;
    for (unsigned int i = bases_of_exp.size(); i-- > 0; ) {
      Expr gosper_sum;
      if (!gosper_algorithm(Recurrence::n,
			    pwr(coefficient, -Recurrence::n) * tmp,
			    index, Recurrence::n, gosper_sum))
	return false;
      else
	sum += gosper_sum;
    }
  }
  sum *= pwr(coefficient, Recurrence::n);
  sum = simplify_ex_for_output(sum, false);
  D_VAR(sum);
  return true;
}

/*!
  Computes two different sums:
  \f$ G(n) = sum_{k = 1}^n a^{n-k} p(b^k) \f$ and
  \f$ H(n) = sum_{k = 2}^n a^{n-k} p(b^k-1) \f$,
  where \f$ a \f$ is stored in \p coefficient, \f$ b \f$ is stored
  in \p divisor and \f$ p(n) \f$ is stored in \p summand.
*/
void
try_to_compute_sum(Type_Bound type, const Expr& summand,
		   const Number& coefficient, const Number& divisor,
		   Expr& sum) {
  D_VAR(summand);
  if (type == LOWER_BOUND) {
    // Compute `G(n)' used for the lower bound.
    if (!compute_sum(LOWER_BOUND, summand, coefficient, divisor, sum))
      // We try to compute `H(n)' because `H(q) <= G(q)'.
      if (!compute_sum(UPPER_BOUND, summand, coefficient, divisor, sum)) {
	Symbol h;
	sum = PURRS::sum(h, 1, Recurrence::n, pwr(coefficient, Recurrence::n-h)
			 * summand.substitute(Recurrence::n, pwr(divisor, h)));
      }
  }
  else
    // Compute `H(n)' used for the upper bound.
    if (!compute_sum(UPPER_BOUND, summand, coefficient, divisor, sum))
      // We try to compute `G(n)' because `H(q+1) <= G(q+1)'.
      if (!compute_sum(LOWER_BOUND, summand, coefficient, divisor, sum)) {
	Symbol h;
	sum = PURRS::sum(h, 2, Recurrence::n, pwr(coefficient, Recurrence::n-h)
			 * summand.substitute(Recurrence::n,
					      pwr(divisor, h) - 1));
      }
}

//! \brief
//! Returns <CODE>true</CODE> if the system is able to approximate the
//! functional equation \f$ x(n) = a x(n/b) + g(n) \f$ with
//! \ p coefficient as \f$ a \f$, \p divisor_arg as \f$ b \f$ and
//! \p inhomogeneous as \f$ g(n) \f$; returns <CODE>false</CODE> otherwise.
bool
known_class_of_functional_eq_rank_1(const Expr& coefficient,
				    const Number& divisor_arg,
				    const Expr& inhomogeneous) {
  assert(divisor_arg.is_rational() && divisor_arg.is_positive());
  // We want that `coefficient' is a positive number and
  // `divisor_arg' a rational number bigger than `1'.
  Number coeff;
  if (!coefficient.is_a_number(coeff) || !coeff.is_positive()
      || divisor_arg < 1)
    return false;

  // FIXME: temporary!
  if (!divisor_arg.is_integer())
    return false;

  if (has_parameters(inhomogeneous))
    return false;

  return true;
}

/*!
  Computes the bound relative to the \f$ i \f$-th addend of the
  form \f$ a b^n \f$ of the non homogeneous term of the functional equation,
  where \f$ a \f$ is contained in \p e, while \f$ b \f$ is contained
  in \p base.

  \param poly         <CODE>true</CODE> if \p e is a (syntactically) a
                      polynomial in the variable \p Recurrence::n;
		      <CODE>false</CODE> otherwise.
  \param e            The coefficient of the exponential in the variable
                      \p Recurrence::n the whose base is \p base.
  \param base         The base of the exponential in the variable
                      \p Recurrence::n.
  \param type         hold <CODE>LOWER_BOUND</CODE> if the system is
                      computing the lower bound, <CODE>UPPER_BOUND</CODE>
		      if it is computing the upper bound.
  \param coeff        The coefficient \f$ \alpha \f$ of the functional
                      equation \f$ x(n) = \alpha x(n / \beta) + g(n) \f$.
  \param divisor_arg  The positive number \f$ \beta \f$ of the functional
                      equation \f$ x(n) = \alpha x(n / \beta) + g(n) \f$.
  \param bound        The part of the bound computed by this function
                      that is that relative to the non homogeneous term
		      \f$ g(n) \f$ of the functional equation
		      \f$ x(n) = \alpha x(n / \beta) + g(n) \f$.
  \param condition    The positive integer starting from which the
                      inhomogeneous term \f$ g(n) \f$ of the functional
		      equation \f$ x(n) = \alpha x(n / \beta) + g(n) \f$
		      is a non negative, non decreasing function.

  \return             <CODE>true</CODE> if the function is able to
                      compute the bound; <CODE>false</CODE> otherwise.
*/
bool
compute_poly_or_no_poly(bool poly, const Expr& e, const Number& base,
			Type_Bound type,
			const Number& coeff, const Number& divisor_arg,
			Expr& bound, Number& condition) {
  if (e.is_a_add())
    for (unsigned int i = e.nops(); i-- > 0; ) {
      const Expr& addend = e.op(i);
      if ((poly && !sharper_bounds_for_polynomial_function(type, addend,
							   coeff, divisor_arg,
							   bound))
	  || (!poly
	      && !sharper_bounds_for_no_polynomial_function(type, addend,
							    coeff, divisor_arg,
							    bound))
	  || base != 1) {
	// Check if the i-th addend of `e' is a non-negative,
	// non-decreasing function.
	if (!is_non_negative_non_decreasing(base, addend, poly,
					    Recurrence::n, true,
					    condition))
	  return false ;
	Expr sum;
	try_to_compute_sum(type, addend * pwr(base, Recurrence::n),
			   coeff, divisor_arg, sum);
	if (type == LOWER_BOUND)
	  bound += sum.substitute(Recurrence::n,
				  log(Recurrence::n) / log(divisor_arg) - 1);
	else
	  bound += sum.substitute(Recurrence::n,
				  log(Recurrence::n) / log(divisor_arg)  + 1);
      }
    }
  else
    if (!e.is_zero()
	&& ((poly && !sharper_bounds_for_polynomial_function(type, e, coeff,
							     divisor_arg,
							     bound))
	    || (!poly
		&& !sharper_bounds_for_no_polynomial_function(type, e, coeff,
							      divisor_arg,
							      bound)))
	|| base != 1) {
      // Check if `e' is a non-negative, non-decreasing function.
      if (!is_non_negative_non_decreasing(base, e, poly, Recurrence::n, true,
					  condition))
	return false;
      Expr sum;
      try_to_compute_sum(type, e, coeff, divisor_arg, sum);
      if (type == LOWER_BOUND)
	bound += sum.substitute(Recurrence::n,
				log(Recurrence::n) / log(divisor_arg) - 1);
      else
	bound += sum.substitute(Recurrence::n,
				log(Recurrence::n) / log(divisor_arg) + 1);
    }
  return true;
}

/*!
  Let \f$ e(n) \f$ be the inhomogeneous term of the functional equation
  contained in \p inhomogeneous.
  This function decomposes the inhomogeneous in the following way:
  \f$ e(n) = \sum_{i=0}^k \alpha_i^n \bigl(p_i(n) + q_i(n)\bigr) \f$, where
  - \f$ \alpha_i \f$ is a expression valid for to be an exponential's base.
    (syntactically different from \p 0);
  - \f$ \alpha_i \neq \alpha_j \f$ if \f$ i \neq j \f$;
  - \f$ p_i(n) \f$ is (syntactically) a polynomial in \f$ n \f$.
  Each expression \f$ \alpha_i^n p_i(n) \f$ and \f$ \alpha_i^n q_i(n) \f$
  is ulteriorly divided in many expression how many are the addends
  of \f$ p_i(n) \f$ and \f$ q_i(n) \f$.
  For each expression computes the lower and the upper bound:
  - try to apply the sharper bounds;
  - try to apply the theorem 2.1 of ...
  - if the previous tests are failed returns the symbolic sum.

  Returns <CODE>true</CODE> if the function is able to compute the bound
  relative to the non homogeneous term of the functional equation;
  <CODE>false</CODE> otherwise.
*/
bool
compute_non_homogeneous_part(Type_Bound type,
			     const Number& coeff, const Number& divisor_arg,
			     const Expr& inhomogeneous,
			     Expr& bound, Number& condition) {
  std::vector<Expr> bases_of_exp;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous.expand(), Recurrence::n,
			 bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
  for (unsigned int i = bases_of_exp.size(); i-- > 0; ) {
    const Expr& base = bases_of_exp[i];
    const Expr& poly_coeff = exp_poly_coeff[i];
    const Expr& no_poly_coeff = exp_no_poly_coeff[i];
    Number num_base;
    if (base.is_a_number(num_base) && num_base.is_positive()) {
      // Base of exponential is equal to `1'.
      if (num_base == 1) {
	// Consider the eventual polynomial part.
	if (!compute_poly_or_no_poly(true, poly_coeff, 1, type,
				     coeff, divisor_arg, bound, condition))
	  return false;
	// Consider the eventual non-polynomial part.
	if (!compute_poly_or_no_poly(false, no_poly_coeff, 1, type,
				     coeff, divisor_arg, bound, condition))
	  return false;
      }
      // Base of exponential is a number different from `1'.
      else {
	if (!poly_coeff.is_zero())
	  // Apply, if it possible, the sharper bounds for the exponential.
	  if (!sharper_bounds_for_exponential(type, num_base, poly_coeff,
					      coeff, divisor_arg, bound))
	    // Consider the eventual polynomial part.
	    if (!compute_poly_or_no_poly(true, poly_coeff, num_base, type,
					 coeff, divisor_arg, bound, condition))
	      return false;
	if (!no_poly_coeff.is_zero())
	  // In this case is not possible to apply the sharper
	  // bounds for the exponential.
	  // Consider the eventual non-polynomial part.
	  if (!compute_poly_or_no_poly(false, no_poly_coeff, num_base,
				       type, coeff, divisor_arg,
				       bound, condition))
	    return false;
      }
    }
    // The base of the exponential is not a number.
    else
      return false;
  }
  return true;
}

/*!
  Computes
  \f[
    a^q x( \lfloor \frac {n}{b^q} \rfloor )
  \f]
  where \f$ a \f$ is stored in \p coefficient, \f$ b \f$
  is stored \p divisor.
*/
void
compute_term_about_initial_condition(Type_Bound type, unsigned int condition,
				     const Number& divisor,
				     const Expr& coefficient,
				     const Expr& q, Expr& bound) {
  // `index_initial_condition' is the argument of the `x' function:
  // `floor(n/b^(log(n)/log(b)))'.
  Expr index_initial_condition;
  if (divisor == 2) {
    index_initial_condition = 1;
    if (index_initial_condition.ex_to_number() < condition)
      index_initial_condition = condition;
  }
  else
    index_initial_condition = Sc(Recurrence::n, divisor);
  
  // Compute `a^q' removing the function `floor' from `q':
  // in the case of "upper bound" consider `a^(log(n)/log(b))';
  // in the case of "lower bound" consider `a^(log(n)/log(b)-1)';
  assert(q.is_the_floor_function());
  Expr tmp = q.arg(0);
  Number num;
  if (coefficient.is_a_number(num)) {
    assert(num > 0);
    if (num < 1 && type == UPPER_BOUND)
      tmp -= 1;
    if (num > 1 && type == LOWER_BOUND)
      tmp -= 1;
  }
  else
    if (type == LOWER_BOUND)
      tmp -= 1;
  bound += pwr(coefficient, tmp) * x(index_initial_condition);
}
					     
} // anonymous namespace

/*!
  If \p kind_of_bound is <CODE>LOWER</CODE> computes the lower bound
  for the functional equation \p *this; otherwise, if \p kind_of_bound
  is <CODE>UPPER</CODE>, computes the upper bound.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::approximate_functional_equation(Bound kind_of_bound) const {
  assert(kind_of_bound == LOWER || kind_of_bound == UPPER);
  if (functional_eq_p->rank() == 1) {
    Number divisor_arg = functional_eq_p->ht_begin()->first;
    Expr coefficient = functional_eq_p->ht_begin()->second;    
    
    if (!known_class_of_functional_eq_rank_1(coefficient, divisor_arg,
					     inhomogeneous_term))
      return TOO_COMPLEX;

    Number coeff = coefficient.ex_to_number();
    const Expr& q = floor(log(n) / log(divisor_arg));

    Expr bound;
    Number condition = -1;
    Type_Bound type = kind_of_bound == LOWER ? LOWER_BOUND : UPPER_BOUND;
    if (!inhomogeneous_term.is_zero()
	&& !compute_non_homogeneous_part(type, coeff, divisor_arg,
					 inhomogeneous_term,
					 bound, condition))
      return TOO_COMPLEX;

    if (condition > 1)
      set_applicability_condition(condition.to_unsigned_int());
    compute_term_about_initial_condition(type, applicability_condition(),
					 divisor_arg, coefficient, q, bound);

    if (kind_of_bound == LOWER) {
      lower_bound_.set_expression(simplify_logarithm(bound));
      if (upper_bound_.has_expression()
	  && upper_bound_.expression() == lower_bound_.expression())
	exact_solution_.set_expression(upper_bound_.has_expression());
    }
    else {
      upper_bound_.set_expression(simplify_logarithm(bound));
      if (lower_bound_.has_expression()
	  && upper_bound_.expression() == lower_bound_.expression())
	exact_solution_.set_expression(lower_bound_.has_expression());
    }
    return SUCCESS;
  }
  else
    return TOO_COMPLEX;
}
