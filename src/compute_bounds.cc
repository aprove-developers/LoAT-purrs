/* Recurrence class implementation: methods (and associated auxiliary
   functions) associated to the computation of lower and upper bounds.
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

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;


void
compute_bounds_for_exp_function(bool lower, const Number& coeff,
				const Number& divisor, const Number& base,
				Expr& bound) {
  // Lower bound.
  if (lower)
    bound += pwr(base, Recurrence::n);
  // Upper bound.
  else {
    Expr constant;
    if (pwr(base, pwr(divisor, 2)) > coeff * pwr(base, divisor)) {
      const Expr& tmp = pwr(base, pwr(divisor, 2));
      constant = coeff * tmp * pwr(tmp - coeff * pwr(base, divisor), -1);
    }
    else
      constant = coeff * pwr(coeff - 1, -1) * pwr(base, divisor)
	* pwr(log(base), -1 - log(coeff) * pwr(log(divisor), -1))
	* gamma(1 + log(coeff) * pwr(log(divisor), -1));
    bound += pwr(base, Recurrence::n)
      + constant * pwr(base, Recurrence::n / divisor);
  }
}

bool
compute_bounds_for_power_of_n(bool lower,
			      const Number& coeff, const Number& divisor,
			      const Number& num, const Number& k,
			      Expr& bound) {
  assert(!k.is_negative());
  // Lower bound.
  // FIXME: check the lower!!
  if (lower) {
    if ((k == 0)
	|| ((coeff == pwr(divisor, k-1) || coeff == pwr(divisor, k))
	    &&  k >= 1)) {
      Expr tmp_lb;
      assert(k == 0 || k >= 1);
      const Expr& log_log = log(Recurrence::n) / log(divisor);
      if (k == 0) {
	const Expr& n_log_log = pwr(Recurrence::n, log(coeff) / log(divisor));
	if (coeff < 1)
	  tmp_lb = (coeff - n_log_log) / (coeff * (1 - coeff));
	else if (coeff == 1)
	  if (divisor.is_integer())
	    tmp_lb = log_log - 1;
	  else // `divisor' is a rational (non-integer) number.  
	    tmp_lb = log_log
	      + log((divisor - 1) / (2 * divisor - 1)) / log(divisor); 
	else
	  tmp_lb = (n_log_log - coeff) / (coeff * (coeff - 1));
      }
      else {
	assert(k >= 1);
	const Expr& inv_div_minus_one = pwr(divisor - 1, -1);
	if (coeff == pwr(divisor, k-1))
	  tmp_lb = divisor * inv_div_minus_one * pwr(Recurrence::n, k)
	    - (pwr(divisor, 2) * inv_div_minus_one + k * log_log - k)
	    * pwr(Recurrence::n, k-1);
	else {
	  assert(coeff == pwr(divisor, k));
	  tmp_lb = (log_log - 1 - k * inv_div_minus_one)
	    * pwr(Recurrence::n, k)
	    + k * pwr(Recurrence::n, k-1) * divisor * inv_div_minus_one;
	}
      }
      bound += num * tmp_lb;
      return true;
    }
    else
      return false;
  }
  // Upper bound.
  else {
    // FIXME: come si puo' fare per evitare che saltino fuori i numeri in
    // virgola mobile? Usare tutte Expr non mi sembra una buona idea.
    const Expr& divisor_ex = divisor;
    const Expr& k_ex = k;
    const Expr& div_k = pwr(divisor_ex, k_ex);
    Number div_k_numeric = pwr(divisor, k);
    const Expr& tmp = pwr(Recurrence::n, k);
    Expr tmp_ub;
    if (coeff < div_k_numeric)
      tmp_ub = tmp * div_k / (div_k - coeff);
    else if (coeff == div_k_numeric)
      tmp_ub = tmp * log(Recurrence::n) / log(divisor);
    else
      tmp_ub = div_k
	* (pwr(Recurrence::n, log(coeff) / log(divisor)) - tmp)
	/ (coeff - div_k);
    bound += num * tmp_ub;
    return true;
  }
}

bool
compute_bounds_for_logarithm_function(bool lower, const Number& coeff,
				      const Number& divisor, const Number& num,
				      Expr& bound) {
  // Lower bound.
  if (lower)
    return false;
  // Upper bound.
  else {
    Expr tmp_bound;
    if (coeff < 1)
      return false;
    if (coeff == 1)
      tmp_bound = Number(1, 2) * log(Recurrence::n)
	* (log(Recurrence::n) / log(divisor) + 1);
    else {
      assert(coeff > 1);
      const Expr& ratio_log = log(coeff) / log(divisor);
      tmp_bound = (1 / (coeff - 1)) * log(Recurrence::n)
	* (pwr(Recurrence::n, ratio_log) - 1)
	+ log(divisor) * pwr(coeff - 1, -2)
	* ((coeff + 1) * pwr(Recurrence::n, ratio_log) - coeff);
    }
    bound =+ num * tmp_bound;
    return true;
  }
}

bool
compute_bounds_for_power_times_logarithm_function(bool lower,
						  const Number& coeff,
						  const Number& divisor,
						  const Number& num,
						  const Number& k,
						  Expr& bound) {
  // Lower bound.
  if (lower)
    return false;
  // Upper bound.
  else {
    Expr tmp_bound;
    Number div_k = pwr(divisor, k);
    const Expr& tmp = pwr(Recurrence::n, k) * log(Recurrence::n);
    if (coeff < div_k)
      tmp_bound = (div_k * tmp) / (div_k - coeff);
    else if (coeff == div_k)
      tmp_bound = tmp * log(Recurrence::n) / log(divisor);
    else
      tmp_bound = (div_k / (coeff - div_k)) 
	* (pwr(Recurrence::n, log(coeff) / log(divisor)) * log(Recurrence::n)
	   - tmp);
    bound += num * tmp_bound;
    return true;
  }
}

/*!
  g(n) = a * n^k.
*/
bool
sharper_bounds_for_polynomial_function(bool lower, const Expr& poly_coeff,
				       const Number& coeff,
				       const Number& divisor,
				       Expr& bound) {
  // Case `g(n) = a n^0' with `a' positive number.
  Number num;
  if (poly_coeff.is_a_number(num) && num.is_positive()
      && compute_bounds_for_power_of_n(lower, coeff, divisor, num, 0, bound))
    return true;
  // Case `g(n) = n^1'.
  else if (poly_coeff == Recurrence::n
	   && compute_bounds_for_power_of_n(lower, coeff, divisor, 1, 1,
					    bound))
    return true;
  // Case `g(n) = n^k' with `k > 1' integer.
  else if (poly_coeff.is_a_power()
	   && poly_coeff.arg(0) == Recurrence::n
	   && poly_coeff.arg(1).is_a_number(num)
	   && num.is_positive()
	   && compute_bounds_for_power_of_n(lower, coeff, divisor, 1, num,
					    bound))
    return true;
  else if (poly_coeff.is_a_mul() && poly_coeff.nops() == 2) {
    const Expr& first = poly_coeff.op(0);
    const Expr& second = poly_coeff.op(1);
    // Case `g(n) = a * n^k', with `a' and `k' positive integer numbers.
    if (first.is_a_number(num) && num.is_positive()) {
      if (second == Recurrence::n
	  && compute_bounds_for_power_of_n(lower, coeff, divisor, num, 1,
					   bound))
	return true;
      Number k;
      if (second.is_a_power()
	  && second.arg(0) == Recurrence::n
	  && second.arg(1).is_a_number(k)
	  && k.is_positive()
	  && compute_bounds_for_power_of_n(lower, coeff, divisor, num, k,
					   bound))
	return true;
    }
    else if (second.is_a_number(num) && num.is_positive()) {
      if (first == Recurrence::n
	  && compute_bounds_for_power_of_n(lower, coeff, divisor, num, 1,
					   bound))
	return true;
      Number k;
      if (first.is_a_power()
	  && first.arg(0) == Recurrence::n
	  && first.arg(1).is_a_number(k)
	  && k.is_positive()
	  && compute_bounds_for_power_of_n(lower, coeff, divisor, num, k,
					   bound))
	return true;
    }
  }
  return false;
}

/*!
  g(n) = a n^k or g(n) = a log(n) or g(n) = a * n^k log(n).
*/
bool
sharper_bounds_for_no_polynomial_function(bool lower,
					  const Expr& no_poly_coeff,
					  const Number& coeff,
					  const Number& divisor,
					  Expr& bound) {
  Number k;
  // Case `g(n) = log(n)'.
  if (no_poly_coeff.is_the_log_function()
      && no_poly_coeff.arg(0) == Recurrence::n
      && compute_bounds_for_logarithm_function(lower, coeff, divisor, 1,
					       bound))
    return true;
  // Case `g(n) = n^k' with `k' rational number.
  else if (no_poly_coeff.is_a_power()
	   && no_poly_coeff.arg(0) == Recurrence::n
	   && no_poly_coeff.arg(1).is_a_number(k)
	   && k.is_positive()
	   && compute_bounds_for_power_of_n(lower, coeff, divisor, 1, k,
					    bound))
    return true;
  else if (no_poly_coeff.is_a_mul()) {
    // 3 possibilities: `a log(n)' or `n^k log(n)' or `a n^k'.
    if (no_poly_coeff.nops() == 2) {
      Number num;
      const Expr& first = no_poly_coeff.op(0);
      const Expr& second = no_poly_coeff.op(1);
      if (first == log(Recurrence::n)) {
	// Case `g(n) = a log(n)' with `a' positive number.
	if (second.is_a_number(num) && num.is_positive()
	    && compute_bounds_for_logarithm_function(lower, coeff, divisor,
						     num, bound))
	  return true;
	// Case `g(n) = n^k log(n)' with `k > 0'.
	// `k = 1'.
	if (second == Recurrence::n
	    && compute_bounds_for_power_times_logarithm_function(lower, coeff,
								 divisor,
								 1, 1, bound))
	  return true;
	// `k > 0 && k != 1'.
	else if (second.is_a_power() && second.arg(0) == Recurrence::n
		 && second.arg(1).is_a_number(k) && k.is_positive()
		 && compute_bounds_for_power_times_logarithm_function(lower,
								      coeff,
								      divisor,
								      1, k,
								      bound))
	  return true;
      }
      else if (second == log(Recurrence::n)) {
	// Case `g(n) = a log(n)' with `a' positive number.
	if (first.is_a_number(num) && num.is_positive()
	    && compute_bounds_for_logarithm_function(lower, coeff, divisor,
						     num, bound))
	  return true;
	// Case `g(n) = n^k log(n)' with `k > 0'.
	// `k = 1'.
	if (first == Recurrence::n
	    && compute_bounds_for_power_times_logarithm_function(lower, coeff,
								 divisor,
								 1, 1, bound))
	    return true;
	// `k > 0 && k != 1'.
	else if (first.is_a_power() && first.arg(0) == Recurrence::n
		 && first.arg(1).is_a_number(k) && k.is_positive()
		 && compute_bounds_for_power_times_logarithm_function(lower,
								      coeff,
								      divisor,
								      1, k,
								      bound))
	  return true;
      }
      // Case `g(n) = a n^k' with `a' positive number and `k' positive
      // rational non-integer number.
      else if (((first.is_a_number(num) && num.is_positive()
	   && second.is_a_power()&& second.arg(0) == Recurrence::n
	   && second.arg(1).is_a_number(k) && k.is_positive())
	  || (second.is_a_number(num) && num.is_positive()
	      && first.is_a_power()&& first.arg(0) == Recurrence::n
	      && first.arg(1).is_a_number(k) && k.is_positive()))
	  && compute_bounds_for_power_of_n(lower, coeff, divisor, num, k,
					   bound))
	return true;
    }
    // Case `g(n) = a n^k log(n)' with `a' positive number and `k >= 1'.
    if (no_poly_coeff.nops() == 3) {
      Expr log_part;
      unsigned j;
      for (unsigned i = 0; i <= 2; ++i)
	if (no_poly_coeff.op(i) == log(Recurrence::n)) {
	  j = i;
	  log_part = no_poly_coeff.op(i);
	}
      Number num_part;
      unsigned h;
      for (unsigned i = 0; i <= 2; ++i)
	if (i != j && no_poly_coeff.op(i).is_a_number(num_part)
	    && num_part.is_positive())
	  h = i;
      Expr power_part;
      for (unsigned i = 0; i <= 2; ++i)
	if (i != j && i != h)
	  power_part = no_poly_coeff.op(i);
      // `k == 1'.
      if (power_part == Recurrence::n
	  && compute_bounds_for_power_times_logarithm_function(lower, coeff,
							       divisor,
							       num_part, 1,
							       bound))
	return true;
      // `k != 1'.
      Number k;
      if (power_part.is_a_power() && power_part.arg(0) == Recurrence::n
	  && power_part.op(1).is_a_number(k) && k.is_positive()
	  && compute_bounds_for_power_times_logarithm_function(lower, coeff,
							       divisor,
							       num_part, k,
							       bound))
	return true;
    }
  }
  return false;
}

/*!  
  g(n) = c^n, c > 1
*/
bool
sharper_bounds_for_exponential(bool lower,
			       const Number& base, const Expr& poly_coeff,
			       const Number& coeff, const Number& divisor,
			       Expr& bound) {
  Number num;
  if (poly_coeff.is_a_number(num) && num.is_positive())
    if (coeff >= 1) {
      Expr tmp_bound;
      compute_bounds_for_exp_function(lower, coeff, divisor, base, tmp_bound);
      bound += num * tmp_bound;
      return true;
    }
  return false;
}

/*!
  If \p for_lower is true then this function compute the sum
  \f$ G(n) = sum_{k = 1}^n a^{n-k} p(b^k) \f$, otherwise
  compute the sum \f$ H(n) = sum_{k = 2}^n a^{n-k} p(b^k-1) \f$,
  where \f$ a \f$ is stored in \p coefficient, \f$ b \f$ is stored
  in \p divisor and \f$ p(n) \f$ is stored in \p summand.
*/
bool
compute_sum(bool for_lower, const Expr& summand,
	    const Number& coefficient, const Number& divisor,
	    Expr& sum) {
  std::vector<Expr> bases_of_exp;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  Expr tmp;
  if (for_lower)
    tmp = summand.substitute(Recurrence::n,
			     pwr(divisor, Recurrence::n));
  else
    tmp = summand.substitute(Recurrence::n,
			     pwr(divisor, Recurrence::n) - 1);
  exp_poly_decomposition(tmp, Recurrence::n,
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
    for (unsigned i = bases_of_exp.size(); i-- > 0; ) {
      Symbol k("k");
      sum += sum_poly_times_exponentials(exp_poly_coeff[i]
					 .substitute(Recurrence::n, k), k,
					 Recurrence::n,
					 bases_of_exp[i] / coefficient);
      // `sum_poly_times_exponentials' computes the sum from 0, whereas
      // we want that the sum starts from `1' if `for_lower' is true
      // or from `2' if `for_lower' is false.
      sum -= exp_poly_coeff[i].substitute(Recurrence::n, 0);
      if (!for_lower)
	sum -= exp_poly_coeff[i].substitute(Recurrence::n, 1)
	  * bases_of_exp[i] / coefficient;
    }
  else {
    Number index = for_lower ? 1 : 2;
    for (unsigned i = bases_of_exp.size(); i-- > 0; ) {
      Expr gosper_sum;
      if (!full_gosper(Recurrence::n, pwr(coefficient, -Recurrence::n) * tmp,
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
try_to_compute_sum(bool lower, const Expr& summand,
		   const Number& coefficient, const Number& divisor,
		   Expr& sum) {
  D_VAR(summand);
  if (lower) {
    // Compute `G(n)' used for the lower bound.
    if (!compute_sum(true, summand, coefficient, divisor, sum))
      // We try to compute `H(n)' because `H(q) <= G(q)'.
      if (!compute_sum(false, summand, coefficient, divisor, sum)) {
	Symbol h;
	sum = PURRS::sum(h, 1, Recurrence::n, pwr(coefficient, Recurrence::n-h)
			 * summand.substitute(Recurrence::n, pwr(divisor, h)));
      }
  }
  else
    // Compute `H(n)' used for the upper bound.
    if (!compute_sum(false, summand, coefficient, divisor, sum))
      // We try to compute `G(n)' because `H(q+1) <= G(q+1)'.
      if (!compute_sum(true, summand, coefficient, divisor, sum)) {
	Symbol h;
	sum = PURRS::sum(h, 2, Recurrence::n, pwr(coefficient, Recurrence::n-h)
			 * summand.substitute(Recurrence::n,
					      pwr(divisor, h) - 1));
      }
}

//! \brief
//! Returns <CODE>true<CODE> if the system is able to approximate the
//! functional equation \f$ x(n) = a x(n/b) + g(n) \f$ with
//! \ p coefficient as \f$ a \f$, \p divisor_arg as \f$ b \f$ and
//! \p inhomogeneous as \f$ g(n) \f$; returns <CODE>false<CODE> otherwise.
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

  if (has_parameters(inhomogeneous))
    return false;

  return true;
}

bool
compute_non_homogeneous_part(bool lower,
			     const Number& coeff, const Number& divisor_arg,
			     const Expr& inhomogeneous, const Expr& q,
			     Expr& bound, Number& condition) {
  std::vector<Expr> bases_of_exp;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous, Recurrence::n,
			 bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
  for (unsigned i = bases_of_exp.size(); i-- > 0; ) {
    const Expr& base = bases_of_exp[i];
    const Expr& poly_coeff = exp_poly_coeff[i];
    const Expr& no_poly_coeff = exp_no_poly_coeff[i];
    Number num_base;
    if (base.is_a_number(num_base) && num_base.is_positive_integer()) {
      // Base of exponential is equal to `1'.
      if (num_base == 1) {
	// Consider the eventual polynomial part.
	if (!poly_coeff.is_zero()
	    && !sharper_bounds_for_polynomial_function(lower, poly_coeff, 
						       coeff, divisor_arg,
						       bound)) {
	  // Check if the polynomial part is a non-negative,
	  // non-decreasing function.
	  if (!is_non_negative_non_decreasing(1, poly_coeff, true,
					      Recurrence::n, true, condition))
	    return false ;
	  Expr sum;
	  try_to_compute_sum(lower, poly_coeff, coeff, divisor_arg, sum);
	  if (lower)
	    bound += sum.substitute(Recurrence::n, q);
	  else
	    bound += sum.substitute(Recurrence::n, q + 1);
	}
	// Consider the eventual non-polynomial part.
	if (!no_poly_coeff.is_zero()
	    && !sharper_bounds_for_no_polynomial_function(lower, no_poly_coeff,
							  coeff, divisor_arg,
							  bound)) {
	  // Check if the non-polynomial part is a non-negative,
	  // non-decreasing function.
	  if (!is_non_negative_non_decreasing(1, no_poly_coeff, false,
					      Recurrence::n, true, condition))
	    return false;
	  Expr sum;
	  try_to_compute_sum(lower, no_poly_coeff, coeff, divisor_arg, sum);
	  if (lower)
	    bound += sum.substitute(Recurrence::n, q);
	  else
	    bound += sum.substitute(Recurrence::n, q + 1);
	}
      }
      // Base of exponential is a number different from `1'.
      else {
	if (!poly_coeff.is_zero())
	  // Apply, if it possible, the sharper bounds for the exponential.
	  if (!sharper_bounds_for_exponential(lower, num_base, poly_coeff,
					      coeff, divisor_arg, bound)) {
	    // Check if the polynomial part times esponential is a
	    // non-negative, non-decreasing function.
	    if (!is_non_negative_non_decreasing(num_base, poly_coeff, true,
						Recurrence::n, true,
						condition))
	      return false ;
	    Expr sum;
	    try_to_compute_sum(lower, pwr(base, Recurrence::n) * poly_coeff,
			       coeff, divisor_arg, sum);
	    if (lower)
	      bound += sum.substitute(Recurrence::n, q);
	    else
	      bound += sum.substitute(Recurrence::n, q + 1);
	  }
	if (!no_poly_coeff.is_zero()) {
	  // In this case is not possible to apply the sharper
	  // bounds for the exponential.
	  // Check if the polynomial part times exponential is a
	  // non-negative, non-decreasing function.
	  if (!is_non_negative_non_decreasing(num_base, no_poly_coeff, false,
					      Recurrence::n, true, condition))
	    return false ;
	  Expr sum;
	  try_to_compute_sum(lower, pwr(base, Recurrence::n) * no_poly_coeff,
			     coeff, divisor_arg, sum);
	  if (lower)
	    bound += sum.substitute(Recurrence::n, q);
	  else
	    bound += sum.substitute(Recurrence::n, q + 1);
	}
      }
    }
    // The base of the exponential is not a number.
    else
      return false;
  }
  return true;
}
					     
} // anonymous namespace

/*!
  Computes
  \f[
    a^q x( \frac {n}{b^q} )
  \f]
  where \f$ a \f$ is stored in
  <CODE>functional_eq_p->ht_begin()->second</CODE>, \f$ b \f$
  is stored <CODE>functional_eq_p->ht_begin()->first</CODE>.
*/
void
PURRS::Recurrence::add_term_with_initial_condition(bool lower, const Expr& q,
						   Expr& bound) const {
  Expr index_initial_condition;
  if (lower)
    index_initial_condition 
      = simplify_logarithm(n / pwr(functional_eq_p->ht_begin()->first, q + 1));
  else
    index_initial_condition 
      = simplify_logarithm(n / pwr(functional_eq_p->ht_begin()->first, q));
  
  if (index_initial_condition.is_a_number()
      && index_initial_condition.ex_to_number() < applicability_condition())
    index_initial_condition = applicability_condition();
  
  bound += pwr(functional_eq_p->ht_begin()->second, q)
    * x(index_initial_condition);
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::approximate_functional_equation_lower() const {
  if (functional_eq_p->rank() == 1) {
    Number divisor_arg = functional_eq_p->ht_begin()->first;
    Expr coefficient = functional_eq_p->ht_begin()->second;    
    
    if (!known_class_of_functional_eq_rank_1(coefficient, divisor_arg,
					     inhomogeneous_term))
      return TOO_COMPLEX;

    Number coeff = coefficient.ex_to_number();
    Expr q_lower = log(n) / log(divisor_arg) - 1;
    Expr lower_bound;

    Number condition = -1;
    if (!inhomogeneous_term.is_zero()
	&& !compute_non_homogeneous_part(true, coeff, divisor_arg,
					 inhomogeneous_term, q_lower,
					 lower_bound, condition))
      return TOO_COMPLEX;

    if (condition > 1)
      set_applicability_condition(condition.to_unsigned());
    add_term_with_initial_condition(true, q_lower, lower_bound);
    
    lower_bound_.set_expression(simplify_logarithm(lower_bound));
    if (upper_bound_.has_expression()
	&& upper_bound_.expression() == lower_bound_.expression())
      exact_solution_.set_expression(upper_bound_.has_expression());
    return SUCCESS;
  }
  else
    return TOO_COMPLEX;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::approximate_functional_equation_upper() const {
  if (functional_eq_p->rank() == 1) {
    Number divisor_arg = functional_eq_p->ht_begin()->first;
    Expr coefficient = functional_eq_p->ht_begin()->second;    
    
    if (!known_class_of_functional_eq_rank_1(coefficient, divisor_arg,
					     inhomogeneous_term))
      return TOO_COMPLEX;

    Number coeff = coefficient.ex_to_number();
    Expr q_upper = log(n) / log(divisor_arg);
    Expr upper_bound;

    Number condition = -1;    
    if (!inhomogeneous_term.is_zero()
	&& !compute_non_homogeneous_part(false, coeff, divisor_arg,
					 inhomogeneous_term, q_upper,
					 upper_bound, condition))
      return TOO_COMPLEX;

    if (condition > 1)
      set_applicability_condition(condition.to_unsigned());
    add_term_with_initial_condition(false, q_upper, upper_bound);
    
    upper_bound_.set_expression(simplify_logarithm(upper_bound));
    if (lower_bound_.has_expression()
	&& upper_bound_.expression() == lower_bound_.expression())
      exact_solution_.set_expression(lower_bound_.has_expression());
    return SUCCESS;
  }
  else
    return TOO_COMPLEX;
}
