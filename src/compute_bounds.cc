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

/*!
  Let \f$ x(n) \f$ be a sequence satisfying the generalized recurrenceg
  \f$ x_n = a x_{n/b} + g(n) \f$ with \f$ g(n) = c^n \f$ and \f$ c > 1 \f$.
  ...
*/
void
bounds_for_exp_function(const Number& coeff, const Number& divisor,
			const Number& base, Expr& lb, Expr& ub) {
  // Lower bound.
  lb += pwr(base, Recurrence::n);

  // Upper bound.
  Expr tmp = pwr(base, pwr(divisor, 2));
  Expr constant = coeff * tmp * pwr(tmp - pwr(coeff, divisor), -1);
  ub += pwr(base, Recurrence::n)
    + constant * pwr(base, Recurrence::n / divisor);
}

/*!
  Let \f$ x(n) \f$ be a sequence satisfying the generalized recurrence
  \f$ x_n = a x_{n/b} + g(n) \f$ with \f$ g(n) = log(n) \f$ and
  \f$ a >= 1 \f$. Then for all \f$ n >= 1 \f$ we have
  \f[
    x(n) \le \begin{cases}
               \dfrac12 \log n \,
	       \Bigl( \dfrac{\log n}{\log(b)} + 1 \Bigr)
	         & \text{if \f$ a = 1 \f$;} \\
	       \dfrac{1}{a - 1} \log(n) \,
	         \bigl( n^{(\log(a))/\log(b)} - 1 \bigr)
	         + \\
	       \qquad +
	         \dfrac{\log(b)}{(a - 1)^2}
	         \bigl( (a + 1) n^{(\log(a))/\log(b)} - a \bigr)
	         & \text{if \f$ a > 1 \f$.} \\
	     \end{cases}
  \f]
*/
void
upper_bound_for_log_function(const Number& coeff, const Number& divisor,
			     Expr& ub) {
  if (coeff == 1)
    ub += Number(1, 2) * log(Recurrence::n)
      * (log(Recurrence::n) / log(divisor) + 1);
  else {
    assert(coeff > 1);
    Expr ratio_log = log(coeff) / log(divisor);
    ub += (1 / (coeff - 1)) * log(Recurrence::n)
      * (pwr(Recurrence::n, ratio_log) - 1)
      + log(divisor) * pwr(coeff - 1, -2)
      * ((coeff + 1) * pwr(Recurrence::n, ratio_log) - coeff);
  }
}

/*!
  Let \f$ x(n) \f$ be a sequence satisfying the generalized recurrence
  \f$ x_n = a x_{n/b} + g(n) \f$ with \f$ g(n) = n^k log(n) \f$,
  where \f$ k \in \Rset \f$, with \f$ k \ge 0 \f$.
  Then for all \f$ n \ge 1 \f$ have
  \f[
    x(n) \le \begin{cases}
               \dfrac{b^k}{b^k - a} \, n^k \log(n)
                 & \text{if \f$ a < b^k \f$;} \\
	       n^k \, \dfrac{(\log n)^2}{\log(b)}
                 & \text{if \f$ a = b^k \f$;} \\
	       \dfrac{b^k}{a - b^k} \, n^k
	         \biggl( \Bigl( \dfrac{a}{b^k}
		         \Bigr)^{(\log(n))/\log(b)} - 1
		 \biggl) \log n & \text{if \f$ a > b^k \f$.} \\
	     \end{cases}
  \f]
*/
void
upper_bound_for_log_times_power_function(const Number& coeff,
					 const Number& divisor,
					 const Number& k, Expr& ub) {
  Number div_k = pwr(divisor, k);
  Expr tmp = pwr(Recurrence::n, k) * log(Recurrence::n);
  if (coeff < div_k)
    ub += (div_k * tmp) / (div_k - coeff);
  else if (coeff == div_k)
    ub += tmp * log(Recurrence::n) / log(divisor);
  else
    ub += tmp * (pwr(coeff / div_k, log(Recurrence::n) / log(divisor)) - 1)
      * div_k / (coeff - div_k);
}

/*!
  Let \f$ x(n) \f$ be a sequence satisfying the generalized recurrence
  \f$ x_n = a x_{n/b} + g(n) \f$ with \f$ g(n) = n^k \f$,
  where \f$ k \in \Rset \f$, with \f$ k \ge 0 \f$.
  Then for all $n \ge 1$ have
  \f[
    x(n) \le \begin{cases}
               \dfrac{b^k}{b^k - a} \, n^k
	         & \text{if \f$ a < b^k \f$;} \\
	       n^k \, \dfrac{\log n}{\log b}
	         & \text{if \f$ a = b^k \f$;} \\
	       \dfrac{b^k}{a - b^k} \, n^k
	         \biggl( \Bigl( \dfrac{a}{b^k}
		         \Bigr)^{(\log n)/\log b} - 1
		 \biggl)      & \text{if \f$ a > b^k \f$.} \\
	     \end{cases}
  \f]
*/
void
upper_bound_for_power_function(const Number& coeff, const Number& divisor,
			       const Number& k, Expr& ub) {
  Number div_k = pwr(divisor, k);
  Expr tmp = pwr(Recurrence::n, k);
  if (coeff < div_k)
    ub += tmp * div_k / (div_k - coeff);
  else if (coeff == div_k)
    ub += tmp * log(Recurrence::n) / log(divisor);
  else
    ub += div_k * (pwr(Recurrence::n, log(coeff) / log(divisor)) - tmp)
      / (coeff - div_k);
}

/*!
  Let \f$ x(n) \f$ be a sequence satisfying the generalized recurrence
  \f$ x_n = a x_{n/b} + g(n) \f$ with \f$ g(n) = n^k \f$,
  where \f$ k \in \Nset \f$.
  Then for all \f$ n \ge 1 \f$ have
  \f[
    x(n) \ge \begin{cases}
               \dfrac{n - a + 1}{(a-1)^2}
                 & \text{if \f$ k = 0 \f$ and \f$ a = b \f$;} \\
               n \Bigl(\dfrac{\log n}{\log a} - 1 \Bigr)
                 & \text{if \f$ k = 1 \f$ and \f$ a = b \f$;} \\
               \dfrac{a^{k-1}}{a^{k-1} - 1}
                 \Bigl( n^k - (a - 1) n^{k-1} \Bigr)
                   & \text{if \f$ k > 1 \f$ and \f$ a = b \f$.}
	     \end{cases}
  \f]
*/
void
lower_bound_for_power_function(const Number& coeff, const Number& k,
			       Expr& lb) {
  assert(k == 0 || k == 1 || k > 1);
  if (k == 0)
    lb += (Recurrence::n - coeff + 1) * pwr(coeff - 1, -2);
  else if (k == 1)
    lb += Recurrence::n * (log(Recurrence::n) / log(coeff) - 1);
  else if (k > 1)
    lb += pwr(coeff, k - 1)
      * (pwr(Recurrence::n, k) - (coeff - 1) * pwr(Recurrence::n, k - 1))
      / (pwr(coeff, k-1) - 1);
}

void
compute_bounds_for_power_of_n(const Number& coeff, const Number& divisor,
			      const Number& num, const Number& k,
			      bool& computed_only_upper,
			      Expr& lb, Expr& ub) {
  Expr tmp_ub = 0;
  upper_bound_for_power_function(coeff, divisor, k, tmp_ub);
  ub += num * tmp_ub;
  if (coeff == divisor && (k == 0 || k >= 1)) {
    Expr tmp_lb = 0;
    lower_bound_for_power_function(coeff, k, tmp_lb);
    lb += num * tmp_lb; 
  }
  else
    computed_only_upper = true;
}

/*!
  a * n^k
*/
bool
sharper_bounds_for_power_of_n(const Expr& poly_coeff,
			      const Number& coeff, const Number& divisor,
			      bool& computed_only_upper, Expr& lb, Expr& ub) {
  bool computed_at_least_one_bound = false;
  Number k;
  // Case `g(n) = a n^0' with `a' positive number.
  Number num;
  if (poly_coeff.is_a_number(num) && num.is_positive()) {
    compute_bounds_for_power_of_n(coeff, divisor, num, 0,
				  computed_only_upper, lb, ub);
    computed_at_least_one_bound = true;
  }
  // Case `g(n) = n^1'.
  else if (poly_coeff == Recurrence::n) {
    compute_bounds_for_power_of_n(coeff, divisor, 1, 1,
				  computed_only_upper, lb, ub);
    computed_at_least_one_bound = true;
  }
  // Case `g(n) = n^k' with `k != 1'.
  else if (poly_coeff.is_a_power()
	   && poly_coeff.arg(0) == Recurrence::n
	   && poly_coeff.arg(1).is_a_number(k)
	   && k.is_positive()) {
    compute_bounds_for_power_of_n(coeff, divisor, 1, k,
				  computed_only_upper, lb, ub);
    computed_at_least_one_bound = true;
  }
  else if (poly_coeff.is_a_mul() && poly_coeff.nops() == 2) {
    const Expr& first = poly_coeff.op(0);
    const Expr& second = poly_coeff.op(1);
    // Case `g(n) = a * n^k', with `a' and `k' positive numbers.
    Number num;
    if (first.is_a_number(num) && num.is_positive()) {
      if (second == Recurrence::n) {
	compute_bounds_for_power_of_n(coeff, divisor, num, 1,
				      computed_only_upper, lb, ub);
	computed_at_least_one_bound = true;
      }
      Number k;
      if (second.is_a_power()
	  && second.arg(0) == Recurrence::n
	  && second.arg(1).is_a_number(k)
	  && k.is_positive()) {
	compute_bounds_for_power_of_n(coeff, divisor, num, k,
				      computed_only_upper, lb, ub);
	computed_at_least_one_bound = true;
      }
    }
    else if (second.is_a_number(num) && num.is_positive()) {
      if (first == Recurrence::n) {
	compute_bounds_for_power_of_n(coeff, divisor, num, 1,
				      computed_only_upper, lb, ub);
	computed_at_least_one_bound = true;
      }
      Number k;
      if (first.is_a_power()
	  && first.arg(0) == Recurrence::n
	  && first.arg(1).is_a_number(k)
	  && k.is_positive()) {
	compute_bounds_for_power_of_n(coeff, divisor, num, k,
				      computed_only_upper, lb, ub);
	computed_at_least_one_bound = true;
      }
    }
  }
  return computed_at_least_one_bound;
}

void
compute_bounds_for_logarithm_function(const Number& coeff,
				      const Number& divisor,
				      const Number& num,
				      Expr& ub) {
  Expr tmp_ub = 0;
  upper_bound_for_log_function(coeff, divisor, tmp_ub);
  ub =+ num * tmp_ub;
}

void
compute_bounds_for_power_times_logarithm_function(const Number& coeff,
						  const Number& divisor,
						  const Number& num,
						  const Number& k,
						  Expr& ub) {
  Expr tmp_ub = 0;
  upper_bound_for_log_times_power_function(coeff, divisor, k, tmp_ub);
  ub += num * tmp_ub;
}

/*!
  g(n) = a log(n) or g(n) = a * n^k log(n).
*/
bool
sharper_bounds_for_no_polynomial_function(const Expr& no_poly_coeff,
					  const Number& coeff,
					  const Number& divisor,
					  bool& computed_only_upper,
					  Expr& /*lb*/, Expr& ub) {
  bool computed_at_least_one_bound = false;
  // Case `g(n) = log(n)'.
  if (no_poly_coeff.is_the_log_function()
      && no_poly_coeff.arg(0) == Recurrence::n) {
    compute_bounds_for_logarithm_function(coeff, divisor, 1, ub);
    computed_only_upper = true;
    computed_at_least_one_bound = true;
  }
  else if (no_poly_coeff.is_a_mul()) {
    if (no_poly_coeff.nops() == 2) {
      const Expr& first = no_poly_coeff.op(0);
      const Expr& second = no_poly_coeff.op(1);
      if (first == log(Recurrence::n)) {
	// Case `g(n) = a log(n)' with `a' positive number.
	Number num;
	if (second.is_a_number(num) && num.is_positive()) {
	  compute_bounds_for_logarithm_function(coeff, divisor, num, ub);
	  computed_only_upper = true;
	  computed_at_least_one_bound = true;
	}
	// Case `g(n) = n^k log(n)' with `k >= 1'.
	Number k;
	if (second == Recurrence::n) {
	  compute_bounds_for_power_times_logarithm_function(coeff, divisor,
							    1, 1, ub);
	  computed_only_upper = true;
	  computed_at_least_one_bound = true;
	}
	else if (second.is_a_power() && second.arg(0) == Recurrence::n
		 && second.arg(1).is_a_number(k) && k.is_positive()) {
	  compute_bounds_for_power_times_logarithm_function(coeff, divisor,
							    1, k, ub);
	  computed_only_upper = true;
	  computed_at_least_one_bound = true;
	}
      }
      if (second == log(Recurrence::n)) {
	// Case `g(n) = a log(n)' with `a' positive number.
	Number num;
	if (first.is_a_number(num) && num.is_positive()) {
	  compute_bounds_for_logarithm_function(coeff, divisor, num, ub);
	  computed_only_upper = true;
	  computed_at_least_one_bound = true;
	}
	// Case `g(n) = n^k log(n)' with `k >= 1'.
	Number k;
	if (first == Recurrence::n) {
	  compute_bounds_for_power_times_logarithm_function(coeff, divisor,
							    1, 1, ub);
	  computed_only_upper = true;
	  computed_at_least_one_bound = true;
	}
	else if (first.is_a_power() && first.arg(0) == Recurrence::n
		 && first.arg(1).is_a_number(k) && k.is_positive()) {
	  compute_bounds_for_power_times_logarithm_function(coeff, divisor,
							    1, k, ub);
	  computed_only_upper = true;
	  computed_at_least_one_bound = true;
	}
      }
    }
    // Case `g(n) = a n^k log(n)' with `a' positive number and `k >= 1'.
    if (no_poly_coeff.nops() == 3) {
      D_MSGVAR("qui ", no_poly_coeff);
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
      if (power_part == Recurrence::n) {
	compute_bounds_for_power_times_logarithm_function(coeff, divisor,
							  num_part, 1, ub);
	computed_only_upper = true;
	computed_at_least_one_bound = true;
      }
      // `k != 1'.
      Number k;
      if (power_part.is_a_power() && power_part.arg(0) == Recurrence::n
	  && power_part.op(1).is_a_number(k) && k.is_positive()) {
	compute_bounds_for_power_times_logarithm_function(coeff, divisor,
							  num_part, k, ub);
	computed_only_upper = true;
	computed_at_least_one_bound = true;
      }
    }
  }
  return computed_at_least_one_bound;
}

/*!  
  g(n) = c^n, c > 1
*/
bool
sharper_bounds_for_exponential(const Number& base, const Expr& poly_coeff,
			       const Number& coeff, const Number& divisor,
			       Expr& lb, Expr& ub) {
  bool computed_at_least_one_bound = false;
  Number num;
  if (poly_coeff.is_a_number(num) && num.is_positive())
    if (coeff >= 1
	&& pwr(base, pwr(divisor, 2)) > coeff * pwr(base, divisor)) {
      Expr tmp_ub = 0;
      Expr tmp_lb = 0;
      bounds_for_exp_function(coeff, divisor, base, tmp_lb, tmp_ub);
      ub += num * tmp_ub;
      lb += num * tmp_lb;
      computed_at_least_one_bound = true;
    }
  return computed_at_least_one_bound;
}

/*!
  If \p for_lower is true then this function compute the sum
  \f$ G(n) = sum_{k = 1}^n a^{n-k} p(b^k) \f$, otherwise
  compute the sum \f$ H(n) = sum_{k = 2}^n a^{n-k} p(b^k-1) \f$,
  where \f$ a \f$ is stored in \p coefficient, \f$ b \f$ is stored
  in \p divisor and \f$ p(n) \f$ is stored in \p summand.
*/
bool
compute_sum(const Expr& summand,
	    const Number& coefficient, const Number& divisor,
	    bool for_lower, Expr& sum) {
  D_MSG("***compute sum");
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
  exp_poly_decomposition(tmp, bases_of_exp,
			 exp_poly_coeff, exp_no_poly_coeff);
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
      if (!full_gosper(pwr(coefficient, -Recurrence::n) * tmp,
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
  If we have already computed the upper bound computes only
  \f$ G(n) \f$.
*/
void
try_to_compute_sum(const Expr& summand,
		   const Number& coefficient, const Number& divisor,
		   bool computed_only_upper,
		   Expr& sum_lower, Expr& sum_upper) {
  D_VAR(summand);
  // Compute `G(n)' used for the lower bound.
  if (!compute_sum(summand, coefficient, divisor, true, sum_lower))
    // We try to compute `H(n)' because `H(q) <= G(q)'.
    if (!compute_sum(summand, coefficient, divisor, false, sum_lower)) {
      Symbol h;
      sum_lower
	= PURRS::sum(h, 1, Recurrence::n, pwr(coefficient, Recurrence::n-h)
		     * summand.substitute(Recurrence::n,
					     pwr(divisor, h)));
    }
  
  if (!computed_only_upper) {
    // Compute `H(n)' used for the upper bound.
    if (!compute_sum(summand, coefficient, divisor, false, sum_upper))
      // We try to compute `G(n)' because `H(q+1) <= G(q+1)'.
      if (!compute_sum(summand, coefficient, divisor, true, sum_upper)) {
	Symbol h;
	sum_upper = PURRS::sum(h, 2, Recurrence::n,
			       pwr(coefficient, Recurrence::n-h)
			       * summand.substitute(Recurrence::n,
						    pwr(divisor, h) - 1));
      }
  }
}
					     
} // anonymous namespace

/*!
  Computes
  \f[
    \a^q x \left( \frac n{\b^q} \right)
  \f]
  where \f$ a \f$ is stored in <CODE>coefficient()</CODE>, \f$ b \f$
  is stored <CODE>divisor_arg()</CODE>.
  \f$ q \f$ is stored in \p q_upper if the value computed will be
  added to \p ub or it is stored in \p q_lower if the value computed
  will be added to \p lb.
*/
void
PURRS::Recurrence::add_term_with_initial_condition(const Expr& q_upper,
						   const Expr& q_lower,
						   Expr& ub, Expr& lb) const {
  Expr index_initial_condition
    = simplify_logarithm(n / pwr(divisor_arg(), q_upper));
  Number index = 0;
  if (index_initial_condition.is_a_number()) {
    index = index_initial_condition.ex_to_number();
    D_VAR(index);
    D_VAR(applicability_condition());
    if (index < applicability_condition())
      index = applicability_condition();
    assert(index.is_positive_integer());
  }
  // The index of the initial condition is not a number.
  if (index == 0) {
    ub += pwr(coefficient(), q_upper) * x(index_initial_condition);
    lb += pwr(coefficient(), q_lower) * x(index_initial_condition);
  }
  // The index of the initial condition is a number.
  else {
    ub += pwr(coefficient(), q_upper)
      * get_initial_condition(index.to_unsigned());
    lb += pwr(coefficient(), q_lower)
      * get_initial_condition(index.to_unsigned());
  }
}

/*
  This function find a lower bound and an upper bound for special recurrences
  that we call <EM>functional equations</EM>.
  For the moment we consider a limited part of these recurrences, i. e.
  equations of the form
  \f[
    x_n = a x_{n/b} + p(n)
  \f]
  where \f$ a > 0 \f$, \f$ b \in \Nset \ setminus \{0, 1\} \f$ and
  \f$ p : \Nset \setminus \{ 0 \} \rightarrow \Rset \f$.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::approximate_functional_equation() const {
  D_VAR(coefficient());
  D_VAR(divisor_arg());
  assert(divisor_arg().is_rational() && divisor_arg().is_positive());
  // We want that `coefficient()' is a positive number and
  // `divisor_arg()' a rational number bigger than `1'.
  Number coeff;
  if (!coefficient().is_a_number(coeff) || !coeff.is_positive()
      || divisor_arg() < 1)
    return TOO_COMPLEX;
  if (has_parameters(inhomogeneous_term)) {
    D_MSG("Functional equation with parameters");
    return TOO_COMPLEX;
  }

  // Consider an upper bound and a lower bound for `q = [log n / log b]'.
  Expr q_upper = log(n) / log(divisor_arg());
  Expr q_lower = q_upper - 1;
  Expr ub = 0;
  Expr lb = 0;
  D_VAR(inhomogeneous_term);
  if (!inhomogeneous_term.is_zero()) {
    std::vector<Expr> bases_of_exp;
    std::vector<Expr> exp_poly_coeff;
    std::vector<Expr> exp_no_poly_coeff;
    exp_poly_decomposition(inhomogeneous_term,
			   bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
    D_VEC(bases_of_exp, 0, bases_of_exp.size()-1);
    D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
    D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);
    for (unsigned i = bases_of_exp.size(); i-- > 0; ) {
      bool computed_only_upper = false;
      const Expr& base = bases_of_exp[i];
      const Expr& poly_coeff = exp_poly_coeff[i];
      const Expr& no_poly_coeff = exp_no_poly_coeff[i];
      Number num_base;
      if (base.is_a_number(num_base) && num_base.is_positive_integer()) {
	// Base of exponential is equal to `1'.
	if (num_base == 1) {
	  // Consider the eventual polynomial part.
	  if (!poly_coeff.is_zero()
	      && (!sharper_bounds_for_power_of_n(poly_coeff, 
						 coeff, divisor_arg(),
						 computed_only_upper, lb, ub)
	      || computed_only_upper)) {
	    // Check if the polynomial part is a non-negative,
	    // non-decreasing function.
	    Number condition = -1;
	    if (!is_non_negative_non_decreasing(poly_coeff, 1,
						poly_coeff, true,
						n, condition))
	      return TOO_COMPLEX ;
	    if (condition > 1)
	      set_applicability_condition(condition.to_unsigned());
	    Expr sum_lower = 0;
	    Expr sum_upper = 0;
	    try_to_compute_sum(poly_coeff, coeff, divisor_arg(),
			       computed_only_upper, sum_lower, sum_upper);
	    lb += sum_lower.substitute(n, q_lower);
	    ub += sum_upper.substitute(n, q_upper + 1);
	  }
	  // Consider the eventual non-polynomial part.
	  if (!no_poly_coeff.is_zero()
	      && (!sharper_bounds_for_no_polynomial_function(no_poly_coeff,
							     coeff,
							     divisor_arg(),
							     computed_only_upper,
							     lb, ub)
		  || computed_only_upper)) {
	    // Check if the non-polynomial part is a non-negative,
	    // non-decreasing function.
	    Number condition = -1;
	    if (!is_non_negative_non_decreasing(no_poly_coeff, 1,
						no_poly_coeff, false,
						n, condition))
	      return TOO_COMPLEX ;
	    if (condition > 1)
	      set_applicability_condition(condition.to_unsigned());
	    Expr sum_lower = 0;
	    Expr sum_upper = 0;
	    try_to_compute_sum(no_poly_coeff, coeff, divisor_arg(),
			       computed_only_upper, sum_lower, sum_upper);
	    lb += sum_lower.substitute(n, q_lower);
	    ub += sum_upper.substitute(n, q_upper + 1);
	  }
	}
	// Base of exponential is a number different from `1'.
	else {
	  if (!poly_coeff.is_zero())
	    // Apply, if it possible, the sharper bounds for the exponential.
	    if (!sharper_bounds_for_exponential(num_base, poly_coeff,
						coeff, divisor_arg(),
						lb, ub)) {
	      // Check if the polynomial part times esponential is a
	      // non-negative, non-decreasing function.
	      Number condition = -1;
	      if (!is_non_negative_non_decreasing(pwr(base, n) * poly_coeff,
						  num_base, poly_coeff, true,
						  n, condition))
		return TOO_COMPLEX ;
	      if (condition > 1)
		set_applicability_condition(condition.to_unsigned());
	      Expr sum_lower = 0;
	      Expr sum_upper = 0;
	      try_to_compute_sum(pwr(base, n) * poly_coeff,
				 coeff, divisor_arg(),
				 computed_only_upper, sum_lower, sum_upper);
	      lb += sum_lower.substitute(n, q_lower);
	      ub += sum_upper.substitute(n, q_upper + 1);
	    }
	  if (!no_poly_coeff.is_zero()) {
	    // In this case is not possible to apply the sharper
	    // bounds for the exponential.
	    // Check if the polynomial part times esponential is a
	    // non-negative, non-decreasing function.
	    Number condition = -1;
	    if (!is_non_negative_non_decreasing(pwr(base, n) * no_poly_coeff,
						num_base, no_poly_coeff, false,
						n, condition))
	      return TOO_COMPLEX ;
	    if (condition > 1)
	      set_applicability_condition(condition.to_unsigned());
	    Expr sum_lower = 0;
	    Expr sum_upper = 0;
	    try_to_compute_sum(pwr(base, n) * no_poly_coeff,
			       coeff, divisor_arg(),
			       computed_only_upper, sum_lower, sum_upper);
	    lb += sum_lower.substitute(n, q_lower);
	    ub += sum_upper.substitute(n, q_upper + 1);
	  }
	}
      }
      // The base of the exponential is not a number.
      else
	return TOO_COMPLEX;
    }  
  }
  add_term_with_initial_condition(q_upper, q_lower, ub, lb);

  upper_bound_.set_expression(simplify_logarithm(ub));
  lower_bound_.set_expression(simplify_logarithm(lb));
  return SUCCESS;
}
