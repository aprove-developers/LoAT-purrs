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
#include "functional_equation.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Cached_Expr.defs.hh"
#include "Recurrence.defs.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

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
  // We already know that `b' is a positive integer and we want also
  // that `a' is a positive number.
  Number coeff;
  if (!coefficient().is_a_number(coeff) || !coeff.is_positive())
    return TOO_COMPLEX;
  // Check if the inhomogeneous term `p(n)' is a non-negative,
  // non-decreasing function. For to do this, the parameters are
  // not allowed.
  if (find_parameters(inhomogeneous_term)
      || !is_non_negative_non_decreasing(inhomogeneous_term, n))
    return TOO_COMPLEX;
  // Compute `sum_{k = 1}^n a^{n-k} p(b^k)'.
  std::vector<Expr> bases_of_exp;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  Expr tmp = inhomogeneous_term.substitute(n, pwr(divisor_arg(), n));
  exp_poly_decomposition(simplify_ex_for_input(tmp, true),
			 bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
  assert(coefficient().is_a_number());
  Expr sum = 0;
  if (vector_not_all_zero(exp_poly_coeff)
      && vector_not_all_zero(exp_no_poly_coeff))
    return TOO_COMPLEX;
  else if (vector_not_all_zero(exp_poly_coeff))
    for (unsigned i = bases_of_exp.size(); i-- > 0; ) {
      Symbol k("k");
      sum += sum_poly_times_exponentials(exp_poly_coeff[i].substitute(n, k), k,
					 bases_of_exp[i] / coefficient());
      // `sum_poly_times_exponentials' computes the sum from 0, whereas
      // we want that the sum start from `1'.
      sum -= exp_poly_coeff[i].substitute(n, 0);
    }
  else
    if (bases_of_exp.size() == 1) {
      Expr& no_poly = exp_no_poly_coeff[0];
      assert(!no_poly.is_a_add());
      if (no_poly.is_a_mul()) {
	for (unsigned i = no_poly.nops(); i-- > 0; )
	  if (!no_poly.op(i).is_polynomial(n)
	      && !no_poly.op(i).is_the_log_function())
	    return TOO_COMPLEX;
      }
      else
	if (!no_poly.is_polynomial(n) && !no_poly.is_the_log_function())
	  return TOO_COMPLEX;
      Symbol k("k");
      no_poly = simplify_logarithm(no_poly);
      sum += sum_poly_times_exponentials(no_poly.substitute(n, k), k,
					 bases_of_exp[0] / coefficient());
      // `sum_poly_times_exponentials' computes the sum from 0, whereas
      // we want that the sum start from `1'.
      sum -= no_poly.substitute(n, 0); 
    }
    else
      return TOO_COMPLEX;
  sum *= pwr(coefficient(), n);
  sum = simplify_ex_for_output(sum, false);
  D_VAR(sum);
  // Consider an upper bound and a lower bound for `q = [log n / log b]'.
  Expr q_upper = log(n) / log(divisor_arg());
  Expr q_lower = q_upper - 1;

  Expr index_initial_condition
    = simplify_logarithm(n / pwr(divisor_arg(), q_upper));
  Number index = 0;
  if (index_initial_condition.is_a_number()) {
    index = index_initial_condition.ex_to_number();
    assert(index.is_positive_integer());
  }

  Expr ub;
  Expr lb;
  if (index == 0) {
    ub = pwr(coefficient(), q_upper) * index_initial_condition
      + sum.substitute(n, q_upper + 1);
    lb = pwr(coefficient(), q_lower) * index_initial_condition
      + sum.substitute(n, q_lower);
  }
  else {
    ub = pwr(coefficient(), q_upper) * get_initial_condition(index.to_int())
      + sum.substitute(n, q_upper + 1);
    lb = pwr(coefficient(), q_lower) * get_initial_condition(index.to_int())
      + sum.substitute(n, q_lower);
  }
  upper_bound_.set_expression(simplify_logarithm(ub));
  lower_bound_.set_expression(simplify_logarithm(lb));

  return SUCCESS;
}
