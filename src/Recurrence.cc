/* Recurrence class implementation (non-inline functions).
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

#include <config.h>

#include "Recurrence.defs.hh"
#include "util.hh"
#include "simplify.hh"
#include "factorize.hh"
#include "finite_order.hh"
#include "Expr.defs.hh"
#include "Cached_Expr.defs.hh"
#include "Non_Linear_Info.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Blackboard.defs.hh"
#include <algorithm>
#include <iostream>
#include <fstream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

void
split(const Expr& e, const Expr& d, Expr& term_with_d, Expr& other_terms) {
  assert(e.is_a_add());
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& term = e.op(i);
    if (term.has(d))
      term_with_d += term;
    else
      other_terms += term;
  }
}

bool
ok_inequalities(const Expr& e, unsigned condition) {
  assert(e.is_a_add());
  Expr term_with_n = 0;
  Expr other_terms = 0;
  split(e, Recurrence::n, term_with_n, other_terms);
  D_VAR(other_terms);
  if (term_with_n == Recurrence::n || term_with_n.is_a_mul()) {
    Expr coeff_n = 1;
    if (term_with_n.is_a_mul()) {
      for (unsigned i = term_with_n.nops(); i-- > 0; ) {
	const Expr& factor = term_with_n.op(i);
	Number num;
	if (!(factor == Recurrence::n || (factor.is_a_power()
					  && factor.arg(0) == Recurrence::n
					  && factor.arg(1).is_a_number(num)
					  && num.is_positive_integer())))
	  coeff_n *= factor;
      }
      D_VAR(coeff_n); 
    }
    Number numer;
    Number denom;
    if (coeff_n.is_a_number(denom) && denom.is_positive()
	&& other_terms.is_a_number(numer) && numer.is_negative())
      if (-numer/denom <= condition)
	return true;
    //       // Devo verificare le altre condizioni iniziali...
    // 	else if () {
    // 	}
  }
  return false;
}

bool
validation_initial_conditions_in_bound(bool upper, const Expr& bound,
				       unsigned index) {
  Expr bound_valuated = bound.substitute(Recurrence::n, index);
  D_VAR(bound_valuated);
  if (bound_valuated != x(index))
    if (bound_valuated.is_a_mul()) {
      Expr coeff_ic = 1;
      for (unsigned i = bound_valuated.nops(); i-- > 0; ) {
	const Expr& factor = bound_valuated.op(i);
	if (factor != x(index))
	  coeff_ic *= factor;
      }
      D_VAR(coeff_ic);
      Number num;
      if (upper) {
	if (coeff_ic.is_a_number(num) && num < 1)
	  return false;
      }
      else
	if (coeff_ic.is_a_number(num) && num > 1)
	  return false;
    }
    else if (bound_valuated.is_a_add()) {
      Expr term_with_ic = 0;
      Expr other_terms = 0;
      split(bound_valuated, x(index),
	    term_with_ic, other_terms);
      Expr coeff_ic = 1;
      if (term_with_ic.is_a_mul()) {
	for (unsigned i = term_with_ic.nops(); i-- > 0; ) {
	  const Expr& factor = term_with_ic.op(i);
	  if (factor != x(index))
	    coeff_ic *= factor;
	}
      }
      D_VAR(coeff_ic);
      D_VAR(other_terms);
      Number num_coeff;
      Number num_other;
      if (coeff_ic.is_a_number(num_coeff)
	  && other_terms.is_a_number(num_other)) {
	if (upper) {
	  if (!(num_coeff >= 1 && num_other.is_positive()))
	    return false;
	}
	else
	  if (!(num_coeff <= 1 && num_other.is_negative()))
	    return false;
      }
      else
	return false;
    }
  return true;
}

} // anonymous namespace

const PURRS::Symbol&
PURRS::Recurrence::n = Symbol("n");

/*!
  Consider the right hand side \p rhs of the order \f$ k \f$ recurrence
  relation
  \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$.
  Let \f$ i \f$ the index for the first initial condition starting from
  which the recurrence is well-defined.
  We try to check that the solution is correct.
  - Validation of initial conditions.
    If <CODE>recurrence_rhs</CODE> is equal to \f$ x(i), \cdots, x(i+k) \f$ for
    \f$ n = i, \cdots, i+k-1 \f$ respectively then
    the initial conditions are verified and we continue to check; otherwise
    return <CODE>false</CODE> because the solution can be wrong or it is not
    simplified enough.
  - Since the initial conditions are verified, we erase from
    <CODE>solution</CODE> all terms containing an initial condition.
    In other words, we check that ther remainder of the solution
    satisfies the same recurrence relation, but with the initial conditions
    all equal to \f$ 0 \f$.
    Starting from the partial solution just computed, we substitute
    \f$ x(n-1) \f$, \f$ x(n-2) \f$, \f$ \dots \f$, \f$ x_{n-k} \f$ into
    <CODE>recurrence_rhs</CODE>.
    We next consider the difference between the partial solution
    and the new right hand side:
    - if it is equal to zero -> returns <CODE>PROVABLY_CORRECT</CODE>:
                                the solution is certainly right.
    - if it is not equal to zero (in a syntactical sense)
                             -> returns <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
			        the solution can be wrong or
				we failed to simplify it.
  FIXME: In the latter case, we will need more powerful tools to
  decide whether the solution is right or it is really wrong and, in this last
  case, to return <CODE>PROVABLY_INCORRECT</CODE>.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_exact_solution(const Recurrence& rec) {
  assert (rec.is_linear_finite_order() || rec.is_non_linear_finite_order());
  unsigned int order_rec;
  unsigned int first_i_c;
  if (rec.is_non_linear_finite_order()) {
    order_rec = rec.order_if_linear();
    first_i_c = rec.non_linear_to_linear_fwdr();
  }
  else {
    order_rec = rec.order();
    first_i_c = rec.first_well_defined_rhs_linear();
  }
  
  D_VAR(rec.recurrence_rhs);
  D_VAR(rec.exact_solution_.expression());
  D_VAR(order_rec);
  D_VAR(first_i_c);
  
  if (order_rec == 0)
    return PROVABLY_CORRECT;
  else {
    // Step 1: validation of initial conditions.
    for (unsigned i = order_rec; i-- > 0; ) {
      Expr solution_evaluated
	= rec.exact_solution_.expression().substitute(n, first_i_c + i);
      solution_evaluated = rec.blackboard.rewrite(solution_evaluated);
      solution_evaluated = simplify_all(solution_evaluated);
      D_VAR(solution_evaluated);
      D_VAR(x(first_i_c + i));
      if (solution_evaluated != x(first_i_c + i))
	return INCONCLUSIVE_VERIFICATION;
    }
    // Step 2: find `partial_solution'.
    // The initial conditions are verified. Build the expression
    // `partial_solution' that has all terms of `solution' minus those
    // containing an initial condition.
    Expr partial_solution = 0;
    if (rec.exact_solution_.expression().is_a_add())
      for (unsigned i = rec.exact_solution_.expression().nops(); i-- > 0; ) {
	if (!rec.exact_solution_.expression().op(i).has_x_function(true))
	  partial_solution += rec.exact_solution_.expression().op(i);
      }
    else
      if (!rec.exact_solution_.expression().has_x_function(true))
	partial_solution = rec.exact_solution_.expression();
    D_VAR(partial_solution);
    // The recurrence is homogeneous.
    if (partial_solution == 0)
      return PROVABLY_CORRECT;
    // Step 3: construct the vector `terms_to_sub': each element of it
    // contains `partial_solution' with `n' substituted by `n - d'
    // (the `d' are the decrements of the terms `x(n - d)').
    // These new expressions contained in the vector `terms_to_sub' are
    // substituted to the correspondenting values in `recurrence_rhs'.
    Expr substituted_rhs;
    std::vector<Expr> terms_to_sub(order_rec);
    
    substituted_rhs = rec.recurrence_rhs;
    for (unsigned i = order_rec; i-- > 0; ) {
      terms_to_sub[i] = simplify_all(partial_solution.substitute
				     (n, n - (i + 1)));
      terms_to_sub[i] = simplify_sum(terms_to_sub[i], true);
      substituted_rhs = substituted_rhs
	.substitute(x(n - (i + 1)), terms_to_sub[i]);
    }
    D_VEC(terms_to_sub, 0, terms_to_sub.size()-1);
    D_VAR(substituted_rhs);
    Expr diff = rec.blackboard.rewrite(partial_solution - substituted_rhs);
    diff = simplify_all(diff);
    D_VAR(diff);
    if (!diff.is_zero())
      if (rec.applied_order_reduction) {
	rec.applied_order_reduction = false;
	// If we have applied the order reduction and we do not have
	// success in the verification of the original recurrence, then
	// we please ourselves if is verified the reduced recurrence.
	Symbol r
	  = rec.insert_auxiliary_definition(mod(n,
						rec.gcd_among_decrements()));
	unsigned dim = rec.coefficients().size()
	  / rec.gcd_among_decrements() + 1;
	std::vector<Expr> new_coefficients(dim);
	Expr inhomogeneous = 0;
	Recurrence rec_rewritten
	  (rewrite_reduced_order_recurrence(rec.recurrence_rhs, r,
					    rec.gcd_among_decrements(),
					    rec.coefficients(),
					    new_coefficients,
					    inhomogeneous));
	rec_rewritten.finite_order_p
	  = new Finite_Order_Info(dim - 1, new_coefficients, 1);
	rec_rewritten.set_type(rec.type());
	rec_rewritten.set_inhomogeneous_term(inhomogeneous);
	rec_rewritten.solve_linear_finite_order();
	D_VAR(rec.exact_solution_.expression());
	return verify_exact_solution(rec_rewritten);
      }
      else
	return INCONCLUSIVE_VERIFICATION;
    return PROVABLY_CORRECT;
  }
}

/*!
  Consider the right hand side \p rhs of the functional equation
  \f$ a x_{n/b} + p(n) \f$.
  If \p upper is <CODE>true</CODE> we try to check that the upper bound
  is correct;
  If \p upper is <CODE>false</CODE> we try to check that the lower bound
  is correct.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_bound(const Recurrence& rec, bool upper) {
  assert(rec.is_functional_equation());
  Expr bound;
  if (upper)
    bound = simplify_sum(rec.upper_bound_.expression(), true);
  else
    bound = simplify_sum(rec.lower_bound_.expression(), true);
  
  // Step 1: validation of initial conditions.
  if (!validation_initial_conditions_in_bound(upper, bound,
					      rec.applicability_condition()))
    return INCONCLUSIVE_VERIFICATION;
  
  // Step 2: find `partial_bound'.
  // We not consider the terms containing the initial conditions:
  // `partial_bound' will contain all the other terms.
  Expr partial_bound = 0;
  if (bound.is_a_add())
    for (unsigned i = bound.nops(); i-- > 0; ) {
      if (!bound.op(i).has_x_function(true))
	partial_bound += bound.op(i);
    }
  else
    if (!bound.has_x_function(true))
      partial_bound = bound;
  D_VAR(partial_bound);
  // The recurrence is homogeneous.
  if (partial_bound == 0)
    return PROVABLY_CORRECT;
  
  // Step 3: verification of the inductive base.
  Number num;
  if (upper && partial_bound
      .substitute(n, rec.applicability_condition()).is_a_number(num)
      && num.is_negative())
    return INCONCLUSIVE_VERIFICATION;
  if (!upper && partial_bound
      .substitute(n, rec.applicability_condition()).is_a_number(num)
      && num.is_positive())
    return INCONCLUSIVE_VERIFICATION;
  
  // Step 4: verification of the inductive step.
  Expr partial_bound_sub
    = partial_bound.substitute(n, n / rec.functional_eq_p->ht_begin()->first);
  Expr approx = rec
    .recurrence_rhs.substitute(x(n / rec.functional_eq_p->ht_begin()->first),
			       partial_bound_sub);
  D_VAR(approx);
  approx = simplify_ex_for_input(approx, true);
  approx = simplify_logarithm(approx);
  D_VAR(approx);
  
  Expr diff;
  if (upper)
    diff = partial_bound - approx;
  else
    diff = approx - partial_bound;
  D_VAR(diff);
  
  if (diff.is_a_number(num)) {
    if (num == 0 || num.is_positive())
      return PROVABLY_CORRECT;
  }
  else if (diff.is_a_mul()) {
    Expr coeff_n = 1;
    for (unsigned i = diff.nops(); i-- > 0; ) {
      const Expr& factor = diff.op(i);
      if (!(factor == n || (factor.is_a_power()
			    && factor.arg(0) == n
			    && factor.arg(1).is_a_number(num)
			    && num.is_positive_integer())))
	coeff_n *= factor;
    }
    D_VAR(coeff_n);
    if (coeff_n.is_a_number(num) && num.is_positive())
      return PROVABLY_CORRECT;
  }
  else if (diff.is_a_add())
    if (ok_inequalities(diff, rec.applicability_condition()))
      return PROVABLY_CORRECT;
  return INCONCLUSIVE_VERIFICATION;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::apply_order_reduction() const {
  applied_order_reduction = true;
  // Build the new recurrence substituting `n' not contained in the
  // `x' functions with `gcd_among_decrements * n + r' and `x(n-k)' with
  // `x(n - k / gcd_among_decrements)'.
  Symbol r = insert_auxiliary_definition(mod(n, gcd_among_decrements()));
  unsigned dim = coefficients().size() / gcd_among_decrements() + 1;
  std::vector<Expr> new_coefficients(dim);
  Expr inhomogeneous = 0;
  Recurrence rec_rewritten
    (rewrite_reduced_order_recurrence(recurrence_rhs, r,
				      gcd_among_decrements(),
				      coefficients(), new_coefficients,
				      inhomogeneous));
  rec_rewritten.finite_order_p
    = new Finite_Order_Info(dim - 1, new_coefficients, 1);
  rec_rewritten.set_type(type());
  rec_rewritten.set_inhomogeneous_term(inhomogeneous);
  Solver_Status status;
  if ((status = rec_rewritten.solve_linear_finite_order()) == SUCCESS) {
    // Now we must compute the solution for the original recurrence.
    // Perform three substitutions:
    // - r                      -> mod(n, gcd_among_decrements);
    // - n                      -> 1 / gcd_among_decrements
    //                             * (n - mod(n, gcd_among_decrements));
    // - x(k), k non-negative integer -> x(mod(n, gcd_among_decrements))
    //                                   + k * gcd_among_decrements.
    exact_solution_.set_expression(come_back_to_original_variable
				   (rec_rewritten
				    .exact_solution_.expression(), r,
				    get_auxiliary_definition(r),
				    gcd_among_decrements()));
    // We must copy the values in the blackboard of the reduced recurrence
    // in the blackboard of the original recurrences.
    blackboard = rec_rewritten.blackboard;
    // If there are defined initial conditions in the map
    // `initial_conditions' then we must to deduce the exact form
    // of the solution and substitute to it the defined initial
    // conditions.
    if (!initial_conditions.empty()) {
      exact_solution_.set_expression
	(simplify_ex_for_output(exact_solution_.expression(), false));
      Expr expanded_solution
	= write_expanded_solution(*this, gcd_among_decrements());
      exact_solution_.set_expression(expanded_solution);
    }
    exact_solution_.set_expression
      (simplify_ex_for_output(exact_solution_.expression(), false));
    recurrence_rhs_rewritten = true;
    return SUCCESS;
  }
  else
    return status;
}

/*!
  Builds a new object recurrence containing a linear recurrence and
  classify it.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::
compute_non_linear_recurrence(Expr& solution_or_bound, unsigned type) const {
  D_MSG("compute_non_linear_recurrence");
  // Build a new object recurrence with a linear recurrence.
  Recurrence rec_rewritten(rhs_transformed_in_linear());
  rec_rewritten.come_from_non_linear_rec = true;
  D_VAR(rec_rewritten.recurrence_rhs);

  Solver_Status status;
  // Classify the linear recurrence `rec_rewritten'.
  if (rec_rewritten.is_classified
      || (status = rec_rewritten.classify_and_catch_special_cases())
      == SUCCESS) {
    assert(rec_rewritten.is_linear_finite_order()
	   || rec_rewritten.is_functional_equation());

    // Linear finite order.
    if (rec_rewritten.is_linear_finite_order())
      if ((status = rec_rewritten.solve_linear_finite_order())
	  == SUCCESS) {
	set_order_if_linear(rec_rewritten.order());
	set_non_linear_to_linear_fwdr
	  (rec_rewritten.first_well_defined_rhs_linear());
	// Transform the solution of the linear recurrence in the solution
	// of the non linear recurrence.
	if (rec_rewritten.exact_solution_.expression() == 0)
	  solution_or_bound = 0;
	else {
	  solution_or_bound = pwr(base_exp_log(),
				  rec_rewritten.exact_solution_.expression());
	  solution_or_bound = substitute_x_function(solution_or_bound,
						    base_exp_log(), false);
	  solution_or_bound = simplify_ex_for_input(solution_or_bound, true);
	  solution_or_bound = simplify_logarithm(solution_or_bound);
	  // Resubstitute eventual auxiliary symbols with the respective
	  // negative number.
	  for (unsigned i = auxiliary_symbols().size(); i-- > 0; )
	    solution_or_bound
	      = solution_or_bound.substitute(auxiliary_symbols()[i],
					     get_auxiliary_definition
					     (auxiliary_symbols()[i]));
	}
	// We must copy the values in the blackboard of the linear recurrence
	// in the blackboard of the original recurrences: they could be
	// necessary in the validation's process of the non-linear recurrence.
	blackboard = rec_rewritten.blackboard;
	recurrence_rhs_rewritten = true;
	return SUCCESS;
      }
      else
	return status;

    // Functional equation.
    else
      if ((status = rec_rewritten.approximate_functional_equation())
	  == SUCCESS)
	// Transform the solution of the linear recurrence in the solution
	// of the non linear recurrence.
	if (type != 0 || (type == 0
			  && rec_rewritten.lower_bound_.expression()
			  == rec_rewritten.upper_bound_.expression())) {
	  if (type == 1)
	    solution_or_bound = pwr(base_exp_log(),
				    rec_rewritten.lower_bound_.expression());
	  else
	    solution_or_bound = pwr(base_exp_log(),
				    rec_rewritten.upper_bound_.expression());
	  solution_or_bound = substitute_x_function(solution_or_bound,
						    base_exp_log(), false);
	  solution_or_bound = simplify_ex_for_input(solution_or_bound, true);
	  solution_or_bound = simplify_logarithm(solution_or_bound);
	  // Resubstitute eventual auxiliary symbols with the respective
	  // negative number.
	  for (unsigned i = auxiliary_symbols().size(); i-- > 0; )
	    solution_or_bound
	      = solution_or_bound.substitute(auxiliary_symbols()[i],
					     get_auxiliary_definition
					     (auxiliary_symbols()[i]));
	  D_VAR(solution_or_bound);
	  recurrence_rhs_rewritten = true;
	  return SUCCESS;
	}
	else
	  return TOO_COMPLEX;
      else
	return status;
  }
  else
    return status;
}

//! \brief
//! Let \p solution_or_bound be the expression that represent the
//! solution or the bound computed for the recurrence \p *this.
//! This function substitutes eventual initial conditions specified
//! by the user shifting the solution or the bound if necessary.
/*!
  \param linear             <CODE>true</CODE> if the system has solved
                            a linear recurrence of finite order;
                            <CODE>false</CODE> if the system has solved
			    a functional equation.
  \param solution_or_bound  Contains the solution or the bound computed
                            for the recurrence \p *this in function of
			    arbitrary initial conditions.

  \return                   The solution or the bounds shifted in agreement
                            with the initial conditions inserted by the user
			    and with the arbitrary initial conditions
			    substituted with the respective values.

  We know the smallest positive integer \f$ s \f$ starting from which the
  recurrence is well-defined. This function check if in the map
  \p initial_conditions therea are some initial conditions
  of the form \f$ x(i) = k \f$ with \f$ k > s \f$: in this case
  the function shifted the solution or the bounds.
  Finally substitute to the arbitrary initial conditions in the solution or
  in the bound the eventual values specified by the user.
*/
Expr
PURRS::Recurrence::
substitute_i_c_shifting(bool linear, const Expr& solution_or_bound) const {
  assert(!initial_conditions.empty());
  Expr sol_or_bound = solution_or_bound;
  unsigned first_well_defined_rhs;
  unsigned order_or_rank;
  if (linear) {
    first_well_defined_rhs = first_well_defined_rhs_linear();
    order_or_rank = order();
  }
  else {
    first_well_defined_rhs = applicability_condition();
    order_or_rank = rank();
  }

  // If the order of linear recurrences is zero than it does not go made
  // no shifts (the rank of functional equations can not to be zero, it is
  // greater or equal to one).
  if (order_or_rank != 0) {
    // Consider the maximum index of `x' function in the map
    // `initial_conditions'.
    unsigned max_i_c = 0;
    for (std::map<unsigned, Expr>::const_iterator i
	   = initial_conditions.begin(),
	   iend = initial_conditions.end(); i != iend; ++i)
      if (i->first > max_i_c)
	max_i_c = i->first;
    
    // Shift initial conditions and the index of the recurrence `n'.
    if (first_well_defined_rhs + order_or_rank - 1 < max_i_c) {
      unsigned shift_forward = max_i_c - first_well_defined_rhs;
      sol_or_bound
	= sol_or_bound.substitute(n, n - (shift_forward - order_or_rank + 1));
      for (unsigned i = order_or_rank; i-- > 0; )
	sol_or_bound
	  = sol_or_bound.substitute(x(i + first_well_defined_rhs),
				    x(i + first_well_defined_rhs
				      + shift_forward - order_or_rank + 1));
    }
  }    

  // Substitute initial conditions with the values in the map
  // `initial_conditions'.
  for (std::map<unsigned, Expr>::const_iterator i = initial_conditions.begin(),
	 iend = initial_conditions.end(); i != iend; ++i)
    sol_or_bound = sol_or_bound.substitute(x(i->first),
					   get_initial_condition(i->first));
  return sol_or_bound;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution() const {
  D_MSG("compute_exact_solution");
  tested_exact_solution = true;
  // See if we have the exact solution already.
  if (exact_solution_.has_expression())
    return SUCCESS;

  // We may not have the exact solution explicitely, yet we may have
  // the lower and the upper bounds equal among themselves.
  if (lower_bound_.has_expression() && upper_bound_.has_expression()
      && lower_bound_.expression() == upper_bound_.expression()) {
    exact_solution_.set_expression(lower_bound_.expression());
    return SUCCESS;
  }

  Solver_Status status;
  if (is_classified
      || (status = classify_and_catch_special_cases()) == SUCCESS) {
    assert(is_linear_finite_order() || is_functional_equation()
	   || is_non_linear_finite_order() || is_linear_infinite_order());

    // Linear finite order.
    if (is_linear_finite_order()) {
      // If the greatest common divisor among the decrements is greater
      // than one, the order reduction is applicable.
      // FIXME: the order reduction is for the moment applied only to
      // recurrences with constant coefficients because the recurrences
      // with variable coefficients are not allowed with parameters.
      if (gcd_among_decrements() > 1 && is_linear_finite_order_const_coeff()) {
	if ((status = apply_order_reduction()) != SUCCESS)
	  return status;
      }
      // We do not have applied the order reduction.
      else
	if ((status = solve_linear_finite_order()) != SUCCESS)
	  return status;

      // Check if there are specified initial conditions and in this case
      // eventually shift the solution in according with them before to
      // substitute the values of the initial conditions to the
      // generic `x(i)'.
      if (!initial_conditions.empty())
	exact_solution_.set_expression
	  (substitute_i_c_shifting(true, exact_solution_.expression()));
      lower_bound_.set_expression(exact_solution_.expression());
      upper_bound_.set_expression(exact_solution_.expression());
      return SUCCESS;
    }
    // Functional equation.
    else if (is_functional_equation())
      if ((status = approximate_functional_equation()) == SUCCESS
	  && lower_bound_.expression() == upper_bound_.expression()) {
	if (!initial_conditions.empty())
	  lower_bound_.set_expression
	    (substitute_i_c_shifting(false, lower_bound_.expression()));
	upper_bound_.set_expression(lower_bound_.expression());
	exact_solution_.set_expression(lower_bound_.expression());
	return SUCCESS;
      }
      else
	return TOO_COMPLEX;
    // Non linear finite order.
    else if (is_non_linear_finite_order()){
      Expr solution;
      if ((status = compute_non_linear_recurrence(solution, 0))
	  == SUCCESS) {
	exact_solution_.set_expression(solution);
	lower_bound_.set_expression(solution);
	upper_bound_.set_expression(solution);
	return SUCCESS;
      }
      else
	return status;
    }
    // Linear infinite order.
    else {
      D_MSG("compute_linear_infinite_order_recurrence");
      std::vector<Expr> coefficients(2);
      coefficients[1] = coeff_first_order();
      Recurrence rec_rewritten(rhs_transformed_in_first_order());
      rec_rewritten.finite_order_p = new Finite_Order_Info(1, coefficients, 1);
      rec_rewritten.set_type(LINEAR_FINITE_ORDER_VAR_COEFF);
      rec_rewritten.set_inhomogeneous_term(inhomog_first_order());
      if ((status = rec_rewritten.solve_linear_finite_order())
	  == SUCCESS) {
	Expr solution = rec_rewritten.exact_solution_.expression();
	// Shift backward: n -> n - 1 and substitution of initial
	// condition `x(1) = 2*x(0)+1'
	// (there are 2 steps in 1: x(0) -> x(1) -> value_of_first_element()).
	solution = solution.substitute(n, n - 1);
	solution = solution
	  .substitute(x(rec_rewritten.first_well_defined_rhs_linear()),
		      value_of_first_element());
	//	solution = simplify_ex_for_output(solution, false);
	exact_solution_.set_expression(solution);
	lower_bound_.set_expression(solution);
	upper_bound_.set_expression(solution);
	return SUCCESS;
      }
      else
	return status;
    }
  }
  else
    return status;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_lower_bound() const {
  D_MSG("compute_lower_bound");
  // See if we have the lower bound already.
  if (lower_bound_.has_expression())
    return SUCCESS;

  // We may not have the lower bound explicitely, yet we may have
  // the exact solution.
  if (exact_solution_.has_expression()) {
    lower_bound_.set_expression(exact_solution_.expression());
    return SUCCESS;
  }

  Solver_Status status;
  if (is_classified
      || (status = classify_and_catch_special_cases()) == SUCCESS) {
    assert(is_linear_finite_order() || is_functional_equation()
	   || is_non_linear_finite_order() || is_linear_infinite_order());

    if (is_linear_finite_order() || is_linear_infinite_order())
      if (!tested_exact_solution) {
	// There is an exact solution.
	if ((status = compute_exact_solution()) == SUCCESS) {
	  lower_bound_.set_expression(exact_solution_.expression());
	  return SUCCESS;
	}
	else
	  return status;
      }
      else
	return TOO_COMPLEX;
    // Functional equation.
    else if (is_functional_equation()) {
      if ((status = approximate_functional_equation()) != SUCCESS)
	return status;
      if (!initial_conditions.empty()) {
	lower_bound_.set_expression
	  (substitute_i_c_shifting(false, lower_bound_.expression()));
	upper_bound_.set_expression
	  (substitute_i_c_shifting(false, upper_bound_.expression()));
      }
      return SUCCESS;
    }
    // Non linear finite order.
    else {
      Expr bound;
      if ((status = compute_non_linear_recurrence(bound, 1))
	  == SUCCESS) {
	lower_bound_.set_expression(bound);
	return SUCCESS;
      }
      else
	return status;
    }
  }
  else
    return status;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_upper_bound() const {
  D_MSG("compute_upper_bound");
  // See if we have the upper bound already.
  if (upper_bound_.has_expression())
    return SUCCESS;

  // We may not have the upper bound explicitely, yet we may have
  // the exact solution.
  if (exact_solution_.has_expression()) {
    upper_bound_.set_expression(exact_solution_.expression());
    return SUCCESS;
  }

  Solver_Status status;
  if (is_classified
      || (status = classify_and_catch_special_cases()) == SUCCESS) {
    assert(is_linear_finite_order() || is_functional_equation()
	   || is_non_linear_finite_order() || is_linear_infinite_order());
    
    if (is_linear_finite_order() || is_linear_infinite_order())
      if (!tested_exact_solution)
	// There is an exact solution.
	if ((status = compute_exact_solution()) == SUCCESS) {
	  upper_bound_.set_expression(exact_solution_.expression());
	  return SUCCESS;
	}
	else
	  return status;
      else
	return TOO_COMPLEX;
    // Functional equation.
    else if (is_functional_equation()) {
      if ((status = approximate_functional_equation()) != SUCCESS)
	return status;
      if (!initial_conditions.empty()) {
	upper_bound_.set_expression
	  (substitute_i_c_shifting(false, upper_bound_.expression()));
	lower_bound_.set_expression
	  (substitute_i_c_shifting(false, lower_bound_.expression()));
      }
      return SUCCESS;
    }
    // Non linear finite order.
    else {
      Expr bound;
      if ((status = compute_non_linear_recurrence(bound, 2))
	  == SUCCESS) {
	upper_bound_.set_expression(bound);
	return SUCCESS;
      }
      else
	return status;
    }
  }
  else
    return status;
}

bool
PURRS::Recurrence::OK() const {
#ifndef NDEBUG
  using std::endl;
  using std::cerr;
#endif

  switch(type_) {
  case ORDER_ZERO:
    if (finite_order_p != 0) {
#ifndef NDEBUG
      cerr << "Recurrence with type unknown or of order zero!" << endl;
#endif
      return false;
    }
  case LINEAR_FINITE_ORDER_CONST_COEFF:
  case LINEAR_FINITE_ORDER_VAR_COEFF:
    if (finite_order_p == 0) {
#ifndef NDEBUG
      cerr << "Recurrence of finite order!" << endl;
#endif
      return false;
    }
    else {
      //      if (! || !)
    }
  default:
    return true;
  }

return true;
}

void
PURRS::Recurrence::dump(std::ostream& s) const {
  s << "solved = "
    << (exact_solution_.has_expression() ? "true" : "false") << std::endl;
  s << "approximated = "
    << ((lower_bound_.has_expression() || upper_bound_.has_expression()) 
	? "true" : "false") << std::endl;
  s << "recurrence_rhs_rewritten = "
    << (recurrence_rhs_rewritten ? "true" : "false") << std::endl;
  s << "recurrence_rhs = " << recurrence_rhs << std::endl;

  s << "auxiliary_definitions:" << std::endl;
  blackboard.dump(s);
  
  if (!initial_conditions.empty()) {
    s << "Initial conditions:" << std::endl;
    for (std::map<unsigned, Expr>::const_iterator i = initial_conditions.begin(),
	   initial_conditions_end = initial_conditions.end();
	 i != initial_conditions_end; ++i)
      s << "  x(" << i->first << ")"
	<< " = " << i->second << std::endl;
  }
  
  //(*functional_eq_p).dump_homogeneous_terms(s);
  s << std::endl;
}
