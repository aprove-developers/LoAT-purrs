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
  We try to check that the solution is correct.
  - Validation of initial conditions.
    If <CODE>recurrence_rhs</CODE> is equal to \f$ x(0), \cdots, x(k) \f$ for
    \f$ n = 0, \cdots, k-1 \f$ respectively then
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
PURRS::Recurrence::verify_solution() const {
  if ((is_linear_finite_order() || is_non_linear_finite_order())
       && exact_solution_.has_expression()) {

    unsigned int order_rec;
    unsigned int first_i_c;
    if (is_non_linear_finite_order()) {
      order_rec = order_if_linear();
      first_i_c = first_i_c_if_linear();
    }
    else {
      order_rec = order();
      first_i_c = first_i_c_for_linear();
    }

    D_VAR(recurrence_rhs);
    D_VAR(order_rec);
    D_VAR(first_i_c);

    if (order_rec == 0)
      return PROVABLY_CORRECT;
    else {
      // Step 1: validation of initial conditions.
      for (unsigned i = order_rec; i-- > 0; ) {
	Expr solution_valuated = exact_solution_.expression()
	  .substitute(n, first_i_c + i);
	solution_valuated = blackboard.rewrite(solution_valuated);
	solution_valuated = simplify_all(solution_valuated);
	D_VAR(solution_valuated);
	unsigned i_c = first_i_c + i;
	if (applied_order_reduction)
	  // FIXME: !!!
	  i_c = mod(Number(i_c), Number(gcd_among_decrements())).to_unsigned();
	if (solution_valuated != x(i_c))
	  return INCONCLUSIVE_VERIFICATION;
      }
      // Step 2: find `partial_solution'.
      // The initial conditions are verified. Build the expression
      // `partial_solution' that has all terms of `solution' minus those
      // containing an initial condition.
      Expr partial_solution = 0;
      if (exact_solution_.expression().is_a_add())
	for (unsigned i = exact_solution_.expression().nops(); i-- > 0; ) {
	  if (!exact_solution_.expression().op(i).has_x_function(true))
	    partial_solution += exact_solution_.expression().op(i);
	}
      else
	if (!exact_solution_.expression().has_x_function(true))
	  partial_solution = exact_solution_.expression();
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

      substituted_rhs = recurrence_rhs;
      for (unsigned i = order_rec; i-- > 0; ) {
	terms_to_sub[i] = simplify_all(partial_solution.substitute
				       (n, n - (i + 1)));
	substituted_rhs = substituted_rhs
	  .substitute(x(n - (i + 1)), terms_to_sub[i]);
      }
      D_VEC(terms_to_sub, 0, terms_to_sub.size()-1);
      D_VAR(substituted_rhs);
      Expr diff = blackboard.rewrite(partial_solution - substituted_rhs);
      diff = simplify_all(diff);
      D_VAR(diff);
      if (!diff.is_zero())
	if (applied_order_reduction) {
	  applied_order_reduction = false;
	  // If we have applied the order reduction and we do not have success
	  // in the verification of the original recurrence, then we please
	  // ourselves if is verified the reduced recurrence.
	  Symbol r = insert_auxiliary_definition(mod(n,
						     gcd_among_decrements()));
	  unsigned dim = coefficients_lfo().size()
	    / gcd_among_decrements() + 1;
	  std::vector<Expr> new_coefficients(dim);
	  Expr inhomogeneous = 0;
	  Recurrence rec_rewritten
	    (rewrite_reduced_order_recurrence(recurrence_rhs, r,
					      gcd_among_decrements(),
					      coefficients_lfo(),
					      new_coefficients,
					      inhomogeneous));
	  rec_rewritten.finite_order_p
	    = new Finite_Order_Info(dim - 1, 0, new_coefficients, 1);
	  rec_rewritten.set_type(type());
	  rec_rewritten.set_inhomogeneous_term(inhomogeneous);
	  rec_rewritten.solve_linear_finite_order();
	  D_VAR(exact_solution_.expression());
	  return rec_rewritten.verify_solution();
	}
	else
	  return INCONCLUSIVE_VERIFICATION;
      return PROVABLY_CORRECT;
    }
  }
  // We failed to solve the recurrence.
  // If the client still insists in asking for the verification...
  return INCONCLUSIVE_VERIFICATION;
}

/*!
  Consider the right hand side \p rhs of the functional equation
  \f$ a x_{n/b} + p(n) \f$.
  If \p upper is true we try to check that the upper bound is correct;
  If \p lower is true we try to check that the lower bound is correct.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_bound(bool upper) const {
  if (is_functional_equation())
    if ((upper && upper_bound_.has_expression())
	|| (!upper && lower_bound_.has_expression())) {
      D_VAR(applicability_condition()); 
      Expr bound;
      if (upper)
	bound = upper_bound_.expression();
      else
	bound = lower_bound_.expression();
      
      // Step 1: validation of initial conditions.
      if (!validation_initial_conditions_in_bound(upper, bound,
						  applicability_condition()))
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
	  .substitute(n, applicability_condition()).is_a_number(num)
	  && num.is_negative())
	return INCONCLUSIVE_VERIFICATION;
      if (!upper && partial_bound
	  .substitute(n, applicability_condition()).is_a_number(num)
	  && num.is_positive())
	return INCONCLUSIVE_VERIFICATION;
      
      // Step 4: verification of the inductive step.
      Expr partial_bound_sub = partial_bound.substitute(n,
							n / divisors_arg()[0]);
      Expr approx = recurrence_rhs.substitute(x(n / divisors_arg()[0]),
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
	if (ok_inequalities(diff, applicability_condition()))
	  return PROVABLY_CORRECT;
    }
  return INCONCLUSIVE_VERIFICATION;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::apply_order_reduction() const {
  applied_order_reduction = true;
  // Build the new recurrence substituting `n' not contained in the
  // `x' functions with `gcd_among_decrements * n + r' and `x(n-k)' with
  // `x(n - k / gcd_among_decrements)'.
  Symbol r = insert_auxiliary_definition(mod(n, gcd_among_decrements()));
  unsigned dim = coefficients_lfo().size() / gcd_among_decrements() + 1;
  std::vector<Expr> new_coefficients(dim);
  Expr inhomogeneous = 0;
  Recurrence rec_rewritten
    (rewrite_reduced_order_recurrence(recurrence_rhs, r,
				      gcd_among_decrements(),
				      coefficients_lfo(), new_coefficients,
				      inhomogeneous));
  rec_rewritten.finite_order_p
    = new Finite_Order_Info(dim - 1, 0, new_coefficients, 1);
  rec_rewritten.set_type(type());
  rec_rewritten.set_inhomogeneous_term(inhomogeneous);
  Solver_Status status;
  if ((status = rec_rewritten.solve_linear_finite_order()) == SUCCESS) {
    // Now we must compute the solution for the original recurrence.
    // Perform three substitutions:
    // - r                      -> mod(n, gcd_among_decrements);
    // - n                      -> 1 / gcd_among_decrements
    //                             * (n - mod(n, gcd_among_decrements));
    // - x(k), k non-negative integer -> x(mod(n, gcd_among_decrements)).
    exact_solution_.set_expression(come_back_to_original_variable
				   (rec_rewritten
				    .exact_solution_.expression(), r,
				    get_auxiliary_definition(r),
				    gcd_among_decrements()));
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
    lower_bound_.set_expression(exact_solution_.expression());
    upper_bound_.set_expression(exact_solution_.expression());
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
  D_MSG(rec_rewritten.recurrence_rhs);

  Solver_Status status;
  // Classify the linear recurrence `rec_rewritten'.
  if ((status == rec_rewritten.classify_and_catch_special_cases())
      == SUCCESS) {
    assert(rec_rewritten.is_linear_finite_order()
	   || rec_rewritten.is_functional_equation());

    // Linear finite order.
    if (rec_rewritten.is_linear_finite_order())
      if ((status = rec_rewritten.solve_linear_finite_order())
	  == SUCCESS) {
	set_order_if_linear(rec_rewritten.order());
	set_first_i_c_if_linear(rec_rewritten.first_i_c_for_linear());
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
  if ((status = classify_and_catch_special_cases()) == SUCCESS) {
    assert(is_linear_finite_order() || is_functional_equation()
	   || is_non_linear_finite_order());

    // Linear finite order.
    if (is_linear_finite_order()) {
      // If the greatest common divisor among the decrements is greater
      // than one, the order reduction is applicable.
      // FIXME: the order reduction is for the moment applied only to
      // recurrences with constant coefficients because the recurrences
      // with variable coefficients are not allowed with parameters.
      if (gcd_among_decrements() > 1 && is_linear_finite_order_const_coeff())
	return apply_order_reduction();

      // We do not have applied the order reduction.
      if ((status = solve_linear_finite_order()) == SUCCESS) {
	lower_bound_.set_expression(exact_solution_.expression());
	upper_bound_.set_expression(exact_solution_.expression());
	return SUCCESS;
      }
      else
	return status;
    }
    // Functional equation.
    else if (is_functional_equation())
      if ((status = approximate_functional_equation()) == SUCCESS
	  && lower_bound_.expression() == upper_bound_.expression()) {
	exact_solution_.set_expression(lower_bound_.expression());
	return SUCCESS;
      }
      else
	return TOO_COMPLEX;
    // Non linear finite order.
    else {
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
  if ((status = classify_and_catch_special_cases()) == SUCCESS) {
    assert(is_linear_finite_order() || is_functional_equation()
	   || is_non_linear_finite_order());

    if (is_linear_finite_order())
      if (!tested_exact_solution)
	// There is an exact solution.
	if ((status = compute_exact_solution()) == SUCCESS) {
	  lower_bound_.set_expression(exact_solution_.expression());
	  return SUCCESS;
	}
	else
	  return status;
      else
	return TOO_COMPLEX;
    // Functional equation.
    else if (is_functional_equation())
      if ((status = approximate_functional_equation()) == SUCCESS)
	return SUCCESS;
      else
	return status;
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
  if ((status = classify_and_catch_special_cases()) == SUCCESS) {
    assert(is_linear_finite_order() || is_functional_equation()
	   || is_non_linear_finite_order());
    if (is_linear_finite_order())
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
    else if (is_functional_equation())
      if ((status = approximate_functional_equation()) == SUCCESS)
	return SUCCESS;
      else
	return status;
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
  case UNKNOWN:
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
  
  s << std::endl;
}
