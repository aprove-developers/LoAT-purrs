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
#include "Order_Reduction_Info.defs.hh"
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
  if (is_linear_finite_order() && exact_solution_.has_expression()) {
    D_VAR(recurrence_rhs);
    D_VAR(order());
    D_VAR(first_initial_condition());
    if (order_reduction_p) {
      D_VAR(old_recurrence_rhs());
      D_VAR(gcd_decrements_old_rhs());
    }

    if (non_linear_p)
      recurrence_rhs = original_recurrence_rhs();
    if (order() == 0)
      return PROVABLY_CORRECT;
    else {
      // Step 1: validation of initial conditions.
      for (unsigned i = order(); i-- > 0; ) {
	Expr solution_valuated = exact_solution_.expression()
	  .substitute(n, first_initial_condition() + i);
	solution_valuated = blackboard.rewrite(solution_valuated);
	solution_valuated = simplify_numer_denom(solution_valuated);
	D_VAR(solution_valuated);
	// We have to substitute `first_initial_condition() + i'
	// in `x(mod(first_initial_condition() + i, gcd))' and, when
	// `gcd_decrements_old_rhs <= i'
	// `x(mod(first_initial_condition() + i, gcd))' is equal to
	// `x(mod(first_initial_condition() + i - gcd, gcd))'.
	unsigned i_c = first_initial_condition() + i;
	if (order_reduction_p && gcd_decrements_old_rhs() <= i)
	  i_c -= gcd_decrements_old_rhs();
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
      std::vector<Expr> terms_to_sub(order());

      unsigned gcd_decrements_old;
      if (order_reduction_p && verified_one_time()) {
	// We have applied the order reduction in order to solve the
	// recurrence.
	substituted_rhs = old_recurrence_rhs();
	gcd_decrements_old = gcd_decrements_old_rhs();
      }
      else {
	// We have not applied the order reduction in order to solve the
	// recurrence.
	substituted_rhs = recurrence_rhs;
	gcd_decrements_old = 1;
      }
	
      for (unsigned i = order(); i-- > 0; ) {
	terms_to_sub[i] = simplify_all(partial_solution.substitute
				       (n, n - (i + 1)
					* gcd_decrements_old));
	substituted_rhs = substituted_rhs
	  .substitute(x(n - (i + 1) * gcd_decrements_old),
		      terms_to_sub[i]);
      }
      D_VEC(terms_to_sub, 0, terms_to_sub.size()-1);
      D_VAR(substituted_rhs);
      Expr diff = blackboard.rewrite(partial_solution - substituted_rhs);
      diff = simplify_all(diff);
      D_VAR(diff);
      if (!diff.is_zero())
	if (order_reduction_p && verified_one_time()) {
	  not_verified_one_time();
	  // If we have applied the order reduction and we do not have success
	  // in the verification of the original recurrence, then we please
	  // ourselves if is verified the reduced redurrence.
	  exact_solution_.set_expression(solution_order_reduced());
	  set_gcd_decrements_old_rhs(0);
	  return verify_solution();
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
      Expr partial_bound_sub = partial_bound.substitute(n, n / divisor_arg());
      Expr approx = recurrence_rhs.substitute(x(n / divisor_arg()),
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
PURRS::Recurrence::compute_exact_solution() const {
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
    assert(is_linear_finite_order() || is_functional_equation());
    if (is_linear_finite_order()) {
      if ((status = solve_linear_finite_order()) == SUCCESS) {
	lower_bound_.set_expression(exact_solution_.expression());
	upper_bound_.set_expression(exact_solution_.expression());
	return SUCCESS;
      }
      else
	return status;
    }
    else // case functional equation
      if ((status = approximate_functional_equation()) == SUCCESS
	  && lower_bound_.expression() == upper_bound_.expression()) {
	exact_solution_.set_expression(lower_bound_.expression());
	return SUCCESS;
      }
      else
	return TOO_COMPLEX;
  }
  else
    return status;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_lower_bound() const {
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
    assert(is_linear_finite_order() || is_functional_equation());
    if (!tested_exact_solution && compute_exact_solution() == SUCCESS) {
      lower_bound_.set_expression(exact_solution_.expression());
      return SUCCESS;
    }
    if (is_functional_equation()
	&& (status = approximate_functional_equation()) == SUCCESS)
      return SUCCESS;
    else
      return TOO_COMPLEX;
  }
  else
    return status;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_upper_bound() const {
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
    assert(is_linear_finite_order() || is_functional_equation());
    if (!tested_exact_solution && compute_exact_solution() == SUCCESS) {
      upper_bound_.set_expression(exact_solution_.expression());
      return SUCCESS;
    }
    if (is_functional_equation()
	&& (status = approximate_functional_equation()) == SUCCESS)
      return SUCCESS;
    else
      return TOO_COMPLEX;
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

  switch(type) {
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
