/* Recurrence class implementation (non-inline functions).
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

#include "Recurrence.defs.hh"
#include "Recurrence.inlines.hh"
#include "ep_decomp.hh"
#include "finite_order.hh"
#include "simplify.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include "Cached_Expr.defs.hh"
#include "Non_Linear_Info.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Blackboard.defs.hh"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

const PURRS::Symbol&
PURRS::Recurrence::n = Symbol("n");

void
PURRS::Recurrence::throw_invalid_argument(const char* method,
					  const char* reason) const {
  std::ostringstream s;
  s << method << ":" << std::endl
    << reason;
  throw std::invalid_argument(s.str());
}

PURRS::index_type
PURRS::Recurrence::first_valid_initial_condition() const {
  // FIXME: Temporary. Add support for exceptions. Actually compute
  // the first valid index for initial conditions instead of just
  // returning `first_valid_index'.
  return first_valid_index;
}

void
PURRS::Recurrence::
set_initial_conditions(const std::map<index_type, Expr>& initial_conditions) {
  if (!initial_conditions_.empty())
    reset_initial_conditions();

  if (classifier_status_ == NOT_CLASSIFIED_YET)
    classify_and_catch_special_cases();

  if (classifier_status_ != CL_SUCCESS)
    throw std::runtime_error("PURRS::Recurrence::set_initial_conditions():\n"
			     "sorry, this recurrence is too difficult.");
  const char* method = "PURRS::Recurrence::set_initial_conditions()";
  std::ostringstream s;
  switch (type_) {
  case ORDER_ZERO:
    // If the recurrence has order zero, the solution is simply the
    // right-hand side.
    break;
  case LINEAR_FINITE_ORDER_CONST_COEFF:
  case LINEAR_FINITE_ORDER_VAR_COEFF:
  case NON_LINEAR_FINITE_ORDER:
    {
      index_type k;
      if (is_non_linear_finite_order()) {
	if (associated_linear_rec().recurrence_rhs == 0)
	  k = 1;
	else {
	  associated_linear_rec().classify_and_catch_special_cases();
	  k = associated_linear_rec().order();
	}
      }
      else
	k = order();
      index_type first = first_valid_index;
      // Check the number of initial conditions.
      if (initial_conditions.size() < k) {
	s << "this is a recurrence of order " << k
	  << " and at least " << k << " initial conditions\n"
	  << "are necessary to uniquely determine it.";
	throw_invalid_argument(method, s.str().c_str());
      }
      // Check if the largest index of the initial conditions in the
      // map `initial_conditions' is not smaller than `first_valid_index'.
      assert(!initial_conditions.empty());
      index_type largest = initial_conditions.rbegin()->first;
      if (largest < first) {
	s << "this is a recurrence of order " << k
	  << " and it is well-defined only for `n >= " << first << "'.\n"
	  << "Hence, the largest index can not to be smaller than "
	  << first << ".";
	throw_invalid_argument(method, s.str().c_str());
      }
      // Check if the indices of the last `k' initial conditions are
      // all consecutive (the `k'-th condition is already ok for the
      // previous check).
      unsigned int j = 1;
      for (std::map<index_type, Expr>::const_reverse_iterator i 
	     = initial_conditions.rbegin(); j < k; j++) {
	// Skip the check on the largest initial condition.
	i++;
	if (i->first != largest - j) {
	  s << "this is a recurrence of order " << k
	    << " and the " << k << " largest initial conditions\n"
	    << " must be consecutive.";
	  throw_invalid_argument(method, s.str().c_str());
	}
      }
      // In the case of non-linear recurrences we at the moment accept
      // only non-negative values for the initial conditions.
      if (is_non_linear_finite_order()) {
	for (std::map<index_type, Expr>::const_iterator i
	       = initial_conditions.begin(),
	       initial_conditions_end = initial_conditions.end();
	     i != initial_conditions_end; ++i) {
	  Number index;
	  if (i->second.is_a_number(index) && index.is_negative())
	    throw_invalid_argument(method,
				   "this is a non-linear recurrence "
				   "and at the moment we accept only\n"
				   "non-negative values for the initial "
				   "conditions.");
	}
      }
    }
    break;
  case WEIGHTED_AVERAGE:
    {
      // The user is in this case not with the driver or with the demo
      // but if call this method passing an empty map of initial
      // conditions.
      if (initial_conditions.empty())
	throw_invalid_argument(method,
			       "this is a weighted-average recurrence "
			       "and we need one initial condition\n"
			       "to uniquely determine it.");
      if (initial_conditions.rbegin()->first < first_valid_index) {
	s << "this recurrence is well-defined for `n >= " << first_valid_index
	  << "'.\nAt least an initial condition with index "
	  << "larger than or equal to " << first_valid_index
	  << "\nmust be present.";
	throw_invalid_argument(method, s.str().c_str());
      }
    }
    break;
  case FUNCTIONAL_EQUATION:
    {
      unsigned int number_i_c
	= (functional_eq_p->ht_begin()->first - 1).to_unsigned_int();;
      // Check the number of initial conditions.
      if (initial_conditions.size() != number_i_c) {
	if (number_i_c == 1)
	  s << "this functional equation needs only "
	    << "the initial condition x(1).";
	else
	  s << "for this functional equation " << number_i_c
	    << " initial conditions are necessary\n"
	    << "to uniquely determine it.";
	throw_invalid_argument(method, s.str().c_str());
      }
      // Check if the indexes of the last `number_i_c' initial conditions
      // are all consecutive in the interval `[1, number_ic]' and the values
      // are numeric.
      if (number_i_c == 1) {
	if (initial_conditions.begin()->first != 1) {
	  s << "this functional equation needs only "
	    << "the initial condition x(1).";
	  throw_invalid_argument(method, s.str().c_str());
	}
	if (!initial_conditions.begin()->second.is_a_number()) {
	  s << "at the moment are allowed only numerical values\n"
	    << "for the initial conditions.";
	  throw_invalid_argument(method, s.str().c_str());
	}
	// FIXME: very temporary `else'!!
	// If the number is negative the homogeneous parts of lower bound
	// and upper bound must be swapped!
	else {
	  Number num = initial_conditions.begin()->second.ex_to_number();
	  if (num.is_negative()) {
	    s << "at the moment are allowed only non-negative values\n"
	      << "for the initial conditions of functional equations.";
	    throw_invalid_argument(method, s.str().c_str());
	  }
	}
      }
      else {
	unsigned int j = 1;
	for (std::map<index_type, Expr>::const_iterator i 
	       = initial_conditions.begin(); j <= number_i_c; ++j) {
	  if (i->first != j) {
	    s << "the indices of the " << number_i_c << " initial conditions\n"
	      << "must be in the interval [1, " << number_i_c << "].";
	    throw_invalid_argument(method, s.str().c_str());
	  }
	  if (!i->second.is_a_number()) {
	    s << "at the moment are allowed only numerical values\n"
	      << "for the initial conditions.";
	    throw_invalid_argument(method, s.str().c_str());
	  }
	  // FIXME: very temporary `else'!!
	  // If the number is negative the homogeneous parts of lower bound
	  // and upper bound must be swapped!
	  else {
	    Number num = i->second.ex_to_number();
	    if (num.is_negative()) {
	      s << "at the moment are allowed only non-negative values\n"
		<< "for the initial conditions of functional equations.";
	      throw_invalid_argument(method, s.str().c_str());
	    }
	  }
	  ++i;
	}
      }
    }
    break;
  default:
    throw std::runtime_error("PURRS internal error: "
			     "set_initial_conditions().");
  }

  // Set the private data `initial_conditions' with the map given.
  initial_conditions_ = initial_conditions;
}

void
PURRS::Recurrence::set_first_valid_index_for_solution() const {
  assert(exact_solution_.has_expression()
	 || lower_bound_.has_expression() || upper_bound_.has_expression());
  index_type index = first_valid_index;
  if (initial_conditions_.empty()) {
    if (is_weighted_average())
      ++index;
  }
  else {
    if (is_linear_finite_order_const_coeff()
	|| is_linear_finite_order_var_coeff()
	|| is_non_linear_finite_order()) {
      index_type order_rec;
      if (is_non_linear_finite_order())
	// This is the case of non-linear recurrences of the form
	// `x(n) = c x(n-1)^a', where `c' and `a' are constants (`a != 1').
	if (coeff_simple_non_linear_rec() != 0)
	  order_rec = 1;
	else
	  order_rec = associated_linear_rec().order();
      else
	order_rec = order();
      std::map<index_type, Expr>::const_reverse_iterator i
	= initial_conditions_.rbegin();
      for (index_type j = 0; j < order_rec-1; j++)
	i++;
      index = std::max(index, i->first);
    }
    else if (is_weighted_average())
      index = std::max(index, initial_conditions_.rbegin()->first) + 1;
  }
  first_valid_index_for_solution_ = index;
}

void
PURRS::Recurrence::check_number_for_evaluation(const char* method,
					       const char* name,
					       const Number& x) const {
  // Check that `x' is non-negative.
  std::ostringstream s;
  if (x.is_negative()) {
    s << "the " << name << " can be evaluated on non-negative numbers";
    throw_invalid_argument(method, s.str().c_str());
  }
  
  switch (type_) {
  case LINEAR_FINITE_ORDER_CONST_COEFF:
  case LINEAR_FINITE_ORDER_VAR_COEFF:
    {
      // Check the number `x' is bigger than
      // `first_valid_index - order() + 1' or
      // is equal to the index of a initial condition.
      unsigned int min_index = first_valid_index + 1 >= order()
	? first_valid_index - order() + 1 : 0;
      if (x < min_index) {
	std::map<index_type, Expr>::const_iterator i
	  = initial_conditions_.find(x.to_unsigned_int());
	if (i == initial_conditions_.end()) {
	  s << "*this is a recurrence of order " << order() << ";\n"
	    << "the least non-negative integer `j' such that the\n"
	    << "recurrence is well-defined for `n >= j' is "
	    << first_valid_index
	    << ".\nThe solution (or the bound) can be evaluated\n"
	    << "for n >= " << min_index << " or must exist an initial\n"
	    << "condition with index equal to " << x;
	  throw_invalid_argument(method, s.str().c_str());
	}
      }
    }
    break;
  case ORDER_ZERO:
  case WEIGHTED_AVERAGE:
    // Check the number `x' is bigger than
    // `first_valid_index' or is equal to the
    // index of a initial condition.
    if (x <= first_valid_index) {
      std::map<index_type, Expr>::const_iterator i
	= initial_conditions_.find(x.to_unsigned_int());
      if (i == initial_conditions_.end()) {
	s << "*this is a weighted-average recurrence and\n"
	  << "the solution (or the bound) is valid from n = "
	  << first_valid_index + 1 << ", \n"
	  << "i.e. from the least non-negative integer `j'\n"
	  << "such that the recurrence is well-defined for `n >= j' plus 1.\n"
	  << "The solution (or the bound) can be evaluated\n"
	  << "for n > " << first_valid_index << " or must exist an initial\n"
	  << "condition with index equal to " << x;
	throw_invalid_argument(method, s.str().c_str());
      }
    }
    break;
  case NON_LINEAR_FINITE_ORDER:
  case FUNCTIONAL_EQUATION:
    // FIXME: to do!
    throw
      "PURRS error: today the evaluation of the solution\n"  
      "is allowed only for linear finite order recurrences and\n"  
      "weighted-average recurrences.";
    break;
  default:
    throw std::runtime_error("PURRS internal error: "
			     "check_number_for_evaluation().");
  }
}

PURRS::Expr
PURRS::Recurrence::evaluate(unsigned int kind, const Number& x) const {
  const char* method;
  const char* name;
  Expr evaluated;
  switch (kind) {
  case 0:
    method = "PURRS::Recurrence::evaluate_exact_solution()";
    name = "solution";
    exact_solution(evaluated);
    break;
  case 1:
    method = "PURRS::Recurrence::evaluate_lower_bound()";
    name = "lower bound";
    lower_bound(evaluated);
    break;
  case 2:
    method = "PURRS::Recurrence::evaluate_upper_bound()";
    name = "upper bound";
    upper_bound(evaluated);
    break;
  default:
    throw std::runtime_error("PURRS internal error: "
			     "evaluate().");    
  }
  // Check that the number `x' is in agreement with the type of the
  // recurrence.
  check_number_for_evaluation(method, name, x);

  // The solution or the bound is valid for the non-negative integer
  // `n' such that `n > first_valid_index'; for `n <= first_valid_index'
  //  the solution is represented by the initial condition with index
  // equal to `x'.
  unsigned int max_index = get_max_index_initial_condition();
  switch (type_) {
  case ORDER_ZERO:
    // If the recurrence has order zero, the solution is simply the
    // right-hand side.
    evaluated = evaluated.substitute(n, x);
    break;
  case LINEAR_FINITE_ORDER_CONST_COEFF:
  case LINEAR_FINITE_ORDER_VAR_COEFF:
  case WEIGHTED_AVERAGE:
    if (x <= max_index)
      evaluated = get_initial_condition(x.to_unsigned_int());
    else
      evaluated = evaluated.substitute(n, x);
    break;
  case NON_LINEAR_FINITE_ORDER:
  case FUNCTIONAL_EQUATION:
    // FIXME: to do!
    throw
      "PURRS error: today the evaluation of the solution\n"  
      "is allowed only for linear finite order recurrences and\n"  
      "weighted-average recurrences.";
    break;
  default:
    throw std::runtime_error("PURRS internal error: "
			     "evaluate_exact_solution().");
  }
  return evaluated;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::apply_order_reduction() const {
  set_order_reduction();
  // Build the new recurrence substituting `n' not contained in the
  // `x' functions with `gcd_among_decrements * n + r' and `x(n-k)' with
  // `x(n - k / gcd_among_decrements)'.
  Symbol r = insert_auxiliary_definition(mod(n, gcd_among_decrements()));
  unsigned int dim = coefficients().size() / gcd_among_decrements() + 1;
  std::vector<Expr> new_coefficients(dim);
  Expr inhomogeneous = 0;
  Recurrence rec_rewritten
    (write_reduced_order_recurrence(recurrence_rhs, r, gcd_among_decrements(),
				    coefficients(), new_coefficients,
				    inhomogeneous));
  rec_rewritten.finite_order_p
    = new Finite_Order_Info(dim - 1, new_coefficients, 1);
  rec_rewritten.set_type(type());
  rec_rewritten.set_inhomogeneous_term(inhomogeneous);
  rec_rewritten.set_first_valid_index(first_valid_index);
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
    exact_solution_.set_expression
      (simplify_ex_for_output(exact_solution_.expression(), false));

    // If there are defined initial conditions in the map
    // `initial_conditions' then we must to deduce the exact form
    // of the solution and substitute to it the defined initial
    // conditions.
    if (!initial_conditions_.empty()) {
      Expr expanded_solution
	= write_expanded_solution(*this, gcd_among_decrements());
      exact_solution_.set_expression(expanded_solution);
    }
    
    bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
    rec_rewritten = true;
    return SUCCESS;
  }
  else
    return status;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::map_status(Classifier_Status classifier_status) {
  switch (classifier_status) {
  case CL_SUCCESS:
    return SUCCESS;
  case CL_INDETERMINATE_RECURRENCE:
    return INDETERMINATE_RECURRENCE;
  case CL_UNSOLVABLE_RECURRENCE:
    return UNSOLVABLE_RECURRENCE;
  case CL_HAS_NON_INTEGER_DECREMENT:
    // Intentionally fall through.
  case CL_MALFORMED_RECURRENCE:
    return MALFORMED_RECURRENCE;
  case CL_DOMAIN_ERROR:
    return DOMAIN_ERROR;
  case CL_HAS_HUGE_DECREMENT:
    // Intentionally fall through.
  case CL_TOO_COMPLEX:
    return TOO_COMPLEX;
  default:
    throw std::runtime_error("PURRS internal error: map_status().");
    break;
  }
}

/*!
  Compute the solution of non-linear recurrence.
  If the recurrence is of the form \f$ x(n) = c x(n-1)^{\alpha} \f$
  then we already know the solution which will store in \p solution_or_bound:
  \f$ x(n) = x(0)^{\alpha^n} c^{\frac{\alpha^n-1}{\alpha-1}} \f$.
  In all the other cases this function builds a new object recurrence
  containing a linear recurrence obtained from that one non-linear, classify
  and solve it; at the end transform the solution of the linear recurrence
  in the solution of the non-linear recurrence.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::
compute_non_linear_recurrence(Expr& solution_or_bound,
			      unsigned int type) const {
  // We consider the simple case of non-linear recurrence of the form
  // `x(n) = c x(n-1)^a', where `c' and `a' are constants (`a != 1').
  // In this case we already know the solution:
  // `x(n) = c^((1-a^n)/(1-a)) x(0)^(a^n)'.
  // FIXME: the system will must consider the following conditions:
  // a > 0, a not natural number -> x(0) >= 0, c >= 0;
  // a < 0, a not natural number -> x(0) > 0, c > 0.
  if (coeff_simple_non_linear_rec() != 0) {
    assert(base_exp_log().is_a_number());
    // Even if in this case we already know the solution,
    // in order to apply the validation's process of the solution,
    // is necessary to know the order of the recurrence.
    associated_linear_rec().classify_and_catch_special_cases();
    solution_or_bound = pwr(coeff_simple_non_linear_rec(),
			    (pwr(base_exp_log(), n)-1)/(base_exp_log()-1))
      * pwr(x(0), pwr(base_exp_log(), n));
    return SUCCESS;
  }
  Classifier_Status classifier_status
    = associated_linear_rec().classify_and_catch_special_cases();
  // Classify the linear recurrence `associated_linear_rec()'.
  if (classifier_status == CL_SUCCESS) {
    assert(associated_linear_rec().is_linear_finite_order()
	   || associated_linear_rec().is_functional_equation());
    
    Solver_Status status;
    // Linear finite order.
    if (associated_linear_rec().is_linear_finite_order())
      // To call `solve_linear_finite_order()' is not enough because
      // it could be possible to apply the order reduction method
      // on `associated_linear_rec()'. So it is necessary to use
      // `compute_exact_solution()'.
      if ((status
	   = associated_linear_rec().compute_exact_solution_finite_order())
	  == SUCCESS) {
	// Transform the solution of the linear recurrence in the solution
	// of the non linear recurrence.
	if (associated_linear_rec().exact_solution_.expression() == 0)
	  solution_or_bound = 0;
	else {
	  solution_or_bound
	    = pwr(base_exp_log(),
		  associated_linear_rec().exact_solution_.expression());
	  solution_or_bound = substitute_x_function(solution_or_bound,
						    base_exp_log(),
						    ARGUMENT_LOG);
	  solution_or_bound = simplify_ex_for_input(solution_or_bound, true);
	  solution_or_bound = simplify_logarithm(solution_or_bound);
	  // Resubstitute eventual auxiliary symbols with the respective
	  // negative number.
	  for (unsigned int i = auxiliary_symbols().size(); i-- > 0; )
	    solution_or_bound
	      = solution_or_bound.substitute(auxiliary_symbols()[i],
					     get_auxiliary_definition
					     (auxiliary_symbols()[i]));
	}
	// We must copy the values in the blackboard of the linear recurrence
	// in the blackboard of the original recurrences: they could be
	// necessary in the validation's process of the non-linear
	// recurrence.
	blackboard = associated_linear_rec().blackboard;
	bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	rec_rewritten = true;
	return SUCCESS;
      }
      else
	return status;

    // Functional equation.
    else {
      if (type == 1) {
	status
	  = associated_linear_rec().approximate_functional_equation(LOWER);
	if (status == SUCCESS)
	  solution_or_bound
	    = pwr(base_exp_log(),
		  associated_linear_rec().lower_bound_.expression());
	else
	  return status;
      }
      else if (type == 2) {
	status
	  = associated_linear_rec().approximate_functional_equation(UPPER);
	if (status == SUCCESS)
	  solution_or_bound
	    = pwr(base_exp_log(),
		  associated_linear_rec().upper_bound_.expression());
	else
	  return status;
      }
      else {// type == 0
	if ((status
	     = associated_linear_rec().approximate_functional_equation(LOWER))
	    != SUCCESS
	    || (status
		= associated_linear_rec().approximate_functional_equation(UPPER))
	    != SUCCESS)
	  return status;
	if (associated_linear_rec().lower_bound_.expression()
	    != associated_linear_rec().upper_bound_.expression())
	  return TOO_COMPLEX;
      }
      solution_or_bound = substitute_x_function(solution_or_bound,
						base_exp_log(), ARGUMENT_LOG);
      solution_or_bound = simplify_ex_for_input(solution_or_bound, true);
      solution_or_bound = simplify_logarithm(solution_or_bound);
      // Resubstitute eventual auxiliary symbols with the respective
      // negative number.
      for (unsigned int i = auxiliary_symbols().size(); i-- > 0; )
	solution_or_bound
	  = solution_or_bound.substitute(auxiliary_symbols()[i],
					 get_auxiliary_definition
					 (auxiliary_symbols()[i]));
      bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
      rec_rewritten = true;
      return SUCCESS;
    }
  }
  else
    return map_status(classifier_status);
}

namespace {
using namespace PURRS;

Expr
increase_argument_x_function(const Expr& e, unsigned int num) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += increase_argument_x_function(e.op(i), num);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten *= increase_argument_x_function(e.op(i), num);
  }
  else if (e.is_a_power())
    return pwr(increase_argument_x_function(e.arg(0), num),
	       increase_argument_x_function(e.arg(1), num));
  else if (e.is_a_function()) {
    if (e.is_the_x_function())
      return x(e.arg(0)+num);
    else if (e.nops() == 1)
      return PURRS::apply(e.functor(), increase_argument_x_function(e.arg(0), num));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int j = 0; j < num_argument; ++j)
	argument[j] = increase_argument_x_function(e.arg(j), num);
      return PURRS::apply(e.functor(), argument);
    }
  }
  else
    e_rewritten = e;
  return e_rewritten; 
}

} // anonymous namespace

//! \brief
//! Solve the weighted-average recurrence in \ref normal_form "normal form"
//! \f[
//!   x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n)
//! \f]
//! using the associated recurrence of the first order whose solution
//! is valid for \f$ n > 1 \f$
//! \f[
//!   x(n) = \frac{f(n)}{1-f(n)} \sum_{k=n_0}^{n-1} x(k) + \frac{g(n)}{1-f(n)};
//! \f]
//! For \f$ n = 1 \f$ it must consider \f$ x(1) = f(1) x(0) + g(1) \f$.
/*!
  In the classification's process have been considered the following steps:
  - eventual rewriting in the normal form of the recurrence;
  - computation of the right hand side of the associated first order
    recurrence;
  - shift forward of the first order recurrence: \f$ n = n + 1 \f$.
  This function performs these other steps:
  - computation of the solution of the first order recurrence;
  - shift backward of the solution: \f$ n = n - 1 \f$;
  - substitution of the initail condition \f$ x(1) = f(1) x(0) + g(1) \f$.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::
compute_weighted_average_recurrence(Expr& solution) const {
  if (weight() == -1) {
    // Special case: `f(n) = -1'.
    // In this case the solution of the recurrence is simply
    // `x(n) = g(n) - g(n-1)' and is not necessary all the normal procedure.
    solution
      = inhomogeneous_term - inhomogeneous_term.substitute(n, n-1).expand();
    return SUCCESS;
  }
  else {
    // Classify the first order recurrence relation associated to that
    // one of type weighted-average.
    associated_first_order_rec().classify_and_catch_special_cases();
    Solver_Status status;
    if ((status = associated_first_order_rec().solve_linear_finite_order())
	== SUCCESS) {
      solution = associated_first_order_rec().exact_solution_.expression();

      unsigned int lower;
      Expr weight_rec;
      Expr inhomogeneous;
      if (recurrence_rewritten) {
	weight_rec = original_weight();
	inhomogeneous = original_inhomogeneous();
	lower = lower_limit();
      }
      else {
	weight_rec = weight();
	inhomogeneous = inhomogeneous_term;
	lower = 0;
      }
      
      // Shift backward: n -> n - 1 - lower.
      solution = solution.substitute(n, n - 1 - lower);
      solution = solution.substitute(x(0),
				     (weight_rec*x(lower)+inhomogeneous)
				     .substitute(n, lower+1));
      //        solution = simplify_ex_for_output(solution, false);
      return SUCCESS;
    }
    else
      return status;
  }
}

namespace {
using namespace PURRS;

void
get_max_index_symbolic_i_c(const Expr& e, unsigned int& max_index) {
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      get_max_index_symbolic_i_c(e.op(i), max_index);
  }
  else if (e.is_a_power()) {
    get_max_index_symbolic_i_c(e.arg(0), max_index);
    get_max_index_symbolic_i_c(e.arg(1), max_index);
  }
  else if (e.is_a_function())
    if (e.is_the_x_function()) {
      const Expr& argument = e.arg(0);
      Number num;
      if (argument.is_a_number(num) && num.is_integer() && num > max_index)
	max_index = num.to_unsigned_int();
    }
    else
      for (unsigned int i = e.nops(); i-- > 0; )
	get_max_index_symbolic_i_c(e.arg(i), max_index);
}

} // anonymous namespace

#if 0
// FIXME: this is an old version that substitutes initial conditions
// in the symbolic solution of the recurrence.
Expr
PURRS::Recurrence::
subs_i_c_finite_order_and_functional_eq(const Expr& solution_or_bound) const {
  Expr sol_or_bound = solution_or_bound;

  index_type first_valid_index_rhs;
  index_type order_or_rank;
  if (is_linear_finite_order()) {
    first_valid_index_rhs = first_valid_index;
    order_or_rank = order();
  }
  else if (is_functional_equation()) {
    first_valid_index_rhs = applicability_condition();
    order_or_rank = rank();
  }
  else {
    assert(is_non_linear_finite_order());
    // FIXME: we must again understand how to work in these cases.
    return sol_or_bound;
  }

  // If the order of linear recurrences is zero than it does not perform
  // shifts (the rank of functional equations can not to be zero, it is
  // greater or equal to one).
  if (order_or_rank != 0) {
    // `max_index_symb_i_c' represents the largest
    // index among the symbolic initial conditions
    // in the solution or in the bound.
    unsigned int max_index_symb_i_c = 0;
      //      = first_valid_index_rhs + order_or_rank - 1;
    get_max_index_symbolic_i_c(solution_or_bound, max_index_symb_i_c);

    // Consider the maximum index of `x' function in the map
    // `initial_conditions', i.e. the largest index of the initial
    // conditions inserted by the user.
    unsigned int max_index_user_i_c = get_max_index_initial_condition();

    // If `max_index_user_i_c' is bigger than `max_index_symb_i_c',
    // then we must shift the solution or the bound.
    // There are two different steps:
    // - 1. shift the index of the symbolic initial conditions:
    //      add `max_index_user_i_c' to the index of the symbolic initial
    //      conditions;
    // - 2. shift the index of the recurrence `n':
    //      `n = n - max_index_user_i_c'.
    if (max_index_user_i_c > max_index_symb_i_c) {
      // Step 1.
      for (index_type i = 0; i < order_or_rank; ++i) {
	sol_or_bound = sol_or_bound.substitute(x(i + first_valid_index_rhs),
					       x(i + max_index_user_i_c
						 - order_or_rank + 1));
      }
      // Step 2.
      if (is_linear_finite_order_var_coeff()) {
	// The solution of `x(n) = a(n) x(n-1) + p(n)' is of the form
	// `x(n) = prod(k, i+1, n, a(k)) x(i)
	//         + prod(k, i+1, n, a(k)) sum(k, i, n, p(k)/a!(k))', where
	// `i' is the least non-negative integer such that the recurrence
	// is well-defined for `n >= i'.
	// If `max_index_user_i_c', `m' for short, is bigger
	// than `i', then the solution is:
	// `x(n) = prod(k, i+1, n, a(k)) / prod(k, i+1, m, a(k)) x(i)
	//         + prod(k, i+1, n, a(k))
	//         * [sum(k, i+1, n, p(k) / prod(j, i+1, k, a(j)))
	//            - sum(k, i+1, m, p(k) / prod(j, i+1, k, a(j)))]'.
	const Expr& homogeneous_term
	  = product_factor() * x(max_index_user_i_c);
	const Expr& non_homogeneous_term = sol_or_bound - homogeneous_term;
	Symbol index;
	sol_or_bound = homogeneous_term
	  / PURRS::prod(index, first_valid_index_rhs+1, max_index_user_i_c,
			coefficients()[1].substitute(n, index)).ex_to_number()
	  + non_homogeneous_term - product_factor()
	  * PURRS::sum(index, first_valid_index_rhs + 1, max_index_user_i_c,
		       (inhomogeneous_term / product_factor())
		       .substitute(n, index));
      }
      else
	// FIXME: this technique is surely valid in the case of
	// linear recurrences with constant coefficients of finite order.
	// To check if it is valid also in the case of functional equations
	// and weighted-average recurrences.
	sol_or_bound
	  = sol_or_bound.substitute(n, n - (max_index_user_i_c
					    - first_valid_index_rhs
					    - order_or_rank + 1));
    }
  }
  // Substitute symbolic initial conditions with the values in the map
  // `initial_conditions'.
  for (std::map<index_type, Expr>::const_iterator i
	 = initial_conditions_.begin(),
	 iend = initial_conditions_.end(); i != iend; ++i)
    sol_or_bound = sol_or_bound.substitute(x(i->first),
					   get_initial_condition(i->first));
  sol_or_bound = simplify_numer_denom(sol_or_bound);
  sol_or_bound = simplify_ex_for_output(sol_or_bound, false);

  return sol_or_bound;  
}
#endif

Expr
PURRS::Recurrence::
compute_solution_linear_finite_order_on_i_c(const Expr& solution) const {
  Expr solution_on_i_c = solution;
  // We will store here the order of the recurrence.
  index_type order_rec = order();
  
  if (initial_conditions_.rbegin()->first >= first_valid_index + order_rec) {
    // Substitute possibly symbolic initial conditions occurring in
    // `solution' (at most `k', where `k' is the order of the recurrence)
    // with arbitrary symbols, which will be also the unknowns of the linear
    // system.
    Expr_List unknowns;
    for (unsigned int i = first_valid_index;
	 i < first_valid_index + order_rec; ++i) {
      Symbol y = Symbol();
      unknowns.append(y);
      solution_on_i_c = solution_on_i_c.substitute(x(i), y);
    }

    // Build the equations of the linear system to solve:
    // substitute to `solution_on_i_c' the `k' biggest index
    // among the indices of the initial conditions contained
    // in `initial_conditions_'.
    assert(!initial_conditions_.empty());
    Expr_List equations;
    unsigned int j = 0;
    for (std::map<index_type, Expr>::const_reverse_iterator i 
	   = initial_conditions_.rbegin(); j < order_rec; j++) {
      Expr lhs = solution_on_i_c.substitute(n, i->first);
      equations.prepend(Expr(lhs, i->second));
      i++;
    }
    
    // Solve the linear system and put the results in the
    // expression `sol_system' (which is a list of equations).
    Expr sol_system = lsolve(equations, unknowns);

    for (unsigned i = sol_system.nops(); i-- > 0; )
      solution_on_i_c = solution_on_i_c.substitute(unknowns.op(i),
						   sol_system.op(i).op(1));
  }
  else {
    // Substitute symbolic initial conditions with the values in the map
    // `initial_conditions'.
    for (std::map<index_type, Expr>::const_iterator i
	   = initial_conditions_.begin(),
	   iend = initial_conditions_.end(); i != iend; ++i)
      solution_on_i_c
	= solution_on_i_c.substitute(x(i->first),
				     get_initial_condition(i->first));
  }
  solution_on_i_c = simplify_numer_denom(solution_on_i_c);
  solution_on_i_c = simplify_ex_for_output(solution_on_i_c, false);
  return solution_on_i_c;
}

Expr
PURRS::Recurrence::
compute_solution_non_linear_finite_order_on_i_c(const Expr& solution) const {
  Expr solution_on_i_c = solution;
  // We will store here the order of the recurrence.
  index_type order_rec = associated_linear_rec().order();

  // If one of the last `k' (`k' is the order of the recurrence)
  // initial conditions is 0 then the solution of the recurrence
  // is 0.
  unsigned int j = 0;
  for (std::map<index_type, Expr>::const_reverse_iterator i 
	 = initial_conditions_.rbegin(); j < order_rec; j++)
    if (i->second == 0)
      return 0;

  if (initial_conditions_.rbegin()->first >= first_valid_index + order_rec) {
    // In the simple case of non-linear recurrence of the form
    // `x(n) = c x(n-1)^a', where `c' and `a' are constants (`a != 1').
    // the solution is `x(n) = c^((1-a^n)/(1-a)) x(0)^(a^n)'.
    // The solution for the recurrence withe the initial condition
    // `x(k)=h', `k > first_valid_index + order_rec' is
    // `x(n) = c^((1-a^(n-k))/(1-a)) h^(a^(n-k))'.
    if (coeff_simple_non_linear_rec() != 0) {
      solution_on_i_c
	= solution_on_i_c.substitute(x(0),
				     initial_conditions_.rbegin()->second);
      solution_on_i_c
	= solution_on_i_c.substitute(n,
				     n - initial_conditions_.rbegin()->first);
      D_VAR(solution_on_i_c);
      return solution_on_i_c;
    }

    Expr solution_of_linear_rec
      = associated_linear_rec().exact_solution_.expression();
    // Substitute possibly symbolic initial conditions occurring in
    // `solution_of_linear_rec' (at most `k', where `k' is the order
    // of the recurrence) with arbitrary symbols, which will be also
    // the unknowns of the linear system.
    Expr_List unknowns;
    for (unsigned int i = first_valid_index;
	 i < first_valid_index + order_rec; ++i) {
      Symbol y = Symbol();
      unknowns.append(y);
      solution_of_linear_rec = solution_of_linear_rec.substitute(x(i), y);
    }

    // Build the equations of the linear system to solve:
    // substitute to `solution_of_linear_rec' the `k' biggest index
    // among the indices of the initial conditions contained
    // in `initial_conditions_'.
    assert(!initial_conditions_.empty());
    Expr_List equations;
    unsigned int j = 0;
    for (std::map<index_type, Expr>::const_reverse_iterator i 
	   = initial_conditions_.rbegin(); j < order_rec; j++) {
      Expr lhs = solution_of_linear_rec.substitute(n, i->first);
      equations.prepend(Expr(lhs, log(i->second)/log(base_exp_log())));
      i++;
    }
    
    // Solve the linear system and put the results in the
    // expression `sol_system' (which is a list of equations).
    Expr sol_system = lsolve(equations, unknowns);

    for (unsigned i = sol_system.nops(); i-- > 0; )
      solution_of_linear_rec
	= solution_of_linear_rec.substitute(unknowns.op(i),
					    sol_system.op(i).op(1));

    // Transform the solution of the linear recurrence in the solution
    // of the non linear recurrence.
    solution_on_i_c = pwr(base_exp_log(), solution_of_linear_rec);
    solution_on_i_c = substitute_x_function(solution_on_i_c,
					    base_exp_log(),
					    ARGUMENT_LOG);
    solution_on_i_c = simplify_ex_for_input(solution_on_i_c, true);
    solution_on_i_c = simplify_logarithm(solution_on_i_c);
  }
  else {
    // Substitute symbolic initial conditions with the values in the map
    // `initial_conditions'.
    for (std::map<index_type, Expr>::const_iterator i
	   = initial_conditions_.begin(),
	   iend = initial_conditions_.end(); i != iend; ++i)
      solution_on_i_c
	= solution_on_i_c.substitute(x(i->first),
				     get_initial_condition(i->first));
  }
  solution_on_i_c = simplify_numer_denom(solution_on_i_c);
  solution_on_i_c = simplify_ex_for_output(solution_on_i_c, false);
  return solution_on_i_c;
}

Expr
PURRS::Recurrence::
compute_solution_weighted_average_on_i_c(const Expr& solution) const {
  Expr weight_rec;
  Expr inhomogeneous;
  if (recurrence_rewritten) {
    weight_rec = original_weight();
    inhomogeneous = original_inhomogeneous();
  }
  else {
    weight_rec = weight();
    inhomogeneous = inhomogeneous_term;
  }
  Expr sol = solution;

  // Consider the maximum index of `x' function in the map
  // `initial_conditions', i.e. the largest index of the initial
  // conditions inserted by the user.
  unsigned int max_index_user_i_c = get_max_index_initial_condition();
  if (max_index_user_i_c > first_valid_index) {
    // To solve `x(n) = f(n)*sum(k,n_0,n-1,x(k))+g(n)' with the
    // initial condition `x(m) = h', where `m > n_0',
    // is like to solve `x(n) = f(n)*sum(k,m,n-1,x(k))+g(n)'.
    Symbol k;
    Recurrence rec(weight_rec * PURRS::sum(k, max_index_user_i_c, n-1, x(k))
		   + inhomogeneous);
    // The following two methods called on the recurrence `rec'
    // surely will return SUCCESS because they work always on the
    // same recurrence (shifted backward).
    rec.classify_and_catch_special_cases();
    rec.compute_weighted_average_recurrence(sol);
    sol = sol.substitute(x(max_index_user_i_c),
			 get_initial_condition(max_index_user_i_c));
  }
  else
    sol = sol.substitute(x(first_valid_index),
			 get_initial_condition(first_valid_index));
  
  return sol;
}

Expr
PURRS::Recurrence::
compute_bound_functional_equation_on_i_c(const Expr& bound) const {
  Expr bound_evaluated = bound;
  const Number& divisor = functional_eq_p->ht_begin()->first;
  if (divisor == 2)
    bound_evaluated
      = bound_evaluated.substitute(x(1), get_initial_condition(1));
  return bound_evaluated;
}

/*!
  FIXME: update the comment!!
  \param solution_or_bound  Contains the solution or the bound computed
                            for the recurrence \p *this in function of
			    arbitrary symbolic initial conditions.

  \return                   The solution or the bound shifted in agreement
                            with the initial conditions inserted by the user
			    and with the eventual remaining symbolic initial
			    conditions.

  We know the least non-negative integer \f$ s \f$ such that
  the recurrence is well-defined for \f$ n \geq s \f$.
  This function checks if in the map \p initial_conditions there are
  some initial conditions of the form \f$ x(i) = k \f$ with \f$ k > s \f$:
  in this case the function shifts the solution or the bound.
  Finally replaces the possible symbolic initial conditions in the
  solution or in the bound with the values specified by the user.
*/
Expr
PURRS::Recurrence::
compute_solution_or_bound_on_i_c(const Expr& solution_or_bound) const {
  assert(!initial_conditions_.empty());
  if (!has_at_least_a_symbolic_initial_condition(solution_or_bound))
    return solution_or_bound;
  
  Expr solution_or_bound_on_i_c;
  switch (type_) {
  case ORDER_ZERO:
    break;
  case LINEAR_FINITE_ORDER_CONST_COEFF:
  case LINEAR_FINITE_ORDER_VAR_COEFF:
    solution_or_bound_on_i_c
      = compute_solution_linear_finite_order_on_i_c(solution_or_bound);
    break;
  case NON_LINEAR_FINITE_ORDER:
    solution_or_bound_on_i_c
      = compute_solution_non_linear_finite_order_on_i_c(solution_or_bound);
    break;
  case WEIGHTED_AVERAGE:
    solution_or_bound_on_i_c
      = compute_solution_weighted_average_on_i_c(solution_or_bound);
    break;
  case FUNCTIONAL_EQUATION:
    solution_or_bound_on_i_c
      = compute_bound_functional_equation_on_i_c(solution_or_bound);
    break;
  default:
    throw std::runtime_error("PURRS internal error: "
			     "compute_solution_or_bound_on_i_c().");
  }
  return solution_or_bound_on_i_c;
}

/*!

*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution_finite_order() const {
  Solver_Status status;
  // If the greatest common divisor among the decrements is greater
  // than one, the order reduction is applicable.
  if (gcd_among_decrements() > 1)
    // FIXME: the order reduction is for the moment applied only to
    // recurrences with constant coefficients because the recurrences
    // with variable coefficients are not allowed with parameters.
    if (is_linear_finite_order_const_coeff()) {
      if ((status = apply_order_reduction()) != SUCCESS)
	return status;
    }
    else {
      if ((status = solve_linear_finite_order()) != SUCCESS)
	return status;
    }
  // We have not applied the order reduction.
  else
    if ((status = solve_linear_finite_order()) != SUCCESS)
      return status;

  lower_bound_.set_expression(exact_solution_.expression());
  upper_bound_.set_expression(exact_solution_.expression());
  
  // Check if there are specified initial conditions and in this case
  // shift the solution according to them if necessary before
  // substituting the values of the initial conditions to the
  // symbolic initial conditions `x(i)'.
  // If the order of the recurrence is `0' no initial conditions
  // must be substituted.
  if (!initial_conditions_.empty() && order() > 0) {
    evaluated_exact_solution_.set_expression
      (compute_solution_or_bound_on_i_c(exact_solution_.expression()));
    // FIXME: This ought to be done more generally in 
    // exact_solution(Expr& e) below.
    evaluated_exact_solution_
      .set_expression(evaluated_exact_solution_
		      .replace_system_generated_symbols(*this));
    evaluated_lower_bound_.set_expression
      (evaluated_exact_solution_.expression());
    evaluated_upper_bound_.set_expression
      (evaluated_exact_solution_.expression());
  }
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution_functional_equation() const {
  Solver_Status status;
  if ((status = approximate_functional_equation(LOWER)) == SUCCESS
      && (status = approximate_functional_equation(UPPER)) == SUCCESS
      && lower_bound_.expression() == upper_bound_.expression()) {
    exact_solution_.set_expression(lower_bound_.expression());

    // Check if there are specified initial conditions and in this case
    // eventually shift the solution in according with them before to
    // substitute the values of the initial conditions to the
    // symbolic initial condition `x(i)'.
    if (!initial_conditions_.empty()) {
      evaluated_lower_bound_.set_expression
	(compute_solution_or_bound_on_i_c(exact_solution_.expression()));
      evaluated_upper_bound_.set_expression
	(evaluated_lower_bound_.expression());
      evaluated_exact_solution_.set_expression
	(evaluated_lower_bound_.expression());
    }
    return SUCCESS;
  }
  else
    return TOO_COMPLEX;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution_non_linear() const {
  Expr solution;
  Solver_Status status  = compute_non_linear_recurrence(solution, 0);
  if (status != SUCCESS)
    return status;  

  exact_solution_.set_expression(solution);
  lower_bound_.set_expression(solution);
  upper_bound_.set_expression(solution);
  // FIXME: before to substitute the initial conditions in the
  // non-linear recurrences we must be sure to have sufficiently
  // simplified the solution.
  // At the moment we accept non-negative values for the initial conditions:
  // if there is the value 0 in one of the last `k' (k' is the order) then the
  // solution is 0; otherwise thre are not problems.

  // Check if there are specified initial conditions and in this case
  // eventually shift the solution in according with them before to
  // substitute the values of the initial conditions to the
  // symbolic initial condition `x(i)'.
  if (!initial_conditions_.empty()) {
    evaluated_exact_solution_.set_expression
      (compute_solution_or_bound_on_i_c(exact_solution_.expression()));
    evaluated_lower_bound_.set_expression
      (evaluated_exact_solution_.expression());
    evaluated_upper_bound_.set_expression
      (evaluated_exact_solution_.expression());
  }
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution_weighted_average() const {
  Expr solution;
  Solver_Status status = compute_weighted_average_recurrence(solution);
  if (status != SUCCESS)
    return status;

  exact_solution_.set_expression(solution);
  lower_bound_.set_expression(solution);
  upper_bound_.set_expression(solution);

  if (!initial_conditions_.empty()) {
    solution = compute_solution_or_bound_on_i_c(solution);
    evaluated_exact_solution_.set_expression(solution);
    evaluated_lower_bound_.set_expression(solution);
    evaluated_upper_bound_.set_expression(solution);
  }
  
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution_max() const {
  Expr solution;
  assert(recurrence_rhs.is_the_max_function());
  // Any well-formed expression containing a max is expressed as
  //   x(n) = max(f(x(0), ..., x(n-1)), a)
  // We interpret it as follows:
  //   x(0) = a
  //   x(n) = max(f(x(0), ..., x(n-1)), a) for all n >= 1.
  Expr numeric_or_symbolic_arg;
  Expr functional_arg;
  if (recurrence_rhs.arg(0).is_a_number() || recurrence_rhs.arg(0).is_a_symbol()) {
    numeric_or_symbolic_arg = recurrence_rhs.arg(0);
    functional_arg = recurrence_rhs.arg(1);
  }
  else if (recurrence_rhs.arg(1).is_a_number() || recurrence_rhs.arg(1).is_a_symbol()) {
    numeric_or_symbolic_arg = recurrence_rhs.arg(1);
    functional_arg = recurrence_rhs.arg(0);
  }
  else return TOO_COMPLEX;

  Recurrence aux_rec(functional_arg);

  aux_rec.classify();
  if (!aux_rec.is_linear_finite_order())
    return TOO_COMPLEX;
  
  // FIXME: check that the coefficients are positive.
  Solver_Status status = aux_rec.compute_exact_solution();
  if (status != SUCCESS)
    return status;


  if (!initial_conditions_.empty()) {
    solution = compute_solution_or_bound_on_i_c(solution);
    evaluated_exact_solution_.set_expression(solution);
    evaluated_lower_bound_.set_expression(solution);
    evaluated_upper_bound_.set_expression(solution);
  }
  
  std::map<index_type, Expr> initial_conditions;
  initial_conditions[0] = numeric_or_symbolic_arg;
  aux_rec.set_initial_conditions(initial_conditions);
  aux_rec.exact_solution(solution);

  exact_solution_.set_expression(solution);
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution() const {
  // It can happen that there is not the exact solution although
  // the system tried to compute it (for example the recurrence is
  // too complex): in order to avoid to repeat the attempt of
  // computation of the solution, we set to `true' the data
  // `tried_to_compute_exact_solution'.
  tried_to_compute_exact_solution = true;

  // See if we already have the exact solution.
  if (exact_solution_.has_expression()) {
    // Check if there are specified initial conditions and in this case
    // eventually shift the solution in according with them before to
    // substitute the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    if (!initial_conditions_.empty()
	&& !evaluated_exact_solution_.has_expression()) {
      evaluated_exact_solution_.set_expression
	(compute_solution_or_bound_on_i_c(exact_solution_.expression()));
      evaluated_lower_bound_.set_expression
	(evaluated_exact_solution_.expression());
      evaluated_upper_bound_.set_expression
	(evaluated_exact_solution_.expression());
    }
    return SUCCESS;
  }

  // We may not have the exact solution explicitely, yet we may have
  // the lower and the upper bounds equal among themselves.
  // FIXME: invece di == usare quella funzione che torna tre valori;
  // sono uguali, sono diversi, non lo so.
  if (lower_bound_.has_expression() && upper_bound_.has_expression()
      && lower_bound_.expression() == upper_bound_.expression()) {
    exact_solution_.set_expression(lower_bound_.expression());
    return SUCCESS;
  }

  Classifier_Status classifier_status = classify_and_catch_special_cases();
  if (classifier_status == CL_SUCCESS)
    switch (type_) {
    case ORDER_ZERO:
    case LINEAR_FINITE_ORDER_CONST_COEFF:
    case LINEAR_FINITE_ORDER_VAR_COEFF:
      return compute_exact_solution_finite_order();
    case FUNCTIONAL_EQUATION:
      return compute_exact_solution_functional_equation();
    case NON_LINEAR_FINITE_ORDER:
      return compute_exact_solution_non_linear();
    case WEIGHTED_AVERAGE:
      return compute_exact_solution_weighted_average();
    case MAX_FUNCTION:
      return compute_exact_solution_max();
    default:
      throw std::runtime_error("PURRS internal error: "
			       "compute_exact_solution().");
    }
  else
    // return the `Solver_Status' associated to `classifier_status'.
    return map_status(classifier_status);
}

void
PURRS::Recurrence::exact_solution(Expr& e) const {
  if (!exact_solution_.has_expression())
    throw std::logic_error("PURRS::Recurrence::exact_solution() called, "
			   "but no exact solution was computed");
  if (initial_conditions_.empty()) {
    // Substitutes all symbols generated by the system contained in
    // \p (*this).expression() with new symbols with shorter
    // name that are not yet used.
    exact_solution_
      .set_expression(exact_solution_.replace_system_generated_symbols(*this));
    e = exact_solution_.expression();
  }
  else {
    // If `evaluated_exact_solution_' is set means that the user
    // has specified the initial conditions. In this case we
    // possibly shift the solution in according with the initial conditions
    // before replacing the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    // FIXME: Always replace system generated symbols and consider that
    // evaluated_exact_solution_ might have been given an expression above.
    if (!evaluated_exact_solution_.has_expression()) {
      evaluated_exact_solution_.set_expression
	(compute_solution_or_bound_on_i_c(exact_solution_.expression()));
      evaluated_exact_solution_
      .set_expression(evaluated_exact_solution_
		      .replace_system_generated_symbols(*this));
      evaluated_lower_bound_.set_expression
	(evaluated_exact_solution_.expression());
      evaluated_upper_bound_.set_expression
	(evaluated_exact_solution_.expression());
    }
    e = evaluated_exact_solution_.expression();
  }
  set_first_valid_index_for_solution();
  assert(has_only_symbolic_initial_conditions(e));
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::
compute_bound_functional_equation(Bound kind_of_bound) const {
  Solver_Status status = approximate_functional_equation(kind_of_bound);
  if (status != SUCCESS)
    return status;

  // Check if there are specified initial conditions and in this case
  // eventually shift the solution in according with them before to
  // substitute the values of the initial conditions to the
  // symbolic initial condition `x(i)'.
  if (!initial_conditions_.empty())
    if (kind_of_bound == LOWER)
      evaluated_lower_bound_.set_expression
	(compute_solution_or_bound_on_i_c(lower_bound_.expression()));
    else
      evaluated_upper_bound_.set_expression
	(compute_solution_or_bound_on_i_c(upper_bound_.expression()));
  
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::
compute_bound_non_linear(Bound kind_of_bound) const {
  assert(kind_of_bound == LOWER || kind_of_bound == UPPER);
  Expr bound;
  unsigned int type = kind_of_bound == LOWER ? 1 : 2;
  Solver_Status status = compute_non_linear_recurrence(bound, type);
  if (status != SUCCESS)
    return status;

  if (kind_of_bound == LOWER)
    lower_bound_.set_expression(bound);
  else
    upper_bound_.set_expression(bound);
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_lower_bound() const {
  // See if we already have the lower bound.
  if (lower_bound_.has_expression()) {
    // Check if there are specified initial conditions and in this case
    // eventually shift the solution in according with them before to
    // substitute the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    if (!initial_conditions_.empty()
	&& !evaluated_lower_bound_.has_expression())
      evaluated_lower_bound_.set_expression
	(compute_solution_or_bound_on_i_c(lower_bound_.expression()));
    return SUCCESS;
  }

  // We may not have the lower bound explicitely, yet we may have
  // the exact solution.
  if (exact_solution_.has_expression()) {
    // Check if there are specified initial conditions and in this case
    // eventually shift the solution in according with them before to
    // substitute the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    if (!initial_conditions_.empty()
	&& !evaluated_lower_bound_.has_expression()) {
      evaluated_exact_solution_.set_expression
	(compute_solution_or_bound_on_i_c(exact_solution_.expression()));
      evaluated_lower_bound_.set_expression
	(evaluated_exact_solution_.expression());
    }
    lower_bound_.set_expression(exact_solution_.expression());
    return SUCCESS;
  }
  
  Classifier_Status classifier_status = classify_and_catch_special_cases();
  if (classifier_status == CL_SUCCESS)
    switch (type_) {
    case ORDER_ZERO:
    case LINEAR_FINITE_ORDER_CONST_COEFF:
    case LINEAR_FINITE_ORDER_VAR_COEFF:
    case WEIGHTED_AVERAGE:
      if (!tried_to_compute_exact_solution)
	return compute_exact_solution();
      else
	return TOO_COMPLEX;
      break;
    case FUNCTIONAL_EQUATION:
      return compute_bound_functional_equation(LOWER); 
      break;
    case NON_LINEAR_FINITE_ORDER:
      return compute_bound_non_linear(LOWER); 
      break;
    case MAX_FUNCTION:
      return TOO_COMPLEX;
      break;
    default:
      throw std::runtime_error("PURRS internal error: "
			       "compute_lower_bound().");
    }
  else
    // return the `Solver_Status' associated to `classifier_status'.
    return map_status(classifier_status);
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_upper_bound() const {
  // See if we already have the upper bound.
  if (upper_bound_.has_expression()) {
    // Check if there are specified initial conditions and in this case
    // eventually shift the solution in according with them before to
    // substitute the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    if (!initial_conditions_.empty()
	&& !evaluated_upper_bound_.has_expression())
      evaluated_upper_bound_.set_expression
	(compute_solution_or_bound_on_i_c(upper_bound_.expression()));
    return SUCCESS;
  }

  // We may not have the upper bound explicitely, yet we may have
  // the exact solution.
  if (exact_solution_.has_expression()) {
    // Check if there are specified initial conditions and in this case
    // eventually shift the solution in according with them before to
    // substitute the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    if (!initial_conditions_.empty()
	&& !evaluated_upper_bound_.has_expression()) {
      evaluated_exact_solution_.set_expression
	(compute_solution_or_bound_on_i_c(exact_solution_.expression()));
      evaluated_upper_bound_.set_expression
	(evaluated_exact_solution_.expression());
    }
    upper_bound_.set_expression(exact_solution_.expression());
    return SUCCESS;
  }

  Classifier_Status classifier_status = classify_and_catch_special_cases();
  if (classifier_status == CL_SUCCESS)
    switch (type_) {
    case ORDER_ZERO:
    case LINEAR_FINITE_ORDER_CONST_COEFF:
    case LINEAR_FINITE_ORDER_VAR_COEFF:
    case WEIGHTED_AVERAGE:
      if (!tried_to_compute_exact_solution)
	return compute_exact_solution();
      else
	return TOO_COMPLEX;
      break;
    case FUNCTIONAL_EQUATION:
      return compute_bound_functional_equation(UPPER); 
      break;
    case NON_LINEAR_FINITE_ORDER:
      return compute_bound_non_linear(UPPER);
      break;
    case MAX_FUNCTION:
      return TOO_COMPLEX;
      break;
    default:
      throw std::runtime_error("PURRS internal error: "
			       "compute_upper_bound().");
    }
  else
    // return the `Solver_Status' associated to `classifier_status'.
    return map_status(classifier_status);
}

void
PURRS::Recurrence::lower_bound(Expr& e) const {
  if (!lower_bound_.has_expression())
    throw std::logic_error("PURRS::Recurrence::lower_bound() called, "
			   "but no lower bounds was computed");
  if (initial_conditions_.empty()) {
    // Substitutes all symbols generated by the system contained in
    // \p (*this).expression() with new symbols with shorter
    // name that are not yet used.
    lower_bound_
      .set_expression(lower_bound_.replace_system_generated_symbols(*this));
    e = lower_bound_.expression();
  }
  else {
    // If `evaluated_lower_bound_' is set means that the user
    // has specified the initial conditions. In this case we
    // possibly shift the lower bound in according with the initial conditions
    // before replacing the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    if (!evaluated_lower_bound_.has_expression()) {
      evaluated_lower_bound_.set_expression
        (compute_solution_or_bound_on_i_c(lower_bound_.expression()));
      evaluated_lower_bound_
      .set_expression(evaluated_lower_bound_
		      .replace_system_generated_symbols(*this));
    }
    e = evaluated_lower_bound_.expression();
  }
  set_first_valid_index_for_solution();
  if (is_functional_equation()
      || (is_non_linear_finite_order()
	  && associated_linear_rec().is_functional_equation()))
    set_definition_Sc();
  assert(has_only_symbolic_initial_conditions(e));
}

void
PURRS::Recurrence::upper_bound(Expr& e) const {
  if (!upper_bound_.has_expression())
    throw std::logic_error("PURRS::Recurrence::upper_bound() called, "
			   "but no upper bounds was computed");
  if (initial_conditions_.empty()) {
    // Substitutes all symbols generated by the system contained in
    // \p (*this).expression() with new symbols with shorter
    // name that are not yet used.
    upper_bound_
      .set_expression(upper_bound_.replace_system_generated_symbols(*this));
    e = upper_bound_.expression();
  }
  else {
    // If `evaluated_upper_bound_' is set means that the user
    // has specified the initial conditions. In this case we
    // possibly shift the upper bound in according with the initial conditions
    // before replacing the values of the initial conditions to the
    // symbolic initial conditions `x(i)'.
    if (!evaluated_upper_bound_.has_expression()) {
      evaluated_upper_bound_.set_expression
	(compute_solution_or_bound_on_i_c(upper_bound_.expression()));
      evaluated_upper_bound_
	.set_expression(evaluated_upper_bound_
			.replace_system_generated_symbols(*this));
    }
    e = evaluated_upper_bound_.expression();
  }
  set_first_valid_index_for_solution();
  if (is_functional_equation()
      || (is_non_linear_finite_order()
	  && associated_linear_rec().is_functional_equation()))
    set_definition_Sc();
  assert(has_only_symbolic_initial_conditions(e));
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
  s << "recurrence_rewritten = "
    << (recurrence_rewritten ? "true" : "false") << std::endl;
  s << "recurrence_rhs = " << recurrence_rhs << std::endl;
  s << "first_valid_index = " << first_valid_index << std::endl;
  if (exact_solution_.has_expression()
      || lower_bound_.has_expression() || upper_bound_.has_expression())
    s << "first_valid_index_for_solution = "
      << first_valid_index_for_solution() << std::endl;
  if ((lower_bound_.has_expression() || upper_bound_.has_expression())
      && !exact_solution_.has_expression()) {
    std::string Sc_function = definition_Sc();
    if (!Sc_function.empty())
      s << Sc_function << std::endl;
  }
  s << "auxiliary_definitions:" << std::endl;
  blackboard.dump(s);
  
  if (!initial_conditions_.empty()) {
    s << "Initial conditions:" << std::endl;
    for (std::map<index_type, Expr>::const_iterator i
	   = initial_conditions_.begin(),
	   initial_conditions_end = initial_conditions_.end();
	 i != initial_conditions_end; ++i)
      s << "  x(" << i->first << ")"
	<< " = " << i->second << std::endl;
  }
  
  //(*functional_eq_p).dump_homogeneous_terms(s);
  s << std::endl;
}
