/* Recurrence class implementation: inline functions.
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

#ifndef PURRS_Recurrence_inlines_hh
#define PURRS_Recurrence_inlines_hh

#include "Non_Linear_Info.defs.hh"
#include "Weighted_Average_Info.defs.hh"
#include "Expr.types.hh"

#include <iostream>
#include <utility>
#include <stdexcept>
#include <algorithm>

namespace Parma_Recurrence_Relation_Solver {

inline
Recurrence::Recurrence()
  : classifier_status_(CL_SUCCESS),
    recurrence_rhs(0),
    recurrence_rewritten(false),
    inhomogeneous_term(0),
    type_(ORDER_ZERO),
    finite_order_p(0),
    functional_eq_p(0),
    non_linear_p(0),
    weighted_average_p(0),
    tried_to_compute_exact_solution(false) {
}

inline
Recurrence::Recurrence(const Expr& e)
  : classifier_status_(NOT_CLASSIFIED_YET),
    recurrence_rhs(e),
    recurrence_rewritten(false),
    inhomogeneous_term(0),
    type_(UNKNOWN),
    finite_order_p(0),
    functional_eq_p(0),
    non_linear_p(0),
    weighted_average_p(0),
    tried_to_compute_exact_solution(false) {
}

inline
Recurrence::Recurrence(const Recurrence& y)
  : classifier_status_(y.classifier_status_),
    recurrence_rhs(y.recurrence_rhs),
    recurrence_rewritten(y.recurrence_rewritten),
    inhomogeneous_term(y.inhomogeneous_term),
    system_rhs(y.system_rhs),
    type_(y.type_),
    finite_order_p(y.finite_order_p),
    functional_eq_p(y.functional_eq_p),    
    non_linear_p(y.non_linear_p),
    weighted_average_p(y.weighted_average_p),
    exact_solution_(y.exact_solution_),
    lower_bound_(y.lower_bound_),
    upper_bound_(y.upper_bound_),
    tried_to_compute_exact_solution(y.tried_to_compute_exact_solution),
    blackboard(y.blackboard),
    initial_conditions_(y.initial_conditions_) {
}

inline Recurrence&
Recurrence::operator=(const Recurrence& y) {
  classifier_status_ = y.classifier_status_;
  recurrence_rhs = y.recurrence_rhs;
  recurrence_rewritten = y.recurrence_rewritten;
  inhomogeneous_term = y.inhomogeneous_term;
  system_rhs = y.system_rhs;
  type_ = y.type_;
  finite_order_p = y.finite_order_p;
  functional_eq_p = y.functional_eq_p;
  non_linear_p = y.non_linear_p;
  weighted_average_p = y.weighted_average_p;
  exact_solution_ = y.exact_solution_;
  lower_bound_ = y.lower_bound_;
  upper_bound_ = y.upper_bound_;
  tried_to_compute_exact_solution = y.tried_to_compute_exact_solution;
  blackboard = y.blackboard;
  initial_conditions_ = y.initial_conditions_;
  return *this;
}

inline
Recurrence::~Recurrence() {
  delete finite_order_p;
  delete functional_eq_p;
  delete non_linear_p;
  delete weighted_average_p;
}

inline void
Recurrence::replace_recurrence(const Expr& e) {
  recurrence_rhs = e;
  classifier_status_ = NOT_CLASSIFIED_YET;
}

inline void
Recurrence::replace_recurrence(unsigned int k, const Expr& e) {
  std::pair<std::map<unsigned int, Expr>::iterator, bool> stat
    = system_rhs.insert(std::map<unsigned int, Expr>::value_type(k, e));
  if (!stat.second)
    // There was already something associated to `k': overwrite it.
    stat.first->second = e;
}

/*!
  The \f$ 3 \$ <CODE>Cached_Expr</CODE> that contain the exact solution,
  the lower bound and the upper bound evaluated on the initial conditions
  stored in the map <CODE>initial_conditions_</CODE> are unset.
  This technique allows to know if the user wants the solution (or the bound)
  with symbolic initial conditions or the solution (or the bounds)
  evaluated simply checking if the expression containing the solution
  (or the bound) evaluated is set.
*/
inline void
Recurrence::reset_initial_conditions() {
  initial_conditions_.clear();
  evaluated_exact_solution_.unset_expression();
  evaluated_lower_bound_.unset_expression();
  evaluated_upper_bound_.unset_expression();
}

inline const std::map<index_type, Expr>&
Recurrence::get_initial_conditions() const {
  return initial_conditions_;
}

inline Expr
Recurrence::get_initial_condition(unsigned int k) const {
  std::map<index_type, Expr>::const_iterator i = initial_conditions_.find(k);
  if (i != initial_conditions_.end())
    return (*i).second;
  else
    return x(k);
}

inline unsigned int
Recurrence::get_max_index_initial_condition() const {
  unsigned int max_index;
  if (initial_conditions_.empty())
    max_index = 0;
  else
    max_index = initial_conditions_.rbegin()->first;
  return max_index;
}

inline Expr
Recurrence::evaluate_exact_solution(const Number& num) const {
  // `0' means that is the case of exact solution.
  return evaluate(0, num);
}

template <class OutputIterator>
void
Recurrence::evaluate_exact_solution(const Number& l, const Number& r,
				    OutputIterator oi) const {
  if (l > r)
    throw_invalid_argument("PURRS::Recurrence::evaluate_exact_solution()",
			   "the numbers represent an interval and the first\n"
			   "must be smaller or equal to the second");
  for (Number i = l; i <= r; ++i) {
    // `0' means that is the case of exact solution.
    *oi = evaluate(0, i);
    ++oi;
  }
}

inline Expr
Recurrence::evaluate_lower_bound(const Number& num) const {
  // `1' means that is the case of lower bound.
  return evaluate(1, num);
}

template <class OutputIterator>
void
Recurrence::evaluate_lower_bound(const Number& l, const Number& r,
				 OutputIterator oi) const {
  if (l > r)
    throw_invalid_argument("PURRS::Recurrence::evaluate_lower_bound()",
			   "the numbers represent an interval and the first\n"
			   "must be smaller or equal to the second");
  for (Number i = l; i <= r; ++i) {
    // `1' means that is the case of lower bound.
    *oi = evaluate(1, i);
    ++oi;
  }
}

inline Expr
Recurrence::evaluate_upper_bound(const Number& num) const {
  // `2' means that is the case of upper bound.
  return evaluate(2, num);
}

template <class OutputIterator>
void
Recurrence::evaluate_upper_bound(const Number& l, const Number& r,
				 OutputIterator oi) const {
  if (l > r)
    throw_invalid_argument("PURRS::Recurrence::evaluate_upper_bound()",
			   "the numbers represent an interval and the first\n"
			   "must be smaller or equal to the second");
  for (Number i = l; i <= r; ++i) {
    // `2' means that is the case of upper bound.
    *oi = evaluate(2, i);
    ++oi;
  }
}

inline void
Recurrence::set_inhomogeneous_term(const Expr& e) const {
  inhomogeneous_term = e;
}

inline const Recurrence::Type&
Recurrence::type() const {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return type_;
}

inline Recurrence::Type&
Recurrence::type() {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return type_;
}

inline void
Recurrence::set_type(const Type& t) const {
  type_ = t;
  classifier_status_ = CL_SUCCESS;
}

inline void
Recurrence::set_order_zero() const {
  type_ = ORDER_ZERO; 
  classifier_status_ = CL_SUCCESS;
}

inline bool
Recurrence::is_linear_finite_order_const_coeff() const {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return type_ == LINEAR_FINITE_ORDER_CONST_COEFF;
}

inline void
Recurrence::set_linear_finite_order_const_coeff() const {
  type_ = LINEAR_FINITE_ORDER_CONST_COEFF;
  classifier_status_ = CL_SUCCESS;
}

inline bool
Recurrence::is_linear_finite_order_var_coeff() const {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return type_ == LINEAR_FINITE_ORDER_VAR_COEFF;
}

inline void
Recurrence::set_linear_finite_order_var_coeff() const {
  type_ = LINEAR_FINITE_ORDER_VAR_COEFF;
  classifier_status_ = CL_SUCCESS;
}

inline bool
Recurrence::is_linear_finite_order() const {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return (type_ == ORDER_ZERO
	  || type_ == LINEAR_FINITE_ORDER_CONST_COEFF
	  || type_ == LINEAR_FINITE_ORDER_VAR_COEFF);
}

inline bool
Recurrence::is_non_linear_finite_order() const {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return type_ == NON_LINEAR_FINITE_ORDER;
}

inline void
Recurrence::set_non_linear_finite_order() const {
  type_ = NON_LINEAR_FINITE_ORDER;
  classifier_status_ = CL_SUCCESS;
}

inline bool
Recurrence::is_functional_equation() const {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return type_ == FUNCTIONAL_EQUATION;
}

inline void
Recurrence::set_functional_equation() const {
  type_ = FUNCTIONAL_EQUATION;
  classifier_status_ = CL_SUCCESS;
}

inline bool
Recurrence::is_weighted_average() const {
  assert(classifier_status_ != NOT_CLASSIFIED_YET);
  return type_ == WEIGHTED_AVERAGE;
}

inline void
Recurrence::set_weighted_average() const {
  type_ = WEIGHTED_AVERAGE;
  classifier_status_ = CL_SUCCESS;
}

inline index_type
Recurrence::order() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->order();
}

inline index_type
Recurrence::first_valid_index() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->first_valid_index();
}

inline void
Recurrence::set_first_valid_index(index_type i_c) const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  finite_order_p->set_first_valid_index(i_c);
}

inline const std::vector<Expr>&
Recurrence::coefficients() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->coefficients();
}

inline std::vector<Expr>&
Recurrence::coefficients() {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->coefficients();
}

inline unsigned int
Recurrence::gcd_among_decrements() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->gcd_among_decrements();
}

inline const Expr&
Recurrence::product_factor() const {
  assert(is_linear_finite_order_var_coeff());
  assert(finite_order_p);
  return finite_order_p->product_factor();
}

inline Expr&
Recurrence::product_factor() {
  assert(is_linear_finite_order_var_coeff());
  assert(finite_order_p);
  return finite_order_p->product_factor();
}

inline void
Recurrence::set_product_factor(const Expr& x) const {
  assert(is_linear_finite_order_var_coeff());
  assert(finite_order_p);
  return finite_order_p->set_product_factor(x);
}

inline bool 
Recurrence::applied_order_reduction() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->applied_order_reduction();
}

inline void 
Recurrence::set_order_reduction() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->set_order_reduction();
}

inline void 
Recurrence::unset_order_reduction() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p->unset_order_reduction();
}

inline index_type
Recurrence::applicability_condition() const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p->applicability_condition();
}

inline void
Recurrence::set_applicability_condition(index_type c) const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p->set_applicability_condition(c);
}

inline index_type
Recurrence::rank() const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p->rank();
}

inline const Recurrence&
Recurrence::associated_linear_rec() const {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->associated_linear_rec();
}

inline Recurrence&
Recurrence::associated_linear_rec() {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->associated_linear_rec();
}

inline const Number&
Recurrence::coeff_simple_non_linear_rec() const {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->coeff_simple_non_linear_rec();
}

inline Number&
Recurrence::coeff_simple_non_linear_rec() {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->coeff_simple_non_linear_rec();
}

inline const Expr&
Recurrence::base_exp_log() const {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->base_exp_log();
}

inline Expr&
Recurrence::base_exp_log() {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->base_exp_log();
}

inline const std::vector<Symbol>&
Recurrence::auxiliary_symbols() const {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->auxiliary_symbols();
}

inline std::vector<Symbol>&
Recurrence::auxiliary_symbols() {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p->auxiliary_symbols();
}

inline const Recurrence&
Recurrence::associated_first_order_rec() const {
  assert(is_weighted_average());
  return weighted_average_p->associated_first_order_rec();
}

inline Recurrence&
Recurrence::associated_first_order_rec() {
  assert(is_weighted_average());
  return weighted_average_p->associated_first_order_rec();
}

inline void
Recurrence::set_original_rhs(const Expr& original_rhs) const {
  assert(is_weighted_average());
  weighted_average_p->set_original_rhs(original_rhs);
}

inline unsigned int
Recurrence::lower_limit() const {
  assert(is_weighted_average());
  assert(weighted_average_p);
  return weighted_average_p->lower_limit();
}

inline void
Recurrence::set_lower_limit(unsigned int lower) const {
  assert(is_weighted_average());
  weighted_average_p->set_lower_limit(lower);
}

inline const Expr&
Recurrence::weight() const {
  assert(is_weighted_average());
  assert(weighted_average_p);
  return weighted_average_p->weight();
}

inline Expr&
Recurrence::weight() {
  assert(is_weighted_average());
  assert(weighted_average_p);
  return weighted_average_p->weight();
}

inline Expr
Recurrence::approximated_solution() const {
  if (exact_solution_.has_expression())
    return blackboard.approximate(exact_solution_.expression());
  else
    // Well, if the client insists...
    return blackboard.approximate(recurrence_rhs);
}

inline Symbol
Recurrence::insert_auxiliary_definition(const Expr& e) const {
  return blackboard.insert_definition(e);
}

inline Expr
Recurrence::get_auxiliary_definition(const Symbol& z) const {
  return blackboard.get_definition(z);
}

inline Expr
Recurrence::substitute_auxiliary_definitions(const Expr& e) const {
  return blackboard.rewrite(e);
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Recurrence_inlines_hh)
