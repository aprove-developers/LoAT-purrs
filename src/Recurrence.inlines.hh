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

#include <iostream>
#include <utility>
#include <stdexcept>
#include <algorithm>

namespace Parma_Recurrence_Relation_Solver {

inline
Recurrence::Recurrence()
  : recurrence_rhs(0),
    recurrence_rewritten(false),
    come_from_non_linear_rec(false),
    inhomogeneous_term(0),
    type_(ORDER_ZERO),
    is_classified(false),
    finite_order_p(0),
    functional_eq_p(0),
    non_linear_p(0),
    infinite_order_p(0),
    tested_exact_solution(false) {
}

inline
Recurrence::Recurrence(const Expr& e)
  : recurrence_rhs(e),
    recurrence_rewritten(false),
    come_from_non_linear_rec(false),
    inhomogeneous_term(0),
    type_(ORDER_ZERO),
    is_classified(false),
    finite_order_p(0),
    functional_eq_p(0),
    non_linear_p(0),
    infinite_order_p(0),
    tested_exact_solution(false) {
}

inline
Recurrence::Recurrence(const Recurrence& y)
  : recurrence_rhs(y.recurrence_rhs),
    recurrence_rewritten(y.recurrence_rewritten),
    come_from_non_linear_rec(y.come_from_non_linear_rec),
    inhomogeneous_term(y.inhomogeneous_term),
    system_rhs(y.system_rhs),
    type_(y.type_),
    is_classified(y.is_classified),
    finite_order_p(y.finite_order_p),
    functional_eq_p(y.functional_eq_p),    
    non_linear_p(y.non_linear_p),
    infinite_order_p(y.infinite_order_p),
    exact_solution_(y.exact_solution_),
    lower_bound_(y.lower_bound_),
    upper_bound_(y.upper_bound_),
    tested_exact_solution(y.tested_exact_solution),
    blackboard(y.blackboard),
    initial_conditions(y.initial_conditions) {
}

inline
Recurrence::~Recurrence() {
  delete finite_order_p;
  delete functional_eq_p;
  delete non_linear_p;
  delete infinite_order_p;
}

inline Recurrence&
Recurrence::operator=(const Recurrence& y) {
  recurrence_rhs = y.recurrence_rhs;
  recurrence_rewritten = y.recurrence_rewritten;
  come_from_non_linear_rec = y.come_from_non_linear_rec;
  inhomogeneous_term = y.inhomogeneous_term;
  system_rhs = y.system_rhs;
  type_ = y.type_;
  is_classified = y.is_classified;
  finite_order_p = y.finite_order_p;
  functional_eq_p = y.functional_eq_p;
  non_linear_p = y.non_linear_p;
  infinite_order_p = y.infinite_order_p;
  exact_solution_ = y.exact_solution_;
  lower_bound_ = y.lower_bound_;
  upper_bound_ = y.upper_bound_;
  tested_exact_solution = y.tested_exact_solution;
  blackboard = y.blackboard;
  initial_conditions = y.initial_conditions;
  return *this;
}

inline void
Recurrence::replace_recurrence(const Expr& e) {
  recurrence_rhs = e;
}

inline void
Recurrence::replace_recurrence(unsigned k, const Expr& e) {
  std::pair<std::map<unsigned, Expr>::iterator, bool> stat
    = system_rhs.insert(std::map<unsigned, Expr>::value_type(k, e));
  if (!stat.second)
    // There was already something associated to `k': overwrite it.
    stat.first->second = e;
}

inline void
Recurrence::replace_initial_condition(unsigned k, const Expr& e) {
  std::pair<std::map<unsigned, Expr>::iterator, bool> stat
    = initial_conditions.insert(std::map<unsigned, Expr>::value_type(k, e));
  if (!stat.second)
    // There was already something associated to `i': overwrite it.
    stat.first->second = e;
}

inline bool
Recurrence::
undefined_initial_conditions(std::set<unsigned int>& undefined) const {
  return false; 
  /*
  assert(is_linear_finite_order() || is_functional_equation());
  if (is_linear_finite_order()) {
  }
  else {
  }
  */
}

inline Expr
Recurrence::get_initial_condition(unsigned k) const {
  std::map<unsigned, Expr>::const_iterator i = initial_conditions.find(k);
  if (i != initial_conditions.end())
    return (*i).second;
  else
    return x(k);
}

inline void
Recurrence::set_inhomogeneous_term(const Expr& e) const {
  inhomogeneous_term = e;
}

inline Recurrence::Type
Recurrence::type() const {
  return type_;
}

inline Recurrence::Type&
Recurrence::type() {
  return type_;
}

inline void
Recurrence::set_type(const Type& t) const {
  type_ = t;
}

inline bool
Recurrence::is_order_zero() const {
  return type_ == ORDER_ZERO; 
}

inline void
Recurrence::set_order_zero() const {
  type_ = ORDER_ZERO; 
}

inline bool
Recurrence::is_linear_finite_order_const_coeff() const {
  return type_ == LINEAR_FINITE_ORDER_CONST_COEFF;
}

inline void
Recurrence::set_linear_finite_order_const_coeff() const {
  type_ = LINEAR_FINITE_ORDER_CONST_COEFF;
}

inline bool
Recurrence::is_linear_finite_order_var_coeff() const {
  return type_ == LINEAR_FINITE_ORDER_VAR_COEFF;
}

inline void
Recurrence::set_linear_finite_order_var_coeff() const {
  type_ = LINEAR_FINITE_ORDER_VAR_COEFF;
}

inline bool
Recurrence::is_linear_finite_order() const {
  return (type_ == ORDER_ZERO
	  || type_ == LINEAR_FINITE_ORDER_CONST_COEFF
	  || type_ == LINEAR_FINITE_ORDER_VAR_COEFF);
}

inline bool
Recurrence::is_non_linear_finite_order() const {
  return type_ == NON_LINEAR_FINITE_ORDER;
}

inline void
Recurrence::set_non_linear_finite_order() const {
  type_ = NON_LINEAR_FINITE_ORDER;
}

inline bool
Recurrence::is_functional_equation() const {
  return type_ == FUNCTIONAL_EQUATION;
}

inline void
Recurrence::set_functional_equation() const {
  type_ = FUNCTIONAL_EQUATION;
}

inline bool
Recurrence::is_linear_infinite_order() const {
  return type_ == LINEAR_INFINITE_ORDER;
}

inline void
Recurrence::set_linear_infinite_order() const {
  type_ = LINEAR_INFINITE_ORDER;
}

inline unsigned int
Recurrence::order() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> order();
}

inline unsigned
Recurrence::first_well_defined_rhs_linear() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> first_well_defined_rhs_linear();
}

inline void
Recurrence::set_first_well_defined_rhs_linear(unsigned i_c) const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  finite_order_p -> set_first_well_defined_rhs_linear(i_c);
}

inline const std::vector<Expr>&
Recurrence::coefficients() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> coefficients();
}

inline std::vector<Expr>&
Recurrence::coefficients() {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> coefficients();
}

inline unsigned
Recurrence::gcd_among_decrements() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> gcd_among_decrements();
}

inline Expr&
Recurrence::product_factor() {
  assert(is_linear_finite_order_var_coeff());
  assert(finite_order_p);
  return finite_order_p -> product_factor();
}

inline Expr
Recurrence::product_factor() const {
  assert(is_linear_finite_order_var_coeff());
  assert(finite_order_p);
  return finite_order_p -> product_factor();
}

inline void
Recurrence::set_product_factor(const Expr& x) const {
  assert(is_linear_finite_order_var_coeff());
  assert(finite_order_p);
  return finite_order_p -> set_product_factor(x);
}

inline bool 
Recurrence::applied_order_reduction() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> applied_order_reduction();
}

inline void 
Recurrence::set_order_reduction() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> set_order_reduction();
}

inline void 
Recurrence::unset_order_reduction() const {
  assert(is_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> unset_order_reduction();
}

inline unsigned
Recurrence::applicability_condition() const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p -> applicability_condition();
}

inline void
Recurrence::set_applicability_condition(unsigned c) const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p -> set_applicability_condition(c);
}

inline size_t
Recurrence::rank() const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p -> rank();
}

inline Expr
Recurrence::rhs_transformed_in_linear() const {
  assert(non_linear_p);
  return non_linear_p -> rhs_transformed_in_linear();
}

inline Expr&
Recurrence::rhs_transformed_in_linear() {
  assert(non_linear_p);
  return non_linear_p -> rhs_transformed_in_linear();
}

inline Expr
Recurrence::base_exp_log() const {
  assert(non_linear_p);
  return non_linear_p -> base_exp_log();
}

inline Expr&
Recurrence::base_exp_log() {
  assert(non_linear_p);
  return non_linear_p -> base_exp_log();
}

inline const std::vector<Symbol>&
Recurrence::auxiliary_symbols() const {
  assert(non_linear_p);
  return non_linear_p -> auxiliary_symbols();
}

inline std::vector<Symbol>&
Recurrence::auxiliary_symbols() {
  assert(non_linear_p);
  return non_linear_p -> auxiliary_symbols();
}

inline unsigned int
Recurrence::order_if_linear() const {
  assert(non_linear_p);
  return non_linear_p -> order_if_linear();
}

inline void
Recurrence::set_order_if_linear(unsigned int x) const {
  assert(non_linear_p);
  return non_linear_p -> set_order_if_linear(x);
}

inline unsigned
Recurrence::non_linear_to_linear_fwdr() const {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  return non_linear_p -> non_linear_to_linear_fwdr();
}

inline void
Recurrence::set_non_linear_to_linear_fwdr(unsigned i_c) const {
  assert(is_non_linear_finite_order());
  assert(non_linear_p);
  non_linear_p -> set_non_linear_to_linear_fwdr(i_c);
}

inline Expr
Recurrence::rhs_transformed_in_first_order() const {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> rhs_transformed_in_first_order();
}

inline Expr&
Recurrence::rhs_transformed_in_first_order() {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> rhs_transformed_in_first_order();
}

inline Expr
Recurrence::coeff_first_order() const {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> coeff_first_order();
}

inline Expr&
Recurrence::coeff_first_order() {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> coeff_first_order();
}

inline Expr
Recurrence::inhomog_first_order() const {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> inhomog_first_order();
}

inline Expr&
Recurrence::inhomog_first_order() {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> inhomog_first_order();
}

inline Expr
Recurrence::weight_inf_order() const {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> weight_inf_order();
}

inline Expr&
Recurrence::weight_inf_order() {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> weight_inf_order();
}

inline unsigned
Recurrence::infinite_order_fwdr() const {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  return infinite_order_p -> infinite_order_fwdr();
}

inline void
Recurrence::set_infinite_order_fwdr(unsigned i_c) const {
  assert(is_linear_infinite_order());
  assert(infinite_order_p);
  infinite_order_p -> set_infinite_order_fwdr(i_c);
}

inline void
Recurrence::exact_solution(Expr& e) const {
  if (!exact_solution_.has_expression())
    throw std::logic_error("PURRS::Recurrence::exact_solution() called, "
			   "but no exact solution were computed");
  exact_solution_
    .set_expression(exact_solution_.replace_system_generated_symbols(*this));
  e = exact_solution_.expression();
}

inline void
Recurrence::lower_bound(Expr& e) const {
  if (!lower_bound_.has_expression())
    throw std::logic_error("PURRS::Recurrence::lower_bound() called, "
			   "but no lower bounds were computed");
  lower_bound_
    .set_expression(lower_bound_.replace_system_generated_symbols(*this));
  e = lower_bound_.expression();
}

inline void
Recurrence::upper_bound(Expr& e) const {
  if (!upper_bound_.has_expression())
    throw std::logic_error("PURRS::Recurrence::upper_bound() called, "
			   "but no upper bounds were computed");
  upper_bound_
    .set_expression(upper_bound_.replace_system_generated_symbols(*this));
  e = upper_bound_.expression();
}

inline Recurrence::Verify_Status
Recurrence::verify_exact_solution() const {
  if (exact_solution_.has_expression()) {
    if (is_linear_finite_order() || is_non_linear_finite_order())
      return verify_exact_solution(*this);
    else if (is_functional_equation())
      // Special case: functional equation of the form `x(n) = x(n/b)'. 
      return PROVABLY_CORRECT;
    // FIXME: temporary! At the moment there is not implemented
    // the verification in the case of infinite order recurrences.
    else { // linear infinite order
      assert(is_linear_infinite_order());
      return INCONCLUSIVE_VERIFICATION;
    }
  }
  else
    throw std::logic_error("PURRS::Recurrence::verify_exact_solution() "
			   "called, but no exact solution were computed");
}

inline Recurrence::Verify_Status
Recurrence::verify_lower_bound() const {
  if (exact_solution_.has_expression())
    return verify_exact_solution();
  else if (lower_bound_.has_expression())
    if (is_functional_equation())
      return verify_bound(*this, false);
    else
      // At the moment this case is impossible!
      return INCONCLUSIVE_VERIFICATION;
  else
    throw std::logic_error("PURRS::Recurrence::verify_lower_bound() "
			   "called, but no lower bound were computed");
}

inline Recurrence::Verify_Status
Recurrence::verify_upper_bound() const {
  if (exact_solution_.has_expression())
    return verify_exact_solution();
  else if (upper_bound_.has_expression())
    if (is_functional_equation())
      return verify_bound(*this, true);
    else
      // At the moment this case is impossible!
      return INCONCLUSIVE_VERIFICATION;
  else
    throw std::logic_error("PURRS::Recurrence::verify_upper_bound() "
			   "called, but no upper bound were computed");
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
