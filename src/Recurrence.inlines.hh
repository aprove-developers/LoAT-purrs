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

namespace Parma_Recurrence_Relation_Solver {

inline
Recurrence::Recurrence()
  : recurrence_rhs(0),
    old_recurrence_rhs(0),
    gcd_decrements_old_rhs(0),
    solution_order_reduced(0),
    recurrence_rhs_rewritten(false),
    inhomogeneous_term(0),
    type(ORDER_ZERO),
    finite_order_p(0),
    functional_eq_p(0),
    solved(false) {
}

inline
Recurrence::Recurrence(const Expr& e)
  : recurrence_rhs(e),
    old_recurrence_rhs(0),
    gcd_decrements_old_rhs(0),
    solution_order_reduced(0),
    recurrence_rhs_rewritten(false),
    inhomogeneous_term(0),
    type(UNKNOWN),
    finite_order_p(0),
    functional_eq_p(0),
    solved(false) {
}

inline
Recurrence::Recurrence(const Recurrence& y)
  : recurrence_rhs(y.recurrence_rhs),
    old_recurrence_rhs(y.old_recurrence_rhs),
    gcd_decrements_old_rhs(y.gcd_decrements_old_rhs),
    solution_order_reduced(y.solution_order_reduced),
    recurrence_rhs_rewritten(y.recurrence_rhs_rewritten),
    inhomogeneous_term(y.inhomogeneous_term),
    system_rhs(y.system_rhs),
    type(y.type),
    finite_order_p(y.finite_order_p),
    functional_eq_p(y.functional_eq_p),    
    solved(y.solved),
    solution(y.solution),
    lower_bound(y.lower_bound),
    upper_bound(y.upper_bound) {
}

inline
Recurrence::~Recurrence() {
  delete finite_order_p;
  delete functional_eq_p;
}

inline Recurrence&
Recurrence::operator=(const Recurrence& y) {
  recurrence_rhs = y.recurrence_rhs;
  old_recurrence_rhs = y.old_recurrence_rhs;
  gcd_decrements_old_rhs = y.gcd_decrements_old_rhs;
  solution_order_reduced = y.solution_order_reduced;
  recurrence_rhs_rewritten = y.recurrence_rhs_rewritten;
  inhomogeneous_term = y.inhomogeneous_term;
  system_rhs = y.system_rhs;
  type = y.type;
  finite_order_p = y.finite_order_p;
  functional_eq_p = y.functional_eq_p;
  solved = y.solved;
  solution = y.solution;
  lower_bound = y.lower_bound;
  upper_bound = y.upper_bound;
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
Recurrence::set_inhomogeneous_term(const Expr& e) const {
  inhomogeneous_term = e;
}

inline bool
Recurrence::is_unknown() const {
  return type == UNKNOWN;
}

inline bool
Recurrence::is_order_zero() const {
  return type == ORDER_ZERO; 
}

inline void
Recurrence::set_order_zero() const {
  type = ORDER_ZERO; 
}

inline bool
Recurrence::is_linear_finite_order_const_coeff() const {
  return type == LINEAR_FINITE_ORDER_CONST_COEFF;
}

inline void
Recurrence::set_linear_finite_order_const_coeff() const {
  type = LINEAR_FINITE_ORDER_CONST_COEFF;
}

inline bool
Recurrence::is_linear_finite_order_var_coeff() const {
  return type == LINEAR_FINITE_ORDER_VAR_COEFF;
}

inline void
Recurrence::set_linear_finite_order_var_coeff() const {
  type = LINEAR_FINITE_ORDER_VAR_COEFF;
}

inline bool
Recurrence::is_linear_finite_order() const {
  return (type == ORDER_ZERO
	  || type == LINEAR_FINITE_ORDER_CONST_COEFF
	  || type == LINEAR_FINITE_ORDER_VAR_COEFF);
}

inline bool
Recurrence::is_non_linear_finite_order() const {
  return type == NON_LINEAR_FINITE_ORDER;
}

inline void
Recurrence::set_non_linear_finite_order() const {
  type = NON_LINEAR_FINITE_ORDER;
}

inline bool
Recurrence::is_functional_equation() const {
  return type == FUNCTIONAL_EQUATION;
}

inline void
Recurrence::set_functional_equation() const {
  type = FUNCTIONAL_EQUATION;
}

inline unsigned int
Recurrence::order() const {
  assert(is_order_zero()
	 || is_linear_finite_order_const_coeff()
	 || is_linear_finite_order_var_coeff()
	 || is_non_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> order();
}

inline unsigned int&
Recurrence::order() {
  assert(is_order_zero()
	 || is_linear_finite_order_const_coeff()
	 || is_linear_finite_order_var_coeff()
	 || is_non_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> order();
}

inline unsigned
Recurrence::first_initial_condition() const {
  assert(is_order_zero()
	 || is_linear_finite_order_const_coeff()
	 || is_linear_finite_order_var_coeff()
	 || is_non_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> first_initial_condition();
}

inline unsigned&
Recurrence::first_initial_condition() {
  assert(is_order_zero()
	 || is_linear_finite_order_const_coeff()
	 || is_linear_finite_order_var_coeff()
	 || is_non_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> first_initial_condition();
}

inline void
Recurrence::set_first_initial_condition(unsigned i_c) const {
  assert(is_order_zero()
	 || is_linear_finite_order_const_coeff()
	 || is_linear_finite_order_var_coeff()
	 || is_non_linear_finite_order());
  assert(finite_order_p);
  finite_order_p -> set_first_initial_condition(i_c);
}

inline const std::vector<Expr>&
Recurrence::coefficients() const {
  assert(is_order_zero()
	 || is_linear_finite_order_const_coeff()
	 || is_linear_finite_order_var_coeff()
	 || is_non_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> coefficients();
}

inline std::vector<Expr>&
Recurrence::coefficients() {
  assert(is_order_zero()
	 || is_linear_finite_order_const_coeff()
	 || is_linear_finite_order_var_coeff()
	 || is_non_linear_finite_order());
  assert(finite_order_p);
  return finite_order_p -> coefficients();
}

inline Expr
Recurrence::coefficient() const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p -> coefficient();
}

inline Expr&
Recurrence::coefficient() {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p -> coefficient();
}

inline unsigned
Recurrence::divisor_arg() const {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p -> divisor_arg();
}

inline unsigned&
Recurrence::divisor_arg() {
  assert(is_functional_equation());
  assert(functional_eq_p);
  return functional_eq_p -> divisor_arg();
}

inline Recurrence::Solver_Status
Recurrence::solve() const {
  Solver_Status status = SUCCESS;
  if (!solved && (status = solve_try_hard()) == SUCCESS)
    solved = true;
  return status;
}

inline Expr
Recurrence::lower_bound_solution() const {
  if (solved || solve()) {
    assert(is_linear_finite_order() || is_functional_equation());
    if (is_linear_finite_order())
      return solution;
    else
      return lower_bound;
  }
  else
    // Well, if the client insists...
    return recurrence_rhs;
}

inline Expr
Recurrence::upper_bound_solution() const {
  if (solved || solve()) {
    assert(is_linear_finite_order() || is_functional_equation());
    if (is_linear_finite_order())
      return solution;
    else
      return upper_bound;
  }
  else
    // Well, if the client insists...
    return recurrence_rhs;
}

inline bool
Recurrence::exact_solution(Expr& exact_sol) const {
  if (solved || solve()) {
    assert(is_linear_finite_order() || is_functional_equation());
    if (is_linear_finite_order())
      exact_sol = solution;
    else
      if (lower_bound == upper_bound)
	exact_sol = lower_bound;
      else
	return false;
  }
  else
    // Well, if the client insists...
    exact_sol = recurrence_rhs;
  return true;
}

inline Expr
Recurrence::approximated_solution() const {
  if (solved || solve())
    return blackboard.approximate(solution);
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
