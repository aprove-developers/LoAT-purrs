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
#include "ep_decomp.hh"
#include "factorize.hh"
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

namespace PURRS = Parma_Recurrence_Relation_Solver;

const PURRS::Symbol&
PURRS::Recurrence::n = Symbol("n");

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

//! \brief
//! Verifies the exact solution for linear recurrences of finite order
//! and non-linear recurrences of finite order.
/*!
  Case 1: linear recurrences of finite order.
  Consider the right hand side \p rhs of the order \f$ k \f$ recurrence
  relation
  \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$.
  Let \f$ i \f$ the non-negative integer starting from which the recurrence
  is well-defined.
  The verification's process is divided in 4 steps:
  -  Validation of initial conditions.
     If <CODE>recurrence_rhs</CODE> is equal to \f$ x(i), \cdots, x(i+k) \f$
     for \f$ n = i, \cdots, i+k-1 \f$ respectively, then
     the initial conditions are verified and we can continue the verification;
     otherwise return <CODE>false</CODE> because the solution can be wrong or
     it is not enough simplified.
  -  Splits the solution in 2 parts: \p homogeneous_part contains the terms
     of the solution with the initial conditions \f$ x(i), \cdots, x(i+k) \f$,
     \p non_homogeneous_part all the other terms.
  -  By substitution, verifies that \p homogeneous_part satisfies the
     homogeneous part of the recurrence
     \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} \f$.
     Considers the difference, called \f$ d1 \f$, between
     \p homogeneous_part and the new right hand side obtained by
     substitution in the omogeneous part of the recurrence:
     - if \f$ d1 = 0 \f$     -> the verification can continue;
     - if \f$ d1 \neq 0 \f$  -> returns <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
 			        the solution can be wrong or we failed to
			        simplify it.
  -  By substitution, verifies that \p non_homogeneous_part satisfies
     the recurrence (in other words, we are considering all initial
     conditions equal to \f$ 0 \f$).
     Considers the difference, called \f$ d2 \f$, between
     \p non_homogeneous_part and the new right hand side obtained by
     substitution:
     - if \f$ d2 = 0 \f$     -> returns <CODE>PROVABLY_CORRECT</CODE>:
                                the solution is certainly right.
     - if \f$ d2 \neq 0 \f$  -> returns <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
 			        the solution can be wrong or we failed to
			        simplify it.
   FIXME: In the latter case, we will need more powerful tools to
   decide whether the solution is right or it is really wrong and, in this
   last case, to return <CODE>PROVABLY_INCORRECT</CODE>.

   Case 2: non-linear recurrences of finite order.
   Considers the order of the linear recurrence associated to that one
   non-linear, since the initial conditions are the same.
   Applies the 4 steps of the previous case.

   Case 3: linear recurrence of infinite order.
   At the moment the system solve infinite order recurrence transforming it
   in linear recurrence of first order, therefore in the solution will be
   only one initial condition.
   The verification's process is divided in 2 steps:
   -  Validation of the initial condition by substitution of \f$ n \f$
      with the value in <CODE>rec.lower_bound_sum()</CODE>, i.e., the first
      value starting from which the solution is valid.
   -  Validation of the solution using mathematic induction.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_exact_solution(const Recurrence& rec) {
  assert(rec.is_linear_finite_order() || rec.is_non_linear_finite_order()
	 || rec.is_linear_infinite_order());
  if (rec.is_linear_finite_order() || rec.is_non_linear_finite_order()) {
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
    
    if (order_rec == 0)
      return PROVABLY_CORRECT;
    else {
      // Step 1: validation of initial conditions.
      for (unsigned i = 0; i < order_rec; ++i) {
	Expr solution_evaluated
	  = rec.exact_solution_.expression().substitute(n, first_i_c + i);
	solution_evaluated = rec.blackboard.rewrite(solution_evaluated);
	solution_evaluated = simplify_all(solution_evaluated);
	if (solution_evaluated != x(first_i_c + i))
	  return INCONCLUSIVE_VERIFICATION;
      }
      
      // Step 2: split the solution in 2 parts: terms with initial conditions
      // are stored in \p homogeneous_part, all the other terms are stored in
      // \p non_homogeneous_part.
      Expr homogeneous_part = 0;
      Expr non_homogeneous_part = 0;
      if (rec.exact_solution_.expression().is_a_add())
	for (unsigned i = rec.exact_solution_.expression().nops(); i-- > 0; ) {
	  if (rec.exact_solution_.expression().op(i).has_x_function_only_ic())
	    homogeneous_part += rec.exact_solution_.expression().op(i);
	  else
	    non_homogeneous_part += rec.exact_solution_.expression().op(i);
	}
      else
	if (rec.exact_solution_.expression().has_x_function_only_ic())
	  homogeneous_part = rec.exact_solution_.expression();
	else
	  non_homogeneous_part = rec.exact_solution_.expression();
      
#if 0
      // Step 3: by substitution, verifies that `homogeneous_part'
      // satisfies the homogeneous part of the recurrence.
      // `substituted_homogeneous_rhs' is the homogeneous part of the
      // recurrence where `n' is substituted by `n - d' (where `d' is
      // the decrement of the i-th term `a_i(n) x(n - d)').
      Expr substituted_homogeneous_rhs
	= rec.recurrence_rhs - rec.inhmogeneous_term;
      // Substitutes in the homogeneous part of the recurrence the terms
      // of the form `x(n-i)'.
      for (unsigned i = 0; i < order_rec; ++i) {
	Expr shifted_solution
	  = simplify_all(homogeneous_part.substitute(n, n - (i + 1)));
	shifted_solution = simplify_sum(shifted_solution, true);
	substituted_homogeneous_rhs = substituted_homogeneous_rhs
	  .substitute(x(n - (i + 1)), shifted_solution);
      }
      Expr diff = rec.blackboard.rewrite(homogeneous_part
					 - substituted_homogeneous_rhs);
      diff = simplify_all(diff);
      if (!diff.is_zero()) {
	diff = simplify_all(diff);
	if (!diff.is_zero())
	  return INCONCLUSIVE_VERIFICATION;
      }
#endif
      
      // The recurrence is homogeneous.
      if (non_homogeneous_part == 0)
	return PROVABLY_CORRECT;
      
#if 0
      std::vector<Expr> bases_of_exp;
      std::vector<Expr> exp_poly_coeff;
      std::vector<Expr> exp_no_poly_coeff;
      exp_poly_decomposition(rec.inhomogeneous_term, Recurrence::n,
			     bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
      
      assert(bases_of_exp.size() == exp_poly_coeff.size()
	     && exp_poly_coeff.size() == exp_no_poly_coeff.size()
	     && exp_no_poly_coeff.size() >= 1);
      
      unsigned num_of_exponentials = bases_of_exp.size();
      D_VEC(bases_of_exp, 0, num_of_exponentials-1);
      D_VEC(exp_poly_coeff, 0, num_of_exponentials-1);
      D_VEC(exp_no_poly_coeff, 0, num_of_exponentials-1);
      
      unsigned max_polynomial_degree = 0;
      for (unsigned i = 0; i < num_of_exponentials; ++i) {
	if (!exp_no_poly_coeff[i].is_zero()) {
	  DD_MSGVAR("No poly: ", exp_no_poly_coeff[i]);
	  goto traditional;
	}
	max_polynomial_degree = std::max(max_polynomial_degree,
					 exp_poly_coeff[i].degree(n));
      }
      
      {
	std::vector<Polynomial_Root> roots;
	std::vector<Number> num_coefficients(order_rec + 1);
	bool all_distinct = true;
	if (rec.is_linear_finite_order_const_coeff()) {
	  Expr characteristic_eq;
	  if (!characteristic_equation_and_its_roots(order_rec,
						     rec.coefficients(),
						     num_coefficients,
						     characteristic_eq, roots,
						     all_distinct))
	    abort();
	}
	// FIXME: if this is the case of linear recurrence with variable
	// coefficient the vector of the roots is empty! 
	
	// Find the maximum degree of a polynomial that may occur in the
	// solution.
	for (unsigned i = 0, nroots = roots.size(); i < nroots; ++i) {
	  max_polynomial_degree += roots[i].multiplicity() - 1;
	  // FIXME: this may be inefficient!
	  for (unsigned j = 0; j < num_of_exponentials; ++j)
	    if (roots[i].value() == bases_of_exp[j])
	      ++max_polynomial_degree;
	}
	
	Expr substituted_rhs = rec.recurrence_rhs;
	for (unsigned i = order_rec; i-- > 0; ) {
	  Expr shifted_solution
	    = non_homogeneous_part.substitute(n, n - (i + 1));
	  //shifted_solution = simplify_sum(shifted_solution, true);
	  substituted_rhs = substituted_rhs
	    .substitute(x(n - (i + 1)), shifted_solution);
	}
	Expr diff
	  = rec.blackboard.rewrite(non_homogeneous_part - substituted_rhs);
	diff = diff.expand();
	
	std::vector<Expr> coefficients_of_exponentials(max_polynomial_degree+1);
	if (diff.is_a_add()) {
	  for (unsigned i = 0; i < diff.nops(); ++i) {
	    Expr summand = diff.op(i);
#if 0
	    if (summand.is_a_mul()) {
	      // Summand has the form `n^k * a^n * b' (with `k' possibly 0 and
	      // `a' and `b' possibly 1).
	      bool done = false;
	      for (unsigned j = 0; j < summand.nops(); ++j) {
		Expr factor = summand.op(j);
		unsigned k;
		if (factor == n)
		  k = 1;
		else if (factor.is_a_power() && factor.arg(0) == n) {
		  assert(factor.arg(1).is_a_number());
		  k = factor.arg(1).ex_to_number().to_unsigned();
		}
		else
		  continue;
		// FIXME: maybe can be implemented in a better way.
		// FIXME: k uninitialized here!
		coefficients_of_exponentials[k] += summand/factor;
		done = true;
		break;
	      }
	      if (!done)
		coefficients_of_exponentials[0] += summand;
	    }
#else
	    if (summand.is_a_mul()) {
	      // Summand has the form `n^k * a^n * b' (with `k' possibly 0 and
	      // `a' and `b' possibly 1).
	      bool done = false;
	      for (unsigned j = 0; (j < summand.nops()) && !done; ++j) {
		Expr factor = summand.op(j);
		if (factor == n) {
		  coefficients_of_exponentials[1] += summand/factor;
		  done = true;
		}
		else if (factor.is_a_power() && factor.arg(0) == n) {
		  assert(factor.arg(1).is_a_number());
		  unsigned k = factor.arg(1).ex_to_number().to_unsigned();
		  assert(k < coefficients_of_exponentials.size());
		  coefficients_of_exponentials[k] += summand/factor;
		  done = true;
		}
	      }
	      // `done' is false if `factor' contains neither `n' nor `n^k'
	      if (!done)
		coefficients_of_exponentials[0] += summand;
	    }
#endif
	    else if (summand.is_a_power()) {
	      // Summand has the form `n^k' or `a^n'
	      if (summand.arg(0) == n) {
		unsigned k = summand.arg(1).ex_to_number().to_unsigned();
		coefficients_of_exponentials[k] += 1;
	      }
	      else
		coefficients_of_exponentials[0] += summand;
	    }
	    else {
	      // Summand is a constant or the symbol `n'.
	      assert(summand.is_a_number() || summand == n);
	      Number k;
	      if (summand.is_a_number(k))
		coefficients_of_exponentials[0] += k;
	      else
		coefficients_of_exponentials[1] += 1;
	    }
	  }
	}
	else {
	  Expr summand = diff;
	  if (summand.is_a_mul()) {
	    // Summand has the form `n^k * a^n * b' (with `k' possibly 0 and
	    // `a' and `b' possibly 1).
	    bool done = false;
	    for (unsigned j = 0; j < summand.nops(); ++j) {
	      Expr factor = summand.op(j);
	      unsigned k;
	      if (factor == n)
		k = 1;
	      else if (factor.is_a_power() && factor.arg(0) == n) {
		assert(factor.arg(1).is_a_number());
		k = factor.arg(1).ex_to_number().to_unsigned();
	      }
	      else
		continue;
	      // FIXME: maybe can be implemented in a better way.
	      coefficients_of_exponentials[k] += summand/factor;
	      done = true;
	      break;
	    }
	    if (!done)
	      coefficients_of_exponentials[0] += summand;
	  }      
	  else if (summand.is_a_power()) {
	    // Summand has the form `n^k' or `a^n'
	    if (summand.arg(0) == n) {
	      unsigned k = summand.arg(1).ex_to_number().to_unsigned();
	      coefficients_of_exponentials[k] += 1;
	    }
	    else
	      coefficients_of_exponentials[0] += summand;
	  }
	  else {
	    // Summand is a constant or the symbol `n'.
	    assert(summand.is_a_number() || summand == n);
	    Number k;
	    if (summand.is_a_number(k))
	      coefficients_of_exponentials[0] += k;
	    else
	      coefficients_of_exponentials[1] += 1;
	  }
	}
	
	D_VEC(coefficients_of_exponentials, 0, max_polynomial_degree);
	
	Number num_tests = num_of_exponentials + order_rec;
	for (unsigned i = 0; i < max_polynomial_degree; ++i) {
	  if (!coefficients_of_exponentials[i].is_zero()) {
	    // Not syntactically 0: try to prove that is it semantically 0.
	    Expr c = coefficients_of_exponentials[i];
	    for (Number r = 0; r < num_tests; ++r) {
	      Expr c_r = simplify_all(c.substitute(n, r));
	      if (!c_r.is_zero()) {
		DD_MSGVAR("Argh!!! ", i);
		DD_MSGVAR("Argh!!! ", r);
		DD_MSGVAR("Argh!!! ", c_r);
		goto traditional;
	      }
	    }
	  }
	}
	return PROVABLY_CORRECT;
      }
      
      
    traditional:
#endif
      // Step 4: by substitution, verifies that `non_homogeneous_part'
      // satisfies the recurrence.
      // Computes `substituted_rhs' by substituting, in the rhs
      // of the recurrence, `n' by `n - d' (where `d' is the decrement
      // of the i-th term `a(n) x(n - d)').
      Expr substituted_rhs = rec.recurrence_rhs;
      for (unsigned i = 0; i < order_rec; ++i) {
	Expr shifted_solution
	  = simplify_all(non_homogeneous_part.substitute(n, n - (i + 1)));
	shifted_solution = simplify_sum(shifted_solution, true);
	substituted_rhs = substituted_rhs
	  .substitute(x(n - (i + 1)), shifted_solution);
      }
      Expr diff = rec.blackboard.rewrite(non_homogeneous_part-substituted_rhs);
      diff = simplify_all(diff);
      if (!diff.is_zero())
	if (rec.applied_order_reduction()) {
	  rec.unset_order_reduction();
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
	    (write_reduced_order_recurrence(rec.recurrence_rhs, r,
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
	else {
	  diff = simplify_all(diff);
	  if (!diff.is_zero())
	    return INCONCLUSIVE_VERIFICATION;
	}
      return PROVABLY_CORRECT;
    }
  }
  else {
    assert(rec.is_linear_infinite_order());
    // The case `f(n) = -1', i.e. the recurrence has the form
    // `x(n) = - sum(k, 0, n-1, x(k)) + g(n)', is special:
    // the solution is simply `x(n) = g(n) - g(n-1)'.
    // FIXME: the traditional validation' process does not work,
    // it is true?
    if (rec.weight_inf_order() == -1)
      return PROVABLY_CORRECT;

    Expr weight;
    Expr inhomogeneous;
    if (rec.upper_bound_sum() == n-1) {
      weight = rec.weight_inf_order();
      inhomogeneous = rec.inhomogeneous_term;
    }
    else {
      // In the case of `upper_bound_sum() == n' the recurrence is
      // transformed in another equivalent recurrence with the upper
      // bound of the sum equal to `n-1': we verify the solution
      // considering this last recurrence.
      assert(rec.upper_bound_sum() == n);
      weight =  rec.weight_inf_order() / (1 - rec.weight_inf_order());
      inhomogeneous = rec.inhomogeneous_term / (1 - rec.weight_inf_order());
    }

    // Note: the solution is valid only for `n > lower_bound_sum()'.

    // Step 1: validation of the initial condition.
    Expr solution_evaluated = rec.exact_solution_.expression()
      .substitute(n, rec.lower_bound_sum()+1);
    solution_evaluated = simplify_all(solution_evaluated);
    if (solution_evaluated
	!= weight.substitute(n, rec.lower_bound_sum()+1)
	* x(rec.lower_bound_sum())
	+ inhomogeneous.substitute(n, rec.lower_bound_sum()+1))
      return INCONCLUSIVE_VERIFICATION;
    
    // Step 2: validation of the solution using mathematic induction.
    // Note: with recurrences of infinite order is not true that if it is
    // homogeneous then is sufficient to verify the initial condition.
    // In fact, the recurrences could be in a non-standard form like
    // `x(n) = f(n) sum(k, 0, n-1, x(k)+h(n))'. For standard form
    // we consider recurrences with the argument of the sum equal to `x(k)'. 

    // Let `x(n) = f(n) sum(k, n_0, n-1, x(k)) + g(n)' be the infinite
    // order recurrence. Now we must consider the expression
    // `x(n) - f(n) x(n_0) - f(n) sum(k, n_0 + 1, n - 1, x(k)) - g(n)',
    // where `x(n)' and `x(k)' are substituted with the solution.
    // If we are able to demonstrate that this expression is syntactically
    // zero, then the solution is correct.
    Symbol h("h");
    Expr diff = PURRS::sum(h, rec.lower_bound_sum() + 1, n - 1,
			   rec.exact_solution_.expression()
			   .substitute(n, h));
    diff = simplify_sum(diff, false, true);
    diff = rec.exact_solution_.expression()
      - diff * weight - x(rec.lower_bound_sum()) * weight - inhomogeneous;
    diff = simplify_all(diff);
    if (diff == 0)
      return PROVABLY_CORRECT;
    else
      return INCONCLUSIVE_VERIFICATION;
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
      if (!bound.op(i).has_x_function_only_ic())
	partial_bound += bound.op(i);
    }
  else
    if (!bound.has_x_function_only_ic())
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
  set_order_reduction();
  // Build the new recurrence substituting `n' not contained in the
  // `x' functions with `gcd_among_decrements * n + r' and `x(n-k)' with
  // `x(n - k / gcd_among_decrements)'.
  Symbol r = insert_auxiliary_definition(mod(n, gcd_among_decrements()));
  unsigned dim = coefficients().size() / gcd_among_decrements() + 1;
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
    bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
    rec_rewritten = true;
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
	bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	rec_rewritten = true;
	return SUCCESS;
      }
      else
	return status;

    // Functional equation.
    else {
      if (type == 1)
	if ((status = rec_rewritten.approximate_functional_equation_lower())
	    == SUCCESS)
	  solution_or_bound = pwr(base_exp_log(),
				  rec_rewritten.lower_bound_.expression());
	else
	  return status;
      else if (type == 2)
	if ((status = rec_rewritten.approximate_functional_equation_upper())
	    == SUCCESS)
	  solution_or_bound = pwr(base_exp_log(),
				  rec_rewritten.upper_bound_.expression());
	else
	  return status;
      else {// type == 0
	if ((status = rec_rewritten.approximate_functional_equation_lower())
	    != SUCCESS
	    || (status = rec_rewritten.approximate_functional_equation_upper())
	    != SUCCESS)
	  return status;
	if (rec_rewritten.lower_bound_.expression()
	    != rec_rewritten.upper_bound_.expression())
	  return TOO_COMPLEX;
      }
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
      bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
      rec_rewritten = true;
      return SUCCESS;
    }
  }
  else
    return status;
}

/*!
  Builds the recurrence of infinite order
  \f$ x(n) = f(n) \sum_{k=n_0}^{n-1} x(k) + g(n) \f$, where
  \f$ f(n) \f$ is stored in \p weight; \f$ g(n) \f$ is stored
  in \p inhomogeneous and \p first_well_defined contains the smallest
  positive integer starting from which the recurrence is well-defined.
  If the system is able to solve the recurrence, then this function
  returns <CODE>SUCCESS</CODE> and the solution is stored in \p solution.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::solve_new_infinite_order_rec(const Expr& weight,
						const Expr& inhomogeneous,
						unsigned first_well_defined,
						Expr& solution) const {
  Symbol h;
  Recurrence rec_rewritten(weight * PURRS::sum(h, 0, n-1, x(h))
			   + inhomogeneous);
  const Expr& coeff_first_order
    = weight  / weight.substitute(n, n-1) * (1 + weight.substitute(n, n-1));
  const Expr& inhomog_first_order = weight
    * (inhomogeneous / weight - (inhomogeneous / weight).substitute(n, n-1));
  rec_rewritten.infinite_order_p
    = new Infinite_Order_Info(coeff_first_order*x(n-1)+inhomog_first_order,
			      coeff_first_order, inhomog_first_order,
			      weight, 0, n-1);
  rec_rewritten.set_linear_infinite_order();
  rec_rewritten.set_infinite_order_fwdr(first_well_defined);
  rec_rewritten.set_inhomogeneous_term(inhomogeneous);
  return rec_rewritten.compute_infinite_order_recurrence(solution);
}

/*!
  Builds a new object <CODE>Recurrence</CODE> containing a linear
  recurrence of finite order built starting from the recurrence
  \p *this of infinite order.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::
compute_infinite_order_recurrence(Expr& solution) const {
  // At the moment we consider only recurrences of infinite order
  // transformable in first order recurrences.
  if (weight_inf_order() == -1) {
    // Special case: the recurrence of infinite order starting
    // from which we have found the components of the linear recurrence
    // of finite order was of the form `x(n) = f(n) sum(k,0,n-1,x(k)) + g(n)'
    // with `f(n) = -1'. In this case the solution of the recurrence
    // is simply `x(n) = g(n) - g(n-1)' and is not necessary all the normal
    // procedure.
    solution
      = inhomogeneous_term - inhomogeneous_term.substitute(n, n-1).expand();
    return SUCCESS;
  }
  else {
    unsigned lower = lower_bound_sum();
    const Expr& upper = upper_bound_sum();
    if (lower == 0 && upper == n-1) {
      std::vector<Expr> coefficients(2);
      // Shift forward `n -> n + 1' the coefficient.
      coefficients[1] = coeff_first_order().substitute(n, n+1);
      Recurrence rec_rewritten(rhs_transformed_in_first_order());
      // It is not necessary to repeat the classification's process
      // because we already know the order, the coefficients and the
      // inhomogeneous term.
      rec_rewritten.finite_order_p = new Finite_Order_Info(1, coefficients, 1);
      if (coeff_first_order().has(n))
	rec_rewritten.set_type(LINEAR_FINITE_ORDER_VAR_COEFF);
      else
	rec_rewritten.set_type(LINEAR_FINITE_ORDER_CONST_COEFF);
      // Shift forward `n -> n + 1' the inhomogeneous term.
      rec_rewritten.set_inhomogeneous_term(inhomog_first_order()
					   .substitute(n, n+1));
      Solver_Status status;
      if ((status = rec_rewritten.solve_linear_finite_order())
	  == SUCCESS) {
	solution = rec_rewritten.exact_solution_.expression();
	// Shift backward: n -> n - 1.
	solution = solution.substitute(n, n - 1);
	solution = solution
	  .substitute(x(rec_rewritten.first_well_defined_rhs_linear()),
		      x(rec_rewritten.first_well_defined_rhs_linear()+1));
	// If there is the initial condition `x(1)' specified then
	// the system substitute it with the respective value; otherwise
	// the system does the substitution `x(1) = 2*x(0)+1'.
	// FIXME: At the moment we substitute here only the initial
	// condition `x(1)'.
	std::map<unsigned, Expr>::const_iterator i
	  = initial_conditions.find(1);
	if (i != initial_conditions.end())
	  solution = solution.substitute(x(1), get_initial_condition(1));
	else
	  solution = solution
	    .substitute(x(rec_rewritten.first_well_defined_rhs_linear()+1),
			(weight_inf_order()*x(0)+inhomogeneous_term)
			.substitute(n, 1));
	//	solution = simplify_ex_for_output(solution, false);
	return SUCCESS;
      }
      else
	return status;
    }
    else
      if (upper == n)
	// Transform the recurrence `x(n) = f(n) sum(k, n_0, n, x(k)) + g(n)'
	// in the equivalent recurrence
	// `x(n) = f(n) / (1-f(n)) sum(k, n_0, n-1, x(k)) + g(n) / (1-f(n))'.
	if (weight_inf_order() != -1) {
	  Expr tmp = 1 - weight_inf_order();
	  Expr weight_rewritten = weight_inf_order() / tmp;
	  Expr inhomogeneous_rewritten = inhomogeneous_term / tmp;
	  // Consider the "shifted" recurrence
	  // `y(n) = f(n+n_0) / (1-f(n+n_0)) sum(k, 0, n-1, x(k))
	  //           + g(n+n_0) / (1-f(n+n_0))'; later on, the solution
	  // `y(n)' will be shifted in order to find `x(n)'.
	  if (lower > 0) {
	    weight_rewritten = weight_rewritten.substitute(n, n + lower);
	    inhomogeneous_rewritten
	      = inhomogeneous_rewritten.substitute(n, n + lower);
	    tmp = tmp.substitute(n, n + lower);
	  }
	  Solver_Status status;
	  if ((status = solve_new_infinite_order_rec(weight_rewritten,
						     inhomogeneous_rewritten,
						     infinite_order_fwdr(),
						     solution))
	      != SUCCESS)
	    return status;
	  // We must shift the solution: n    -> n - lower,
	  //                             x(a) -> x(a + lower).
	  if (lower > 0) {
	    solution = solution.substitute(n, n - lower);
	    solution = solution.substitute(x(infinite_order_fwdr()),
					   x(infinite_order_fwdr() + lower));
	  }
	  return SUCCESS;
	}
	else
	  // FIXME: how we can do?
	  return TOO_COMPLEX;
    // Consider the "shifted" recurrence
    // `y(n) = f(n+n_0) sum(k, 0, n-1, x(k)) + g(n+n_0)';
    // later on, the solution `y(n)' will be shifted in order to find `x(n)'.
    // Note: if `infinite_order_fwdr() > 0' means that the recurrence does
    // not have any sense for positive integer less than
    // `infinite_order_fwdr()', so it is not possible to consider the
    // aforesaid method.
      else if (lower > 0 && infinite_order_fwdr() == 0) {
	const Expr& weight_rewritten
	  = weight_inf_order().substitute(n, n + lower);
	const Expr& inhomogeneous_rewritten
	  = inhomogeneous_term.substitute(n, n + lower);
	Solver_Status status;
	if ((status = solve_new_infinite_order_rec(weight_rewritten,
						   inhomogeneous_rewritten,
						   infinite_order_fwdr(),
						   solution))
	    != SUCCESS)
	  return status;
	// We must shift the solution: n    -> n - lower,
	//                             x(a) -> x(a + lower).
	solution = solution.substitute(n, n - lower);
	solution = solution.substitute(x(infinite_order_fwdr()),
				       x(infinite_order_fwdr() + lower));
	return SUCCESS;
      }
      else
	return TOO_COMPLEX;
  }
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
  recurrence is well-defined. This function checks if in the map
  \p initial_conditions there are some initial conditions
  of the form \f$ x(i) = k \f$ with \f$ k > s \f$: in this case
  the function shifts the solution or the bound.
  Finally substitutes to the arbitrary initial conditions in the solution or
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
    
    // Shift initial conditions.
    if (first_well_defined_rhs + order_or_rank - 1 < max_i_c) {
      unsigned shift_forward = max_i_c - first_well_defined_rhs;
      for (unsigned i = order_or_rank; i-- > 0; )
	sol_or_bound
	  = sol_or_bound.substitute(x(i + first_well_defined_rhs),
				    x(i + first_well_defined_rhs
				      + shift_forward - order_or_rank + 1));
      // The solution of `x(n) = a(n) x(n-1) + p(n)' is of the form
      // `x(n) = prod(k,i+1,n,a(k)) x(i)
      //         + prod(k,i+1,n,a(k)) sum(k,i,n,p(k)/a!(k))', where
      // `i' is the smallest positive integer starting from which the
      // recurrence is well-defined. If `max_i_c', `m' for short, is bigger
      // than `i', then the solution is:
      // `x(n) = prod(k,i+1,n,a(k)) / prod(k,i+1,m,a(k)) x(i)
      //         +prod(k,i+1,n,a(k))*[sum(k,i+1,n,p(k)/prod(j,i+1,k,a(j)))
      //                              - sum(k,i+1,m,p(k)/prod(j,i+1,k,a(j)))]'.
      if (is_linear_finite_order_var_coeff()) {
	const Expr& homogeneous_term = product_factor()
	  * x(first_well_defined_rhs + shift_forward);
	const Expr& non_homogeneous_term = sol_or_bound - homogeneous_term;
	Symbol index;
	sol_or_bound = homogeneous_term
	  / PURRS::prod(index, first_well_defined_rhs+1, max_i_c,
			coefficients()[1].substitute(n, index)).ex_to_number()
	  + non_homogeneous_term - product_factor()
	  * PURRS::sum(index, first_well_defined_rhs + 1, max_i_c,
		       (inhomogeneous_term / product_factor())
		       .substitute(n, index));
      }
      else
	// Shift the index of the recurrence `n'.
	// FIXME: this technique is surely valid in the case of
	// linear recurrences with constant coefficients of finite order.
	// To check if it is valid also in the case of functional equations
	// and infinite order recurrences.
	sol_or_bound
	  = sol_or_bound.substitute(n,
				    n - (shift_forward - order_or_rank + 1));
    }
  }

  // Substitute initial conditions with the values in the map
  // `initial_conditions'.
  for (std::map<unsigned, Expr>::const_iterator i = initial_conditions.begin(),
	 iend = initial_conditions.end(); i != iend; ++i)
    sol_or_bound = sol_or_bound.substitute(x(i->first),
					   get_initial_condition(i->first));

  sol_or_bound = simplify_numer_denom(sol_or_bound);
  sol_or_bound = simplify_ex_for_output(sol_or_bound, false);
  return sol_or_bound;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution() const {
  bool& tried_to_compute_exact_solution
    = const_cast<bool&>(tested_exact_solution);
  tried_to_compute_exact_solution = true;
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
      if ((status = approximate_functional_equation_lower()) == SUCCESS
	  && (status = approximate_functional_equation_upper()) == SUCCESS
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
      Expr solution;
      if ((status = compute_infinite_order_recurrence(solution))
	  == SUCCESS) {
	// FIXME: At the moment we substitute here only the initial
	// condition `x(0)'.
	std::map<unsigned, Expr>::const_iterator i
	  = initial_conditions.find(0);
	if (i != initial_conditions.end())
	  solution = solution.substitute(x(0), get_initial_condition(0));
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
      if (!tested_exact_solution)
	// There is an exact solution.
	return compute_exact_solution();
      else
	return TOO_COMPLEX;
    // Functional equation.
    else if (is_functional_equation()) {
      if ((status = approximate_functional_equation_lower()) != SUCCESS)
	return status;
      if (!initial_conditions.empty())
	lower_bound_.set_expression
	  (substitute_i_c_shifting(false, lower_bound_.expression()));
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
	return compute_exact_solution();
      else
	return TOO_COMPLEX;
    // Functional equation.
    else if (is_functional_equation()) {
      if ((status = approximate_functional_equation_upper()) != SUCCESS)
	return status;
      if (!initial_conditions.empty())
	upper_bound_.set_expression
	  (substitute_i_c_shifting(false, upper_bound_.expression()));
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
  s << "recurrence_rewritten = "
    << (recurrence_rewritten ? "true" : "false") << std::endl;
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
