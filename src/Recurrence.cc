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

namespace PURRS = Parma_Recurrence_Relation_Solver;

const PURRS::Symbol&
PURRS::Recurrence::n = Symbol("n");

namespace {
using namespace PURRS;

void
split(const Expr& e, const Expr& d, Expr& term_with_d, Expr& other_terms) {
  assert(e.is_a_add());
  for (unsigned int i = e.nops(); i-- > 0; ) {
    const Expr& term = e.op(i);
    if (term.has(d))
      term_with_d += term;
    else
      other_terms += term;
  }
}

bool
ok_inequalities(const Expr& e, unsigned int condition) {
  assert(e.is_a_add());
  Expr term_with_n = 0;
  Expr other_terms = 0;
  split(e, Recurrence::n, term_with_n, other_terms);
  D_VAR(other_terms);
  if (term_with_n == Recurrence::n || term_with_n.is_a_mul()) {
    Expr coeff_n = 1;
    if (term_with_n.is_a_mul()) {
      for (unsigned int i = term_with_n.nops(); i-- > 0; ) {
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

} // anonymous namespace


PURRS::Recurrence::Verify_Status
PURRS::Recurrence::
validate_initial_conditions(index_type order,
			    index_type first_valid_index) const {
  for (index_type i = 0; i < order; ++i) {
    Expr e = exact_solution_.expression().substitute(n, first_valid_index + i);
    // Expand blackboard's definitions in order to increase the
    // opportunities for simplifications.
    e = blackboard.rewrite(e);
    e = simplify_all(e);
    // If in the map `initial_conditions' there is an expression
    // `e' correspondent to the integer `first_valid_index + i'
    // (the index of the initial condition), then returns `e';
    // returns `x(first_valid_index + i)' otherwise.
    if (e != get_initial_condition(first_valid_index + i))
      // FIXME: pravably_incorrect nei casi semplici.
      return INCONCLUSIVE_VERIFICATION;
  }
  return PROVABLY_CORRECT;
}

//! \brief
//! Verify the exact solution of recurrence \p *this, where the recurrence
//! can be linear of finite order or non-linear of finite order.
/*!
  Case 1: linear recurrences of finite order.
  Consider the right hand side \p rhs of the order \f$ k \f$ recurrence
  relation
  \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$,
  which is stored in the expression \p recurrence_rhs.
  Let \f$ i \f$ be \p first_valid_index, that is, the least non-negative
  integer \f$ j \f$ such that the recurrence is well-defined for
  \f$ n \geq j \f$.
  Assume that the system has produced the expression
  \p exact_solution_.expression(), which is a function of the
  variable \f$ n \f$.

  The verification's process is divided in 4 steps:
  -  Validation of \ref initial_conditions "symbolic initial conditions".
     Evaluate the expression \p exact_solution_.expression() for
     \f$ n = i, \cdots, i+k-1 \f$ and simplify as much as possible.
     If the final result is <EM>synctactically</EM> equal to \f$ x(n) \f$
     for \f$ n = i, \cdots, i+k-1 \f$, the symbolic initial conditions are
     verified and we can proceed to step 2; otherwise return
     <CODE>INCONCLUSIVE_VERIFICATION</CODE> because the
     solution can be wrong or it is not enough simplified.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
  -  Split \p exact_solution_.expression() in 2 expressions:
     \p summands_with_i_c contains the summands with an occurrence of a
     symbolic initial conditions \f$ x(i), \cdots, x(i+k) \f$;
     \p summands_without_i_c contains all the other summands.
  -  Verify that \p summands_with_i_c satisfies the homogeneous part of
     the recurrence
     \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} \f$.
     Replace \f$ x(n-i) \f$ by \p exact_solution_.expression()
     evaluated at \f$ n-i \f$ (for \f$ i = 1, \cdots, k \f$) in the
     above expression, and store the result in
     \p substituted_homogeneous_rhs.
     Consider the difference
     \f$ d1 = summands_with_i_c - substituted_homogeneous_rhs \f$:
     - if \f$ d1 = 0 \f$     -> proceed to step 3;
     - if \f$ d1 \neq 0 \f$  -> returns <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
 			        the solution can be wrong or we failed to
			        simplify it.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
  -  Verify that \p summands_without_i_c satisfies the recurrence
     (in other words, we are considering all initial conditions equal
     to \f$ 0 \f$).
     Replace \f$ x(n-i) \f$ by \p exact_solution_.expression()
     evaluated at \f$ n-i \f$ (for \f$ i = 1, \cdots, k \f$) in the
     right hand side of the recurrence, and store the result in
     \p substituted_rhs.
     Consider the difference
     \f$ d2 = summands_without_i_c - substituted_rhs \f$:
     - if \f$ d2 = 0 \f$     -> returns <CODE>PROVABLY_CORRECT</CODE>:
                                the solution is certainly right.
     - if \f$ d2 \neq 0 \f$  -> returns <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
 			        the solution can be wrong or we failed to
			        simplify it.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.

   FIXME: In the latter case, we will need more powerful tools to
   decide whether the solution is right or it is really wrong and, in this
   last case, to return <CODE>PROVABLY_INCORRECT</CODE>.
     
   Case 2: non-linear recurrences of finite order.
   Considers the order of the linear recurrence associated to that one
   non-linear, since the initial conditions are the same.
   Applies the 4 steps of the previous case.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_finite_order() const {
  // We will store here the order of the recurrence.
  index_type order_rec;
  // We will store here the the least non-negative integer `j'
  // such that the recurrence is well-defined for `n >= j':
  // the index `k' of the symbolic initial condition `x(k)'
  // will start from it.
  index_type first_i_c;
  // FIXME: chiarire discorso dell'ordine con le non lineari.
  if (is_non_linear_finite_order()) {
    order_rec = associated_linear_rec().order();
    first_i_c = associated_linear_rec().first_valid_index();
  }
  else {
    order_rec = order();
    first_i_c = first_valid_index();
  }
  
  if (order_rec == 0)
    // We call recurrence of `order zero' special recurrences of the
    // form `x(n) = rhs', where `rhs' contains only functions of `n',
    // parameters and `x(k_1),...,x(k_m)' where `m >= 0' and
    // `k_1,...,k_m' are non-negative integers.
    // In these cases are not mathematical relationship expressing `x(n)'
    // as some combination of `x(i)', with i < n and the solution
    // is simply `rhs'.
    return PROVABLY_CORRECT;
  
  // Step 1: validation of symbolic initial conditions.
  Verify_Status status = validate_initial_conditions(order_rec, first_i_c);
  if (status == INCONCLUSIVE_VERIFICATION || status == PROVABLY_INCORRECT)
    return status;
  
  // Step 2: split the solution in 2 parts: terms with initial conditions
  // are stored in `summands_with_i_c', all the other terms are stored in
  // `summands_without_i_c'.
  Expr summands_with_i_c = 0;
  Expr summands_without_i_c = 0;
  const Expr& exact_solution = exact_solution_.expression();
  if (exact_solution.is_a_add())
    for (unsigned int i = exact_solution.nops(); i-- > 0; ) {
      const Expr& addend_exact_solution = exact_solution.op(i);
      if (has_symbolic_initial_conditions(addend_exact_solution))
	summands_with_i_c += addend_exact_solution;
      else
	summands_without_i_c += addend_exact_solution;
    }
  else
    if (has_symbolic_initial_conditions(exact_solution))
      summands_with_i_c = exact_solution;
    else
      summands_without_i_c = exact_solution;
  
#if 0
  // Step 3: by substitution, verifies that `summands_with_i_c'
  // satisfies the homogeneous part of the recurrence.
  // `substituted_homogeneous_rhs' is the homogeneous part of the
  // recurrence where `n' is substituted by `n - d' (where `d' is
  // the decrement of the i-th term `a_i(n) x(n - d)').
  Expr substituted_homogeneous_rhs = recurrence_rhs - inhmogeneous_term;
  // Substitutes in the homogeneous part of the recurrence the terms
  // of the form `x(n-i)'.
  for (index_type i = 0; i < order_rec; ++i) {
    Expr shifted_solution
      = simplify_all(summands_with_i_c.substitute(n, n - (i + 1)));
    shifted_solution = simplify_sum(shifted_solution, REWRITE_UPPER_LIMIT);
    substituted_homogeneous_rhs
      = substituted_homogeneous_rhs
      .substitute(x(n - (i + 1)), shifted_solution);
  }
  Expr diff
    = blackboard.rewrite(summands_with_i_c - substituted_homogeneous_rhs);
  diff = simplify_all(diff);
  if (!diff.is_zero()) {
    diff = simplify_all(diff);
    if (!diff.is_zero())
      return INCONCLUSIVE_VERIFICATION;
  }
#endif
  
  // The recurrence is homogeneous.
  if (summands_without_i_c == 0)
    return PROVABLY_CORRECT;
  
#if 0
  std::vector<Expr> bases_of_exp;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous_term.expand(), Recurrence::n,
			 bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
  
  assert(bases_of_exp.size() == exp_poly_coeff.size()
	 && exp_poly_coeff.size() == exp_no_poly_coeff.size()
	 && exp_no_poly_coeff.size() >= 1);
  
  unsigned int num_of_exponentials = bases_of_exp.size();
  D_VEC(bases_of_exp, 0, num_of_exponentials-1);
  D_VEC(exp_poly_coeff, 0, num_of_exponentials-1);
  D_VEC(exp_no_poly_coeff, 0, num_of_exponentials-1);
  
  unsigned int max_polynomial_degree = 0;
  for (unsigned int i = 0; i < num_of_exponentials; ++i) {
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
    if (is_linear_finite_order_const_coeff()) {
      Expr characteristic_eq;
      if (!characteristic_equation_and_its_roots(order_rec,
						 coefficients(),
						 num_coefficients,
						 characteristic_eq, roots,
						 all_distinct))
	abort();
    }
    // FIXME: if this is the case of linear recurrence with variable
    // coefficient the vector of the roots is empty! 
    
    // Find the maximum degree of a polynomial that may occur in the
    // solution.
    for (unsigned int i = 0, nroots = roots.size(); i < nroots; ++i) {
      max_polynomial_degree += roots[i].multiplicity() - 1;
      // FIXME: this may be inefficient!
      for (unsigned int j = 0; j < num_of_exponentials; ++j)
	if (roots[i].value() == bases_of_exp[j])
	  ++max_polynomial_degree;
    }
    
    Expr substituted_rhs = recurrence_rhs;
    for (index_type i = order_rec; i-- > 0; ) {
      Expr shifted_solution = summands_without_i_c.substitute(n, n - (i + 1));
      //shifted_solution = simplify_sum(shifted_solution, REWRITE_UPPER_LIMIT);
      substituted_rhs = substituted_rhs
	.substitute(x(n - (i + 1)), shifted_solution);
    }
    Expr diff = blackboard.rewrite(summands_without_i_c - substituted_rhs);
    diff = diff.expand();
    
    std::vector<Expr> coefficients_of_exponentials(max_polynomial_degree+1);
    if (diff.is_a_add()) {
      for (unsigned int i = 0; i < diff.nops(); ++i) {
	Expr summand = diff.op(i);
#if 0
	if (summand.is_a_mul()) {
	  // Summand has the form `n^k * a^n * b' (with `k' possibly 0 and
	  // `a' and `b' possibly 1).
	  bool done = false;
	  for (unsigned int j = 0; j < summand.nops(); ++j) {
	    Expr factor = summand.op(j);
	    unsigned int k;
	    if (factor == n)
	      k = 1;
	    else if (factor.is_a_power() && factor.arg(0) == n) {
	      assert(factor.arg(1).is_a_number());
	      k = factor.arg(1).ex_to_number().to_unsigned int();
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
	  for (unsigned int j = 0; (j < summand.nops()) && !done; ++j) {
	    Expr factor = summand.op(j);
	    if (factor == n) {
	      coefficients_of_exponentials[1] += summand/factor;
	      done = true;
	    }
	    else if (factor.is_a_power() && factor.arg(0) == n) {
	      assert(factor.arg(1).is_a_number());
	      unsigned int k = factor.arg(1).ex_to_number().to_unsigned int();
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
	    unsigned int k = summand.arg(1).ex_to_number().to_unsigned_int();
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
	for (unsigned int j = 0; j < summand.nops(); ++j) {
	  Expr factor = summand.op(j);
	  unsigned int k;
	  if (factor == n)
	    k = 1;
	  else if (factor.is_a_power() && factor.arg(0) == n) {
	    assert(factor.arg(1).is_a_number());
	    k = factor.arg(1).ex_to_number().to_unsigned_int();
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
	  unsigned int k = summand.arg(1).ex_to_number().to_unsigned_int();
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
    for (unsigned int i = 0; i < max_polynomial_degree; ++i) {
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

  // Step 4: by substitution, verifies that `summands_without_i_c'
  // satisfies the recurrence.
  // Computes `substituted_rhs' by substituting, in the rhs
  // of the recurrence, `n' by `n - d' (where `d' is the decrement
  // of the i-th term `a(n)*x(n - d)').
  Expr substituted_rhs = recurrence_rhs;
  for (index_type d = 1; d <= order_rec; ++d) {
    Expr shifted_solution
      = simplify_all(summands_without_i_c.substitute(n, n - d));
    shifted_solution = simplify_sum(shifted_solution, REWRITE_UPPER_LIMIT);
    substituted_rhs = substituted_rhs.substitute(x(n - d), shifted_solution);
  }
  Expr diff = summands_without_i_c - substituted_rhs;
  // Differently from the step 1 (validation of symbolic initial condition)
  // the expression `diff' now contains `n' and is more difficult
  // to simplify it. For this motive we performed simplification also
  // before to expand blackboard's definitions.
  diff = blackboard.rewrite(diff);
  diff = simplify_all(diff);
  if (diff.is_zero())
    return PROVABLY_CORRECT;
  
  // If we have applied the order reduction and we did not succeed
  // in the verification of the original recurrence, then
  // we try to verify the solution of the reduced recurrence.
  // Note that "verify the original recurrence" is the first tentative
  // because the solution given to the user is relative to the non-reduced
  // recurrence. The verification of the reduced recurrence, tried only
  // if the previous verification fails, does not guarantee the correctness
  // of the solution supplied to the user in the case of errors in the
  // computation of the solution of the original recurrence starting
  // from the solution of the reduced recurrence.
  if (applied_order_reduction()) {
    unset_order_reduction();
    Symbol r = insert_auxiliary_definition(mod(n, gcd_among_decrements()));
    unsigned int dim = coefficients().size() / gcd_among_decrements() + 1;
    std::vector<Expr> new_coefficients(dim);
    Expr inhomogeneous = 0;
    Recurrence rec_rewritten
      (write_reduced_order_recurrence(recurrence_rhs, r,
				      gcd_among_decrements(),
				      coefficients(), new_coefficients,
				      inhomogeneous));
    rec_rewritten.finite_order_p
      = new Finite_Order_Info(dim - 1, new_coefficients, 1);
    rec_rewritten.set_type(type());
    rec_rewritten.set_inhomogeneous_term(inhomogeneous);
    rec_rewritten.solve_linear_finite_order();
    return rec_rewritten.verify_exact_solution();
  }
  else {
    diff = simplify_all(diff);
    if (diff.is_zero())
      return PROVABLY_CORRECT;
    else
      return INCONCLUSIVE_VERIFICATION;
  }
}

//! Verify the exact solution of the weighted-average recurrence \p *this.
/*!
  Consider the right hand side of a weighted-average recurrence in
  \f$ f(n) \sum_{k=0}^{n-1} x(k) + g(n) \f$,
  which is stored in the expression \p recurrence_rhs.
  Assume that the system has produced the expression
  \p exact_solution_.expression(), which is a function of the
  variable \f$ n \f$.

  The verification's process is divided in 2 steps:
  -  Validation of the \ref initial_condition
     "symbolic initial condition x(0)".
     Evaluate the expression \p exact_solution_.expression() for
     \f$ n = 0 \f$ and simplify as much as possible.
     If the final result is <EM>synctactically</EM> equal to \f$ x(0) \f$
     the symbolic initial conditions is verified and we can proceed to
     step 2; otherwise return <CODE>INCONCLUSIVE_VERIFICATION</CODE>
     because the solution can be wrong or it is not enough simplified.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
  -  Validation of the solution using substitution.
     Consider the difference
     \f$ d = x(n) - (f(n) sum(k, 0, n-1, x(k)) + g(n)) \f$, with
     \f$ x(n) \f$ and \f$ x(k) \f$ replaced with the solution.
     Since the closed formula for \f$ x(n) \f$ is guaranteed to hold for
     \f$ n >= 1 \f$ only, the lower limit of the sum must start from
     \f$ 1 \f$ and the term for \f$ n = 0 \f$ must be considered before.
     Hence, the difference that we consider is
     \f$ d = x(n) - (f(n) x(0) + f(n) sum(k, 1, n - 1, x(k)) + g(n)) \f$:
     - if \f$ d = 0 \f$     -> returns <CODE>PROVABLY_CORRECT</CODE>:
                               the solution is certainly right.
     - if \f$ d \neq 0 \f$  -> returns <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
 			       the solution can be wrong or we failed to
			       simplify it.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_weighted_average() const {
  Expr weight_rec = weight();
  
  // The case `f(n) = -1', i.e. the recurrence has the form
  // `x(n) = - sum(k, 0, n-1, x(k)) + g(n)', is special:
  // the solution is simply `x(n) = g(n) - g(n-1)'.
  // FIXME: the traditional validation' process does not work,
  // is it true?
  if (weight_rec == -1)
    // FIXME: verify!!!
    return PROVABLY_CORRECT;
  
  // Note: the solution is valid only for `n > 0'.
 
  // Step 1: validation of the initial condition.
  Expr e = exact_solution_.expression().substitute(n, 1);
  e = simplify_all(e);
  if (e != weight_rec.substitute(n, 1) * x(0)
      + inhomogeneous_term.substitute(n, 1))
    // FIXME: provably_incorrect...
    return INCONCLUSIVE_VERIFICATION;
  
  // Step 2: validation of the solution.
  // Consider the `sum(k, 0, n-1, x(k)', with `x(k)' replaced by
  // the solution, and tries to semplify it.
  Symbol h;
  Expr diff
    = PURRS::sum(h, 1, n - 1, exact_solution_.expression().substitute(n, h));
  diff = simplify_sum(diff, COMPUTE_SUM);
  // Consider the difference
  // `x(n) - (f(n) x(0) + f(n) sum(k, 1, n - 1, x(k)) + g(n))'
  // and tries to simplify it.
  diff = exact_solution_.expression()
    - diff * weight_rec - x(0) * weight_rec - inhomogeneous_term;
  diff = simplify_all(diff);
  if (diff == 0)
    return PROVABLY_CORRECT;
  else
    return INCONCLUSIVE_VERIFICATION;
}

PURRS::Recurrence::Verify_Status
PURRS::Recurrence::validate_initial_conditions(Bound kind_of_bound,
					       const Expr& bound,
					       unsigned int index) const {
  Expr bound_valuated = bound.substitute(n, index);
  D_VAR(bound_valuated);
  if (bound_valuated != x(index))
    if (bound_valuated.is_a_mul()) {
      Expr coeff_ic = 1;
      for (unsigned int i = bound_valuated.nops(); i-- > 0; ) {
	const Expr& factor = bound_valuated.op(i);
	if (factor != x(index))
	  coeff_ic *= factor;
      }
      D_VAR(coeff_ic);
      Number num;
      if (kind_of_bound == UPPER) {
	if (coeff_ic.is_a_number(num) && num < 1)
	  return INCONCLUSIVE_VERIFICATION;
      }
      else
	if (coeff_ic.is_a_number(num) && num > 1)
	  return INCONCLUSIVE_VERIFICATION;
    }
    else if (bound_valuated.is_a_add()) {
      Expr term_with_ic = 0;
      Expr other_terms = 0;
      split(bound_valuated, x(index),
	    term_with_ic, other_terms);
      Expr coeff_ic = 1;
      if (term_with_ic.is_a_mul()) {
	for (unsigned int i = term_with_ic.nops(); i-- > 0; ) {
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
	if (kind_of_bound == UPPER) {
	  if (!(num_coeff >= 1 && num_other.is_positive()))
	    return INCONCLUSIVE_VERIFICATION;
	}
	else
	  if (!(num_coeff <= 1 && num_other.is_negative()))
	    return INCONCLUSIVE_VERIFICATION;
      }
      else
	return INCONCLUSIVE_VERIFICATION;
    }
  return PROVABLY_CORRECT;
}

/*!
  Consider the right hand side \p rhs of the functional equation
  \f$ a x_{n/b} + p(n) \f$.
  If \p kind_of_bound is <CODE>UPPER</CODE> we try to check the
  correctness of the upper bound
  If \p kind_of_bound is <CODE>LOWER</CODE> we try to check the
  correctness of the lower bound
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_bound(Bound kind_of_bound) const{
  assert(is_functional_equation());
  assert(kind_of_bound == UPPER || kind_of_bound == LOWER);
  Expr bound;
  if (kind_of_bound == UPPER)
    bound = simplify_sum(upper_bound_.expression(), REWRITE_UPPER_LIMIT);
  else
    bound = simplify_sum(lower_bound_.expression(), REWRITE_UPPER_LIMIT);
  
  // Step 1: validation of initial conditions.
  Verify_Status status
    = validate_initial_conditions(kind_of_bound, bound,
				  applicability_condition());
  if (status == INCONCLUSIVE_VERIFICATION || status == PROVABLY_INCORRECT)
    return status;
  
  // Step 2: find `partial_bound'.
  // We not consider the terms containing the initial conditions:
  // `partial_bound' will contain all the other terms.
  Expr partial_bound = 0;
  if (bound.is_a_add())
    for (unsigned int i = bound.nops(); i-- > 0; ) {
      if (!has_symbolic_initial_conditions(bound.op(i)))
	partial_bound += bound.op(i);
    }
  else
    if (!has_symbolic_initial_conditions(bound))
      partial_bound = bound;
  D_VAR(partial_bound);
  // The recurrence is homogeneous.
  if (partial_bound == 0)
    return PROVABLY_CORRECT;
  
  // Step 3: verification of the inductive base.
  Number num;
  if (kind_of_bound == UPPER
      && partial_bound
      .substitute(n, applicability_condition()).is_a_number(num)
      && num.is_negative())
    return INCONCLUSIVE_VERIFICATION;
  if (kind_of_bound == LOWER
      && partial_bound
      .substitute(n, applicability_condition()).is_a_number(num)
      && num.is_positive())
    return INCONCLUSIVE_VERIFICATION;
  
  // Step 4: verification of the inductive step.
  Expr partial_bound_sub
    = partial_bound.substitute(n, n / functional_eq_p->ht_begin()->first);
  Expr approx
    = recurrence_rhs.substitute(x(n / functional_eq_p->ht_begin()->first),
				partial_bound_sub);
  D_VAR(approx);
  approx = simplify_ex_for_input(approx, true);
  approx = simplify_logarithm(approx);
  D_VAR(approx);
  
  Expr diff;
  if (kind_of_bound == UPPER)
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
    for (unsigned int i = diff.nops(); i-- > 0; ) {
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
  return INCONCLUSIVE_VERIFICATION;
}

//! \brief
//! Verify the exact solution of recurrence \p *this, where the recurrence
//! can be linear of finite order, non-linear of finite order,
//! weighted-average or a functional equation.
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_exact_solution() const {
  if (!exact_solution_.has_expression())
    throw std::logic_error("PURRS::Recurrence::verify_exact_solution() "
			   "called, but no exact solution was computed");
  // Case 1 and case 2: linear recurrences of finite order
  // and non-linear recurrences of finite order. 
  if (is_linear_finite_order() || is_non_linear_finite_order())
    return verify_finite_order();
  // Case 3: weighted-average recurrence.
  else if (is_weighted_average())
    return verify_weighted_average();
  else {
    assert(is_functional_equation());
    // Case 4 (special case): functional equation of the form
    // `x(n) = x(n/b)'. 
    return INCONCLUSIVE_VERIFICATION;
  }
}

PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_lower_bound() const {
  if (exact_solution_.has_expression())
    return verify_exact_solution();
  else if (lower_bound_.has_expression())
    if (is_functional_equation())
      return verify_bound(LOWER);
    else
      // At the moment this case is impossible!
      return INCONCLUSIVE_VERIFICATION;
  else
    throw std::logic_error("PURRS::Recurrence::verify_lower_bound() "
			   "called, but no lower bound was computed");
}

PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_upper_bound() const {
  if (exact_solution_.has_expression())
    return verify_exact_solution();
  else if (upper_bound_.has_expression())
    if (is_functional_equation())
      return verify_bound(UPPER);
    else
      // At the moment this case is impossible!
      return INCONCLUSIVE_VERIFICATION;
  else
    throw std::logic_error("PURRS::Recurrence::verify_upper_bound() "
			   "called, but no upper bound was computed");
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
    // is necessary to classify the linear recurrence relation
    // associated to that one non-linear.
    // FIXME: rivedere la verifica per capire se questo passaggio e'
    // davvero necessario.
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
      if ((status = associated_linear_rec().solve_linear_finite_order())
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
						    base_exp_log(), false);
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
						base_exp_log(), false);
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

/*!
  Builds the weighted-average recurrence
  \f$ x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n) \f$, where
  \f$ f(n) \f$ is stored in \p weight; \f$ g(n) \f$ is stored
  in \p inhomogeneous and \p first_valid_index contains the least
  non-negative integer \f$ j \f$ such that the recurrence is
  well-defined for \f$ n \geq j \f$.
  If the system is able to solve the recurrence, then this function
  returns <CODE>SUCCESS</CODE> and the solution is stored in \p solution.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::solve_new_weighted_average_rec(const Expr& weight,
						  const Expr& inhomogeneous,
						  index_type first_valid_index,
						  Expr& solution) const {
  Symbol h;
  Recurrence rec_rewritten(weight * PURRS::sum(h, 0, n-1, x(h))
			   + inhomogeneous);
  const Expr& coeff_first_order
    = simplify_all(weight / weight.substitute(n, n-1)
		   * (1 + weight.substitute(n, n-1)));
  const Expr& inhomog_first_order
    = simplify_all(weight * (inhomogeneous / weight
			     - (inhomogeneous / weight).substitute(n, n-1)));
  rec_rewritten.weighted_average_p
    = new Weighted_Average_Info(Recurrence(coeff_first_order*x(n-1)
					   +inhomog_first_order), weight);
  rec_rewritten.set_weighted_average();
  associated_first_order_rec().set_first_valid_index(first_valid_index);
  rec_rewritten.set_inhomogeneous_term(inhomogeneous);
  return rec_rewritten.compute_weighted_average_recurrence(solution);
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
      return apply(e.functor(), increase_argument_x_function(e.arg(0), num));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int j = 0; j < num_argument; ++j)
	argument[j] = increase_argument_x_function(e.arg(j), num);
      return apply(e.functor(), argument);
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
  - shift forward of the first order recurence: \f$ n = n + 1 \f$.
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
      // Shift backward: n -> n - 1.
      solution = solution.substitute(n, n - 1);
      solution = solution
	.substitute(x(associated_first_order_rec().first_valid_index()),
		    x(associated_first_order_rec().first_valid_index()+1));
      // If there is the initial condition `x(1)' specified then
      // the system substitute it with the respective value; otherwise
      // the system performs the substitution `x(1) = 2*x(0)+1'.
      std::map<unsigned int, Expr>::const_iterator i
	= initial_conditions.find(1);
      if (i != initial_conditions.end())
	solution = solution.substitute(x(1), get_initial_condition(1));
      else
	solution = solution
	  .substitute(x(associated_first_order_rec().first_valid_index()+1),
		      (weight()*x(0)+inhomogeneous_term).substitute(n, 1));
      //	solution = simplify_ex_for_output(solution, false);
      return SUCCESS;
    }
    else
      return status;
  }
}

//! \brief
//! Let \p solution_or_bound be the expression that represent the
//! solution or the bound computed for the recurrence \p *this.
//! This function substitutes eventual symbolic initial conditions
//! specified by the user shifting the solution or the bound if necessary.
/*!
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
  Finally substitutes to the arbitrary initial conditions in the solution or
  in the bound the eventual values specified by the user.
*/
Expr
PURRS::Recurrence::
substitute_i_c_shifting(const Expr& solution_or_bound) const {
  assert(!initial_conditions.empty());
  Expr sol_or_bound = solution_or_bound;
  index_type first_valid_index_rhs;
  index_type order_or_rank;
  if (is_linear_finite_order()) {
    first_valid_index_rhs = first_valid_index();
    order_or_rank = order();
  }
  else if (is_functional_equation()) {
    first_valid_index_rhs = applicability_condition();
    order_or_rank = rank();
  }
  else {
    assert(is_non_linear_finite_order() || is_weighted_average());
    // FIXME: we must again understand how to work in these cases.
    return sol_or_bound;
  }

  // If the order of linear recurrences is zero than it does not perform
  // shifts (the rank of functional equations can not to be zero, it is
  // greater or equal to one).
  if (order_or_rank != 0) {
    // Consider the maximum index of `x' function in the map
    // `initial_conditions', i.e. the largest index of the initial
    // conditions inserted by the user.
    unsigned int max_index_user_i_c = 0;
    for (std::map<unsigned int, Expr>::const_iterator i
	   = initial_conditions.begin(),
	   iend = initial_conditions.end(); i != iend; ++i)
      if (i->first > max_index_user_i_c)
	max_index_user_i_c = i->first;
    
    // `max_index_symb_i_c' represents the largest
    // index of the symbolic initial conditions.
    unsigned int max_index_symb_i_c
      = first_valid_index_rhs + order_or_rank - 1;
    // If `max_index_user_i_c' is bigger than `max_index_symb_i_c',
    // then we must shift the solution or the bound.
    // There are two different steps:
    // - 1. shift the index of the symbolic initial conditions;
    // - 2. shift the index of the recurrence `n'.
    if (max_index_user_i_c > max_index_symb_i_c) {
      // Step 1.
      for (index_type i = 0; i < order_or_rank; ++i)
	sol_or_bound = sol_or_bound.substitute(x(i + first_valid_index_rhs),
					       x(i + max_index_user_i_c
						 - order_or_rank + 1));
      // Step 2..
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
  for (std::map<unsigned int, Expr>::const_iterator i
	 = initial_conditions.begin(),
	 iend = initial_conditions.end(); i != iend; ++i)
    sol_or_bound = sol_or_bound.substitute(x(i->first),
					   get_initial_condition(i->first));
  sol_or_bound = simplify_numer_denom(sol_or_bound);
  sol_or_bound = simplify_ex_for_output(sol_or_bound, false);
  return sol_or_bound;
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
  
  // Check if there are specified initial conditions and in this case
  // eventually shift the solution in according with them before to
  // substitute the values of the initial conditions to the
  // symbolic initial condition `x(i)'.
  if (!initial_conditions.empty())
    exact_solution_.set_expression
      (substitute_i_c_shifting(exact_solution_.expression()));
  lower_bound_.set_expression(exact_solution_.expression());
  upper_bound_.set_expression(exact_solution_.expression());
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution_functional_equation() const {
  Solver_Status status;
  if ((status = approximate_functional_equation(LOWER)) == SUCCESS
      && (status = approximate_functional_equation(UPPER)) == SUCCESS
      && lower_bound_.expression() == upper_bound_.expression()) {
    // Check if there are specified initial conditions and in this case
    // eventually shift the solution in according with them before to
    // substitute the values of the initial conditions to the
    // symbolic initial condition `x(i)'.
    if (!initial_conditions.empty())
      lower_bound_.set_expression
	(substitute_i_c_shifting(lower_bound_.expression()));
    upper_bound_.set_expression(lower_bound_.expression());
    exact_solution_.set_expression(lower_bound_.expression());
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
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_exact_solution_weighted_average() const {
  Expr solution;
  Solver_Status status = compute_weighted_average_recurrence(solution);
  if (status != SUCCESS)
    return status;

#if 0
  // FIXME: to improve the function `substitute_i_c_shifting()'
  // so that to accept also weighted-average recurrences.
  // Check if there are specified initial conditions and in this case
  // eventually shift the solution in according with them before to
  // substitute the values of the initial conditions to the
  // symbolic initial condition `x(i)'.
  if (!initial_conditions.empty())
    exact_solution_.set_expression
      (substitute_i_c_shifting(exact_solution_.expression()));
  lower_bound_.set_expression(exact_solution_.expression());
  upper_bound_.set_expression(exact_solution_.expression());
#else
  // FIXME: At the moment we substitute here only the initial
  // condition `x(0)'.
  std::map<unsigned int, Expr>::const_iterator i
    = initial_conditions.find(0);
  if (i != initial_conditions.end())
    solution = solution.substitute(x(0), get_initial_condition(0));
  exact_solution_.set_expression(solution);
  lower_bound_.set_expression(solution);
  upper_bound_.set_expression(solution);
#endif
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
  if (exact_solution_.has_expression())
    return SUCCESS;

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
    default:
      throw std::runtime_error("PURRS internal error: "
			       "compute_exact_solution().");
    }
  else
    // return the `Solver_Status' associated to `classifier_status'.
    return map_status(classifier_status);
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
  if (!initial_conditions.empty())
    if (kind_of_bound == LOWER)
      lower_bound_.set_expression
	(substitute_i_c_shifting(lower_bound_.expression()));
    else
      upper_bound_.set_expression
	(substitute_i_c_shifting(upper_bound_.expression()));

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
  if (lower_bound_.has_expression())
    return SUCCESS;

  // We may not have the lower bound explicitely, yet we may have
  // the exact solution.
  if (exact_solution_.has_expression()) {
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
  if (upper_bound_.has_expression())
    return SUCCESS;

  // We may not have the upper bound explicitely, yet we may have
  // the exact solution.
  if (exact_solution_.has_expression()) {
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
    default:
      throw std::runtime_error("PURRS internal error: "
			       "compute_upper_bound().");
    }
  else
    // return the `Solver_Status' associated to `classifier_status'.
    return map_status(classifier_status);
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
    for (std::map<unsigned int, Expr>::const_iterator i
	   = initial_conditions.begin(),
	   initial_conditions_end = initial_conditions.end();
	 i != initial_conditions_end; ++i)
      s << "  x(" << i->first << ")"
	<< " = " << i->second << std::endl;
  }
  
  //(*functional_eq_p).dump_homogeneous_terms(s);
  s << std::endl;
}
