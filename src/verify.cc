/* To be written.
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
http://www.cs.unip.it/purrs/ . */

#include <config.h>

#define NEW_VERIFICATION 1

#include "ep_decomp.hh"
#include "fact_decomp.hh"
#include "finite_order.hh"
#include "simplify.hh"
#include "Expr.defs.hh"
#include "Cached_Expr.defs.hh"
#include "Blackboard.defs.hh"
#include "Recurrence.defs.hh"
#include "Recurrence.inlines.hh"
#include <vector>
#include <algorithm>

namespace PURRS = Parma_Recurrence_Relation_Solver;

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
validate_initial_conditions(bool is_symbolic_solution, index_type order,
			    const std::vector<Expr>& coefficients_i_c,
			    const Expr& summands_with_i_c,
			    const Expr& summands_without_i_c) const {
  if (!is_symbolic_solution) {
    for (index_type i = 0; i < order; ++i) {
      index_type index =
	get_max_index_initial_condition() > first_valid_index + i
	? get_max_index_initial_condition() : first_valid_index + i;
      Expr e = simplify_all(summands_without_i_c.substitute(n, index));
      e = blackboard.rewrite(e);
      e = simplify_all(e);
      // Get from the map `initial_conditions_' the value associated
      // to the symbolic initial conditions `x(index)'.
      if (e != get_initial_condition(index))
	// FIXME: provably_incorrect nei casi semplici.
	return INCONCLUSIVE_VERIFICATION;
    }
    return PROVABLY_CORRECT;
  }

  for (index_type i = 0; i < order; ++i) {
    index_type index = first_valid_index + i;
    // In the case of non-linear recurrences the homogeneous part
    // of the solution can contains more than one symbolic initial
    // condition: the vector `coefficients_i_c' is in this case empty.
    if (is_non_linear_finite_order()) {
      // The expression `e' can be more difficult to simplify.
      // For this motive we performed simplification also
      // before to expand blackboard's definitions.
      Expr e = simplify_all(summands_with_i_c.substitute(n, index));
      // Expand blackboard's definitions in order to increase the
      // opportunities for simplifications.
      e = blackboard.rewrite(e);
      e = simplify_all(e);
      if (e != x(index))
	// FIXME: provably_incorrect nei casi semplici.
	return INCONCLUSIVE_VERIFICATION;
    }
    else {
      unsigned int coefficients_i_c_size = coefficients_i_c.size();
      // coefficients of initial conditions.
      for (unsigned j = 0; j < coefficients_i_c_size; ++j) {
	D_VAR(coefficients_i_c[j]);
	Expr e = simplify_all(coefficients_i_c[j].substitute(n, index));
	// Expand blackboard's definitions in order to increase the
	// opportunities for simplifications.
	e = blackboard.rewrite(e);
	e = simplify_all(e);
	D_VAR(index);
	if (index == first_valid_index + j) {
	  if (e != 1)
	    // FIXME: provably_incorrect nei casi semplici.
	    return INCONCLUSIVE_VERIFICATION;
	}
	else
	  if (!e.is_zero())
	    // FIXME: provably_incorrect nei casi semplici.
	    return INCONCLUSIVE_VERIFICATION;
      }
    }
    // The non-homogeneous part of the solution.
    Expr e = simplify_all(summands_without_i_c.substitute(n, index));
    e = blackboard.rewrite(e);
    e = simplify_all(e);
    if (!e.is_zero())
      // FIXME: provably_incorrect nei casi semplici.
      return INCONCLUSIVE_VERIFICATION;
  }
  return PROVABLY_CORRECT;
}

PURRS::Recurrence::Verify_Status
PURRS::Recurrence::
traditional_step_3(index_type order_rec, const Expr& e) const {
  // By substitution, verifies that `e' satisfies the homogeneous part
  // of the recurrence.
  // Computes `substituted_homogeneous_rhs' by substituting, in the
  // hoomogeneous part of the recurrence, `n' by `n - d' (where `d' is
  // the decrement of the i-th term `a_i(n)*x(n - d)').
  Expr substituted_homogeneous_rhs = recurrence_rhs - inhomogeneous_term;
  // Substitutes in the homogeneous part of the recurrence the terms
  // of the form `x(n-i)'.
  // FIXME: the code in the "else" is, from a theoretical point of view,
  // better for efficiency because it performs only 1 symbolic substitution
  // and d numeric substitutions, but in practice is not so.
#if 1
  for (index_type d = 1; d <= order_rec; ++d)
    if (substituted_homogeneous_rhs.has(x(n-d))) {
      Expr shifted_solution
	= simplify_all(e.substitute(n, n - d));
      shifted_solution = simplify_sum(shifted_solution, REWRITE_UPPER_LIMIT);
      substituted_homogeneous_rhs
	= substituted_homogeneous_rhs.substitute(x(n - d), shifted_solution);
    }
#else
  Symbol h;
  Expr shifted_solution = simplify_all(e.substitute(n, n - h));
  for (index_type d = 1; d <= order_rec; ++d)
    if (substituted_homogeneous_rhs.has(x(n-d))) {
      Expr shifted_solution_num
	= shifted_solution.substitute(h, d);
      shifted_solution_num = simplify_sum(shifted_solution_num,
					  REWRITE_UPPER_LIMIT);
      substituted_homogeneous_rhs
	= substituted_homogeneous_rhs.substitute(x(n - d),
						 shifted_solution_num);
    }
#endif
  Expr diff = e - substituted_homogeneous_rhs;
  // The expression `diff' can be more difficult to simplify.
  // For this motive we performed simplification also
  // before to expand blackboard's definitions.
  diff = blackboard.rewrite(diff);
  diff = simplify_all(diff);
  if (!diff.is_zero()) {
    diff = simplify_all(diff);
    if (!diff.is_zero())
      return INCONCLUSIVE_VERIFICATION;
  }
  return PROVABLY_CORRECT;
}

// FIXME: to finish comment!!!
/*!
  Find the maximum degree of a polynomial that may occur in the
  solution in the following way:
  we assume that \f$ \lambda_1 \f$, ..., \f$ \lambda_m \f$ are the
  roots of the characteristic equation, of multiplicity
  \f$ r_1 \f$, ..., \f$ r_m \f$ respectively, where \f$ \lambda_i \f$
  are complex numbers all distict and  the multiplicity are positive
  integer (\f$ r_1 +...+ r_m = k \f$ where \f$ k \f$ is the order
  of the recurrence).
  Moreover, we assume that the non-homogeneous part of the recurrence
  has the following form
  \f[
    p_1(n) \alpha_1^n + ... + p_j(n) \alpha_j^n
  \f]
  where \f$ \alpha_i \f$ are complex numbers all distinct and the
  \f$ p_i \f$ are polynomials (different from \f$ 0 \f$).
  We set
  \f[
    q_1 =
    \begin{cases}
      r_1 - 1,
         \quad \text{if } \lambda_1 \notin \{ \alpha_1, ..., \alpha_j \}\\
      r_1 - 1, + \deg(p_i)
         \quad \text{if } \lambda_1 = \alpha_i.\\
    \end{cases}
  \f]
  Analogously for \f$ q_2 \f$, ..., \f$ q_m \f$.
  The maximum degree of a polynomial that may occur in the
  solution is the maximum of the set
  \f[
    \{ q_1, ..., q_m, \deg(p_1), ..., \deg(p_j) \}.
  \f]
  Since \f$ r_i \le k \f$ is possible to consider the MAGGIORAZIONE
  \f[
    k - 1 + \max \{ \deg(p_1), ..., \deg(p_j) \}.
  \f]
  This method can seem inefficient because builds memory's locations
  that will not be used, but is not necessary the comparison between
  the roots \f$ \lambda_i \f$ and the bases of the exponentials
  \f$ \alpha_j \f$ (we also observe that in the memory's locations
  unused there is \f$ 0 \f$ and so we can avoid useless checks.
  In the following code we use the MAGGIORAZIONE
  \f[
    max_polynomial_degree
      =\max\{ r_1, ..., r_m \} - 1 + \max \{ \deg(p_1), ..., \deg(p_j) \}.
  \f]
*/
bool
PURRS::Recurrence::
verify_new_method_const_coeff(index_type order_rec, const Expr& e,
			      const std::vector<Polynomial_Root> roots,
			      bool inhomogeneous_part) const {
  // The maximum degree of a polynomial that may occur in the solution.
  unsigned int max_polynomial_degree = 0;
  unsigned int num_of_exponentials = 0;
  if (inhomogeneous_part) {
    std::vector<Expr> bases_of_exp;
    std::vector<Expr> exp_poly_coeff;
    std::vector<Expr> exp_no_poly_coeff;
    exp_poly_decomposition(inhomogeneous_term.expand(), Recurrence::n,
			   bases_of_exp, exp_poly_coeff, exp_no_poly_coeff);
    
    assert(bases_of_exp.size() == exp_poly_coeff.size()
	   && exp_poly_coeff.size() == exp_no_poly_coeff.size()
	   && exp_no_poly_coeff.size() >= 1);
    
    num_of_exponentials = bases_of_exp.size();
    for (unsigned int i = 0; i < num_of_exponentials; ++i) {
      if (!exp_no_poly_coeff[i].is_zero()) {
	D_MSGVAR("No poly: ", exp_no_poly_coeff[i]);
	return false;
      }
      max_polynomial_degree = std::max(max_polynomial_degree,
				       exp_poly_coeff[i].degree(n));
    }
  }
  unsigned int max_multiplicity = 1;
  for (unsigned int i = 0, nroots = roots.size(); i < nroots; ++i)
    if (roots[i].multiplicity() > max_multiplicity)
      max_multiplicity = roots[i].multiplicity();
  max_polynomial_degree += max_multiplicity - 1;

  Expr substituted_rhs = recurrence_rhs;
  if (!inhomogeneous_part)
    substituted_rhs -= inhomogeneous_term;
  
  // FIXME: the code in the "else" is, from a theoretical point of view,
  // better for efficiency because it performs only 1 symbolic substitution
  // and d numeric substitutions, but in practice is not so.
#if 1
  for (index_type i = order_rec; i-- > 0; )
    if (substituted_rhs.has(x(n - (i + 1)))) {
      const Expr& shifted_solution = e.substitute(n, n - (i + 1));
      //shifted_solution = simplify_sum(shifted_solution, REWRITE_UPPER_LIMIT);
      substituted_rhs = substituted_rhs
	.substitute(x(n - (i + 1)), shifted_solution);
    }
#else
  Symbol h;
  Expr shifted_solution = simplify_all(e.substitute(n, n - h));
  for (index_type i = order_rec; i-- > 0; )
    if (substituted_rhs.has(x(n - (i + 1)))) {
      const Expr& shifted_solution_num
	= shifted_solution.substitute(h, i + 1);
      //shifted_solution_num = simplify_num(shifted_solution_num,
      //                                    REWRITE_UPPER_LIMIT);
      substituted_rhs = substituted_rhs
	.substitute(x(n - (i + 1)), shifted_solution_num);
    }
#endif
  // FIXME: ritardare la riscrittura.
  Expr diff = blackboard.rewrite(e - substituted_rhs);
  diff = diff.expand();
  D_VAR(diff);
  if (diff == 0)
    return true;

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
	    k = factor.arg(1).ex_to_number().to_unsigned_int();
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
	    unsigned int k = factor.arg(1).ex_to_number().to_unsigned_int();
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
	assert(summand.is_a_number() || summand == n
	       // FIXME: added this condition requested, for example,
	       // by the recurrence number 904 (12-12-03) of the file heap:
	       // x(n) = a*x(n-1)+b+c*(n-1). Is it right?
	       || summand.is_a_symbol());
	if (summand.is_a_number() || summand.is_a_symbol())
	  coefficients_of_exponentials[0] += summand;
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
      assert(summand.is_a_number() || summand == n
	     // FIXME: added this condition, is it right?
	     || summand.is_a_symbol());
      Number k;
      if (summand.is_a_number(k))
	coefficients_of_exponentials[0] += k;
      else
	coefficients_of_exponentials[1] += 1;
    }
  }
  
  D_VEC(coefficients_of_exponentials, 0, max_polynomial_degree);
  Number num_tests = num_of_exponentials + order_rec;
  for (unsigned int i = 0; i <= max_polynomial_degree; ++i) {
    if (!coefficients_of_exponentials[i].is_zero()) {
      // Not syntactically 0: try to prove that is it semantically 0.
      Expr c = coefficients_of_exponentials[i];
      for (Number r = 0; r < num_tests; ++r) {
	Expr c_r = simplify_all(c.substitute(n, r));
	if (!c_r.is_zero()) {
	  D_MSGVAR("Argh!!! ", i);
	  D_MSGVAR("Argh!!! ", r);
	  D_MSGVAR("Argh!!! ", c_r);
	  return false;
	}
      }
    }
  }
  return true; 
}

bool
PURRS::Recurrence::
verify_new_method_var_coeff(index_type order_rec,
			    const Expr& summands_without_i_c) const {
  // If `summands_without_i_c' contains `sum()', `prod()'
  // we do not apply this method.
  if (summands_without_i_c.has_sum_or_prod_function())
    return false;

  Expr substituted_rhs = recurrence_rhs;
  for (index_type i = order_rec; i-- > 0; ) {
    const Expr& shifted_solution
      = summands_without_i_c.substitute(n, n - (i + 1));
    //shifted_solution = simplify_sum(shifted_solution, REWRITE_UPPER_LIMIT);
    substituted_rhs = substituted_rhs
      .substitute(x(n - (i + 1)), shifted_solution);
  }
  Expr diff = blackboard.rewrite(summands_without_i_c - substituted_rhs);
  diff = diff.expand();
  D_VAR(diff);
  D_VAR(diff.numerator());
  diff = diff.numerator();
  if (diff == 0)
    return true;

  // FIXME: to be finished!
  std::vector<Number> argument_factorials;
  std::vector<Expr> coeff_factorials;
  factorial_decomposition(diff, n, argument_factorials, coeff_factorials);
  int position = 0;
  for (unsigned int i = argument_factorials.size(); i-- > 0; )
    if (argument_factorials[i] > position)
      position = i;
  Expr coeff_dominant = coeff_factorials[position];
  // FIXME:
  // To call the method `verify_new_method_const_coeff()' (only the
  // necessary part) on `coeff_dominant'.
  return false;
}

namespace {
using namespace PURRS;

void
split_solution(bool is_symbolic_solution, const Expr& exact_solution,
	       Expr& summands_with_i_c, Expr& summands_without_i_c) {
  if (is_symbolic_solution)
    if (exact_solution.is_a_add())
      for (unsigned int i = exact_solution.nops(); i-- > 0; ) {
	const Expr& addend_exact_solution = exact_solution.op(i);
	if (has_at_least_a_symbolic_initial_condition(addend_exact_solution))
	  summands_with_i_c += addend_exact_solution;
	else
	  summands_without_i_c += addend_exact_solution;
      }
    else
      if (has_only_symbolic_initial_conditions(exact_solution))
	summands_with_i_c = exact_solution;
      else
	summands_without_i_c = exact_solution;
  else
    summands_without_i_c = exact_solution;
  D_VAR(summands_without_i_c);
}

void
fill_vector_coefficients_i_c(const Expr& summands_with_i_c, unsigned int gcd,
			     unsigned int first_valid_index,
			     std::vector<Expr>& coefficients_i_c) {
  D_VAR(summands_with_i_c);
  if (summands_with_i_c.is_a_add()) {
    std::vector<Expr> index_i_c;
    for (unsigned int i = summands_with_i_c.nops(); i-- > 0; ) {
      const Expr& addend = summands_with_i_c.op(i);
      if (addend.is_a_mul()) {
	Expr arg = 0;
	Expr coeff = 1;
	for (unsigned int j = addend.nops(); j-- > 0; ) {
	  const Expr& factor = addend.op(j);
	  if (factor.is_the_x_function()) {
	    assert(arg == 0);
	    arg = factor.arg(0);
	  }
	  else
	    coeff *= factor;
	}
	Number numeric_arg;
	if (arg.is_a_number(numeric_arg)) {
	  assert(numeric_arg.is_nonnegative_integer());
	  assert(numeric_arg.to_unsigned_int()-first_valid_index
		 < coefficients_i_c.size());
	  coefficients_i_c[numeric_arg.to_unsigned_int()-first_valid_index]
	    += coeff;
	}
	else if (arg.is_the_mod_function())
	  coefficients_i_c[0] += coeff;
	else {
	  assert(arg.is_a_add() && arg.nops() == 2);
	  assert(arg.op(0).is_a_number() || arg.op(1).is_a_number());
	  unsigned int num = arg.op(0).is_a_number()
	    ? arg.op(0).ex_to_number().to_unsigned_int()
	    : arg.op(1).ex_to_number().to_unsigned_int();
	  assert(num/gcd-first_valid_index < coefficients_i_c.size());
	  coefficients_i_c[num/gcd-first_valid_index] += coeff;
	}
      }
      else {
	assert(addend.is_the_x_function());
	const Expr& arg = addend.arg(0);
	Number numeric_arg;
	if (arg.is_a_number(numeric_arg)) {
	  assert(numeric_arg.is_nonnegative_integer());
	  assert(numeric_arg.to_unsigned_int()-first_valid_index
		 < coefficients_i_c.size());
	  coefficients_i_c[numeric_arg.to_unsigned_int()-first_valid_index]
	    += 1;
	}
	else if (arg.is_the_mod_function())
	  coefficients_i_c[0] += 1;
	else {
	  assert(arg.is_a_add() && arg.nops() == 2
		 && (arg.op(0).is_a_number() || arg.op(1).is_a_number()));
	  unsigned int num = arg.op(0).is_a_number()
	    ? arg.op(0).ex_to_number().to_unsigned_int()
	    : arg.op(1).ex_to_number().to_unsigned_int();
	  assert(num/gcd-first_valid_index < coefficients_i_c.size());
	  coefficients_i_c[num/gcd-first_valid_index] += 1;
	}
      }
    }
  }
  else if (summands_with_i_c.is_a_mul()) {
    coefficients_i_c[0] = 1;
    for (unsigned int j = summands_with_i_c.nops(); j-- > 0; ) {
      const Expr& factor = summands_with_i_c.op(j);
      if (!factor.is_the_x_function())
	coefficients_i_c[0] *= factor;
    }
  }
  else {
    // FIXME: chiamare metodo `is_a_symbolic_initial_condition()'.
    assert(summands_with_i_c.is_the_x_function());
    coefficients_i_c[0] = 1;
  }
}

} // anonymous namespace

/*!
  Case 1: linear recurrences of finite order.
  Consider the right hand side \p rhs of the order \f$ k \f$ recurrence
  relation
  \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dots + a_k * x_{n-k} + p(n) \f$,
  which is stored in the expression \p recurrence_rhs.
  Let \f$ i \f$ be the \ref first_valid_index "first_valid_index" and
  assume that the system has produced the expression
  \p exact_solution_.expression(), which is a function of the
  variable \f$ n \f$.

  The verification's process is divided in 4 steps:
  -  Split \p exact_solution_.expression() in 2 expressions:
     \p summands_with_i_c contains the summands with an occurrence of a
     symbolic initial conditions \f$ x(i), \cdots, x(i+k) \f$;
     \p summands_without_i_c contains all the other summands.
     Consider the homogeneous part of the solution: the coefficients of
     the symbolic initial conditions \f$ x(i), \cdots, x(i+k) \f$ are put
     in the vector \p coefficients_i_c.
  -  Validation of \ref initial_conditions "symbolic initial conditions".
     Evaluate the elements of the vector \p coefficients_i_c for
     \f$ n = i, \cdots, i+k-1 \f$ and simplify as much as possible.
     Only the element corresponding to the \f$ i \f$-th symbolic initial
     condition must be <EM>synctactically</EM> equal to \f$ 1 \f$, while
     all the other elements must be <EM>synctactically</EM> equal to
     \f$ 0 \f$: if it is not so the function returns
     <CODE>INCONCLUSIVE_VERIFICATION</CODE> because
     the solution can be wrong or it is not enough simplified.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
     If the previous verification on the coefficients of symbolic
     initial conditions is successfully, then evaluate \p summands_without_i_c
     for \f$ n = i, \cdots, i+k-1 \f$ and simplify as much as possible.
     If the final result is not <EM>synctactically</EM> equal to \f$ 0 \f$
     the function returns <CODE>INCONCLUSIVE_VERIFICATION</CODE> because the
     solution can be wrong or it is not enough simplified.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
  -  Verify that the elements of the vector \p coefficients_i_c satisfy
     the homogeneous part of the recurrence
     \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} \f$.
     There are two different way:
     - Replace \f$ x(n-i) \f$ by the \f$ i \f$-th elements of the vector
       evaluated at \f$ n-i \f$ (for \f$ i = 1, \cdots, k \f$) in the
       above expression, and store the result in
       \p substituted_homogeneous_rhs.
       For each element of the vector consider the difference
       \f$ d_i = summands_with_i_c - substituted_homogeneous_rhs \f$:
       - if \f$ d_i = 0 \f$     -> proceed with the other elements or,
                                   if they are finished, proceed to step 4;
       - if \f$ d_i \neq 0 \f$  -> returns
                                   <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
				   the solution can be wrong or we failed to
				   simplify it.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
     - A new method explained in the paper
       "Checking and Confining the Solutions of Recurrence Relations".
       FIXME: to be written
  -  There are two different way:
     - Verify that \p summands_without_i_c satisfies the recurrence
       (in other words, we are considering all initial conditions equal
       to \f$ 0 \f$).
       Replace \f$ x(n-i) \f$ by \p exact_solution_.expression()
       evaluated at \f$ n-i \f$ (for \f$ i = 1, \cdots, k \f$) in the
       right hand side of the recurrence, and store the result in
       \p substituted_rhs.
       Consider the difference
       \f$ d = summands_without_i_c - substituted_rhs \f$:
       - if \f$ d = 0 \f$     -> returns <CODE>PROVABLY_CORRECT</CODE>:
                                 the solution is certainly right.
       - if \f$ d \neq 0 \f$  -> returns
                                 <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
				 the solution can be wrong or we failed to
				 simplify it.
	 FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.

     - A new method explained in the paper
       "Checking and Confining the Solutions of Recurrence Relations".
       FIXME: to be written

   FIXME: In the latter case, we will need more powerful tools to
   decide whether the solution is right or it is really wrong and, in this
   last case, to return <CODE>PROVABLY_INCORRECT</CODE>.
     
   Case 2: non-linear recurrences of finite order.
   Considers the order of the linear recurrence associated to that one
   non-linear, since the initial conditions are the same.
   Applies the 4 steps of the previous case.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_finite_order(bool partial_verification) const {
  // If `is_symbolic_solution' is true then the solution contains
  // symbolic initial conditions (e.g. x(i), with i non-negative integer);
  // otherwise the solution has been evaluated on the values contained in
  // the map `initial_conditions_'.
  bool is_symbolic_solution = false;
  if (initial_conditions_.empty())
    is_symbolic_solution = true;

  // We will store here the order of the recurrence.
  index_type order_rec;
  // We will store here the greatest common divisor among the
  // decrements `d' of the terms `x(n-d)' occurring in the right
  // hand side of the recurrence.
  unsigned int gcd;
  // The order of the non-linear recurrence is equal to the order of the
  // associated first_order linear recurrence.
  // In order to compute the solution we need the order of the linear
  // recurrence and so we can use it and to avoid to store in a variable
  // the order of the non-linear recurrence.
  // The same thing holds also for the gcd. 
  if (is_non_linear_finite_order()) {
    order_rec = associated_linear_rec().order();
    gcd = associated_linear_rec().gcd_among_decrements();
  }
  else {
    order_rec = order();
    gcd = gcd_among_decrements();
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

  Expr exact_solution;
  // Expand the solution in the case `*this' has been solved with the
  // order reduction method (if we are in the case of solution NOT symbolic
  // then it has already been expanded before substituting the values
  // of the map `initial_conditions_').
  if (is_linear_finite_order_const_coeff() && applied_order_reduction()
      && is_symbolic_solution)
    // FIXME: with `2' means that the solution is never expanded!
    // At the moment we always verify the solution of the "reduced"
    // recurrence.
    if (order_rec < 2)
      exact_solution = write_expanded_solution(*this, gcd_among_decrements());
    else {
      // In this case the expanded solution is too much complicated, then
      // we try to verify the solution of the reduced recurrence.
      // The verification of the reduced recurrence, does not guarantee
      // the correctness of the solution supplied to the user in the case
      // of errors in the computation of the solution of the original
      // recurrence starting from the solution of the reduced recurrence.
      unset_order_reduction();
      Symbol r = insert_auxiliary_definition(mod(n, gcd));
      unsigned int dim = coefficients().size() / gcd + 1;
      std::vector<Expr> new_coefficients(dim);
      Expr inhomogeneous = 0;
      Recurrence rec_rewritten
	(write_reduced_order_recurrence(recurrence_rhs, r, gcd, 
					coefficients(), new_coefficients,
					inhomogeneous));
      rec_rewritten.finite_order_p
	= new Finite_Order_Info(dim - 1, new_coefficients, 1);
      rec_rewritten.set_type(type());
      rec_rewritten.set_inhomogeneous_term(inhomogeneous);
      rec_rewritten.set_first_valid_index(first_valid_index);
      rec_rewritten.solve_linear_finite_order();
      return rec_rewritten.verify_exact_solution(true);
    }
  else
    if (is_symbolic_solution)
      exact_solution = exact_solution_.expression();
    else
      exact_solution = evaluated_exact_solution_.expression();
  D_VAR(exact_solution);

  // Step 1: split the solution in 2 parts: terms with initial conditions
  // are stored in `summands_with_i_c', all the other terms are stored in
  // `summands_without_i_c'.
  // In the case of not symbolic solution there are not symbolic initial
  // conditions: all the solution is stored in `summands_without_i_c'.
  Expr summands_with_i_c = 0;
  Expr summands_without_i_c = 0;
  split_solution(is_symbolic_solution, exact_solution,
		 summands_with_i_c, summands_without_i_c);

  // Prepare a vector containing the coefficients of the symbolic
  // initial conditions occurring in the solution.
  std::vector<Expr> coefficients_i_c(order_rec);
  if (is_linear_finite_order() && is_symbolic_solution)
    fill_vector_coefficients_i_c(summands_with_i_c, gcd, first_valid_index,
				 coefficients_i_c);
  D_VEC(coefficients_i_c, 0, coefficients_i_c.size()-1);
  
  // Step 2: validation of symbolic initial conditions.
  Verify_Status status = validate_initial_conditions(is_symbolic_solution,
						     order_rec,
						     coefficients_i_c,
						     summands_with_i_c,
						     summands_without_i_c);
  if (status != PROVABLY_CORRECT)
    return status;

  // Step 3.
#if NEW_VERIFICATION
  // In order to apply the method explained in the paper
  // "Checking and Confining the Solutions of Recurrence Relations"
  // are necessary the roots of the characteristic equation associated
  // to linear finite order with constant coefficients.
  std::vector<Polynomial_Root> roots;
  if (is_linear_finite_order_const_coeff()) {
    std::vector<Number> num_coefficients(order_rec + 1);
    bool all_distinct = true;
    Expr characteristic_eq;
    if (!characteristic_equation_and_its_roots(order_rec, coefficients(),
					       num_coefficients,
					       characteristic_eq, roots,
					       all_distinct)) {
#if 0
      abort();
#else
      for (unsigned int i = 0; i < order_rec/gcd; ++i) {
	Verify_Status verify_status = traditional_step_3(order_rec,
							 coefficients_i_c[i]);
	if (verify_status != PROVABLY_CORRECT)
	  return verify_status;
      }
      goto continue_with_step_4;
#endif
    }
    else {
      // Step 3: new method.
      bool new_step_3_failed = false;
      for (unsigned int i = 0; i < order_rec/gcd; ++i)
	if (!verify_new_method_const_coeff(order_rec, coefficients_i_c[i],
					   roots, false)) {
	  new_step_3_failed = true;
	  break;
	}
      if (!new_step_3_failed)
	goto continue_with_step_4;
    }
  }
#endif
   
  {
    // The traditional step 3 is applied in the case of non-linear recurrence
    // or if the new method applied on the homogeneous part of the
    // recurrence is failed.
    if (is_linear_finite_order())
      for (unsigned int i = 0; i < order_rec/gcd; ++i) {
	Verify_Status verify_status = traditional_step_3(order_rec,
							 coefficients_i_c[i]);
	if (verify_status != PROVABLY_CORRECT)
	  return verify_status;
      }
    // In the case of non-linear recurrences the homogeneous part
    // of the solution can contains more than one symbolic initial
    // condition: the vector `coefficients_i_c' is in this case empty.
    else {
      assert(is_non_linear_finite_order());
      Verify_Status verify_status = traditional_step_3(order_rec,
						       summands_with_i_c);
      if (verify_status != PROVABLY_CORRECT)
	return verify_status;
    }
  }

 continue_with_step_4:
  // The recurrence is homogeneous.
  if (summands_without_i_c == 0)
    if (partial_verification)
      return PARTIAL_PROVABLY_CORRECT;
    else
      return PROVABLY_CORRECT;
  
#if NEW_VERIFICATION
  // Step 4: the method of the paper
  // "Checking and Confining the Solutions of Recurrence Relations"
  // for linear finite order with constant coefficients.
  if (is_linear_finite_order_const_coeff())
    if (verify_new_method_const_coeff(order_rec, summands_without_i_c, roots,
				      true))
      if (partial_verification)
	return PARTIAL_PROVABLY_CORRECT;
      else
	return PROVABLY_CORRECT;
#endif
  
#if NEW_VERIFICATION
  // FIXME: to be finished!
//    // Step 4: the method of the paper
//    // "Checking and Confining the Solutions of Recurrence Relations"
//    // for linear finite order with variable coefficients.
//    if (is_linear_finite_order_var_coeff())
//      if (verify_new_method_var_coeff(order_rec, summands_without_i_c))
//        return PROVABLY_CORRECT;
#endif

  // Traditional way in order to verify exact solution of `*this'..

  // Step 4: by substitution, verifies that `summands_without_i_c'
  // satisfies the recurrence.
  // Computes `substituted_rhs' by substituting, in the rhs
  // of the recurrence, `n' by `n - d' (where `d' is the decrement
  // of the i-th term `a_i(n)*x(n - d)').
  Expr substituted_rhs = recurrence_rhs;
  // FIXME: the code in the "else" is, from a theoretical point of view,
  // better for efficiency because it performs only 1 symbolic substitution
  // and d numeric substitutions, but in practice is not so.
#if 1
  for (index_type d = 1; d <= order_rec; ++d)
    if (substituted_rhs.has(x(n - d))) {
      Expr shifted_solution
	= simplify_all(summands_without_i_c.substitute(n, n - d));
      shifted_solution = simplify_sum(shifted_solution, REWRITE_UPPER_LIMIT);
      substituted_rhs = substituted_rhs.substitute(x(n - d), shifted_solution);
    }
#else
  Symbol h;
  Expr shifted_solution
    = simplify_all(summands_without_i_c.substitute(n, n - h));
  for (index_type d = 1; d <= order_rec; ++d)
    if (substituted_rhs.has(x(n - d))) {
      Expr shifted_solution_num
	= shifted_solution.substitute(h, d);
      shifted_solution_num = simplify_sum(shifted_solution_num,
					  REWRITE_UPPER_LIMIT);
      substituted_rhs = substituted_rhs.substitute(x(n - d),
						   shifted_solution_num);
    }
#endif  

  Expr diff = summands_without_i_c - substituted_rhs;
  // The expression `diff' can be more difficult to simplify.
  // For this motive we performed simplification also
  // before to expand blackboard's definitions.
  diff = blackboard.rewrite(diff);
  diff = simplify_all(diff);
  if (diff.is_zero())
    if (partial_verification)
      return PARTIAL_PROVABLY_CORRECT;
    else
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
    Symbol r = insert_auxiliary_definition(mod(n, gcd));
    unsigned int dim = coefficients().size() / gcd + 1;
    std::vector<Expr> new_coefficients(dim);
    Expr inhomogeneous = 0;
    Recurrence rec_rewritten
      (write_reduced_order_recurrence(recurrence_rhs, r, gcd, 
				      coefficients(), new_coefficients,
				      inhomogeneous));
    rec_rewritten.finite_order_p
      = new Finite_Order_Info(dim - 1, new_coefficients, 1);
    rec_rewritten.set_type(type());
    rec_rewritten.set_inhomogeneous_term(inhomogeneous);
    rec_rewritten.set_first_valid_index(first_valid_index);
    rec_rewritten.solve_linear_finite_order();
    return rec_rewritten.verify_exact_solution();
  }
  else {
    diff = simplify_all(diff);
    if (diff.is_zero())
      if (partial_verification)
	return PARTIAL_PROVABLY_CORRECT;
      else
	return PROVABLY_CORRECT;
    else
      return INCONCLUSIVE_VERIFICATION;
  }
}

/*!
  Consider the right hand side of a weighted-average recurrence
  (or of a recurrence transformable in a weighted-average recurrence)
  \f[
    x(n) = f(n) \sum_{k=n_0}^{n-1} x(k) + g(n),
  \f]
  where \f$ n_0 \in \Nset \cup \{ 0 \}.
  Let \p recurrence_rhs be the variable containing the right-hand side
  of the recurrence.
  Assume that the system has produced the expression
  \p exact_solution_.expression(), which is a function of the
  variable \f$ n \f$.

  The verification's process is divided in 2 steps:
  -  Validation of the \ref initial_condition
     "symbolic initial condition x(n_0)".
     Evaluate the expression \p exact_solution_.expression() for
     \f$ n = n_0 \f$ simplify as much as possible.
     If the final result is <EM>synctactically</EM> equal to \f$ x(n_0) \f$
     the symbolic initial conditions is verified and we can proceed to
     step 2; otherwise return <CODE>INCONCLUSIVE_VERIFICATION</CODE>
     because the solution can be wrong or it is not enough simplified.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
  -  Validation of the solution using substitution.
     Consider the difference
     \f[
       d = x(n) - \left( f(n) \sum_{k=n_0}^{n-1} x(k) + g(n) \right)
     \f],
     with \f$ x(n) \f$ and \f$ x(k) \f$ replaced with the solution.
     Since the closed formula for \f$ x(n) \f$ is guaranteed to hold for
     \f$ n >= n_0 + 1 \f$ only, the lower limit of the sum must start from
     \f$ n_0 + 1 \f$ and the term for \f$ n = n_0 \f$ must be considered
     before. Hence, the difference that we consider is
     \f[
       d = x(n) - \left( f(n) x(n_0) + f(n) \sum_{k=1}^{n-1} x(k)
                   + g(n) \right).
     \f]
     There are 2 possibilities:
     - if \f$ d = 0 \f$     -> returns <CODE>PROVABLY_CORRECT</CODE>:
                               the solution is certainly right.
     - if \f$ d \neq 0 \f$  -> returns <CODE>INCONCLUSIVE_VERIFICATION</CODE>:
 			       the solution can be wrong or we failed to
			       simplify it.
     FIXME: in some cases it can return <CODE>PROVABLY_INCORRECT</CODE>.
*/
PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_weighted_average() const {
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
  
  if (!initial_conditions_.empty())
    if (initial_conditions_.rbegin()->first > lower)
      // To solve `x(n) = f(n)*sum(k,n_0,n-1,x(k))+g(n)' with the
      // initial condition `x(m) = h', where `m > n_0',
      // is like to solve `x(n) = f(n)*sum(k,m,n-1,x(k))+g(n)'.
      lower = initial_conditions_.rbegin()->first;
  
  // The case `f(n) = -1', i.e. the recurrence has the form
  // `x(n) = - sum(k, n_0, n-1, x(k)) + g(n)', is special:
  // the solution is simply `x(n) = g(n) - g(n-1)'.
  // FIXME: the traditional validation' process does not work,
  // is it true?
  if (weight_rec == -1)
    // FIXME: verify!!!
    return PROVABLY_CORRECT;
  
  // Note: the solution is valid only for `n > lower'.
  Expr exact_solution;
  if (evaluated_exact_solution_.has_expression())
    exact_solution = evaluated_exact_solution_.expression();
  else
    exact_solution = exact_solution_.expression();

  // Step 1: validation of the initial condition.
  Expr e = exact_solution.substitute(n, lower+1);
  e = simplify_all(e);
  if (e != (weight_rec * get_initial_condition(lower)
	    + inhomogeneous).substitute(n, lower+1))
    // FIXME: provably_incorrect...
    return INCONCLUSIVE_VERIFICATION;
  
  // Step 2: validation of the solution.
  // Consider the `sum(k, n_0, n-1, x(k)', with `x(k)' replaced by
  // the solution, and tries to semplify it.
  Symbol h;
  Expr diff = PURRS::sum(h, lower+1, n - 1, exact_solution.substitute(n, h));
  diff = simplify_sum(diff, COMPUTE_SUM);
  // Consider the difference
  // `x(n) - (f(n) x(n_0) + f(n) sum(k, n_0+1, n-1, x(k)) + g(n))'
  // and tries to simplify it.
  diff = exact_solution - weight_rec * get_initial_condition(lower)
    - weight_rec * diff - inhomogeneous;
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
  
  // FIXME: temporary!!!!!!!!
  // With the introduction of the function `Sc()' the validation
  // of the initial condition must change!!!!!
#if 0
  // Step 1: validation of initial conditions.
  Verify_Status status
    = validate_initial_conditions(kind_of_bound, bound,
				  applicability_condition());
  if (status == INCONCLUSIVE_VERIFICATION || status == PROVABLY_INCORRECT)
    return status;
#endif  

  // Step 2: find `partial_bound'.
  // We not consider the terms containing the initial conditions:
  // `partial_bound' will contain all the other terms.
  Expr partial_bound = 0;
  if (bound.is_a_add())
    for (unsigned int i = bound.nops(); i-- > 0; ) {
      if (!has_only_symbolic_initial_conditions(bound.op(i)))
	partial_bound += bound.op(i);
    }
  else
    if (!has_only_symbolic_initial_conditions(bound))
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

PURRS::Recurrence::Verify_Status
PURRS::Recurrence::verify_exact_solution(bool partial_verification) const {
  if (!exact_solution_.has_expression())
    throw std::logic_error("PURRS::Recurrence::verify_exact_solution() "
			   "called, but no exact solution was computed");

  switch (type_) {
  case ORDER_ZERO:
  case LINEAR_FINITE_ORDER_CONST_COEFF:
  case LINEAR_FINITE_ORDER_VAR_COEFF:
    return verify_finite_order(partial_verification);
  case NON_LINEAR_FINITE_ORDER:
    {
      Verify_Status status = verify_finite_order(partial_verification);
      if (status == INCONCLUSIVE_VERIFICATION && !partial_verification) {
	// In the simple case of non-linear recurrence of the form
	// `x(n) = c x(n-1)^a', where `c' and `a' are constants (`a != 1'),
	// we know the solution without computing the associated linear
	// recurrence. Now we want to know the solution of the associated
	// linear recurrence.
	if (coeff_simple_non_linear_rec() != 0)
	  associated_linear_rec().compute_exact_solution_finite_order();
	return associated_linear_rec().verify_exact_solution(true);
      }
      else
	return status;
    }
    break;
  case WEIGHTED_AVERAGE:
    return verify_weighted_average();
  case FUNCTIONAL_EQUATION:
    // FIXME: see again what to make in the case of the
    // functional equations.
    // Special case: functional equation of the form `x(n) = x(n/b)'.
    return INCONCLUSIVE_VERIFICATION;
  default:
    throw std::runtime_error("PURRS internal error: "
			     "verify_exact_solution().");
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
