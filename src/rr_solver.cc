/* Definition of the main recurrence relation solver.
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

#include <config.h>

#define NOISY 0

#include "globals.hh"
#include "util.hh"
#include "sum_poly.hh"
#include "rr_solver.hh"
#include "simplify.hh"
#include "alg_eq_solver.hh"
#include <climits>
#include <algorithm>
#include <string>

// TEMPORARY
#include <iostream>
#include <fstream>

using namespace GiNaC;

/*!
  Returns <CODE>true</CODE> if and only if \p e is of the form
  \f$ n - d \f$ with \f$ d \f$ an integer; in this case the opposite
  of \f$ d \f$ to \p decrement.
*/
static bool
get_constant_decrement(const GExpr& e, const GSymbol& n, GNumber& decrement) {
  static GExpr n_plus_d = n + wild(0);
  GList substitution;
  if (match(e, n_plus_d, substitution)) {
    GExpr d = get_binding(substitution, 0);
    if (is_a<numeric>(d)) {
      GNumber i = ex_to<numeric>(d);
      if (i.is_integer())
	decrement = -i;
      return true;
    }
  }
  return false;
}

/*!
  Let
  \f$ x_n = a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$
  the recurrence relation of order \f$ k \f$ and the coefficients \f$ a_j \f$
  (for \f$ i = 1, \dotsc, k \f$) are stored in the vector \p coefficients.
  This function returns the characteristic equation
  \f$ x^k - ( a_1 * x^{k-1} + \dotsb + a^{k-1} * x + a_k ) \f$.
  Since <CODE>find_roots()</CODE> solve equations only with integer
  coefficients the \p coefficients' elements must be rationals and,
  eventually, if are not integer, builds an other vector,
  \p int_coefficients, with the element of \p coefficients multiplied for
  the lcm among their denominators.
*/
static GExpr
build_characteristic_equation(const GSymbol& x,
			      const std::vector<GNumber>& coefficients) {
  for (unsigned i = coefficients.size(); i-- > 0; )
    if (!coefficients[i].is_rational())
      throw
	"PURRS error: today the algebraic equation solver works\n"
	"only with integer coefficients.\n"
	"Please come back tomorrow.";
  std::vector<GNumber> denominators;
  // Find the least common multiple among the denominators of the
  // rational elements of `coefficients'.
  for (unsigned i = coefficients.size(); i-- > 0; )
    if (!coefficients[i].is_integer())
      denominators.push_back(coefficients[i].denom());
  GExpr p = 0;
  // Build the vector `int_coefficients' with the elements of
  // `coefficients' multiplied for the least common multiple
  // `least_com_mul'.
  // There is no need to `order' because the order of the recurrence relation
  // is equal to `coefficients.size() - 1'.
  if (denominators.size() != 0) {
    GNumber least_com_mul = lcm(denominators);
    std::vector<GNumber> int_coefficients(coefficients);
    // In the first position of `coefficients' there is 0.
    for (unsigned i = coefficients.size(); i-- > 1; )
      int_coefficients[i] *= least_com_mul;
    for (unsigned i = coefficients.size() - 1; i-- > 0; )
      p += pow(x, i) * (-int_coefficients[coefficients.size() - 1 - i]);
    p += least_com_mul * pow(x, coefficients.size() - 1);
  }
  else {
    for (unsigned i = coefficients.size() - 1; i-- > 0; )
      p += pow(x, i) * (-coefficients[coefficients.size() - 1 - i]);
    p += pow(x, coefficients.size() - 1);
  }
  return p;
}

//! Returns <CODE>true</CODE> if at least one element of the vector \p v
//! is different from \f$ 0 \f$.
//! Returns <CODE>false</CODE> otherwise, i. e., all the elements of \p v
//! are equal to \f$ 0 \f$. 
static bool
vector_not_all_zero(const std::vector<GExpr>& v) {
  for (unsigned i = v.size(); i-- > 0; )
    if (!v[i].is_zero())
      return true;
  return false;
}

static GExpr
return_sum(const bool distinct, const GSymbol& n, const GNumber& order,
	   const GExpr& coeff, const GSymbol& alpha, const GSymbol& lambda) {
  GSymbol k("k");
  GSymbol x("x");
  GExpr q_k = coeff.subs(n == k);
  GExpr symbolic_sum;  
  if (distinct)
    symbolic_sum = sum_poly_times_exponentials(q_k, k, n, x);
  else
    symbolic_sum = sum_poly_times_exponentials(q_k, k, n, 1);
  // `sum_poly_times_exponentials' calculates the sum from 0 while
  // we want to start from `order'.
  symbolic_sum -= q_k.subs(k == 0);
  for (int j = 1; j < order; ++j)
    symbolic_sum -= q_k.subs(k == j) * pow(alpha, j) * pow(lambda, -j);
  if (distinct)
    symbolic_sum = symbolic_sum.subs(x == alpha/lambda);
  symbolic_sum *= pow(lambda, n);
  symbolic_sum = simplify_on_output_ex(symbolic_sum.expand(), n, false);

  return symbolic_sum;
}

/*!
  Consider the inhomogeneous term \f$ e(n) \f$ in the form polynomials
  times exponentials:
  \f$ e(n) = \sum_{i=0}^k \alpha_i^{n} p_i(n) \f$ (where \f$ k \f$ is
  the number of exponentials).
  We call \p \lambda the generic root of the characteristic equation
  and \p alpha the generic base of an exponential.

  This function fills the two vectors of <CODE>GExpr</CODE>
  \p symbolic_sum_distinct and \p symbolic_sum_no_distinct,
  with two different sums:
  fixed tha base \p alpha_i, for each root \p lambda 
  - if \f$ \alpha_i \neq \lambda \f$ then
    \f[
      symbolic\_sum\_distinct[i]
        = \lambda^n * f_i(\alpha_i / \lambda)
        = \lambda^n * \sum_{k=order}^n (\alpha_i / \lambda)^k \cdot p_i(k);
    \f]
  - if \f$ \alpha_i = \lambda \f$ then
    \f[
      symbolic\_sum\_no\_distinct[i]
        = \lambda^n * f_i(1)
        = \lambda^n * \sum_{k=order}^n p_i(k).
    \f]
    In the i-th position of \p symbolic_sum_no_distinct, in the first case,
    and of \p symbolic_sum_distinct, in the second case, is put \f$ 0 \f$
    so that the two vector have always the same dimensions.
*/
static void
compute_symbolic_sum(const GSymbol& n,
		     const GSymbol& alpha, const GSymbol& lambda,
		     const std::vector<Polynomial_Root>& roots,
		     const std::vector<GExpr>& base_of_exps,
		     const std::vector<GExpr>& exp_poly_coeff,
		     std::vector<GExpr>& symbolic_sum_distinct,
		     std::vector<GExpr>& symbolic_sum_no_distinct) {
  // Compute the order of the recurrence relation.
  GNumber order = 0;
  for (unsigned i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  unsigned r = 0;
  for (unsigned i = base_of_exps.size(); i-- > 0; )
    for (unsigned j = roots.size(); j-- > 0; ) {
      bool distinct = true;
      if (roots[j].value().is_equal(base_of_exps[i]))
	distinct = false;
      
      // The root is different from the exponential's base.
      if (distinct) {
	if (order <= 2)
	  symbolic_sum_distinct.push_back(return_sum(true, n, order,
						     exp_poly_coeff[i],
						     alpha, lambda));
	else
	  symbolic_sum_distinct.push_back(return_sum(true, n, order,
						     exp_poly_coeff[r],
						     alpha, lambda));
	symbolic_sum_no_distinct.push_back(0);
      }
      // The root is equal to the exponential's base.
      else {
	if (order <= 2)
	  symbolic_sum_no_distinct.push_back(return_sum(false, n, order,
							exp_poly_coeff[i],
							alpha, lambda));
	else
	  symbolic_sum_no_distinct.push_back(return_sum(false, n, order,
							exp_poly_coeff[r],
							alpha, lambda));
	symbolic_sum_distinct.push_back(0);
      }
      ++r;
    }
#if NOISY
  D_VEC(symbolic_sum_distinct, 0, symbolic_sum_distinct.size()-1);
  D_VEC(symbolic_sum_no_distinct, 0, symbolic_sum_no_distinct.size()-1);
#endif
}

/*!
  Consider the vectors \p symbolic_sum_distinct and \p symbolic_sum_no_distinct
  that contain all the symbolic sums of the inhomogeneous term's terms that
  are polynomial or the product of a polynomial and an exponential,
  For each sum this function
  - substitutes to \p lambda the corresponding value of the characteristic
    equation's root;
  - substitutes to \p alpha the corresponding base of the eventual
    exponential.  Returns a <CODE>GExpr</CODE> \p solution with the
    sum of all sums of the vectors.
*/
static GExpr
subs_to_sum_roots_and_bases(const GSymbol& alpha, const GSymbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    const std::vector<GExpr>& base_of_exps,
			    std::vector<GExpr>& symbolic_sum_distinct,
			    std::vector<GExpr>& symbolic_sum_no_distinct) {
  // Compute the order of the recurrence relation.
  GNumber order = 0;
  for (unsigned i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  GExpr solution = 0;
  unsigned r = 0;
  for (unsigned i = base_of_exps.size(); i-- > 0; )
    for (unsigned j = roots.size(); j-- > 0; ) {
      GExpr base_exp = base_of_exps[i];
      GExpr tmp;
      if (!base_exp.is_equal(roots[j].value()))
	tmp = symbolic_sum_distinct[r].subs(lst(alpha == base_exp,
						lambda == roots[j].value()));
      else
	tmp
	  = symbolic_sum_no_distinct[r].subs(lst(alpha == base_exp,
						 lambda == roots[j].value()));
      if (order == 2)
	solution += tmp * pow(-1, j);
      else
	solution += tmp;
      ++r;
    }
  return solution;
}

static void
exp_poly_decomposition(const GExpr& e, const GSymbol& n,
		       std::vector<GExpr>& alpha,
		       std::vector<GExpr>& p,
		       std::vector<GExpr>& q);

static void
add_initial_conditions(const GExpr& g_n, const GSymbol& n,
		       const std::vector<GNumber>& coefficients,
		       const std::vector<GExpr>& initial_conditions,
		       GExpr& solution);

static GExpr
solve_constant_coeff_order_1(const GSymbol& n,
			     const std::vector<GExpr>& base_of_exps,
			     const std::vector<GExpr>& exp_poly_coeff,
			     const std::vector<GExpr>& exp_no_poly_coeff,
			     const std::vector<Polynomial_Root>& roots,
			     const std::vector<GExpr>& initial_conditions);

static GExpr
solve_constant_coeff_order_2(const GSymbol& n, GExpr& g_n, const int order,
			     const bool all_distinct,
			     const std::vector<GExpr>& base_of_exps,
			     const std::vector<GExpr>& exp_poly_coeff,
			     const std::vector<GExpr>& exp_no_poly_coeff,
			     const std::vector<GNumber>& coefficients,
			     const std::vector<Polynomial_Root>& roots);

static GExpr
solve_constant_coeff_order_k(const GSymbol& n, GExpr& g_n, const int order,
			     const bool all_distinct,
			     const std::vector<GExpr>& base_of_exps,
			     const std::vector<GExpr>& exp_poly_coeff,
			     const std::vector<GExpr>& exp_no_poly_coeff,
			     const std::vector<GNumber>& coefficients,
			     const std::vector<Polynomial_Root>& roots);

static bool
verify_solution(const GExpr& solution, const int& order, const GExpr& rhs,
		const GSymbol& n);

/*!
  Returns an expression that is equivalent to \p e and that is
  "maximally expanded" with respect to addition.  This amounts, among
  other things, to distribute multiplication over addition.
*/
GExpr
additive_form(const GExpr& e) {
  return e.expand();
}

/*!
  This function solves recurrences of SOME TYPE provided they are
  supplied in SOME FORM. (Explain.)
*/
Solver_Status
solve(const GExpr& rhs, const GSymbol& n, GExpr& solution) {
  // The following code depends on the possibility of recovering
  // the various parts of `rhs' as summands of an additive expression.
  GExpr e = additive_form(rhs);
  
  // Initialize the computation of the order of the linear part of the
  // recurrence.  This works like the computation of a maximum: it is
  // the maximum `k' such that `rhs = a*x(n-k) + b' where `a' is not
  // syntactically 0. 
  int order = -1;

  // We will store here the coefficients of linear part of the recurrence.
  std::vector<GExpr> coefficients;

  // Will be set to true if at least one element of coefficients is
  // non-constant.
  bool has_non_constant_coefficients = false;

  do {
    // These patterns are used repeatedly for pattern matching.
    // We avoid recreating them over and over again by declaring
    // them static.
    static GExpr x_i = x(wild(0));
    static GExpr x_i_plus_r = x_i + wild(1);
    static GExpr a_times_x_i = wild(1)*x_i;
    static GExpr a_times_x_i_plus_r = a_times_x_i + wild(2);

    // This will hold the substitutions produced by the various match
    // operations.
    GList substitution;
    
    // This will hold the index `i' in contexts of the form `x(i)'.
    GExpr i;
    
    // This will hold the coefficient `a' in contexts of the form
    // `a*x(i)'.
    GExpr a;
    
    // The following matches are attempted starting from the most
    // common, then the second most common and so forth.  The check
    // `if (!i.has(n))' is necessary because otherwise do not accept
    // `x(i)' with `i' numeric in a general recurrence relation
    // (es. x(n) = x(n-1)+x(0)).
    if (clear(substitution), match(e, x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = 1;
      e = get_binding(substitution, 1);
    }
    else if (clear(substitution), match(e, a_times_x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = get_binding(substitution, 1);
      e = get_binding(substitution, 2);
    }
    else if (clear(substitution), match(e, a_times_x_i, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = get_binding(substitution, 1);
      e = 0;
    }
    else if (clear(substitution), match(e, x_i, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = 1;
      e = 0;
    }
    else
      break;

    GNumber decrement;
    if (!get_constant_decrement(i, n, decrement))
      return HAS_NON_INTEGER_DECREMENT;
    if (decrement == 0)
      return HAS_NULL_DECREMENT;
    if (decrement < 0)
      return HAS_NEGATIVE_DECREMENT;
    // Make sure that (1) we can represent `decrement' as a long, and
    // (2) we will be able to store the coefficient into the
    // appropriate position of the `coefficients' vector.
    if (decrement >= LONG_MAX || decrement >= coefficients.max_size())
      return HAS_HUGE_DECREMENT;

    // Detect non-constant coefficients, i.e., those with occurrences of `n'.
    if (a.has(n))
      has_non_constant_coefficients = true;

    // The `order' is defined as the maximum value of `index'.
    unsigned long index = decrement.to_long();
    if (order < 0 || index > unsigned(order))
      order = index;

    // The vector `coefficients' contains in the `i'-th position the
    // coefficient of `x(n-i)'.  The first position always contains 0.
    if (index > coefficients.size())
      coefficients.insert(coefficients.end(),
			  index - coefficients.size(),
			  GNumber(0));
    if (index == coefficients.size())
      coefficients.push_back(a);
    else
      coefficients[index] += a;
  } while (e != 0);
  
  // Special case: `e' is a function of `n', the parameters and of
  // `x(k_1)', ..., `x(k_m)' where `m >= 0' and `k_1', ..., `k_m' are
  // non-negative integers. 
  if (order < 0) {
    solution = e;
    return OK;
  }

  D_VAR(order);
  D_VEC(coefficients, 1, order);
  D_MSGVAR("Inhomogeneous term: ", e);

  // Simplifies expanded expressions, in particular rewrites nested powers.
  e = simplify_on_input_ex(e, n, true);

  // We search exponentials in `n' (for this the expression `e'
  // must be expanded).
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `e'. In the `i'-th position of the vectors
  // `exp_poly_coeff' and `exp_no_poly_coeff' there are respectively
  // the polynomial part and possibly non polynomial part of the coefficient
  // of the exponential with the base in `i'-th position of `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of base_of_exps[i]^n.
  std::vector<GExpr> base_of_exps;
  std::vector<GExpr> exp_poly_coeff;
  std::vector<GExpr> exp_no_poly_coeff;
  exp_poly_decomposition(e, n,
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);

  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  // Create the vector of initial conditions.
  std::vector<GExpr> initial_conditions(order);
  for (int i = 0; i < order; ++i)
    initial_conditions[i] = x(i);

  // Compute the characteristic equation and its roots.
  GExpr characteristic_eq;
  GSymbol y("y");
  std::vector<Polynomial_Root> roots;
  bool all_distinct = true;
  // FIXME: il seguente if sull'ordine e' temporaneo perche' cosi' si
  // riescono a fare le parametriche del primo ordine almeno.
  // Temporaneo fino a che `find_roots()' accettera' i parametri anche
  // per equazioni di grado superiore al primo.
  std::vector<GNumber> num_coefficients(order+1);
  if (order == 1) {
    characteristic_eq = y - coefficients[1];
    roots.push_back(coefficients[1]);
  }
  else {
  // Check if the vector `coefficients' contains only numeric
  // elements and in this case use a vector of GNumber.
    for (unsigned i = coefficients.size(); i--> 0; )
      if (is_a<numeric>(coefficients[i]))
	num_coefficients[i] = ex_to<numeric>(coefficients[i]);
      else
	throw
	  "PURRS error: today the recurrence relation\n"
	  "does not support irrationals coefficients.\n"
	  "Please come back tomorrow.";
    characteristic_eq
      = build_characteristic_equation(y, num_coefficients);
    if (!find_roots(characteristic_eq, y, roots, all_distinct))
      return TOO_COMPLEX;
  }

  D_VAR(characteristic_eq);
  D_VEC(roots, 0, roots.size()-1);
  D_MSG("");

  GExpr g_n;
  switch (order) {
  case 1:
    if (!has_non_constant_coefficients)
      solution = solve_constant_coeff_order_1(n, base_of_exps, exp_poly_coeff,
					      exp_no_poly_coeff, roots,
					      initial_conditions);
    else
      throw
	"PURRS error: today we only solve recurrence"
	"relations with constant coefficients.\n"
	"Please come back tomorrow.";
    break;
  case 2:
    if (!has_non_constant_coefficients)
      solution = solve_constant_coeff_order_2(n, g_n, order, all_distinct,
					      base_of_exps, exp_poly_coeff,
					      exp_no_poly_coeff, 
					      num_coefficients, roots);
    else
      throw
	"PURRS error: today we only solve recurrence"
	"relations with constant coefficients.\n"
	"Please come back tomorrow.";
    break;
  default:
    if (!has_non_constant_coefficients)
      solution = solve_constant_coeff_order_k(n, g_n, order, all_distinct,
					      base_of_exps, exp_poly_coeff,
					      exp_no_poly_coeff,
					      num_coefficients, roots);
    else
      throw
	"PURRS error: today we only solve recurrence"
	"relations with constant coefficients.\n"
	"Please come back tomorrow.";
    break;
  }
  
  if (order > 1)
    add_initial_conditions(g_n, n, num_coefficients, initial_conditions,
			   solution);
  
  D_MSGVAR("Before calling simplify: ", solution);
  solution = simplify_on_output_ex(solution.expand(), n, false);
  
  // Only for the output.
  GList conditions;
  for (unsigned i = order; i-- > 0; )
    conditions.append(initial_conditions[i]);
  solution = solution.collect(conditions);
  
  if (!verify_solution(solution, order, rhs, n))
    D_MSG("ATTENTION: the following solution can be wrong\n"
	  "or not enough simplified.");

  return OK;
}


void
impose_condition(const std::string&) {
}

/*!
  Assuming that \p rhs contains occurrences of <CODE>x(n-k)</CODE>
  where <CODE>k</CODE> is a negative integer, this function tries to
  perform suitable changes of variables that preserve the meaning of
  the recurrence relation, but transforms it into its <EM>standard
  form</EM> <CODE>x(n) = new_rhs</CODE>, where <CODE>new_rhs</CODE>
  does not contain any instance of <CODE>x(n-k)</CODE>, with a
  negative integer <CODE>k</CODE>.
  Returns <CODE>true</CODE> if the transformation is successful;
  returns <CODE>false</CODE> otherwise.
*/
static bool
eliminate_negative_decrements(const GExpr& /* rhs */, GExpr& /* new_rhs */) {
  // Let `j' be the largest positive integer such that `x(n+j)' occurs
  // in `rhs' with a coefficient `a' which is not syntactically 0.
  // Then the changes of variables include replacing `n' by `n-j',
  // changing sign, and division by `a'.

  return false;
}

/*!
  Here we assume that \p rhs contains occurrences of <CODE>x(n)</CODE> itself.
  Therefore the recurrence may be impossible.  This function decides
  if this is the case and, if so, it returns <CODE>false</CODE>.  If the
  recurrence is solvable, it is rewritten into its normal form, which
  is then written in <CODE>new_rhs</CODE>, and the function returns
  <CODE>true</CODE>.
*/
static bool
eliminate_null_decrements(const GExpr& /* rhs */, GExpr& /* new_rhs */) {

  //   Assume that `rhs = a*x(n) + b' and that `b' does not contain
  //   `x(n)'.  The following cases are possible:
  //   - If `a = 1' and `b' does not contain any occurrence of `x(n-k)'
  //      where `k' is a positive integer, the recurrence is impossible.
  //   - If `a = 1' and `b' contains `x(n-k)' for some positive integer `k'
  //     and with a coefficient that is not syntactically 0, we remove
  //     `x(n)' from both sides of `x(n) = rhs', and then rewrite the
  //     recurrence into its standard form.
  //   - If `a != 1' we move `a*x(n)' to the left-hand side, and divide
  //     through by `1 - a', obtaining the standard form, which is 
  //     `(rhs - a*x(n)) / (1-a)'.

  return false;
}

/*!
  This function solves recurrences of SOME TYPE provided they
  are supplied in SOME FORM. (Explain.)
  It does that by repeatedly calling solve() and handling
  the errors that may arise.
*/
Solver_Status
solve_try_hard(const GExpr& rhs, const GSymbol& n, GExpr& solution) {
  bool exit_anyway = false;
  Solver_Status status;
  do {
    status = solve(rhs, n, solution);
    switch (status) {
    case OK:
      break;
    case HAS_NON_INTEGER_DECREMENT:
    case HAS_HUGE_DECREMENT:
      exit_anyway = true;
      break;
    case HAS_NEGATIVE_DECREMENT:
      {
	GExpr new_rhs;
	if (eliminate_negative_decrements(rhs, new_rhs))
	    status = solve(new_rhs, n, solution);
	else
	  exit_anyway = true;
      }
      break;
    case HAS_NULL_DECREMENT:
      {
	GExpr new_rhs;
	if (eliminate_null_decrements(rhs, new_rhs))
	    status = solve(new_rhs, n, solution);
	else
	  exit_anyway = true;
      }
      exit_anyway = true;
      break;
    default:
      throw std::runtime_error("PURRS internal error: "
			       "solve_try_hard.");
      break;
    }
  } while (!exit_anyway && status != OK);
  return status;
}

static void
exp_poly_decomposition_factor(const GExpr& base,
			      const GExpr& e, const GSymbol& n,
			      std::vector<GExpr>& alpha,
			      std::vector<GExpr>& p,
			      std::vector<GExpr>& q) {
  unsigned alpha_size = alpha.size();
  unsigned position;
  bool found = false;
  for (unsigned i = alpha_size; i-- > 0; )
    if (base == alpha[i]) {
      position = i;
      found = true;
      break;
    }
  if (!found) {
    alpha.push_back(base);
    p.push_back(0);
    q.push_back(0);
    position = alpha_size;
  }
  // Here `alpha[position]' contains `base' and the polynomial and
  // possibly not polynomial parts of `e' can be added to
  // `p[position]' and `q[position]', respectively.
  GExpr polynomial;
  GExpr possibly_not_polynomial;
  isolate_polynomial_part(e, n, polynomial, possibly_not_polynomial);
  p[position] += polynomial;
  q[position] += possibly_not_polynomial;
}

static void
exp_poly_decomposition_summand(const GExpr& e, const GSymbol& n,
			       std::vector<GExpr>& alpha,
			       std::vector<GExpr>& p,
			       std::vector<GExpr>& q) {
  static GExpr exponential = pow(wild(0), n);
  GList substitution;
  unsigned num_factors = is_a<mul>(e) ? e.nops() : 1;
  if (num_factors == 1) {
    if (match(e, exponential, substitution)) {
      // We have found something of the form `pow(base, n)'.
      GExpr base = get_binding(substitution, 0);
      assert(base != 0);
      if (is_scalar_representation(base, n)) {
	// We have found something of the form `pow(base, n)'
	// and `base' is good for the decomposition.
	exp_poly_decomposition_factor(base, 1, n, alpha, p, q);
	return;
      }
    }
  }
  else
    for (unsigned i = num_factors; i-- > 0; ) {
      if (clear(substitution), match(e.op(i), exponential, substitution)) {
	// We have found something of the form `pow(base, n)'.
	GExpr base = get_binding(substitution, 0);
	assert(base != 0);
	if (is_scalar_representation(base, n)) {
	  // We have found something of the form `pow(base, n)'
	  // and `base' is good for the decomposition: determine
	  // `r = e/pow(base, n)'.
	  GExpr r = 1;
	  for (unsigned j = num_factors; j-- > 0; )
	    if (i != j)
	      r *= e.op(j);
	  exp_poly_decomposition_factor(base, r, n, alpha, p, q);
	  return;
	}
      }
    }
  // No proper exponential found: this is treated like `pow(1, n)*e'.
  exp_poly_decomposition_factor(1, e, n, alpha, p, q);
}

/*!
  Let \f$ e(n) \f$ be the expression in \p n contained in \p e,
  which is assumed to be already expanded.
  This function computes a decomposition
  \f$ e(n) = \sum_{i=0}^k \alpha_i^n \bigl(p_i(n) + q_i(n)\bigr) \f$, where
  - \f$ \alpha_i \f$ is a expression valid for to be an exponential's base.
    (syntactically different from \p 0);
  - \f$ \alpha_i \neq \alpha_j \f$ if \f$ i \neq j \f$;
  - \f$ p_i(n) \f$ is (syntactically) a polynomial in \f$ n \f$.

  The expressions corresponding to \f$ \alpha_i \f$, \f$ p_i \f$ and
  \f$ q_i \f$ are stored in the \f$ i \f$-th position of the vectors
  \p alpha, \p p and \p q, respectively.
*/
static void
exp_poly_decomposition(const GExpr& e, const GSymbol& n,
		       std::vector<GExpr>& alpha,
		       std::vector<GExpr>& p,
		       std::vector<GExpr>& q) {
  unsigned num_summands = is_a<add>(e) ? e.nops() : 1;
  // An upper bound to the number of exponentials is the number of
  // summands in `e': reserve space in the output vectors so that
  // no reallocations will be required.
  alpha.reserve(num_summands);
  p.reserve(num_summands);
  q.reserve(num_summands);
  if (num_summands > 1)
    for (unsigned i = num_summands; i-- > 0; )
      exp_poly_decomposition_summand(e.op(i), n, alpha, p, q);
  else
    exp_poly_decomposition_summand(e, n, alpha, p, q);
}

/*!
  Adds to the sum already computed those corresponding to the initial
  conditions:
  \f[
    \sum_{i=0}^{order - 1} g_{n-i}
      \bigl( x_i - \sum_{j=1}^i a_j x_{i-j} \bigr).
  \f]
*/
// FIXME: il vettore `coefficients' dovra' diventare di `GExpr' quando
// sapremo risolvere anche le eq. di grado superiore al primo con i
// parametri.
static void
add_initial_conditions(const GExpr& g_n, const GSymbol& n,
                       const std::vector<GNumber>& coefficients,
		       const std::vector<GExpr>& initial_conditions,
		       GExpr& solution) {
  // `coefficients.size()' has `order + 1' elements because in the first
  // position there is the value 0. 
  for (unsigned i = coefficients.size() - 1; i-- > 0; ) {
    GExpr g_n_i = g_n.subs(n == n - i);
    GExpr tmp = initial_conditions[i];
    for (unsigned j = i; j > 0; j--)
      tmp -= coefficients[j] * initial_conditions[i-j];
    solution += tmp * g_n_i;
  }
}

/*!
  Consider the linear recurrence relation of first order with
  constant coefficients
  \f[
    x_n = \lambda x_{n-1} + p(n),
  \f]
  where <CODE>p(n)</CODE> is a function defined over the natural numbers.
  In this case we know the final formula that give the solution (we observe
  that the coefficient coincides with the root of the characteristic
  equation):
  \f[
    x_n = \lambda^n * x_0
          + \sum_{k=1}^n \lambda^{n-k} p(k).
  \f]
*/
static GExpr
solve_constant_coeff_order_1(const GSymbol& n,
			     const std::vector<GExpr>& base_of_exps,
			     const std::vector<GExpr>& exp_poly_coeff,
			     const std::vector<GExpr>& exp_no_poly_coeff,
			     const std::vector<Polynomial_Root>& roots,
			     const std::vector<GExpr>& initial_conditions) {
  GExpr solution;
  // Calculates the solution of the first order recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential.
  if (!vector_not_all_zero(exp_no_poly_coeff)) {
    GSymbol alpha("alpha");
    GSymbol lambda("lambda");
    std::vector<GExpr> symbolic_sum_distinct;
    std::vector<GExpr> symbolic_sum_no_distinct;
    compute_symbolic_sum(n, alpha, lambda, roots,
			 base_of_exps, exp_poly_coeff,
			 symbolic_sum_distinct, symbolic_sum_no_distinct);
    // Substitutes to the sums in the vectors `symbolic_sum_distinct' or
    // `symbolic_sum_no_distinct' the value of the characteristic equation's
    // root and of the bases of the eventual exponentials.
    // In `solution' put the sum of all sums of the vectors after the
    // substitution.
    solution = subs_to_sum_roots_and_bases(alpha, lambda, roots,
					   base_of_exps,
					   symbolic_sum_distinct,
					   symbolic_sum_no_distinct);
    // FIXME: per ora non si puo' usare la funzione
    // `add_initial_conditions' perche' richiede un vettore di
    // `GNumber' come `coefficients' e voglio risolvere anche le
    // parametriche.
    // add_initial_conditions(1, n, coefficients, initial_conditions,
    //		              solution);
    solution += initial_conditions[0] * pow(roots[0].value(), n);
  }
  else
    throw
      "PURRS error: today we only allow inhomogeneous terms\n"
      "in the form of polynomials or product of exponentials\n"
      "and polynomials.\n"
      "Please come back tomorrow.";
  return solution;
}

/*!
  Consider the <EM>fundamental</EM> solution of the associated
  homogeneous equation, which is
  \f[
    \begin{cases}
      g_n = a_1 g_{n-1} + a_2 g_{n-2} + \cdots + a_k g_{n-k}, \\
      g_0 = 1, \\
      g_n = a_1 g_{n-1} + a_2 g_{n-2} + \cdots + a_{n-1} g_1 + a_n g_0
        & \text{for $1 \le n < k$,} \\
    \end{cases}
  \f]
  and the general solution of the homogeneous recurrence <CODE>g_n</CODE>:
  - if the roots are simple, i. e., they are all distinct, then
    \f$ g_n = \alpha_1 \lambda_1^n + \cdots + \alpha_k \lambda_k^n \f$,
  - if there are multiple roots then
    \f[
      g_n = \sum_{j=1}^r (\alpha_{j,0} + \alpha_{j,1}n
            + \cdots + \alpha_{j,\mu_j-1}n^{\mu_j-1}) \lambda_j^n,
    \f]
  where the roots <CODE>lambda</CODE> of the recurrence are knowed.
  This function defines and solves the system in <CODE>k</CODE>
  equations and <CODE>k</CODE> unknowns to found the <CODE>alpha</CODE>
  to insert in order to determine <CODE>g_n</CODE>.
  Returns in the matrix \p solution the solution of the system. 
*/
static GMatrix
solve_system(const bool all_distinct,
	     const std::vector<GNumber>& coefficients,
	     const std::vector<Polynomial_Root>& roots) {
  unsigned coefficients_size = coefficients.size();
  // Prepare a list with the elments for the right hand side of the system
  // to solve.
  // The elements of the list are `g_0, g_1, ..., g_{k-1}'.
  // Note that `tmp[i]' is built on `tmp[i-1]'.
  std::vector<GExpr> tmp(coefficients_size - 1);
  tmp[0] = 1;
  for (unsigned i = 1; i < coefficients_size - 1; ++i)
    for (unsigned j = 0; j < i; ++j)
      tmp[i] += coefficients[j+1] * tmp[i-j-1];
  GList g_i;
  for (unsigned i = coefficients_size - 1; i-- > 0; )
    g_i.prepend(tmp[i]);
  
  // Prepare a list with the coefficients of the equations of the system
  // to solve. This calculus is based on different form of `g_n'
  // in according to the roots' multiplicity.
  GList coeff_equations;
  if (all_distinct)
    for (unsigned i = coefficients_size - 1; i-- > 0; )
      for (unsigned j = coefficients_size - 1; j-- > 0; )
	coeff_equations.prepend(pow(roots[j].value(), i));
  else
    for (unsigned h = 0; h < coefficients_size - 1; ++h)
      for (unsigned i = roots.size(); i-- > 0; ) {
	for (GNumber j = roots[i].multiplicity(); j-- > 1; )
	  coeff_equations.append(pow(h, j) * pow(roots[i].value(), h));
	coeff_equations.append(pow(roots[i].value(), h));
      }

  // Define the matrices and solve the system.
  GMatrix coeff_alpha(coefficients_size - 1, coefficients_size - 1,
		      coeff_equations);
  GMatrix rhs(coefficients_size - 1, 1, g_i);
  GMatrix vars(coefficients_size - 1, 1);
  for (unsigned i = coefficients_size - 1; i-- > 0; )
    vars(i, 0) = GSymbol();
  GMatrix solution = coeff_alpha.solve(vars, rhs);
#if NOISY
  std::cout << "g_i: " << g_i << std::endl;
  std::cout << "equations: " << coeff_equations << std::endl;
  std::cout << "system: " << coeff_alpha << std::endl;
  std::cout << "rhs: " << rhs << std::endl;
  std::cout << "alpha_i: " << sol << std::endl;
#endif
  return solution;
}

static GExpr
find_g_n(const GSymbol& n, const bool all_distinct, const GMatrix sol,
	 const std::vector<Polynomial_Root>& roots) {
  // Compute the order of the recurrence relation.
  GNumber order = 0;
  for (unsigned i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  GExpr g_n = 0;
  if (all_distinct)
    for (int i = 0; i < order; ++i)
      g_n += sol(i, 0) * pow(roots[i].value(), n);
  else
    for (unsigned i = roots.size(); i-- > 0; ) {
      int h = 0;
      for (GNumber j = roots[i].multiplicity(); j-- > 0 && h < order; ) {
	g_n += sol(h, 0) * pow(n, j) * pow(roots[i].value(), n);
	++h;
      }
    }
  return g_n;
}

/*!
  Consider \f$ p(n) \f$ the non homogeneous part of the recurrence relation
  \f$ p(n) = c_1 p_1^n + c_2 p_2^n + \dotsb + c_k p_k^n \f$ and
  \f$ g(n) = d_1 g_1^n + d_2 g_2^n + \dotsb + d_k g_h^n \f$.
  We consider the decompositions of both the previous sums:
  \p bases_of_exp and \p exp_poly_coeff will contain respectively 
  \f$ p_1, p_2, \cdots, p_k \f$ and \f$ c_1, c_2, \dots, c_k \f$;
  \p bases_exp_g_n and \p g_n_poly_coeff instead will contain respectively
  \f$ g_1, g_2, \cdots, g_h \f$ and \f$ d_1, d_2, \cdots, d_h \f$.
  We build the new vector \p poly_coeff_tot multiplting the elements
  of \p exp_poly_coeff and \p g_n_poly_coeff, i. e. it will contain
  \f$ c_1 d_1, c_1 d_2, \cdots, c_1 d_h, c_2 d_1, \cdots, c_k d_h \f$.
  The function <CODE>compute_symbolic_sum()</CODE> will fill the vectors
  \p symbolic_sum_distinct and \p symbolic sum_distinct from the position
  \f$ 0 \f$ and it will have like parameters the elements of \p poly_coeff_tot.
*/
static void
prepare_for_symbolic_sum(const GSymbol& n, const GExpr& g_n,
			 const std::vector<Polynomial_Root>& roots,
			 const std::vector<GExpr>& exp_poly_coeff,
			 std::vector<GExpr>& poly_coeff_tot) {
  std::vector<GExpr> bases_exp_g_n;
  std::vector<GExpr> g_n_poly_coeff;
  std::vector<GExpr> g_n_no_poly_coeff;
  exp_poly_decomposition(g_n, n, bases_exp_g_n,
			 g_n_poly_coeff, g_n_no_poly_coeff);
  // `bases_of_exp_g_n' must have same elements of `roots' in the same order.
  bool equal = true;
  for (unsigned i = roots.size(); i-- > 0; )
    if (!roots[i].value().is_equal(bases_exp_g_n[i]))
      equal = false;
  if (!equal) {
    std::vector<GExpr> tmp_exp(roots.size());
    std::vector<GExpr> tmp_coeff_poly(roots.size());
    std::vector<GExpr> tmp_coeff_no_poly(roots.size());
    for (unsigned i = roots.size(); i-- > 0; )
      tmp_exp[i] = roots[i].value();
    for (unsigned i = tmp_exp.size(); i-- > 0; )
      for (unsigned j = bases_exp_g_n.size(); j-- > 0; )
	if (tmp_exp[i].is_equal(bases_exp_g_n[j])) {
	  tmp_coeff_poly[i] = g_n_poly_coeff[j];
	  tmp_coeff_no_poly[i] = g_n_no_poly_coeff[j];
	}
    copy(tmp_exp.begin(), tmp_exp.end(), bases_exp_g_n.begin());
    copy(tmp_coeff_poly.begin(), tmp_coeff_poly.end(), g_n_poly_coeff.begin());
    copy(tmp_coeff_no_poly.begin(), tmp_coeff_no_poly.end(),
	 g_n_no_poly_coeff.begin());
  }
  // The roots are simple, i. e., their multiplicity is 1.
  for (unsigned i = exp_poly_coeff.size(); i-- > 0; )
    for (unsigned j = g_n_poly_coeff.size(); j-- > 0; )
      poly_coeff_tot.push_back(exp_poly_coeff[i] * g_n_poly_coeff[j]);
#if NOISY
  D_VEC(bases_exp_g_n, 0, bases_exp_g_n.size()-1);
  D_VEC(g_n_poly_coeff, 0, g_n_poly_coeff.size()-1);
  D_VEC(g_n_no_poly_coeff, 0, g_n_no_poly_coeff.size()-1);
  D_VEC(poly_coeff_tot, 0, poly_coeff_tot.size()-1);
#endif
}

static GExpr
compute_non_homogeneous_part(const GSymbol& n, const GExpr& g_n,
			     const int order,
			     const std::vector<GExpr>& base_of_exps,
			     const std::vector<GExpr>& exp_poly_coeff) {
  GExpr solution_tot = 0;
  std::vector<GExpr> bases_exp_g_n;
  std::vector<GExpr> g_n_poly_coeff;
  std::vector<GExpr> g_n_no_poly_coeff;
  exp_poly_decomposition(g_n, n, bases_exp_g_n,
			 g_n_poly_coeff, g_n_no_poly_coeff);
#if NOISY
  D_VEC(bases_exp_g_n, 0, bases_exp_g_n.size()-1);
  D_VEC(g_n_poly_coeff, 0, g_n_poly_coeff.size()-1);
  D_VEC(g_n_no_poly_coeff, 0, g_n_no_poly_coeff.size()-1);
#endif  
  for (unsigned i = bases_exp_g_n.size(); i-- > 0; )
    for (unsigned j = base_of_exps.size(); j-- > 0; ) {
      GExpr solution = 0;
      GSymbol k("k");
      GExpr g_n_coeff_k = g_n_poly_coeff[i].subs(n == n - k);
      GExpr exp_poly_coeff_k = exp_poly_coeff[j].subs(n == k);
      solution = sum_poly_times_exponentials(g_n_coeff_k * exp_poly_coeff_k,
					     k, n, 1/bases_exp_g_n[i]
					     * base_of_exps[j]);
      // `sum_poly_times_exponentials' calculates the sum from 0 while
      // we want to start from `order'.
      solution -= (g_n_coeff_k * exp_poly_coeff_k).subs(k == 0);
      for (int h = 1; h < order; ++h)
	solution -= (g_n_coeff_k * exp_poly_coeff_k).subs(k == h)
	  * pow(1/bases_exp_g_n[i] * base_of_exps[j], h);
      solution *= pow(bases_exp_g_n[i], n);
      solution_tot += solution;
    }
  return solution_tot;
}

/*!
  Consider the linear recurrence relation of the second order with
  constant coefficients
  \f[
    x_n = a_1 x_{n-1} + a_2 + p(n),
  \f]
  where <CODE>p(n)</CODE> is a function defined over the natural numbers.
  If the roots of the characteristic equation \f$ \lambda_1 \f$ and
  \f$ \lambda_2 \f$ are distinct then we know the final formula that
  give the solution:
  \f[
    x_n = g_{n-1} x_1 + a_2 g_{n-2} x_0
          + \frac{\lambda_1^{n+1}}{\lambda_1-\lambda_2}
	    \sum_{k=2}^n \lambda_1^{-k} p(k)
          - \frac{\lambda_2^{n+1}}{\lambda_1-\lambda_2}
	    \sum_{k=2}^n \lambda_2^{-k} p(k),
  \f]
  where
  \f[
    g(n) = \frac{\lambda_1^{n+1}-\lambda_2^{n+1}}{\lambda_1-\lambda_2}.
  \f]
  If there is one root, \f$ \lambda \f$, with double multiplicity,
  computes \f$ g_n = (\alpha_1 + \alpha_2) * \lambda^n \f$,
  where \f$ \alpha_1 \f$ and \f$ \alpha_2 \f$ are complex numbers.
  Introduced the <EM>fundamental</EM> solution of the associated
  homogeneous equation, which is
  \f[
    \begin{cases}
      g_n = a_1 g_{n-1} + a_2 g_{n-2},\\
      g_0 = 1, \\
      g_1 = a_1 g_0,
    \end{cases}
  \f]
  this function returns the general solution of recurrence relation
  which is calculated by the formula
  \f[
    x_n = \sum_{i=k}^n g_{n-i} p(i)
          + \sum_{i=0}^{k-1} g_{n-i}
	    \Bigl( x_i - \sum_{j=1}^i a_j x_{i-j} \Bigr).
  \f]
  The two sums in the previous formula correspond to the non-homogeneous
  part \f$ p(n) \f$ and to the initial conditions (computed afterwards by
  the function <CODE>add_initial_conditions()</CODE>), respectively.
*/
static GExpr
solve_constant_coeff_order_2(const GSymbol& n, GExpr& g_n, const int order,
			     const bool all_distinct, 
			     const std::vector<GExpr>& base_of_exps,
			     const std::vector<GExpr>& exp_poly_coeff,
			     const std::vector<GExpr>& exp_no_poly_coeff,
			     const std::vector<GNumber>& coefficients,
			     const std::vector<Polynomial_Root>& roots) {
  GExpr solution;
  // Calculates the solution of the second order recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential.
  if (!vector_not_all_zero(exp_no_poly_coeff))
    if (all_distinct) {
      GExpr root_1 = roots[0].value();
      GExpr root_2 = roots[1].value();
      GExpr diff_roots = root_1 - root_2;
      GSymbol alpha("alpha");
      GSymbol lambda("lambda");
      std::vector<GExpr> symbolic_sum_distinct;
      std::vector<GExpr> symbolic_sum_no_distinct;
      compute_symbolic_sum(n, alpha, lambda, roots,
			   base_of_exps, exp_poly_coeff,
			   symbolic_sum_distinct, symbolic_sum_no_distinct);
      for (unsigned j = symbolic_sum_distinct.size(); j-- > 0; ) {
	symbolic_sum_no_distinct[j] *= lambda / diff_roots;
	symbolic_sum_distinct[j] *= lambda / diff_roots;
      }
      // Substitutes to the sums in the vector `symbolic_sum_distinct'
      // or `symbolic_sum_no_distinct' the corresponding values of the
      // characteristic equation's roots and of the bases of the
      // eventual exponentials and in `solution' put the sum of all
      // sums of the vector after the substitution.
      solution = subs_to_sum_roots_and_bases(alpha, lambda, roots,
					     base_of_exps,
					     symbolic_sum_distinct,
					     symbolic_sum_no_distinct);
      g_n = (pow(root_1, n+1) - pow(root_2, n+1)) / diff_roots;
      // FIXME: forse conviene semplificare g_n
      D_VAR(g_n);
    }
    else {
      // The characteristic equation
      // x^2 + a_1 * x + a_2 = 0 has a double root.
      assert(roots[0].multiplicity() == 2);      
      
      // Solve system in order to finds `alpha_i' (i = 1,...,order).
      GMatrix sol = solve_system(all_distinct, coefficients, roots);
      
      // Finds `g_n', always taking into account the root's multiplicity
      g_n = find_g_n(n, all_distinct, sol, roots);
      D_VAR(g_n);
      
      solution = compute_non_homogeneous_part(n, g_n, order, base_of_exps,
					      exp_poly_coeff);
    }
  else
    throw
      "PURRS error: today we only allow inhomogeneous terms\n"
      "in the form of polynomials or product of exponentials\n"
      "and polynomials.\n"
      "Please come back tomorrow.";
  return solution;
 }

/*!
  Consider the linear recurrence relation of order <CODE>k</CODE> with
  constant coefficients
  \f[
    x_n = a_1 x_{n-1} + a_2 x_{n-2} + \cdots + a_k x_{n-k} + p(n),
  \f]
  where <CODE>p(n)</CODE> is a function defined over the natural numbers.
  Knowing the roots \f$ \lambda_1, \cdots, \lambda_k \f$ of the
  characteristic equation, builds the general solution of the homogeneous
  recurrence <CODE>g_n</CODE>:
  - if the roots are simple, i. e., they are all distinct, then
    \f$ g_n = \alpha_1 \lambda_1^n + \cdots + \alpha_k \lambda_k^n \f$,
  - if there are multiple roots then
    \f[
      g_n = \sum_{j=1}^r (\alpha_{j,0} + \alpha_{j,1}n
            + \cdots + \alpha_{j,\mu_j-1}n^{\mu_j-1}) \lambda_j^n,
    \f]
  where \f$ \alpha_1, \cdots, \alpha_k \f$ are complex numbers
  (<CODE>g_n</CODE> in the fisrt case is contained in those of the second
  case as special case).
  Introduced the <EM>fundamental</EM> solution of the associated
  homogeneous equation, which is
  \f[
    \begin{cases}
      g_n = a_1 g_{n-1} + a_2 g_{n-2} + \cdots + a_k g_{n-k}, \\
      g_0 = 1, \\
      g_n = a_1 g_{n-1} + a_2 g_{n-2} + \cdots + a_{n-1} g_1 + a_n g_0
        & \text{for $1 \le n < k$,} \\
     \end{cases}
  \f]
  this function returns the general solution of recurrence relation
  which is calculated by the formula
  \f[
    x_n = \sum_{i=k}^n g_{n-i} p(i)
          + \sum_{i=0}^{k-1} g_{n-i}
	    \Bigl( x_i - \sum_{j=1}^i a_j x_{i-j} \Bigr).
  \f]
  The two sums in the previous formula correspond to the non-homogeneous
  part \f$ p(n) \f$ and to the initial conditions (computed afterwards by
  the function <CODE>add_initial_conditions()</CODE>), respectively.
*/
static GExpr
solve_constant_coeff_order_k(const GSymbol& n, GExpr& g_n,
			     const int order, const bool all_distinct,
			     const std::vector<GExpr>& base_of_exps,
			     const std::vector<GExpr>& exp_poly_coeff,
			     const std::vector<GExpr>& exp_no_poly_coeff,
			     const std::vector<GNumber>& coefficients,
			     const std::vector<Polynomial_Root>& roots) {
  GExpr solution;  
  // Calculates the solution of the recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential.
  if (!vector_not_all_zero(exp_no_poly_coeff)) {
    // Solve system in order to finds `alpha_i' (i = 1,...,order).
    GMatrix sol = solve_system(all_distinct, coefficients, roots);
    
    // Finds `g_n', always taking into account the root's multiplicity
    g_n = find_g_n(n, all_distinct, sol, roots);
    D_VAR(g_n);
    
    if (all_distinct) {      
      // Prepare for to compute the symbolic sum.
      std::vector<GExpr> poly_coeff_tot;
      prepare_for_symbolic_sum(n, g_n, roots, exp_poly_coeff, poly_coeff_tot);
      GSymbol alpha("alpha");
      GSymbol lambda("lambda");
      std::vector<GExpr> symbolic_sum_distinct;
      std::vector<GExpr> symbolic_sum_no_distinct;
      compute_symbolic_sum(n, alpha, lambda, roots,
			   base_of_exps, poly_coeff_tot,
			   symbolic_sum_distinct, symbolic_sum_no_distinct);
      // Substitutes to the sums in the vector `symbolic_sum_distinct'
      // or `symbolic_sum_no_distinct' the corresponding values of the
      // characteristic equation's roots and of the bases of the
      // eventual exponentials and in `solution' put the sum of all
      // sums of the vector after the substitution.
      solution = subs_to_sum_roots_and_bases(alpha, lambda, roots,
					     base_of_exps,
					     symbolic_sum_distinct,
					     symbolic_sum_no_distinct);
    }
    else
      // There are roots with multiplicity greater than 1.
      solution = compute_non_homogeneous_part(n, g_n, order, base_of_exps,
					      exp_poly_coeff);
  }
  else
    throw
      "PURRS error: today we only allow inhomogeneous terms\n"
      "in the form of polynomials or product of exponentials\n"
      "and polynomials.\n"
      "Please come back tomorrow.";
  return solution;
}

static void
print_bad_exp(const GExpr& e, const GExpr rhs, bool conditions) {
  std::ofstream outfile("not_verified.out", std::ios_base::app);
  if (conditions)
    outfile << std::endl << "not verified initial conditions in x(n) = "
	    << rhs << std::endl;
  else
    outfile << std::endl << "diff not zero in x(n) = " << rhs << std::endl;
  outfile << e << std::endl;
}

/*!
  Consider the right hand side \p rhs of the order \f$ k \f$ recurrence
  relation
  \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$.
  We try to check that the solution is correct.
  - Validation of initial conditions.
    If \p rhs is equal to \f$ x(0), \cdots, x(k) \f$ for
    \f$ n = 0, \cdots, k-1 \f$ respectively then
    the initial conditions are verified and we continue to check; otherwise
    return false because the solution can be wrong or it is not
    simplified enough.
  - Since the initial conditions are verified, we erase from \p solution
    all terms containing an initial condition.
    In other words, we check that ther remainder of the solution
    satisfies the same recurrence relation, but with the initial conditions
    all equal to \f$ 0 \f$.
    Starting from the partial solution just computed, we substitute
    \f$ x(n-1) \f$, \f$ x(n-2) \f$, \f$ \dots \f$, \f$ x_{n-k} \f$ into \p rhs.
    We next consider the difference \p diff between the partial solution
    and the new right hand side:
    - if \p diff is equal to zero     -> return <CODE>true</CODE>:
                                         the solution is certainly right.
    - if \p diff is not equal to zero (in a syntactical sense)
                                      -> return <CODE>false</CODE>:   
                                         the solution can be wrong or
                                         we failed to simplify it.

  FIXME: In the latter case, we will need more powerful tools to
  decide whether the solution is right or it is really wrong.
*/
static bool
verify_solution(const GExpr& solution, const int& order, const GExpr& rhs,
		const GSymbol& n) {
  // Validation of initial conditions.
  for (int i = order; i-- > 0; ) {
    GExpr g_i = x(i);
    GExpr sol_subs = simplify_numer_denom(solution.subs(n == i));
    if (!g_i.is_equal(sol_subs)) {
      print_bad_exp(sol_subs, rhs, true);
      return false;
    }
  }
  // The initial conditions are verified. Build an other expression
  // that has all terms of `solution' minus those containing an initial
  // condition.
  GExpr partial_solution = 0;
  for (unsigned i = solution.nops(); i-- > 0; )
    if (!solution.op(i).match(x(wild(0)))
	&& !solution.op(i).match(wild(1) * x(wild(0))))
      partial_solution += solution.op(i);

  std::vector<GExpr> terms_to_sub(order);
  for (int i = 0; i < order; ++i)
    terms_to_sub[i] = partial_solution.subs(n == n - i - 1);
  GExpr substituted_rhs = simplify_on_input_ex(rhs.expand(), n, true);
  for (unsigned i = terms_to_sub.size(); i-- > 0; )
    substituted_rhs = substituted_rhs.subs(x(n - i - 1) == terms_to_sub[i]);
  GExpr diff = (partial_solution - substituted_rhs).expand();
  diff = simplify_numer_denom(diff);
  if (diff.is_zero())
    return true;
  else {
    print_bad_exp(diff, rhs, false);
    return false;
  }
}
