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

#ifndef NOISY
#define NOISY 0
#endif

#include "rr_solver.hh"
#include "gosper.hh"
#include "alg_eq_solver.hh"
#include "simplify.hh"
#include "sum_poly.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include "Expr_List.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Matrix.defs.hh"

#include <climits>
#include <algorithm>
#include <string>

// TEMPORARY
#include <iostream>
#include <fstream>

namespace Parma_Recurrence_Relation_Solver {

/*!
  Returns <CODE>true</CODE> in two cases:
  - if \p e is of the form \f$ n - d \f$ with \f$ d \f$ an integer:
    in this case assign the opposite of \f$ d \f$ to \p decrement.
  - if \p e is equal to \p n: in this case decrement is zero.
*/
static bool
get_constant_decrement(const Expr& e, const Symbol& n, Number& decrement) {
  static Expr n_plus_d = n + wild(0);
  Expr_List substitution;
  if (e.match(n_plus_d, substitution)) {
    Expr d = get_binding(substitution, 0);
    if (d.is_a_number()) {
      Number i = d.ex_to_number();
      if (i.is_integer())
	decrement = -i;
      return true;
    }    
  }
  else if (e == n) {
    decrement = 0;
     return true;
  }
  return false;
}

/*!
  Let
  \f$ x_n = a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$
  be a recurrence relation of order \f$ k \f$ whose coefficients \f$ a_j \f$
  (for \f$ i = 1, \dotsc, k \f$) are stored in the vector \p coefficients.
  This function returns the left hand side of the characteristic equation
  \f$ x^k - ( a_1 * x^{k-1} + \dotsb + a^{k-1} * x + a_k ) \f$.
  Since <CODE>find_roots()</CODE> only solves equations with integer
  coefficients, the elements in the vector \p coefficients must be 
  rational numbers.
  Tf they are not integers, this function builds another vector,
  \p int_coefficients, with the element of \p coefficients multiplied by
  the least common multiple of their denominators.
*/
static Expr
build_characteristic_equation(const Symbol& x,
			      const std::vector<Number>& coefficients) {
  for (unsigned i = coefficients.size(); i-- > 0; )
    if (!coefficients[i].is_rational())
      throw
	"PURRS error: today the algebraic equation solver works\n"
	"only with integer coefficients.\n"
	"Please come back tomorrow.";
  std::vector<Number> denominators;
  // Find the least common multiple of the denominators of the
  // rational elements of `coefficients'.
  for (unsigned i = coefficients.size(); i-- > 0; )
    if (!coefficients[i].is_integer())
      denominators.push_back(coefficients[i].denominator());
  Expr p = 0;
  // Build the vector `int_coefficients' with the elements of
  // `coefficients' multiplied by the least common multiple
  // `least_com_mul'.
  // There is no need to know the `order' because the order of the recurrence relation
  // is equal to `coefficients.size() - 1'.
  if (denominators.size() != 0) {
    Number least_com_mul = lcm(denominators);
    std::vector<Number> int_coefficients(coefficients);
    // In the first position of `coefficients' there is 0.
    for (unsigned i = coefficients.size(); i-- > 1; )
      int_coefficients[i] *= least_com_mul;
    for (unsigned i = coefficients.size() - 1; i-- > 0; )
      p += pwr(x, i) * (-int_coefficients[coefficients.size() - 1 - i]);
    p += least_com_mul * pwr(x, coefficients.size() - 1);
  }
  else {
    for (unsigned i = coefficients.size() - 1; i-- > 0; )
      p += pwr(x, i) * (-coefficients[coefficients.size() - 1 - i]);
    p += pwr(x, coefficients.size() - 1);
  }
#if NOISY
  std::cout << "characteristic equation = " << p << std::endl;
#endif
  return p;
}

/*
  Builds the characteristic equation and computes its roots.
  The boolean \p all_distinct is true if all roots are distinct, i.e., every
  roots has multiplicity equal to one, false otherwise.
  Returns <CODE>false<CODE> if it does not succeed to find roots of the
  characteristic equations; returns <CODE>true<CODE> otherwise. 
*/
static bool
characteristic_equation_and_its_roots(int order,
				      const std::vector<Expr>& coefficients,
				      std::vector<Number>& num_coefficients,
				      Expr& characteristic_eq,
				      std::vector<Polynomial_Root>& roots,
				      bool& all_distinct) {
  Symbol y("y");
  // FIXME: il seguente if sull'ordine e' temporaneo perche' cosi' si
  // riescono a fare le parametriche del primo ordine almeno.
  // Temporaneo fino a che `find_roots()' accettera' i parametri anche
  // per equazioni di grado superiore al primo.
  if (order == 1) {
    characteristic_eq = y - coefficients[1];
    roots.push_back(coefficients[1]);
  }
  else {
    // Check if the vector `coefficients' contains only numeric
    // elements and in this case use a vector of Number.
    for (unsigned i = coefficients.size(); i--> 0; )
      if (coefficients[i].is_a_number())
	num_coefficients[i] = coefficients[i].ex_to_number();
      else
	throw
	  "PURRS error: today the recurrence relation of order\n"
	  "greater than one, does not support parametric coefficients.\n"
	  "Please come back tomorrow.";
    characteristic_eq = build_characteristic_equation(y, num_coefficients);
    if (!find_roots(characteristic_eq, y, roots, all_distinct))
      return false;
  }
#if NOISY
  D_VEC(roots, 0, roots.size()-1);
  D_MSG("");
#endif
  return true;
}

//! \brief
//! Returns <CODE>true</CODE> if at least one element of the vector \p v
//! is different from \f$ 0 \f$.
//! Returns <CODE>false</CODE> otherwise, i. e., all the elements of \p v
//! are equal to \f$ 0 \f$. 
static bool
vector_not_all_zero(const std::vector<Expr>& v) {
  for (unsigned i = v.size(); i-- > 0; )
    if (!v[i].is_zero())
      return true;
  return false;
}

static Expr
return_sum(bool distinct, const Symbol& n, const Number& order,
	   const Expr& coeff, const Symbol& alpha, const Symbol& lambda) {
  Symbol k("k");
  Symbol x("x");
  Expr q_k = coeff.subs(n, k);
  Expr symbolic_sum;  
  if (distinct)
    symbolic_sum = sum_poly_times_exponentials(q_k, k, n, x);
  else
    symbolic_sum = sum_poly_times_exponentials(q_k, k, n, 1);
  // `sum_poly_times_exponentials' computes the sum from 0, while
  // we want it to start from `order'.
  symbolic_sum -= q_k.subs(k, 0);
  for (Number j = 1; j < order; ++j)
    symbolic_sum -= q_k.subs(k, j) * pwr(alpha, j) * pwr(lambda, -j);
  if (distinct)
    symbolic_sum = symbolic_sum.subs(x, alpha/lambda);
  symbolic_sum *= pwr(lambda, n);
  symbolic_sum = simplify_on_output_ex(symbolic_sum.expand(), n, false);

  return symbolic_sum;
}

/*!
  Consider an inhomogeneous term \f$ e(n) \f$ which is a sum of products of 
  polynomials and exponentials:
  \f$ e(n) = \sum_{i=0}^k \alpha_i^{n} p_i(n) \f$ (where \f$ k \f$ is
  the number of exponentials, and we assume that the \f$ \alpha_j \f$'s
  are distinct complex numbers).
  We let \f$ \lambda \f$ denote the generic root of the characteristic 
  equation and \$f \alpha \f$ the generic base of an exponential.

  This function fills the two vectors of <CODE>Expr</CODE>
  \p symbolic_sum_distinct and \p symbolic_sum_no_distinct,
  with two different sums: first fix the base \f$ \alpha_i \f$.
  For each root \f$ \lambda \f$
  - if \f$ \alpha_i \neq \lambda \f$ then
    \f[
      symbolic{\_}sum{\_}distinct[i]
        = \lambda^n * f_i(\alpha_i / \lambda)
        = \lambda^n * \sum_{k=order}^n (\alpha_i / \lambda)^k \cdot p_i(k);
    \f]
  - if \f$ \alpha_i = \lambda \f$ then
    \f[
      symbolic{\_}sum{\_}no{\_}distinct[i]
        = \lambda^n * f_i(1)
        = \lambda^n * \sum_{k=order}^n p_i(k).
    \f]
    We put \f$ 0 \f$ in the i-th position of \p symbolic_sum_no_distinct,
    in the first case, and of \p symbolic_sum_distinct, in the second case,
    so that the two vectors have always the same dimensions.
*/
static void
compute_symbolic_sum(const Symbol& n,
		     const Symbol& alpha, const Symbol& lambda,
		     const std::vector<Polynomial_Root>& roots,
		     const std::vector<Expr>& base_of_exps,
		     const std::vector<Expr>& exp_poly_coeff,
		     std::vector<Expr>& symbolic_sum_distinct,
		     std::vector<Expr>& symbolic_sum_no_distinct) {
  // Compute the order of the recurrence relation.
  Number order = 0;
  for (unsigned i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  unsigned r = 0;
  for (unsigned i = base_of_exps.size(); i-- > 0; )
    for (unsigned j = roots.size(); j-- > 0; ) {
      bool distinct = true;
      if (roots[j].value() == base_of_exps[i])
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
}

/*!
  Consider the vectors \p symbolic_sum_distinct and \p symbolic_sum_no_distinct
  that contain all the symbolic sums of the inhomogeneous term's terms that
  are polynomial or the product of a polynomial and an exponential,
  For each sum this function
  - substitutes to \p lambda the corresponding value of the characteristic
    equation's root;
  - substitutes to \p alpha the corresponding base of the eventual
    exponential.  Returns a <CODE>Expr</CODE> \p solution with the
    sum of all sums of the vectors.
*/
static Expr
subs_to_sum_roots_and_bases(const Symbol& alpha, const Symbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    const std::vector<Expr>& base_of_exps,
			    std::vector<Expr>& symbolic_sum_distinct,
			    std::vector<Expr>& symbolic_sum_no_distinct) {
  // Compute the order of the recurrence relation.
  Number order = 0;
  for (unsigned i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  Expr solution = 0;
  unsigned r = 0;
  for (unsigned i = base_of_exps.size(); i-- > 0; )
    for (unsigned j = roots.size(); j-- > 0; ) {
      Expr base_exp = base_of_exps[i];
      Expr tmp;
      if (base_exp != roots[j].value())
	tmp = symbolic_sum_distinct[r]
	  .subs(Expr_List(alpha, lambda),
		Expr_List(base_exp, roots[j].value()));
      else
	tmp
	  = symbolic_sum_no_distinct[r]
	  .subs(Expr_List(alpha, lambda),
		Expr_List(base_exp, roots[j].value()));
      if (order == 2 && (j & 1) == 1)
	solution -= tmp;
      else
	// order != 2 or j even.
	solution += tmp;
      ++r;
    }
  return solution;
}

static void
exp_poly_decomposition(const Expr& e, const Symbol& n,
		       std::vector<Expr>& alpha,
		       std::vector<Expr>& p,
		       std::vector<Expr>& q);

static void
add_initial_conditions(const Expr& g_n, const Symbol& n,
		       const std::vector<Number>& coefficients,
		       const std::vector<Expr>& initial_conditions,
		       Expr& solution);

static Expr
solve_constant_coeff_order_1(const Symbol& n, int order,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff,
			     const std::vector<Expr>& exp_no_poly_coeff,
			     const std::vector<Polynomial_Root>& roots,
			     const std::vector<Expr>& initial_conditions);

static Expr
solve_constant_coeff_order_2(const Symbol& n, Expr& g_n, int order,
			     bool all_distinct,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff,
			     const std::vector<Expr>& exp_no_poly_coeff,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots);

static Expr
solve_constant_coeff_order_k(const Symbol& n, Expr& g_n, int order,
			     bool all_distinct,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff,
			     const std::vector<Expr>& exp_no_poly_coeff,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots);

static Expr
solve_variable_coeff_order_1(const Symbol& n, const Expr& p_n,
			     const Expr& coefficient);

static Expr
compute_alpha_factorial(const Expr& alpha, const Symbol& n);

static bool
verify_solution(const Expr& solution, int order, const Expr& rhs,
		const Symbol& n);

/*!
  Returns an expression that is equivalent to \p e and that is
  "maximally expanded" with respect to addition.  This amounts, among
  other things, to distribute multiplication over addition.
*/
Expr
additive_form(const Expr& e) {
  return e.expand();
}

/*!
  This function solves recurrences of SOME TYPE provided they are
  supplied in SOME FORM. (Explain.)
*/
Solver_Status
solve(const Expr& rhs, const Symbol& n, Expr& solution) {
#if NOISY
  D_VAR(rhs);
#endif
  // The following code depends on the possibility of recovering
  // the various parts of `rhs' as summands of an additive expression.
  Expr e = additive_form(rhs);
  
  // Initialize the computation of the order of the linear part of the
  // recurrence.  This works like the computation of a maximum: it is
  // the maximum `k' such that `rhs = a*x(n-k) + b' where `a' is not
  // syntactically 0. 
  int order = -1;

  // We will store here the coefficients of linear part of the recurrence.
  std::vector<Expr> coefficients;

  // Will be set to true if at least one element of coefficients is
  // non-constant.
  bool has_non_constant_coefficients = false;

  do {
    // These patterns are used repeatedly for pattern matching.
    // We avoid recreating them over and over again by declaring
    // them static.
    static Expr x_i = x(wild(0));
    static Expr x_i_plus_r = x_i + wild(1);
    static Expr a_times_x_i = wild(1)*x_i;
    static Expr a_times_x_i_plus_r = a_times_x_i + wild(2);

    // This will hold the substitutions produced by the various match
    // operations.
    Expr_List substitution;
    
    // This will hold the index `i' in contexts of the form `x(i)'.
    Expr i;
    
    // This will hold the coefficient `a' in contexts of the form
    // `a*x(i)'.
    Expr a;
    
    // The following matches are attempted starting from the most
    // common, then the second most common and so forth.  The check
    // `if (!i.has(n))' is necessary because otherwise do not accept
    // `x(i)' with `i' numeric in a general recurrence relation
    // (es. x(n) = x(n-1)+x(0)).
    if (clear(substitution), e.match(x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = 1;
      e = get_binding(substitution, 1);
    }
    else if (clear(substitution), e.match(a_times_x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = get_binding(substitution, 1);
      e = get_binding(substitution, 2);
    }
    else if (clear(substitution), e.match(a_times_x_i, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = get_binding(substitution, 1);
      e = 0;
    }
    else if (clear(substitution), e.match(x_i, substitution)) {
      i = get_binding(substitution, 0);
      if (!i.has(n))
	break;
      a = 1;
      e = 0;
    }
    else
      break;

    Number decrement;
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
			  Number(0));
    if (index == coefficients.size())
      coefficients.push_back(a);
    else
      coefficients[index] += a;
  } while (!e.is_zero());
  
  // `order' is negative in two cases. 
  if (order < 0)
    if (e.has(x(wild(0))))
      // 1. the recurrence is not linear, i.e. there is a non-linear term
      // containing `x(n-k)', where `k' non-negative-integer.
      return NOT_LINEAR_RECURRENCE;
    else { 
      // 2. `e' is a function of `n', the parameters and of
      // `x(k_1)', ..., `x(k_m)' where `m >= 0' and `k_1', ..., `k_m' are
      // non-negative integers. 
      solution = e;
      return OK;
    }
#if NOISY
  D_VAR(order);
  D_VEC(coefficients, 1, order);
  D_MSGVAR("Inhomogeneous term: ", e);
#endif
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
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(e, n,
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
#if NOISY
  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);
#endif
  // FIXME: the initial conditions can not start always from 0:
  // make a function for this check.
  // Create the vector of initial conditions.
  std::vector<Expr> initial_conditions(order);
  for (int i = 0; i < order; ++i)
    initial_conditions[i] = x(i);

  // `num_coefficients' and `g_n' are defined here because they are
  // necessary in the function `add_initial_conditions()' (at the end
  // of function `solve()').
  std::vector<Number> num_coefficients(order + 1);
  Expr g_n;
  switch (order) {
  case 1:
    if (!has_non_constant_coefficients) {
      Expr characteristic_eq;
      std::vector<Polynomial_Root> roots;
      bool all_distinct = true;
      if (!characteristic_equation_and_its_roots(order, coefficients,
						 num_coefficients,
						 characteristic_eq, roots,
						 all_distinct))
	return TOO_COMPLEX;
      solution = solve_constant_coeff_order_1(n, order, base_of_exps,
					      exp_poly_coeff,
					      exp_no_poly_coeff, roots,
					      initial_conditions);
    }
    else {
#if 1
      throw
	"PURRS error: today we only solve recurrence"
	"relations with constant coefficients.\n"
	"Please come back tomorrow.";
#else
      solution = solve_variable_coeff_order_1(n, e, coefficients[1]);
#endif
    }
    break;
  case 2:
    if (!has_non_constant_coefficients) {
      Expr characteristic_eq;
      std::vector<Polynomial_Root> roots;
      bool all_distinct = true;
      if (!characteristic_equation_and_its_roots(order, coefficients,
						 num_coefficients,
						 characteristic_eq, roots,
						 all_distinct))
	return TOO_COMPLEX;
      solution = solve_constant_coeff_order_2(n, g_n, order, all_distinct,
					      base_of_exps, exp_poly_coeff,
					      exp_no_poly_coeff, 
					      num_coefficients, roots);
    }
    else
      throw
	"PURRS error: today we only solve second order "
	"recurrence relations with constant coefficients.\n"
	"Please come back tomorrow.";
    break;
  default:
    if (!has_non_constant_coefficients) {
      Expr characteristic_eq;
      std::vector<Polynomial_Root> roots;
      bool all_distinct = true;
      if (!characteristic_equation_and_its_roots(order, coefficients,
						 num_coefficients,
						 characteristic_eq, roots,
						 all_distinct))
	return TOO_COMPLEX;
      solution = solve_constant_coeff_order_k(n, g_n, order, all_distinct,
					      base_of_exps, exp_poly_coeff,
					      exp_no_poly_coeff,
					      num_coefficients, roots);
    }
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
#if NOISY
  D_MSGVAR("Before calling simplify: ", solution);
#endif
  solution = simplify_on_output_ex(solution.expand(), n, false);
  
  // Only for the output.
  Expr_List conditions;
  for (unsigned i = order; i-- > 0; )
    conditions.append(initial_conditions[i]);
  solution = solution.collect(conditions);

#if 0  
  if (!verify_solution(solution, order, rhs, n)) {
    std::cout << "x(n) = " << rhs << std::endl;
    std::cout << " -> solution wrong or not enough simplified." << std::endl;
    std::cout << std::endl;
  }
#endif

  return OK;
}

void
impose_condition(const std::string&) {
}

/*!
  If \p possibly_dec is greater than \p max_decrement then \p possibly_dec
  becomes the new maximum decrement and it is assigned to \p max_decrement,
  in this case \p possibly_coeff becomes the new \p coefficient.
*/
void
assign_max_decrement_and_coeff(const Expr& possibly_dec,
			       const Expr& possibly_coeff, const Symbol& n,
			       int& max_decrement, Expr& coefficient) {
  Number decrement;
  get_constant_decrement(possibly_dec, n, decrement);
  int dec = -decrement.to_int();
  if (dec > max_decrement) {
    max_decrement = dec;
    coefficient = possibly_coeff;
  }
}

/*!
  Let \p e be the right hand side of a linear recurrence.
  This functions seeks the largest positive integer \f$ j \f$ such that
  \f$ x(n+j) \f$ occurs in \p e and its coefficient.
  These two values are stored in \p max_decrement and \p coefficient,
  respectively.
*/
void
find_max_decrement_and_coeff(const Expr& e, const Symbol& n,
			     const Expr& x_i, const Expr& a_times_x_i,
			     int& max_decrement, Expr& coefficient) {
  Expr_List substitution;
  if (e.is_a_add()) {
    for (unsigned j = e.nops(); j-- > 0; )
      if (clear(substitution), e.op(j).match(a_times_x_i, substitution))
	assign_max_decrement_and_coeff(get_binding(substitution, 0),
				       get_binding(substitution, 1),
				       n, max_decrement, coefficient);
      else if (clear(substitution), e.op(j).match(x_i, substitution))
	assign_max_decrement_and_coeff(get_binding(substitution, 0), 1,
				       n, max_decrement, coefficient);
  }
  else if (e.is_a_mul()) {
    if (clear(substitution), e.match(a_times_x_i, substitution))
      assign_max_decrement_and_coeff(get_binding(substitution, 0),
				     get_binding(substitution, 1),
				     n, max_decrement, coefficient);
  }
  else
    if (clear(substitution), e.match(x_i, substitution))
      assign_max_decrement_and_coeff(get_binding(substitution, 0), 1,
				     n, max_decrement, coefficient);
}

/*!
  Assuming that \p rhs contains occurrences of \f$ x(n-k) \f$
  where \f$ k \f$ is a negative integer, this function
  performs suitable changes of variables that preserve the meaning of
  the recurrence relation, but transforms it into its <EM>standard
  form</EM> \f$ x(n) = new_rhs \f$, where \f$ new_rhs \f$
  does not contain any instance of \f$ x(n-k) \f$, with a
  negative integer \f$ k \f$.
*/
static void
eliminate_negative_decrements(const Expr& rhs, Expr& new_rhs,
			      const Symbol& n) {
  // Seeks `max_decrement', i.e., the largest positive integer `j' such that
  // `x(n+j)' occurs in `rhs' with a coefficient `coefficient' which is not
  // syntactically 0.
  Expr x_i = x(wild(0));
  Expr a_times_x_i = x_i * wild(1);
  int max_decrement = INT_MIN;
  Expr coefficient;
  find_max_decrement_and_coeff(rhs, n, x_i, a_times_x_i,
			       max_decrement, coefficient);
  // The changes of variables includes replacing `n' by `n-max_decrement',
  // changing sign, and division by `coefficient'.
  new_rhs = rhs.subs(n, n-max_decrement);
  new_rhs *= -1;
  new_rhs = new_rhs.subs(x(n), - x(n-max_decrement)
			 * pwr(coefficient, -1));
  new_rhs /= coefficient;
}

/*!
  Here we assume that \p rhs contains occurrences of \f$ x(n) \f$ itself.
  Therefore the recurrence may be impossible.  This function decides
  if this is the case and, if so, it returns <CODE>false</CODE>.  If the
  recurrence is solvable, it is rewritten into its normal form, which
  is then written in \f$ new_rhs \f$, and the function returns
  <CODE>true</CODE>.
*/
static bool
eliminate_null_decrements(const Expr& rhs, Expr& new_rhs,
			  const Symbol& n) {
  Expr_List substitution;
  // Let `rhs = a*x(n) + b' and that `b' does different to zero
  // and does not contain `x(n)'.  The following cases are possible:
  // 1. If `a = 1' and `b' does not contain any occurrence of `x(n-k)'
  //    where `k' is a positive integer, the recurrence is impossible.
  // 2. If `a = 1' and `b' contains `x(n-k)' for some positive integer `k'
  //    and with a coefficient that is not syntactically 0, we remove
  //    `x(n)' from both sides of `x(n) = rhs', and then rewrite the
  //    recurrence into its standard form.
  // 3. If `a != 1' we move `a*x(n)' to the left-hand side, and divide
  //    through by `1 - a', obtaining the standard form, which is 
  //    `(rhs - a*x(n)) / (1-a)'.
  if (rhs.is_a_add()) {
    // Tries `a' and `b'.
    Expr a;
    Expr b = rhs;
    for (unsigned j = rhs.nops(); j-- > 0; )
      if (clear(substitution), rhs.op(j).match(x(n)*wild(0), substitution)) {
	a = get_binding(substitution, 0);
	b -= rhs.op(j);
      }
      else if (clear(substitution), rhs.op(j).match(x(n), substitution)) {
	a = 1;
	b -= rhs.op(j);
      }
    Expr x_i = x(wild(0));
    Expr a_times_x_i = x_i * wild(1);
    if (a == 1)
      // Case 1.
      if (!b.has(x_i) && !b.has(a_times_x_i))
	return false;
    // Case 2.
      else {
	// Seeks `max_decrement', i.e., the largest integer `j' (it may be
	// non positive) such that `x(n+j)' occurs in `b' with a coefficient
	// `coefficient' which is not syntactically 0.
	int max_decrement = INT_MIN;
	Expr coefficient;
	find_max_decrement_and_coeff(b, n, x_i, a_times_x_i,
				     max_decrement, coefficient);
	// Rewrites the recurrence into its standard form:
	// removes from `b' the term that will be the right hand side of the
	// recurrence, i.e. `x(n+max_decrement)'; changes variable replacing
	// `n+max_decrement' by `n', changes sign and divides for the
	// coefficient of `x(n+max_decrement)'.
	new_rhs = b - coefficient * x(n+max_decrement);
	new_rhs = new_rhs.subs(n, n-max_decrement);
	new_rhs *= -1;
	new_rhs /= coefficient;
      }
    else {
      // Case 3.
      new_rhs = b * pwr(1 - a, -1);
    }
  }
  // Let `rhs = a*x(n)'.
  else if (rhs.is_a_mul()) {
    if (clear(substitution), rhs.match(x(n)*wild(1), substitution))
      new_rhs = 0;
  }
  // Let `rhs = x(n)'.
  else if (rhs == x(n))
    return false;
  
  return true;
}

/*!
  This function solves recurrences of SOME TYPE provided they
  are supplied in SOME FORM. (Explain.)
  It does that by repeatedly calling solve() and handling
  the errors that may arise.
*/
Solver_Status
solve_try_hard(const Expr& rhs, const Symbol& n, Expr& solution) {
  bool exit_anyway = false;
  Solver_Status status;
  do {
    status = solve(rhs, n, solution);
    switch (status) {
    case OK:
      break;
    case HAS_NON_INTEGER_DECREMENT:
    case HAS_HUGE_DECREMENT:
    case TOO_COMPLEX:
      exit_anyway = true;
      break;
    case HAS_NEGATIVE_DECREMENT:
      {
	Expr new_rhs;
	eliminate_negative_decrements(rhs, new_rhs, n);
	//std::cout << "Recurrence tranformed x(n) = " << new_rhs << std::endl;
	status = solve(new_rhs, n, solution);
      }
      break;
    case HAS_NULL_DECREMENT:
      {
	Expr new_rhs;
	if (eliminate_null_decrements(rhs, new_rhs, n)) {
	  //std::cout << "Recurrence tranformed x(n) = " << new_rhs << std::endl;
	  status = solve(new_rhs, n, solution);
	}
	else
	  throw std::logic_error("PURRS error: "
				 "impossible recurrence");
      }
      exit_anyway = true;
      break;
    case NOT_LINEAR_RECURRENCE:
      // Next we will solve also non-linear recurrence and we will
      // call an opportune function.
      throw
	"PURRS error: today only solve linear recurrences. \n"
	"Please come back tomorrow.";
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
exp_poly_decomposition_factor(const Expr& base,
			      const Expr& e, const Symbol& n,
			      std::vector<Expr>& alpha,
			      std::vector<Expr>& p,
			      std::vector<Expr>& q) {
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
  Expr polynomial;
  Expr possibly_not_polynomial;
  isolate_polynomial_part(e, n, polynomial, possibly_not_polynomial);
  p[position] += polynomial;
  q[position] += possibly_not_polynomial;
}

static void
exp_poly_decomposition_summand(const Expr& e, const Symbol& n,
			       std::vector<Expr>& alpha,
			       std::vector<Expr>& p,
			       std::vector<Expr>& q) {
  static Expr exponential = pwr(wild(0), n);
  Expr_List substitution;
  unsigned num_factors = e.is_a_mul() ? e.nops() : 1;
  if (num_factors == 1) {
    if (e.match(exponential, substitution)) {
      // We have found something of the form `power(base, n)'.
      Expr base = get_binding(substitution, 0);
      assert(!base.is_zero());
      if (is_scalar_representation(base, n)) {
	// We have found something of the form `power(base, n)'
	// and `base' is good for the decomposition.
	exp_poly_decomposition_factor(base, 1, n, alpha, p, q);
	return;
      }
    }
  }
  else
    for (unsigned i = num_factors; i-- > 0; ) {
      if (clear(substitution), e.op(i).match(exponential, substitution)) {
	// We have found something of the form `power(base, n)'.
	Expr base = get_binding(substitution, 0);
	assert(!base.is_zero());
	if (is_scalar_representation(base, n)) {
	  // We have found something of the form `power(base, n)'
	  // and `base' is good for the decomposition: determine
	  // `r = e/power(base, n)'.
	  Expr r = 1;
	  for (unsigned j = num_factors; j-- > 0; )
	    if (i != j)
	      r *= e.op(j);
	  exp_poly_decomposition_factor(base, r, n, alpha, p, q);
	  return;
	}
      }
    }
  // No proper exponential found: this is treated like `power(1, n)*e'.
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
exp_poly_decomposition(const Expr& e, const Symbol& n,
		       std::vector<Expr>& alpha,
		       std::vector<Expr>& p,
		       std::vector<Expr>& q) {
  unsigned num_summands = e.is_a_add() ? e.nops() : 1;
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
// FIXME: il vettore `coefficients' dovra' diventare di `Expr' quando
// sapremo risolvere anche le eq. di grado superiore al primo con i
// parametri.
static void
add_initial_conditions(const Expr& g_n, const Symbol& n,
                       const std::vector<Number>& coefficients,
		       const std::vector<Expr>& initial_conditions,
		       Expr& solution) {
  // `coefficients.size()' has `order + 1' elements because in the first
  // position there is the value 0. 
  for (unsigned i = coefficients.size() - 1; i-- > 0; ) {
    Expr g_n_i = g_n.subs(n, n - i);
    Expr tmp = initial_conditions[i];
    for (unsigned j = i; j > 0; j--)
      tmp -= coefficients[j] * initial_conditions[i-j];
    solution += tmp * g_n_i;
  }
}

/*!
  Applies the Gosper's algorithm to express in closed form, if it is
  possible, sum with the summand an hypergeometric term not polynomials
  or polynomials times exponentials (these last sums are computed by the
  functions <CODE>compute_symbolic_sum()</CODE> and
  <CODE>subs_to_sum_roots_and_bases()</CODE>).  This function returns
  <CODE>true</CODE> if the summand is an hypergeometric term,
  independently if it is possible or not to express the sum in closed form.
  Returns <CODE>false</CODE> otherwise.
*/
static bool
gosper_algorithm(int order, const Symbol& n,
		 const std::vector<Expr>& base_of_exps,
		 const std::vector<Expr>& exp_no_poly_coeff,
		 const std::vector<Polynomial_Root>& roots,
       Expr& solution) {
  solution = 0;
  for (unsigned i = exp_no_poly_coeff.size(); i-- > 0; ) {
    Expr tmp;
    if (!exp_no_poly_coeff[i].is_zero()) {
      // FIXME: the lower bound for the sum is not `order'
      // ex. \sum 1/(n^2-1) must be start from 2.
      // FIXME: this is a temporary assert untile he generalization of
      // this function.
      assert(order == 1);
      Expr t = pwr(base_of_exps[i], n) * exp_no_poly_coeff[i]
	* pwr(roots[0].value(), -n);
      //std::cout << "t(n) = " << t << std::endl;
      if (!gosper(t, n, Number(order), n, tmp))
	return false;
    }
    solution += tmp;
  }
  return true;
}

/*!
  Consider the linear recurrence relation of first order with
  constant coefficients
  \f[
    x_n = \lambda x_{n-1} + p(n),
  \f]
  where \f$ p(n) \f$ is a function defined over the natural numbers.
  In this case we know the final formula that give the solution (we observe
  that the coefficient coincides with the root of the characteristic
  equation):
  \f[
    x_n = \lambda^n * x_0
          + \sum_{k=1}^n \lambda^{n-k} p(k).
  \f]
*/
static Expr
solve_constant_coeff_order_1(const Symbol& n, int order,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff,
			     const std::vector<Expr>& exp_no_poly_coeff,
			     const std::vector<Polynomial_Root>& roots,
			     const std::vector<Expr>& initial_conditions) {
  Expr solution = 0;
  // Computes the sum when `\lambda^{n-k} p(k)' is a polynomial or
  // a product of a polynomial times an exponential.
  if (vector_not_all_zero(exp_poly_coeff)) {
    Symbol alpha("alpha");
    Symbol lambda("lambda");
    std::vector<Expr> symbolic_sum_distinct;
    std::vector<Expr> symbolic_sum_no_distinct;
    compute_symbolic_sum(n, alpha, lambda, roots,
			 base_of_exps, exp_poly_coeff,
			 symbolic_sum_distinct, symbolic_sum_no_distinct);
    // Substitutes to the sums in the vectors `symbolic_sum_distinct' or
    // `symbolic_sum_no_distinct' the value of the characteristic equation's
    // root and of the bases of the eventual exponentials.
    // In `solution' put the sum of all sums of the vectors after the
    // substitution.
    solution += subs_to_sum_roots_and_bases(alpha, lambda, roots,
					    base_of_exps,
					    symbolic_sum_distinct,
					    symbolic_sum_no_distinct);
  }
  // Computes the sum when `\lambda^{n-k} p(k)' is not a polynomial or
  // a product of a polynomial times an exponential.
  // The summand must be an hypergeometric term
  if (vector_not_all_zero(exp_no_poly_coeff)) {
    Expr gosper_solution;
    if (gosper_algorithm(order, n, base_of_exps, exp_no_poly_coeff, roots,
			 gosper_solution))
      solution += gosper_solution;
    else {
      // FIXME: the summand is not hypergeometric:
      // no chance of using Gosper's algorithm.
    }
  }
  // FIXME: per ora non si puo' usare la funzione
  // `add_initial_conditions' perche' richiede un vettore di
  // `Number' come `coefficients' e voglio risolvere anche le
  // parametriche (g_n pu' essere posta uguale ad 1 in questo caso).
  // add_initial_conditions(g_n, n, coefficients, initial_conditions,
  //		              solution);
  solution += initial_conditions[0] * pwr(roots[0].value(), n);

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
  and the general solution of the homogeneous recurrence \f$ g_n \f$:
  - if the roots are simple, i. e., they are all distinct, then
    \f$ g_n = \alpha_1 \lambda_1^n + \cdots + \alpha_k \lambda_k^n \f$,
  - if there are multiple roots then
    \f[
      g_n = \sum_{j=1}^r (\alpha_{j,0} + \alpha_{j,1}n
            + \cdots + \alpha_{j,\mu_j-1}n^{\mu_j-1}) \lambda_j^n,
    \f]
  where the roots \f$ lambda \f$ of the recurrence are knowed.
  This function defines and solves the system in \f$ k \f$
  equations and \f$ k \f$ unknowns to found the \f$ alpha \f$
  to insert in order to determine \f$ g_n \f$.
  Returns in the matrix \p solution the solution of the system. 
*/
static Matrix
solve_system(bool all_distinct,
	     const std::vector<Number>& coefficients,
	     const std::vector<Polynomial_Root>& roots) {
  unsigned coefficients_size = coefficients.size();
  // Prepare a list with the elments for the right hand side of the system
  // to solve.
  // The elements of the list are `g_0, g_1, ..., g_{k-1}'.
  // Note that `tmp[i]' is built on `tmp[i-1]'.
  std::vector<Expr> tmp(coefficients_size - 1);
  tmp[0] = 1;
  for (unsigned i = 1; i < coefficients_size - 1; ++i)
    for (unsigned j = 0; j < i; ++j)
      tmp[i] += coefficients[j+1] * tmp[i-j-1];
  Expr_List g_i;
  for (unsigned i = coefficients_size - 1; i-- > 0; )
    g_i.prepend(tmp[i]);
  
  // Prepare a list with the coefficients of the equations of the system
  // to solve. This calculus is based on different form of `g_n'
  // in according to the roots' multiplicity.
  Expr_List coeff_equations;
  if (all_distinct)
    for (unsigned i = coefficients_size - 1; i-- > 0; )
      for (unsigned j = coefficients_size - 1; j-- > 0; )
	coeff_equations.prepend(pwr(roots[j].value(), i));
  else
    for (unsigned h = 0; h < coefficients_size - 1; ++h)
      for (unsigned i = roots.size(); i-- > 0; ) {
	for (Number j = roots[i].multiplicity(); j-- > 1; )
	  coeff_equations.append(pwr(h, j) * pwr(roots[i].value(), h));
	coeff_equations.append(pwr(roots[i].value(), h));
      }

  // Define the matrices and solve the system.
  Matrix coeff_alpha(coefficients_size - 1, coefficients_size - 1,
		      coeff_equations);
  Matrix rhs(coefficients_size - 1, 1, g_i);
  Matrix vars(coefficients_size - 1, 1);
  for (unsigned i = coefficients_size - 1; i-- > 0; )
    vars(i, 0) = Symbol();
  Matrix solution = coeff_alpha.solve(vars, rhs);

  return solution;
}

static Expr
find_g_n(const Symbol& n, bool all_distinct, const Matrix sol,
	 const std::vector<Polynomial_Root>& roots) {
  // Compute the order of the recurrence relation.
  Number order = 0;
  for (unsigned i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  Expr g_n = 0;
  if (all_distinct)
    for (unsigned i = 0; i < order; ++i)
      g_n += sol(i, 0) * pwr(roots[i].value(), n);
  else
    for (unsigned i = roots.size(); i-- > 0; ) {
      unsigned h = 0;
      for (Number j = roots[i].multiplicity(); j-- > 0 && h < order; ) {
	g_n += sol(h, 0) * pwr(n, j) * pwr(roots[i].value(), n);
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
prepare_for_symbolic_sum(const Symbol& n, const Expr& g_n,
			 const std::vector<Polynomial_Root>& roots,
			 const std::vector<Expr>& exp_poly_coeff,
			 std::vector<Expr>& poly_coeff_tot) {
  std::vector<Expr> bases_exp_g_n;
  std::vector<Expr> g_n_poly_coeff;
  std::vector<Expr> g_n_no_poly_coeff;
  exp_poly_decomposition(g_n, n, bases_exp_g_n,
			 g_n_poly_coeff, g_n_no_poly_coeff);
  // `bases_of_exp_g_n' must have same elements of `roots' in the same order.
  bool equal = true;
  for (unsigned i = roots.size(); i-- > 0; )
    if (roots[i].value() != bases_exp_g_n[i])
      equal = false;
  if (!equal) {
    std::vector<Expr> tmp_exp(roots.size());
    std::vector<Expr> tmp_coeff_poly(roots.size());
    std::vector<Expr> tmp_coeff_no_poly(roots.size());
    for (unsigned i = roots.size(); i-- > 0; )
      tmp_exp[i] = roots[i].value();
    for (unsigned i = tmp_exp.size(); i-- > 0; )
      for (unsigned j = bases_exp_g_n.size(); j-- > 0; )
	if (tmp_exp[i] == bases_exp_g_n[j]) {
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
}

static Expr
compute_non_homogeneous_part(const Symbol& n, const Expr& g_n, int order,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff) {
  Expr solution_tot = 0;
  std::vector<Expr> bases_exp_g_n;
  std::vector<Expr> g_n_poly_coeff;
  std::vector<Expr> g_n_no_poly_coeff;
  exp_poly_decomposition(g_n, n, bases_exp_g_n,
			 g_n_poly_coeff, g_n_no_poly_coeff);
  for (unsigned i = bases_exp_g_n.size(); i-- > 0; )
    for (unsigned j = base_of_exps.size(); j-- > 0; ) {
      Expr solution = 0;
      Symbol k("k");
      Expr g_n_coeff_k = g_n_poly_coeff[i].subs(n, n - k);
      Expr exp_poly_coeff_k = exp_poly_coeff[j].subs(n, k);
      solution = sum_poly_times_exponentials(g_n_coeff_k * exp_poly_coeff_k,
					     k, n, 1/bases_exp_g_n[i]
					     * base_of_exps[j]);
      // `sum_poly_times_exponentials' calculates the sum from 0 while
      // we want to start from `order'.
      solution -= (g_n_coeff_k * exp_poly_coeff_k).subs(k, 0);
      for (int h = 1; h < order; ++h)
	solution -= (g_n_coeff_k * exp_poly_coeff_k).subs(k, h)
	  * pwr(1/bases_exp_g_n[i] * base_of_exps[j], h);
      solution *= pwr(bases_exp_g_n[i], n);
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
  where \f$ p(n) \f$ is a function defined over the natural numbers.
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
static Expr
solve_constant_coeff_order_2(const Symbol& n, Expr& g_n, int order,
			     bool all_distinct, 
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff,
			     const std::vector<Expr>& exp_no_poly_coeff,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots) {
  Expr solution;
  // Calculates the solution of the second order recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential.
  if (!vector_not_all_zero(exp_no_poly_coeff))
    if (all_distinct) {
      Expr root_1 = roots[0].value();
      Expr root_2 = roots[1].value();
      Expr diff_roots = root_1 - root_2;
      Symbol alpha("alpha");
      Symbol lambda("lambda");
      std::vector<Expr> symbolic_sum_distinct;
      std::vector<Expr> symbolic_sum_no_distinct;
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
      g_n = (pwr(root_1, n+1) - pwr(root_2, n+1)) / diff_roots;
      // FIXME: forse conviene semplificare g_n
#if NOISY
      D_VAR(g_n);
#endif
    }
    else {
      // The characteristic equation
      // x^2 + a_1 * x + a_2 = 0 has a double root.
      assert(roots[0].multiplicity() == 2);      
      
      // Solve system in order to finds `alpha_i' (i = 1,...,order).
      Matrix sol = solve_system(all_distinct, coefficients, roots);
      
      // Finds `g_n', always taking into account the root's multiplicity
      g_n = find_g_n(n, all_distinct, sol, roots);
#if NOISY
      D_VAR(g_n);
#endif
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
  Consider the linear recurrence relation of order \f$ k \f$ with
  constant coefficients
  \f[
    x_n = a_1 x_{n-1} + a_2 x_{n-2} + \cdots + a_k x_{n-k} + p(n),
  \f]
  where \f$ p(n) \f$ is a function defined over the natural numbers.
  Knowing the roots \f$ \lambda_1, \cdots, \lambda_k \f$ of the
  characteristic equation, builds the general solution of the homogeneous
  recurrence \f$ g_n \f$:
  - if the roots are simple, i. e., they are all distinct, then
    \f$ g_n = \alpha_1 \lambda_1^n + \cdots + \alpha_k \lambda_k^n \f$,
  - if there are multiple roots then
    \f[
      g_n = \sum_{j=1}^r (\alpha_{j,0} + \alpha_{j,1}n
            + \cdots + \alpha_{j,\mu_j-1}n^{\mu_j-1}) \lambda_j^n,
    \f]
  where \f$ \alpha_1, \cdots, \alpha_k \f$ are complex numbers
  (\f$ g_n \f$ in the fisrt case is contained in those of the second
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
static Expr
solve_constant_coeff_order_k(const Symbol& n, Expr& g_n,
			     int order, bool all_distinct,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff,
			     const std::vector<Expr>& exp_no_poly_coeff,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots) {
  Expr solution;  
  // Calculates the solution of the recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential.
  if (!vector_not_all_zero(exp_no_poly_coeff)) {
    // Solve system in order to finds `alpha_i' (i = 1,...,order).
    Matrix sol = solve_system(all_distinct, coefficients, roots);
    
    // Finds `g_n', always taking into account the root's multiplicity
    g_n = find_g_n(n, all_distinct, sol, roots);
#if NOISY
    D_VAR(g_n);
#endif
    if (all_distinct) {      
      // Prepare for to compute the symbolic sum.
      std::vector<Expr> poly_coeff_tot;
      prepare_for_symbolic_sum(n, g_n, roots, exp_poly_coeff, poly_coeff_tot);
      Symbol alpha("alpha");
      Symbol lambda("lambda");
      std::vector<Expr> symbolic_sum_distinct;
      std::vector<Expr> symbolic_sum_no_distinct;
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

//! \brief
//! When possible, computes \f$ alpha!(n) \f$ if \f$ alpha \f$ is a sum
//! or terms, otherwise returns the symbolic product.
static Expr
alpha_factorial_if_add(const Expr& alpha, const Symbol& n) {
  Expr alpha_factorial;
  Expr_List substitution;
  bool alpha_factorial_computed = false;
  if (alpha.match(n + wild(0), substitution)) {
    Expr tmp = get_binding(substitution, 0);
    if (tmp.is_a_number()) {
      Number num = tmp.ex_to_number();
      if (num.is_positive_integer()) {
	alpha_factorial = factorial(alpha) / factorial(num);
	alpha_factorial_computed = true;
      }
    }
  }
  else if (alpha == 2*n+1) {
    alpha_factorial = factorial(2*n+1) / (pwr(2, n) * factorial(n));
    alpha_factorial_computed = true;
  }
  else {
    // Allows to compute `alpha!(n)' for function as `a*n+a*b'
    // (`a' not rational).
    Expr a = alpha.content(n);
    if (a != 1) {
      alpha_factorial = compute_alpha_factorial(alpha.primpart(n), n)
	* compute_alpha_factorial(a, n);
      alpha_factorial_computed = true;
    }
    // To compute numerator and denominator is useful because allows
    // to solve cases as `a/b * n + c/d': infact consider separately
    // `a*n + c*d' (that we are able to solve if `a = 1 && c/d is
    // positive integer' or `a = 2 && c*d = 1) and `b*d'.
    Expr num;
    Expr den;
    alpha.numerator_denominator(num, den);
    if (den != 1) {
      alpha_factorial = compute_alpha_factorial(num, n)
	* pwr(compute_alpha_factorial(den, n), -1);
      alpha_factorial_computed = true;
    }
  }
  if (!alpha_factorial_computed) {
    Symbol h("h");
    Expr alpha_h = alpha.subs(n, h);
    alpha_factorial = prod(Expr(h), Expr(1), Expr(n), alpha_h);
  }
  return alpha_factorial;
}

//! \brief
//! When possible, computes \f$ alpha!(n) \f$ if \f$ alpha \f$ is a
//! power, otherwise returns the symbolic product.
static Expr
alpha_factorial_if_power(const Expr& alpha, const Symbol& n) {
  Expr alpha_factorial;
  bool alpha_factorial_computed = false;
  if (alpha.op(0).has(n)) {
    if (alpha.op(1).is_a_number()) {
      alpha_factorial = pwr(compute_alpha_factorial(alpha.op(0), n),
			    alpha.op(1));
      alpha_factorial_computed = true;
    }
  }
  // In this case `alpha!(n) = k^{\sum_{h=1}^n f(h)}'.
  else {
    std::vector<Expr> base_of_exps;
    std::vector<Expr> exp_poly_coeff;
    std::vector<Expr> exp_no_poly_coeff;
    exp_poly_decomposition(alpha.op(1), n,
			   base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
    Expr new_exponent = 0;
    // `f(h)' is a polynomial or a product of a polynomial times an
    // exponential.
    if (vector_not_all_zero(exp_poly_coeff)) {
      Symbol k("k");
      for (unsigned i = base_of_exps.size(); i-- > 0; ) {
	Expr coeff_k = exp_poly_coeff[i].subs(n, k);
	new_exponent += sum_poly_times_exponentials(coeff_k, k, n,
						    base_of_exps[i]);
	// `sum_poly_times_exponentials' computes the sum from 0, while
	// we want it to start from `1'.
	new_exponent -= coeff_k.subs(k, 0);
      }
      alpha_factorial = pwr(alpha.op(0), new_exponent);
      alpha_factorial_computed = true;
    }
    // FIXME: aggiungere anche 
    // if (vector_not_all_zero(exp_no_poly_coeff)) {...}
    // per risolvere altre sommatorie risolvibili solo con gosper.
  }
  if (!alpha_factorial_computed) {
    Symbol h("h");
    Expr alpha_h = alpha.subs(n, h);
    alpha_factorial = prod(Expr(h), Expr(1), Expr(n), alpha_h);
  }
  return alpha_factorial;
}

//! \brief
//! Let \f$ alpha(n) \f$ be a not constant expression in the
//! variable \f$ n \f$.
//! This functions computes \f$ \alpha!(n) \f$ defined as follows:
//! \f[
//!   \alpha!(0) \defeq 1,
//!   \qquad
//!   \alpha!(n) \defeq \prod_{k=1}^n \alpha(k).
//! \f]
/*!
  When possible to find the closed form for \f$ \prod_{k=1}^n \alpha(k) \f$,
  we compute it; when it is not possible we returns the symbolic function
  for the product.
  We defined inductively \f$ \alpha!(n) \f$ as follows:
  - if \f$ alpha \f$ is a constant, i.e. it not contains \f$ n \f$,
    then \f$ alpha!(n) = alpha^n \f$;
  - if \f$ alpha = n \f$, then \f$ alpha!(n) = n! \f$;
  - if \f$ alpha = n + k \f$, where \f$ k \in \Nset \setminus \{0\} \f$,
    then \f$ alpha!(n) = \frac{(n + k)!}{k!}  \f$;
  - if \f$ alpha = 2*n+1 \f$,
    then \f$ alpha!(n) = \frac{(2*n + 1)!}{2^n * n} \f$;
  - if \f$ alpha \f$ is a power there are two cases.
    We consider \f$ a \f$ and \f$ b \f$ so that \f$ alpha = a^b \f$, 
    - if \f$ a \f$ contains \f$ n \f$ and \f$ b \f$ is a number,
      then \f$ alpha!(n) = a!(n)^b;
    - if \f$ a \f$ not contains \f$ n, i.e. \f$ a \f$ is a constant,
      then \f$ alpha!(n) = k^{\sum_{h=1}^n f(h)} \f$;
  - if \f$ alpha = alpha_1 \cdots alpha_k \f$, where \f$ alpha_i \f$,
    for \f$ i = 1, \dots, k \f$, is one of the previous case,
    then \f$ alpha!(n) =  alpha_1!(n) \cdots alpha_k!(n) \f$.  
*/
static Expr
compute_alpha_factorial(const Expr& alpha, const Symbol& n) {
  Expr alpha_factorial;
  if (!alpha.has(n))
    alpha_factorial = pwr(alpha, n);
  else if (alpha == n)
    alpha_factorial = factorial(alpha);
  else if (alpha.is_a_add())
    alpha_factorial = alpha_factorial_if_add(alpha, n);
  else if (alpha.is_a_power())
    alpha_factorial = alpha_factorial_if_power(alpha, n);
  else if (alpha.is_a_mul()) {
    alpha_factorial = 1;
    for (unsigned i = alpha.nops(); i-- > 0; )
      alpha_factorial *= compute_alpha_factorial(alpha.op(i), n);
  }
  else {
    Symbol h("h");
    Expr alpha_h = alpha.subs(n, h);
    alpha_factorial = prod(Expr(h), Expr(1), Expr(n), alpha_h);
  }
  return alpha_factorial;
}

/*!
  Consider the linear recurrence relation of first order with
  constant coefficients
  \f[
    x_n = \alpha(n) x_{n-1} + p(n),
  \f]
  where \f$ p(n) \f$ is a function defined over the natural numbers and
  \f$ \alpha \f$ is not constant.
  We set \f$ y_n \defeq x_n/\alpha!(n) \f$, where
  \f[
    \alpha!(0) \defeq 1,
    \qquad
    \alpha!(n) \defeq \prod_{k=1}^n \alpha(k).
  \f]
  We find that \f$ y_n \f$ satisfies the recurrence
  \f[
    y_n = y_{n-1} + \frac{p(n)}{\alpha!(n)}
  \f]
  whose solution is
  \f[
    y_n = x_0 + \sum_{k=1}^n \frac{p(k)}{k!},
  \f]
  so that our problem has been brought to the computation of a finite sum.
  At the end we find the formula for \f$ x_n \f$:
  \f[
    x_n = y_n \cdot \alpha!(n).
  \f]
*/
static Expr
solve_variable_coeff_order_1(const Symbol& n, const Expr& p_n,
			     const Expr& coefficient) {
  // Compute `\alpha!(n)' when possible.
  Expr alpha_factorial = compute_alpha_factorial(coefficient, n);
#if NOISY
  D_VAR(alpha_factorial);
#endif
  // Compute the non-homogeneous term for the recurrence
  // `y_n = y_{n-1} + \frac{p(n)}{\alpha!(n)}'.
  // In this case is better to jump a part of Gosper's step one:
  // `r(n) = \frac{t(n+1)}{t(n)}
  //       = \frac{p(n+1)}{\alpha!(n+1)} * \frac{\alpha!(n)}{p(n)}
  //       = \frac{p(n+1)}{p(n) * \alpha(n+1)}'.
  Expr new_p_n;
  if (!p_n.is_zero()) {
    Expr p_n_plus_one = p_n.subs(n, n+1);
    D_VAR(p_n_plus_one);
    Expr alpha_plus_one = coefficient.subs(n, n+1);
    D_VAR(alpha_plus_one);
    new_p_n = p_n_plus_one / (p_n * alpha_plus_one);
    D_VAR(new_p_n);
  }
  else
    new_p_n = 0;
  Expr new_rhs = x(n-1) + new_p_n;
//    D_VAR(new_rhs);
//    std::vector<Expr> base_of_exps;
//    std::vector<Expr> exp_poly_coeff;
//    std::vector<Expr> exp_no_poly_coeff;
//    exp_poly_decomposition(new_p_n, n,
//  			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
//    std::vector<Polynomial_Root> new_roots;
//    new_roots.push_back(Expr(1));

//    Expr solution;
//  #if 1
//    int order = 1;
//    Expr t_n = p_n/alpha_factorial;
//    if (!gosper_algorithm(t_n, order, n, base_of_exps,
//  			exp_poly_coeff, exp_no_poly_coeff,
//  			new_roots, solution)) {
//      // FIXME: the summand is not hypergeometric:
//      // no chance of using Gosper's algorithm.
//    }
//    solution += x(0);
//  #else
//    Expr solution = solve_constant_coeff_order_1(n, 1, base_of_exps,
//  					       exp_poly_coeff,
//  					       exp_no_poly_coeff, new_roots,
//  					       initial_conditions);
//  #endif
  //solution *= alpha_factorial;
  Expr solution;
#if NOISY
  D_VAR(new_p_n);
  D_VAR(new_rhs);
  D_MSGVAR("Solution new_rhs = ", solution);
#endif
  return solution;
}


static void
print_bad_exp(const Expr& e, const Expr rhs, bool conditions) {
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
verify_solution(const Expr& solution, int order, const Expr& rhs,
		const Symbol& n) {
  // FIXME: the initial conditions can not start always from 0:
  // `order' is temporary until we will consider a method in order to
  // know the right initial conditions.

  // Validation of initial conditions.
  for (int i = order; i-- > 0; ) {
    Expr g_i = x(i);
    Expr sol_subs = simplify_numer_denom(solution.subs(n, i));
    if (g_i != sol_subs) {
      print_bad_exp(sol_subs, rhs, true);
      return false;
    }
  }
  // The initial conditions are verified. Build an other expression
  // that has all terms of `solution' minus those containing an initial
  // condition.
  Expr partial_solution = 0;
  for (unsigned i = solution.nops(); i-- > 0; )
    if (!solution.op(i).match(x(wild(0)))
	&& !solution.op(i).match(wild(1) * x(wild(0))))
      partial_solution += solution.op(i);

  std::vector<Expr> terms_to_sub(order);
  for (int i = 0; i < order; ++i)
    terms_to_sub[i] = partial_solution.subs(n, n - i - 1);
  Expr substituted_rhs = simplify_on_input_ex(rhs.expand(), n, true);
  for (unsigned i = terms_to_sub.size(); i-- > 0; )
    substituted_rhs = substituted_rhs.subs(x(n - i - 1), terms_to_sub[i]);
  Expr diff = (partial_solution - substituted_rhs);
  // `simplify_factorials_and_exponentials()' must be call on not
  // expanded expression.
  diff = simplify_factorials_and_exponentials(diff, n).expand();
  diff = simplify_numer_denom(diff);
  if (!diff.is_zero()) {
    diff = simplify_factorials_and_exponentials(diff, n).expand();
    if (!diff.is_zero()) {
      print_bad_exp(diff, rhs, false);
      return false;
    }
  }
  return true;
}

} // namespace Parma_Recurrence_Relation_Solver
