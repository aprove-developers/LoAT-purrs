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
http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#ifndef NOISY
#define NOISY 0
#endif

#include "finite_order.hh"

#include "util.hh"
#include "gosper.hh"
#include "simplify.hh"
#include "sum_poly.hh"
#include "ep_decomp.hh"
#include "alg_eq_solver.hh"
#include "Expr.defs.hh"
#include "Expr_List.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Matrix.defs.hh"
#include "Recurrence.defs.hh"

#include <algorithm>
#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

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
Expr
build_characteristic_equation(const Symbol& x,
			      const std::vector<Number>& coefficients) {
  for (unsigned int i = coefficients.size(); i-- > 0; )
    assert(coefficients[i].is_rational());

  std::vector<Number> denominators;
  // Find the least common multiple of the denominators of the
  // rational elements of `coefficients'.
  for (unsigned int i = coefficients.size(); i-- > 0; )
    if (!coefficients[i].is_integer())
      denominators.push_back(coefficients[i].denominator());
  Expr p = 0;
  // Build the vector `int_coefficients' with the elements of
  // `coefficients' multiplied by the least common multiple
  // `least_com_mul'.
  // There is no need to know the `order' because the order of the
  // recurrence relation is equal to `coefficients.size() - 1'.
  if (denominators.size() != 0) {
    Number least_com_mul = lcm(denominators);
    std::vector<Number> int_coefficients(coefficients);
    // In the first position of `coefficients' there is 0.
    for (unsigned int i = coefficients.size(); i-- > 1; )
      int_coefficients[i] *= least_com_mul;
    for (unsigned int i = coefficients.size() - 1; i-- > 0; )
      p += pwr(x, i) * (-int_coefficients[coefficients.size() - 1 - i]);
    p += least_com_mul * pwr(x, coefficients.size() - 1);
  }
  else {
    for (unsigned int i = coefficients.size() - 1; i-- > 0; )
      p += pwr(x, i) * (-coefficients[coefficients.size() - 1 - i]);
    p += pwr(x, coefficients.size() - 1);
  }
  return p;
}

Expr
return_sum(bool distinct, unsigned int lower_bound_sum, const Expr& coeff,
	   const Symbol& alpha, const Symbol& lambda) {
  Symbol k("k");
  Symbol x("x");
  Expr q_k = coeff.substitute(Recurrence::n, k);
  Expr symbolic_sum;  
  if (distinct)
    symbolic_sum = sum_poly_times_exponentials(q_k, k, Recurrence::n, x);
  else
    symbolic_sum = sum_poly_times_exponentials(q_k, k, Recurrence::n, 1);
  // `sum_poly_times_exponentials' computes the sum from 0, whereas
  // we want that the sum start from `lower_bound_sum'.
  for (Number j = 0; j < lower_bound_sum; ++j)
    symbolic_sum -= q_k.substitute(k, j) * pwr(alpha, j) * pwr(lambda, -j);
  if (distinct)
    symbolic_sum = symbolic_sum.substitute(x, alpha/lambda);
  symbolic_sum *= pwr(lambda, Recurrence::n);
  symbolic_sum = simplify_ex_for_output(symbolic_sum, false);
  return symbolic_sum;
}

Expr
rewrite_factor(const Expr& e, const Symbol& r, unsigned int gcd_among_decrements) {
  if (e.is_a_power())
    return pwr(rewrite_factor(e.arg(0), r, gcd_among_decrements),
	       rewrite_factor(e.arg(1), r, gcd_among_decrements));
  else if (e.is_a_function())
    if (e.is_the_x_function()) {
      Expr argument = e.arg(0);
      assert(argument.is_a_add() && argument.nops() == 2);
      assert(argument.op(0).is_a_number() || argument.op(1).is_a_number());
      Number decrement;
      argument.op(0).is_a_number(decrement)
	|| argument.op(1).is_a_number(decrement);
      return x(Recurrence::n + decrement / gcd_among_decrements);
    }
    else if (e.nops() == 1)
      return apply(e.functor(),
		   rewrite_factor(e.arg(0), r, gcd_among_decrements));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
	argument[i] = rewrite_factor(e.arg(i), r, gcd_among_decrements);
      return apply(e.functor(), argument);
    }
  else if (e == Recurrence::n)
    return gcd_among_decrements * Recurrence::n + r;
  return e;
}

/*!
  ...
*/
Expr
rewrite_term(const Expr& e, const Symbol& r, unsigned int gcd_among_decrements) {
  unsigned int num_factors = e.is_a_mul() ? e.nops() : 1;
  Expr e_rewritten = 1;
  if (num_factors > 1)
    for (unsigned int i = num_factors; i-- > 0; )
      e_rewritten *= rewrite_factor(e.op(i), r, gcd_among_decrements);
  else
    e_rewritten = rewrite_factor(e, r, gcd_among_decrements);
  return e_rewritten;
}

} // anonymous namespace

/*
  Builds the characteristic equation and computes its roots.
  The boolean \p all_distinct is true if all roots are distinct, i.e., every
  roots has multiplicity equal to one, false otherwise.
  Returns <CODE>false</CODE> if it does not succeed to find roots of the
  characteristic equations; returns <CODE>true<CODE> otherwise. 
*/
bool
PURRS::
characteristic_equation_and_its_roots(index_type order,
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
    // `NON_RATIONAL' because the coefficient could be a parameter. 
    roots.push_back(Polynomial_Root(coefficients[1], NON_RATIONAL));
  }
  else {
    // Check if the vector `coefficients' contains only numeric
    // elements and in this case use a vector of Number.
    for (unsigned int i = coefficients.size(); i--> 0; )
      if (coefficients[i].is_a_number())
	num_coefficients[i] = coefficients[i].ex_to_number();
      else
	return false;
    characteristic_eq = build_characteristic_equation(y, num_coefficients);
    D_VAR(characteristic_eq);
    if (!find_roots(characteristic_eq, y, roots, all_distinct))
      return false;
  }
  for (unsigned int i = roots.size(); i-- > 0; ) {
    D_VAR(roots[i].value());
    D_VAR(roots[i].multiplicity());
  }
  return true;
}

/*!
  Consider an inhomogeneous term \f$ e(n) \f$ which is a sum of products of 
  polynomials and exponentials:
  \f$ e(n) = \sum_{i=0}^k \alpha_i^{n} p_i(n) \f$ (where \f$ k \f$ is
  the number of exponentials, and we assume that the \f$ \alpha_j \f$'s
  are distinct complex numbers).
  \f$ \lambda \f$ denote the generic root of the characteristic 
  equation and \f$ \alpha \f$ the generic base of an exponential.

  This function fills the two vectors of <CODE>Expr</CODE>,
  \p symbolic_sum_distinct and \p symbolic_sum_no_distinct,
  with two different sums: first fix the base \f$ \alpha_i \f$.
  For each root \f$ \lambda \f$
  - if \f$ \alpha_i \neq \lambda \f$ then
    \f[
      symbolic{\_}sum{\_}distinct[i]
        = \lambda^n * f_i(\alpha_i / \lambda)
        = \lambda^n * \sum_{k=h}^n (\alpha_i / \lambda)^k \cdot p_i(k);
    \f]
  - if \f$ \alpha_i = \lambda \f$ then
    \f[
      symbolic{\_}sum{\_}no{\_}distinct[i]
        = \lambda^n * f_i(1)
        = \lambda^n * \sum_{k=h}^n p_i(k).
    \f]
  The lower bound `h' of the sums is the order of the recurrence
  shifted forward of a positive integer that indicates the least
  non-negative integer \f$ j \f$ such that the recurrence is
  well-defined for \f$ n \geq j \f$.
  We put \f$ 0 \f$ in the i-th position of \p symbolic_sum_no_distinct,
  in the first case, and of \p symbolic_sum_distinct, in the second case,
  so that the two vectors have always the same dimensions.
*/
void
PURRS::compute_symbolic_sum(const Symbol& alpha, const Symbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    const std::vector<Expr>& base_of_exps,
			    const std::vector<Expr>& exp_poly_coeff,
			    unsigned int lower_bound_sum,
			    std::vector<Expr>& symbolic_sum_distinct,
			    std::vector<Expr>& symbolic_sum_no_distinct) {
  // Compute the order of the recurrence relation.
  Number order = 0;
  for (unsigned int i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();

  unsigned int r = 0;
  for (unsigned int i = base_of_exps.size(); i-- > 0; )
    for (unsigned int j = roots.size(); j-- > 0; ) {
      bool distinct = true;
      if (roots[j].value() == base_of_exps[i])
	distinct = false;
      
      // The root is different from the exponential's base.
      if (distinct) {
	if (order <= 2)
	  symbolic_sum_distinct.push_back(return_sum(true, lower_bound_sum,
						     exp_poly_coeff[i],
						     alpha, lambda));
	else
	  symbolic_sum_distinct.push_back(return_sum(true, lower_bound_sum,
						     exp_poly_coeff[r],
						     alpha, lambda));
	symbolic_sum_no_distinct.push_back(0);
      }
      // The root is equal to the exponential's base.
      else {
	if (order <= 2)
	  symbolic_sum_no_distinct.push_back(return_sum(false, lower_bound_sum,
							exp_poly_coeff[i],
							alpha, lambda));
	else
	  symbolic_sum_no_distinct.push_back(return_sum(false, lower_bound_sum,
							exp_poly_coeff[r],
							alpha, lambda));
	symbolic_sum_distinct.push_back(0);
      }
      ++r;
    }
  D_VEC(symbolic_sum_distinct, 0, symbolic_sum_distinct.size()-1);
  D_VEC(symbolic_sum_no_distinct, 0, symbolic_sum_no_distinct.size()-1);
}

/*!
  Consider the vectors \p symbolic_sum_distinct and \p symbolic_sum_no_distinct
  that contain all the symbolic sums of the inhomogeneous term's terms that
  are polynomial or the product of a polynomial and an exponential,
  For each sum this function
  - substitutes to \p lambda the corresponding value of the characteristic
    equation's root;
  - substitutes to \p alpha the corresponding base of the eventual
    exponential.
  Returns a new expression containing the sum of all sums of the vectors
  after the substitutions.
*/
PURRS::Expr
PURRS::subs_to_sum_roots_and_bases(const Symbol& alpha, const Symbol& lambda,
				   const std::vector<Polynomial_Root>& roots,
				   const std::vector<Expr>& base_of_exps,
				   std::vector<Expr>& symbolic_sum_distinct,
				   std::vector<Expr>& symbolic_sum_no_distinct) {
  // Compute the order of the recurrence relation.
  Number order = 0;
  for (unsigned int i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  Expr solution = 0;
  unsigned int r = 0;
  for (unsigned int i = base_of_exps.size(); i-- > 0; )
    for (unsigned int j = roots.size(); j-- > 0; ) {
      const Expr& base_exp = base_of_exps[i];
      Expr tmp;
      if (base_exp != roots[j].value())
	tmp = symbolic_sum_distinct[r];
      else
	tmp = symbolic_sum_no_distinct[r];
      tmp = tmp.substitute(alpha, base_exp);
      tmp = tmp.substitute(lambda, roots[j].value());
      if (order == 2 && (j & 1) == 1)
	solution -= tmp;
      else
	// order != 2 or j even.
	solution += tmp;
      ++r;
    }
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
        \quad \text{for } 1 \le n < k, \\
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
PURRS::Matrix
PURRS::solve_system(bool all_distinct,
		    const std::vector<Number>& coefficients,
		    const std::vector<Polynomial_Root>& roots) {
  unsigned int coefficients_size = coefficients.size();
  // Prepare a list with the elments for the right hand side of the system
  // to solve.
  // The elements of the list are `g_0, g_1, ..., g_{k-1}'.
  // Note that `tmp[i]' is built on `tmp[i-1]'.
  std::vector<Expr> tmp(coefficients_size - 1);
  tmp[0] = 1;
  for (unsigned int i = 1; i < coefficients_size - 1; ++i)
    for (unsigned int j = 0; j < i; ++j)
      tmp[i] += coefficients[j+1] * tmp[i-j-1];
  Expr_List g_i;
  for (unsigned int i = coefficients_size - 1; i-- > 0; )
    g_i.prepend(tmp[i]);
  
  // Prepare a list with the coefficients of the equations of the system
  // to solve. This calculus is based on different form of `g_n'
  // in according to the roots' multiplicity.
  Expr_List coeff_equations;
  if (all_distinct)
    for (unsigned int i = coefficients_size - 1; i-- > 0; )
      for (unsigned int j = coefficients_size - 1; j-- > 0; )
	coeff_equations.prepend(pwr(roots[j].value(), i));
  else
    for (unsigned int h = 0; h < coefficients_size - 1; ++h)
      for (unsigned int i = roots.size(); i-- > 0; ) {
	for (unsigned int j = roots[i].multiplicity(); j-- > 1; )
	  coeff_equations.append(exact_pwr(h, j) * pwr(roots[i].value(), h));
	coeff_equations.append(pwr(roots[i].value(), h));
      }

  // Define the matrices and solve the system.
  Matrix coeff_alpha(coefficients_size - 1, coefficients_size - 1,
		      coeff_equations);
  Matrix rhs(coefficients_size - 1, 1, g_i);
  Matrix vars(coefficients_size - 1, 1);
  for (unsigned int i = coefficients_size - 1; i-- > 0; )
    vars(i, 0) = Symbol();
  // FIXME: in the case of `all_distinct = true' we have a Vandermonde's
  // matrix.
  // It is more efficient to solve the system with the method of inverse
  // matrix, where the closed form in order to compute the inverse of
  // Vandermonde's matrix is given in "Knuth, Fundamental Algorithms,
  // Addison-Wesley Publishing Company" (second edition pag.36).
  Matrix solution = coeff_alpha.solve(vars, rhs);

  return solution;
}

/*!
  ...
*/
PURRS::Expr
PURRS::find_g_n(bool all_distinct, const Matrix& sol,
		const std::vector<Polynomial_Root>& roots) {
  // Compute the order of the recurrence relation.
  Number order = 0;
  for (unsigned int i = roots.size(); i-- > 0; )
    order += roots[i].multiplicity();
  Expr g_n = 0;
  if (all_distinct)
    for (index_type i = 0; i < order; ++i)
      g_n += sol(i, 0) * pwr(roots[i].value(), Recurrence::n);
  else
    for (unsigned int i = roots.size(), h = 0; i-- > 0; )
      for (Number j = roots[i].multiplicity(); j-- > 0; h++)
	g_n += sol(h, 0) * pwr(Recurrence::n, j)
	  * pwr(roots[i].value(), Recurrence::n);
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
void
PURRS::prepare_for_symbolic_sum(const Expr& g_n,
				const std::vector<Polynomial_Root>& roots,
				const std::vector<Expr>& exp_poly_coeff,
				std::vector<Expr>& poly_coeff_tot) {
  std::vector<Expr> bases_exp_g_n;
  std::vector<Expr> g_n_poly_coeff;
  std::vector<Expr> g_n_no_poly_coeff;
  exp_poly_decomposition(g_n.expand(), Recurrence::n, bases_exp_g_n,
			 g_n_poly_coeff, g_n_no_poly_coeff);
  // `bases_of_exp_g_n' must have same elements of `roots' in the same order.
  bool equal = true;
  for (unsigned int i = roots.size(); i-- > 0; )
    if (roots[i].value() != bases_exp_g_n[i])
      equal = false;
  if (!equal) {
    std::vector<Expr> tmp_exp(roots.size());
    std::vector<Expr> tmp_coeff_poly(roots.size());
    std::vector<Expr> tmp_coeff_no_poly(roots.size());
    for (unsigned int i = roots.size(); i-- > 0; )
      tmp_exp[i] = roots[i].value();
    for (unsigned int i = tmp_exp.size(); i-- > 0; )
      for (unsigned int j = bases_exp_g_n.size(); j-- > 0; )
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
  for (unsigned int i = exp_poly_coeff.size(); i-- > 0; )
    for (unsigned int j = g_n_poly_coeff.size(); j-- > 0; )
      poly_coeff_tot.push_back(exp_poly_coeff[i] * g_n_poly_coeff[j]);
}

/*!
  ...
*/
PURRS::Expr
PURRS::compute_non_homogeneous_part(const Expr& g_n, index_type order,
				    const std::vector<Expr>& base_of_exps,
				    const std::vector<Expr>& exp_poly_coeff) {
  Expr solution_tot = 0;
  std::vector<Expr> bases_exp_g_n;
  std::vector<Expr> g_n_poly_coeff;
  std::vector<Expr> g_n_no_poly_coeff;
  exp_poly_decomposition(g_n.expand(), Recurrence::n, bases_exp_g_n,
			 g_n_poly_coeff, g_n_no_poly_coeff);
  for (unsigned int i = bases_exp_g_n.size(); i-- > 0; )
    for (unsigned int j = base_of_exps.size(); j-- > 0; ) {
      Expr solution = 0;
      Symbol k("k");
      Expr g_n_coeff_k = g_n_poly_coeff[i].substitute(Recurrence::n,
						      Recurrence::n - k);
      Expr exp_poly_coeff_k = exp_poly_coeff[j].substitute(Recurrence::n, k);
      solution
	= sum_poly_times_exponentials(g_n_coeff_k * exp_poly_coeff_k, k,
				      Recurrence::n,
				      1/bases_exp_g_n[i] * base_of_exps[j]);
      // `sum_poly_times_exponentials' computes the sum from 0, whereas
      // we want that the sum start from `order'.
      for (index_type h = 0; h < order; ++h)
	solution -= (g_n_coeff_k * exp_poly_coeff_k).substitute(k, h)
	  * pwr(1/bases_exp_g_n[i] * base_of_exps[j], h);
      solution *= pwr(bases_exp_g_n[i], Recurrence::n);
      solution_tot += solution;
    }
  return solution_tot;
}

/*!
  Applies the Gosper's algorithm to express in closed form, if it is
  possible, sum with the summand an hypergeometric term not polynomials
  or polynomials times exponentials (these last sums are computed by the
  functions <CODE>compute_symbolic_sum()</CODE> and
  <CODE>subs_to_sum_roots_and_bases()</CODE>).
*/
PURRS::Expr
PURRS::
compute_sum_with_transcendental_method(const Number& lower, const Expr& upper,
				       const std::vector<Expr>& base_of_exps,
				       const std::vector<Expr>&
				       exp_no_poly_coeff,
				       const std::vector<Polynomial_Root>&
				       roots) {
  Expr solution = 0;
  for (unsigned int i = exp_no_poly_coeff.size(); i-- > 0; ) {
    Expr gosper_solution;
    if (!exp_no_poly_coeff[i].is_zero()) {
      // FIXME: for the moment use this function only when the `order'
      // is one, then `roots' has only one elements.
      Expr t_n = pwr(base_of_exps[i], Recurrence::n) * exp_no_poly_coeff[i]
	* pwr(roots[0].value(), -Recurrence::n);
      D_VAR(t_n);
      if (!gosper_algorithm(Recurrence::n, t_n, lower, upper,
			    gosper_solution)) {
	// The summand is not hypergeometric:
	// no chance of using Gosper's algorithm.
	Symbol h;
	gosper_solution += PURRS::sum(h, lower, Recurrence::n,
				      pwr(roots[0].value(), -h)
				      * exp_no_poly_coeff[i]
				      .substitute(Recurrence::n, h)
				      * pwr(base_of_exps[i], Recurrence::n)
				      .substitute(Recurrence::n, h));
      }
    }
    solution += pwr(roots[0].value(), Recurrence::n) * gosper_solution;
  }
  return solution;
}

//! \brief
//! Let \p old_rhs be the right hand side of a recurrence to which is
//! possible to reduce the order. This function computes and returns
//! the right hand side of the reduced order recurrence.
/*!
  Let \f$ x(n) = a_1 x(n-k_1) + \dotsb + a_h x(n-k_h) + p(n) \f$ be
  a recurrence such that \f$ g = gcd(k_1, \dotsc, k_h) > 1 \f$. Its
  right hand side is contained in \p old_rhs, while \p r is a symbol
  that represent the quantity \f$ mod(n, g) \f$.
  For these recurrences is possible to reduce the order.
  This function computes and returns the right hand side of the reduced
  order recurrence.
  Note: the expression \p old_rhs must be simplified so that not to have
  nested power (using the function <CODE>simplify_ex_for_input()</CODE>).
*/
PURRS::Expr
PURRS::write_reduced_order_recurrence(const Expr& old_rhs, const Symbol& r,
				      unsigned int gcd_among_decrements,
				      const std::vector<Expr>& coefficients,
				      std::vector<Expr>& new_coefficients,
				      Expr& inhomogeneous) {
  unsigned int num_summands = old_rhs.is_a_add() ? old_rhs.nops() : 1;
  // Build the right hand side of the reduced order recurrences and,
  // at the same time, build its inhomogeneous term.
  Expr new_rhs = 0;
  if (num_summands > 1)
    for (unsigned int i = num_summands; i-- > 0; ) {
      const Expr& tmp = rewrite_term(old_rhs.op(i), r, gcd_among_decrements);
      new_rhs += tmp;
      if (!tmp.has_x_function(Recurrence::n))
	inhomogeneous += tmp;
    }
  else {
    const Expr& tmp = rewrite_term(old_rhs, r, gcd_among_decrements);
    new_rhs += tmp;
    if (!tmp.has_x_function(Recurrence::n))
      inhomogeneous += tmp;
  }
  // Find the coefficients of the reduced order recurrence.
  for (unsigned int i = coefficients.size(); i-- > 0; )
    if (i % gcd_among_decrements == 0)
      new_coefficients[i / gcd_among_decrements] = coefficients[i];
  return new_rhs;
}

//! Let \p e be the expression that represent the solution of the reduced
//! order recurrence: starting from it, this function computes and returns
//! a new expression containing the solution of the original recurrence.
/*!
  \param e                      The solution of the reduced order recurrence.
  \param r                      The symbol used instead of the expression
                                \f$ mod(n, gcd_among_decrements) \f$.
  \param m                      The expression
                                \f$ mod(n, gcd_among_decrements) \f$.
  \param gcd_among_decrements   The greatest common divisor among the
                                decrements \f$ d \f$ of the terms
				\f$ x(n-d) \f$ in the right hand side
				of the original recurrence.

  \return                       A new expression containing the solution
                                of the original recurrence.

  This function crosses inductively the expression \p e and performs
  the following substitutions:
  - r                            -> mod(n, gcd_among_decrements);
  - n                            -> 1 / gcd_among_decrements
                                    * (n - mod(n, gcd_among_decrements));
  - x(k), k non-negative integer -> x(mod(n, gcd_among_decrements))
                                     + k * gcd_among_decrements.
*/
PURRS::Expr 
PURRS::come_back_to_original_variable(const Expr& e, const Symbol& r,
				      const Expr& m,
				      unsigned int gcd_among_decrements) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += come_back_to_original_variable(e.op(i), r,
						    m, gcd_among_decrements);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten *= come_back_to_original_variable(e.op(i), r,
						    m, gcd_among_decrements);
  }
  else if (e.is_a_power())
    return pwr(come_back_to_original_variable(e.arg(0), r, m,
					      gcd_among_decrements),
	       come_back_to_original_variable(e.arg(1), r, m,
					      gcd_among_decrements));
  else if (e.is_a_function())
    if (e.is_the_x_function()) {
      assert(e.arg(0).is_a_number());
      Number arg = e.arg(0).ex_to_number();
      return x(gcd_among_decrements * arg + m);
    }
    else if (e.nops() == 1)
      return apply(e.functor(),
		   come_back_to_original_variable(e.arg(0), r,
						  m, gcd_among_decrements));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
	argument[i] = come_back_to_original_variable(e.arg(i), r,
						     m, gcd_among_decrements);
      return apply(e.functor(), argument);
    }
  else if (e == Recurrence::n)
    return Number(1, gcd_among_decrements) * (Recurrence::n - m);
  else if (e == r)
    return m;
  else
    return e;
  return e_rewritten;
}

//! \brief
//! Returns the expanded solution of the recurrence \p rec
//! to which we have applied the order reduction.
/*!
  There are two different steps:
  1)  the problem of expressing the part depending on the initial conditions;
  2)  to give an expression for \f$ r = n \bmod k \f$.
  
  Let \f$ \omega_h \f$, for \f$h, 1 \dots, k \f$ the \f$ k \f$-th roots
  of unity.
  1) Let
     \f[
       \eta_{r,k}(n) \defeq \frac 1k \sum_{h=1}^k \omega_h^{n-r},
     \f]
     so that \f$ \eta_{r,k}(n) = 1 \f$ if \f$ n \equiv r \bmod k \f$,
     and is~\f$ 0 \f$ otherwise.
     Therefore we have
     \f[
       x_r = \sum_{s=0}^{k-1} \eta_{s,k}(n) x_s.
     \f]

  2) \f[
       n \bmod k
       =
       \frac12 (k-1) + \sum_{h=2}^k \frac{\omega_h}{1-\omega_h} \, \omega_h^n.
     \f]
*/
PURRS::Expr
PURRS::Recurrence::write_expanded_solution(const Recurrence& rec,
					   unsigned int gcd_among_decrements) {
  // `term_with_ic' will contain the term of the solution relative to
  // the initial condition; `remainder_solution' will contain the
  // rest of the solution.
  Expr term_with_ic = 0;
  Expr remainder_solution = 0;
  if (rec.exact_solution_.expression().is_a_add())
    for (unsigned int i = rec.exact_solution_.expression().nops(); i-- > 0; )
      if (rec.exact_solution_.expression()
	  .op(i).has_symbolic_initial_conditions())
	term_with_ic += rec.exact_solution_.expression().op(i);
      else
	remainder_solution += rec.exact_solution_.expression().op(i);
  else {
    assert(rec.exact_solution_.expression().is_the_x_function());
    term_with_ic = rec.exact_solution_.expression();
  }

  // Compute `to_sub_in_solution': it is the value to substitute to the
  // function `mod()'.
  Expr to_sub_in_solution = 0;
  const Expr& theta = 2*Constant::Pi/gcd_among_decrements;
  for (unsigned int j = 1; j <= gcd_among_decrements; ++j) {
    const Expr& root_of_unity = cos(j*theta) + Number::I*sin(j*theta);
    // Skip the contribution of the root of unity equal to `1'.
    if (root_of_unity != 1)
      to_sub_in_solution += root_of_unity * pwr(1 - root_of_unity, -1)
	* pwr(root_of_unity, Recurrence::n);
  }
  // Add the contribution of the root of unity equal to `1'.
  to_sub_in_solution += Number(1, 2) * (gcd_among_decrements - 1);

  // 1. rewrite the part of the solution depending on the
  // initial conditions.
  for (unsigned int h = rec.order()/gcd_among_decrements; h-- > 0; ) {
    Expr initial_condition = 0;
    for (unsigned int i = 0; i < gcd_among_decrements; ++i) {
      Expr tmp = 0;
      for (unsigned int j = 1; j <= gcd_among_decrements; ++j) {	  
	const Expr& root_of_unity = cos(j*theta) + Number::I*sin(j*theta);
	tmp += pwr(root_of_unity, Recurrence::n - (i + h));
      }
      D_VAR(i + h * gcd_among_decrements);
      initial_condition += tmp * x(i + h * gcd_among_decrements)
	/ gcd_among_decrements;
    }
    D_VAR(x(mod(n, gcd_among_decrements) + h * gcd_among_decrements));
    D_VAR(initial_condition);
    term_with_ic = term_with_ic.substitute(x(mod(n, gcd_among_decrements)
					     + h * gcd_among_decrements),
					   initial_condition);
  }
  term_with_ic = term_with_ic.substitute(mod(n, gcd_among_decrements),
					 to_sub_in_solution);

  // 2. rewrite the part of the solution not depending on the
  // initial conditions.
  remainder_solution
    = remainder_solution.substitute(mod(n, gcd_among_decrements),
				    to_sub_in_solution);
  D_VAR(remainder_solution);
  
  return simplify_ex_for_output(term_with_ic + remainder_solution, false);
}

/*!
  ...
*/
void
PURRS::substitute_non_rational_roots(const Recurrence& rec,
				     std::vector<Polynomial_Root>& roots) {
  for (unsigned int i = roots.size(); i-- > 0; )
    if (roots[i].is_non_rational())
      roots[i].value() = rec.insert_auxiliary_definition(roots[i].value());
  for (unsigned int i = roots.size(); i-- > 0; )
    D_VAR(roots[i].value());
}

