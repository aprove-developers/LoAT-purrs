/* Recurrence class implementation: methods (and associated auxiliary
   functions) associated to the computation of the exact solution.
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

#include "util.hh"
#include "numerator_denominator.hh"
#include "simplify.hh"
#include "compute_prod.hh"
#include "ep_decomp.hh"
#include "finite_order.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Expr_List.defs.hh"
#include "Matrix.defs.hh"
#include "Finite_Order_Info.defs.hh"
#include "Cached_Expr.defs.hh"
#include "Recurrence.defs.hh"
#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

/*!
  This function adds the part of the solution of linear recurrences
  with constant coefficients of finite order corresponding to the initial
  conditions:
  \f[
    \sum_{i=0}^{order - 1} g_{n-i}
      \bigl( x_i - \sum_{j=1}^i a_j x_{i-j} \bigr).
  \f]
*/
// FIXME: il vettore `coefficients' dovra' diventare di `Expr' quando
// sapremo risolvere anche le eq. di grado superiore al primo con i
// parametri.
void
PURRS::Recurrence::
add_term_with_initial_conditions(const Expr& g_n,
				 const std::vector<Number>& coefficients) const {
  // `coefficients.size()' has `order + 1' elements because in the first
  // position there is the value 0.
  Expr solution = exact_solution_.expression();
  for (unsigned i = coefficients.size() - 1; i-- > 0; ) {
    Expr g_n_i = g_n.substitute(n, n - i);
    Expr tmp = x(first_well_defined_rhs_linear() + i);
    for (unsigned j = i; j > 0; j--)
      tmp -= coefficients[j] * x(first_well_defined_rhs_linear() + i - j);
    solution += tmp * g_n_i;
  }
  exact_solution_.set_expression(solution);
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
  In this function we compute, when possible, the sum of the previous
  formula; when this is not possible, we return the symbolic sum.
*/
PURRS::Expr
PURRS::Recurrence::
solve_constant_coeff_order_1(const std::vector<Polynomial_Root>& roots) const {
  // We search exponentials in `n' (for this the expression
  // `inhomogeneous_term' must be expanded).
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `inhomogeneous_term'. In the `i'-th position of
  // the vectors `exp_poly_coeff' and `exp_no_poly_coeff' there are
  // respectively the polynomial part and possibly non polynomial part of
  // the coefficient of the exponential with the base in `i'-th position of
  // `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of base_of_exps[i]^n.
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous_term, Recurrence::n,
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  Expr solution = 0;
  // Computes the sum when `\lambda^{n-k} p(k)' is a polynomial or
  // a product of a polynomial times an exponential.
  if (vector_not_all_zero(exp_poly_coeff)) {
    Symbol alpha("alpha");
    Symbol lambda("lambda");
    std::vector<Expr> symbolic_sum_distinct;
    std::vector<Expr> symbolic_sum_no_distinct;
    compute_symbolic_sum(alpha, lambda, roots,
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
  }
  // Computes the sum when `\lambda^{n-k} p(k)' is not a polynomial or
  // a product of a polynomial times an exponential.
  // The summand must be an hypergeometric term.
  if (vector_not_all_zero(exp_no_poly_coeff)) {
    Expr gosper_solution;
    if (compute_sum_with_gosper_algorithm(first_well_defined_rhs_linear() + 1,
					  n, base_of_exps, exp_no_poly_coeff,
					  roots, gosper_solution))
      solution += gosper_solution;
    else {
      // FIXME: the summand is not hypergeometric:
      // no chance of using Gosper's algorithm.
      Symbol h;
      unsigned lower = first_well_defined_rhs_linear() > 0
	? first_well_defined_rhs_linear() + 1 : 1;
      solution += PURRS::sum(h, lower, n,
			     pwr(roots[0].value(), n - h)
			     * inhomogeneous_term.substitute(n, h));
    }
  }
  return solution;
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
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::
solve_variable_coeff_order_1(const Expr& coefficient) const {
  if (has_parameters(coefficient)) {
    D_MSG("Variable coefficient with parameters");
    return TOO_COMPLEX;
  }
  // `z' will contain the largest positive or null integer, if it exists,
  // that cancel the denominator of the coefficient.
  // If this integer does not exist then `z' is left to 0.
  Number z = 0;
  if (!largest_positive_int_zero(coefficient, n, z))
    return TOO_COMPLEX;
  // Find the largest positive or null integer that cancel the denominator of
  // `inhomogeneous_term' and store it in `z' if it is bigger than the
  // current `z'.
  if (!inhomogeneous_term.is_zero())
    if (!largest_positive_int_zero(denominator(inhomogeneous_term), n, z))
      return TOO_COMPLEX;
  // The initial conditions will start from `z'.
  set_first_well_defined_rhs_linear(z.to_unsigned());
  Symbol index;
  Expr alpha_factorial
    = compute_product(index, z + 1,
		      transform_in_single_fraction(coefficient
						   .substitute(n, index)));
  // FIXME: this simplification simplifies the value of `alpha_factorial'
  // but not the solution because we need to the simplification about
  // factorials and exponentials for the output.
  //alpha_factorial = simplify_factorials_and_exponentials(alpha_factorial);
  D_VAR(alpha_factorial);
  // Compute the non-homogeneous term for the recurrence
  // `y_n = y_{n-1} + \frac{p(n)}{\alpha!(n)}'.
  // In this case is better to jump a part of Gosper's step one:
  // `r(n) = \frac{t(n+1)}{t(n)}
  //       = \frac{p(n+1)}{\alpha!(n+1)} * \frac{\alpha!(n)}{p(n)}
  //       = \frac{p(n+1)}{p(n) * \alpha(n+1)}'.
  Expr new_inhomogeneous_term;
  Expr solution = 0;
  if (!inhomogeneous_term.is_zero()) {
    new_inhomogeneous_term = inhomogeneous_term.substitute(n, n+1) 
      / (inhomogeneous_term * coefficient.substitute(n, n+1));
    new_inhomogeneous_term = simplify_all(new_inhomogeneous_term);
    D_VAR(new_inhomogeneous_term);
    std::vector<Expr> base_of_exps;
    std::vector<Expr> exp_poly_coeff;
    std::vector<Expr> exp_no_poly_coeff;
    exp_poly_decomposition(new_inhomogeneous_term, Recurrence::n,
			   base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
    std::vector<Polynomial_Root> new_roots;
    new_roots.push_back(Polynomial_Root(Expr(1), RATIONAL));
    if (!compute_sum_with_gosper_algorithm(first_well_defined_rhs_linear() + 1,
					   n, base_of_exps, exp_poly_coeff,
					   exp_no_poly_coeff, new_roots,
					   inhomogeneous_term/alpha_factorial,
					   solution)) {
      // FIXME: the summand is not hypergeometric:
      // no chance of using Gosper's algorithm.
      // vedere direttamente il rapporto p(k)/alpha!(k) se e' sommabile
      // (forse prima di vedere gosper)
      Symbol h;
      solution
	+= alpha_factorial * x(z)
	+ alpha_factorial * PURRS::sum(h, z + 1, n,
				       inhomogeneous_term.substitute(n, h) 
				       / alpha_factorial.substitute(n, h));
      exact_solution_.set_expression(solution);
      return SUCCESS;
    }
  }
  solution += x(z);
  solution *= alpha_factorial;
  exact_solution_.set_expression(solution);
  return SUCCESS;
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
  the function <CODE>add_term_with_initial_conditions()</CODE>), respectively.
*/
PURRS::Expr
PURRS::Recurrence::
solve_constant_coeff_order_2(Expr& g_n, bool all_distinct,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots) const {
  // We search exponentials in `n' (for this the expression
  // `inhomogeneous_term' must be expanded).
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `inhomogeneous_term'.
  // In the `i'-th position of the vectors
  // `exp_poly_coeff' and `exp_no_poly_coeff' there are respectively
  // the polynomial part and possibly non polynomial part of the coefficient
  // of the exponential with the base in `i'-th position of `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of base_of_exps[i]^n.
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous_term, Recurrence::n,
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  Expr solution;
  // Calculates the solution of the second order recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential; return the symbolic solution with
  // the object `sum' otherwise.
  if (all_distinct) {
    const Expr& root_1 = roots[0].value();
    const Expr& root_2 = roots[1].value();
    Expr diff_roots = root_1 - root_2;
    g_n = (pwr(root_1, Recurrence::n+1) - pwr(root_2, Recurrence::n+1))
      / diff_roots;
    if (!vector_not_all_zero(exp_no_poly_coeff)) {
      Symbol alpha("alpha");
      Symbol lambda("lambda");
      std::vector<Expr> symbolic_sum_distinct;
      std::vector<Expr> symbolic_sum_no_distinct;
      compute_symbolic_sum(alpha, lambda, roots,
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
    }
    else {
      Symbol h;
      unsigned lower = first_well_defined_rhs_linear() > 1
	? first_well_defined_rhs_linear() + 1 : 2;
      solution = 1 / diff_roots
	* (pwr(root_1, Recurrence::n+1)
	   * PURRS::sum(h, lower, Recurrence::n, pwr(root_1, -h)
			* inhomogeneous_term.substitute(Recurrence::n, h))
	   - (pwr(root_2, Recurrence::n+1)
	      * PURRS::sum(h, lower, Recurrence::n, pwr(root_2, -h)
			   * inhomogeneous_term.substitute(Recurrence::n, h))));
    }
  }
  else {
    // The characteristic equation
    // x^2 + a_1 * x + a_2 = 0 has a double root.
    assert(roots[0].multiplicity() == 2);      
    
    // Solve system in order to finds `alpha_i' (i = 1,...,order).
    Matrix sol = solve_system(all_distinct, coefficients, roots);
    
    // Finds `g_n', always taking into account the root's multiplicity.
    g_n = find_g_n(all_distinct, sol, roots);
    if (!vector_not_all_zero(exp_no_poly_coeff))
      solution = compute_non_homogeneous_part(g_n, order(), base_of_exps,
					      exp_poly_coeff);
    else {
      Symbol h;
      unsigned lower = first_well_defined_rhs_linear() > 1
	? first_well_defined_rhs_linear() + 1 : 2;
      solution
	= PURRS::sum(h, lower, Recurrence::n,
		     g_n.substitute(Recurrence::n, Recurrence::n - h)
		     * inhomogeneous_term.substitute(Recurrence::n, h));
    }
  }
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
        \quad \text{for } 1 \le n < k, \\
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
  the function <CODE>add_term_with_initial_conditions()</CODE>), respectively.
*/
PURRS::Expr
PURRS::Recurrence::
solve_constant_coeff_order_k(Expr& g_n, bool all_distinct,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots) const {
  // We search exponentials in `n' (for this the expression
  // `inhomogeneous_term' must be expanded).
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `inhomogeneous_term'.
  // In the `i'-th position of the vectors
  // `exp_poly_coeff' and `exp_no_poly_coeff' there are respectively
  // the polynomial part and possibly non polynomial part of the coefficient
  // of the exponential with the base in `i'-th position of `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of base_of_exps[i]^n.
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous_term, Recurrence::n,
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  Expr solution;
  // Calculates the solution of the recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential; return the symbolic solution with
  // the object `sum' otherwise.

  // Solve system in order to finds `alpha_i' (i = 1,...,order).
  Matrix sol = solve_system(all_distinct, coefficients, roots);

  // Finds `g_n', always taking into account the root's multiplicity
  g_n = find_g_n(all_distinct, sol, roots);
  D_VAR(g_n);

  if (all_distinct)
    if (!vector_not_all_zero(exp_no_poly_coeff)) {
      // Prepare for to compute the symbolic sum.
      std::vector<Expr> poly_coeff_tot;
      prepare_for_symbolic_sum(g_n, roots, exp_poly_coeff, poly_coeff_tot);
      Symbol alpha("alpha");
      Symbol lambda("lambda");
      std::vector<Expr> symbolic_sum_distinct;
      std::vector<Expr> symbolic_sum_no_distinct;
      compute_symbolic_sum(alpha, lambda, roots,
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
    else {
      Symbol h;
      unsigned lower = first_well_defined_rhs_linear() > order() - 1
	? first_well_defined_rhs_linear() + 1 : order();
      solution
	= PURRS::sum(h, lower, Recurrence::n,
		     g_n.substitute(Recurrence::n, Recurrence::n - h)
		     * inhomogeneous_term.substitute(Recurrence::n, h));
    }
  else
    if (!vector_not_all_zero(exp_no_poly_coeff))
      // There are roots with multiplicity greater than 1.
      solution = compute_non_homogeneous_part(g_n, order(), base_of_exps,
					      exp_poly_coeff);
    else {
      Symbol h;
      unsigned lower = first_well_defined_rhs_linear() > order() - 1
	? first_well_defined_rhs_linear() + 1 : order();
      solution
	= PURRS::sum(h, lower, Recurrence::n,
		     g_n.substitute(Recurrence::n, Recurrence::n - h)
		     * inhomogeneous_term.substitute(Recurrence::n, h));
    }
  return solution;
}

/*!
  Solves the linear recurrence of finite order (we have already considered
  the banal case of order zero recurrence), i. e., linear finite
  order recurrence with constant or variable coefficients.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::solve_linear_finite_order() const {
  D_VAR(order());
  D_VEC(coefficients(), 1, order());
  D_VAR(inhomogeneous_term);

  if (is_order_zero()) {
    exact_solution_.set_expression(inhomogeneous_term);
    exact_solution_.set_expression
      (simplify_ex_for_output(exact_solution_.expression(), false));
    exact_solution_.set_expression
      (simplify_factorials_and_exponentials(exact_solution_.expression()));
    exact_solution_.set_expression
      (simplify_logarithm(exact_solution_.expression()));
    return SUCCESS;
  }
  // Find the largest positive or null integer that cancel the denominator of
  // `inhomogeneous_term' and store it in `z' if it is bigger than `0'.
  Number z = 0;
  if (!inhomogeneous_term.is_zero()) {
    if (has_parameters(denominator(inhomogeneous_term))) {
      D_MSG("Linear finite order with parameters in the denominator of "
	    "the inhomogeneous term.");
      return TOO_COMPLEX;
    }
    // The system not finds an integer that cancel `inhomogeneous_term' or
    // starting from which `inhomogeneous_term' is well-defined
    // (polynomials are always well-defined).
    if (!largest_positive_int_zero(denominator(inhomogeneous_term), n, z))
      return TOO_COMPLEX;
  }
  // The initial conditions will start from `z'.
  set_first_well_defined_rhs_linear(z.to_unsigned());
  D_VAR(first_well_defined_rhs_linear());

  std::vector<Number> num_coefficients(order() + 1);
  std::vector<Polynomial_Root> roots;
  bool all_distinct = true;
  if (is_linear_finite_order_const_coeff()) {
    Expr characteristic_eq;
    if (!characteristic_equation_and_its_roots(order(), coefficients(),
					       num_coefficients,
					       characteristic_eq, roots,
					       all_distinct))
      return TOO_COMPLEX;
  }

  // `g_n' is defined here because it is necessary in the function
  // `add_term_with_initial_conditions()' (at the end of function
  // `solve_linear_finite_order()').
  Expr g_n;
  switch (order()) {
  case 1:
    {
      if (is_linear_finite_order_const_coeff())
	exact_solution_
	  .set_expression(solve_constant_coeff_order_1(roots));
      else {
	Solver_Status status;
	if ((status = solve_variable_coeff_order_1(coefficients()[1]))
	    != SUCCESS)
	  return status;
      }
    }
    break;

  case 2:
    if (is_linear_finite_order_const_coeff()) {
      // If there is some root not rational then, for efficiency, we substitute
      // it with an arbitrary symbol.
      substitute_non_rational_roots(*this, roots);
      exact_solution_
	.set_expression(solve_constant_coeff_order_2(g_n, all_distinct,
						     num_coefficients, roots));
    }
    else
      // For the time being, we only solve second order
      // recurrence relations with constant coefficients.
      return TOO_COMPLEX;
    break;

  default:
    if (is_linear_finite_order_const_coeff()) {
      // If there is some root not rational then, for efficiency, we substitute
      // it with an arbitrary symbol.
      substitute_non_rational_roots(*this, roots);
      exact_solution_.
	set_expression(solve_constant_coeff_order_k(g_n, all_distinct,
						    num_coefficients, roots));
    }
    else
      // For the time being, we only solve recurrence relations
      // of order 3 and more only if they have constant coefficients.
      return TOO_COMPLEX;
    break;
  }

  if (is_linear_finite_order_const_coeff())
    if (order() == 1 ) {
      // FIXME: per ora non si puo' usare la funzione
      // `add_term_with_initial_conditions' perche' richiede un vettore di
      // `Number' come `coefficients' e voglio risolvere anche le
      // parametriche (g_n pu' essere posta uguale ad 1 in questo caso).
      // add_term_with_initial_conditions(g_n, coefficients, solution);
      Expr solution = exact_solution_.expression();
      // The substitution of the initial condition will be later:
      // -  in the case of order recuction on the expanded solution 
      //    (and not on the reduced solution);
      // -  in the case of non linear recurrence on solution after
      //    the change of the variable.
      if (applied_order_reduction || come_from_non_linear_rec)
	solution
	  += x(first_well_defined_rhs_linear()) * pwr(coefficients()[1], n);
      else
	solution
	  += x(first_well_defined_rhs_linear()) * pwr(coefficients()[1], n);
      exact_solution_.set_expression(solution);
    }
    else
      add_term_with_initial_conditions(g_n, num_coefficients);

  // Resubstitutes eventually auxiliary definitions contained in
  // the solution with their original values.
  //exact_solution_.set_expression(blackboard.rewrite(solution));

  D_MSGVAR("Before calling simplify: ", exact_solution_.expression());
  exact_solution_.set_expression
    (simplify_ex_for_output(exact_solution_.expression(), false));
  exact_solution_.set_expression
    (simplify_factorials_and_exponentials(exact_solution_.expression()));
  exact_solution_.set_expression
    (simplify_logarithm(exact_solution_.expression()));

  // Only for the output.
  if (exact_solution_.expression().is_a_add()) {
    Expr_List conditions;
    for (unsigned i = order(); i-- > 0; )
      conditions.append(x(first_well_defined_rhs_linear() + i));
    // FIXME: `collect' throws an exception if the object to collect has
    // non-integer exponent. 
    exact_solution_
      .set_expression(exact_solution_.expression().collect(conditions));
  }

  return SUCCESS;
}
