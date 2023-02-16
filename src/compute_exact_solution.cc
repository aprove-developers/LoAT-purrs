/* Recurrence class implementation: methods (and associated auxiliary
   functions) associated to the computation of the exact solution.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#ifndef NOISY
#define NOISY 0
#endif

#include "util.hh"
#include "simplify.hh"
#include "sum_poly.hh"
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
#include "Recurrence.inlines.hh"

#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

/*!
  This function returns the part of the solution of linear recurrences
  with constant coefficients of finite order corresponding to the initial
  conditions:
  \f[
    \sum_{i=0}^{order - 1} g_{n-i}
      \bigl( x_i - \sum_{j=1}^i a_j x_{i-j} \bigr).
  \f]
  Suitable modifications are considered when the recurrence is not
  well-defined for every natural number, i.e., when
  \ref first_valid_index "first_valid_index" is not \f$ 0 \f$.
*/
// FIXME: il vettore `coefficients' dovra' diventare di `Expr' quando
// sapremo risolvere anche le eq. di grado superiore al primo con i
// parametri.
Expr
compute_term_about_initial_conditions(const Expr& g_n,
				      const std::vector<Number>& coefficients,
				      index_type first_valid_index) {
  Expr term = 0;
  // `coefficients.size() - 1' is equal to the order of the recurrence.
  for (unsigned int i = coefficients.size() - 1; i-- > 0; ) {
    const Expr& g_n_i
      = g_n.substitute(Recurrence::n, Recurrence::n - i - first_valid_index);
    Expr tmp = x(first_valid_index + i);
    for (unsigned int j = i; j > 0; j--)
      tmp -= coefficients[j] * x(first_valid_index + i - j);
    term += tmp * g_n_i;
  }
  return term;
}

} // anonymous namespace


/*!
  Consider the linear recurrence relation of first order with
  constant coefficients
  \f[
    x_n = \lambda x_{n-1} + p(n),
  \f]
  where \f$ p(n) \f$ is a function defined over the natural numbers.
  In this case we know the final formula that gives the solution (we remark
  that the coefficient coincides with the root of the characteristic
  equation):
  \f[
    x_n = \lambda^n * x_0
          + \sum_{k=1}^n \lambda^{n-k} p(k).
  \f]
  This function computes, when possible, the closed form for the sum
  of the previous formula; when this is not possible, the symbolic sum
  is returned.
*/
PURRS::Expr
PURRS::Recurrence::
solve_constant_coeff_order_1(const std::vector<Polynomial_Root>& roots) const {
  assert(is_linear_finite_order_const_coeff());
  assert(order() == 1);
  // We search exponentials in `Recurrence::n'.
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `inhomogeneous_term'. In the `i'-th position of
  // the vectors `exp_poly_coeff' and `exp_no_poly_coeff' there are
  // respectively the polynomial part and possibly non polynomial part of
  // the coefficient of the exponential with the base in `i'-th position of
  // `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of `base_of_exps[i]^n'.
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous_term.expand(), Recurrence::n,
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  Expr solution = 0;
  // Computes the sum when `\lambda^{n-k} p(k)' is a polynomial or
  // a product of a polynomial times an exponential.
  if (vector_not_all_zero(exp_poly_coeff)) {
    Symbol alpha;
    Symbol lambda;
    std::vector<Expr> symbolic_sum_distinct;
    std::vector<Expr> symbolic_sum_no_distinct;
    compute_symbolic_sum(alpha, lambda, roots, base_of_exps, exp_poly_coeff,
			 order() + first_valid_index, 
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
  // a product of a polynomial times an exponential with transcendental
  // methods.
  if (vector_not_all_zero(exp_no_poly_coeff))
    // If the summand is an hypergeometric term then is applied the
    // Gosper's algorithm; otherwise is returned the symbolic sum.
    solution += compute_sum_with_transcendental_method(first_valid_index + 1,
						       n, base_of_exps,
						       exp_no_poly_coeff,
						       roots);
  return solution;
}

/*!
  Consider the linear recurrence relation of first order with
  variable coefficient
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
solve_variable_coeff_order_1(const std::vector<Expr>& coefficients,
			     Expr& solution) const {
  assert(is_linear_finite_order_var_coeff());
  assert(order() == 1);
  // `product_factor' is `alpha!(n)'.
  Expr numerator;
  Expr denominator;
  Symbol index;
  coefficients[1].substitute(n, index).numerator_denominator(numerator,
							     denominator);
  const Expr& product_factor
    = simplify_binomials_factorials_exponentials
    (compute_product(index, first_valid_index + 1, numerator/denominator));
  set_product_factor(product_factor);

  // Build the recurrence with constant coefficient of the first order
  // `y_n = y_{n-1} + \frac{p(n)}{\alpha!(n)}'.
  Recurrence rec_const_coeff(x(n-1)+inhomogeneous_term/product_factor);
  std::vector<Polynomial_Root> new_roots;
  new_roots.push_back(Polynomial_Root(Expr(1), RATIONAL));
  rec_const_coeff.finite_order_p
    = new Finite_Order_Info(1, coefficients, 1);
  rec_const_coeff.set_type(LINEAR_FINITE_ORDER_CONST_COEFF);
  rec_const_coeff.set_inhomogeneous_term(inhomogeneous_term/product_factor);
  rec_const_coeff.set_first_valid_index(first_valid_index); 
  solution = product_factor
    * (rec_const_coeff.solve_constant_coeff_order_1(new_roots)
       + x(first_valid_index));
  return SUCCESS;
}

namespace {
using namespace PURRS;

Expr
substitute_pwr_diff_symbols(const Expr& e,
			    const Symbol& alpha, const Symbol& lambda,
			    const Expr& subs_to_diff) {
  Expr e_substituted;
  if (e.is_a_add()) {
    e_substituted = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_substituted += substitute_pwr_diff_symbols(e.op(i), alpha, lambda,
						   subs_to_diff);
  }
  else if (e.is_a_mul()) {
    e_substituted = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_substituted *= substitute_pwr_diff_symbols(e.op(i), alpha, lambda,
						   subs_to_diff);
  }
  else if (e.is_a_power()) {
    const Expr& base = e.arg(0);
    const Expr& exponent = e.arg(1);
    Number exp;
    if (exponent.is_a_number(exp) && exp.is_negative()) {
      if (base == lambda - alpha)
	return pwr(subs_to_diff, -exp);
      if (base == alpha - lambda)
	return pwr(-subs_to_diff, -exp);
    }
    return pwr(substitute_pwr_diff_symbols(e.arg(0), alpha, lambda,
					   subs_to_diff),
	       substitute_pwr_diff_symbols(e.arg(1), alpha, lambda,
					   subs_to_diff));
  }
  else if (e.is_a_function())
    if (e.nops() == 1)
      return PURRS::apply(e.functor(), substitute_pwr_diff_symbols(e.arg(0),
							    alpha, lambda,
							    subs_to_diff));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int j = 0; j < num_argument; ++j)
	argument[j] = substitute_pwr_diff_symbols(e.arg(j), alpha, lambda,
						  subs_to_diff);
      return PURRS::apply(e.functor(), argument);
    }
  else
    e_substituted = e;
  return e_substituted;
}

void
find_max_exponent_of_lambda(const Expr& e, const Symbol& lambda,
			    unsigned int& max) {
  if (e.is_a_add() || e.is_a_mul())
    for (unsigned int i = e.nops(); i-- > 0; )
      find_max_exponent_of_lambda(e.op(i), lambda, max);
  else if (e.is_a_power()) {
    Number exp;
    if (e.arg(0) == lambda && e.arg(1).is_a_number(exp)) {
      assert(exp.is_positive_integer());
      unsigned int candidate_max = exp.to_unsigned_int();
      if (candidate_max > max)
	max = candidate_max;
    }
    else {
      find_max_exponent_of_lambda(e.arg(0), lambda, max);
      find_max_exponent_of_lambda(e.arg(1), lambda, max);
    }
  }
  else if (e.is_a_function()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      find_max_exponent_of_lambda(e.arg(i), lambda, max);
  }
}

} // anonymous namespace

/*!
  Consider the linear recurrence relation of the second order with
  constant coefficients
  \f[
    x_n = a_1 x_{n-1} + a_2 + p(n),
  \f]
  where \f$ p(n) \f$ is a function defined over the natural numbers.
  If the roots of the characteristic equation \f$ \lambda_1 \f$ and
  \f$ \lambda_2 \f$ are distinct then we know the final formula that
  gives the solution:
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
  the function <CODE>compute_term_about_initial_conditions()</CODE>),
  respectively.
*/
PURRS::Expr
PURRS::Recurrence::
solve_constant_coeff_order_2(Expr& g_n, bool all_distinct,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots) const {
  assert(is_linear_finite_order_const_coeff());
  assert(order() == 2);
  // We search exponentials in `Recurrence::n'.
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `inhomogeneous_term'.
  // In the `i'-th position of the vectors
  // `exp_poly_coeff' and `exp_no_poly_coeff' there are respectively
  // the polynomial part and possibly non polynomial part of the coefficient
  // of the exponential with the base in `i'-th position of `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of `base_of_exps[i]^n'.
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous_term.expand(), Recurrence::n,
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  Expr solution;
  // Calculates the solution of the second order recurrences when
  // the inhomogeneous term is a polynomial or the product of a
  // polynomial and an exponential; return the symbolic solution otherwise.
  if (all_distinct) {
    const Expr& root_1 = roots[0].value();
    const Expr& root_2 = roots[1].value();
    // We use the actual values of the roots of the characteristic equation
    // in order to compute `diff_roots' so that it is possible to simplify it.
    Expr diff_roots = root_1 - root_2;
    diff_roots = blackboard.rewrite(diff_roots).expand();
    D_VAR(diff_roots);
    g_n = (pwr(root_1, Recurrence::n+1) - pwr(root_2, Recurrence::n+1))
      / diff_roots;
    if (!vector_not_all_zero(exp_no_poly_coeff)) {
      Symbol alpha;
      Symbol lambda;
      std::vector<Expr> symbolic_sum_distinct;
      std::vector<Expr> symbolic_sum_no_distinct;
      compute_symbolic_sum(alpha, lambda, roots,
			   base_of_exps, exp_poly_coeff,
			   order() + first_valid_index,
			   symbolic_sum_distinct, symbolic_sum_no_distinct);
      for (unsigned int j = symbolic_sum_distinct.size(); j-- > 0; ) {
	symbolic_sum_no_distinct[j] *= lambda / diff_roots;
	symbolic_sum_distinct[j] *= lambda / diff_roots;
      }
      // Substitute all the occurrences of `(lambda-alpha)^(-k)'
      // with `((lambda+alpha-c)/(-alpha^2+alpha*c+d))^k',
      // where `c' is the coefficient of `x(n-1)',
      // `d' is the coefficient of `x(n-2)', `k' is a positive integer.
      // Notice that GiNaC can rewrite the same expression in different
      // ways: it is possible to have `-(alpha-lambda)^(-k)' and
      // then it must be substituted with
      // `(-(lambda+alpha-c)/(-alpha^2+alpha*c+d))^k'
      // (more details are available in the paper
      // "The Automatic Solution of Recurrence Relations
      // V. Transforming and Simplifying the Solutions").
      if (!vector_not_all_zero(symbolic_sum_no_distinct))
	for (unsigned int j = symbolic_sum_distinct.size(); j-- > 0; ) {
	  const Expr& subs_to_diff = (lambda + alpha - coefficients[1])
	    *pwr(-pwr(alpha, 2) + alpha*coefficients[1] + coefficients[2], -1);
	  symbolic_sum_distinct[j]
	    = substitute_pwr_diff_symbols(symbolic_sum_distinct[j],
					  alpha, lambda, subs_to_diff);
	  // To expand is necessary for finding the powers of `lambda'.
	  symbolic_sum_distinct[j] = symbolic_sum_distinct[j].expand();
#if 0
	  // Find the maximum exponent of `lambda'.
	  unsigned int max_exponent = 1;
	  find_max_exponent_of_lambda(symbolic_sum_distinct[j], lambda,
				      max_exponent);
	  // Substitute every occurrence of powers of lambda (with
	  // positive integer exponent):
	  // `lambda^2 = c*lambda +d' (from the characteristic equation);
	  // `lambda^3 = (c*lambda + d) * lambda
	  //           = ... = c^2*lambda + c*d + d*lambda'
	  // and so forth. In this way we no longer have powers  (with
	  // positive integer exponent) of lambda.

// 	  // First possibility.
// 	  const Expr& expr_to_subs = coefficients[1]*lambda + coefficients[2];
// 	  for (unsigned int i = max_exponent; i > 1; --i)
// 	    symbolic_sum_distinct[j]
// 	      = symbolic_sum_distinct[j].substitute(pwr(lambda, i),
// 						    pwr(lambda, i - 2)
// 						    * expr_to_subs).expand();

	  // Second possibility.
	  std::vector<Expr> powers_of_lambda(max_exponent+1);
	  powers_of_lambda[2] = coefficients[1]*lambda + coefficients[2];
	  for (unsigned i = 3; i <= max_exponent; ++i) {
	    powers_of_lambda[i] = (powers_of_lambda[i-1]*lambda).expand();
	    powers_of_lambda[i]
	      = powers_of_lambda[i].substitute(pwr(lambda,2),
					       powers_of_lambda[2]).expand();
	  }
	  for (unsigned i = 2; i <= max_exponent; ++i)
	    symbolic_sum_distinct[j]
	      = symbolic_sum_distinct[j].substitute(pwr(lambda, i),
						    powers_of_lambda[i]);

// 	  // Third possibility.
// 	  for (unsigned int i = max_exponent; i > 1; --i) {
// 	    D_MSG("");
// 	    D_VAR(symbolic_sum_distinct[j]);
// 	    const Expr& tmp = coefficients[1]*pwr(lambda, i-1)
// 	      + coefficients[2]*pwr(lambda, i-2);
// 	    symbolic_sum_distinct[j]
// 	      = symbolic_sum_distinct[j].substitute(pwr(lambda, i), tmp);
// 	    D_MSG("");
// 	    D_VAR(symbolic_sum_distinct[j]);
// 	  }
#endif
	}
      // Substitutes to the sums in the vector `symbolic_sum_distinct'
      // or `symbolic_sum_no_distinct' the corresponding values of the
      // characteristic equation's roots and of the bases of the
      // possible exponentials and in `solution' put the sum of all
      // sums of the vector after the substitution.
      solution = subs_to_sum_roots_and_bases(alpha, lambda, roots,
					     base_of_exps,
					     symbolic_sum_distinct,
					     symbolic_sum_no_distinct);
      D_VAR(solution);
    }
    else {
      Symbol h;
      // In this function we know the order: 2.
      unsigned int lower = first_valid_index + 2;
      solution = 1 / diff_roots
	* (pwr(root_1, Recurrence::n+1)
	   * PURRS::sum(h, lower, Recurrence::n, pwr(root_1, -h)
			* inhomogeneous_term.substitute(Recurrence::n, h))
	   - (pwr(root_2, Recurrence::n+1)
	      * PURRS::sum(h, lower, Recurrence::n, pwr(root_2, -h)
			   * inhomogeneous_term
			   .substitute(Recurrence::n, h))));
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
      // In this function we know the order: 2.
      unsigned int lower = first_valid_index + 2;
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
  the function <CODE>compute_term_about_initial_conditions()</CODE>),
  respectively.
*/
PURRS::Expr
PURRS::Recurrence::
solve_constant_coeff_order_k(Expr& g_n, bool all_distinct,
			     const std::vector<Number>& coefficients,
			     const std::vector<Polynomial_Root>& roots) const {
  assert(is_linear_finite_order_const_coeff());
  assert(order() > 2);
  // We search exponentials in `n'.
  // The vector `base_of_exps' contains the exponential's bases
  // of all exponentials in `inhomogeneous_term'.
  // In the `i'-th position of the vectors
  // `exp_poly_coeff' and `exp_no_poly_coeff' there are respectively
  // the polynomial part and possibly non polynomial part of the coefficient
  // of the exponential with the base in `i'-th position of `base_of_exp'.
  // `exp_poly_coeff[i] + exp_no_poly_coeff[i]' represents the
  // coefficient of `base_of_exps[i]^n'.
  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(inhomogeneous_term.expand(), Recurrence::n,
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
      Symbol alpha;
      Symbol lambda;
      std::vector<Expr> symbolic_sum_distinct;
      std::vector<Expr> symbolic_sum_no_distinct;
      compute_symbolic_sum(alpha, lambda, roots,
			   base_of_exps, poly_coeff_tot,
			   order() + first_valid_index,
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
      unsigned int lower = first_valid_index + order();
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
      unsigned int lower = first_valid_index + order();
      solution
	= PURRS::sum(h, lower, Recurrence::n,
		     g_n.substitute(Recurrence::n, Recurrence::n - h)
		     * inhomogeneous_term.substitute(Recurrence::n, h));
    }
  return solution;
}

/*!
  Compute a special solution of the linear recurrence
  of finite order with constant coefficients
  \f$ x(n) = homo_rhs + poly * base^n \f$
  where
    `homo_rhs' is the homogeneous part of the rhs of the recurrence
    `poly' is a polynomial in the variable `n'
    `base' is some complex number, which is a root of the characteristic
       polynomial with multiplicity `mult' (possibly zero)
  The solution has the form \f$ p * base^n \f$
  where `p' is a polynomial of degree deg(poly) + mult
*/
PURRS::Expr
compute_special_solution(const Expr& homo_rhs, const index_type order_rec,
			 const Expr& poly, const Expr& base,
			 const unsigned int mult) {
  D_VAR(homo_rhs);
  D_VAR(order_rec);
  D_VAR(poly);
  D_VAR(base);
  D_VAR(mult);

  // Use fast, tailor-made method for summing polynomials when possible
  if ((base == 1) && (homo_rhs == x(Recurrence::n-1))) {
    Expr q = sum_poly_alt(poly, Recurrence::n, Recurrence::n);
    return(q);
  }

  // Build the generic polynomial `q' of the correct degree
  // with unknown coefficients.
  unsigned int deg = poly.degree(Recurrence::n) + mult;
  Expr_List unknowns;
  for (unsigned int i = 0; i < deg + 1; ++i)
    unknowns.append(Symbol());
  Expr q = 0;
  for (unsigned int i = 0; i < deg + 1; ++i)
    q += pwr(Recurrence::n, i) * unknowns.op(i);
  q *= pwr(base, Recurrence::n);
  // At this point, the expression `q' has the desired form poly * base^n.

  Expr diff = x(Recurrence::n) - homo_rhs
    - poly * pwr(base, Recurrence::n);

  // Substitute the expected solution `q' into the original recurrence.
  for (index_type d = 0; d <= order_rec; ++d) {
    const Expr& shifted_poly = q.substitute(Recurrence::n, Recurrence::n - d);
    diff = diff.substitute(x(Recurrence::n - d), shifted_poly);
  }

  // Remove a factor `base^n' from all terms of the expression `diff'
  diff *= pwr(base, -Recurrence::n);
  // After expansion, `diff' is a polynomial with unknown coefficients
  // with the property that the expression `q' built above satisfies
  // the original recurrence if, and only if, `diff' is the zero polynomial.
  diff = diff.expand();

  // FIXME: We need to explicitly simplify diff otherwise it will not be
  // recognized as a polynomial. But remove the unneeded simplifications
  // from below.
  diff = simplify_ex_for_output(diff, false);
  diff = simplify_binomials_factorials_exponentials(diff);
  diff = simplify_logarithm(diff);
  diff = diff.expand();
  D_VAR(diff);

  // Set up a system of `deg + 1' unknowns, forcing all coefficients
  // of the polynomial `diff' to vanish.
  Expr_List equations;
  for (unsigned int i = 0; i < deg + 1; ++i) {
    const Expr& coeff = diff.coeff(Recurrence::n, i);
    equations.prepend(Expr(coeff, 0));
  }

  // Solve the system.
  Expr sol_system = lsolve(equations, unknowns);

  // Compute the actual form of the function `q' by substituting the
  // correct values of the coefficients found above.
  // It is expected that exactly `mult' unknowns, corresponding to the
  // coefficients of degree 0, ..., mult-1, will be undetermined. We
  // can safely substitute 0 for the undetermined unknowns.
  unsigned int removed_unknowns = 0;
  for (unsigned int i = sol_system.nops(); i-- > 0; ) {
    if (sol_system.op(i).rhs() != unknowns.op(i))
      q = q.substitute(unknowns.op(i), sol_system.op(i).op(1));
    else {
      q = q.substitute(unknowns.op(i), 0);
      removed_unknowns++;
    }
  }

  // FIXME: Check that the removed unknowns are the right ones.
  assert(removed_unknowns == mult);

  D_VAR(q);
  return q;
}


/*!
  Solves the linear recurrence of finite order with constant or variable
  coefficients.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::solve_linear_finite_order() const {
  assert(is_linear_finite_order());
  D_VAR(order());
  D_VEC(coefficients(), 1, order());
  D_VAR(inhomogeneous_term);

  Expr solution;
  // We call recurrences of `order zero' special recurrences of the
  // form `x(n) = rhs', where `rhs' contains only functions of `n',
  // parameters and `x(k_1)', ..., `x(k_m)' where `m >= 0' and
  // `k_1', ..., `k_m' are non-negative integers.
  // In this case `*this' is not a proper recurrence and the solution
  // is simply `rhs' (now equal to `inhomogeneous_term').
  if (order() == 0) {
    solution = inhomogeneous_term;
    solution = simplify_ex_for_output(solution, false);
    solution = simplify_binomials_factorials_exponentials(solution);
    solution = simplify_logarithm(solution);
    exact_solution_.set_expression(solution);
    return SUCCESS;
  }

  // FIXME: usare i vettori solo dall'ordine 2 in su.
  // L'ordine 1 gestirlo a parte senza equazione caratteristica.
  // numeric_coefficients()
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

#define OLD_METHOD 0
#if OLD_METHOD
  // `g_n' is defined here because it is necessary in the function
  // `compute_term_about_initial_conditions()' (at the end of function
  // `solve_linear_finite_order()').
  Expr g_n;
  switch (order()) {
  case 1:
    {
      if (is_linear_finite_order_const_coeff())
	// Compute the non-homogeneous part of the solution.
	solution = solve_constant_coeff_order_1(roots);
      else {
	Solver_Status status;
	// Compute the complete solution.
	if ((status = solve_variable_coeff_order_1(coefficients(), solution))
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
      // Compute the non-homogeneous part of the solution.
      solution = solve_constant_coeff_order_2(g_n, all_distinct,
					      num_coefficients, roots);
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
      // Compute the non-homogeneous part of the solution.
      solution = solve_constant_coeff_order_k(g_n, all_distinct,
					      num_coefficients, roots);
    }
    else
      // For the time being, we only solve recurrence relations
      // of order 3 and more only if they have constant coefficients.
      return TOO_COMPLEX;
    break;
  }
  D_MSGVAR("Before calling simplify: ", exact_solution_.expression());
  solution = simplify_ex_for_output(solution, false);
  solution = simplify_binomials_factorials_exponentials(solution);
  solution = simplify_logarithm(solution);

  if (is_linear_finite_order_const_coeff())
    if (order() == 1 )
      // FIXME: per ora non si puo' usare la funzione
      // `compute_term_about_initial_conditions' perche' richiede un
      // vettore di `Number' come `coefficients' e voglio risolvere anche
      // le parametriche (g_n pu' essere posta uguale ad 1 in questo caso).
      // compute_term_about_initial_conditions(g_n, coefficients, solution);
      exact_solution_.set_expression(x(first_valid_index)
				     * pwr(coefficients()[1],
					   n-(first_valid_index-order()+1))
				     + solution);
    else
      exact_solution_.set_expression(compute_term_about_initial_conditions
				     (g_n, num_coefficients, first_valid_index)
				     + solution);
  else
    // In the case of variable coefficients the expression contained in
    // `solution' is already the sum of the homogeneous part and the
    // non-homogeneous part of the solution of the recurrence.
    exact_solution_.set_expression(solution);

  // Resubstitute auxiliary definitions possibly appearing in
  // the solution with their original values.
  //exact_solution_.set_expression(blackboard.rewrite(solution));
  
#else
  /*
  // `g_n' is defined here because it is necessary in the function
  // `compute_term_about_initial_conditions()' (at the end of function
  // `solve_linear_finite_order()').
  Expr g_n;
  */
  if (is_linear_finite_order_const_coeff()) {
#if DEBUG
    for (size_t i = 0; i < roots.size(); ++i)
      std::cerr << i << " - " << roots[i].value()  << " - " << roots[i].multiplicity() << std::endl;
#endif
    Expr homo_rhs = recurrence_rhs - inhomogeneous_term;
    D_MSGVAR("rhs: ", recurrence_rhs);
    D_MSGVAR("Homogeneous term: ", homo_rhs);
    const index_type order_rec = order();
    Expr sol = 0;
    Expr poly;
    Expr base;
    std::vector<Expr> base_of_exps;
    std::vector<Expr> exp_poly_coeff;
    std::vector<Expr> exp_no_poly_coeff;
    exp_poly_decomposition(inhomogeneous_term.expand(), Recurrence::n,
			   base_of_exps, exp_poly_coeff, exp_no_poly_coeff);
    D_VEC(base_of_exps, 0, base_of_exps.size()-1);
    D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
    D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

    /*
    // FIXME: If we are able to deal with non-zero extra terms in some cases,
    // do not give up.
    if (vector_not_all_zero(exp_no_poly_coeff))
      return TOO_COMPLEX;
    */

    // FIXME: This is a copy-and-paste of the old resolution method. It is
    // supposed to be used, temporarily, until the new method solves all the
    // recurrences the old method attempted to solve. At the moment only the
    // recurrences whose inhomogeneous term cannot be expressed as sum of
    // polynomials*exponentials require the old method.

    if (vector_not_all_zero(exp_no_poly_coeff)) {

      // This is the old method. It will be eventually be removed from here.

      // `g_n' is defined here because it is necessary in the function
      // `compute_term_about_initial_conditions()' (at the end of function
      // `solve_linear_finite_order()').
      Expr g_n;
      switch (order()) {
      case 1:
	{
	  if (is_linear_finite_order_const_coeff())
	    // Compute the non-homogeneous part of the solution.
	    solution = solve_constant_coeff_order_1(roots);
	  else {
	    Solver_Status status;
	    // Compute the complete solution.
	    if ((status = solve_variable_coeff_order_1(coefficients(), solution))
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
	  // Compute the non-homogeneous part of the solution.
	  solution = solve_constant_coeff_order_2(g_n, all_distinct,
						  num_coefficients, roots);
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
	  // Compute the non-homogeneous part of the solution.
	  solution = solve_constant_coeff_order_k(g_n, all_distinct,
						  num_coefficients, roots);
	}
	else
	  // For the time being, we only solve recurrence relations
	  // of order 3 and more only if they have constant coefficients.
	  return TOO_COMPLEX;
	break;
      }
      D_MSGVAR("Before calling simplify: ", exact_solution_.expression());
      solution = simplify_ex_for_output(solution, false);
      solution = simplify_binomials_factorials_exponentials(solution);
      solution = simplify_logarithm(solution);
      
      if (is_linear_finite_order_const_coeff())
	if (order() == 1 )
	  // FIXME: per ora non si puo' usare la funzione
	  // `compute_term_about_initial_conditions' perche' richiede un
	  // vettore di `Number' come `coefficients' e voglio risolvere anche
	  // le parametriche (g_n pu' essere posta uguale ad 1 in questo caso).
	  // compute_term_about_initial_conditions(g_n, coefficients, solution);
	  exact_solution_.set_expression(x(first_valid_index)
					 * pwr(coefficients()[1],
					       n-(first_valid_index-order()+1))
					 + solution);
	else {
	  /*
	  Expr tmp = compute_term_about_initial_conditions(g_n, num_coefficients, first_valid_index) + solution;
	  tmp = blackboard.rewrite(tmp.expand());
	  D_MSGVAR("Solution with old method: ", tmp);
	  exact_solution_.set_expression(tmp);
	  */
	  exact_solution_.set_expression(blackboard.rewrite
					 (compute_term_about_initial_conditions
					  (g_n, num_coefficients, first_valid_index) + solution).expand());

	}
      else
	// In the case of variable coefficients the expression contained in
	// `solution' is already the sum of the homogeneous part and the
	// non-homogeneous part of the solution of the recurrence.
	exact_solution_.set_expression(solution);
      
      // Resubstitute auxiliary definitions possibly appearing in
      // the solution with their original values.
      //exact_solution_.set_expression(blackboard.rewrite(solution));
      
      return SUCCESS;
    }


    // If there are non-rational roots, substitute them with symbols for
    // efficiency.
    substitute_non_rational_roots(*this, roots);

    for (size_t j = 0; j < base_of_exps.size(); ++j) {
      poly = exp_poly_coeff[j];
      base = base_of_exps[j];
      size_t mult = 0;
      for (size_t k = roots.size(); k-- > 0; )
	if (blackboard.rewrite(roots[k].value()) == base)
	  mult = roots[k].multiplicity();
    
      sol += compute_special_solution(homo_rhs, order_rec, poly, base, mult);
    }
    solution = sol; 
  }
  else 
    if (order()==1) {
      Solver_Status status;
      // Compute the complete solution.
      if ((status = solve_variable_coeff_order_1(coefficients(), solution))
	  != SUCCESS)
	return status;
    }
    else
      return TOO_COMPLEX;

  D_MSGVAR("Before calling simplify: ", exact_solution_.expression());
  solution = simplify_ex_for_output(solution, false);
  solution = simplify_binomials_factorials_exponentials(solution);
  solution = simplify_logarithm(solution);

  if (is_linear_finite_order_const_coeff())
    if (order() == 1 )
      // FIXME: This could be integrated in the method below if 
      // parameters are not a problem.
      exact_solution_.set_expression(blackboard.rewrite((x(first_valid_index) - solution.substitute(n, first_valid_index))
							* pwr(coefficients()[1],
							      n-first_valid_index)
							+ solution));
    else {
      unsigned int num_roots = roots.size();
      Expr_List unknowns;
      std::vector<unsigned int> root_for;
      std::vector<unsigned int> multiplicity;
      Expr general_solution = solution;
      for (unsigned int i = 0; i < num_roots; ++i) {
	unsigned int root_multiplicity = roots[i].multiplicity();
	for (unsigned int j = 1; j <= root_multiplicity; ++j) {
	  unknowns.append(Symbol());
	  root_for.push_back(i);
	  multiplicity.push_back(j);
	  general_solution += unknowns.op(unknowns.nops()-1) * pwr(roots[i].value(),n) * pwr(n, j-1);
	}
      }
      D_MSGVAR("General solution: ", general_solution);
      unsigned int num_unknowns = unknowns.nops();
      Expr_List equations;
      for (unsigned int i = 0; i < num_unknowns; ++i) {
	//	Expr tmp = x(first_valid_index + i) - general_solution.substitute(n, first_valid_index + i);
	//	D_MSGVAR("Equation: ", tmp);
	equations.prepend(Expr(x(first_valid_index + i) - general_solution.substitute(n, first_valid_index + i), 0));
	D_MSGVAR("Equation: ", equations.op(0));
      }

      // Solve the system.
      Expr sol_system = lsolve(equations, unknowns);

      // FIXME: Would throwing an exception be better?
      if (sol_system.nops() == 0)
	return TOO_COMPLEX;

      // Substitute the correct values of the coefficients found above.
      for (unsigned int i = sol_system.nops(); i-- > 0; ) {
	D_MSGVAR("Solution: ", sol_system.op(i).rhs());
	general_solution = general_solution.substitute(unknowns.op(i), sol_system.op(i).rhs());
     }

      D_MSGVAR("Complete solution: ", general_solution);
     
      // Resubstitute auxiliary definitions possibly appearing in
      // the solution with their original values.
      exact_solution_.set_expression(blackboard.rewrite(general_solution));
    }
  else
    // In the case of variable coefficients the expression contained in
    // `solution' is already the sum of the homogeneous part and the
    // non-homogeneous part of the solution of the recurrence.
    exact_solution_.set_expression(blackboard.rewrite(solution));

#endif

  // Only for the output.
  if (exact_solution_.expression().is_a_add()) {
    Expr_List conditions;
    for (index_type i = order(); i-- > 0; )
      conditions.append(x(first_valid_index + i));
    // FIXME: `collect' throws an exception if the object to collect has
    // non-integer exponent. 
    exact_solution_
      .set_expression(exact_solution_.expression().collect(conditions));
  }

  return SUCCESS;
}
