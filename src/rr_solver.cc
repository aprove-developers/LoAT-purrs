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

#include "globals.hh"
#include "util.hh"
#include "sum_poly.hh"
#include "alg_eq_solver.hh"
#include <climits>

// TEMPORARY
#include <iostream>

using namespace GiNaC;

#define NOISY 1

static GExpr
get_binding(const GList& l, unsigned wild_index) {
  assert(wild_index < l.nops());
  assert(l.op(wild_index).info(GiNaC::info_flags::relation_equal));
  assert(l.op(wild_index).lhs() == GiNaC::wild(wild_index));
  return l.op(wild_index).rhs();
}

static bool
get_constant_decrement(const GExpr& e, const GSymbol& n, GNumber& decrement) {
  static GExpr n_plus_d = n + GiNaC::wild(0);
  GList substitution;
  if (match(e, n_plus_d, substitution)) {
    GExpr d = get_binding(substitution, 0);
    if (GiNaC::is_a<GiNaC::numeric>(d)) {
      decrement = -GiNaC::ex_to<GiNaC::numeric>(d);
      return true;
    }
  }
  return false;
}

/*!
  Transforms expressions of the form \f$(a^x)^b)\f$ in the expression \p p
  into \f$a^(b*x)\f$.
*/
static GExpr
mul_exponents_exp(GExpr& p, const GList lst_of_exp) {
  GExpr q = 0;

  for (size_t i = lst_of_exp.nops(); i-- > 0; ) {
    GExpr tmp = lst_of_exp.op(i);
    tmp = tmp.subs(pow(pow(wild(0), wild(1)), wild(2)) == 
		   pow(wild(0), wild(1)*wild(2)));
    // Finds coefficient of exponential 'lst_of_exp.op(i)'.
    GList substitution;
    clear(substitution);
    if (match(p, lst_of_exp.op(i)*wild(0), substitution) ||
	match(p, lst_of_exp.op(i)*wild(0) + wild(1), substitution)) {
      GExpr coeff = get_binding(substitution, 0);
      q += tmp*coeff;
      p -= coeff*lst_of_exp.op(i);
    }
    // The exponential's coefficient is 1.
    else {
      q += tmp;
      p -= lst_of_exp.op(i);
    }
  }
  // Now p does not contain other exponential of the form (a^x)^b. 
  q += p;

  return q;
}

/*!
  Transforms expressions of the form \f$a^(b*x)\f$ in the expression \p p
  into \f$(a^b)^x\f$.
*/
static GExpr
split_exp(GExpr& p, const GSymbol& x, const GList lst_of_exp) {
  GExpr q = 0;

  for (size_t i = lst_of_exp.nops(); i-- > 0; ) {
    GExpr tmp = lst_of_exp.op(i);
    tmp = tmp.subs(pow(wild(0),x*wild(1)) == pow(pow(wild(0),wild(1)),x));
    // Finds coefficient of exponential 'lst_of_exp.op(i)'.
    GList substitution;
    clear(substitution);
    if (match(p, lst_of_exp.op(i)*wild(0), substitution) ||
	match(p, lst_of_exp.op(i)*wild(0) + wild(1), substitution)) {
      GExpr coeff = get_binding(substitution, 0);
      q += tmp*coeff;
      p -= coeff*lst_of_exp.op(i);
    }
    // The exponential's coefficient is 1.
    else {
      q += tmp;
      p -= lst_of_exp.op(i);
    }
  }
  // Now p does not contain other exponential of the form a^(b*x). 
  q += p;

  return q;
}

/*!
  Transforms expressions of the form \f$a^x*b^x\f$ in the expression \p p
  into \f$(a^b)^x\f$.
*/
static GExpr 
union_exp(GExpr& p, const GSymbol& x, const GList lst_of_exp,
	  const unsigned n) {
  GExpr q = 0;

  for (size_t i = lst_of_exp.nops(); i-- > 0; ) {
    GExpr tmp = lst_of_exp.op(i);
    if (n == 1) {
      tmp = tmp.subs(pow(wild(0),x)*pow(wild(1),x) == pow(wild(0)*wild(1),x));
      q += tmp;
    } 
    else {
      tmp = tmp.subs(pow(wild(0),x)*pow(wild(1),x)*wild(2) == 
		     pow(wild(0)*wild(1),x));
      // Finds coefficient of exponential 'lst_of_exp.op(i)'.  
      GList l;
      GExpr coeff;
      if (p.find(lst_of_exp.op(i) * wild(), l)) {
	coeff = l.op(0);
	coeff = coeff.subs(wild(0)*pow(wild(1), x)*pow(wild(2), x) == wild(0));
      }
      q += tmp*coeff;
    }
    p -= lst_of_exp.op(i);
  }
  // Now p does not contain other exponential of the form a^x*b^x.
  q += p;

  return q;
}

static void
transform_exponentials(GExpr& e, const GSymbol& n) {
  GList lst_of_exp;

  // Transforms (a^n)^b into a^(b*n).
  static GExpr a_x_b = pow(pow(wild(0), wild(1)), wild(2));
  while (e.find(a_x_b, lst_of_exp)) {
    e = mul_exponents_exp(e, lst_of_exp);
    clear(lst_of_exp);
  }

  // Transforms a^(b*n+c) into (a^b)^n*a^c.
  static GExpr a_bn = pow(wild(0), n*wild(1));
  static GExpr a_n_c = pow(wild(0), n+wild(1));
  static GExpr a_bn_c = pow(wild(0), n*wild(1)+wild(2));
  if (e.find(a_bn_c, lst_of_exp) || e.find(a_n_c, lst_of_exp)) 
    e = e.expand();
  if (e.find(a_bn, lst_of_exp))
    e = split_exp(e, n, lst_of_exp);

  // Transforms a^n*b^n into (a*b)^n.
  static GExpr a_n_b_n = pow(wild(0), n)*pow(wild(1), n);
  static GExpr c_a_n_b_n = pow(wild(0), n)*pow(wild(1), n)*wild(2);
  while (e.find(a_n_b_n, lst_of_exp) || e.find(c_a_n_b_n, lst_of_exp)) { 
    if (e.find(a_n_b_n, lst_of_exp))
	e = union_exp(e, n, lst_of_exp, 1);
    else
      e = union_exp(e, n, lst_of_exp, 2);
    clear(lst_of_exp);
  }
}

static GMatrix
decomposition_inhomogeneous_term(const GExpr& e, const GSymbol& n);

static bool
solution_1_poly_times_exponentials(const GSymbol& x_0, const GSymbol& n,
				   const GMatrix& decomposition,
				   const std::vector<GNumber>& coefficients,
				   GExpr& solution);

static bool
solution_2_poly_times_exponentials(const GSymbol& x_0, const GSymbol& x_1,
				   const GSymbol& n,
				   const GMatrix& decomposition,
				   const std::vector<GNumber>& coefficients,
				   GExpr& solution);

bool
solve(const GExpr& rhs, const GSymbol& n) {
  static GExpr x_i = x(GiNaC::wild(0));
  static GExpr x_i_plus_r = x_i + GiNaC::wild(1);
  static GExpr a_times_x_i = GiNaC::wild(1)*x_i;
  static GExpr a_times_x_i_plus_r = a_times_x_i + GiNaC::wild(2);

  static GList substitution;

  int order = -1;
  std::vector<GNumber> coefficients;
  GExpr e = rhs;
  bool failed = false;
  do {
    GExpr i;
    GExpr a;
    // The following matches are attempted starting from the most common,
    // then the second most common and so forth.
    if (clear(substitution), match(e, x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = get_binding(substitution, 1);
    }
    else if (clear(substitution), match(e, a_times_x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = get_binding(substitution, 2);
    }
    else if (clear(substitution), match(e, a_times_x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = 0;
    }
    else if (clear(substitution), match(e, x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = 0;
    }
    else
      break;

    GNumber decrement;
    if (!get_constant_decrement(i, n, decrement)) {
      failed = true;
      break;
    }
#if NOISY
    std::cout << "decrement = " << decrement << std::endl;
#endif
    if (!decrement.is_integer()
	|| decrement < 0
	|| decrement >= coefficients.max_size()) {
      failed = true;
      break;
    }
    if (!GiNaC::is_a<GiNaC::numeric>(a)) {
      failed = true;
      break;
    }
    GNumber coefficient = GiNaC::ex_to<GiNaC::numeric>(a);
    // FIXME: turn this assertion into something more appropriate.
    assert(decrement >= LONG_MIN && decrement <= LONG_MAX);
    unsigned long index = decrement.to_long();
    if (order < 0 || index > unsigned(order))
      order = index;
    // FIXME: why 'coefficients' is created with one useless element, 
    // the first, that is always zero?
    if (index > coefficients.size())
      coefficients.insert(coefficients.end(),
			  index - coefficients.size(),
			  GNumber(0));
    if (index == coefficients.size())
      coefficients.insert(coefficients.end(), coefficient);
    else
      coefficients[index] += coefficient;
  } while (e != 0);

  // See if what is left is the inhomogeneous term,
  // i.e., if all the occurrences of `x(e)' are such that
  // `e' is numeric.
  GList occurrences;
  if (e.find(x_i, occurrences)) {
    for (unsigned i = 0, n = occurrences.nops(); i < n; ++i)
      if (!is_a<numeric>(occurrences.op(i))) {
	failed = true;
	break;
      }
  }

  if (failed)
    return false;
#if NOISY
  std::cout << "Order = " << order << std::endl;
  std::cout << "Coefficients = ";
  for (int i = 1; i <= order; ++i)
    std::cout << coefficients[i] << " ";
  std::cout << std::endl;
  std::cout << "Inhomogeneous term = " << e << std::endl;
#endif
  // The factors of the form a^(bn+c) (a,b,c numeric) must be transformed
  // into (a^b)^n*a^c. GiNaC tranforms only a^(bn+c) in a^c*a^(bn) but not
  // a^(bn) into (a^b)^n.
  // Besides transform factor of the form (a^n)^b into a a^(b*n)
  // (GiNaC do not make this).
  transform_exponentials(e, n);

  // Now certainly the inhomogeneous term contains only exponential
  // with exponent n (or a power of n).
  // 'decomposition' is a matrix with two rows and a number of columns
  // which is at most the number of exponentials in the inhomogeneous term
  // plus one.
  // In every column there is a exponential in the first row and its
  // coefficient in the second row. In the last column there is the
  // constant exponential with its coefficients. 
  GMatrix decomposition = decomposition_inhomogeneous_term(e, n);
#if NOISY
  std::cout << "Inhomogeneous term's decomposition"
	    << decomposition << std::endl;
#endif
  // Calculates the number of columns of the matrix 'decomposition'.
  unsigned num_columns = decomposition.cols();
  GExpr solution;
  switch (order) {
  case 1:
    {
      bool poly_times_exp = true;
      for (size_t i = 0; i < num_columns; ++i) {
	GExpr exponential = decomposition(0, i);
	GExpr coeff_of_exp = decomposition(1, i);
	// Check if the inhomogeneous term is a polynomial or a the
	// product of a polynomial and an exponential.
	if (is_a<power>(exponential) &&
	    !GiNaC::is_a<GiNaC::numeric>(exponential.op(0)))
	  poly_times_exp = false;
	if (!is_a<power>(exponential) && exponential != 1) {
	  poly_times_exp = false;
	}
	if (!coeff_of_exp.info(info_flags::polynomial))
	  poly_times_exp = false;
      }
      if (poly_times_exp) {
	// Calculates the solution of the first order recurrences when
	// the inhomogeneous term is a polynomial or the product of a
	// polynomial and an exponential.      
	GSymbol x_0("x_0");
	solution_1_poly_times_exponentials(x_0, n, decomposition,
					   coefficients, solution);
      }
      else 
	throw ("PURRS error: at the moment the recurrence "
	       "relation is solved only when the inhomogeneous term "
	       "is polynomials or product of exponentials times "
	       "polynomials.");
      break;
    }
  case 2:
    {
      bool poly_times_exp = true;
      for (size_t i = 0; i < num_columns; ++i) {
	GExpr exponential = decomposition(0, i);
	GExpr coeff_of_exp = decomposition(1, i);
	// Check if the inhomogeneous term is a polynomial or a the
	// product of a polynomial and an exponential.
	if (is_a<power>(exponential) &&
	    !GiNaC::is_a<GiNaC::numeric>(exponential.op(0)))
	  poly_times_exp = false;
	if (!is_a<power>(exponential) && exponential != 1)
	  poly_times_exp = false;
	if (!coeff_of_exp.info(info_flags::polynomial))
	  poly_times_exp = false;
      }
      if (poly_times_exp) {
	// Calculates the solution of the second order recurrences when
	// the inhomogeneous term is a polynomial or the product of a
	// polynomial and an exponential.      
	GSymbol x_0("x_0"), x_1("x_1");
	solution_2_poly_times_exponentials(x_0, x_1, n, decomposition,
					   coefficients, solution);
      }	
      else 
	throw ("PURRS error: at the moment the recurrence "
	       "relation is solved only when the inhomogeneous term "
	       "is polynomials or product of exponentials times "
	       "polynomials."); 
      break;
    }
  default:
    {
      throw ("PURRS error: order too large"); 
    } 
  } 
  
#if NOISY
  std::cout << "Solution  " << solution << std::endl << std::endl;
#endif
  
  return true;
}

/*!
  This function makes a matrix with two rows and a number of columns
  does not exceed the number of exponentials in the inhomogeneous term
  plus one.
  The function gives the decomposition 
  \f$ e(n) = sum_{i=0}^k \alpha_i^j \cdot p(n)_i \f$ with
  - \f$ \alpha_i \ne \alpha_j \f$ if \f$ i \ne j \f$
  - \p p does not contains exponentials.
  It returns the matrix whose \f$ i\f$-th column contains 
  \f$ \alpha_i^n \f$ and \f$ p(n)_i \f$
  for \f$ i = 1, \ldots, k \f$.
*/
static GMatrix
decomposition_inhomogeneous_term(const GExpr& e, const GSymbol& n) {
  GExpr p, q;
  GList(lst_of_exp);

  p = e;
  GExpr pattern = pow(wild(), n);
  p.find(pattern, lst_of_exp);
  p = p.collect(lst_of_exp);

  if (lst_of_exp.nops() == 0)
    // In this case there are no axponentials.
    return GMatrix (2, 1, lst(1, p));
  else {
    // In this case there is at least one exponential in 'p'
    GList row_exp;
    GList row_coeff;
    q = p;
    for (size_t i = lst_of_exp.nops(); i-- > 0; ) {
      GList addendum;
      // 'addendum' has always only one element because we have collected
      // 'p' with respect to the exponentials.
      p.find(lst_of_exp.op(i) * wild(), addendum);
      row_exp.append(lst_of_exp.op(i));
      if (addendum.nops() != 0) {
	// The exponential's coefficient is not 1.
	GExpr coeff = addendum.op(0);
	coeff = coeff.subs(wild(1)*pow(wild(2), n) == wild(1));
	row_coeff.append(coeff);
	q -= addendum.op(0);
      }
      else {
	// The exponential's coefficient is 1.
	row_coeff.append(1);
	q -= lst_of_exp.op(i);
      }
    }
    // Now 'q' does not contains any exponentials or product of exponential 
    // times other expressions.
    unsigned int num_columns;
    if (!q.is_zero()) {
      row_exp.append(1);
      row_coeff.append(q);
      num_columns = lst_of_exp.nops() + 1;
    }
    else
      num_columns = lst_of_exp.nops();

    // The next 'for' is necessary for to use the constructor matrices
    //      matrix::matrix(unsigned r, unsigned c, const lst& l);
    unsigned row_coeff_nops = row_coeff.nops();
    for (size_t i = 0; i < row_coeff_nops; ++i)
      row_exp.append(row_coeff.op(i));

    GMatrix terms_divided(2, num_columns, row_exp);
      
    return terms_divided;
  }
}

/*!
  Calculates the solution of the first order recurrence
  \f$x(n) = \alpha * x(n-1) + p(n)\f$ with \p x_0 as
  initial condition.
  The solution is given by the closed formula
  \f$x(n) = \alpha^n * x_0 + sum_{k=1}^n \alpha^(n-k)*p(k)\f$.
  \p decomposition is a matrix that contains a decomposition of
  the inhomogeneous term \p p(n) and \p coefficients a vector that
  contains the coefficients of the recurrence.
 */
static bool
solution_1_poly_times_exponentials(const GSymbol& x_0, const GSymbol& n,
				   const GMatrix& decomposition,
				   const std::vector<GNumber>& coefficients,
				   GExpr& solution_tot) {
  // The closed formula of the solution can be rewritten as
  // x(n) = \alpha^n * ( x_0 - p(0) + sum_{k=0}^n \alpha^(-k)*p(k) ).
  solution_tot = x_0;
  GSymbol k("k");
  // Calculates the number of columns of the matrix 'decomposition'.
  unsigned num_columns = decomposition.cols();
  GExpr alpha_n = pow(coefficients[1],n);
  for (size_t i = 0; i < num_columns; ++i) {
    GExpr solution = 0;
    GExpr exponential = decomposition(0, i);
    GExpr coeff_of_exp = decomposition(1, i);
    GExpr coeff_of_exp_k = coeff_of_exp.subs(n == k);
    if (is_a<power>(exponential)) 
      solution = sum_poly_times_exponentials(coeff_of_exp_k, k, n, 
					     exponential.op(0) *
					     pow(coefficients[1], -1));
    else
      // This is the case of the constant exponential.
      solution = sum_poly_times_exponentials(coeff_of_exp_k, k, n,
					     pow(coefficients[1], -1)); 
    // 'sum_poly_times_exponentials' calculates the sum from 0 while
    // we want to start from 1. 
    solution -= coeff_of_exp.subs(n == 0);	
    solution_tot += solution;
  }
  solution_tot *= alpha_n;
  solution_tot = solution_tot.expand();
  transform_exponentials(solution_tot, n);

  return true;
} 

static void
build_characteristic_equation(GExpr& p, const GSymbol x,
			      const std::vector<GNumber>& coefficients) {
  // The function 'find_roots' wants only integer coefficients
  // for to calculates the roots of the equation.
  GNumber c = 0;
  unsigned coefficients_size = coefficients.size();
  for (unsigned i = 1; i < coefficients_size; ++i)
    if (!coefficients[i].is_integer() && coefficients[i].denom() > c) 
      c = coefficients[i].denom();
  if (c != 0) {
    // Build the vector 'int_coefficients' with the elements of
    // 'coefficients' multiplies for a factor 'c' that is the
    // bigger among the denominator of the elements of 'coefficients'.
    std::vector<GNumber> int_coefficients(coefficients);
    for (unsigned i = 1; i < coefficients_size; ++i)
	int_coefficients[i] *= c;
    p = c * pow(x, 2) + (-int_coefficients[1])*x + (-int_coefficients[2]); 
  }
  else
    p = pow(x, 2) + (-coefficients[1])*x + (-coefficients[2]);
}

/*!
  Calculates the solution of the second order recurrence
  \f$x(n) = \alpha * x(n-1) + \beta * x(n-2) + p(n)\f$, (\f$ beta \neq 0\f$),
  with \p x_0 and \p x_1 as initials conditions.
  We consider the homogeneous recurrence
  \f$g(n) = \alpha * g(n-1) + \beta * g(n-2)\f$
  with \f$g_0 = 1\f$ and \f$g_1 = \alpha\f$ as initials conditions for the
  fundamental solution.
  The complete solution is given by the formula
  \f$x(n) = g(n-1) * x_1 + \beta g(n-2) * x_0 +
         \sum_{k=2}^n g(n-k) * p(k)\f$,   for \f$n \ge 2\f$. 
  \p decomposition is a matrix that contains a decomposition of
  the inhomogeneous term \p p(n) and \p coefficients a vector that
  contains the coefficients of the recurrence.
 */
static bool
solution_2_poly_times_exponentials(const GSymbol& x_0, const GSymbol& x_1,
				   const GSymbol& n, 
				   const GMatrix& decomposition,
				   const std::vector<GNumber>& coefficients,
				   GExpr& solution_tot) {
    solution_tot = 0;
    // Calculates the number of columns of the matrix.
    unsigned num_columns = decomposition.cols();

    // We first solve the equation with g(n) by means of the
    // characteristic equation wich, by definition, is the polynomial
    // equation x^2 + \alpha * x + \beta = 0.
    // Build the characteristic equation.
    GExpr p = 0;
    static GSymbol x("x");
    build_characteristic_equation(p, x, coefficients);
#if NOISY
    std::cout << "p " << p << std::endl;
#endif
    std::vector<Polynomial_Root> roots;
    bool all_distinct;
    if (!find_roots(p, x, roots, all_distinct))
      return false;
    // We proceed in two different ways in according to the value of
    // the boolean 'all_distinct', i. e., if all roots of the
    // characteristic equation are distinct or not.
    if (all_distinct) {
      // In this case g(n) has the explicit expression
      // g(n) = \frac{\lambda_1^{n+1} - \lambda_2^{n+1}}
      // {\lambda_1 - \lambda_2}.   
      // Hence for n \ge 2 we have
      // x_n = g(n-1) * x_1 + \beta g(n-2) * x_0 +
      // \frac{\lambda_1^{n+1}}{\lambda_1 - \lambda_2}
      //   \sum_{k=2}^n \lambda_1^{-k} \, p(k) -
      //\frac{\lambda_2^{n+1}}{\lambda_1 - \lambda_2}
      //   \sum_{k=2}^n \lambda_2^{-k} \, p(k).
      
      // Construct g(n).
      GExpr root_1 = roots[0].value();
      GExpr root_2 = roots[1].value();
#if NOISY
      std::cout << "root_1 " << root_1 << std::endl;
      std::cout << "root_2 " << root_2 << std::endl;
#endif
      GExpr diff_roots = root_1 - root_2;
      GExpr g_n = (pow(root_1, n+1) - pow(root_2, n+1)) / diff_roots;
#if NOISY
      std::cout << "g_n " << g_n << std::endl;
#endif
      // FIXME: this is necessary for a bug of GiNaC about subs
      // (in the next release of GiNaC it must be fix and then
      // the GSymbol k will be useless).
      GSymbol k("k");
      GExpr g_n_1 = g_n.subs(n == k-1);
      GExpr g_n_2 = g_n.subs(n == k-2);
      g_n_1 = g_n_1.subs(k == n);
      g_n_2 = g_n_2.subs(k == n);

      solution_tot = x_1 * g_n_1 + x_0 * g_n_2 * coefficients[2];
      for (size_t i = 0; i < num_columns; ++i) {
	GExpr solution_1 = 0;
	GExpr solution_2 = 0;
	GExpr exponential = decomposition(0, i);
	GExpr coeff_of_exp = decomposition(1, i);
	GExpr coeff_of_exp_k = coeff_of_exp.subs(n == k);
	if (is_a<power>(exponential)) {
	  solution_1 = sum_poly_times_exponentials(coeff_of_exp_k, k, n, 
						   exponential.op(0)*
						   pow(root_1, -1));
	  solution_2 = sum_poly_times_exponentials(coeff_of_exp_k, k, n, 
						   exponential.op(0)*
						   pow(root_2, -1));
	}     
	else {
	  // This is the case of the constant exponential.
	  solution_1 = sum_poly_times_exponentials(coeff_of_exp_k, k, n, 
						   pow(root_1, -1));
	  solution_2 = sum_poly_times_exponentials(coeff_of_exp_k, k, n, 
						   pow(root_2, -1));
	}
	// 'sum_poly_times_exponentials' calculates the sum from 0 while
	// we want to start from 2.
	solution_1 -= coeff_of_exp.subs(n == 0) + coeff_of_exp.subs(n == 1) *
	  (1/root_1);
	solution_2 -= coeff_of_exp.subs(n == 0) + coeff_of_exp.subs(n == 1) *
	  (1/root_2);
	solution_1 *= pow(root_1, n+1) / diff_roots;
	solution_2 *= pow(root_2, n+1) / diff_roots;
	solution_tot += solution_1 - solution_2;
      }
      solution_tot = solution_tot.expand();
      solution_tot = solution_tot.collect(lst(x_0, x_1));
      transform_exponentials(solution_tot, n);
    }
    else {
      // 'all_distinct == false': the characteristic equation
      // x^2 + \alpha * x + \beta = 0 has a double root.
      assert(roots[0].multiplicity() == 2);
      
      GSymbol a("a"), b("b");
      GExpr root = roots[0].value();
#if NOISY
      std::cout << "root " << root << std::endl;
#endif
      
      if (decomposition(1, 0).is_zero()) {
	// The binary recurrence is homogeneous.
	// The solution of the homogeneous recurrence is of the form
	// x_n = (a + b * n) \lambda^n where
	// \lambda is the root with multiplicity 2.
	// In this case we must to solve a linear system:
	//             a = x_0
	//             (a + b) * \lambda = x_1.
	solution_tot = (a + b * n) * pow(root, n);
	// Solved the system with the inverse matrix'method.
	GMatrix vars(2, 2, lst(1, 0, root, root));
	GMatrix rhs(2, 1, lst(x_0, x_1));
	GMatrix sol(2, 1);
	sol = vars.inverse();
	sol = sol.mul(rhs);
#if NOISY
	std::cout << "matrix solution " << sol << std::endl;
#endif
	solution_tot = solution_tot.subs(a == sol(0,0));
	solution_tot = solution_tot.subs(b == sol(1,0));
      }
      else {
	// The binary recurrence is non-homogeneous.
	// We set
	// g(n) = \alpha * g(n-1) + \beta * g(n-2)
	// with g_0 = 1 and g_1 = \alpha as initials conditions for the
	// fundamental solution.
	// g_n = (a + b * n) \lambda^n where
	// \lambda is the root with multiplicity 2.
	// In this case we must to solve a linear system:
	//             a = 1
	//             (a + b) * \lambda = \alpha.
	GExpr g_n = (a + b * n) * pow(root, n);
	// Solved the system with the inverse matrix'method.
	GMatrix vars(2, 2, lst(1, 0, root, root));
	GMatrix rhs(2, 1, lst(1, coefficients[1]));
	GMatrix sol(2, 1);
	sol = vars.inverse();
	sol = sol.mul(rhs);
#if NOISY
	std::cout << "matrix solution " << sol << std::endl;
#endif
	g_n = g_n.subs(a == sol(0,0));
	g_n = g_n.subs(b == sol(1,0));
	
	// FIXME: this is necessary for a bug of GiNaC about subs
	// (in the next release of GiNaC it must be fix and then
	// the GSymbol k will be useless).
	GSymbol k("k");
	GExpr g_n_1 = g_n.subs(n == k-1);
	GExpr g_n_2 = g_n.subs(n == k-2);
	g_n_1 = g_n_1.subs(k == n);
	g_n_2 = g_n_2.subs(k == n);
	
	solution_tot = x_1 * g_n_1 + x_0 * g_n_2 * coefficients[2];
	for (size_t i = 0; i < num_columns; ++i) {
	  GExpr solution = 0;
	  GExpr exponential = decomposition(0, i);
	  GExpr coeff_of_exp = decomposition(1, i);
	  GExpr coeff_of_exp_k = coeff_of_exp.subs(n == k);
	  GExpr g_n_k = g_n.subs(n == n - k);
	  if (is_a<power>(exponential))
	    solution = sum_poly_times_exponentials(coeff_of_exp_k * g_n_k,
						   k, n, exponential.op(0));
	  else
	    // This is the case of the constant exponential.
	    solution = sum_poly_times_exponentials(coeff_of_exp_k * g_n_k,
						   k, n, 1);
	  // 'sum_poly_times_exponentials' calculates the sum from 0 while
	  // we want to start from 2.
	  solution = solution.expand();
	  solution -= (coeff_of_exp_k * g_n_k).subs(k == 0) +
	              (coeff_of_exp_k * g_n_k).subs(k == 1) * (1/root);
	  solution = solution.expand();
	  solution_tot += solution;
	}
	solution_tot = solution_tot.expand();
	solution_tot = solution_tot.collect(lst(x_0, x_1));
	transform_exponentials(solution_tot, n);
      }
    }

  return true;
}
