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
#include "simplify.hh"
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

/*!
  Returns <CODE>true</CODE> if and only if the <CODE>GExpr</CODE> \p e
  is of the form \f$ n + d \f$ with \f$ d \f$ a <CODE>GiNaC::numeric</CODE>;
  in this case \p decrement contains the opposite of \f$ d \f$.
*/
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
  Give a vector of <CODE>GNumber</CODE> with a number of elements like
  the order \p order of the recurrence relation, builds the characteristic
  equation in the variable \p x.
  Whence, if we have the recurrence relation of the order \f$ k \f$
  \f$ x_n = a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$
  then the coefficients \f$ a_j \f$ are in the vector \p coefficients and the
  function returns the characteristic equation \p
  \f$ x^k = a_1 * x^{k-1} + \dotsb + a^{k-1} * x + a_k \f$.
*/
static GExpr
build_characteristic_equation(const int order, const GSymbol& x,
			      const std::vector<GNumber>& coefficients) {
  GExpr p;
  p = 0;
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
    for (int i = 0; i < order; ++i)
      p += pow(x, i) * (-int_coefficients[order - i]);
    p += c * pow(x, order);
  }
  else {
    for (int i = 0; i < order; ++i)
      p += pow(x, i) * (-coefficients[order - i]);
    p += pow(x, order);
  }
  return p;
}

/*!
  Returns <CODE>true</CODE> if and only if the inhomogeneous term (that is
  decomposed in a matrix \p decomposition is a polynomial or the
  product of a polynomial and an exponential; <CODE>false</CODE> otherwise.
*/
static bool
check_poly_times_exponential(const GMatrix& decomposition) {
  // Calculates the number of columns of the matrix 'decomposition'.
  unsigned num_columns = decomposition.cols();
  bool poly_times_exp = true;
  for (size_t i = 0; i < num_columns; ++i)
    if (!decomposition(2, i).is_zero())
      poly_times_exp = false;
  return poly_times_exp;
}

/*!
  We consider the recurrence relation's inhomogeneous term 
  \f$ p(n) = sum_{i=0}^decomposition.cols() q(n)_i alpha_i^{n} \f$ and the
  vector \p roots of the characteristic equation's roots and we call
  \p \lambda the generic root.
  This function fills the two <CODE>vector<GExpr></CODE>, with dimension
  equal to \p decomposition.cols(), \p symbolic_sum_distinct and
  \p symbolic_sum_no_distinct with two different sum:
  for \f$ j = 0, \dotsc, roots.size() \f$ and for
  \f$ i = 0, \dotsc, decomposition.cols() \f$
  -  if \f$ alpha_i \neq \lambda_j \f$ then
     \f$ symbolic_sum_distinct[i] = f_n_i(\alpha / \lambda)
           = \lambda^n * sum_{k=order}^n {\alpha / \lambda}^k \cdot q(k)_i \f$;
  -  if \f$ alpha_i = \lambda_j \f$ then
     \f$ symbolic_sum_no_distinct[i] = f_n_i(\lambda)
           = \lambda^n * sum_{k=order}^n q(k)_i \f$.
*/
static void
compute_symbolic_sum(const int order, const GSymbol& n, const GSymbol& alpha,
		     const GSymbol& lambda,
		     const std::vector<Polynomial_Root> roots,
		     const GMatrix& decomposition,
		     std::vector<GExpr>& symbolic_sum_distinct,
		     std::vector<GExpr>& symbolic_sum_no_distinct) {
  unsigned num_columns = decomposition.cols();
  for (size_t i = num_columns; i-- > 0; ) {
    bool distinct = true;
    for (size_t j = roots.size(); j-- > 0; ) {
      if (is_a<power>(decomposition(0,i))) {
	if (roots[j].value().is_equal(decomposition(0,i).op(0)))
	  distinct = false;
      }
      else
	if (roots[j].value().is_equal(decomposition(0,i)))
	  distinct = false;
      GSymbol k("k");
      GExpr q_k = decomposition(1, i).subs(n == k);
      // The root is deifferent from the exponential's base.
      if (distinct) {
	GSymbol x("x");
	symbolic_sum_distinct[i] = sum_poly_times_exponentials(q_k, k, n, x);
	// 'sum_poly_times_exponentials' calculates the sum from 0 while
	// we want to start from 'order'.
	symbolic_sum_distinct[i] -= q_k.subs(k == 0);
	for (int j = 1; j < order; ++j)
	  symbolic_sum_distinct[i] -= q_k.subs(k == j) * pow(lambda, -1);
	
	symbolic_sum_distinct[i] =
	  symbolic_sum_distinct[i].subs(x == alpha/lambda);
	symbolic_sum_distinct[i] *= pow(lambda, n);
	symbolic_sum_distinct[i] =
	  simplify_on_output_ex(symbolic_sum_distinct[i].expand());
      }
      // The root is equal to the exponential's base.
      else {
	symbolic_sum_no_distinct[i] = sum_poly_times_exponentials(q_k, k, n,
								  1);
	// 'sum_poly_times_exponentials' calculates the sum from 0 while
	// we want to start from 'order'.
	symbolic_sum_no_distinct[i] -= q_k.subs(k == 0);
	for (int j = 1; j < order; ++j)
	  symbolic_sum_no_distinct[i] -= q_k.subs(k == j) * pow(lambda, -1);
	symbolic_sum_no_distinct[i] *= pow(lambda, n);
	symbolic_sum_no_distinct[i] =
	  simplify_on_output_ex(symbolic_sum_no_distinct[i].expand());
      }
    }     
  }
}

static GMatrix
decomposition_inhomogeneous_term(const GExpr& e, const GSymbol& n);

static void
subs_to_sum_roots_and_bases(const GSymbol& alpha, const GSymbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    const GMatrix& decomposition,
			    std::vector<GExpr>& symbolic_sum_distinct,
			    std::vector<GExpr>& symbolic_sum_no_distinct,
			    GExpr& solution);

static GExpr
order_2_sol_roots_no_distinct(const GSymbol& n, const GMatrix& decomposition,
			      const std::vector<GExpr> initial_conditions,
			      const std::vector<GNumber>& coefficients,
			      const std::vector<Polynomial_Root>& roots);

/*!
...
*/
bool
solve(const GExpr& rhs, const GSymbol& n, GExpr& solution) {
  static GExpr x_i = x(GiNaC::wild(0));
  static GExpr x_i_plus_r = x_i + GiNaC::wild(1);
  static GExpr a_times_x_i = GiNaC::wild(1)*x_i;
  static GExpr a_times_x_i_plus_r = a_times_x_i + GiNaC::wild(2);

  int order = -1;
  std::vector<GExpr> coefficients;
  GExpr e = rhs;

  // Special case: 'e' is only a function in n or a constant.
  // If there is not in 'e' 'x(argument)' or there is but with
  // 'argument' that not contains the variable 'n', then the solution
  // is obviously 'e'. 
  GList occurrences;
  bool finished = true;
  if (e.find(x_i, occurrences)) {
    int occurrences_nops = occurrences.nops();
    if (occurrences_nops != 0)
      for (int i = 0; i < occurrences_nops; ++i) {
	GExpr argument = occurrences.op(i).subs(x(wild(0)) == wild(0));
	if (argument.has(n)) {
	  finished = false;
	  break;
	}
      }
  }
  if (finished) {
    solution = e; 
#if NOISY 
    std::cout << "Solution " << solution << std::endl << std::endl;
#endif
    return true;
  }

  static GList substitution;
  bool failed = false;
  do {
    GExpr i;
    GExpr a;
    // The following matches are attempted starting from the most common,
    // then the second most common and so forth.
    // The check 'if (!i.has(n))' is necessary because otherwise do not
    // accept 'x(i)' with 'i' numeric.
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
    // For the moment the coefficients of recurrence relation
    // must be constants, i. e., it does not contains the variable n.
    if (a.has(n))
      throw ("PURRS error: for the moment we only solve recurrence "
	     "relations with constant coefficients. ");
    GExpr coefficient = a;

    // FIXME: turn this assertion into something more appropriate.
    assert(decrement >= LONG_MIN && decrement <= LONG_MAX);
    unsigned long index = decrement.to_long();
    if (order < 0 || index > unsigned(order))
      order = index;

    // The vector 'coefficients' contains in the i-th position
    // the coefficient of x(n-i); in the first position and in the
    // positions corresponding to x(n-i) absent, there is '0'.
    if (index > coefficients.size())
      coefficients.insert(coefficients.end(),
			  index - coefficients.size(),
			  GNumber(0));
    if (index == coefficients.size())
      coefficients.insert(coefficients.end(), coefficient);
    else
      coefficients[index] += coefficient;
  } while (e != 0);
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

  // Simplifies expressions, in particular rewrites nested powers.
  e = simplify_on_input_ex(e);

  // 'decomposition' is a matrix with three rows and a number of columns
  // which is at most the number of exponentials in the inhomogeneous term
  // plus one.
  // In every column there is a exponential in the first row and its
  // coefficient in the second and third row: in the second row there is
  // the polynomial part of the coefficient and in the third row there is
  // non polynomial part of the coefficients. If the coefficient is not
  // polynomial then in the second row there is zero and, similarly, in the
  // third row there is zero if the coefficient is a polynomial.
  // In the last column there is the constant exponential with its
  // coefficients. 
  GMatrix decomposition = decomposition_inhomogeneous_term(e, n);
#if NOISY
  std::cout << "Inhomogeneous term's decomposition"
	    << decomposition << std::endl;
#endif
  // Creates the vector of initials conditions.
  // FIXME: fare controllo sulle condizioni iniziali
  // (numero minore dell'ordine).
  std::vector<GExpr> initial_conditions(order);
  for (int i = 0; i < order; ++i)
    initial_conditions[i] = x(i);

  // Computes the characteristic equation and its roots.
  GExpr characteristic_eq;
  GSymbol y("y");
  std::vector<Polynomial_Root> roots;
  bool all_distinct;
  // FIXME: il seguente if sull'ordine e' temporaneo perche' cosi' si
  // riescono a fare le parametriche del primo ordine almeno.
  // Temporaneo fino a che 'find_roots' accettera' i parametri anche
  // per equazioni di grado superiore al primo. 
  std::vector<GNumber> num_coefficients(order+1);
  if (order == 1) {
    characteristic_eq = y - coefficients[1];
    roots.push_back(coefficients[1]);
  } 
  else {
    // Check that the vector 'coefficients' does not contains
    // parameters and in this case uses a vector of GNumber.
    for (int i = 1; i <= order; ++i)
      if (is_a<numeric>(coefficients[i]))
	num_coefficients[i] = GiNaC::ex_to<GiNaC::numeric>(coefficients[i]);
      else
	throw("PURRS error: the second order recurrence relations does not "
	      "support the case of the coefficients with "
	      "parameters. ");
    characteristic_eq = build_characteristic_equation(order, y,
						      num_coefficients);
    if (!find_roots(characteristic_eq, y, roots, all_distinct))
      return false;
  }
#if NOISY
  std::cout << "characteristic equation " << characteristic_eq << std::endl;
  for (size_t i = roots.size(); i-- > 0; )
    std::cout << "root_" << i << ": " << roots[i].value() << "    ";
  std::cout << std::endl;
#endif

  GSymbol alpha("alpha");
  GSymbol lambda("lambda");
  std::vector<GExpr> symbolic_sum_distinct(decomposition.cols());
  std::vector<GExpr> symbolic_sum_no_distinct(decomposition.cols());
  if (all_distinct)
    // Fills the two vector 'symbolic_sum_distinct' and
    // 'symbolic_sum_no_distinct' with two different symbolic sum
    // A SECONDA CHE the roots are equal or not to the exponential's bases.
    compute_symbolic_sum(order, n, alpha, lambda, roots, decomposition,
			 symbolic_sum_distinct, symbolic_sum_no_distinct);
  switch (order) {
  case 1:
    {
      // Calculates the solution of the first order recurrences when
      // the inhomogeneous term is a polynomial or the product of a
      // polynomial and an exponential.
      if (check_poly_times_exponential(decomposition)) {
	// Substitutes to the sums in the vectors 'symbolic_sum_distinct' or
	// 'symbolic_sum_no_distinct' the corresponding values of the
	// characteristic equation's roots and of the bases of the eventual
	// exponentials and in 'solution' put the sum of all sums of the
	// vectors after the substitution.
	subs_to_sum_roots_and_bases(alpha, lambda, roots, decomposition,
				    symbolic_sum_distinct,
				    symbolic_sum_no_distinct, solution);
	solution += initial_conditions[0] * pow(coefficients[1] ,n);
#if NOISY
	std::cout << "Solutions before calling simplify: " << std::endl;
	std::cout << solution << std::endl;
#endif        
	solution = simplify_on_output_ex(solution.expand());
      }
      else 
	throw ("PURRS error: for the moment the recurrence "
	       "relation is solved only when the inhomogeneous term "
	       "is polynomials or product of exponentials times "
	       "polynomials.");
      break;
    }
  case 2:
    {
      if (all_distinct)
	// Calculates the solution of the second order recurrences when
	// the inhomogeneous term is a polynomial or the product of a
	// polynomial and an exponential.
	if (check_poly_times_exponential(decomposition)) {
	  GExpr root_1 = roots[0].value();
	  GExpr root_2 = roots[1].value();
	  GExpr diff_roots = root_1 - root_2;
	  for (size_t j = symbolic_sum_distinct.size(); j-- > 0; ) {
	    symbolic_sum_no_distinct[j] *= lambda / diff_roots;
	    symbolic_sum_distinct[j] *= lambda / diff_roots;
	  }
	  // Substitutes to the sums in the vector 'symbolic_sum_distinct' or
	  // 'symbolic_sum_no_distinct' the corresponding values of the
	  // characteristic equation's roots and of the bases of the eventual
	  // exponentials and in 'solution' put the sum of all sums of the
	  // vector after the substitution.
	  subs_to_sum_roots_and_bases(alpha, lambda, roots, decomposition,
				      symbolic_sum_distinct,
				      symbolic_sum_no_distinct, solution);
	  GExpr g_n = (pow(root_1, n+1) - pow(root_2, n+1)) / diff_roots;
	  // FIXME: forse conviene semplificare g_n
#if NOISY
	  std::cout << "g_n " << g_n << std::endl;
#endif
	  GExpr g_n_1 = g_n.subs(n == n - 1);
	  GExpr g_n_2 = g_n.subs(n == n - 2);
	  solution += initial_conditions[1] * g_n_1 + 
	    initial_conditions[0] * g_n_2 * num_coefficients[2];
	  solution = simplify_on_output_ex(solution.expand());
	  solution = solution.collect(lst(initial_conditions[0],
					  initial_conditions[1]));
	}
	else 
	  throw ("PURRS error: for the moment the recurrence "
		 "relation is only solved when the inhomogeneous term "
		 "is polynomials or product of exponentials times "
		 "polynomials."); 
      else {
	// The characteristic equation
	// x^2 + a_1 * x + a_2 = 0 has a double root.
	assert(roots[0].multiplicity() == 2);      
	solution = order_2_sol_roots_no_distinct(n, decomposition,
						 initial_conditions,
						 num_coefficients, roots);
      }
      break;
    }
  default:
    throw ("PURRS error: order too large"); 
  }
  return true;
}

/*!
  This function makes a matrix with three rows and a number of columns
  does not exceed the number of exponentials in the inhomogeneous term
  plus one.
  The function gives the decomposition 
  \f$ e(n) = sum_{i=0}^k \alpha_i^j \cdot p(n)_i \f$ with
  - \f$ \alpha_i \ne \alpha_j \f$ if \f$ i \ne j \f$
  - \p p does not contains exponentials.
  It returns the matrix whose \f$ i\f$-th column, for \f$ i = 1, \ldots, k \f$,
  contains \f$ \alpha_i^n \f$ in the first row and \f$ p(n)_i \f$
  in the second and third row: the polynomial part of \f$ p(n)_i \f$ in the
  second row and the non polynomial part in the third row.
*/
static GMatrix
decomposition_inhomogeneous_term(const GExpr& e, const GSymbol& n) {
  GExpr p, q;
  GList(lst_of_exp);
  p = e;
  GExpr pattern = pow(wild(), n);
  p.find(pattern, lst_of_exp);
  // 'lst_of_exp' contains all exponentials.
  p = p.collect(lst_of_exp);

  if (lst_of_exp.nops() == 0) {
    // In this case there are not any exponentials.
    GExpr p_poly; 
    GExpr p_no_poly;
    // In 'p' there are not nested powers.
    assign_poly_part_and_no_poly_part(p, n, p_poly, p_no_poly);
    return GMatrix (3, 1, lst(1, p_poly, p_no_poly));
  }
  else {
    // In this case there is at least one exponential in 'p'.
    // 'row_exp' will contains all exponentials; 'row_coeff_poly' and
    // 'row_coeff_no_poly' will contain polynomial part and non
    // polynomial part of exponential's coefficients, rispectively.
    GList row_exp;
    GList row_coeff_poly;
    GList row_coeff_no_poly;
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
	GExpr coeff_poly;
	GExpr coeff_no_poly;
	// In 'coeff' there are not nested powers.
	assign_poly_part_and_no_poly_part(coeff, n, coeff_poly, coeff_no_poly);
	row_coeff_poly.append(coeff_poly);
	row_coeff_no_poly.append(coeff_no_poly);
	q -= addendum.op(0);
      }
      else {
	// The exponential's coefficient is 1.
	row_coeff_poly.append(1);
	row_coeff_no_poly.append(0);
	q -= lst_of_exp.op(i);
      }
    }
    // Now 'q' does not contains any exponentials or product of exponential 
    // times other expressions.
    unsigned int num_columns;
    if (!q.is_zero()) {
      row_exp.append(1);
      GExpr q_poly;
      GExpr q_no_poly;
      assign_poly_part_and_no_poly_part(q, n, q_poly, q_no_poly);
      row_coeff_poly.append(q_poly);
      row_coeff_no_poly.append(q_no_poly);

      num_columns = lst_of_exp.nops() + 1;
    }
    else
      num_columns = lst_of_exp.nops();
 
    // The next 'for' is necessary for to use the constructor matrices
    //      matrix::matrix(unsigned r, unsigned c, const lst& l);
    unsigned row_coeff_poly_nops = row_coeff_poly.nops();
    for (size_t i = 0; i < row_coeff_poly_nops; ++i)
      row_exp.append(row_coeff_poly.op(i));
    for (size_t i = 0; i < row_coeff_poly_nops; ++i)
      row_exp.append(row_coeff_no_poly.op(i));
    
    GMatrix terms_divided(3, num_columns, row_exp);
    
    return terms_divided;
  }
}

/*!
  Give a vector \p symbolic_sum that contains all the symbolic sums of the
  inhomogeneous term's terms that are polynomial or the product of a
  polynomial and an exponential, this function substitutes to the sums the
  corresponding values of the characteristic equation's roots and of the bases
  of the eventual exponentials.
  Returns a <CODE>GExpr</CODE> \p solution with the sum of all sums of the
  vector.
 */
static void
subs_to_sum_roots_and_bases(const GSymbol& alpha, const GSymbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    const GMatrix& decomposition,
			    std::vector<GExpr>& symbolic_sum_distinct,
			    std::vector<GExpr>& symbolic_sum_no_distinct,
			    GExpr& solution) {
  solution = 0;
  for (size_t i = roots.size(); i-- > 0; )
    for (size_t j = symbolic_sum_distinct.size(); j-- > 0; ) {
      GExpr base_exp;
      if (is_a<power>(decomposition(0, j)))
	base_exp = decomposition(0, j).op(0);
      else
	base_exp = decomposition(0, j);
      GExpr tmp;
      if (base_exp.is_equal(roots[i].value()))
      tmp = symbolic_sum_no_distinct[j].subs(lst(alpha == base_exp,
						 lambda == roots[i].value()));
      else
	tmp = symbolic_sum_distinct[j].subs(lst(alpha == base_exp,
						lambda == roots[i].value()));
      solution += tmp * pow(-1, i);
    }
}

static GExpr
order_2_sol_roots_no_distinct(const GSymbol& n, const GMatrix& decomposition,
			      const std::vector<GExpr> initial_conditions,
			      const std::vector<GNumber>& coefficients,
			      const std::vector<Polynomial_Root>& roots) {
  GExpr solution_tot;
  // Calculates the number of columns of the matrix.
  unsigned num_columns = decomposition.cols();

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
    // In this case we must solve a linear system:
    //             a = x(0)
    //             (a + b) * \lambda = x(1).
    solution_tot = (a + b * n) * pow(root, n);
    // Solved the system with the inverse matrix'method.
    GMatrix vars(2, 2, lst(1, 0, root, root));
    GMatrix rhs(2, 1, lst(initial_conditions[0], initial_conditions[1]));
    GMatrix sol(2, 1);
    sol = vars.inverse();
    sol = sol.mul(rhs);
    
    solution_tot = solution_tot.subs(lst(a == sol(0,0), b == sol(1,0)));
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
    g_n = g_n.subs(lst(a == sol(0,0), b == sol(1,0)));
    GExpr g_n_1 = g_n.subs(n == n - 1);
    GExpr g_n_2 = g_n.subs(n == n - 2);
	
    solution_tot = initial_conditions[1]*g_n_1 + 
                   initial_conditions[0]*g_n_2*coefficients[2];
    std::cout << "SOL TOT " << solution_tot << std::endl;

    for (size_t i = 0; i < num_columns; ++i) {
      GExpr solution = 0;
      GExpr exponential = decomposition(0, i);
      GExpr coeff_of_exp = decomposition(1, i);
      GSymbol k("k");
      GExpr coeff_of_exp_k = coeff_of_exp.subs(n == k);
      GExpr g_n_k = g_n.subs(n == n - k);
      // In this case g_n_k always contains root^(n-k).
      // We pass to the function 'sum_poly_exponentials' only
      // (1/root)^k that it is multiplied for exponential.op(0).
      // Hence we pass only the polynomial part of g_n with
      // the substitution n == n - k.
      GExpr poly_g_n = g_n.subs(wild(0)*power(wild(1), n) == wild(0));
      poly_g_n = poly_g_n.subs(n == n - k);
      if (is_a<power>(exponential)) {
	solution = sum_poly_times_exponentials(coeff_of_exp_k * poly_g_n,
					       k, n, exponential.op(0) *
					       pow(root, -1));
	// 'sum_poly_times_exponentials' calculates the sum from 0 while
	// we want to start from 2.
	solution -= (coeff_of_exp_k * poly_g_n).subs(k == 0) +
	  (coeff_of_exp_k * poly_g_n).subs(k == 1) * 
	  (1/root) * exponential.op(0);
      }
      else {
	// This is the case of the constant exponential.
	solution = sum_poly_times_exponentials(coeff_of_exp_k * poly_g_n,
					       k, n, pow(root, -1));
	// 'sum_poly_times_exponentials' calculates the sum from 0 while
	// we want to start from 2.
	solution -= (coeff_of_exp_k * poly_g_n).subs(k == 0) +
	  (coeff_of_exp_k * poly_g_n).subs(k == 1) * (1/root);
      }
      // We have passed to the function 'sum_poly_times_exponentials'
      // only (1/root)^k then now we must multiply for (root)^n.
      solution *= power(root, n);
      solution = solution.expand();
      solution_tot += solution;
    }
    solution_tot = simplify_on_output_ex(solution_tot.expand());
    solution_tot = solution_tot.collect(lst(initial_conditions[0],
					    initial_conditions[1]));
  }
  return solution_tot;
}
