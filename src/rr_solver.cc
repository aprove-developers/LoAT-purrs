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
#include "simplify.hh"
#include "alg_eq_solver.hh"
#include <climits>

// TEMPORARY
#include <iostream>

using namespace GiNaC;

/*!
  Returns <CODE>true</CODE> if and only if \p e is of the form
  \f$ n + d \f$ with \f$ d \f$ a <CODE>GiNaC::numeric</CODE>;
  in this case the opposite of \f$ d \f$ is assigned to \p decrement.
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
  Given a vector of <CODE>GNumber</CODE> with a number of elements like
  the order \p order of the recurrence relation, builds the characteristic
  equation in the variable \p x.
  Therefore, if we have the order \f$ k \f$ recurrence relation
  \f$ x_n = a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$,
  the coefficients \f$ a_j \f$ are in the vector \p coefficients and the
  function returns the characteristic equation
  \f$ x^k - ( a_1 * x^{k-1} + \dotsb + a^{k-1} * x + a_k ) = 0 \f$.
*/
static GExpr
build_characteristic_equation(const int order, const GSymbol& x,
			      const std::vector<GNumber>& coefficients) {
  GExpr p;
  p = 0;
  unsigned coefficients_size = coefficients.size();
  
  // We know taht the coefficients are 'numerics', but we want only
  // rationals.
  for (unsigned i = 1; i < coefficients_size; ++i)
    if (!coefficients[i].is_rational())
      throw("PURRS error: today the algebraic equation solver works\n"
            "only with integer coefficients.\n"
            "Please come back tomorrow.");

  // The function 'find_roots' wants only integer coefficients
  // for to calculates the roots of the equation.
  std::vector<GNumber> denominators;
  // Find the least common multiple among the denominators of the
  // rational elements of 'coefficients'.
  for (unsigned i = 1; i < coefficients_size; ++i)
    if (!coefficients[i].is_integer())
      denominators.push_back(coefficients[i].denom());
  if (denominators.size() != 0) {
    GNumber least_com_mul = lcm(denominators);
    // Build the vector 'int_coefficients' with the elements of
    // 'coefficients' multiplies for the least common multiple
    // 'least_com_mul'.
    std::vector<GNumber> int_coefficients(coefficients);
    for (unsigned i = 1; i < coefficients_size; ++i)
      int_coefficients[i] *= least_com_mul;
    for (int i = 0; i < order; ++i)
      p += pow(x, i) * (-int_coefficients[order - i]);
    p += least_com_mul * pow(x, order);
  }
  else {
    for (int i = 0; i < order; ++i)
      p += pow(x, i) * (-coefficients[order - i]);
    p += pow(x, order);
  }
  return p;

}

/*!
  Returns <CODE>true</CODE> if and only if the inhomogeneous term
  is a polynomial or the product of a polynomial and an exponential;
  <CODE>false</CODE> otherwise.
  The vector \p exp_no_poly_coeff contains the non polynomial part
  (if it exists) of the inhomogeneous term.
*/
static bool
check_poly_times_exponential(const std::vector<GExpr>& exp_no_poly_coeff) {
  for (unsigned i = exp_no_poly_coeff.size(); i-- > 0; )
    if (!exp_no_poly_coeff[i].is_zero())
      return false;
  return true;
}

/*!
  We consider
  - the recurrence relation's inhomogeneous term when is in the form
    polynomials times exponentials:
    \f$ e(n) = \sum_{i=0}^k \alpha_i^{n} p_i(n) \f$ (where \f$ k \f$ is
    the number of exponentials);
  - the vector \p roots of the characteristic equation's roots. 
  We call \f$ \lambda \f$ the generic root.

  This function fills the two vectors of <CODE>GExpr</CODE>
  \p symbolic_sum_distinct and \p symbolic_sum_no_distinct, with dimension
  equal to \p base_of_exps.size(), with two different sums:
  for \f$ j = 0, \dotsc, roots.size() \f$ and for
  \f$ i = 0, \dotsc, base_of_exps.size() \f$
  -  if \f$ \alpha_i \neq \lambda_j \f$ then
     \p symbolic_sum_distinct[i] = \f$ f_i(\alpha_i / \lambda)
        = \lambda^n * \sum_{k=order}^n (\alpha_i / \lambda)^k \cdot q_i(k) \f$;
  -  if \f$ \alpha_i = \lambda_j \f$ then
     \p symbolic_sum_no_distinct[i] = \f$ f_i(\lambda)
        = \lambda^n * \sum_{k=order}^n q_i(k) \f$.
*/
static void
compute_symbolic_sum(const int order, const GSymbol& n,
		     const GSymbol& alpha, const GSymbol& lambda,
		     const std::vector<Polynomial_Root> roots,
		     std::vector<GExpr>& base_of_exps,
		     std::vector<GExpr>& exp_poly_coeff,
		     std::vector<GExpr>& symbolic_sum_distinct,
		     std::vector<GExpr>& symbolic_sum_no_distinct) {
  for (size_t i = base_of_exps.size(); i-- > 0; ) {
    bool distinct = true;
    for (size_t j = roots.size(); j-- > 0; ) {
      if (roots[j].value().is_equal(base_of_exps[i]))
	distinct = false;

      GSymbol k("k");
      GExpr q_k = exp_poly_coeff[i].subs(n == k);

      // The root is different from the exponential's base.
      if (distinct) {
	GSymbol x("x");
	symbolic_sum_distinct[i] = sum_poly_times_exponentials(q_k, k, n,
                                                               x);
	// 'sum_poly_times_exponentials' calculates the sum from 0 while
        // we want to start from 'order'.
        symbolic_sum_distinct[i] -= q_k.subs(k == 0);
        for (int j = 1; j < order; ++j)
          symbolic_sum_distinct[i] -= q_k.subs(k == j) * pow(lambda, -1);
        
        symbolic_sum_distinct[i] =
          symbolic_sum_distinct[i].subs(x == alpha/lambda);
        symbolic_sum_distinct[i] *= pow(lambda, n);
        symbolic_sum_distinct[i] =
          simplify_on_output_ex(symbolic_sum_distinct[i].expand(), n, false);
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
          simplify_on_output_ex(symbolic_sum_no_distinct[i].expand(), n,
                                false);
      }
    }
  }
}

static void
exp_poly_decomposition(const GExpr& e, const GSymbol& n,
		       std::vector<GExpr>& alpha,
		       std::vector<GExpr>& p,
		       std::vector<GExpr>& q);

static void
subs_to_sum_roots_and_bases(const GSymbol& alpha, const GSymbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    std::vector<GExpr>& base_of_exps,
			    std::vector<GExpr>& symbolic_sum_distinct,
			    std::vector<GExpr>& symbolic_sum_no_distinct,
			    GExpr& solution);

static void
add_initial_conditions(GExpr& g_n, const GSymbol& n,
		       std::vector<GNumber>& coefficients,
		       const std::vector<GExpr> initial_conditions,
		       GExpr& solution);

static GExpr
order_2_sol_roots_no_distinct(const GSymbol& n,
			      const std::vector<GExpr> base_of_exps,
			      const std::vector<GExpr> exp_poly_coeff,
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
  // 'expand()' is necessary to simplify expressions, to find coefficients
  // and to decompose the inhomogeneous term of the recurrence relation.  
  GExpr e = rhs.expand();

  // Special case: 'e' is only a function in \p n or a constant.
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

    D_VAR(solution);

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
    // accept 'x(i)' with 'i' numeric in a general recurrence relation
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
    if (!get_constant_decrement(i, n, decrement)) {
      failed = true;
      break;
    }

    // FIXME: fare controllo sulle condizioni iniziali, cioe' verificare
    // che se c'e' 'x(i)' con 'i' numerico allora 'i' sia un numero
    // minore dell'ordine.


    D_VAR(decrement);

    if (!decrement.is_integer()
	|| decrement < 0
	|| decrement >= coefficients.max_size()) {
      failed = true;
      break;
    }
    // For the moment the coefficients of recurrence relation
    // must be constants, i. e., it does not contains the variable n.
    if (a.has(n))
      throw ("PURRS error: today we only solve recurrence "
	     "relations with constant coefficients.\n"
             "Please come back tomorrow.");
    GExpr coefficient = a;

    // FIXME: turn this assertion into something more appropriate.
    assert(decrement >= LONG_MIN && decrement <= LONG_MAX);
    unsigned long index = decrement.to_long();
    if (order < 0 || index > unsigned(order))
      order = index;

    // The vector 'coefficients' contains in the i-th position
    // the coefficient of x(n-i); in the first position and in the
    // positions corresponding to x(n-i) absent there is '0'.
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


  D_VAR(order);
  D_VEC(coefficients, 1, order);
  D_MSGVAR("Inhomogeneous term: ", e);


  // Simplifies expanded expressions, in particular rewrites nested powers.
  e = simplify_on_input_ex(e, n, true);

  // We search exponentials in the variable 'n', so the expression 'e'
  // must be expanded because otherwise the function not recognizes the
  // the exponentials in the variable 'n' (i.e. 2^(n-1) is not considerated
  // while 1/2*2^n, obtained from previous by 'expand()', yes).
  // In every column there is a exponential in the first row and its
  // coefficient in the second and third row: in the second row there is
  // the polynomial part of the coefficient and in the third row there is
  // non polynomial part of the coefficients. If the coefficient is not
  // polynomial then in the second row there is zero and, similarly, in the
  // third row there is zero if the coefficient is a polynomial.
  // In the last column there is the constant exponential with its
  // coefficients.
  std::vector<GExpr> base_of_exps;
  std::vector<GExpr> exp_poly_coeff;
  std::vector<GExpr> exp_no_poly_coeff;
  exp_poly_decomposition(e, n, base_of_exps, exp_poly_coeff,
			 exp_no_poly_coeff);

  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  // Creates the vector of initials conditions.
  std::vector<GExpr> initial_conditions(order);
  for (int i = 0; i < order; ++i)
    initial_conditions[i] = x(i);

  // Computes the characteristic equation and its roots.
  GExpr characteristic_eq;
  GSymbol y("y");
  std::vector<Polynomial_Root> roots;
  bool all_distinct = true;
 
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
	throw("PURRS error: today the second order recurrence relations\n"
	      "does not support parametric coefficients.\n"
	      "Please come back tomorrow.");
    characteristic_eq = build_characteristic_equation(order, y,
						      num_coefficients);
    if (!find_roots(characteristic_eq, y, roots, all_distinct))
      return false;
  }


  D_VAR(characteristic_eq);
  D_VEC(roots, 0, roots.size()-1);
  D_MSG("");


  GSymbol alpha("alpha");
  GSymbol lambda("lambda");
  std::vector<GExpr> symbolic_sum_distinct(base_of_exps.size());
  std::vector<GExpr> symbolic_sum_no_distinct(base_of_exps.size());

  if (all_distinct)
    // Fills the two vector 'symbolic_sum_distinct' and
    // 'symbolic_sum_no_distinct' with two different symbolic sum
    // in according to the roots are equal or not to the exponential's bases.
    compute_symbolic_sum(order, n, alpha, lambda, roots,
			 base_of_exps, exp_poly_coeff,
			 symbolic_sum_distinct, symbolic_sum_no_distinct);
  switch (order) {
  case 1:
    {
      // Calculates the solution of the first order recurrences when
      // the inhomogeneous term is a polynomial or the product of a
      // polynomial and an exponential.
      if (check_poly_times_exponential(exp_no_poly_coeff)) {
	// Substitutes to the sums in the vectors 'symbolic_sum_distinct' or
	// 'symbolic_sum_no_distinct' the corresponding values of the
	// characteristic equation's roots and of the bases of the eventual
	// exponentials and in 'solution' put the sum of all sums of the
	// vectors after the substitution.
	subs_to_sum_roots_and_bases(alpha, lambda, roots, base_of_exps,
				    symbolic_sum_distinct,
				    symbolic_sum_no_distinct, solution);
	// FIXME: per ora non si puo' usare la funzione
	// 'add_initial_conditions' perche' richiede un vettore di
	// 'GNumber' come 'coefficients' e voglio risolvere anche le
	// parametriche. 
//  	GExpr g_n = pow(roots[0].value(), n);
//  	add_initial_conditions(g_n, n, coefficients,
//  			       initial_conditions, solution);
	solution += initial_conditions[0] * pow(roots[0].value() ,n);

	D_MSGVAR("Before calling simplify: ", solution);

	solution = simplify_on_output_ex(solution.expand(), n, false);

	//int r = verify_solution(solution, order, e,coefficients);
      }
      else 
	throw ("PURRS error: today we only allow inhomogeneous terms\n"
	       "in the form of polynomials or product of exponentials\n"
	       "and polynomials.\n"
	       "Please come back tomorrow.");
      break;
    }
  case 2:
    {
      // There are two roots with multiplicity equal to 1.
      if (all_distinct)
	// Calculates the solution of the second order recurrences when
	// the inhomogeneous term is a polynomial or the product of a
	// polynomial and an exponential.
	if (check_poly_times_exponential(exp_no_poly_coeff)) {
	  GExpr root_1 = roots[0].value();
	  GExpr root_2 = roots[1].value();
	  GExpr diff_roots = root_1 - root_2;
	  for (size_t j = symbolic_sum_distinct.size(); j-- > 0; ) {
	    symbolic_sum_no_distinct[j] *= lambda / diff_roots;
	    symbolic_sum_distinct[j] *= lambda / diff_roots;
	  }
	  // Substitutes to the sums in the vector 'symbolic_sum_distinct'
	  // or 'symbolic_sum_no_distinct' the corresponding values of the
	  // characteristic equation's roots and of the bases of the
	  // eventual exponentials and in 'solution' put the sum of all
	  // sums of the vector after the substitution.
	  subs_to_sum_roots_and_bases(alpha, lambda, roots, base_of_exps,
				      symbolic_sum_distinct,
				      symbolic_sum_no_distinct, solution);
	  GExpr g_n = (pow(root_1, n+1) - pow(root_2, n+1)) / diff_roots;
	  // FIXME: forse conviene semplificare g_n

	  D_VAR(g_n);

	  add_initial_conditions(g_n, n, num_coefficients,
				 initial_conditions, solution);

	  D_MSGVAR("Before calling simplify: ", solution);

	  solution = simplify_on_output_ex(solution.expand(), n, false);
	  solution = solution.collect(lst(initial_conditions[0],
					  initial_conditions[1]));
	}
	else 
	  throw ("PURRS error: today we only allow inhomogeneous terms\n"
		 "in the form of polynomials or product of exponentials\n"
		 "and polynomials.\n"
		 "Please come back tomorrow.");
      else {
	// The characteristic equation
	// x^2 + a_1 * x + a_2 = 0 has a double root.
	assert(roots[0].multiplicity() == 2);      
	solution = order_2_sol_roots_no_distinct(n, base_of_exps,
						 exp_poly_coeff,
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
  Let \f$ e(n) \f$ be the expression in \p n contained in \p e,
  which is assumed to be already expanded.
  This function computes a decomposition
  \f$ e(n) = \sum_{i=0}^k \alpha_i^n \bigl(p_i(n) + q_i(n)\bigr) \f$, where
  - \f$ \alpha_i \f$ is a ground expression
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
  Given a vector \p symbolic_sum that contains all the symbolic sums of the
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
			    std::vector<GExpr>& base_of_exps,
			    std::vector<GExpr>& symbolic_sum_distinct,
			    std::vector<GExpr>& symbolic_sum_no_distinct,
			    GExpr& solution) {
  solution = 0;
  for (size_t i = roots.size(); i-- > 0; )
    for (size_t j = symbolic_sum_distinct.size(); j-- > 0; ) {
      GExpr base_exp = base_of_exps[j];
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

/*!
  Adds to sum related to the inhmogeneous terms the sum
  \f$ \sum_{i=0}^{order - 1} g_{n-i}
    \bigl( x_i - \sum_{j=1}^i a_j x_{i-j} \bigr) \f$
  corresponding to the initial conditions'part.
*/
// FIXME: il vettore 'coefficients' dovra' diventare di 'GExpr' quando
// sapremo risolvere anche le eq. di grado superiore al primo con i
// parametri.
static void
add_initial_conditions(GExpr& g_n, const GSymbol& n,
		       std::vector<GNumber>& coefficients,
		       const std::vector<GExpr> initial_conditions,
		       GExpr& solution) {
  // 'coefficients.size()' has 'order + 1' elements because in the first
  // position there is the value 0. 
  std::vector<GExpr> g_n_i(coefficients.size() - 1);
  for (size_t i = g_n_i.size(); i-- > 0; )
    g_n_i[i] = g_n.subs(n == n - i);
  for (size_t i = g_n_i.size(); i-- > 0; ) {
    GExpr tmp = initial_conditions[i];
    for (size_t j = i; j > 0; j--)
      tmp -= coefficients[j]*initial_conditions[i-j];
    solution += tmp * g_n_i[i];
  }
}

static GExpr
order_2_sol_roots_no_distinct(const GSymbol& n,
			      const std::vector<GExpr> base_of_exps,
			      const std::vector<GExpr> exp_poly_coeff,
			      const std::vector<GNumber>& coefficients,
			      const std::vector<Polynomial_Root>& roots) {
  GExpr solution_tot;

  unsigned num_exponentials = base_of_exps.size();

  GSymbol a("a"), b("b");
  GExpr root = roots[0].value();

  D_VAR(root);

  if (num_exponentials == 0) {
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
    // Solved the system with the inverse matrix method.
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

    for (size_t i = 0; i < num_exponentials; ++i) {
      GExpr solution = 0;
      GExpr base_of_exp = base_of_exps[i];
      GExpr coeff_of_exp = exp_poly_coeff[i];
      GSymbol k("k");
      GExpr coeff_of_exp_k = coeff_of_exp.subs(n == k);
      GExpr g_n_k = g_n.subs(n == n - k);
      // In this case g_n_k always contains root^(n-k).
      // We pass to the function 'sum_poly_exponentials' only
      // (1/root)^k that it is multiplied for 'base_of_exp'.
      // Hence we pass only the polynomial part of g_n with
      // the substitution n == n - k.
      GExpr poly_g_n = g_n.subs(wild(0) * power(wild(1), n) == wild(0));
      poly_g_n = poly_g_n.subs(n == n - k);
      solution = sum_poly_times_exponentials(coeff_of_exp_k * poly_g_n,
					     k, n, base_of_exp *
					     pow(root, -1));
      // 'sum_poly_times_exponentials' calculates the sum from 0 while
      // we want to start from 2.
      solution -= (coeff_of_exp_k * poly_g_n).subs(k == 0) +
	(coeff_of_exp_k * poly_g_n).subs(k == 1) * 
	(1/root) * base_of_exp;
      // We have passed to the function 'sum_poly_times_exponentials'
      // only (1/root)^k then now we must multiply for (root)^n.
      solution *= power(root, n);
      solution = solution.expand();
      solution_tot += solution;
    }
    solution_tot = simplify_on_output_ex(solution_tot.expand(), n, false);
    solution_tot = solution_tot.collect(lst(initial_conditions[0],
					    initial_conditions[1]));
  }
  return solution_tot;
}

/*
static int
verify_solution(const GExpr& solution, const int order, const GExpr& e,
		const std::vector<GNumber>& coefficients) {
  std::vector<GExpr> terms(order);
  for (size_t i = 0; i = order; ++i)
    term[i] = solution.subs(n == n - i - 1);
  int diff = solution;
  for (size_t i = terms.size(); i-- > 0; )
    diff -= coefficients[] 
  return r;
}
*/
