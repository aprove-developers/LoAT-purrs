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
#include <fstream>

using namespace GiNaC;

/*!
  Returns <CODE>true</CODE> if and only if \p e is of the form
  \f$ n + d \f$ with \f$ d \f$ a <CODE>GiNaC::numeric</CODE>;
  in this case the opposite of \f$ d \f$ is assigned to \p decrement.
*/
static bool
get_constant_decrement(const GExpr& e, const GSymbol& n, GNumber& decrement) {
  static GExpr n_plus_d = n + wild(0);
  GList substitution;
  if (match(e, n_plus_d, substitution)) {
    GExpr d = get_binding(substitution, 0);
    if (is_a<numeric>(d)) {
      decrement = -ex_to<numeric>(d);
      return true;
    }
  }
  return false;
}

/*!
  Given a vector of <CODE>GNumber</CODE> with a number of elements like
  the order \p order of the recurrence relation, builds the characteristic
  equation in the variable \p x.
  Therefore, if we have the recurrence relation
  \f$ x_n = a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$
  of order \f$ k \f$ with the coefficients \f$ a_j \f$ in the vector
  \p coefficients, this function returns the characteristic equation
  \f$ x^k - ( a_1 * x^{k-1} + \dotsb + a^{k-1} * x + a_k ) \f$.

  Since <CODE>find_roots()</CODE> solve equations only with integer
  coefficients the \p coefficients' elements must be rationals and,
  if they are not integer, builds an other vector, \p int_coefficients,
  with the element of \p coefficients multiplied for the lcm among their
  denominators.
*/
static GExpr
build_characteristic_equation(const int order, const GSymbol& x,
			      const std::vector<GNumber>& coefficients) {
  unsigned coefficients_size = coefficients.size();  
  for (unsigned i = 1; i < coefficients_size; ++i)
    if (!coefficients[i].is_rational())
      throw("PURRS error: today the algebraic equation solver works\n"
            "only with integer coefficients.\n"
            "Please come back tomorrow.");
  std::vector<GNumber> denominators;
  // Find the least common multiple among the denominators of the
  // rational elements of 'coefficients'.
  for (unsigned i = 1; i < coefficients_size; ++i)
    if (!coefficients[i].is_integer())
      denominators.push_back(coefficients[i].denom());
  GExpr p = 0;
  // Build the vector 'int_coefficients' with the elements of
  // 'coefficients' multiplied for the least common multiple
  // 'least_com_mul'.
  if (denominators.size() != 0) {
    GNumber least_com_mul = lcm(denominators);
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

static GExpr
return_sum(const bool distinct, const GSymbol& n, const int order,
	   const GExpr& coeff, const GSymbol& alpha, const GSymbol& lambda) {
  GSymbol k("k");
  GSymbol x("x");
  GExpr q_k = coeff.subs(n == k);
  GExpr symbolic_sum;  
  if (distinct)
    symbolic_sum = sum_poly_times_exponentials(q_k, k, n, x);
  else
    symbolic_sum = sum_poly_times_exponentials(q_k, k, n, 1);
  // 'sum_poly_times_exponentials' calculates the sum from 0 while
  // we want to start from 'order'.
  symbolic_sum -= q_k.subs(k == 0);
  for (int j = 1; j < order; ++j)
    symbolic_sum -= q_k.subs(k == j) * alpha * pow(lambda, -1);
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
  \p symbolic_sum_distinct and \p symbolic_sum_no_distinct
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
compute_symbolic_sum(const int order, const GSymbol& n,
		     const GSymbol& alpha, const GSymbol& lambda,
		     const std::vector<Polynomial_Root>& roots,
		     std::vector<GExpr>& base_of_exps,
		     std::vector<GExpr>& exp_poly_coeff,
		     std::vector<GExpr>& symbolic_sum_distinct,
		     std::vector<GExpr>& symbolic_sum_no_distinct) {
  for (unsigned i = base_of_exps.size(); i-- > 0; ) {
    bool distinct = true;
    for (unsigned j = roots.size(); j-- > 0; ) {
      if (roots[j].value().is_equal(base_of_exps[i]))
	distinct = false;
      
      // The root is different from the exponential's base.
      if (distinct) {
	symbolic_sum_distinct.push_back(return_sum(true, n, order,
						   exp_poly_coeff[i],
						   alpha, lambda));
	symbolic_sum_no_distinct.push_back(0);
      }
      // The root is equal to the exponential's base.
      else {
	symbolic_sum_no_distinct.push_back(return_sum(false, n, order,
						      exp_poly_coeff[i],
						      alpha, lambda));
	symbolic_sum_distinct.push_back(0);
      }
    }
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
    exponential.
  Returns a <CODE>GExpr</CODE> \p solution with the sum of all sums of the
  vectors.
 */
static void
subs_to_sum_roots_and_bases(const GSymbol& alpha, const GSymbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    std::vector<GExpr>& base_of_exps,
			    std::vector<GExpr>& symbolic_sum_distinct,
			    std::vector<GExpr>& symbolic_sum_no_distinct,
			    GExpr& solution) {
  solution = 0;
  unsigned ind = 0;
  for (unsigned i = base_of_exps.size(); i-- > 0; )
    for (unsigned j = roots.size(); j-- > 0; ) {
      GExpr base_exp = base_of_exps[i];
      GExpr tmp;
      if (!base_exp.is_equal(roots[j].value()))
	tmp = 
	  symbolic_sum_distinct[ind].subs(lst(alpha == base_exp,
					      lambda == roots[j].value()));
      else
	tmp = 
	  symbolic_sum_no_distinct[ind].subs(lst(alpha == base_exp,
						 lambda == roots[j].value()));
      solution += tmp * pow(-1, j);
      ++ind;
    }
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
order_2_sol_roots_no_distinct(const GSymbol& n, GExpr& g_n,
			      const std::vector<GExpr>& base_of_exps,
			      const std::vector<GExpr>& exp_poly_coeff,
			      const std::vector<GExpr>& initial_conditions,
			      const std::vector<GNumber>& coefficients,
			      const std::vector<Polynomial_Root>& roots);

static GExpr
solve_linear_constant_coeff(const GSymbol& n, GExpr& g_n,
			    const int order, const bool all_distinct,
			    const std::vector<GExpr>& base_of_exps,
			    const std::vector<GExpr>& exp_poly_coeff,
			    const std::vector<GNumber>& coefficients,
			    const std::vector<Polynomial_Root>& roots);

static bool
verify_solution(const GExpr& solution, const int& order, const GExpr& rhs,
		const GSymbol& n);

/*!
...
*/
bool
solve(const GExpr& rhs, const GSymbol& n, GExpr& solution) {
  static GExpr x_i = x(wild(0));
  static GExpr x_i_plus_r = x_i + wild(1);
  static GExpr a_times_x_i = wild(1)*x_i;
  static GExpr a_times_x_i_plus_r = a_times_x_i + wild(2);

  int order = -1;
  std::vector<GExpr> coefficients;
  // 'expand()' is necessary to simplify expressions, to find coefficients
  // and to decompose the inhomogeneous term of the recurrence relation.  
  GExpr e = rhs.expand();

  // Special case: 'e' is only a function in 'n' or a constant.
  // If there is not in 'e' 'x(argument)' or there is but with
  // 'argument' that not contains the variable 'n', then the solution
  // is obviously 'e'.
  GList occurrences;
  bool finished = true;
  if (e.find(x_i, occurrences)) {
    unsigned occurrences_nops = occurrences.nops();
    if (occurrences_nops != 0)
      for (unsigned i = 0; i < occurrences_nops; ++i) {
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
  // The vector 'base_of_exps' contains the exponential's bases of all
  // exponentials in 'e'. In the i-th position of the vectors 'exp_poly_coeff'
  // and 'exp_no_poly_coeff' there are respectively the polynomial part and
  // non polynomial part of the coefficient of the exponential with the base
  // in i-th position of 'base_of_exp' so that
  // exp_poly_coeff[i] + exp_no_poly_coeff[i] represents the coefficient of
  // base_of_exps[i]^n. 
  std::vector<GExpr> base_of_exps;
  std::vector<GExpr> exp_poly_coeff;
  std::vector<GExpr> exp_no_poly_coeff;
  exp_poly_decomposition(e, n, base_of_exps, exp_poly_coeff,
			 exp_no_poly_coeff);

  D_VEC(base_of_exps, 0, base_of_exps.size()-1);
  D_VEC(exp_poly_coeff, 0, exp_poly_coeff.size()-1);
  D_VEC(exp_no_poly_coeff, 0, exp_no_poly_coeff.size()-1);

  // Creates the vector of initial conditions.
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
  // Temporaneo fino a che 'find_roots()' accettera' i parametri anche
  // per equazioni di grado superiore al primo. 
  std::vector<GNumber> num_coefficients(order+1);
  if (order == 1) {
    characteristic_eq = y - coefficients[1];
    roots.push_back(coefficients[1]);
  } 
  else {
    // Check that the vector 'coefficients' contains only numeric
    // elements and in this case uses a vector of GNumber.
    for (int i = 1; i <= order; ++i)
      if (is_a<numeric>(coefficients[i]))
	num_coefficients[i] = ex_to<numeric>(coefficients[i]);
      else
	throw("PURRS error: today the recurrence relation\n"
	      "does not support irrationals coefficients.\n"
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
  std::vector<GExpr> symbolic_sum_distinct;
  std::vector<GExpr> symbolic_sum_no_distinct;

  if (all_distinct) {
    // Fills the two vector 'symbolic_sum_distinct' and
    // 'symbolic_sum_no_distinct' with two different symbolic sum
    // in according to the roots are equal or not to the exponential's bases.
    compute_symbolic_sum(order, n, alpha, lambda, roots,
			 base_of_exps, exp_poly_coeff,
			 symbolic_sum_distinct, symbolic_sum_no_distinct);
  }
  GExpr g_n;
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
	// g_n = pow(roots[0].value(), n);
	solution += initial_conditions[0] * pow(roots[0].value() ,n);
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
	  for (unsigned j = symbolic_sum_distinct.size(); j-- > 0; ) {
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
	  g_n = (pow(root_1, n+1) - pow(root_2, n+1)) / diff_roots;
	  // FIXME: forse conviene semplificare g_n
	  D_VAR(g_n);
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
	solution = order_2_sol_roots_no_distinct(n, g_n, base_of_exps,
						 exp_poly_coeff,
						 initial_conditions,
						 num_coefficients, roots);
      }
      break;
    }
  default:
      // Calculates the solution of linear recurrence relations with
      // constant coefficients (of order grater than two) when
      // the inhomogeneous term is a polynomial or the product of a
      // polynomial and an exponential.
      if (check_poly_times_exponential(exp_no_poly_coeff)) {
	solution = solve_linear_constant_coeff(n, g_n, order, all_distinct,
					       base_of_exps, exp_poly_coeff,
					       num_coefficients, roots);
	D_VAR(solution);
      }
      else
	throw ("PURRS error: today we only allow inhomogeneous terms\n"
	       "in the form of polynomials or product of exponentials\n"
	       "and polynomials.\n"
	       "Please come back tomorrow.");
      break;
  }
  
  if (order > 2 || (order == 2 && all_distinct))
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
  
  return true;
  }

/*
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
    assign_polynomial_part(p, n, p_poly, p_no_poly);
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
	assign_polynomial_part(coeff, n, coeff_poly, coeff_no_poly);
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
      assign_polynomial_part(q, n, q_poly, q_no_poly);
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
*/

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
  assign_polynomial_part(e, n, polynomial, possibly_not_polynomial);
  p[position] += polynomial;
  q[position] += possibly_not_polynomial;
}

/*!
  Definition of a <EM>valid_base</EM> for an exponential in inductive way:
  - every <CODE>GiNaC::numeric</CODE> is a <EM>valid_base</EM>;
  - every <CODE>GiNaC::constant</CODE> is a <EM>valid_base</EM>;
  - every <CODE>GiNaC::symbol</CODE> different from \p n is a
    <EM>valid_base</EM>;
  - given \f$ e \f$ a <CODE>GiNaC::power</CODE>,
    if \f$ e.op(0) \f$ and \f$ e.op(1) \f$ are <EM>valid_base</EM>
    then \f$ e \f$ is a <EM>valid_base</EM>;
  - given \f$ f \f$ a <CODE>GiNaC::function</CODE>,
    if \f$ f.op(0) \f$ is <EM>valid_base</EM>
    then \f$ f \f$ is a <EM>valid_base</EM>;
  - given the binary operations sum (\f$ + \f$) and multiplication (\f$ * \f$),
    if \f$ a \f$ and \f$ b \f$ are <EM>valid_base</EM> then
    \f$ a + b \f$ and \f$ a * b \f$ are <EM>valid_base</EM>.
*/
static bool
valid_base(const GExpr& e, const GSymbol& n) {
  bool ok;
  if (is_a<numeric>(e))
    ok = true;
  else if (is_a<constant>(e))
    ok = true;
  else if (is_a<symbol>(e) && !e.is_equal(n))
    ok = true;
  else if (is_a<power>(e))
    ok = valid_base(e.op(0), n) && valid_base(e.op(1), n);
  else if (is_a<function>(e))
    ok = valid_base(e.op(0), n);
  else if (is_a<add>(e) || is_a<mul>(e)) {
    ok = true;
    for (unsigned i = e.nops(); i-- > 0; )
      ok = ok && valid_base(e.op(i), n);
  }
  else
    ok = false;
  
  return ok;
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
      if (valid_base(base, n)) {
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
	if (valid_base(base, n)) {
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
// FIXME: il vettore 'coefficients' dovra' diventare di 'GExpr' quando
// sapremo risolvere anche le eq. di grado superiore al primo con i
// parametri.
static void
add_initial_conditions(const GExpr& g_n, const GSymbol& n,
		       const std::vector<GNumber>& coefficients,
		       const std::vector<GExpr>& initial_conditions,
		       GExpr& solution) {
  // 'coefficients.size()' has 'order + 1' elements because in the first
  // position there is the value 0. 
  for (unsigned i = coefficients.size() - 1; i-- > 0; ) {
    GExpr g_n_i = g_n.subs(n == n - i);
    GExpr tmp = initial_conditions[i];
    for (unsigned j = i; j > 0; j--)
      tmp -= coefficients[j] * initial_conditions[i-j];
    solution += tmp * g_n_i;
  }
}

static GExpr
order_2_sol_roots_no_distinct(const GSymbol& n, GExpr& g_n,
			      const std::vector<GExpr>& base_of_exps,
			      const std::vector<GExpr>& exp_poly_coeff,
			      const std::vector<GExpr>& initial_conditions,
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
    g_n = (a + b * n) * pow(root, n);
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

    for (unsigned i = 0; i < num_exponentials; ++i) {
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
  }
  return solution_tot;
}

static void
prepare_system(const int order, const bool all_distinct,
	       const std::vector<GNumber>& coefficients,
	       const std::vector<Polynomial_Root>& roots,
	       GList& g_i, GList& coeff_equations) {
  // Prepares a list with the elments for the 'rhs' of the system
  // to will find the numbers 'alpha_i' (for i = 1,...,k where k is the
  // order of the recurrence).
  // The first element is g_0, then g_1, ..., g_k.
  // Attenzione: nel caso bool distinct == false cambia! NON E' VERO!!!!
  std::vector<GExpr> tmp(order);
  tmp[0] = 1;
  for (int i = 1; i < order; ++i)
    for (int j = 0; j < i; ++j)
      tmp[i] += coefficients[j+1] * tmp[i-j-1];
  for (int i = 0; i < order; ++i)
    g_i.append(tmp[i]);
  
  // Prepares a list with the coefficients to insert in matrix 'system'
  // to will find the numbers 'alpha_i' (for i = 1,...,k).
  // This calculus is based on
  // 'g_n = \alpha_1 \lambda_1^n + \cdots + \alpha_k \lambda_k^n' if the roots
  // are all distinct; on
  // 'g_n = \sum_{j=1}^r (\alpha_{j,0} + \alpha_{j,1}n
  //        + \cdots + \alpha_{j,\mu_j-1}n^{\mu_j-1}) \lambda_j^n'
  // if there are multiple roots.
  if (all_distinct)
    for (int i = 0; i < order; ++i)
      for (int j = 0; j < order; ++j)
	coeff_equations.append(pow(roots[j].value(), i)); 
  else
    for (int h = 0; h < order; ++h)
      for (unsigned i = roots.size(); i-- > 0; ) {
	for (GNumber j = roots[i].multiplicity(); j-- > 1; )
	  coeff_equations.append(pow(h, j) * pow(roots[i].value(), h));
	coeff_equations.append(pow(roots[i].value(), h));
      }
#if NOISY
  std::cout << "g_i: " << g_i << std::endl;
  std::cout << "equations: " << coeff_equations << std::endl;
#endif
}

static GExpr
compute_non_homogeneous_part_solution(const GSymbol& n, const GExpr& g_n,
				      const int order,
				      const std::vector<GExpr>& base_of_exps,
				      const std::vector<GExpr>&
				      exp_poly_coeff) {
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
      // 'sum_poly_times_exponentials' calculates the sum from 0 while
      // we want to start from 'order'.
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
  Consider the linear recurrence relation of order \f$ k \f$ with constant
  coefficients
  \f$ x_n = a_1 x_{n-1} + a_2 x_{n-2} + \cdots + a_k x_{n-k} + p(n) \f$,
  where \f$ p(n) \f$ is a polynomial or a product of polynomials times
  exponentials.
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
  (\f$ g_n \f$ in the fisrt case is contained those of the second case as
  special case).
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
  part \f$ p(n) \f$ and to the initial conditions, respectively.
*/
static GExpr
solve_linear_constant_coeff(const GSymbol& n, GExpr& g_n,
			    const int order, const bool all_distinct,
			    const std::vector<GExpr>& base_of_exps,
			    const std::vector<GExpr>& exp_poly_coeff,
			    const std::vector<GNumber>& coefficients,
			    const std::vector<Polynomial_Root>& roots) {

  GList g_i;
  GList coeff_equations;
  // Prepare the elements to insert in the system.
  prepare_system(order, all_distinct, coefficients, roots,
		 g_i, coeff_equations);  
  GMatrix system(order, order, coeff_equations);
  GMatrix rhs(order, 1, g_i);
#if NOISY
  std::cout << "system: " << system << std::endl;
  std::cout << "rhs: " << rhs << std::endl;
#endif
  GMatrix vars(order, 1);
  for (int i = 0; i < order; ++i)
    vars(i, 0) = GSymbol();
  // Solve system in order to finds 'alpha_i' (i = 1,...,order).
  GMatrix sol = system.solve(vars, rhs);
#if NOISY
  std::cout << "alpha_i: " << sol << std::endl;
#endif
  // Finds 'g_n', always taking into account the root's multiplicity
  g_n = 0;
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
  D_VAR(g_n);
  // Computes '\sum_{i=k}^n g_{n-i} p(i)'.
  GExpr solution = compute_non_homogeneous_part_solution(n, g_n, order,
							 base_of_exps,
							 exp_poly_coeff);
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
  // that has all terms of 'solution' minus those containing an initial
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
