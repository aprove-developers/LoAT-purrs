/* To be written
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

#ifndef NOISY
#define NOISY 0
#endif

#include "gosper.hh"
#include "alg_eq_solver.hh"
#include "simplify.hh"
#include "util.hh"
#include <vector>
#include <algorithm>

// TEMPORARY
#include <iostream>

namespace Parma_Recurrence_Relation_Solver {

static const unsigned
FACTOR_THRESHOLD = 100;

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/

//! Gosper's algorithm, step 1: see Chapter 5 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
/*!
  By definition, an expression \f$ t(n) \f$ is an
  <EM>hypergeometric term</EM> if \f$ t(n+1) / t(n) \f$ is a rational
  function of \f$ n \f$.
  This function returns <CODE>true</CODE> if \p t is an hypergeometric term
  and in this case \p r_n stores the ratio \f$ t(n+1) / t(n) \f$.
  Returns <CODE>false</CODE> otherwise.
*/
static bool
gosper_step_one(const Expr& t_n, Expr& r_n, const Symbol& n, bool full) {
  // Not is the case of variable coefficient.
  if (full) {
    Expr t_plus_one = t_n.subs(n, n+1); 
    r_n = simplify_factorials_and_exponentials(t_plus_one, n)
      * pwr(simplify_factorials_and_exponentials(t_n, n), -1);
  }
  // The following use of `numerator_denominator()' simplify ulteriorly
  // `r_n' (we can not to call `simplify_numer_denom()' because it expandes
  // the expressions).
  Expr r_n_num;
  Expr r_n_den;
  r_n.numerator_denominator(r_n_num, r_n_den);
  r_n = r_n_num * pwr(r_n_den, -1);
  D_VAR(r_n);
  if (is_rational_function(r_n, n))
    return true;
  else {
    D_MSG("t(n) not hypergeometric term");
    return false;
  }
}

/*!
  The resultant's non-negative integral roots are needed in Step 2.1
  of Gosper's algorithm.
*/
static void
compute_resultant_and_its_roots(const Expr& f, const Expr& g, 
				const Symbol& n,
				std::vector<Number>& integer_roots) {
  Symbol h("h");
  Expr temp_g = g.subs(n, n + h);
  Expr R = resultant(f, temp_g, n);
  R = R.primpart(h);
  if (!R.is_integer_polynomial())
    R = convert_to_integer_polynomial(R, h);
  std::vector<Number> potential_roots;
  Number constant_term = abs(R.tcoeff(h).ex_to_number());
  // If `constant_term == 0', divide `R' by `h', and repeat.
  // The constant `0' is a root of the original resultant `R'
  // and is therefore pushed in the vector `potential_roots'.
  // We ensure that `0' is inserted only once.
  while (constant_term == 0) {
    if (potential_roots.size() == 0)
      potential_roots.push_back(0);
    R = quo(R, h, h);
    constant_term = abs(R.tcoeff(h).ex_to_number());
  }
  find_divisors(constant_term, potential_roots);
  // Find non-negative integral roots of the resultant.
  for(unsigned i = potential_roots.size(); i-- > 0; ) {
    Number temp = R.subs(h, potential_roots[i]).ex_to_number();
    if (temp == 0)
      integer_roots.push_back(potential_roots[i]);
  }

  // It is more efficient to have the roots sorted. 
  sort(integer_roots.begin(), integer_roots.end());
}

//! Gosper's algorithm, step 2: see \S5.3 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
/*!
  We assume that \p f and \p g are relatively prime polynomials, stored
  in the list \p n_d.
*/
static void
gosper_step_two(const Expr& r_n, const Symbol& n,
		Expr& a_n, Expr& b_n, Expr& c_n) {
  // Gosper's algorithm, step 2.1.
  Expr f;
  Expr g;
  r_n.numerator_denominator(f, g);
  // To do `expand()' is necessary in order to have right answers from
  // `lcoeff()'.
  f = f.expand();
  g = g.expand();
  std::vector<Number> integer_roots;
  compute_resultant_and_its_roots(f, g, n, integer_roots);
  // Gosper's algorithm, step 2.2.
  // `a_n' and `b_n' are used below as starting values for a sequence
  // of polynomials.

  // normalize f and g, and store conversion factor in Z
  Expr lead_f = f.lcoeff(n);
  Expr lead_g = g.lcoeff(n);
  Expr Z = lead_f * pwr(lead_g, -1);
  a_n = f * pwr(lead_f, -1);
  b_n = g * pwr(lead_g, -1);
  // Computation of the output polynomials.
  c_n = 1;
  unsigned integer_roots_size = integer_roots.size();
  for (unsigned i = 0; i < integer_roots_size; ++i) {
    Expr temp_b_n = (b_n.subs(n, n + integer_roots[i])).expand();
    Expr s = general_gcd(a_n, temp_b_n, n);
    a_n = quo(a_n, s, n);
    Expr temp_s = s.subs(n, n - integer_roots[i]);
    b_n = quo(b_n, temp_s, n);
    for (Number j = 1; j <= integer_roots[i]; ++j)
      c_n *= s.subs(n, n - j);
  }
  a_n *= Z;
  // The polynomials `a_n' and `b_n' may have rational coefficients.
  // Multiply numerator and denominator of the fraction `a_n/b_n'
  // by suitable integers, so that the output values of `a_n' and `b_n'
  // have integer coefficients. 
  Number a_n_factor;
  a_n = convert_to_integer_polynomial(a_n, n, a_n_factor);
  a_n *= a_n_factor.numerator();
  b_n *= a_n_factor.denominator();
  Number b_n_factor;
  b_n = convert_to_integer_polynomial(b_n, n, b_n_factor);
  a_n *= b_n_factor.denominator();
  b_n *= b_n_factor.numerator();
  D_VAR(a_n);
  D_VAR(b_n);
  D_VAR(c_n);
}

static bool
find_polynomial_solution(const Symbol& n, const Number& deg_x,
			 const Expr& a_n, const Expr& b_n, const Expr& c_n,
			 Expr& x_n) {
  unsigned deg_a = a_n.degree(n);
  unsigned deg_b = b_n.degree(n);
  unsigned deg_c = c_n.degree(n);
  // `number_of_coeffs' is the number of coefficients of 
  // the polynomial `p', that is, 1 + deg_x.
  unsigned number_of_coeffs = (1 + deg_x).to_int();

  // Compute the real number of unknowns of the polynomial equation
  // `a(n) * p(n+1) - b(n-1) * p(n) - c(n) = 0'.
  // In general, this is larger than `number_of_coeffs' because 
  // the polynomial equation above may have large degree.
  unsigned number_of_unknowns = number_of_coeffs;
  number_of_unknowns += deg_a > deg_b ? deg_a : deg_b;
  number_of_unknowns = number_of_unknowns > deg_c ? number_of_unknowns : deg_c;

  Expr_List unknowns;
  for (unsigned i = 0; i < number_of_unknowns; ++i)
    unknowns.append(Symbol());

  // Builds the generic polynomial `p' of degree `deg_x'.
  x_n = 0;
  for (unsigned i = 0; i < number_of_coeffs; ++i)
    x_n += pwr(n, i) * unknowns.op(i);

  Expr x_n_shift = x_n.subs(n, n+1);
  Expr b_shift = b_n.subs(n, n-1);

  // Considers the recurrence relation to solve.
  Expr rr = a_n * x_n_shift - b_shift * x_n - c_n;
  rr = rr.expand();

  // Builds the lists to put in the matrix `rr_coefficients' and `rhs'.
  Expr_List equations;
  for (unsigned i = 0; i < number_of_unknowns; ++i) {
    Expr lhs = rr.coeff(n, i);
    equations.prepend(Expr(lhs, 0));
  }
    
  Expr solution = lsolve(equations, unknowns);  
  if (solution.nops() == 0)
    return false;

  // Builds the solution `x(n)'.
  for (unsigned i = 0; i < number_of_coeffs; ++i)
    x_n = x_n.subs(unknowns.op(i), solution.op(i).op(1));

  D_VAR(x_n);
  return true;
}

//! Gosper's algorithm, step 3: see \S5.4 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
/*!
  Returns <CODE>true</CODE> if it finds, if exists, a non-zero polynomial
  solution \f$ x(n) \f$ of \f$ a(n) * x(n+1) - b(n-1) * x(n) = c(n) \f$,
  where \f$ a(n) \f$, \f$ b(n) \f$ and \f$ c(n) \f$ are polynomials such that
  \f$ gcd(a(n), b(n+h)) = 1 \f$, for all non-negative integers
  \f$ h \f$. The solution \f$ x(n) \f$ is stored in \p x_n. 
  Returns <CODE>false</CODE> otherwise, i. e., it not finds a non-zero
  polynomial solution \f$ x(n) \f$.
*/
static bool
gosper_step_three(const Expr& a_n, const Expr& b_n, const Expr& c_n,
		  const Symbol& n, Expr& x_n) {
  // Gosper's algorithm, step 3.1.
  // Finds the degree of `x(n)'.
  unsigned deg_a = a_n.degree(n);
  unsigned deg_b = b_n.degree(n);
  unsigned deg_c = c_n.degree(n);
  Expr lead_a = a_n.lcoeff(n);
  Expr lead_b = b_n.lcoeff(n);
  Number deg_x = -1;
  // On output, if a possible degree exists,
  // `deg_x' will contain a non-negative number.
  if (deg_a != deg_b || lead_a != lead_b) {
    if (deg_a >= deg_b && deg_c >= deg_a)
      deg_x = deg_c - deg_a;
    if (deg_b > deg_a && deg_c >= deg_b)
      deg_x = deg_c - deg_b;
  }
  else {
    // `deg_a = deg_b' and `lead_a = lead_b'.
    Expr shift_b = b_n.subs(n, n - 1);
    Expr A = a_n.coeff(n, deg_a - 1);
    Expr B = shift_b.coeff(n, deg_a - 1);
    Expr B_A_e = (B - A) * pwr(lead_a, -1);
    Number B_A = B_A_e.ex_to_number();
    Number possible_deg = Number(deg_c) - Number(deg_a) + 1;
    if (B_A.is_nonnegative_integer())
      if (B_A > possible_deg)
	possible_deg = B_A;
    deg_x = possible_deg >= 0 ? possible_deg : -1;
  }
  if (deg_x == -1)
    return false;
  D_MSGVAR("Degree of x(n): ", deg_x);

  // Gosper's algorithm, step 3.2.
  return find_polynomial_solution(n, deg_x, a_n, b_n, c_n, x_n);
}

/*!
  Gosper's algorithm, step 4: see Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/
static Expr
gosper_step_four(const Expr& t, const Expr& b_n, const Expr& c_n,
		 const Expr& x_n, const Symbol& n,
		 const Number& lower, const Expr& upper,
		 Expr solution) {
  Expr shift_b = b_n.subs(n, n-1);
  Expr z_n = shift_b * x_n * t * pwr(c_n, -1);
  z_n = simplify_numer_denom(z_n);
  // The Gosper's algorithm computes summation with the lower limit `0'
  // and the upper limit `n - 1': in this case, once we have `z_n',
  // the sum that we are looking for is `z_n - z_0'.
  // In general the solution will be `z_n - z_{lower}'.
  solution = z_n - z_n.subs(n, lower);
  // We must modify the sum if its upper limit is not `n - 1'.
  if (upper != n-1)
    if (upper == n)
      solution = solution.subs(n, n + 1);
    else {
      Expr n_plus_i = n + wild(0);
      assert(upper == n_plus_i);
      Expr_List substitution;
      upper.match(n_plus_i, substitution);
      Expr i = get_binding(substitution, 0);
      solution = solution.subs(n, n + i + 1);
    }
  solution = simplify_numer_denom(solution);
  return solution;
}

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
  Let
  \f[
    S_n = \sum_{k=0}^{n-1} t_k
  \f]
  with \f$ t_k \f$ an <EM>hypergeometric term</EM> that does not depend on
  \f$ n \f$, i. e., consecutive term ratio
  \f[
    r(k) = \frac{t_{k+1}}{t_k}
  \f]
  is a rational function of \f$ k \f$.
  This function returns <CODE>false</CODE> if \f$ t_k \f$ is not an
  hypergeometric term. 
  Returns <CODE>true</CODE> if \f$ t_k \f$ is an hypergeometric term.
  There are two case:
  -  it is possible to express \f$ S_n \f$ in closed form and the solution
     is stored in \p solution;
  -  it is not possible to express \f$ S_n \f$ in closed form and returns
     in \p solution the symbolic sum \f$ \sum_{k=lower_limit}^{upper} t_k \f$.
*/
bool
full_gosper(const Expr& t_n, const Symbol& n,
	    const Number& lower, const Expr& upper, Expr& solution) {
  Expr r_n;
  if (!gosper_step_one(t_n, r_n, n, true))
    // `t' is not hypergeometric: no chance of using Gosper's algorithm.
    return false;
  Expr a_n;
  Expr b_n;
  Expr c_n;
  gosper_step_two(r_n, n, a_n, b_n, c_n);
  Expr x_n;
  if (gosper_step_three(a_n, b_n, c_n, n, x_n))
    solution = gosper_step_four(t_n, b_n, c_n, x_n, n, lower, upper,
				solution);
  else {
    // `t' is not Gosper-summable, i. e., there is not hypergeometric
    // solution.
    Symbol h("h");
    Expr t_h = t_n.subs(n, h);
    solution = sum(Expr(h), Expr(lower), upper, t_h);
  }
  D_MSGVAR("The sum is: ", solution);

  return true;
}

/*!
  This function not represent the full Gosper's algorithm because
  not computes \f$ r(n) = t(n+1) / t(n) \f$ in the first step but
  \f$ r(n) \f$ is received like argument.
*/
bool
partial_gosper(const Expr& t_n, Expr& r_n, const Symbol& n,
	       const Number& lower, const Expr& upper, Expr& solution) {
  if (!gosper_step_one(t_n, r_n, n, false))
    // `t' is not hypergeometric: no chance of using Gosper's algorithm.
    return false;
  Expr a_n;
  Expr b_n;
  Expr c_n;
  gosper_step_two(r_n, n, a_n, b_n, c_n);
  Expr x_n;
  if (gosper_step_three(a_n, b_n, c_n, n, x_n))
    solution = gosper_step_four(t_n, b_n, c_n, x_n, n, lower, upper,
				solution);
  else {
    // `t' is not Gosper-summable, i. e., there is not hypergeometric
    // solution.
    Symbol h("h");
    Expr t_h = t_n.subs(n, h);
    solution = sum(Expr(h), Expr(lower), upper, t_h);
  }
  D_MSGVAR("The sum is: ", solution);

  return true;
}


} // namespace Parma_Recurrence_Relation_Solver
