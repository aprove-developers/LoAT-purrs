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
#include "simplify.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"

#include <vector>
#include <algorithm>

// TEMPORARY
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

static const unsigned int
FACTOR_THRESHOLD = 100;

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/

namespace {
using namespace PURRS;

//! Gosper's algorithm, step 1: see Chapter 5 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
/*!
  By definition, an expression \f$ t(m) \f$ is a
  <EM>hypergeometric term</EM> if \f$ t(m+1) / t(m) \f$ is a rational
  function of \f$ m \f$.
  This function returns <CODE>true</CODE> if \p t is a hypergeometric term
  and in this case \p r_m stores the ratio \f$ t(m+1) / t(m) \f$.
  Returns <CODE>false</CODE> otherwise.
*/
bool
gosper_step_one(const Symbol& m, const Expr& t_m, Expr& r_m) {
  D_VAR(t_m);
  r_m = simplify_binomials_factorials_exponentials(t_m.substitute(m, m+1))
    * pwr(simplify_binomials_factorials_exponentials(t_m), -1);
  // FIXME: we must understand the better simplification to use in this case.
  r_m = simplify_numer_denom(r_m);
  D_VAR(r_m);
  if (r_m.is_rational_function(m))
    return true;
  else {
    D_MSG("t(n) not hypergeometric term");
    return false;
  }
}

//! Compute the resultant's non-negative integral roots are needed
//! in Step 2.1 of Gosper's algorithm.
/*!
  Compute the resultant's non-negative integral roots are needed in Step 2.1
  of Gosper's algorithm.
  Returns <CODE>true</CODE> if the system is able to compute the resultant
  and its roots, <CODE>false</CODE> otherwise.
*/
bool
compute_resultant_and_its_roots(const Symbol& m, const Expr& f, const Expr& g,
				std::vector<Number>& integer_roots) {
  Symbol h("h");
  Expr temp_g = g.substitute(m, m + h);
  Expr R = resultant(f, temp_g, m);
  R = simplify_all(R);
  assert(R.is_rational_function(h));
  // We must try the roots of the rational function `R': if it is
  // a rational function but not a polynomial, we consider only the zeros
  // of its numerator.
  if (!R.is_polynomial(h))
    R = R.numerator();
  R = R.expand();
  R = R.primpart(h);
  if (!R.is_integer_polynomial(h))
    R = convert_to_integer_polynomial(R, h);
  std::vector<Number> potential_roots;
  if (!R.tcoeff(h).is_a_number())
    return false;
  Number constant_term = abs(R.tcoeff(h).ex_to_number());
  // If `constant_term == 0', divide `R' by `h', and repeat.
  // The constant `0' is a root of the original resultant `R'
  // and is therefore pushed in the vector `potential_roots'.
  // We ensure that `0' is inserted only once.
  while (constant_term == 0) {
    if (potential_roots.size() == 0)
      potential_roots.push_back(0);
    R = quo(R, h, h);
    if (!R.tcoeff(h).is_a_number())
      return false;
    constant_term = abs(R.tcoeff(h).ex_to_number());
  }
  if (!find_divisors(constant_term, potential_roots))
    return false;
  // Find non-negative integral roots of the resultant.
  for(unsigned int i = potential_roots.size(); i-- > 0; ) {
    Number temp = R.substitute(h, potential_roots[i]).ex_to_number();
    if (temp == 0)
      integer_roots.push_back(potential_roots[i]);
    }
  // It is more efficient to have the roots sorted. 
  sort(integer_roots.begin(), integer_roots.end());
  return true;
}

//! Gosper's algorithm, step 2: see \S5.3 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
/*!
  We assume that \p f and \p g are relatively prime polynomials, stored
  in the list \p n_d.
*/
bool
gosper_step_two(const Symbol& m, const Expr& r_m,
		Expr& a_m, Expr& b_m, Expr& c_m) {
  // Gosper's algorithm, step 2.1.
  Expr f;
  Expr g;
  r_m.numerator_denominator(f, g);
  // It is necessary to `expand()' in order to have the right answer from
  // `lcoeff()'.
  f = f.expand();
  g = g.expand();
  std::vector<Number> integer_roots;
  if (!compute_resultant_and_its_roots(m, f, g, integer_roots))
    return false;
  // Gosper's algorithm, step 2.2.
  // `a_m' and `b_m' are used below as starting values for a sequence
  // of polynomials.

  // normalize f and g, and store conversion factor in Z
  Expr lead_f = f.lcoeff(m);
  Expr lead_g = g.lcoeff(m);
  Expr Z = lead_f * pwr(lead_g, -1);
  a_m = f * pwr(lead_f, -1);
  b_m = g * pwr(lead_g, -1);
  // Computation of the output polynomials.
  c_m = 1;
  unsigned int integer_roots_size = integer_roots.size();
  for (unsigned int i = 0; i < integer_roots_size; ++i) {
    Expr temp_b_m = (b_m.substitute(m, m + integer_roots[i])).expand();
    Expr s = general_gcd(a_m, temp_b_m, m);
    a_m = quo(a_m, s, m);
    Expr temp_s = s.substitute(m, m - integer_roots[i]);
    b_m = quo(b_m, temp_s, m);
    for (Number j = 1; j <= integer_roots[i]; ++j)
      c_m *= s.substitute(m, m - j);
  }
  a_m *= Z;
  // The polynomials `a_m' and `b_m' may have rational coefficients.
  // Multiply numerator and denominator of the fraction `a_m/b_m'
  // by suitable integers, so that the output values of `a_m' and `b_m'
  // have integer coefficients. 
  if (!a_m.is_integer_polynomial(m)) {
    Expr a_m_factor;
    a_m = convert_to_integer_polynomial(a_m, m, a_m_factor);
    a_m *= a_m_factor.numerator();
    b_m *= a_m_factor.denominator();
  }
  if (!b_m.is_integer_polynomial(m)) {
    Expr b_m_factor;
    b_m = convert_to_integer_polynomial(b_m, m, b_m_factor);
    a_m *= b_m_factor.denominator();
    b_m *= b_m_factor.numerator();
  }
  D_VAR(a_m);
  D_VAR(b_m);
  D_VAR(c_m);

  return true;
}

bool
find_polynomial_solution(const Symbol& m, const Number& deg_x, const Expr& a_m,
			 const Expr& b_m, const Expr& c_m, Expr& x_m) {
  unsigned int deg_a = a_m.degree(m);
  unsigned int deg_b = b_m.degree(m);
  unsigned int deg_c = c_m.degree(m);
  // `number_of_coeffs' is the number of coefficients of 
  // the polynomial `p', that is, 1 + deg_x.
  unsigned int number_of_coeffs = (1 + deg_x).to_unsigned_int();

  // Compute the real number of unknowns of the polynomial equation
  // `a(m) * p(m+1) - b(m-1) * p(m) - c(m) = 0'.
  // In general, this is larger than `number_of_coeffs' because 
  // the polynomial equation above may have large degree.
  unsigned int number_of_unknowns = number_of_coeffs;
  number_of_unknowns += deg_a > deg_b ? deg_a : deg_b;
  number_of_unknowns = number_of_unknowns > deg_c ? number_of_unknowns : deg_c;

  Expr_List unknowns;
  for (unsigned int i = 0; i < number_of_unknowns; ++i)
    unknowns.append(Symbol());

  // Builds the generic polynomial `p' of degree `deg_x'.
  x_m = 0;
  for (unsigned int i = 0; i < number_of_coeffs; ++i)
    x_m += pwr(m, i) * unknowns.op(i);

  Expr x_m_shift = x_m.substitute(m, m+1);
  Expr b_shift = b_m.substitute(m, m-1);

  // Considers the recurrence relation to solve.
  Expr rr = a_m * x_m_shift - b_shift * x_m - c_m;
  rr = rr.expand();

  // Builds the lists to put in the matrix `rr_coefficients' and `rhs'.
  Expr_List equations;
  for (unsigned int i = 0; i < number_of_unknowns; ++i) {
    Expr lhs = rr.coeff(m, i);
    equations.prepend(Expr(lhs, 0));
  }
    
  Expr solution = lsolve(equations, unknowns);  
  if (solution.nops() == 0)
    return false;

  // Builds the solution `x(n)'.
  for (unsigned int i = 0; i < number_of_coeffs; ++i)
    x_m = x_m.substitute(unknowns.op(i), solution.op(i).op(1));

  D_VAR(x_m);
  return true;
}

//! Gosper's algorithm, step 3: see \S5.4 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
/*!
  Returns <CODE>true</CODE> if it finds, if exists, a non-zero polynomial
  solution \f$ x(m) \f$ of \f$ a(m) * x(m+1) - b(m-1) * x(m) = c(m) \f$,
  where \f$ a(m) \f$, \f$ b(m) \f$ and \f$ c(m) \f$ are polynomials such that
  \f$ gcd(a(m), b(m+h)) = 1 \f$, for all non-negative integers
  \f$ h \f$. The solution \f$ x(m) \f$ is stored in \p x_m. 
  Returns <CODE>false</CODE> otherwise, i. e., it does not find a non-zero
  polynomial solution \f$ x(m) \f$.
*/
bool
gosper_step_three(const Symbol& m, const Expr& a_m, const Expr& b_m,
		  const Expr& c_m, Expr& x_m) {
  // Gosper's algorithm, step 3.1.
  // Finds the degree of `x(n)'.
  unsigned int deg_a = a_m.degree(m);
  unsigned int deg_b = b_m.degree(m);
  unsigned int deg_c = c_m.degree(m);
  Expr lead_a = a_m.lcoeff(m);
  Expr lead_b = b_m.lcoeff(m);
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
    Expr A = a_m.coeff(m, deg_a - 1);
    Expr B = b_m.substitute(m, m - 1).expand().coeff(m, deg_a - 1);
    Number B_A = ((B - A) * pwr(lead_a, -1)).ex_to_number();
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
  return find_polynomial_solution(m, deg_x, a_m, b_m, c_m, x_m);
}

/*!
  Gosper's algorithm, step 4: see Chapter 5 of \f$ A = B \f$, by 
  M.~Petkovsek, H.~Wilf and D.~Zeilberger.
*/
Expr
gosper_step_four(const Symbol& m, const Expr& t,
		 const Expr& b_m, const Expr& c_m, const Expr& x_m,
		 const Number& lower, const Expr& upper, Expr solution) {
  Expr shift_b = b_m.substitute(m, m - 1);
  Expr z_m = shift_b * x_m * t * pwr(c_m, -1);
  z_m = simplify_numer_denom(z_m);
  // The Gosper's algorithm computes summation with the lower limit `0'
  // and the upper limit `m - 1': in this case, once we have `z_m',
  // the sum that we are looking for is `z_m - z_0'.
  // In general the solution will be `z_m - z_{lower}'.
  solution = z_m - z_m.substitute(m, lower);
  // We must modify the sum if its upper limit is not `m - 1'.
  if (upper != m - 1)
    if (upper == m)
      solution = solution.substitute(m, m + 1);
    else {
      assert(upper.is_a_add() && upper.nops() == 2);
      const Expr& a = upper.op(0);
      const Expr& b = upper.op(1);
      assert(a == m || b == m);
      assert(a.is_a_number() || b.is_a_number());
      Number num;
      if (a == m)
	num = b.ex_to_number();
      if (b == m)
	num = a.ex_to_number();
      assert(num.is_integer());
      solution = solution.substitute(m, m + num + 1);
    }
  solution = simplify_numer_denom(solution);
  return solution;
}

} // anonymous namespace

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkovsek, H.~Wilf and D.~Zeilberger.
  Let
  \f[
    S_m = \sum_{k=0}^{m-1} t_k
  \f]
  with \f$ t_k \f$ a <EM>hypergeometric term</EM> that does not depend on
  \f$ n \f$, i. e., consecutive term ratio
  \f[
    r(k) = \frac{t_{k+1}}{t_k}
  \f]
  is a rational function of \f$ k \f$.
  This function returns <CODE>false</CODE> if \f$ t_k \f$ is not a
  hypergeometric term.
  It returns <CODE>true</CODE> if \f$ t_k \f$ is a hypergeometric term.
  There are two cases:
  -  it is possible to express \f$ S_m \f$ in closed form and the solution
     is stored in \p solution;
  -  it is not possible to express \f$ S_m \f$ in closed form and returns
     in \p solution the symbolic sum \f$ \sum_{k=lower}^{upper} t_k \f$.
*/
bool
PURRS::gosper_algorithm(const Symbol& m, const Expr& t_m,
			const Number& lower, const Expr& upper,
			Expr& solution) {
  Expr r_m;
  if (!gosper_step_one(m, t_m, r_m))
    // `t_m' is not hypergeometric: no chance of using Gosper's algorithm.
    return false;
  Expr a_m;
  Expr b_m;
  Expr c_m;
  if (!gosper_step_two(m, r_m, a_m, b_m, c_m))
    // Problem in the computation of the resultant and its roots.
    return false;
  Expr x_m;
  if (gosper_step_three(m, a_m.expand(), b_m.expand(), c_m.expand(), x_m))
    solution = gosper_step_four(m, t_m, b_m, c_m, x_m, lower, upper, solution);
  else {
    // `t_m' is not Gosper-summable, i. e., there is not hypergeometric
    // solution.
    Symbol h;
    solution = PURRS::sum(h, lower, upper, t_m.substitute(m, h));
  }
  D_MSGVAR("The sum is: ", solution);

  return true;
}
