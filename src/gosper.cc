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
#include "numerator_denominator.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Recurrence.defs.hh"

#include <vector>
#include <algorithm>

// TEMPORARY
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

static const unsigned
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
  By definition, an expression \f$ t(n) \f$ is a
  <EM>hypergeometric term</EM> if \f$ t(n+1) / t(n) \f$ is a rational
  function of \f$ n \f$.
  This function returns <CODE>true</CODE> if \p t is a hypergeometric term
  and in this case \p r_n stores the ratio \f$ t(n+1) / t(n) \f$.
  Returns <CODE>false</CODE> otherwise.
*/
bool
gosper_step_one(const Expr& t_n, Expr& r_n, bool full) {
  // Not is the case of variable coefficient.
  if (full) 
    r_n = simplify_factorials_and_exponentials(t_n.substitute(Recurrence::n,
							      Recurrence::n+1))
      * pwr(simplify_factorials_and_exponentials(t_n), -1);
  // FIXME: we must understand the better simplification to use in this case.
  r_n = simplify_numer_denom(r_n);
  D_VAR(r_n);
  if (r_n.is_rational_function(Recurrence::n))
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
compute_resultant_and_its_roots(const Expr& f, const Expr& g,
				std::vector<Number>& integer_roots) {
  Symbol h("h");
  Expr temp_g = g.substitute(Recurrence::n, Recurrence::n + h);
  Expr R = resultant(f, temp_g, Recurrence::n);
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
  if (!find_divisors(constant_term, potential_roots))
    return false;
  // Find non-negative integral roots of the resultant.
  for(unsigned i = potential_roots.size(); i-- > 0; ) {
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
gosper_step_two(const Expr& r_n, Expr& a_n, Expr& b_n, Expr& c_n) {
  // Gosper's algorithm, step 2.1.
  Expr f;
  Expr g;
  numerator_denominator_purrs(r_n, f, g);
  // It is necessary to `expand()' in order to have the right answer from
  // `lcoeff()'.
  f = f.expand();
  g = g.expand();
  std::vector<Number> integer_roots;
  if (!compute_resultant_and_its_roots(f, g, integer_roots))
    return false;
  // Gosper's algorithm, step 2.2.
  // `a_n' and `b_n' are used below as starting values for a sequence
  // of polynomials.

  // normalize f and g, and store conversion factor in Z
  Expr lead_f = f.lcoeff(Recurrence::n);
  Expr lead_g = g.lcoeff(Recurrence::n);
  Expr Z = lead_f * pwr(lead_g, -1);
  a_n = f * pwr(lead_f, -1);
  b_n = g * pwr(lead_g, -1);
  // Computation of the output polynomials.
  c_n = 1;
  unsigned integer_roots_size = integer_roots.size();
  for (unsigned i = 0; i < integer_roots_size; ++i) {
    Expr temp_b_n
      = (b_n.substitute(Recurrence::n,
			Recurrence::n + integer_roots[i])).expand();
    Expr s = general_gcd(a_n, temp_b_n, Recurrence::n);
    a_n = quo(a_n, s, Recurrence::n);
    Expr temp_s = s.substitute(Recurrence::n,
			       Recurrence::n - integer_roots[i]);
    b_n = quo(b_n, temp_s, Recurrence::n);
    for (Number j = 1; j <= integer_roots[i]; ++j)
      c_n *= s.substitute(Recurrence::n, Recurrence::n - j);
  }
  a_n *= Z;
  // The polynomials `a_n' and `b_n' may have rational coefficients.
  // Multiply numerator and denominator of the fraction `a_n/b_n'
  // by suitable integers, so that the output values of `a_n' and `b_n'
  // have integer coefficients. 
  Number a_n_factor;
  a_n = convert_to_integer_polynomial(a_n, Recurrence::n, a_n_factor);
  a_n *= numerator(a_n_factor);
  b_n *= denominator(a_n_factor);
  Number b_n_factor;
  b_n = convert_to_integer_polynomial(b_n, Recurrence::n, b_n_factor);
  a_n *= denominator(b_n_factor);
  b_n *= numerator(b_n_factor);
  D_VAR(a_n);
  D_VAR(b_n);
  D_VAR(c_n);

  return true;
}

bool
find_polynomial_solution(const Number& deg_x, const Expr& a_n,
			 const Expr& b_n, const Expr& c_n, Expr& x_n) {
  unsigned deg_a = a_n.degree(Recurrence::n);
  unsigned deg_b = b_n.degree(Recurrence::n);
  unsigned deg_c = c_n.degree(Recurrence::n);
  // `number_of_coeffs' is the number of coefficients of 
  // the polynomial `p', that is, 1 + deg_x.
  unsigned number_of_coeffs = (1 + deg_x).to_unsigned();

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
    x_n += pwr(Recurrence::n, i) * unknowns.op(i);

  Expr x_n_shift = x_n.substitute(Recurrence::n, Recurrence::n+1);
  Expr b_shift = b_n.substitute(Recurrence::n, Recurrence::n-1);

  // Considers the recurrence relation to solve.
  Expr rr = a_n * x_n_shift - b_shift * x_n - c_n;
  rr = rr.expand();

  // Builds the lists to put in the matrix `rr_coefficients' and `rhs'.
  Expr_List equations;
  for (unsigned i = 0; i < number_of_unknowns; ++i) {
    Expr lhs = rr.coeff(Recurrence::n, i);
    equations.prepend(Expr(lhs, 0));
  }
    
  Expr solution = lsolve(equations, unknowns);  
  if (solution.nops() == 0)
    return false;

  // Builds the solution `x(n)'.
  for (unsigned i = 0; i < number_of_coeffs; ++i)
    x_n = x_n.substitute(unknowns.op(i), solution.op(i).op(1));

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
  Returns <CODE>false</CODE> otherwise, i. e., it does not find a non-zero
  polynomial solution \f$ x(n) \f$.
*/
bool
gosper_step_three(const Expr& a_n, const Expr& b_n, const Expr& c_n,
		  Expr& x_n) {
  // Gosper's algorithm, step 3.1.
  // Finds the degree of `x(n)'.
  unsigned deg_a = a_n.degree(Recurrence::n);
  unsigned deg_b = b_n.degree(Recurrence::n);
  unsigned deg_c = c_n.degree(Recurrence::n);
  Expr lead_a = a_n.lcoeff(Recurrence::n);
  Expr lead_b = b_n.lcoeff(Recurrence::n);
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
    Expr A = a_n.coeff(Recurrence::n, deg_a - 1);
    Expr B = b_n.substitute(Recurrence::n,
			    Recurrence::n - 1).coeff(Recurrence::n, deg_a - 1);
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
  return find_polynomial_solution(deg_x, a_n, b_n, c_n, x_n);
}

/*!
  Gosper's algorithm, step 4: see Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/
Expr
gosper_step_four(const Expr& t, const Expr& b_n, const Expr& c_n,
		 const Expr& x_n, const Number& lower, const Expr& upper,
		 Expr solution) {
  Expr shift_b = b_n.substitute(Recurrence::n, Recurrence::n - 1);
  Expr z_n = shift_b * x_n * t * pwr(c_n, -1);
  z_n = simplify_numer_denom(z_n);
  // The Gosper's algorithm computes summation with the lower limit `0'
  // and the upper limit `n - 1': in this case, once we have `z_n',
  // the sum that we are looking for is `z_n - z_0'.
  // In general the solution will be `z_n - z_{lower}'.
  solution = z_n - z_n.substitute(Recurrence::n, lower);
  // We must modify the sum if its upper limit is not `n - 1'.
  if (upper != Recurrence::n - 1)
    if (upper == Recurrence::n)
      solution = solution.substitute(Recurrence::n, Recurrence::n + 1);
    else {
      assert(upper.is_a_add() && upper.nops() == 2);
      const Expr& a = upper.op(0);
      const Expr& b = upper.op(1);
      assert(a == Recurrence::n || b == Recurrence::n);
      assert(a.is_a_number() || b.is_a_number());
      Number num;
      if (a == Recurrence::n)
	num = b.ex_to_number();
      if (b == Recurrence::n)
	num = a.ex_to_number();
      assert(num.is_integer());
      solution = solution.substitute(Recurrence::n, Recurrence::n + num + 1);
    }
  solution = simplify_numer_denom(solution);
  return solution;
}

} // anonymous namespace

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
  Let
  \f[
    S_n = \sum_{k=0}^{n-1} t_k
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
  -  it is possible to express \f$ S_n \f$ in closed form and the solution
     is stored in \p solution;
  -  it is not possible to express \f$ S_n \f$ in closed form and returns
     in \p solution the symbolic sum \f$ \sum_{k=lower_limit}^{upper} t_k \f$.
*/
bool
PURRS::full_gosper(const Expr& t_n, const Number& lower, const Expr& upper,
		   Expr& solution) {
  Expr r_n;
  if (!gosper_step_one(t_n, r_n, true))
    // `t' is not hypergeometric: no chance of using Gosper's algorithm.
    return false;
  Expr a_n;
  Expr b_n;
  Expr c_n;
  if (!gosper_step_two(r_n, a_n, b_n, c_n))
    // Problem in the computation of the resultant and its roots.
    return false;
  Expr x_n;
  if (gosper_step_three(a_n, b_n, c_n, x_n))
    solution = gosper_step_four(t_n, b_n, c_n, x_n, lower, upper, solution);
  else {
    // `t' is not Gosper-summable, i. e., there is not hypergeometric
    // solution.
    Symbol h;
    solution = PURRS::sum(h, lower, upper, t_n.substitute(Recurrence::n, h));
  }
  D_MSGVAR("The sum is: ", solution);

  return true;
}

/*!
  This function does not represent the full Gosper algorithm because
  it does not compute \f$ r(n) = t(n+1) / t(n) \f$ in the first step.
  Instead, \f$ r(n) \f$ is received as an argument.
  This is useful when solving first-order recurrence relations with
  variable coefficients, because it is easy to compute the expression
  \f$ r(n) \f$ from the coefficients of the recurrence itself.
  In other words, we avoid useless calls to simplification routines.
*/
bool
PURRS::partial_gosper(const Expr& t_n, Expr& r_n,
		      const Number& lower, const Expr& upper, Expr& solution) {
  if (!gosper_step_one(t_n, r_n, false))
    // `t' is not hypergeometric: no chance of using Gosper's algorithm.
    return false;
  Expr a_n;
  Expr b_n;
  Expr c_n;
  if (!gosper_step_two(r_n, a_n, b_n, c_n))
    // Problem in the computation of the resultant and its roots.
    return false;
  Expr x_n;
  if (gosper_step_three(a_n, b_n, c_n, x_n))
    solution = gosper_step_four(t_n, b_n, c_n, x_n, lower, upper, solution);
  else {
    // `t' is not Gosper-summable, i. e., there is not hypergeometric
    // solution.
    Symbol h;
    solution = PURRS::sum(h, lower, upper, t_n.substitute(Recurrence::n, h));
  }
  D_MSGVAR("The sum is: ", solution);

  return true;
}
