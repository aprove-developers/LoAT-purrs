/* ************************
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

#include "gosper.hh"

#include "alg_eq_solver.hh"
#include "simplify.hh"
#include "util.hh"
#include <vector>
#include <algorithm>

// TEMPORARY
#include <iostream>

using namespace GiNaC;

#define NOISY 0

static const unsigned
FACTOR_THRESHOLD = 100;

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/

/*!
  By definition, an expression \p p is a <EM>hypergeometric term</EM>
  if \f$ p(n+1) / p(n) \f$ is a rational function of \f$ n \f$.
  If this is the case, return also numerator and denominator.
*/
// static bool
// is_hypergeometric_term(const GExpr& p, const GSymbol& n,
// 		       GExpr& num_den) {

//   GExpr q = (p.subs(n == n+1)) / p;

//   q.expand(); // FIXME: simplify factorials here!
//   if (!q.info(info_flags::rational_function))
//     return false;
// //    here GiNaC guarantees that numerator and denominator are coprime
//   num_den = q.numer_denom();
//   return true;
// }

//! Gosper's algorithm, step 1: see Chapter 5 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
static bool
gosper_step_one(const GExpr&/*t*/, const GSymbol& /*n*/, GExpr& num_den_r_n) {
  // FIXME: general simplifications to be inserted here.
  num_den_r_n = numer_denom(num_den_r_n);
  return true;
}

/*!
  The resultant's non-negative integral roots are needed in Step 2.1
  of Gosper's algorithm.
*/
static void
compute_resultant_and_its_roots(const GExpr& f, const GExpr& g, 
				const GSymbol& n,
				std::vector<GNumber>& integer_roots) {
  GSymbol h("h");
  GExpr temp_g = g.subs(n == n + h);
  GExpr R = resultant(f, temp_g, n);
  R = R.primpart(h);
  if (!R.info(info_flags::integer_polynomial))
    R = convert_to_integer_polynomial(R, h);
  
  std::vector<GNumber> potential_roots;
  GNumber constant_term = abs(ex_to<GiNaC::numeric>(R.tcoeff(h)));
  // If `constant_term == 0', divide `R' by `h', and repeat.
  // The constant `0' is a root of the original resultant `R'
  // and is therefore pushed in the vector `potential_roots'.
  // We ensure that `0' is inserted only once.
  while (constant_term == 0) {
    if (potential_roots.size() == 0)
      potential_roots.push_back(0);
    R = quo(R, h, h);
    constant_term = abs(ex_to<GiNaC::numeric>(R.tcoeff(h)));
  }
  find_divisors(constant_term, potential_roots);
  // Find non-negative integral roots of the resultant.
  for(unsigned i = potential_roots.size(); i-- > 0; ) {
    GNumber temp = ex_to<GiNaC::numeric>(R.subs(h == potential_roots[i]));
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
gosper_step_two(const GExpr& r_n, const GSymbol& n,
		GExpr& a_n, GExpr& b_n, GExpr& c_n) {
  // Gosper's algorithm, step 2.1.
  GExpr f = expand(r_n.op(0)); // the numerator
  GExpr g = expand(r_n.op(1)); // the denominator

  std::vector<GNumber> integer_roots;
  compute_resultant_and_its_roots(f, g, n, integer_roots);

  // Gosper's algorithm, step 2.2.
  // `a_n' and `b_n' are used below as starting values for a sequence
  // of polynomials.

  // normalize f and g, and store conversion factor in Z
  GExpr lead_f = f.lcoeff(n);
  GExpr lead_g = g.lcoeff(n);
  GExpr Z = lead_f * pow(lead_g, -1);
  a_n = f * pow(lead_f, -1);
  b_n = g * pow(lead_g, -1);

  // Computation of the output polynomials.
  c_n = 1;
  unsigned integer_roots_size = integer_roots.size();
  for (unsigned i = 0; i < integer_roots_size; ++i) {
    GExpr temp_b_n = (b_n.subs(n == n + integer_roots[i])).expand();
    GExpr s = general_gcd(a_n, temp_b_n, n);
    a_n = quo(a_n, s, n);
    GExpr temp_s = s.subs(n == n - integer_roots[i]);
    b_n = quo(b_n, temp_s, n);
    for (GNumber j = 1; j <= integer_roots[i]; ++j)
      c_n *= s.subs(n == n - j);
  }
  a_n *= Z;
  // The polynomials `a_n' and `b_n' may have rational coefficients.
  // Multiply numerator and denominator of the fraction `a_n/b_n'
  // by suitable integers, so that the output values of `a_n' and `b_n'
  // have integer coefficients. 
  GNumber a_n_factor;
  a_n = convert_to_integer_polynomial(a_n, n, a_n_factor);
  GExpr a_n_d = numer_denom(a_n_factor);
  a_n *= a_n_d.op(0);
  b_n *= a_n_d.op(1);
  GNumber b_n_factor;
  b_n = convert_to_integer_polynomial(b_n, n, b_n_factor);
  GExpr b_n_d = numer_denom(b_n_factor);
  a_n *= b_n_d.op(1);
  b_n *= b_n_d.op(0);
#if NOISY
  std::cout << "a(n) = " << a_n << std::endl;
  std::cout << "b(n) = " << b_n << std::endl;
  std::cout << "c(n) = " << c_n << std::endl;
#endif
}

static bool
find_polynomial_solution(const GSymbol& n, const GNumber& deg_x,
			 const GExpr& a_n, const GExpr& b_n, const GExpr& c_n,
			 GExpr& x_n) {
  unsigned deg_a = a_n.degree(n);
  unsigned deg_b = b_n.degree(n);
  unsigned deg_c = c_n.degree(n);
  // `number_of_coeffs' is the number of coefficients of 
  // the polynomial `p', that is, 1 + deg_x.
  // FIXME: is it possible to convert a GiNaC integer to an int or an unsigned?
  // For the moment, we use a dirty trick.
  unsigned number_of_coeffs = 0;
  for (unsigned i = 0; i <= deg_x; ++i)
    ++number_of_coeffs;

  // Compute the real number of unknowns of the polynomial equation
  // `a(n) * p(n+1) - b(n-1) * p(n) - c(n) = 0'.
  // In general, this is larger than `number_of_coeffs' because 
  // the polynomial equation above may have large degree.
  unsigned number_of_unknowns = number_of_coeffs;
  number_of_unknowns += deg_a > deg_b ? deg_a : deg_b;
  number_of_unknowns = number_of_unknowns > deg_c ? number_of_unknowns : deg_c;

  GList unknowns;
  for (unsigned i = 0; i < number_of_unknowns; ++i)
    unknowns.append(GSymbol());

  // Builds the generic polynomial `p' of degree `deg_x'.
  x_n = 0;
  for (unsigned i = 0; i < number_of_coeffs; ++i)
    x_n += pow(n, i) * unknowns.op(i);

  GExpr x_n_shift = x_n.subs(n == n+1);
  GExpr b_shift = b_n.subs(n == n-1);

  // Considers the recurrence relation to solve.
  GExpr rr = a_n * x_n_shift - b_shift * x_n - c_n;
  rr = rr.expand();

  // Builds the lists to put in the matrix `rr_coefficients' and `rhs'.
  GList equations;
  for (unsigned i = 0; i < number_of_unknowns; ++i) {
    GExpr lhs = rr.coeff(n, i);
    equations.prepend(ex(lhs == 0));
  }
    
  GExpr solution = lsolve(equations, unknowns);  
  if (solution.nops() == 0)
    return false;

  // Builds the solution `x(n)'.
  for (unsigned i = 0; i < number_of_coeffs; ++i)
    x_n = x_n.subs(unknowns.op(i) == solution.op(i).op(1));

#if NOISY
  std::cout << "Solution x(n) = " << x_n << std::endl;
#endif  
  return true;
}

//! Gosper's algorithm, step 3: see \S5.4 of \f$ A = B \f$, by 
//! M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
/*!
  Returns <CODE>true</CODE> if it finds, if exists, a non-zero polynomial
  solution \f$ x(n) \f$ of \f$ a(n) * x(n+1) - b(n-1) * x(n) = c(n) \f$,
  where /f$ a(n) /f$, /f$ b(n) /f$ and /f$ c(n) /f$ are polynomials such that
  /f$ gcd(a(n), b(n+h)) = 1 \f$, for all non-negative integers
  /f$ h /f$. The solution /f$ x(n) /f$ is stored in /p x_n. 
  Returns <CODE>false</CODE> otherwise, i. e., it not finds a non-zero
  polynomial solution /f$ x(n) /f$. In this case
  /f$ x(n) = \sum_{k=0}^{n-1} t_k /f$.
*/
static bool
gosper_step_three(const GExpr& a_n, const GExpr& b_n, const GExpr& c_n,
		  const GSymbol& n, GExpr& x_n) {
  // Gosper's algorithm, step 3.1.
  // Finds the degree of `x(n)'.
  unsigned deg_a = a_n.degree(n);
  unsigned deg_b = b_n.degree(n);
  unsigned deg_c = c_n.degree(n);
  GExpr lead_a = a_n.lcoeff(n);
  GExpr lead_b = b_n.lcoeff(n);
  GNumber deg_x = -1;
  // On output, if a possible degree exixts,
  // `deg_x' will contain a non-negative number.
  if (deg_a != deg_b || !lead_a.is_equal(lead_b)) {
    if (deg_a > deg_b && deg_c >= deg_a)
      deg_x = deg_c - deg_a;
    if (deg_b > deg_a && deg_c >= deg_b)
      deg_x = deg_c - deg_b;
  }
  else {
    // `deg_a = deg_b' and `lead_a = lead_b'.
    GExpr shift_b = b_n.subs(n == n - 1);
    GExpr A = a_n.coeff(n, deg_a - 1);
    GExpr B = shift_b.coeff(n, deg_a - 1);
    assert(is_a<numeric>((B - A) * pow(lead_a, -1)));
    GNumber B_A = ex_to<numeric>((B - A) * pow(lead_a, -1));
    if (!B_A.is_nonneg_integer())
      B_A = -1;
    GNumber possible_deg;
    if (deg_c + 1 < deg_a)
      possible_deg = -1;
    if ((possible_deg - B_A).is_pos_integer())
      deg_x = possible_deg;
    else
      deg_x = B_A;
  }
  if (deg_x.is_equal(-1))
    return false;
#if NOISY
  std::cout << "Degree of x(n) = " << deg_x << std::endl;
#endif

  // Gosper's algorithm, step 3.2.
  return find_polynomial_solution(n, deg_x, a_n, b_n, c_n, x_n);
}

/*!
  Gosper's algorithm, step 4: see Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/
static GExpr
gosper_step_four(const GExpr& t, const GExpr& b_n, const GExpr& c_n,
		 const GExpr& x_n, const GSymbol& n,
		 const int lower_bound, const GExpr& upper_bound,
		 GExpr solution) {
  GExpr shift_b = b_n.subs(n == n-1);
  GExpr z_n = shift_b * x_n * t * pow(c_n, -1);
  z_n = simplify_numer_denom(z_n);
  // The Gosper's algorithm computes summation with the lower bound `0'
  // and the upper bound `n - 1': in this case, once we have `z_n',
  // the sum that we are looking for is `z_n - z_0'.
  // In general the solution will be `z_n - z_{lower_bound}'.
  solution = z_n - z_n.subs(n == lower_bound);
  // We must modify the sum if its upper bound is not `n - 1'.
  if (!upper_bound.is_equal(n-1))
    if (upper_bound.is_equal(n))
      solution = solution.subs(n == n + 1);
    else {
      GExpr n_plus_i = n + wild(0);
      assert(upper_bound.is_equal(n_plus_i));
      GList substitution;
      match(upper_bound, n_plus_i, substitution);
      GExpr i = get_binding(substitution, 0);
      solution = solution.subs(n == n + i + 1);
    }
  solution = simplify_numer_denom(solution);
  return solution;
}

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/
// FIXME: `r_n' is temporary until the implementation of step one that
// will build `r_n' starting from `t'.
bool
gosper(const GExpr& t, GExpr& r_n, const GSymbol& n,
       const int lower_bound, const GExpr& upper_bound, GExpr& solution) {
  if (!gosper_step_one(t, n, r_n))
    // `t' is not hypergeometric: no chance of using Gosper's algorithm.
    return false;
  GExpr a_n;
  GExpr b_n;
  GExpr c_n;
  gosper_step_two(r_n, n, a_n, b_n, c_n);
  GExpr x_n;
  if (!gosper_step_three(a_n, b_n, c_n, n, x_n)) {
    std::cout << "No non-zero polynomial solution" << std::endl;
    return false;
  }
  solution = gosper_step_four(t, b_n, c_n, x_n, n, lower_bound, upper_bound,
			      solution);
#if NOISY
  std::cout << "The sum is: " << solution << std::endl;
#endif
  return true;
}
