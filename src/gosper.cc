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

#include "util.hh"
#include "alg_eq_solver.hh"
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
  if \f$ p(x+1) / p(x) \f$ is a rational function of \f$ x \f$.
  If this is the case, return also numerator and denominator.
*/
// static bool
// is_hypergeometric_term(const GExpr& p, const GSymbol& x,
// 		       GExpr& num_den) {

//   GExpr q = (p.subs(x == x+1)) / p;

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
gosper_step_one(const GExpr& t, const GSymbol& /* n */, GExpr& num_den_t) {
  // FIXME: general simplifications to be inserted here.
  num_den_t = numer_denom(t);
  return true;
}

/*!
  The resultant's non-negative integral roots are needed in Step 2.1
  of Gosper's algorithm.
*/
static void
compute_resultant_and_its_roots(const GExpr& f, const GExpr& g, 
				const GSymbol& x,
				std::vector<GNumber>& integer_roots) {
  GSymbol h("h");
  GExpr temp_g = g.subs(x == x + h);
  GExpr R = resultant(f, temp_g, x);
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
gosper_step_two(const GExpr& n_d, const GSymbol& x,
		GExpr& p, GExpr& q, GExpr& r) {
  // Gosper's algorithm, step 2.1.
  GExpr f = expand(n_d.op(0)); // the numerator
  GExpr g = expand(n_d.op(1)); // the denominator

  std::vector<GNumber> integer_roots;
  compute_resultant_and_its_roots(f, g, x, integer_roots);

  // Gosper's algorithm, step 2.2.
  // `p' and `q' are used below as starting values for a sequence
  // of polynomials.

  // normalize f and g, and store conversion factor in Z
  GExpr lead_f = f.lcoeff(x);
  GExpr lead_g = g.lcoeff(x);
  GExpr Z = lead_f * pow(lead_g, -1);
  p = f * pow(lead_f, -1);
  q = g * pow(lead_g, -1);

  // Computation of the output polynomials.
  r = 1;
  unsigned integer_roots_size = integer_roots.size();
  for (unsigned i = 0; i < integer_roots_size; ++i) {
    GExpr temp_q = (q.subs(x == x + integer_roots[i])).expand();
    GExpr s = general_gcd(p, temp_q, x);
    p = quo(p, s, x);
    GExpr temp_s = s.subs(x == x - integer_roots[i]);
    q = quo(q, temp_s, x);
    for (GNumber j = 1; j <= integer_roots[i]; ++j)
      r *= s.subs(x == x - j);
  }
  p *= Z;
  // The polynomials `p' and `q' may have rational coefficients.
  // Multiply numerator and denominator of the fraction `p/q'
  // by suitable integers, so that the output values of `p' and `q'
  // have integer coefficients. 
  GNumber p_factor;
  p = convert_to_integer_polynomial(p, x, p_factor);
  GExpr p_n_d = numer_denom(p_factor);
  p *= p_n_d.op(0);
  q *= p_n_d.op(1);
  GNumber q_factor;
  q = convert_to_integer_polynomial(q, x, q_factor);
  GExpr q_n_d = numer_denom(q_factor);
  p *= q_n_d.op(1);
  q *= q_n_d.op(0);
#if NOISY
  std::cout << "a(x) = " << p << std::endl;
  std::cout << "b(x) = " << q << std::endl;
  std::cout << "c(x) = " << r << std::endl;
#endif
}

static bool
find_polynomial_solution(const GSymbol& x, const GNumber& deg_x,
			 const GExpr& a_x, const GExpr& b_x, const GExpr& c_x,
			 GExpr& x_n) {
  unsigned deg_a = a_x.degree(x);
  unsigned deg_b = b_x.degree(x);
  unsigned deg_c = c_x.degree(x);
  // `number_of_coeffs' is the number of coefficients of 
  // the polynomial `p', that is, 1 + deg_x.
  // FIXME: is it possible to convert a GiNaC integer to an int or an unsigned?
  // For the moment, we use a dirty trick.
  unsigned number_of_coeffs = 0;
  for (unsigned i = 0; i <= deg_x; ++i)
    ++number_of_coeffs;

  // Compute the real number of unknowns of the polynomial equation
  // `a(x) * p(x+1) - b(x-1) * p(x) - c(x) = 0'.
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
    x_n += pow(x, i) * unknowns.op(i);

  GExpr x_n_shift = x_n.subs(x == x+1);
  GExpr b_shift = b_x.subs(x == x-1);

  // Considers the recurrence relation to solve.
  GExpr rr = a_x * x_n_shift - b_shift * x_n - c_x;
  rr = rr.expand();

  // Builds the lists to put in the matrix `rr_coefficients' and `rhs'.
  GList equations;
  for (unsigned i = 0; i < number_of_unknowns; ++i) {
    GExpr lhs = rr.coeff(x, i);
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
gosper_step_three(const GExpr& a_x, const GExpr& b_x, const GExpr& c_x,
		  const GSymbol& x, GExpr& x_n) {
  // Gosper's algorithm, step 3.1.
  // Finds the degree of `x(n)'.
  unsigned deg_a = a_x.degree(x);
  unsigned deg_b = b_x.degree(x);
  unsigned deg_c = c_x.degree(x);
  GExpr lead_a = a_x.lcoeff(x);
  GExpr lead_b = b_x.lcoeff(x);
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
    GExpr shift_b = b_x.subs(x == x - 1);
    GExpr A = a_x.coeff(x, deg_a - 1);
    GExpr B = shift_b.coeff(x, deg_a - 1);
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
  return find_polynomial_solution(x, deg_x, a_x, b_x, c_x, x_n);
}

/*!
  Gosper's algorithm, step 4: see Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/
static GExpr
gosper_step_four(const GExpr& /* b_n */, const GExpr& /* c_n*/,
		 const GExpr& x_n, const GSymbol& /* n */) {
  // FIXME: to do
  return x_n;
}

/*!
  Gosper's algorithm, from Chapter 5 of \f$ A = B \f$, by 
  M.~Petkov\v sek, H.~Wilf and D.~Zeilberger.
*/
bool
gosper(const GExpr& t, const GSymbol& n, GExpr& solution) {
  GExpr n_d;
  if (!gosper_step_one(t, n, n_d))
    // `t' is not hypergeometric: no chance of using Gosper's algorithm.
    return false;
  GExpr a_n;
  GExpr b_n;
  GExpr c_n;
  gosper_step_two(n_d, n, a_n, b_n, c_n);
  GExpr x_n;
  if (!gosper_step_three(a_n, b_n, c_n, n, x_n)) {
    std::cout << "No non-zero polynomial solution" << std::endl;
    return false;
  }
  solution = gosper_step_four(b_n, c_n, x_n, n); 
  return true;
}
