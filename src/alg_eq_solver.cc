/* Definition of the main function of the algebraic equation solver.
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
#include "alg_eq_solver.hh"
#include <cassert>
#include <iostream>

using namespace GiNaC;

/*!
  ...
*/
static const unsigned FIND_DIVISORS_MAX = 11;

static const unsigned
FIND_DIVISORS_THRESHOLD = FIND_DIVISORS_MAX*FIND_DIVISORS_MAX;

static GExpr zero = 0;

/*!
  This routine inserts into \p divisors all the positive divisors of
  the strictly positive integer \p n which is also less than
  <CODE>FIND_DIVISORS_THRESHOLD</CODE>.
*/
static void
find_divisors(GNumber n, std::vector<GNumber>& divisors) {
  assert(n.is_pos_integer());
  assert(n > 0 && n < FIND_DIVISORS_THRESHOLD);
  unsigned m = n.to_int();
  // Once a divisor `i' is found, it is pushed onto the vector divisors,
  // along with its conjugate `j = n/i', provided that `j' is less than `i'.
  if (m == 1)
    divisors.push_back(1);
  else
    for (unsigned i = 1, j = m; i < FIND_DIVISORS_MAX && i < j; ++i) {
      j = m / i;
      unsigned r = m % i;
      if (r == 0) {
	divisors.push_back(i);
	if (i < j)
	  divisors.push_back(j);
      }
    }
}

static bool
find_roots(const GExpr& p, const GSymbol& x,
	   std::vector<Polynomial_Root>& roots, GNumber multiplicity = 1);

static bool
find_power_roots(const GExpr& p, const GSymbol& x,
		 std::vector<Polynomial_Root>& roots);

static void
solve_equation_2(const GExpr& b, const GExpr& c,
		 GExpr& x1, GExpr& x2);

bool
solve_equation_3(const GNumber& a1, const GNumber& a2, const GNumber& a3,
		 GExpr& x1, GExpr& x2, GExpr& x3);

void
solve_equation_4(const GNumber& a1, const GNumber& a2,
		 const GNumber& a3, const GNumber& a4,
		 GExpr& x1, GExpr& x2, GExpr& x3, GExpr& x4);

/*!
  Let \p p be a polynomial with integer coefficients in \p x and
  \p roots be a (possibly non-empty) vector.
  This function appends to \p roots some or all the roots of \p p.
  Let \f$n\f$ be the degree of \p p and let \f$h\f$ and \f$k\f$
  be the value of <CODE>roots.size()</CODE> on entry and on exit
  to and from find_roots(), respectively.
  If \f$n = k-h\f$, then the positions \f$h, h+1, \ldots, k-1\f$
  of \p roots contain <EM>all</EM> the (possibly complex) roots of \p p
  and the function returns <CODE>true</CODE>.
  If \f$n \neq k-h\f$, then the positions \f$h, h+1, \ldots, k-1\f$
  of \p roots contain <EM>some</EM> of the roots of \p p and the function
  returns <CODE>true</CODE> if it was possible to prove that <EM>all</EM>
  the roots of \p p of maximal modulus are among those inserted
  into \p roots;
  the function returns <CODE>false</CODE> otherwise.
  The parameter \p all_distinct is set to <CODE>true</CODE> if it was
  possible to prove that all the roots of \p p are distinct;
  \p all_distinct is set to <CODE>false</CODE> otherwise.
*/
bool
find_roots(const GExpr& p, const GSymbol& x,
	   std::vector<Polynomial_Root>& roots,
	   bool& all_distinct) {
  assert(p.info(info_flags::integer_polynomial));
  assert(!p.info(info_flags::numeric));
  // Compute a square-free decomposition for p.
  GExpr q = sqrfree(p.expand(), lst(x));
  // There are now 4 cases depending on q's principal functor:
  //
  // 1) q is a product of two or more factors: e.g., (1+x)^2*(2+x);
  // 2) q is a factor to the power of its multiplicity (at least 2)
  //    in q: e.g., (1+x)^2;
  // 3) q is a sum of monomials: e.g., 1+x+x^2;
  // 4) q is the symbol x.
  //
  // In cases 1 and 2 there are multiple roots,
  // in cases 3 and 4 there are not.
  if (is_a<add>(q)) {
    all_distinct = true;
    return find_roots(q, x, roots);
  }
  else if (is_a<mul>(q)) {
    all_distinct = false;
    for (unsigned i = 0, n = q.nops(); i < n; ++i) {
      GExpr factor = q.op(i);
      if (is_a<power>(factor)) {
	if (!find_power_roots(factor, x, roots))
	  return false;
      }
      else if (!find_roots(factor, x, roots))
	return false;
    }
    return true;
  }
  else if (is_a<power>(q)) {
    all_distinct = false;
    return find_power_roots(q, x, roots);
  }
  else {
    assert(is_a<symbol>(q));
    // 0 is the only solution of x = 0.
    all_distinct = true;
    roots.push_back(zero);
    return true;
  }
}

static bool
find_power_roots(const GExpr& p, const GSymbol& x,
		 std::vector<Polynomial_Root>& roots) {
  assert(is_a<power>(p));
  GExpr base = p.op(0);
  assert(is_a<numeric>(p.op(1)));
  GNumber exponent = ex_to<numeric>(p.op(1));
  assert(exponent.is_pos_integer() && exponent >= 2);
  if (!find_roots(base, x, roots, exponent))
    // No way: we were unable to solve the base.
    return false;
  return true;
}

static bool
find_roots(const GExpr& p, const GSymbol& x,
	   std::vector<Polynomial_Root>& roots,
	   GNumber multiplicity) {
  int ldegree = p.ldegree(x);
  assert(ldegree <= 1);
  GExpr q;
  if (ldegree == 1) {
    roots.push_back(Polynomial_Root(zero, multiplicity));
    q = quo(p, x, x);
  }
  else
    q = p;

  assert(is_a<numeric>(q.lcoeff(x)));
  assert(is_a<numeric>(q.tcoeff(x)));
  GNumber lc = ex_to<numeric>(q.lcoeff(x));
  GNumber tc = ex_to<numeric>(q.tcoeff(x));
  int degree = q.degree(x);
  if (degree == 1) {
    roots.push_back(Polynomial_Root(-tc/lc, multiplicity));
    return true;
  }

  // Here `q' has degree at least 2 and the least coefficient is non-zero.
  GNumber abs_lc = abs(lc);
  GNumber abs_tc = abs(tc);
  if (abs_lc < FIND_DIVISORS_THRESHOLD && abs_tc < FIND_DIVISORS_THRESHOLD) {
    bool coefficients_changed = false;
    std::vector<GNumber> abs_lc_divisors;
    std::vector<GNumber> abs_tc_divisors;
    find_divisors(abs_lc, abs_lc_divisors);
    find_divisors(abs_tc, abs_tc_divisors);
    for (unsigned l = 0, ml = abs_lc_divisors.size(); l < ml; ++l) 
      for (unsigned t = 0, mt = abs_tc_divisors.size(); t < mt; ++t) {
	GNumber r = abs_tc_divisors[t] / abs_lc_divisors[l];
	if (q.subs(x == r) == 0) {
	  q = quo(q, x-r, x);
	  --degree;
	  roots.push_back(Polynomial_Root(r, multiplicity));
	  coefficients_changed = true;
	}
	r = -r;
	if (q.subs(x == r) == 0) {
	  q = quo(q, x-r, x);
	  --degree;
	  roots.push_back(Polynomial_Root(r, multiplicity));
	  coefficients_changed = true;
	}
      }
    if (degree == 0)
      // FIXME: a decent comment here.
      return true;

    // Here `q' has only simple roots which are either irrational or complex.
    if (coefficients_changed) {
      assert(is_a<numeric>(q.lcoeff(x)));
      assert(is_a<numeric>(q.tcoeff(x)));
      lc = ex_to<numeric>(q.lcoeff(x));
      tc = ex_to<numeric>(q.tcoeff(x));
    }
  }

  assert(degree > 1);
  if (degree <= 4) {
    unsigned position = roots.size();
    // Insert `degree' elements at the end of roots.
    roots.insert(roots.end(), degree, 0);

    switch (degree) {
    case 2:
      {
	assert(is_a<numeric>(q.coeff(x, 1)));
	GNumber b = ex_to<numeric>(q.coeff(x, 1)) / lc;
	GNumber c = tc / lc;
	solve_equation_2(b, c,
			 roots[position].value(),
			 roots[position+1].value());
	return true;
      }
    case 3:
      {
	assert(is_a<numeric>(q.coeff(x, 1)));
	assert(is_a<numeric>(q.coeff(x, 2)));
	GNumber a1 = ex_to<numeric>(q.coeff(x, 2)) / lc;
	GNumber a2 = ex_to<numeric>(q.coeff(x, 1)) / lc;
	GNumber a3 = tc / lc;
	solve_equation_3(a1, a2, a3,
			 roots[position].value(),
			 roots[position+1].value(),
			 roots[position+2].value());
	return true;
      }
    case 4:
      {
	assert(is_a<numeric>(q.coeff(x, 1)));
	assert(is_a<numeric>(q.coeff(x, 2)));
	assert(is_a<numeric>(q.coeff(x, 3)));
	GNumber a1 = ex_to<numeric>(q.coeff(x, 3)) / lc;
	GNumber a2 = ex_to<numeric>(q.coeff(x, 2)) / lc;
	GNumber a3 = ex_to<numeric>(q.coeff(x, 1)) / lc;
	GNumber a4 = tc / lc;
	solve_equation_4(a1, a2, a3, a4,
			 roots[position].value(),
			 roots[position+1].value(),
			 roots[position+2].value(),
			 roots[position+3].value());
	return true;
      }
    default:
      abort();
      break;
    }
  }
  return false;
}


//! Solve the equation \f$x^2 + b x + c = 0\f$.
static void
solve_equation_2(const GExpr& b, const GExpr& c,
		 GExpr& x1, GExpr& x2) {
  GExpr sqrt_d = sqrt(b*b - 4*c);
  x1 = (-b + sqrt_d)/2;
  x2 = (-b - sqrt_d)/2;
}

static GExpr
cubic_root(const GExpr& e) {
  static GExpr one_third = GExpr(1)/3;
  if (is_a<numeric>(e)) {
    GNumber n = ex_to<numeric>(e); 
    if (n.is_rational() && n < 0)
      return (-pow(-n, one_third));
  }
  std::cout << "Computing cubic_root(" << e << ")" << std::endl;
  return pow(e, one_third);
}

/*!
  Solve the equation \f$x^3 + a_1 x^2 + a_2 x + a_3 = 0\f$
  and return <CODE>true</CODE> if and only if all the solutions are real.
  \f$x_1\f$ is guaranteed to be a real solution.
*/
bool
solve_equation_3(const GNumber& a1, const GNumber& a2, const GNumber& a3,
		 GExpr& x1, GExpr& x2, GExpr& x3) {
  GExpr Q = (3*a2 - pow(a1, 2)) / 9;
  GExpr R = (9*a1*a2 - 27*a3 -2*pow(a1, 3)) / 54;
  GExpr d = pow(Q, 3) + pow(R, 2);
  GExpr a1_div_3 = a1/3;
  if (d < 0) {
    GExpr sqrt_minus_Q = sqrt(-Q);
    GExpr theta = acos(-R/(Q*sqrt_minus_Q));
    x1 = -a1_div_3 + 2*sqrt_minus_Q*cos(theta/3);
    x2 = -a1_div_3 + 2*sqrt_minus_Q*cos((theta+2*Pi)/3);
    x3 = -a1_div_3 + 2*sqrt_minus_Q*cos((theta+4*Pi)/3);
  }
  else {
    GExpr sqrt_d = sqrt(d);
    GExpr S = cubic_root(R + sqrt_d);
    GExpr T = cubic_root(R - sqrt_d);
    GExpr S_plus_T = S + T;
    GExpr t1 = -S_plus_T/2 - a1_div_3;
    GExpr t2 = (S - T)*I*sqrt(GExpr(3))/2;
    x1 = S_plus_T - a1_div_3;
    x2 = t1 + t2;
    x3 = t1 - t2;
  }
  // The roots are all real if and only if d <= 0.
  return d <= 0;
}

//! Solve the equation \f$x^4 + a_1 x^3 + a_2 x^2 + a_3 x + a_4 = 0\f$.
void
solve_equation_4(const GNumber& a1, const GNumber& a2,
		 const GNumber& a3, const GNumber& a4,
		 GExpr& x1, GExpr& x2, GExpr& x3, GExpr& x4) {
  GExpr y1;
  GExpr y2;
  GExpr y3;
#if 0
  bool all_real =
#endif
    solve_equation_3(-a2,
		     a1*a3 - 4*a4,
		     4*a2*a4 - a3*a3 - a1*a1*a4,
		     y1, y2, y3);
  GExpr d1 = pow(a1, 2) - 4*a2  + 4*y1;
  GExpr d2 = pow(y1, 2) - 4*a4;
  GExpr sqrt_d1 = sqrt(d1);
  GExpr sqrt_d2 = sqrt(d2);
  solve_equation_2((a1 + sqrt_d1)/2, (y1 - sqrt_d2)/2, x1, x2);
  solve_equation_2((a1 - sqrt_d1)/2, (y1 + sqrt_d2)/2, x3, x4);
}
