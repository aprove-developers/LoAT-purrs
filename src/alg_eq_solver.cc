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

#define NOISY 0
 
#include "globals.hh"
#include "alg_eq_solver.hh"
#include "poly_factor.hh"
#include "simplify.hh"
#include "util.hh"
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
  // Once a divisor `i' is found, it is pushed onto the vector `divisors'
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

/*!
  This function takes a polynomial expression \f$p(x)\f$ and returns 
  the largest integer \f$n\f$ such that there is a polynomial \f$q\f$
  such that \f$p(x) = q(x^n)\f$ and the polynomial \f$q\f$ itself.
*/
static unsigned
is_nested_polynomial(const GExpr& p, const GSymbol& x, GExpr& q) {
  unsigned degree = p.degree(x);  
  if (degree == 0) {
    // The constant polynomial.
    q = p;
    return 0;
  }
  // Here the degree is at least 1.
  if (p.coeff(x, 1) != 0) {
    // The gcd of the coefficients is 1.
    q = p;
    return 1;
  }
  // Here the degree is at least 2 and the polynomial does not have a
  // linear term.  Look for the first non-zero coefficient (apart from
  // the constant term).
  unsigned i = 2;
  while(p.coeff(x, i) == 0)
    ++i;

  // Here i >= 2 and the polynomial has the form a_0 + a_i x^i + ... 

  // Check whether all exponents are multiple of some integer `n'.
  // We first set n = i, and update its value every time that we 
  // find a non-zero coefficient.
  // The routine ends as soon as `n' reaches the value of 1
  // (this means that the gcd of all exponents of non-zero 
  // monomials is 1) or when the polynomial has been entirely processed.
  unsigned n = i;
  for (unsigned j = i+1; j <= degree && n > 1; ++j)
    // If n ==1 there is no need to read the rest of the polynomial.
    if (p.coeff(x, j) !=0)
      n = gcd(n, j);

  // Here n is the largest integer such that there is 
  // a polynomial q such that p(x) == q(x^n). 
  // Now we compute q.
  if (n > 1) {
    q = p.coeff(x, 0);
    // Note that `n' divides `degree'.
    for (unsigned j = 1, m = degree/n; j <= m; ++j)
      q += p.coeff(x, n*j) * pow(x, j); 
  }
  else
    // n == 1, the polynomial q is equal to the polynomial p.
    q = p;
  return n;
}

static bool
find_roots(const GExpr& p, const GSymbol& x,
	   std::vector<Polynomial_Root>& roots, GNumber multiplicity);

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
  // There are now 4 cases depending on the principal functor of `q':
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
    return find_roots(q, x, roots, 1);
  }
  else if (is_a<mul>(q)) {
    all_distinct = false;
    for (unsigned i = 0, n = q.nops(); i < n; ++i) {
      GExpr factor = q.op(i);
      if (is_a<power>(factor)) {
	if (!find_power_roots(factor, x, roots))
	  return false;
      }
      else if (!find_roots(factor, x, roots, 1))
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
      // The polynomial has been factored completely.
      return true;

    // Here `q' has only simple roots that are either
    // - rational with large (>= FIND_DIVISORS_THRESHOLD)
    //   numerator or denominator,
    // - irrational, or
    // - complex.
    if (coefficients_changed) {
      assert(is_a<numeric>(q.lcoeff(x)));
      assert(is_a<numeric>(q.tcoeff(x)));
      lc = ex_to<numeric>(q.lcoeff(x));
      tc = ex_to<numeric>(q.tcoeff(x));
    }
  }

  // Direct solution for polynomials of degree between 1 and 4.
  if (degree <= 4) {
    unsigned position = roots.size();
    // Insert `degree' elements at the end of roots.
    roots.insert(roots.end(), degree, 0);

    switch (degree) {
    case 1:
      {
	roots.push_back(Polynomial_Root(-tc/lc, multiplicity));
	return true;
      }
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
	// FIXME: call `is_nested_polynomial(q, x, r)' here?
	// Consider, for example, 1+x^4+x^2 = 0.
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
    }
  }

  // If we want to solve q(x) = 0 and we know that q(x) = r(x^n), then
  // we need to find the roots y_1, ... y_k of r(y) = 0, and then we
  // need to solve x^n = y_1, x^n = y_2, ..., x^n = y_k.
  // Once a root of x^n = y_1, call it x_1, has been found, all the
  // roots of the equation x^n = y_1 can be found by multiplying x_1
  // by the n-th roots of unity, that is by the complex numbers z such
  // that z^n = 1.
  // We note that the set of roots found in this way does not depend
  // on the particular root x_1 that we have chosen.
  GExpr r;
  int nested_degree = is_nested_polynomial(q, x, r);
  if (nested_degree > 1) {
    // We need a vector to hold the roots of `r'.
    std::vector<Polynomial_Root> roots_r;
    // Avoid any useless reallocation.
    roots_r.reserve(r.degree(x));
    if (find_roots(r, x, roots_r, 1)) {
      size_t num_r_roots = roots_r.size();
      GExpr theta = 2*Pi/nested_degree;
      for (int j = 0; j < nested_degree; ++j) {
	GExpr root_of_unity = cos(j*theta) + I*sin(j*theta);
	for (size_t i = 0; i < num_r_roots; ++i)
	  roots.push_back(Polynomial_Root(pow(roots_r[i].value(),
					      numeric(1)/nested_degree)
					  * root_of_unity, multiplicity));
      }
      return true;
    }
  }

  // Try to factorize the polynomial.
  std::vector<GExpr> factors;
  int num_factors = poly_factor(q, x, factors);
  if (num_factors > 1) {
    for (int i = num_factors-1; i >= 0; --i)
      if (!find_roots(factors[i], x, roots, 1))
	return false;
    return true;
  }
  return false;
}


//! Solve the equation \f$x^2 + b x + c = 0\f$.
static void
solve_equation_2(const GExpr& b, const GExpr& c,
		 GExpr& x1, GExpr& x2) {
  GSymbol n("n");
  GExpr sqrt_d = sqrt(b*b - 4*c);

  D_MSGVAR("Before: ", sqrt_d);

  sqrt_d = simplify_on_output_ex(sqrt_d, n, false);

  D_VAR(sqrt_d);

  x1 = (-b + sqrt_d)/2;
  x2 = (-b - sqrt_d)/2;

  x1 = simplify_on_output_ex(x1, n, false);
  x2 = simplify_on_output_ex(x2, n, false);
}

/*!
  Solve the equation \f$x^3 + a_1 x^2 + a_2 x + a_3 = 0\f$
  and return <CODE>true</CODE> if and only if all the solutions are real.
  \f$x_1\f$ is guaranteed to be a real solution.
*/
bool
solve_equation_3(const GNumber& a1, const GNumber& a2, const GNumber& a3,
		 GExpr& x1, GExpr& x2, GExpr& x3) {
  GSymbol n("n");
  GNumber Q = (3*a2 - a1*a1) / 9;
  GNumber R = (9*a1*a2 - 27*a3 -2*a1*a1*a1) / 54;
  GNumber d = Q*Q*Q + R*R;
  GNumber a1_div_3 = a1/3;
  if (d < 0) {
    GExpr sqrt_minus_Q = sqrt(-Q);
    GExpr theta = acos(-R/(Q*sqrt_minus_Q));
    x1 = -a1_div_3 + 2*sqrt_minus_Q*cos(theta/3);
    x2 = -a1_div_3 + 2*sqrt_minus_Q*cos((theta+2*Pi)/3);
    x3 = -a1_div_3 + 2*sqrt_minus_Q*cos((theta+4*Pi)/3);
  }
  else {
    // When d > 0 there is one real and two complex conjugate roots.
    // In order to avoid taking cube roots of negative numbers we
    // check the signs of the expressions A and B below (which need
    // not be rational) without actually computing them:
    // - if Q >= 0 then d >= R^2 and sqrt_d >= abs(R) so that 
    //   A >= 0 and B < 0;
    // - if Q < 0 and R < 0 then a similar argument shows that 
    //   A < 0 and B < 0; 
    // - if Q < 0 and R >= 0 then A > 0 and B > 0.
    GExpr sqrt_d = sqrt(d);
    GExpr A = R + sqrt_d;
    GExpr B = R - sqrt_d;
    GExpr S;
    GExpr T;
    if (Q >= 0) {
      S = cubic_root(A);
      T = -cubic_root(-B);
    }
    else {
      // Q < 0
      if (R < 0) {
	S = -cubic_root(-A);
	T = -cubic_root(-B);
      }
      else {
	// R >= 0
	S = cubic_root(A);
	T = cubic_root(B);
      }
    }

    D_MSGVAR("Before: ", S); 
    D_MSGVAR("Before: ", T);

    S = simplify_on_output_ex(S, n, false);
    T = simplify_on_output_ex(T, n, false);

    D_VAR(S); 
    D_VAR(T);

    GExpr S_plus_T = S + T;

    // FIXME: S+T are of the form (a+b)^(1/3) + (a-b)^(1/3).
    // Is there a way to simplify this?

    GExpr t1 = -S_plus_T/2 - a1_div_3;
    GExpr t2 = (S - T)*I*sqrt(GExpr(3))/2;
    x1 = S_plus_T - a1_div_3;
    x2 = t1 + t2;
    x3 = t1 - t2;
  }

  x1 = simplify_on_output_ex(x1, n, false);
  x2 = simplify_on_output_ex(x2, n, false);
  x3 = simplify_on_output_ex(x3, n, false);

  // The roots are all real if and only if d <= 0.
  return d <= 0;
}

//! Solve the equation \f$x^4 + a_1 x^3 + a_2 x^2 + a_3 x + a_4 = 0\f$.
void
solve_equation_4(const GNumber& a1, const GNumber& a2,
		 const GNumber& a3, const GNumber& a4,
		 GExpr& x1, GExpr& x2, GExpr& x3, GExpr& x4) {
  GNumber f = a2 - 3*a1*a1*numeric(1)/8;
  GNumber g = a3 + a1*a1*a1/8 - a1*a2/2;
  GNumber h = a4 - 3*a1*a1*a1*a1/256 + a1*a1*a2/16 - a1*a3/4;

  D_VAR(f); 
  D_VAR(g); 
  D_VAR(h);

  GExpr y1;
  GExpr y2;
  GExpr y3;
  // If 'g' is zero then the auxiliary equation
  // y^3 +  f/2*y^2 + (f*f - 4*h)/16*y - g*g/64 = 0
  // has the root 0, in this case we solve the simpler auxiliary equation
  // y^2 +  f/2*y + (f*f - 4*h)/16 = 0.
  if (g == 0) {
    solve_equation_2(f/2, (f*f - 4*h)/16, y1, y2);
    // Both roots are nonzero.
    assert(y1 != 0 && y2 != 0);
  }
  else
    // We are deliberately ignoring the return value.
    (void) solve_equation_3(f/2, (f*f - 4*h)/16, -g*g/64, y1, y2, y3);

  GExpr p, q;
  p = sqrt(y1);
  q = sqrt(y2);

  D_MSGVAR("Before: ", p); 
  D_MSGVAR("Before: ", q);

  // FIXME: the one and only `n' symbol should be global,
  // i.e., created once and for all.
  GSymbol n("n");
  p = simplify_on_output_ex(p, n, false);
  q = simplify_on_output_ex(q, n, false);

  D_VAR(p); 
  D_VAR(q);

  GExpr r = numeric(-g)/(8*p*q);

  D_MSGVAR("Before: ", r); 

  r = simplify_on_output_ex(r, n, false);
  GExpr s = numeric(a1)/4;

  D_VAR(r); 
  D_VAR(s); 

  x1 = p + q + r - s;
  x2 = p - q - r - s;
  x3 = -p + q - r - s;
  x4 = -p - q + r - s;

  D_MSG("Solutions before calling simplify: ");
  D_VAR(x1); 
  D_VAR(x2); 
  D_VAR(x3); 
  D_VAR(x4);
  D_MSG("");

  x1 = simplify_on_output_ex(x1, n, false);
  x2 = simplify_on_output_ex(x2, n, false);
  x3 = simplify_on_output_ex(x3, n, false);
  x4 = simplify_on_output_ex(x4, n, false);
}

std::ostream&
operator<<(std::ostream& s, const Polynomial_Root& r) {
  GNumber multiplicity = r.multiplicity();
  if (multiplicity != 1)
    s << "mult: " << multiplicity << ", val: ";
  s << r.value();
  return s;
}
