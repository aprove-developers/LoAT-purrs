/* Definition of some utility functions.
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

#include "Expr.defs.hh"
#include "Number.defs.hh"
#include "Symbol.defs.hh"
#include "util.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

bool
PURRS::vector_not_all_zero(const std::vector<PURRS::Expr>& v) {
  for (unsigned i = v.size(); i-- > 0; )
    if (!v[i].is_zero())
      return true;
  return false;
}

/*!
  Computes the gcd between the integers \p n and \p m.
*/
int
PURRS::gcd(int n, int m) {
  int r = m;
  while (r != 0){
    r = n % m;
    n = m; 
    m = r;
  }
  return n;  
}

/*!
  Accept rational polynomial as input, and normalize them 
  before calling GiNaC's <CODE>gcd()</CODE>.
*/
PURRS::Expr
PURRS::general_gcd(const Expr& p, const Expr& q, const Symbol& x) {
  Expr f = convert_to_integer_polynomial(p, x);
  Expr g = convert_to_integer_polynomial(q, x);
  return gcd(f,g);
}

/*!
  Computes the LCM among the numbers in the vector \p v.
*/
PURRS::Number
PURRS::lcm(const std::vector<Number>& v) {
  Number n = 1;
  for (unsigned i = v.size(); i-- > 0; )
    n = lcm(n, v[i]);
  return n;
}

PURRS::Expr
PURRS::cubic_root(const Expr& e) {
  static Expr one_third = Expr(1)/3;
  return pwr(e, one_third);
}

void
PURRS::clear(Expr_List& l) {
  for (unsigned n = l.nops(); n-- > 0; )
    l.remove_first();
  assert(l.nops() == 0);
}

/*!
  We assume that \p substitution has been produced by
  <CODE>GiNaC::match()</CODE> and that the binding for the wildcard of
  index \p wild_index is in the position \p wild_index of \p substitution.
*/
PURRS::Expr
PURRS::get_binding(const Expr_List& substitution, unsigned wild_index) {
  assert(wild_index < substitution.nops());
  assert(substitution.op(wild_index).is_relation_equal());
  //assert(substitution.op(wild_index).lhs() == wild(wild_index));
  return substitution.op(wild_index).rhs();
}


//! \brief
//! Isolates a polynomial part of \p e and assigns it to \p polynomial,
//! assigning the corresponding possibly non-polynomial part of \p e
//! to \p rest.
/*!
  Given an expression \f$ e \f$, <EM>isolating a polynomial part of
  \f$ e \f$ in \f$ x \f$</EM> means finding two expressions \f$ p \f$
  and \f$ r \f$ such that \f$ p \f$ is a polynomial in \f$ x \f$ and
  \f$ e = p + r \f$.
*/
void
PURRS::isolate_polynomial_part(const Expr& e, const Symbol& x,
			Expr& polynomial, Expr& rest) {
  if (e.is_a_add()) {
    polynomial = 0;
    rest = 0;
    for (unsigned i = e.nops(); i-- > 0; ) {
      if (e.op(i).is_polynomial(x))
	polynomial += e.op(i);
      else
	rest += e.op(i);
    }
  }
  else if (e.is_polynomial(x)) {
    polynomial = e;
    rest = 0;
  }
  else {
    rest = e;
    polynomial = 0;
  }
}


/*!
  A polynomial with integer coefficients is <EM>primitive</EM> if its
  leading coefficient is positive and its coefficients have no common factor.
  This function converts a polynomial with rational coefficients into the
  associate primitive polynomial divides the input polynomial.
*/
PURRS::Expr
PURRS::convert_to_integer_polynomial(const Expr& p, const Symbol& x) {
  assert(p.is_rational_polynomial());
  unsigned deg_p = p.degree(x);

  // Choose non-zero starting value and compute least common
  // multiple of denominators.
  Number t_lcm = p.coeff(x, deg_p).ex_to_number().denominator();
  for (unsigned i = 0; i <= deg_p; ++i) {
    Expr t_coeff = p.coeff(x, i);
    t_lcm = lcm(t_lcm, t_coeff.ex_to_number().denominator());
  }
  return (p * t_lcm).primpart(x);
}

/*!
  Converts a polynomial with rational coefficients into the 
  associate primitive polynomial with integer coefficients.
  This version also returns a Number containing the factor used to 
  convert.
*/
PURRS::Expr
PURRS::convert_to_integer_polynomial(const Expr& p, const Symbol& x,
                              Number& factor) {
  assert(p.is_rational_polynomial());
  unsigned deg_p = p.degree(x);

  // Choose non-zero starting value and compute least common
  // multiple of denominators.
  Number t_lcm = p.coeff(x, deg_p).ex_to_number().denominator();
  for (unsigned i = 0; i <= deg_p; ++i) {
    Expr t_coeff = p.coeff(x, i);
    t_lcm = lcm(t_lcm, t_coeff.ex_to_number().denominator());
  }

  Expr q = (p * t_lcm).primpart(x);
  factor  = p.lcoeff(x).ex_to_number();
  factor *= pwr(q.lcoeff(x), -1).ex_to_number();
  return q;
}

//! \brief
//! Computes the resultant of the polynomials \p f and \p g with 
//! rational coefficients, using Euclid's algorithm.
//! Returns the solution in \p res.
/*!
  The following properties of the resultant are used:
  - \f$ Res(g, f) = (-1)^{\deg(f)\deg(g)} Res(f, g) \f$; 
  - \f$ Res(f, g) = a^{\deg(g) - \deg(r)} Res(f, r) \f$
    if \f$ g = fq + r \f$ for some polynomials \f$ q \f$  and \f$ r \f$.
    Here \f$ a \f$ is the leading coefficient of the polynomial \f$ f \f$.
  - \f$ Res(f, b) = b^{\deg(f)} \f$ if \f$ b \f$ is a scalar.
*/
PURRS::Expr
PURRS::resultant(const Expr& p, const Expr& q, const Symbol& x) {
  D_VAR(p);
  D_VAR(q);
  assert(p.is_rational_polynomial());
  assert(q.is_rational_polynomial());
  Expr f = p.expand();
  Expr g = q.expand();
  Expr res = 1;
  unsigned deg_f = f.degree(x);
  unsigned deg_g = g.degree(x);

  // Special case: `f' or `g' is a constant polynomial. By definition
  // `Res(f, g) = f.lcoeff(n)^g.degree(n) * g.lcoeff(n)^f.degree(n)'. 
  if (deg_f == 0 || deg_g == 0)
    res = pwr(f.lcoeff(x), deg_g) * pwr(g.lcoeff(x), deg_f);
  else {
    // Modified Euclid's algorithm starts here.
    while (deg_f > 0) {
      // `prem()' computes the pseudo-remainder of `g' and `f' which satisfies
      // `factor * g = factor * f * q + prem(g, f, x)' where `q' is the
      // quozient of `g' and `f' and
      // `factor = f.lcoeff(x)^(g.degree(x) - f.degree(x) + 1)'.
      Expr r = prem(g, f, x);
      Expr factor = pwr(f.lcoeff(x), g.degree(x) - f.degree(x) + 1);
      // The rest of euclidean's division is given by the ratio
      // `pseudo-remainder / factor'.
      r *= pwr(factor, -1);
      unsigned deg_r = r.degree(x);
      Expr a = f.lcoeff(x);
      // Using rule two.
      res *= pwr(a, deg_g - deg_r);
      // Using rule one.
      if ((deg_f * deg_r) & 1 != 0)
	// `deg_f * deg_r' is odd.
	res = -res;
      g = f;
      f = r;
      deg_f = f.degree(x);
      deg_g = g.degree(x);
    }
    // Here `f' is a constant: use rule three.
    res *= pwr(f, deg_g);
  }
  D_MSGVAR("Resultant(f(x), g(x)): ", res);
  return res;
}
