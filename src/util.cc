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

#include "util.hh"

namespace Parma_Recurrence_Relation_Solver {

#ifndef NOISY
#define NOISY 0
#endif

/*!
  Computes the gcd between the integers \p n and \p m.
*/
int
gcd(int n, int m) {
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
Expr
general_gcd(const Expr& p, const Expr& q, const Symbol& x) {
  Expr f = convert_to_integer_polynomial(p, x);
  Expr g = convert_to_integer_polynomial(q, x);
  return gcd(f,g);
}

/*!
  Computes the LCM among the numbers in the vector \p v.
*/
Number
lcm(const std::vector<Number>& v) {
  Number n = 1;
  for (unsigned i = v.size(); i-- > 0; )
    n = lcm(n, v[i]);
  return n;
}

Expr
cubic_root(const Expr& e) {
  static Expr one_third = Expr(1)/3;
  return pwr(e, one_third);
}

void
clear(Expr_List& l) {
  for (unsigned n = l.nops(); n-- > 0; )
    l.remove_first();
  assert(l.nops() == 0);
}

/*!
  We assume that \p substitution has been produced by
  <CODE>GiNaC::match()</CODE> and that the binding for the wildcard of
  index \p wild_index is in the position \p wild_index of \p substitution.
*/
Expr
get_binding(const Expr_List& substitution, unsigned wild_index) {
  assert(wild_index < substitution.nops());
  assert(substitution.op(wild_index).is_relation_equal());
  //assert(substitution.op(wild_index).lhs() == wild(wild_index));
  return substitution.op(wild_index).rhs();
}

//! \brief
//! Returns <CODE>true</CODE> if \p e is a scalar rapresentation for \p x;
//! returns <CODE>false</CODE> otherwise.
/*!
  This function realizes the definition of <EM>scalar representation
  for \f$ x \f$</EM>, where \f$ x \f$ is any symbol.
  This is more briefly written <EM>scalar</EM> and defined inductively
  as follows:
  - every number is a scalar;
  - every symbolic constant is a scalar;
  - every parameter different from \f$ x \f$ is a scalar;
  - if \f$ f \f$ is any unary function and \f$ x \f$ is a
    scalar representation, then \f$ f(x) \f$ is a scalar;
  - if \f$ a \f$ and \f$ b \f$ are scalars then
    \f$ a+b \f$, \f$ a*b \f$, and \f$ a^b \f$ are scalars.
*/
bool
is_scalar_representation(const Expr& e, const Symbol& x) {
  if (e.is_a_number())
    return true;
  else if (e.is_a_constant())
    return true;
  else if (e.is_a_symbol() && e != x)
    return true;
  else if (e.is_a_power())
    return is_scalar_representation(e.op(0), x)
      && is_scalar_representation(e.op(1), x);
  else if (e.is_a_function())
    return is_scalar_representation(e.op(0), x);
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (!is_scalar_representation(e.op(i), x))
	return false;
    return true;
  }
  return false;
}

//! \brief
//! Returns <CODE>true</CODE> if \p e is a polynomial in \p x;
//! returns <CODE>false</CODE> otherwise.
/*!
  This function realizes the definition of <EM>polynomial in \f$ x \f$</EM>,
  where \f$ x \f$ is any symbol.
  This is more briefly written <EM>polynomial</EM> and defined inductively
  as follows:
  - every scalar representation for \f$ x \f$ is a polynomial;
  - \f$ x \f$ is a polynomial;
  - if \f$ a \f$ is a polynomial in \f$ x \f$ and \f$ b \f$ is a positive
    integer, then \f$ a^b \f$ is a polynomial;
  - if \f$ a \f$ and \f$ b \f$ are polynomials then
    \f$ a + b \f$ and \f$ a * b \f$ are polynomials.
*/
static bool
is_polynomial(const Expr& e, const Symbol& x) {
  if (is_scalar_representation(e, x))
    return true;
  else if (e == x)
    return true;
  else if (e.is_a_power()) {
    if (is_polynomial(e.op(0), x))
      if (e.op(1).is_a_number()) {
	Number exponent = e.op(1).ex_to_number();
	if (exponent.is_positive_integer())
	  return true;
      }
  }
  else if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (!is_polynomial(e.op(i), x))
	return false;
    return true;
  }
  return false;
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
isolate_polynomial_part(const Expr& e, const Symbol& x,
			Expr& polynomial, Expr& rest) {
  if (e.is_a_add()) {
    polynomial = 0;
    rest = 0;
    for (unsigned i = e.nops(); i-- > 0; ) {
      if (is_polynomial(e.op(i), x))
	polynomial += e.op(i);
      else
	rest += e.op(i);
    }
  }
  else if (is_polynomial(e, x)) {
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
Expr
convert_to_integer_polynomial(const Expr& p, const Symbol& x) {
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
Expr
convert_to_integer_polynomial(const Expr& p, const Symbol& x,
                              Number& factor) {
  assert(p.is_rational_polynomial());
  unsigned deg_p = p.degree(x);

  // Choose non-zero starting value and compute least common
  // multiple of denominators.
  Number t_lcm = p.coeff(x, deg_p).denominator().ex_to_number();
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
  - \f$ R(g, f) = (-1)^{\deg(f)\deg(g)} R(f, g) \f$; 
  - \f$ R(f, g) = a^{\deg(g) - \deg(r)} R(f, r) \f$
    if \f$ g = fq + r \f$ for some polynomials \f$ q \f$  and \f$ r \f$.
    Here \f$ a \f$ is the leading coefficient of the polynomial \f$ f \f$.
  - \f$ R(f, b) = b^{\deg(f)} \f$ if \f$ b \f$ is a scalar.
*/
Expr
resultant(const Expr& p, const Expr& q, const Symbol& x) {
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

} // namespace Parma_Recurrence_Relation_Solver
