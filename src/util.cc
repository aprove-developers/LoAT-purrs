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

#include "numerator_denominator.hh"
#include "Expr.defs.hh"
#include "Number.defs.hh"
#include "Symbol.defs.hh"
#include "Recurrence.defs.hh"
#include "util.hh"
#include <algorithm>
#include <iterator>

namespace PURRS = Parma_Recurrence_Relation_Solver;

static const unsigned
FACTOR_THRESHOLD = 100;

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

/*!
  Construct a partial factorization of the integer \p n.
  \p n is tested for divisibility by 2 and by odd integers between 3
  and <CODE>FACTOR_THRESHOLD</CODE>.
  The partially factored form is returned in the pair of vectors 
  \p bases and \p exponents, of <CODE>Number</CODE>s
  and <CODE>int</CODE>s respectively.
*/
void 
PURRS::partial_factor(const Number& n,
		      std::vector<Number>& bases,
		      std::vector<int>& exponents) {
  assert(n.is_integer());
  Number m = abs(n);
  assert(m != 0);
  int k = 0;
  while (mod(m, 2) == 0) { // the case 2 is handled separately 
    m /= 2;
    ++k;
  }
  if (k > 0) {
    bases.push_back(2);
    exponents.push_back(k);
  }
  for (unsigned i = 3; (i < FACTOR_THRESHOLD) && (i * i <= m); i += 2) {
    k = 0;
    while (mod(m, i) == 0) { // test for divisibility by the odd integer i
      m /= i;
      ++k;
    }
    if (k > 0) {
      bases.push_back(i);
      exponents.push_back(k);
    }
  }
  if (m > 1) { // here n has not necessarily been factored completely 
    bases.push_back(m);
    exponents.push_back(1);
  }
}

namespace {
using namespace PURRS;

void
split_bases_exponents_factor(const Expr& e,
			     std::vector<Expr>& bases,
			     std::vector<Expr>& exponents) {
  D_MSGVAR("---", e);
  Number e_num;
  if (e.is_a_number(e_num) && e_num.is_integer()) {
    std::vector<Number> e_num_bases;
    std::vector<int> e_num_exponents;
    partial_factor(e_num, e_num_bases, e_num_exponents);
    D_VEC(e_num_bases, 0, e_num_bases.size()-1);
    D_VEC(e_num_exponents, 0, e_num_exponents.size()-1);
    // If `e_num' is a negative integer the sign minus is ignored from
    // `partial_factor()'.
    if (!e_num.is_nonnegative_integer()) {
      bases.push_back(-1);
      exponents.push_back(1);
    }
    // `e_num == 1' is ignored from `partial_factor()'.
    if (e == 1) {
      bases.push_back(1);
      exponents.push_back(1);
    }
    copy(e_num_bases.begin(), e_num_bases.end(),
	 inserter(bases, bases.begin()));
    copy(e_num_exponents.begin(), e_num_exponents.end(),
	 inserter(exponents, exponents.begin()));
  }
  else
    if (e.is_a_power()) {
      bases.push_back(e.arg(0));
      exponents.push_back(e.arg(1));
    }
    else {
      bases.push_back(e);
      exponents.push_back(1);
    }
  D_VEC(bases, 0, bases.size()-1);
  D_VEC(exponents, 0, exponents.size()-1);
}

} // anonymous namespace

/*!
  Given an expression \f$ e \f$, this function returns bases and exponents of
  each factor of \f$ e \f$ in a pair of vectors.
*/
void
PURRS::split_bases_exponents(const Expr& e, 
			     std::vector<Expr>& bases,
			     std::vector<Expr>& exponents) {
  if (e.is_a_mul())
    for (unsigned i = e.nops(); i-- > 0; )
      split_bases_exponents_factor(e.op(i), bases, exponents);
  else
    split_bases_exponents_factor(e, bases, exponents);
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

namespace {
using namespace PURRS;

/*!
  Returns <CODE>true</CODE> if \p e is on the form
  \f$ a n + b \f$ with \f$ a \in \Nset \setminus \{0\} \f$ and
  \f$ b \in \Zset \f$;
  returns <CODE>false</CODE> otherwise. 
*/
bool
ok_argument_factorial(const Expr& argument, const Symbol& n) {
  D_VAR(argument);
  if (argument.is_a_add() && argument.nops() == 2) {
    const Expr& first = argument.op(0);
    const Expr& second = argument.op(1);
    Number a;
    Number b;
    if ((first.is_a_mul() && first.nops() == 2
	 && (first.op(0) == n && first.op(1).is_a_number(a)
	     || first.op(1) == n && first.op(0).is_a_number(a))
	 && second.is_a_number(b))
	||
	(second.is_a_mul() && second.nops() == 2
	 && (second.op(0) == n && second.op(1).is_a_number(a)
	     || second.op(1) == n && second.op(0).is_a_number(a))
	 && first.is_a_number(b)))
      if (a.is_positive_integer() && b.is_integer())
	return true;
    if ((first == n && second.is_a_number(b))
	|| (second == n && first.is_a_number(b)))
      if (b.is_integer())
	return true;
  }
  return false;
}

} // anonymous namespace

/*!
  Let \f$ e(n) \f$ be the expression in \p n contained in \p e,
  which is assumed to be already expanded and with denominator equal to 1.
  This function find the largest positive integer that cancel \f$ e(n) \f$
  and, if it is bigger than \p z, store it in \p z; if do not exist a
  positive integer that cancel \f$ e(n) \f$ or exist but smaller than \p z,
  then \p z is left unchanged.
*/
void
PURRS::largest_positive_int_zero(const Expr& e, Number& z) {
  assert(denominator(e) == 1);
  if (e.is_polynomial(Recurrence::n)) {
    Expr partial_e = e;
    unsigned lower_degree = partial_e.ldegree(Recurrence::n);
    while (lower_degree > 0) {
      partial_e = quo(partial_e, Recurrence::n, Recurrence::n);
      lower_degree = partial_e.ldegree(Recurrence::n);
      if (z < 0)
	z = 0;
    }
    std::vector<Number> potential_roots;
    Number constant_term
      = abs(partial_e.tcoeff(Recurrence::n).ex_to_number());
    // Find the divisors of the constant term.
    if (constant_term.is_positive_integer())
      find_divisors(constant_term, potential_roots);
    // Find non-negative integral roots of the denominator.
    for(unsigned i = potential_roots.size(); i-- > 0; ) {
      Number temp = partial_e.substitute(Recurrence::n,
					 potential_roots[i]).ex_to_number();
      if (temp == 0 &&  potential_roots[i] > z)
	z = potential_roots[i];
    }
  }
  else {
    if (e.is_a_add() || e.is_a_mul())
      for (unsigned i = e.nops(); i-- > 0; )
	largest_positive_int_zero(e.op(i), z);
    else if (e.is_a_function()) {
      if (e.is_the_log_function())
	if (e.arg(0).is_polynomial(Recurrence::n))
	  largest_positive_int_zero(e.arg(0) - 1, z);
	else if (e.arg(0).is_the_log_function()) {
	  largest_positive_int_zero(e.arg(0), z);
	  ++z;
	}  
      else if (e.is_the_factorial_function())
	// If it is in the form `(a n + b)!' with `a' positive integer
	// then we find the minimum `n' such that `a n + b >= 0'.
	if (ok_argument_factorial(e.arg(0), Recurrence::n))
	  largest_positive_int_zero(e.arg(0), z);
    }
  }
}

//! Returns <CODE>true</CODE> if \p e contains parameters;
//! returns <CODE>false</CODE> otherwise.
/*!
  The parameters are all symbols different from \p n and the initial
  conditions \f$ x(k) \f$ with \f$ k \f$ a positive integer.
  Note: \p e does not contain \f$ x(f) \f$ with \f$ f \f$ an expression
  containig \p n.
*/
bool
PURRS::find_parameters(const Expr& e) {
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (find_parameters(e.op(i)))
	return true;
  }
  else if (e.is_a_power()) {
    if (find_parameters(e.arg(0)) || find_parameters(e.arg(1)))
      return true;
  }
  else if (e.is_a_function()) {
    // In this case the function `x' is surely an initial condition.
    if (e.is_the_x_function())
      return true;
    else
      for (unsigned i = e.nops(); i-- > 0; )
	if (find_parameters(e.arg(i)))
	  return true;
  }
  else
    if (e.is_a_symbol() && e != Recurrence::n)
      return true;
  return false;
}
