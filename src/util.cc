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

#include "util.hh"
#include "factorize.hh"
#include "alg_eq_solver.hh"
#include "numerator_denominator.hh"
#include "Expr.defs.hh"
#include "Number.defs.hh"
#include "Symbol.defs.hh"
#include "Recurrence.defs.hh"
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
  if (m == 1 && bases.size() == 0) {
    bases.push_back(1);
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

Expr
sylvester_matrix_resultant(const Expr& /*p*/, const Expr& /*q*/) {
  throw
    "PURRS error: function\n"
    "`compute_resultant_with_determinant()': work in progress.\n"
    "Please come back tomorrow.";
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
      // `prem()' wants only rational polynomials. The expressions
      // `f' and `g' are surely polynomials but in this point they could be
      // not enough simplified so that the system could not recognize them.
      if (!f.is_rational_polynomial()) {
	Expr common_factor;
	Expr rem;
	factorize(f, common_factor, rem);
	f = (common_factor * rem).expand();
      }
      if (!g.is_rational_polynomial()) {
	Expr common_factor;
	Expr rem;
	factorize(g, common_factor, rem);
	g = (common_factor * rem).expand();
      }
      if (!f.is_rational_polynomial() || !g.is_rational_polynomial())
	// The last chanche to compute the resultant is to use the
	// method of the Sylvester matrix
	// (see http://mathworld.wolfram.com/SylvesterMatrix.html).
	return sylvester_matrix_resultant(f.expand(), g.expand());
      Expr r = prem(g, f, x);
      Expr factor = pwr(f.lcoeff(x), deg_g - deg_f + 1);
      // The rest of euclidean's division is given by the ratio
      // `pseudo-remainder / factor'.
      r *= pwr(factor, -1);
      unsigned deg_r = r.expand().degree(x);
      // Using rule two.
      res *= pwr(f.lcoeff(x), deg_g - deg_r);
      // Using rule one.
      if ((deg_f * deg_r) & 1 != 0)
	// `deg_f * deg_r' is odd.
	res = -res;
      g = f.expand();
      f = r.expand();
      deg_f = f.degree(x);
      deg_g = g.degree(x);
    }
    // Here `f' is a constant: use rule three.
    res *= pwr(f, deg_g);
  }
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
  // `a != 1' and `b != 0'.
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
  // `a != 1' and `b == 0'.
  else if (argument.is_a_mul() && argument.nops() == 2) {
    const Expr& first = argument.op(0);
    const Expr& second = argument.op(1);
    Number a;
    if ((first == n && second.is_a_number(a))
	|| (second == n && first.is_a_number(a)))
      if (a.is_positive_integer())
	return true;
  }
  // `a == 1' and `b == 0'.
  else if (argument == n)
    return true;
  return false;
}

/*!
  Let \f$ e(n) \f$ be the expression in \p n contained in \p e,
  which is assumed to be already expanded and with denominator equal to 1.
  This function returns <CODE>true</CODE> if it finds an integer that cancel
  \f$ e(n) \f$ or starting from which \f$ e(n) \f$ is well-defined
  (polynomials are always well-defined); if this integer is bigger than \p z,
  stores it in \p z, otherwise \p z is left unchanged.
  Returns <CODE>false</CODE> if it does not find the integer looked for.
*/
bool
largest_positive_int_zero_on_expanded_ex(const Expr& e, Number& z) {
  assert(e.is_expanded());
  bool ok = false;
  // FIXME: `if (e.has(Recurrence::n))' is necessary because here there is
  // the method `PURRS::is_polynomial()' and there are the methods
  // `GiNaC::ldegree()', `GiNaC::tcoeff()': for example `log(2)' is a
  // PURRS::polynomial but it is not a GiNaC::polynomial.
  if (e.has(Recurrence::n)) {
    // `e' is a polynomial in `n'.
    if (e.is_polynomial(Recurrence::n)) {
      // Polynomials are always well-defined.
      ok = true;
      std::vector<Polynomial_Root> roots;
      bool all_distinct;
      if (find_roots(e, Recurrence::n, roots, all_distinct))
	for (unsigned i = roots.size(); i-- > 0; ) {
	  Number num_root;
	  if (roots[i].value().is_a_number(num_root) && num_root.is_real())
	    while (z < num_root)
	      ++z;
	}
    }
    else {
      // `e' is not a polynomial in `n'.
      if (e.is_a_add() || e.is_a_mul()) {
	for (unsigned i = e.nops(); i-- > 0; )
	  if (!largest_positive_int_zero_on_expanded_ex(e.op(i), z))
	    return false;
	ok = true;
      }
      else if (e.is_a_function()) {
	if (e.is_the_log_function()) {
	  if (e.arg(0).is_polynomial(Recurrence::n)) {
	    ok = true;
	    largest_positive_int_zero_on_expanded_ex(e.arg(0), z);
	    ++z;
	  }
	  else if (e.arg(0).is_the_log_function()) {
	    ok = true;
	    largest_positive_int_zero_on_expanded_ex(e.arg(0), z);
	    // Consider an approximation of the Napier's number.
	    Number tmp = pwr(Number(2718, 1000), z);
	    while (z < tmp)
	      ++z;
	  }
	}
	else if (e.is_the_factorial_function())
	  // If it is in the form `(a n + b)!' with `a' positive integer
	  // then we find the minimum `n' such that `a n + b >= 0'.
	  if (ok_argument_factorial(e.arg(0), Recurrence::n))
	    if (largest_positive_int_zero_on_expanded_ex(e.arg(0), z))
	      ok = true;
      }
      else if (e.is_a_power()) {
	Number base;
	// `if' necessary for the problem of the GiNaC::polynomial:
	// ex. `2^n = e^(n log(2))' and `log(2)' is not a GiNaC::polynomial.
	if (e.arg(0).is_a_number(base) && base.is_positive()
	    && e.arg(1).has(Recurrence::n))
	  ok = true;
	// `p(n)^q(n) = e^(q(n)*log(p(n)))'.
	else if (largest_positive_int_zero(e.arg(1) * log(e.arg(0)), z))
	  ok = true;
      }
    }
  }
  else
    // Consider constant polynomials.
    if (e.is_polynomial(Recurrence::n))
      ok = true;
  return ok;
}

} // anonymous namespace

bool
PURRS::largest_positive_int_zero(const Expr& e, Number& z) {
  if (largest_positive_int_zero_on_expanded_ex(numerator(e).expand(), z)
      && largest_positive_int_zero_on_expanded_ex(denominator(e).expand(), z))
    return true;
  else
    return false;
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
PURRS::has_parameters(const Expr& e) {
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (has_parameters(e.op(i)))
	return true;
  }
  else if (e.is_a_power()) {
    if (has_parameters(e.arg(0)) || has_parameters(e.arg(1)))
      return true;
  }
  else if (e.is_a_function()) {
    // In this case the function `x' is surely an initial condition.
    if (e.is_the_x_function())
      return true;
    else
      for (unsigned i = e.nops(); i-- > 0; )
	if (has_parameters(e.arg(i)))
	  return true;
  }
  else
    if (e.is_a_symbol() && e != Recurrence::n)
      return true;
  return false;
}

/*!
  If \p do_power is true this function substitutes every occurrence
  of the function \f$ x() \f$ with \f$ k^{x()} \f$.
  If \p do_power is false substitutes every occurrence
  of the function \f$ x() \f$ with \f$ log_k x() \f$.
*/
Expr
PURRS::substitute_x_function(const Expr& e, const Expr& k, bool do_power) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten += substitute_x_function(e.op(i), k, do_power);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten *= substitute_x_function(e.op(i), k, do_power);
  }
  else if (e.is_a_power())
    return pwr(substitute_x_function(e.arg(0), k, do_power),
 	       substitute_x_function(e.arg(1), k, do_power));
  else if (e.is_a_function()) {
    if (e.is_the_x_function())
      if (do_power)
	return pwr(k, e);
      else
	return log(e) / log(k);
    else if (e.nops() == 1)
      return apply(e.functor(), substitute_x_function(e.arg(0), k, do_power));
    else {
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned i = 0; i < num_argument; ++i)
	argument[i] = substitute_x_function(e.arg(i), k, do_power);
      return apply(e.functor(), argument);
    }
  }
  else
    return e;
  return e_rewritten;
}
