/* Definition of some utility functions.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#ifndef NOISY
#define NOISY 0
#endif

#include "util.hh"
#include "factorize.hh"
#include "alg_eq_solver.hh"
#include "Matrix.defs.hh"
#include "Expr.defs.hh"
#include "Number.defs.hh"
#include "Symbol.defs.hh"
#include "Recurrence.defs.hh"
#include <algorithm>
#include <iterator>

namespace PURRS = Parma_Recurrence_Relation_Solver;

#define Napier exp(Expr(1))

static const unsigned int
FACTOR_THRESHOLD = 100;

static const unsigned int FIND_DIVISORS_MAX = 11;

static const unsigned int
FIND_DIVISORS_THRESHOLD = FIND_DIVISORS_MAX*FIND_DIVISORS_MAX;

bool
PURRS::vector_not_all_zero(const std::vector<Expr>& v) {
  for (unsigned int i = v.size(); i-- > 0; )
    if (!v[i].is_zero())
      return true;
  return false;
}

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

PURRS::Number
PURRS::lcm(const std::vector<Number>& v) {
  Number n = 1;
  for (unsigned int i = v.size(); i-- > 0; )
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
  for (unsigned int i = 3; (i < FACTOR_THRESHOLD) && (i * i <= m); i += 2) {
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

//! Finds all the positive divisors of the strictly positive integer \p n
//! if it is less than <CODE>FIND_DIVISORS_THRESHOLD</CODE>.
/*!
  This routine inserts into \p divisors all the positive divisors of
  the strictly positive integer \p n if it is less than
  <CODE>FIND_DIVISORS_THRESHOLD</CODE>.
  Returns <CODE>false</CODE> if \p n is bigger than
  <CODE>FIND_DIVISORS_THRESHOLD</CODE> and, in this case, the function not
  find the divisors of \p n; returns <CODE>true</CODE> otherwise.
*/
bool
PURRS::find_divisors(Number n, std::vector<Number>& divisors) {
  assert(n.is_positive_integer());
  if (n < FIND_DIVISORS_THRESHOLD) {
    unsigned int m = n.to_unsigned_int();
    // Once a divisor `i' is found, it is pushed onto the vector `divisors'
    // along with its conjugate `j = n/i', provided that `j' is less than `i'.
    if (m == 1)
      divisors.push_back(1);
    else
      for (unsigned int i = 1, j = m; i < FIND_DIVISORS_MAX && i < j; ++i) {
	j = m / i;
	unsigned int r = m % i;
	if (r == 0) {
	  divisors.push_back(i);
	  if (i < j)
	    divisors.push_back(j);
	}
      }
    return true;
  }
  return false;
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
    for (unsigned int i = e.nops(); i-- > 0; )
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
    for (unsigned int i = e.nops(); i-- > 0; ) {
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
  assert(p.is_rational_polynomial(x));
  unsigned int deg_p = p.degree(x);

  // Choose non-zero starting value and compute least common
  // multiple of denominators.
  Expr t_lcm = p.coeff(x, deg_p).denominator();
  for (unsigned int i = 0; i <= deg_p; ++i)
    t_lcm = lcm(t_lcm, p.coeff(x, i).denominator());
  return (p * t_lcm).expand().primpart(x);
}

/*!
  Converts a polynomial with rational coefficients into the 
  associate primitive polynomial with integer coefficients.
  This version also returns an expression containing the factor
  used to convert.
*/
PURRS::Expr
PURRS::convert_to_integer_polynomial(const Expr& p, const Symbol& x,
				     Expr& factor) {
  assert(p.is_rational_polynomial(x));
  unsigned int deg_p = p.degree(x);

  // Choose non-zero starting value and compute least common
  // multiple of denominators.
  Expr t_lcm = p.coeff(x, deg_p).denominator();
  for (unsigned int i = 0; i <= deg_p; ++i)
    t_lcm = lcm(t_lcm, p.coeff(x, i).denominator());

  const Expr& q = (p * t_lcm).expand().primpart(x);
  factor = p.lcoeff(x) * pwr(q.lcoeff(x), -1);
  return q;
}

//! \brief
//! Returns the resultant between the polynomials \p p and \p q
//! in the variable \p x using the Sylvester matrix method.
/*!
  Let \f$ m \f$ and \f$ n \f$ be the degree of the polynomials \p p
  and \p q, respectively.
  The resultant between \p p and \p q is equal to the determinant
  of the Sylvester matrix.
  The Sylvester matrix is a \f$ m + n \times m + n \f$ matrix built
  like follow: in the first \f$ n \f$ rows there are the coefficients
  of \p p and in the remaining \f$ m \f$ rows there are the coefficients
  of \p q.
  The coefficients of \p p occurs in the first \f$ m + 1 \f$ positions
  of the first row, while in the remaining positions there are all zeros;
  in the second row the coefficients are shifted forward by one position,
  so will be one zero in the first position, and so forth. In the
  \f$ n \f$-th row the coefficients are in the last \f$ m + 1 \f$ positions,
  while in the first positions there are all zeroes.
  The coefficients of \p q in the last \f$ m \f$ rows are arranged in the
  analogous way.
*/
PURRS::Expr
PURRS::sylvester_matrix_resultant(const Expr& p, const Expr& q,
				  const Symbol& x) {
  assert(p.is_polynomial(x) && q.is_polynomial(x));
  unsigned int deg_p = p.degree(x);
  unsigned int deg_q = q.degree(x);
  // The list `elements_sylvester_matrix' will contain all the elements
  // of the matrix `sylvester' in succession starting from the first row.
  Expr_List elements_sylvester_matrix;
  for (unsigned int i = 1; i <= deg_q; ++i) {
    for (unsigned int j = 1; j < i; ++j)
      elements_sylvester_matrix.append(0);
    for (unsigned int j = i + deg_p; j >= i; --j)
      elements_sylvester_matrix.append(p.coeff(x, j - i));
    for (unsigned int j = i + deg_p + 1; j <= deg_p + deg_q; ++j)
      elements_sylvester_matrix.append(0);
  }
  for (unsigned int i = 1; i <= deg_p; ++i) {
    for (unsigned int j = 1; j < i; ++j)
      elements_sylvester_matrix.append(0);
    for (unsigned int j = i + deg_q; j >= i; --j)
      elements_sylvester_matrix.append(q.coeff(x, j - i));
    for (unsigned int j = i + deg_q + 1; j <= deg_p + deg_q; ++j)
      elements_sylvester_matrix.append(0);
  }
  Matrix sylvester(deg_p + deg_q, deg_p + deg_q, elements_sylvester_matrix);
  return sylvester.determinant();
}

//! \brief
//! Computes the resultant of the polynomials \p f and \p g in the
//! variable \p x with rational coefficients, using Euclid's algorithm or,
//! when this last fails, using the Sylvester matrix method.
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
  assert(p.is_rational_polynomial(x));
  assert(q.is_rational_polynomial(x));
  Expr f = p.expand();
  Expr g = q.expand();
  Expr res = 1;
  unsigned int deg_f = f.degree(x);
  unsigned int deg_g = g.degree(x);

  // Special case: `f' or `g' is a constant polynomial. By definition
  // `Res(f, g) = f.lcoeff(n)^g.degree(n) * g.lcoeff(n)^f.degree(n)'. 
  if (deg_f == 0 || deg_g == 0)
    res = pwr(f.lcoeff(x), deg_g) * pwr(g.lcoeff(x), deg_f);
  else {
    // Modified Euclid's algorithm starts here.
    while (deg_f > 0) {
      // `prem()' computes the pseudo-remainder of `g' and `f' which satisfies
      // `factor * g = factor * f * q + prem(g, f, x)' where `q' is the
      // quotient of `g' and `f' and
      // `factor = f.lcoeff(x)^(g.degree(x) - f.degree(x) + 1)'.
      // `prem()' wants only rational polynomials. The expressions
      // `f' and `g' are surely polynomials but in this point they could be
      // not enough simplified so that the system could not recognize them.
      if (!f.is_rational_polynomial(x)) {
	Expr common_factor;
	Expr rem;
	factorize(f, common_factor, rem);
	f = (common_factor * rem).expand();
      }
      if (!g.is_rational_polynomial(x)) {
	Expr common_factor;
	Expr rem;
	factorize(g, common_factor, rem);
	g = (common_factor * rem).expand();
      }
      if (!f.is_rational_polynomial(x) || !g.is_rational_polynomial(x))
	// The last chanche to compute the resultant is to use the
	// Sylvester matrix method.
	return sylvester_matrix_resultant(p.expand(), q.expand(), x);
      Expr r = prem(g, f, x);
      Expr factor = pwr(f.lcoeff(x), deg_g - deg_f + 1);
      // The rest of euclidean's division is given by the ratio
      // `pseudo-remainder / factor'.
      r *= pwr(factor, -1);
      unsigned int deg_r = r.expand().degree(x);
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
  Let \f$ e(x) \f$ be the expression in \p x contained in \p e,
  which is assumed to be already expanded and with denominator equal to 1.
  This function returns <CODE>true</CODE> if \f$ e(x) \f$ is well-defined:
  polynomials are always well-defined while if \f$ e(x) \f$ is not a
  polynomial tries the positive integer starting from which is well-defined
  and, if it is bigger than the actual \p z, stored it in \p z.
  If the function is not able to find the positive integer looked for then
  returns <CODE>false</CODE>.
  Moreover, tries an integer that cancel \f$ e(x) \f$ and if it is bigger
  than \p z, stores it in \p z, otherwise \p z is left unchanged.
*/
bool
find_domain_in_N_on_expanded_ex(const Expr& e, const Symbol& x, Number& z) {
  bool ok = false;
  if (e.has(x)) {
    // `e' is a polynomial in `x'.
    if (e.is_polynomial(x)) {
      // Polynomials are always well-defined.
      ok = true;
      if (e.is_rational_polynomial(x)) {
	Expr e_tmp = convert_to_integer_polynomial(e, x);
	std::vector<Polynomial_Root> roots;
	bool all_distinct;
	if (find_roots(e_tmp, x, roots, all_distinct))
	  for (unsigned int i = roots.size(); i-- > 0; ) {
	    Number num_root;
	    if (roots[i].value().is_a_number(num_root) && num_root.is_real())
	      while (z < num_root)
		++z;
	  }
      }
    }
    else {
      // `e' is not a polynomial in `x'.
      if (e.is_a_add() || e.is_a_mul()) {
	for (unsigned int i = e.nops(); i-- > 0; )
	  if (!find_domain_in_N_on_expanded_ex(e.op(i), x, z))
	    return false;
	ok = true;
      }
      else if (e.is_a_function()) {
	if (e.is_the_log_function()) {
	  if (e.arg(0).is_polynomial(x)) {
	    ok = true;
	    find_domain_in_N_on_expanded_ex(e.arg(0) - 1, x, z);
	  }
	  else if (e.arg(0).is_the_log_function()) {
	    ok = true;
	    find_domain_in_N_on_expanded_ex(e.arg(0), x, z);
	    while (compare(z, Napier) == -1)
	      ++z;
	  }
	}
	else if (e.is_the_factorial_function())
	  // If it is in the form `(a n + b)!' with `a' positive integer
	  // then we find the minimum `n' such that `a n + b >= 0'.
	  if (ok_argument_factorial(e.arg(0), x))
	    if (find_domain_in_N_on_expanded_ex(e.arg(0), x, z))
	      ok = true;
      }
      else if (e.is_a_power()) {
	Number base;
	// `if' necessary for the problem of the GiNaC::polynomial:
	// ex. `2^n = e^(n log(2))' and `log(2)' is not a GiNaC::polynomial.
	if (e.arg(0).is_a_number(base) && base.is_positive()
	    && e.arg(1).has(x))
	  ok = true;
	// `p(n)^q(n) = e^(q(n)*log(p(n)))'.
	else if (find_domain_in_N(e.arg(1) * log(e.arg(0)), x, z))
	  ok = true;
      }
    }
  }
  else
    // The expression does not contain `x' and then is a polynomial in `x'.
    ok = true;
  return ok;
}

} // anonymous namespace

bool
PURRS::find_domain_in_N(const Expr& e, const Symbol& x, Number& z) {
  return find_domain_in_N_on_expanded_ex(e.numerator().expand(), x, z)
    && find_domain_in_N_on_expanded_ex(e.denominator().expand(), x, z);
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
    for (unsigned int i = e.nops(); i-- > 0; )
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
      for (unsigned int i = e.nops(); i-- > 0; )
	if (has_parameters(e.arg(i)))
	  return true;
  }
  else
    if (e.is_a_symbol() && e != Recurrence::n)
      return true;
  return false;
}

Expr
PURRS::substitute_x_function(const Expr& e, const Expr& k,
			     const Mode_Subs mode_subs_x) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += substitute_x_function(e.op(i), k, mode_subs_x);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten *= substitute_x_function(e.op(i), k, mode_subs_x);
  }
  else if (e.is_a_power())
    return pwr(substitute_x_function(e.arg(0), k, mode_subs_x),
 	       substitute_x_function(e.arg(1), k, mode_subs_x));
  else if (e.is_a_function()) {
    if (e.is_the_x_function())
      switch (mode_subs_x) {
      case EXPONENT:
	return pwr(k, e);
      case ARGUMENT_LOG:
	  return log(e) / log(k);
      default:
	throw std::runtime_error("PURRS internal error: "
				 "substitute_x_function().");
      }
    else if (e.nops() == 1)
      return PURRS::apply(e.functor(), substitute_x_function(e.arg(0), k,
						      mode_subs_x));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
	argument[i] = substitute_x_function(e.arg(i), k, mode_subs_x);
      return PURRS::apply(e.functor(), argument);
    }
  }
  else
    return e;
  return e_rewritten;
}

bool
PURRS::has_only_symbolic_initial_conditions(const Expr& e) {
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned int i = e.nops(); i-- > 0; )
      if (!has_only_symbolic_initial_conditions(e.op(i)))
	return false;
  }
  else if (e.is_a_power()) {
    if (!has_only_symbolic_initial_conditions(e.arg(0))
	|| !has_only_symbolic_initial_conditions(e.arg(1)))
      return false;
  }
  else if (e.is_a_function())
    if (e.is_the_x_function()) {
      const Expr& argument = e.arg(0);
      Number arg;
      if (!(argument.is_a_number(arg) && arg.is_nonnegative_integer())
	  && !(argument.is_the_mod_function()
	       && argument.arg(0) == Recurrence::n
	       && argument.arg(1).is_a_number(arg)
	       && arg.is_positive_integer())
	  && !(argument.is_a_add() && argument.nops() == 2
	       && (argument.op(0).is_the_mod_function()
		   || argument.op(1).is_the_mod_function()))
	  && !(argument.is_the_Sc_function()
	       && argument.arg(1).is_a_number(arg)
	       && arg.is_positive_integer()))
	return false;
    }
    else
      for (unsigned int i = e.nops(); i-- > 0; )
	if (!has_only_symbolic_initial_conditions(e.arg(i)))
	  return false;
  return true;
}

bool
PURRS::has_at_least_a_symbolic_initial_condition(const Expr& e) {
  return e.has_x_function() && has_only_symbolic_initial_conditions(e);
}
