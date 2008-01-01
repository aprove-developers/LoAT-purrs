/* Declaration of some utility functions.
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

#ifndef PURRS_util_hh
#define PURRS_util_hh 1

#include "Number.defs.hh"
#include "Symbol.types.hh"
#include "Expr.defs.hh"
#include "Expr_List.types.hh"
#include <vector>
#include <iostream>

namespace Parma_Recurrence_Relation_Solver {

//! \brief
//! Returns <CODE>true</CODE> if at least one element of the vector \p v
//! is different from \f$ 0 \f$.
//! Returns <CODE>false</CODE> otherwise.
bool
vector_not_all_zero(const std::vector<Expr>& v);

//! \brief
//! Computes the gcd, greatest common divisor, between the integers
//! \p n and \p m.
int
gcd(int n, int m);

//! \brief
//! Computes the gcd, greatest common divisor, between polynomials
//! with possibly rational coefficients.
Expr
general_gcd(const Expr& p, const Expr& q, const Symbol& x);

//! \brief
//! If all the elements of \p v are integers returns the lcm,
//! least common multiple, among them; returns the product of
//! them otherwise.
Number
lcm(const std::vector<Number>& v);

//! Construct a partial factorization of the integer \p n.
void 
partial_factor(const Number& n,
	       std::vector<Number>& bases, std::vector<int>& exponents);

//! Finds divisors of positive integer \p n.
bool
find_divisors(Number n, std::vector<Number>& divisors);

//! \brief
//! Returns bases and exponents of each factor of \p e in the pair of vectors
//! \p bases and \p exponents.
void
split_bases_exponents(const Expr& e,
		      std::vector<Expr>& bases, std::vector<Expr>& exponents);

//! Computes the cubic root of \p e.
Expr
cubic_root(const Expr& e);

void
isolate_polynomial_part(const Expr& p, const Symbol& var,
			Expr& poly, Expr& no_poly);

//! Finds associate primitive polynomial with integer coefficients.
Expr
convert_to_integer_polynomial(const Expr& p, const Symbol& x);

//! Finds associate primitive polynomial with integer coefficients.
Expr
convert_to_integer_polynomial(const Expr& p, const Symbol& x,
			      Expr& factor);

//! \brief
//! Returns the resultant between the polynomials \p p and \p q
//! in the variable \p x using the Sylvester matrix method.
Expr
sylvester_matrix_resultant(const Expr& p, const Expr& q, const Symbol& x);

//! \brief
//! Returns the resultant between the polynomials \p p and \p q
//! in the variable \p x, using Euclid's algorithm or,
//! when this last fails, using the Sylvester matrix method.
Expr
resultant(const Expr& p, const Expr& q, const Symbol& x);

bool
find_domain_in_N(const Expr& e, const Symbol& x, Number& z);

bool
has_parameters(const Expr& e);

//! Kinds of substitution of the symbolic initial condition \f$ x() \f$.
enum Mode_Subs {
  //! Substitute the function \f$ x() \f$ with \f$ k^{x()} \f$.
  EXPONENT,

  //! Substitute the function \f$ x() \f$ with \f$ log_k x() \f$.
  ARGUMENT_LOG,
};

//! \brief
//! Substitutes all the symbolic initial conditions occurring in \p e
//! with another expression depending from \p mode_subs_x.
/*!
  If \p mode_subs_x is equal to <CODE>EXPONENT</CODE> then this function
  substitutes every occurrence in \p e of the function \f$ x() \f$ with
  \f$ k^{x()} \f$.
  If \p mode_subs_x is equal to <CODE>ARGUMENT_LOG</CODE> then this function
  substitutes every occurrence in \p e of the function \f$ x() \f$ with
  \f$ \log_k x() \f$.
*/
Expr
substitute_x_function(const Expr& e, const Expr& k,
		      const Mode_Subs mode_subs_x);

//! \brief
//! Returns <CODE>true</CODE> if and only if the functions \f$ x() \f$
//! occurring in \p *this are all symbolic initial conditions; returns
//! <CODE>false</CODE> otherwise, i.e., if there is at least one
//! function \f$ x() \f$ that is not a symbolic initial condition.
/*!
  A function \f$ x() \f$ is a symbolic initial condition in the
  following cases:
  - the argument is a positive integer;
  - the argument is parametric;
  - the argument is equal to \f$ mod(n, k) + h \f$, with
    \f$ k, h \in \Nset \f$;
  - the argument is equal to \f$ Sc(n, b) \f$, whit \f$ b \in \Zset_+ \f$.
*/
bool has_only_symbolic_initial_conditions(const Expr& e);

//! \brief
//! Returns <CODE>true</CODE> if and only there is at least
//! a symbolic initial condition in \p *this; returns
//! <CODE>false</CODE> otherwise.
/*!
  A function \f$ x() \f$ is a symbolic initial condition in the
  following cases:
  - the argument is a positive integer;
  - the argument is parametric;
  - the argument is equal to \f$ mod(n, k) + h \f$, with
    \f$ k, h \in \Nset \f$;
  - the argument is equal to \f$ Sc(n, b) \f$, whit \f$ b \in \Zset_+ \f$.
*/
bool has_at_least_a_symbolic_initial_condition(const Expr& e);

#define DD_MSG(s) std::cout << s << std::endl
#define DD_VAR(x) std::cout << #x " = " << x << std::endl
#define DD_MSGVAR(s, x) std::cout << s << #x " = " << x << std::endl
#define DD_VEC(vec, first, last) \
do { \
  for (int i = (int) first; i <= (int) last; ++i) \
    std::cout << #vec << "[" << i << "] = " << vec[i] << std::endl; \
} while (0)

#if defined(NOISY) && NOISY
#define D_MSG(s) DD_MSG(s)
#define D_VAR(x) DD_VAR(x) 
#define D_MSGVAR(s, x) DD_MSGVAR(s, x) 
#define D_VEC(vec, first, last) DD_VEC(vec, first, last)
#else
#define D_MSG(s)
#define D_VAR(x)
#define D_MSGVAR(s, x)
#define D_VEC(vec, first, last)
#endif

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_util_hh)
