/* Declaration of some utility functions.
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

//! Computes the gcd between \p n and \p m.
int
gcd(int n, int m);

//! Computes the gcd between polynomials with possibly rational coefficients.
Expr
general_gcd(const Expr& p, const Expr& q, const Symbol& x);

//! Computes the lcm among the integers in the vector \p v.
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
largest_positive_int_zero(const Expr& e, const Symbol& x, Number& z);

bool
has_parameters(const Expr& e);

Expr
substitute_x_function(const Expr& e, const Expr& k, bool do_power);

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
