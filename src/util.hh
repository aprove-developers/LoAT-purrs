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

//! Computes the gcd between \f$n\f$ and \f$m\f$.
int
gcd(int n, int m);

//! Computes the gcd between polynomials with possibly rational coefficients.
Expr
general_gcd(const Expr& p, const Expr& q, const Symbol& x);

//! Computes the lcm among the integers in the vector \f$numbers\f$.
Number
lcm(const std::vector<Number>& numbers);

//! Computes the cubic root of \f$e\f$.
Expr
cubic_root(const Expr& e);

//! Removes all the elements from \p l.
void
clear(Expr_List& l);

//! Returns the expression that \p substitution binds
//! to the wildcard of index \p wild_index.
Expr
get_binding(const Expr_List& substitution, unsigned wild_index);

void
isolate_polynomial_part(const Expr& p, const Symbol& var,
			Expr& poly, Expr& no_poly);

//! Finds associate primitive polynomial with integer coefficients.
Expr
convert_to_integer_polynomial(const Expr& p, const Symbol& x);

//! Finds associate primitive polynomial with integer coefficients.
Expr
convert_to_integer_polynomial(const Expr& p, const Symbol& x,
			      Number& factor);

//! Computes the resultant between the polynomials \p p and \p q.
Expr
resultant(const Expr& p, const Expr& q, const Symbol& x);

#if defined(NOISY) && NOISY
#define D_MSG(s) std::cout << s << std::endl
#define D_VAR(x) std::cout << #x " = " << x << std::endl
#define D_MSGVAR(s, x) std::cout << s << #x " = " << x << std::endl
#define D_VEC(vec, first, last) \
do { \
  for (int i = (int) first; i <= (int) last; ++i) \
    std::cout << #vec << "[" << i << "] = " << vec[i] << std::endl; \
} while (0)
#else
#define D_MSG(s)
#define D_VAR(x)
#define D_MSGVAR(s, x)
#define D_VEC(vec, first, last)
#endif

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_util_hh)
