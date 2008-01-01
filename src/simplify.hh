/* To be written.
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

#ifndef PURRS_simplify_hh
#define PURRS_simplify_hh 1

#include "Expr.types.hh"

namespace Parma_Recurrence_Relation_Solver {

//! Kinds of sum's simplification.
enum Sum_Simplification_Kind {
  //! \brief
  //! Rewrite the sum with the upper limit of the form \f$ n + j \f$
  //! (\f$ n \f$ is a symbol and \f$ j \in \Zset \f$), so that the
  //! upper limit is \f$ n \f$.
  REWRITE_UPPER_LIMIT,
  
  //! \brief
  //! Split the sum in many sums how many are the addends of the summand
  //! and compute, when possible, symbolic sums.
  COMPUTE_SUM
};

Expr
simplify_ex_for_input(const Expr& e, bool input);

Expr
simplify_ex_for_output(const Expr& e, bool input);

Expr
simplify_numer_denom(const Expr& e);

Expr
simplify_binomials_factorials_exponentials(const Expr& e);

Expr
simplify_sum(const Expr& e, const Sum_Simplification_Kind simplification);

Expr
simplify_collect_sums(const Expr& e);

Expr
simplify_logarithm(const Expr& e);

Expr
simplify_all(const Expr& e);

} // namespace Parma_Recurrence_Relation_Solver

#endif
