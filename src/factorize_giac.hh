/* To be written.
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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#ifndef PURRS_factorize_giac_hh
#define PURRS_factorize_giac_hh 1

#include "Blackboard.types.hh"
#include "Expr.types.hh"
#include "Symbol.types.hh"
#include <vector>
#include <utility>

namespace Parma_Recurrence_Relation_Solver {

void
factorize_giac(const Expr& e, Expr& num_factors, Expr& not_num_factors);

void
factorize_giac(const Expr& e, const Symbol& x,
	       Expr& num_factors, Expr& not_num_factors);

void
factors_giac(const Expr& e, std::vector<std::pair<Expr, int> >& factors);

void
factors_giac(const Expr& e, const Symbol& x,
	     std::vector<std::pair<Expr, int> >& factors_e);

} // namespace Parma_Recurrence_Relation_Solver

#endif
