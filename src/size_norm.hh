/* Declare the public functions for computing size norms.
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

#ifndef PURRS_approx_expr_hh
#define PURRS_approx_expr_hh 1

#include "Expr.types.hh"
#include "Matrix.types.hh"

namespace Parma_Recurrence_Relation_Solver {

//! Returns the size-norm of \p e.
unsigned int
size_norm(const Expr& e);

//! Returns the size-norm of \p ,.
unsigned int
size_norm(const Matrix& m);

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_approx_expr_hh)
