/* Declarations of global objects.
   Copyright (C) 2001-2003 Roberto Bagnara <bagnara@cs.unipr.it>

This file is part of the Parma Polyhedra Library (PURRS).

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

For the most up-to-date information see the Parma Polyhedra Library
site: http://www.cs.unipr.it/purrs/ . */

#ifndef PURRS_globals_hh
#define PURRS_globals_hh 1

namespace Parma_Recurrence_Relation_Solver {

//! \brief
//! An unsigned integral type for representing different kinds of
//! indices of the recurrence.
/*!
  Are <CODE>index_type</CODE> the following elements:
  - the order of a recurrence;
  - the rank of a functional equation;
  - the least non-negative integer \f$ j \f$ such that
    the recurrence is well-defined for \f$ n \geq j \f$.
*/
typedef unsigned int index_type;

} // namespace Parma_Recurrence_Relation_Solver

#endif
