/* Bool3 class declaration.
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

#ifndef PURRS_Bool3_defs_hh
#define PURRS_Bool3_defs_hh 1

#include "Bool3.types.hh"
#include <iosfwd>

namespace Parma_Recurrence_Relation_Solver {

//! Logical or.
Bool3 lor(Bool3 a, Bool3 b);

//! Logical and.
Bool3 land(Bool3 a, Bool3 b);

//! Logical not.
Bool3 lnot(Bool3 a);

//! Logical and operator.
Bool3 operator&&(Bool3 a, Bool3 b);

//! Logical or operator.
Bool3 operator||(Bool3 a, Bool3 b);

//! Logical not operator.
Bool3 operator!(Bool3 a);

//! True if and only if \p a is a strong assertion.
bool is_strong(Bool3 a);

//! True if and only if \p a is stronger than \p b.
bool stronger(Bool3 a, Bool3 b);

//! True if and only if \p a is stronger than or equal to \p b.
bool stronger_eq(Bool3 a, Bool3 b);

//! Output operator.
std::ostream& operator<<(std::ostream& s, Bool3 a);

} // namespace Parma_Recurrence_Relation_Solver

#include "Bool3.inlines.hh"

#endif // !defined(PURRS_Bool3_defs_hh)
