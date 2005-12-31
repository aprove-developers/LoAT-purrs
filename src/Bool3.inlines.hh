/* Bool3 class implementation: inline functions.
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

#ifndef PURRS_Bool3_inlines_hh
#define PURRS_Bool3_inlines_hh 1

namespace Parma_Recurrence_Relation_Solver {

inline Bool3
lor(Bool3 a, Bool3 b) {
  return
    (a == ALWAYS || b == ALWAYS) ? ALWAYS
    : (a == NEVER && b == NEVER) ? NEVER
    : MAYBE;
}

inline Bool3
operator||(Bool3 a, Bool3 b) {
  return lor(a, b);
}

inline Bool3
land(Bool3 a, Bool3 b) {
  return
    (a == NEVER || b == NEVER) ? NEVER
    : (a == ALWAYS && b == ALWAYS) ? ALWAYS
    : MAYBE;
}

inline Bool3
operator&&(Bool3 a, Bool3 b) {
  return land(a, b);
}

inline Bool3
lnot(Bool3 a) {
  return
    a == NEVER ? ALWAYS
    : a == ALWAYS ? NEVER
    : MAYBE;
}

inline Bool3
operator!(Bool3 a) {
  return lnot(a);
}

inline bool
is_strong(Bool3 a) {
  return a == ALWAYS || a == NEVER;
}

inline bool
stronger(Bool3 a, Bool3 b) {
  return b == MAYBE && (a == ALWAYS || a == NEVER);
}

inline bool
stronger_eq(Bool3 a, Bool3 b) {
  return a == b || (b == MAYBE && (a == ALWAYS || a == NEVER));
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Bool3_inlines_hh)
