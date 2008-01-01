/* Init class implementation: inline functions.
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

#ifndef PURRS_Init_inlines_hh
#define PURRS_Init_inlines_hh 1

#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Init::Init() {
  // Only when the first Init object is constructed...
  if (count++ == 0) {
    // ... create and immediately destroy an instance
    // of each user-defined function.
    Symbol a;
    Expr* dummy;
    dummy = new Expr(x(a));
    delete dummy;
    dummy = new Expr(sum(a, 1, 2, a));
    delete dummy;
    dummy = new Expr(prod(a, 1, 2, a));
    delete dummy;
  }
}

inline
Init::~Init() {
  // Only when the last Init object is destroyed...
  if (--count == 0)
    // ... do nothing for the time being.
    ;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Init_inlines_hh)
