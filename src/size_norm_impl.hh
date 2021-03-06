/* Inline functions implementing the computation of size_norms.
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

#ifndef PURRS_size_norm_impl_hh
#define PURRS_size_norm_impl_hh 1

#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

template <typename SymbolHandler>
inline unsigned int
generic_size_norm(const Expr& e, const SymbolHandler& sh) {
  int count = 1;
  if (e.is_a_add() || e.is_a_mul())
    for (unsigned int i = e.nops(); i-- > 0; )
      count += generic_size_norm(e.op(i), sh);
  else if (e.is_a_power())
    count += generic_size_norm(e.arg(0), sh) + generic_size_norm(e.arg(1), sh);
  else if (e.is_a_function())
    for (unsigned int i = e.nops(); i-- > 0; )
      count += generic_size_norm(e.arg(i), sh);
  else if (e.is_a_symbol())
    count = sh.size_norm(e.ex_to_symbol());
  // Leave count to 1 in case e is a number or a constant.
  return count;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_size_norm_impl_hh)
