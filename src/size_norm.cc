/* Implementation of the public functions for computing size norms.
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

#include <config.h>

#include "size_norm.hh"
#include "size_norm_impl.hh"
#include "Expr.types.hh"
#include "Symbol.types.hh"
#include "Matrix.defs.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {

struct Norm_1 {
  unsigned size_norm(const PURRS::Symbol&) const {
    return 1;
  }
};

}

unsigned
PURRS::size_norm(const Expr& e) {
  return generic_size_norm(e, Norm_1());
}

unsigned
PURRS::size_norm(const Matrix& m) {
  int count = 0;
  for (unsigned i = m.num_rows(); i-- > 0; )
    for (unsigned j = m.num_columns(); j-- > 0; )
      count += size_norm(m(i, j));
  return count;
}
