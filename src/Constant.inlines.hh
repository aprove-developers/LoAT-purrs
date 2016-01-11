/* Constant class implementation: inline functions.
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

#ifndef PURRS_Constant_inlines_hh
#define PURRS_Constant_inlines_hh

#include "Constant.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Constant::Constant() {
}

inline
Constant::Constant(const Constant& k)
  : c(k.c) {
};

inline Constant&
Constant::operator=(const Constant& k) {
  c = k.c;
  return *this;
};

inline
Constant::Constant(const GiNaC::constant& gc)
  : c(gc) {
}

inline
Constant::~Constant() {
}

inline
Constant Constant::Pi() {
  return GiNaC::Pi;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Constant_inlines_hh)
