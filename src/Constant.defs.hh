/* Constant class declaration.
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

#ifndef PURRS_Constant_defs_hh
#define PURRS_Constant_defs_hh 1

#include "Constant.types.hh"
#include "Expr.types.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

//! The symbolic constants (e.g. \f$ \pi \f$).
/*!
  ...
*/
class Constant {
public:
  //! Default constructor.
  Constant();

  //! Copy-constructor.
  Constant(const Constant& k);

  //! Destructor.
  ~Constant();

  //! Assignment operator.
  Constant& operator=(const Constant& k);

  //! The Archimedes'constant \f$ \pi = 3.14159\dots \f$.
  static const Constant Pi;

private:
  friend class Expr;

  friend bool operator==(const Expr& e, const Constant& c);

  GiNaC::constant c;

  //! Builds the constant corresponding to \p gc.
  Constant(const GiNaC::constant& gc);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Constant.inlines.hh"

#endif // !defined(PURRS_Constant_defs_hh)
