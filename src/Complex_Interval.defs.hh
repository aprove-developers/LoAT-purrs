/* Complex_Interval class declaration.
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

#ifndef PURRS_Complex_Interval_defs_hh
#define PURRS_Complex_Interval_defs_hh 1

#include "Complex_Interval.types.hh"
#include "Expr.types.hh"
#include "CInterval.types.hh"
#include "complint.hh"

#include <ginac/ginac.h>
#include <string>

namespace Parma_Recurrence_Relation_Solver {

class Complex_Interval {
public:
  //! Builds a symbol with a unique name.
  Complex_Interval();

  //! Builds a symbol named \p n.
  Complex_Interval(const CInterval& y);

  //! Copy-constructor.
  Complex_Interval(const Complex_Interval& y);

  //! Destructor.
  ~Complex_Interval();

  //! Assignment operator.
  Complex_Interval& operator=(const Complex_Interval& y);

  const CInterval& get_interval() const;

private:
  friend class Expr;

  GiNaC::complint i;

  //! Builds the symbol corresponding to \p gi.
  Complex_Interval(const GiNaC::complint& y);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Complex_Interval.inlines.hh"

#endif // !defined(PURRS_Complex_Interval_defs_hh)
