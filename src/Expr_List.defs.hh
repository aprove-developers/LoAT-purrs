/* *****************
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

#ifndef PURRS_Expr_List_defs_hh
#define PURRS_Expr_List_defs_hh 1

#include "Expr_List.types.hh"
#include "Expr.types.hh"
#include "Symbol.types.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

class Parma_Recurrence_Relation_Solver::Expr_List {
public:
  //! Ordinary copy-constructor.
  Expr_List();

  explicit Expr_List(const Symbol& symb);
  explicit Expr_List(const Expr& exp1, const Expr& exp2);
  explicit Expr_List(const Expr& exp1, const Expr& exp2, const Expr& exp3,
		     const Expr& exp4, const Expr& exp5);

  //! Copy-constructor.
  Expr_List(const Expr_List& lst);

  //! Destructor.
  ~Expr_List();

  //! Assignment operator.
  Expr_List& operator=(const Expr_List& lst);

  unsigned nops() const;
  Expr op(unsigned i) const;

  Expr_List append(const Expr& exp);
  Expr_List prepend(const Expr& exp);
  // FIXME: dovrebbe tornare Expr_List&?
  Expr_List remove_first();

private:
  GiNaC::lst l;

  friend class Expr;
  friend class Matrix;

public:
  Expr_List(const GiNaC::lst& gl);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Expr_List.inlines.hh"

#endif // !defined(PURRS_Expr_List_defs_hh)
