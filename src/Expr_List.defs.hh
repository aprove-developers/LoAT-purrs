/* Expr_List class declaration.
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

class Expr_List {
public:
  //! Default constructor.
  Expr_List();

  //! Builds a list containing \p x.
  explicit Expr_List(const Symbol& x);

  //! Builds a list containing \p x1 and \p x2.
  Expr_List(const Expr& x1, const Expr& x2);

  //! Builds a list containing \p x1, \p x2, \p x3, \p x4 and \p x5.
  Expr_List(const Expr& x1, const Expr& x2, const Expr& x3,
	    const Expr& x4, const Expr& x5);

  //! Builds a list containing \p x1, \p x2, \p x3, \p x4 and \p x5.
  Expr_List(const Expr& x1, const Expr& x2, const Expr& x3,
	    const Expr& x4, const Expr& x5, const Expr& x6);

  //! Copy-constructor.
  Expr_List(const Expr_List& x);

  //! Destructor.
  ~Expr_List();

  //! Assignment operator.
  Expr_List& operator=(const Expr_List& x);

  //! Returns the number of expressions of \p *this.
  unsigned nops() const;

  //! Returns the \f$ i \f$-th element of \p *this.
  /*!
    \exception std::out_of_range thrown if
                                 \f$ i \notin \{0, \dotsc, nops() - 1 \}.
  */
  // FIXME: ginac -> wrong answer!
  Expr op(unsigned i) const;

  //! Appends \p x to \p *this.
  Expr_List& append(const Expr& x);

  //! Prepends \p x to \p *this.
  Expr_List& prepend(const Expr& x);

  //! Removes the first element from \p *this.
  /*!
    \exception std::logic_error thrown if \p *this is empty.
  */
  // FIXME: ginac -> segmentation fault!
  Expr_List& remove_first();

private:
  friend class Expr;
  friend class Matrix;

  friend Expr sqrfree(const Expr& x, const Expr_List& y);
  friend Expr lsolve(const Expr_List& x, const Expr_List& y);

  GiNaC::lst l;

  //! Builds the expression corresponding to \p gl.
  Expr_List(const GiNaC::lst& gl);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Expr_List.inlines.hh"

#endif // !defined(PURRS_Expr_List_defs_hh)
