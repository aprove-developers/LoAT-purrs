/* A class for holding expressions and attributes we may have or may have
   not computed (and that we do not want to compute more than once).
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

#ifndef PURRS_Cached_Expr_defs_hh
#define PURRS_Cached_Expr_defs_hh 1

#include "Recurrence.types.hh"
#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Cached_Expr {
public:
  //! Default constructor.
  Cached_Expr();

  //! Copy-constructor.
  Cached_Expr(const Cached_Expr& y);

  //! Destructor.
  ~Cached_Expr();

  //! Assignment operator.
  Cached_Expr& operator=(const Cached_Expr& y);

  //! \brief
  //! Returns <CODE>true</CODE> if <CODE>expression_</CODE>
  //! has been setted; returns <CODE>false</CODE> otherwise.
  bool has_expression() const;

  //! Returns <CODE>expression_</CODE>.
  const Expr& expression() const;

  //! Sets <CODE>expression_</CODE> with \p e.
  void set_expression(const Expr& e);

  //! \brief
  //! Substitutes all <EM>bad</EM> symbols contained in
  //! \p (*this).expression() with <EM>good</EM> symbols that
  //! are not yet used. Besides, substitutes the name of the bad
  //! symbols also in the blackboard of the rcurrence \p rec.
  Expr remove_bad_symbols(const Recurrence& rec) const;

  void set_size(unsigned long n);
  unsigned long size();
  bool has_size() const;

private:
  //! Indicates if there is a value in the data <CODE>expression_</CODE>.
  bool has_expression_;

  //! The value of the expression.
  Expr expression_;

  bool has_size_;
  unsigned long size_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Cached_Expr.inlines.hh"

#endif // !defined(PURRS_Cached_Expr_defs_hh)


