/* Declaration of the Blackboard class: a class to keep track
   of symbol definitions.
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

#ifndef PURRS_Blackboard_defs_hh
#define PURRS_Blackboard_defs_hh 1

#include "Blackboard.types.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include <map>
#include <deque>
#include <iosfwd>

namespace Parma_Recurrence_Relation_Solver {

/*!
  This class keeps track of a set of symbol definitions.
  A symbol can be associated to an expression that, in turn,
  may contain other symbols.
*/
class Blackboard {
public:
  //! Default constructor: created the empty blackboard.
  Blackboard();

  //! Copy-constructor.
  Blackboard(const Blackboard& y);

  //! Destructor.
  ~Blackboard();

  //! Assignment operator.
  Blackboard& operator=(const Blackboard& y);

  //! Returns a new symbol \f$ z \f$ and records the equation \f$ z = e \f$.
  Symbol insert_definition(const Expr& e);

  //! \brief
  //! Returns the right-hand side of the auxiliary equation \f$ z = e \f$,
  //! if such an auxiliary equation exists;
  //! returns the expression \f$ z \f$ otherwise.
  Expr get_definition(const Symbol& z) const;

  //! Rewrite \f$ e \f$ according to the definitions in \p *this.
  Expr rewrite(const Expr& e) const;

  //! \brief
  //! Compute the size norm of \f$ e \f$ according to the definitions
  //! in \p *this.
  unsigned size_norm(const Expr& e) const;

  //! \brief
  //! Compute the size norm of \f$ e \f$ according to the definitions
  //! in \p *this.
  unsigned size_norm(const Symbol& s) const;

  //! Approximate \f$ e \f$ according to the definitions in \p *this.
  Expr approximate(const Expr& e) const;

  void Blackboard::dump(std::ostream& s) const;

private:
  template <typename T>
  class Cached {
  public:
    Cached();

    unsigned long timestamp;
    T value;
  };

  Expr rewrite(Cached<Expr>& ce) const;
  unsigned size_norm(Cached<unsigned>& ce) const;
  Expr approximate(Cached<Expr>& ce) const;

  class Definition {
  public:
    explicit Definition(const Expr& e);

    Expr rhs;
    Cached<unsigned> size;
    Cached<Expr> approximation;
    Cached<Expr> expansion;
  };

  std::map<Symbol, unsigned> index;

  mutable std::deque<Definition> definitions;

  mutable unsigned long timestamp;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Blackboard.inlines.hh"

#endif // !defined(PURRS_Blackboard_defs_hh)
