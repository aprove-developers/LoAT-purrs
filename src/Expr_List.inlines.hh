/* Expr_List class implementation: inline functions.
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

#ifndef PURRS_Expr_List_inlines_hh
#define PURRS_Expr_List_inlines_hh

#include "Symbol.defs.hh"
#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Expr_List::Expr_List() {
}

inline
Expr_List::Expr_List(const Symbol& symb)
  : l(symb.s) {
};

inline
Expr_List::Expr_List(const Expr& exp1, const Expr& exp2)
  : l(exp1.e, exp2.e) {
};

inline
Expr_List::Expr_List(const Expr& exp1, const Expr& exp2, const Expr& exp3,
		     const Expr& exp4, const Expr& exp5)
  : l(exp1.e, exp2.e, exp3.e, exp4.e, exp5.e) {
};

inline
Expr_List::Expr_List(const Expr_List& lst)
  : l(lst.l) {
};

inline Expr_List&
Expr_List::operator=(const Expr_List& lst) {
  l = lst.l;
  return *this;
};

inline
Expr_List::Expr_List(const GiNaC::lst& gl)
  : l(gl) {
}

inline
Expr_List::~Expr_List() {
}

inline unsigned
Expr_List::nops() const {
  return l.nops();
}

inline Expr
Expr_List::op(unsigned i) const {
  return l.op(i);
}

inline Expr_List&
Expr_List::append(const Expr& exp) {
  l.append(exp.e);
  return *this;
}

inline Expr_List&
Expr_List::prepend(const Expr& exp) {
  l.prepend(exp.e);
  return *this;
}

inline Expr_List&
Expr_List::remove_first() {
  l.remove_first();
  return *this;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Expr_List_inlines_hh)