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
Expr_List::Expr_List(const Symbol& x)
  : l(x.s) {
}

inline
Expr_List::Expr_List(const Expr& e)
  : l(static_cast<const GiNaC::ex>(e)) {
}

inline
Expr_List::Expr_List(const Expr& e1, const Expr& e2)
  : l(static_cast<const GiNaC::ex>(e1), static_cast<const GiNaC::ex>(e2)) {
};

inline
Expr_List::Expr_List(const Expr& e1, const Expr& e2, const Expr& e3,
		     const Expr& e4, const Expr& e5)
  : l(static_cast<const GiNaC::ex>(e1),
      static_cast<const GiNaC::ex>(e2),
      static_cast<const GiNaC::ex>(e3),
      static_cast<const GiNaC::ex>(e4),
      static_cast<const GiNaC::ex>(e5)) {
};

inline
Expr_List::Expr_List(const Expr& e1, const Expr& e2, const Expr& e3,
		     const Expr& e4, const Expr& e5, const Expr& e6)
  : l(static_cast<const GiNaC::ex>(e1),
      static_cast<const GiNaC::ex>(e2),
      static_cast<const GiNaC::ex>(e3),
      static_cast<const GiNaC::ex>(e4),
      static_cast<const GiNaC::ex>(e5),
      static_cast<const GiNaC::ex>(e6)) {
};

inline
Expr_List::Expr_List(const Expr_List& x)
  : l(x.l) {
};

inline Expr_List&
Expr_List::operator=(const Expr_List& x) {
  l = x.l;
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
Expr_List::append(const Expr& x) {
  l.append(static_cast<const GiNaC::ex>(x));
  return *this;
}

inline Expr_List&
Expr_List::prepend(const Expr& x) {
  l.prepend(static_cast<const GiNaC::ex>(x));
  return *this;
}

inline Expr_List&
Expr_List::remove_first() {
  l.remove_first();
  return *this;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Expr_List_inlines_hh)
