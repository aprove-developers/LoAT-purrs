/* Functional_Equation_Info class implementation: inline functions.
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

#ifndef PURRS_Functional_Equation_Info_inlines_hh
#define PURRS_Functional_Equation_Info_inlines_hh

#include "Functional_Equation_Info.types.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Functional_Equation_Info::Functional_Equation_Info(const Expr& a, unsigned b)
  : coefficient_(a),
    divisor_arg_(b) {
}

inline
Functional_Equation_Info::
Functional_Equation_Info(const Functional_Equation_Info& y)
  : coefficient_(y.coefficient_),
    divisor_arg_(y.divisor_arg_) {
}

inline
Functional_Equation_Info::~Functional_Equation_Info() {
}

inline Functional_Equation_Info&
Functional_Equation_Info::operator=(const Functional_Equation_Info& y) {
  coefficient_ = y.coefficient_;
  divisor_arg_ = y.divisor_arg_;
  return *this;
}

inline Expr
Functional_Equation_Info::coefficients() const {
  return coefficient_;
}

inline Expr&
Functional_Equation_Info::coefficients() {
  return coefficient_;
}

inline unsigned
Functional_Equation_Info::divisor_arg() const {
  return divisor_arg_;
}

inline unsigned&
Functional_Equation_Info::divisor_arg() {
  return divisor_arg_;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Functional_Equation_Info_inlines_hh)
