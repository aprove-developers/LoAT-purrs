/* Finite_Order_Info class implementation: inline functions.
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

#ifndef PURRS_Finite_Order_Info_inlines_hh
#define PURRS_Finite_Order_Info_inlines_hh

#include "Finite_Order_Info.types.hh"
#include <vector>
#include <algorithm>
#include <iterator>

namespace Parma_Recurrence_Relation_Solver {

inline
Finite_Order_Info::Finite_Order_Info(int k, unsigned i_c,
				     const std::vector<Expr>& coeffs)
  : order_(k),
    first_initial_condition_(i_c),
    coefficients_(coeffs) {
}

inline
Finite_Order_Info::Finite_Order_Info(const Finite_Order_Info& y) 
  : order_(y.order_),
    first_initial_condition_(y.first_initial_condition_),
    coefficients_(y.coefficients_) {
}

inline
Finite_Order_Info::~Finite_Order_Info() {
}

inline Finite_Order_Info&
Finite_Order_Info::operator=(const Finite_Order_Info& y) { 
  order_ = y.order_;
  first_initial_condition_ = y.first_initial_condition_;
  coefficients_ = y.coefficients_;
  return *this;
}

inline unsigned int
Finite_Order_Info::order() const {
  return order_;
}

inline unsigned int&
Finite_Order_Info::order() {
  return order_;
}

inline unsigned
Finite_Order_Info::first_initial_condition() const {
  return first_initial_condition_;
}

inline unsigned&
Finite_Order_Info::first_initial_condition() {
  return first_initial_condition_;
}

inline const std::vector<Expr>&
Finite_Order_Info::coefficients() const {
  return coefficients_;
}

inline std::vector<Expr>&
Finite_Order_Info::coefficients() {
  return coefficients_;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Finite_Order_Info_inlines_hh)
