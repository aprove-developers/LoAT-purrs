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
Finite_Order_Info::Finite_Order_Info(int k, const std::vector<unsigned>& decs,
				     const std::vector<Expr>& coeffs)
  : order(k),
    decrements(decs),
    coefficients(coeffs) {
}

inline
Finite_Order_Info::Finite_Order_Info(const Finite_Order_Info& y) 
  : order(y.order),
    decrements(y.decrements),
    initial_conditions(y.initial_conditions),
    coefficients(y.coefficients) {
}

inline
Finite_Order_Info::~Finite_Order_Info() {
}

inline Finite_Order_Info&
Finite_Order_Info::operator=(const Finite_Order_Info& y) { 
  order = y.order;
  decrements = y.decrements;
  initial_conditions = y.initial_conditions;
  coefficients = y.coefficients;
  return *this;
}

inline int
Finite_Order_Info::get_order() const {
  return order;
}

inline const std::vector<unsigned>&
Finite_Order_Info::get_decrements() const {
  return decrements; 
}

inline const std::vector<unsigned>&
Finite_Order_Info::get_initial_conditions() const {
  return initial_conditions;
}

inline void
Finite_Order_Info::set_decrements(const std::vector<unsigned> decs) {
  copy(decs.begin(), decs.end(), inserter(decrements, decrements.begin()));
}

inline void
Finite_Order_Info::set_initial_conditions(const std::vector<unsigned> i_c) {
  copy(i_c.begin(), i_c.end(),
       inserter(initial_conditions, initial_conditions.begin()));
}

inline void
Finite_Order_Info::set_coefficients(const std::vector<Expr> coeffs) {
  copy(coeffs.begin(), coeffs.end(),
       inserter(coefficients, coefficients.begin()));
}

inline void
Finite_Order_Info::add_decrement(unsigned d) {
  decrements.push_back(d);
}

inline void
Finite_Order_Info::add_initial_condition(unsigned i) {
  initial_conditions.push_back(i);
}

inline void
Finite_Order_Info::add_coefficient(const Expr& c) {
  coefficients.push_back(c);
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Finite_Order_Info_inlines_hh)
