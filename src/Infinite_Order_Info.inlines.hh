/* Infinite_Order_Info class implementation: inline functions.
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

#ifndef PURRS_Infinite_Order_Info_inlines_hh
#define PURRS_Infinite_Order_Info_inlines_hh

#include "Infinite_Order_Info.types.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Infinite_Order_Info::Infinite_Order_Info(const Expr& new_rhs,
					 const Expr& coefficient,
					 const Expr& inhomogeneous,
					 const Expr& first_element) 
  : rhs_transformed_in_first_order_(new_rhs),
    coeff_first_order_(coefficient),
    inhomog_first_order_(inhomogeneous),
    value_of_first_element_(first_element) {
}

inline
Infinite_Order_Info::Infinite_Order_Info(const Infinite_Order_Info& y)
  : rhs_transformed_in_first_order_
(y.rhs_transformed_in_first_order_),
    coeff_first_order_(y.coeff_first_order_),
    inhomog_first_order_(y.inhomog_first_order_),
    value_of_first_element_(y.value_of_first_element_) {
}

inline
Infinite_Order_Info::~Infinite_Order_Info() {
}

inline Infinite_Order_Info&
Infinite_Order_Info::operator=(const Infinite_Order_Info& y) {
  rhs_transformed_in_first_order_
    = y.rhs_transformed_in_first_order_;
  coeff_first_order_ = y.coeff_first_order_;
  inhomog_first_order_ = y.inhomog_first_order_;
  value_of_first_element_ = y.value_of_first_element_;
  return *this;
}

inline Expr
Infinite_Order_Info::rhs_transformed_in_first_order() const {
  return rhs_transformed_in_first_order_;
}

inline Expr&
Infinite_Order_Info::rhs_transformed_in_first_order() {
  return rhs_transformed_in_first_order_;
}

inline Expr
Infinite_Order_Info::coeff_first_order() const {
  return coeff_first_order_;
}

inline Expr&
Infinite_Order_Info::coeff_first_order() {
  return coeff_first_order_;
}

inline Expr
Infinite_Order_Info::inhomog_first_order() const {
  return inhomog_first_order_;
}

inline Expr&
Infinite_Order_Info::inhomog_first_order() {
  return inhomog_first_order_;
}

inline Expr
Infinite_Order_Info::value_of_first_element() const {
  return value_of_first_element_;
}

inline Expr&
Infinite_Order_Info::value_of_first_element() {
  return value_of_first_element_;
}

// inline unsigned
// Infinite_Order_Info::infinite_order_fwdr() const {
//   return infinite_order_fwdr_;
// }

// inline unsigned&
// Infinite_Order_Info::infinite_order_fwdr() {
//     return infinite_order_fwdr_;
// }

// inline void
// Infinite_Order_Info::set_infinite_order_fwdr(unsigned i_c) {
//   infinite_order_fwdr_ = i_c;
// }

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Infinite_Order_Info_inlines_hh)
