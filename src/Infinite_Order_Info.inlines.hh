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
					 const Expr& weight,
					 unsigned i_c) 
  : rhs_transformed_in_first_order_var_coeffs_(new_rhs),
    weight_(weight),
    infinite_order_fwdr_(i_c) {
}

inline
Infinite_Order_Info::Infinite_Order_Info(const Infinite_Order_Info& y)
  : rhs_transformed_in_first_order_var_coeffs_
(y.rhs_transformed_in_first_order_var_coeffs_),
    weight_(y.weight_),
    infinite_order_fwdr_(y.infinite_order_fwdr_) {
}

inline
Infinite_Order_Info::~Infinite_Order_Info() {
}

inline Infinite_Order_Info&
Infinite_Order_Info::operator=(const Infinite_Order_Info& y) {
  rhs_transformed_in_first_order_var_coeffs_
    = y.rhs_transformed_in_first_order_var_coeffs_;
  weight_ = y.weight_;
  infinite_order_fwdr_ = y.infinite_order_fwdr_;
  return *this;
}

inline Expr
Infinite_Order_Info::rhs_transformed_in_first_order_var_coeffs() const {
  return rhs_transformed_in_first_order_var_coeffs_;
}

inline Expr&
Infinite_Order_Info::rhs_transformed_in_first_order_var_coeffs() {
  return rhs_transformed_in_first_order_var_coeffs_;
}

inline Expr
Infinite_Order_Info::weight() const {
  return weight_;
}

inline Expr&
Infinite_Order_Info::weight() {
  return weight_;
}

inline unsigned
Infinite_Order_Info::infinite_order_fwdr() const {
  return infinite_order_fwdr_;
}

inline unsigned&
Infinite_Order_Info::infinite_order_fwdr() {
    return infinite_order_fwdr_;
}

inline void
Infinite_Order_Info::set_infinite_order_fwdr(unsigned i_c) {
  infinite_order_fwdr_ = i_c;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Infinite_Order_Info_inlines_hh)
