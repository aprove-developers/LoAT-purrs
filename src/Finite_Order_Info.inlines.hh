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

namespace Parma_Recurrence_Relation_Solver {

inline
Finite_Order_Info::Finite_Order_Info(index_type k,
				     const std::vector<Expr>& coeffs,
				     index_type first_valid_index,
				     unsigned int gcd)
  : order_(k),
    coefficients_(coeffs),
    first_valid_index_(first_valid_index),
    gcd_among_decrements_(gcd),
    product_factor_(0),
    applied_order_reduction_(false) {
}

inline
Finite_Order_Info::Finite_Order_Info(const Finite_Order_Info& y) 
  : order_(y.order_),
    coefficients_(y.coefficients_),
    first_valid_index_(y.first_valid_index_),
    gcd_among_decrements_(y.gcd_among_decrements_),
    product_factor_(y.product_factor_),
    applied_order_reduction_(y.applied_order_reduction_) {
}

inline
Finite_Order_Info::~Finite_Order_Info() {
}

inline Finite_Order_Info&
Finite_Order_Info::operator=(const Finite_Order_Info& y) { 
  order_ = y.order_;
  coefficients_ = y.coefficients_;
  first_valid_index_ = y.first_valid_index_;
  gcd_among_decrements_ = y.gcd_among_decrements_;
  product_factor_ = y.product_factor_;
  applied_order_reduction_ = y.applied_order_reduction_;
  return *this;
}

inline index_type
Finite_Order_Info::order() const {
  return order_;
}

inline index_type
Finite_Order_Info::first_valid_index() const {
  return first_valid_index_;
}

inline void
Finite_Order_Info::set_first_valid_index(index_type i_c) {
  first_valid_index_ = i_c;
}

inline const std::vector<Expr>&
Finite_Order_Info::coefficients() const {
  return coefficients_;
}

inline std::vector<Expr>&
Finite_Order_Info::coefficients() {
  return coefficients_;
}

inline unsigned int
Finite_Order_Info::gcd_among_decrements() const {
  return gcd_among_decrements_;
}

inline const Expr&
Finite_Order_Info::product_factor() const {
  return product_factor_;
}

inline Expr&
Finite_Order_Info::product_factor() {
  return product_factor_;
}

inline void
Finite_Order_Info::set_product_factor(const Expr& x) {
  product_factor_ = x;
}

inline bool
Finite_Order_Info::applied_order_reduction() const {
  return applied_order_reduction_;
}

inline void
Finite_Order_Info::set_order_reduction() {
  applied_order_reduction_ = true;
}

inline void
Finite_Order_Info::unset_order_reduction() {
  applied_order_reduction_ = false;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Finite_Order_Info_inlines_hh)
