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
Finite_Order_Info::Finite_Order_Info(int k,const std::vector<Expr>& coeffs,
				     unsigned gcd, unsigned i_c,
				     Expr prod_factor)
  : order_(k),
    coefficients_(coeffs),
    gcd_among_decrements_(gcd),
    first_well_defined_rhs_linear_(i_c),
    product_factor_(prod_factor) {
}

inline
Finite_Order_Info::Finite_Order_Info(const Finite_Order_Info& y) 
  : order_(y.order_),
    coefficients_(y.coefficients_),
    gcd_among_decrements_(y.gcd_among_decrements_),
    first_well_defined_rhs_linear_(y.first_well_defined_rhs_linear_),
    product_factor_(y.product_factor_) {
}

inline
Finite_Order_Info::~Finite_Order_Info() {
}

inline Finite_Order_Info&
Finite_Order_Info::operator=(const Finite_Order_Info& y) { 
  order_ = y.order_;
  coefficients_ = y.coefficients_;
  gcd_among_decrements_ = y.gcd_among_decrements_;
  first_well_defined_rhs_linear_ = y.first_well_defined_rhs_linear_;
  product_factor_ = y.product_factor_;
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
Finite_Order_Info::first_well_defined_rhs_linear() const {
  return first_well_defined_rhs_linear_;
}

inline unsigned&
Finite_Order_Info::first_well_defined_rhs_linear() {
  return first_well_defined_rhs_linear_;
}

inline void
Finite_Order_Info::set_first_well_defined_rhs_linear(unsigned i_c) {
  first_well_defined_rhs_linear_ = i_c;
}

inline const std::vector<Expr>&
Finite_Order_Info::coefficients() const {
  return coefficients_;
}

inline std::vector<Expr>&
Finite_Order_Info::coefficients() {
  return coefficients_;
}

inline unsigned
Finite_Order_Info::gcd_among_decrements() const {
  return gcd_among_decrements_;
}

inline unsigned&
Finite_Order_Info::gcd_among_decrements() {
  return gcd_among_decrements_;
}

inline Expr&
Finite_Order_Info::product_factor() {
  return product_factor_;
}

inline Expr
Finite_Order_Info::product_factor() const {
  return product_factor_;
}

inline void
Finite_Order_Info::set_product_factor(const Expr& x) {
  product_factor_ = x;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Finite_Order_Info_inlines_hh)
