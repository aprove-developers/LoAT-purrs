/* Order_Reduction_Info class implementation: inline functions.
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

#ifndef PURRS_Order_Reduction_Info_inlines_hh
#define PURRS_Order_Reduction_Info_inlines_hh

#include "Order_Reduction_Info.types.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Order_Reduction_Info::Order_Reduction_Info(const Expr& rhs, unsigned gcd,
					   const Symbol& r)
  : old_recurrence_rhs_(rhs),
    gcd_decrements_old_rhs_(gcd),
    symbol_for_mod_(r),
    verified_one_time_(true) {
}

inline
Order_Reduction_Info::Order_Reduction_Info(const Order_Reduction_Info& y)
  : old_recurrence_rhs_(y.old_recurrence_rhs_),
    gcd_decrements_old_rhs_(y.gcd_decrements_old_rhs_),
    symbol_for_mod_(y.symbol_for_mod_),
    solution_order_reduced_(y.solution_order_reduced_),
    verified_one_time_(y.verified_one_time_) {
}

inline
Order_Reduction_Info::~Order_Reduction_Info() {
}

inline Order_Reduction_Info&
Order_Reduction_Info::operator=(const Order_Reduction_Info& y) {
  old_recurrence_rhs_ = y.old_recurrence_rhs_;
  gcd_decrements_old_rhs_ = y.gcd_decrements_old_rhs_;
  symbol_for_mod_ = y.symbol_for_mod_;
  solution_order_reduced_ = y.solution_order_reduced_;
  verified_one_time_ = y.verified_one_time_;
  return *this;
}

inline Expr
Order_Reduction_Info::old_recurrence_rhs() const {
  return old_recurrence_rhs_;
}

inline Expr&
Order_Reduction_Info::old_recurrence_rhs() {
  return old_recurrence_rhs_;
}

inline unsigned
Order_Reduction_Info::gcd_decrements_old_rhs() const {
  return gcd_decrements_old_rhs_;
}

inline unsigned&
Order_Reduction_Info::gcd_decrements_old_rhs() {
  return gcd_decrements_old_rhs_;
}

inline void
Order_Reduction_Info::set_gcd_decrements_old_rhs(unsigned g) {
  gcd_decrements_old_rhs_ = g;
}

inline Symbol
Order_Reduction_Info::symbol_for_mod() const {
  return symbol_for_mod_;
}

inline Symbol&
Order_Reduction_Info::symbol_for_mod() {
  return symbol_for_mod_;
}

inline Expr
Order_Reduction_Info::solution_order_reduced() const {
  return solution_order_reduced_;
}

inline Expr&
Order_Reduction_Info::solution_order_reduced() {
  return solution_order_reduced_;
}

inline void
Order_Reduction_Info::set_solution_order_reduced(const Expr& e) {
  solution_order_reduced_ = e;
}

inline bool
Order_Reduction_Info::verified_one_time() const {
  return verified_one_time_;
}

inline void
Order_Reduction_Info::not_verified_one_time() {
  verified_one_time_ = false;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Order_Reduction_Info_inlines_hh)
