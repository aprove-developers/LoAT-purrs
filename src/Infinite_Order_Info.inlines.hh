/* Weighted_Average_Info class implementation: inline functions.
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

inline Weighted_Average_Info::
Weighted_Average_Info(const Recurrence& associated_first_order_rec,
		      const Expr& infinite_order_weight)
  : associated_first_order_rec_(associated_first_order_rec),
    infinite_order_weight_(infinite_order_weight) {
}

inline
Weighted_Average_Info::Weighted_Average_Info(const Weighted_Average_Info& y)
  : associated_first_order_rec_(y.associated_first_order_rec_),
    original_rhs_(y.original_rhs_),
    infinite_order_weight_(y.infinite_order_weight_),
    first_valid_index_inf_order_(y.first_valid_index_inf_order_) {
}

inline
Weighted_Average_Info::~Weighted_Average_Info() {
}

inline Weighted_Average_Info&
Weighted_Average_Info::operator=(const Weighted_Average_Info& y) {
  associated_first_order_rec_ = y.associated_first_order_rec_;
  original_rhs_ = y.original_rhs_;
  infinite_order_weight_ = y.infinite_order_weight_;
  first_valid_index_inf_order_ = y.first_valid_index_inf_order_;
  return *this;
}

inline const Recurrence&
Weighted_Average_Info::associated_first_order_rec() const {
  return associated_first_order_rec_;
}

inline Recurrence&
Weighted_Average_Info::associated_first_order_rec() {
  return associated_first_order_rec_;
}

inline void
Weighted_Average_Info::set_original_rhs(const Expr& original_rhs) {
  original_rhs_ = original_rhs;
}

inline const Expr&
Weighted_Average_Info::infinite_order_weight() const {
  return infinite_order_weight_;
}

inline Expr&
Weighted_Average_Info::infinite_order_weight() {
  return infinite_order_weight_;
}

inline index_type
Weighted_Average_Info::first_valid_index_inf_order() const {
  return first_valid_index_inf_order_;
}

inline void
Weighted_Average_Info::set_first_valid_index_inf_order(index_type i_c) {
  first_valid_index_inf_order_ = i_c;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Infinite_Order_Info_inlines_hh)
