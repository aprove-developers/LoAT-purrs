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

inline Infinite_Order_Info::
Infinite_Order_Info(const Recurrence& associated_first_order_rec,
		    const Expr& infinite_order_weight)
  : associated_first_order_rec_(associated_first_order_rec),
    infinite_order_weight_(infinite_order_weight) {
}

inline
Infinite_Order_Info::Infinite_Order_Info(const Infinite_Order_Info& y)
  : associated_first_order_rec_(y.associated_first_order_rec_),
    infinite_order_weight_(y.infinite_order_weight_) {
}

inline
Infinite_Order_Info::~Infinite_Order_Info() {
}

inline Infinite_Order_Info&
Infinite_Order_Info::operator=(const Infinite_Order_Info& y) {
  associated_first_order_rec_ = y.associated_first_order_rec_;
  infinite_order_weight_ = y.infinite_order_weight_;
  return *this;
}

inline const Recurrence&
Infinite_Order_Info::associated_first_order_rec() const {
  return associated_first_order_rec_;
}

inline Recurrence&
Infinite_Order_Info::associated_first_order_rec() {
  return associated_first_order_rec_;
}

inline const Expr&
Infinite_Order_Info::infinite_order_weight() const {
  return infinite_order_weight_;
}

inline Expr&
Infinite_Order_Info::infinite_order_weight() {
  return infinite_order_weight_;
}

inline index_type
Infinite_Order_Info::first_valid_index_inf_order() const {
  return first_valid_index_inf_order_;
}

inline void
Infinite_Order_Info::set_first_valid_index_inf_order(index_type i_c) {
  first_valid_index_inf_order_ = i_c;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Infinite_Order_Info_inlines_hh)
