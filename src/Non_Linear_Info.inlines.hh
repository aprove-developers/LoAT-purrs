/* Non_Linear_Info class implementation: inline functions.
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

#ifndef PURRS_Non_Linear_Info_inlines_hh
#define PURRS_Non_Linear_Info_inlines_hh

#include "Non_Linear_Info.types.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Non_Linear_Info::Non_Linear_Info(const Expr& rhs, const Expr& new_rhs,
				 const Expr& base_exp_log,
				 const std::vector<Symbol> auxiliary_symbols)
  : original_recurrence_rhs_(rhs),
    rhs_transformed_in_linear_(new_rhs),
    base_exp_log_(base_exp_log),
    auxiliary_symbols_(auxiliary_symbols),
    order_if_linear_(0),
    first_i_c_if_linear_(0) {
}

inline
Non_Linear_Info::Non_Linear_Info(const Non_Linear_Info& y)
  : original_recurrence_rhs_(y.original_recurrence_rhs_),
    rhs_transformed_in_linear_(y.rhs_transformed_in_linear_),
    base_exp_log_(y.base_exp_log_),
    auxiliary_symbols_(y.auxiliary_symbols_),
    order_if_linear_(y.order_if_linear_),
    first_i_c_if_linear_(y.first_i_c_if_linear_) {
}

inline
Non_Linear_Info::~Non_Linear_Info() {
}

inline Non_Linear_Info&
Non_Linear_Info::operator=(const Non_Linear_Info& y) {
  original_recurrence_rhs_ = y.original_recurrence_rhs_;
  rhs_transformed_in_linear_ = y.rhs_transformed_in_linear_;
  base_exp_log_ = y.base_exp_log_;
  auxiliary_symbols_ = y.auxiliary_symbols_;
  order_if_linear_ = y.order_if_linear_;
  first_i_c_if_linear_ = y.first_i_c_if_linear_;
  return *this;
}

inline Expr
Non_Linear_Info::original_recurrence_rhs() const {
  return original_recurrence_rhs_;
}

inline Expr&
Non_Linear_Info::original_recurrence_rhs() {
  return original_recurrence_rhs_;
}

inline Expr
Non_Linear_Info::rhs_transformed_in_linear() const {
  return rhs_transformed_in_linear_;
}

inline Expr&
Non_Linear_Info::rhs_transformed_in_linear() {
  return rhs_transformed_in_linear_;
}

inline Expr
Non_Linear_Info::base_exp_log() const {
  return base_exp_log_;
}

inline Expr&
Non_Linear_Info::base_exp_log() {
  return base_exp_log_;
}

inline const std::vector<Symbol>&
Non_Linear_Info::auxiliary_symbols() const {
  return auxiliary_symbols_;
}
  
inline std::vector<Symbol>&
Non_Linear_Info::auxiliary_symbols() {
  return auxiliary_symbols_;
}

inline unsigned int
Non_Linear_Info::order_if_linear() const {
  return order_if_linear_;
}

inline unsigned int&
Non_Linear_Info::order_if_linear() {
  return order_if_linear_;
}

inline void
Non_Linear_Info::set_order_if_linear(unsigned int x) {
  order_if_linear_ = x;
}

inline unsigned
Non_Linear_Info::first_i_c_if_linear() const {
  return first_i_c_if_linear_;
}

inline unsigned&
Non_Linear_Info::first_i_c_if_linear() {
  return first_i_c_if_linear_;
}

inline void
Non_Linear_Info::set_first_i_c_if_linear(unsigned i_c) {
  first_i_c_if_linear_ = i_c;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Non_Linear_Info_inlines_hh)
