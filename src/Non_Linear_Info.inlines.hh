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
#include "util.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Non_Linear_Info::Non_Linear_Info(const Recurrence& associated_linear_rec,
				 const std::pair<Number, Expr>& coeff_and_base,
				 const std::vector<Symbol> auxiliary_symbols)
  : associated_linear_rec_(associated_linear_rec),
    coeff_and_base_(coeff_and_base),
    auxiliary_symbols_(auxiliary_symbols) {
}

inline
Non_Linear_Info::Non_Linear_Info(const Non_Linear_Info& y)
  : associated_linear_rec_(y.associated_linear_rec_),
    coeff_and_base_(y.coeff_and_base_),
    auxiliary_symbols_(y.auxiliary_symbols_) {
}

inline
Non_Linear_Info::~Non_Linear_Info() {
}

inline Non_Linear_Info&
Non_Linear_Info::operator=(const Non_Linear_Info& y) {
  associated_linear_rec_ = y.associated_linear_rec_;
  coeff_and_base_ = y.coeff_and_base_;
  auxiliary_symbols_ = y.auxiliary_symbols_;
  return *this;
}

inline const Recurrence&
Non_Linear_Info::associated_linear_rec() const {
  return associated_linear_rec_;
}

inline Recurrence&
Non_Linear_Info::associated_linear_rec() {
  return associated_linear_rec_;
}

inline Number
Non_Linear_Info::coeff_simple_non_linear_rec() const {
  return coeff_and_base_.first;
}

inline Number&
Non_Linear_Info::coeff_simple_non_linear_rec() {
  return coeff_and_base_.first;
}

inline Expr
Non_Linear_Info::base_exp_log() const {
  return coeff_and_base_.second;
}

inline Expr&
Non_Linear_Info::base_exp_log() {
  return coeff_and_base_.second;
}

inline const std::vector<Symbol>&
Non_Linear_Info::auxiliary_symbols() const {
  return auxiliary_symbols_;
}
  
inline std::vector<Symbol>&
Non_Linear_Info::auxiliary_symbols() {
  return auxiliary_symbols_;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Non_Linear_Info_inlines_hh)
