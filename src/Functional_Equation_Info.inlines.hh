/* Functional_Equation_Info class implementation: inline functions.
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

#ifndef PURRS_Functional_Equation_Info_inlines_hh
#define PURRS_Functional_Equation_Info_inlines_hh

#include "Functional_Equation_Info.types.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Functional_Equation_Info::
Functional_Equation_Info(unsigned int k, const std::vector<Expr>& coeffs,
			 const std::vector<Number>& divisors, unsigned c)
  : rank_(k),
    coefficients_fe_(coeffs),
    divisors_arg_(divisors),
    applicability_condition_(c) {
}

inline
Functional_Equation_Info::
Functional_Equation_Info(const Functional_Equation_Info& y)
  : rank_(y.rank_),
    coefficients_fe_(y.coefficients_fe_),
    divisors_arg_(y.divisors_arg_),
    applicability_condition_(y.applicability_condition_) {
}

inline
Functional_Equation_Info::~Functional_Equation_Info() {
}

inline Functional_Equation_Info&
Functional_Equation_Info::operator=(const Functional_Equation_Info& y) {
  rank_ = y.rank_;
  coefficients_fe_ = y.coefficients_fe_;
  divisors_arg_ = y.divisors_arg_;
  applicability_condition_ = y.applicability_condition_;
  return *this;
}

inline unsigned int
Functional_Equation_Info::rank() const {
  return rank_;
}

inline unsigned int&
Functional_Equation_Info::rank() {
  return rank_;
}

inline const std::vector<Expr>&
Functional_Equation_Info::coefficients_fe() const {
  return coefficients_fe_;
}

inline std::vector<Expr>&
Functional_Equation_Info::coefficients_fe() {
  return coefficients_fe_;
}

inline const std::vector<Number>&
Functional_Equation_Info::divisors_arg() const {
  return divisors_arg_;
}

inline std::vector<Number>&
Functional_Equation_Info::divisors_arg() {
  return divisors_arg_;
}

inline unsigned
Functional_Equation_Info::applicability_condition() const {
  return applicability_condition_;
}

inline unsigned&
Functional_Equation_Info::applicability_condition() {
  return applicability_condition_;
}

inline void
Functional_Equation_Info::set_applicability_condition(unsigned c) {
  applicability_condition_ = c;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Functional_Equation_Info_inlines_hh)
