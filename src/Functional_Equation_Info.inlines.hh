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
#include <sstream>

namespace Parma_Recurrence_Relation_Solver {

inline
Functional_Equation_Info::
Functional_Equation_Info(const std::map<Number, Expr>& hom_terms)
  : homogeneous_terms(hom_terms),
    applicability_condition_(1),
    definition_Sc_() {
}

inline
Functional_Equation_Info::
Functional_Equation_Info(const Functional_Equation_Info& y)
  : homogeneous_terms(y.homogeneous_terms),
    applicability_condition_(y.applicability_condition_),
    definition_Sc_(y.definition_Sc_) {
}

inline
Functional_Equation_Info::~Functional_Equation_Info() {
}

inline Functional_Equation_Info&
Functional_Equation_Info::operator=(const Functional_Equation_Info& y) {
  homogeneous_terms = y.homogeneous_terms;
  applicability_condition_ = y.applicability_condition_;
  definition_Sc_ = y.definition_Sc_;
  return *this;
}

inline index_type
Functional_Equation_Info::applicability_condition() const {
  return applicability_condition_;
}

inline void
Functional_Equation_Info::set_applicability_condition(index_type c) {
  applicability_condition_ = c;
}

inline Functional_Equation_Info::ht_iterator
Functional_Equation_Info::ht_begin() {
  return homogeneous_terms.begin();
}

inline Functional_Equation_Info::ht_iterator
Functional_Equation_Info::ht_end() {
  return homogeneous_terms.end();
}

inline Functional_Equation_Info::ht_const_iterator
Functional_Equation_Info::ht_begin() const {
  return homogeneous_terms.begin();
}

inline Functional_Equation_Info::ht_const_iterator
Functional_Equation_Info::ht_end() const {
  return homogeneous_terms.end();
}

inline index_type
Functional_Equation_Info::rank() const {
  return homogeneous_terms.size();
}

inline std::string
Functional_Equation_Info::definition_Sc() const {
  return definition_Sc_;
}

inline void
Functional_Equation_Info::set_definition_Sc() {
  assert(!homogeneous_terms.empty());
  assert(homogeneous_terms.size() == 1);
  const Number& divisor = ht_begin()->first;
  std::ostringstream s;
  if (divisor != 2)
    s << "Sc(n, " << divisor
      << ") = [ n/" << divisor <<"^[log(n)/log(" << divisor << ")] ]";
  definition_Sc_ = s.str();
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Functional_Equation_Info_inlines_hh)
