/* *******************: a recurrence: inline functions.
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

#ifndef PURRS_Number_inlines_hh
#define PURRS_Number_inlines_hh

#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

Number
operator+(const Number& lh, const Number& rh) {
  return lh + rh;
}

Number
operator-(const Number& lh, const Number& rh) {
  return lh - rh;
}

Number
operator*(const Number& lh, const Number& rh) {
  return lh * rh;
}

Number
operator/(const Number& lh, const Number& rh) {
  return lh / rh;
}

Number
operator+(const Number& lh) {
  return +lh;
}

Number
operator-(const Number& lh) {
  return -lh;
}

inline bool
Number::operator==(const Number& num) const {
  return n == num.n;
}

inline bool
Number::operator!=(const Number& num) const {
  return n != num.n;
}

inline bool
Number::operator>(const Number& num) const {
  return n > num.n;
}

inline bool
Number::operator<(const Number& num) const {
  return n < num.n;
}

inline bool
Number::operator>=(const Number& num) const {
  return n >= num.n;
}

inline bool
Number::operator<=(const Number& num) const {
  return n <= num.n;
}

// inline relational
// Number::operator==(const Symbol& lh, const Number& rh) const {
//   return lh == rh;
// }

inline
Number::Number() {
}

inline
Number::Number(int i)
  : n(i) {
}

inline
Number::Number(Expr& exp)
  : n(exp.e) {
}

inline
Number::Number(const Number& num)
  : n(num.n) {
};

inline Number&
Number::operator=(const Number& num) {
  n = num.n;
  return *this;
};

inline
Number::Number(const GiNaC::numeric& gn)
  : n(gn) {
}

inline
Number::~Number() {
}

inline Number
abs(const Number& n) {
  return abs(n);
}

inline Number
sqrt(const Number &n) {
  return sqrt(n);
};

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Number_inlines_hh)
