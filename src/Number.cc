/* Symbol class implementation (non-inline functions).
   Copyright (C) 2001, 2002 Roberto Bagnara <bagnara@cs.unipr.it>

This file is part of the Parma Polyhedra Library (PPL).

The PPL is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The PPL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.

For the most up-to-date information see the Parma Polyhedra Library
site: http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#include "Number.defs.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

const PURRS::Number PURRS::I = GiNaC::I;

int
PURRS::to_int(const Number& n) {
  return to_int(n);
}

long
PURRS::to_long(const Number& n) {
  return to_long(n);
}


bool
PURRS::Number::is_positive() const {
  return n.is_positive();
}

bool
PURRS::Number::is_integer() const {
  return n.is_integer();
}

bool
PURRS::Number::is_pos_integer() const {
  return n.is_pos_integer();
}

bool
PURRS::Number::is_nonneg_integer() const {
  return n.is_nonneg_integer();
}

bool
PURRS::Number::is_even() const {
  return n.is_even();
}

bool
PURRS::Number::is_odd() const {
  return n.is_odd();
}

bool
PURRS::Number::is_prime() const {
  return n.is_prime();
}

bool
PURRS::Number::is_rational() const {
  return n.is_rational();
}

bool
PURRS::Number::is_real() const {
  return n.is_real();
}

bool
PURRS::Number::is_cinteger() const {
  return n.is_cinteger();
}

bool
PURRS::Number::is_crational() const {
  return n.is_crational();
}
