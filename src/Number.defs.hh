/* *****************
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

#ifndef PURRS_Number_defs_hh
#define PURRS_Number_defs_hh 1

#include "Number.types.hh"
#include "Expr.types.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

// binary arithmetic operators Number with Number
Number operator+(const Number& lh, const Number& rh);
Number operator-(const Number& lh, const Number& rh);
Number operator*(const Number& lh, const Number& rh);
Number operator/(const Number& lh, const Number& rh);

// unary operators
Number operator+(const Number &lh);
Number operator-(const Number &lh);

// increment / decrement operators
Number& operator++(Number& rh);
Number& operator--(Number& rh);
// FIXME: int i e' superflua? Ma se la tolgo c'e' ambiguita' con quella sopra
Number operator++(Number& lh, int i);
Number operator--(Number& lh, int i);

class Parma_Recurrence_Relation_Solver::Number {
public:
  //! Ordinary copy-constructor.
  Number();

  //! Builds the integer number \p i.
  Number(int i);

  //! Copy-constructor.
  Number(const Number& s);

  //! Destructor.
  ~Number();

  //! Assignment operator.
  Number& operator=(const Number& s);

  bool operator==(const Number &num) const;
  bool operator!=(const Number &num) const;
  bool operator>(const Number& num) const;
  bool operator<(const Number& num) const;
  bool operator>=(const Number& num) const;
  bool operator<=(const Number& num) const;

  //  relational operator==(const Symbol& lh, const Number& rh);

  bool is_positive() const;
  bool is_integer() const;
  bool is_pos_integer() const;
  bool is_nonneg_integer() const;
  bool is_even() const;
  bool is_odd() const;
  bool is_prime() const;
  bool is_rational() const;
  bool is_real() const;
  bool is_cinteger() const;
  bool is_crational() const;

  Number numer_denom() const;

  Number real() const;
  Number imag() const;
  Number numer() const;
  Number denom() const;

private:
  friend class Expr;

  GiNaC::numeric n;

public:
  // Made public only to initialize the constant I.
  Number(const GiNaC::numeric& gn);
};

extern const Number I;

int to_int(const Number& n);
long to_long(const Number& n);

Number abs(const Number& n);
Number irem(const Number& a, const Number& b);

} // namespace Parma_Recurrence_Relation_Solver

#include "Number.inlines.hh"

#endif // !defined(PURRS_Number_defs_hh)
