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

//! Binary arithmetic operators Number with Number.
Number operator+(const Number& lh, const Number& rh);
Number operator-(const Number& lh, const Number& rh);
Number operator*(const Number& lh, const Number& rh);
Number operator/(const Number& lh, const Number& rh);

//! Unary operators.
Number operator+(const Number &lh);
Number operator-(const Number &lh);

//! Increment / decrement operators.
Number& operator++(Number& rh);
Number& operator--(Number& rh);
Number operator++(Number& lh, int);
Number operator--(Number& lh, int);

bool operator==(long x, const Number& y);
bool operator!=(long x, const Number& y);
bool operator>(long x, const Number& y);
bool operator<(long x, const Number& y);
bool operator>=(long x, const Number& y);
bool operator<=(long x, const Number& y);

Number power(const Number& b, const Number& e);
Number lcm(const Number& x, const Number& y);
Number factorial(const Number& n);

class Number {
public:
  //! Default constructor.
  Number();

  //! Builds the integer number \p i.
  Number(int i);

  //! Builds the integer number \p i.
  Number(unsigned int i);

  //! Builds the integer number \p i.
  Number(long i);

  //! Builds the integer number \p i.
  Number(unsigned long i);

  //! Builds the number \p numer / denom.
  Number(long numer, long denom);

  //! Copy-constructor.
  Number(const Number& s);

  //! Destructor.
  ~Number();

  //! Assignment operator.
  Number& operator=(const Number& s);

  bool operator==(const Number& num) const;
  bool operator!=(const Number& num) const;
  bool operator>(const Number& num) const;
  bool operator<(const Number& num) const;
  bool operator>=(const Number& num) const;
  bool operator<=(const Number& num) const;

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

  int to_int() const;
  long to_long() const;
  Number real() const;
  Number imag() const;
  Number numer() const;
  Number denom() const;

private:
  friend class Expr;
  friend Number power(const Number& b, const Number& e);
  friend Number lcm(const Number& x, const Number& y);
  friend Number factorial(const Number& n);

  GiNaC::numeric n;

public:
  // Made public only to initialize the constant I.
  Number(const GiNaC::numeric& gn);
};

extern const Number I;

Number gcd(const Number& a, const Number& b);
Number abs(const Number& n);
Number irem(const Number& a, const Number& b);

} // namespace Parma_Recurrence_Relation_Solver

#include "Number.inlines.hh"

#endif // !defined(PURRS_Number_defs_hh)
