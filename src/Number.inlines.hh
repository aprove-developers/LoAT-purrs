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

inline Number
operator+(const Number& x, const Number& y) {
  return x.n + y.n;
}

inline Number
operator-(const Number& x, const Number& y) {
  return x.n - y.n;
}

inline Number
operator*(const Number& x, const Number& y) {
  return x.n * y.n;
}

inline Number
operator/(const Number& x, const Number& y) {
  return x.n / y.n;
}

inline Number
operator+(const Number& x) {
  return x;
}

inline Number
operator-(const Number& x) {
  return -x.n;
}

inline Number&
Number::operator++() {
  ++n;
  return *this;
}

inline Number&
Number::operator--() {
  --n;
  return *this;
}

inline Number
Number::operator++(int) {
  return n++;
}

inline Number
Number::operator--(int) {
  return n--;
}

inline bool
operator==(const Number& x, const Number& y) {
  return x.n == y.n;
}

inline bool
operator!=(const Number& x, const Number& y) {
  return x.n != y.n;
}

inline bool
operator>(const Number& x, const Number& y) {
  return x.n > y.n;
}

inline bool
operator<(const Number& x, const Number& y) {
  return x.n < y.n;
}

inline bool
operator>=(const Number& x, const Number& y) {
  return x.n >= y.n;
}

inline bool
operator<=(const Number& x, const Number& y) {
  return x.n <= y.n;
}

inline
Number::Number() {
}

inline
Number::Number(int i)
  : n(i) {
}

inline
Number::Number(unsigned int i)
  : n(i) {
}

inline
Number::Number(long i)
  : n(i) {
}

inline
Number::Number(unsigned long i)
  : n(i) {
}

inline
Number::Number(long n, long d)
  : n(n, d) {
}

inline
Number::Number(const Number& num)
  : n(num.n) {
}

inline Number&
Number::operator=(const Number& y) {
  n = y.n;
  return *this;
}

inline Number&
Number::operator=(int i) {
  n = i;
  return *this;
}

inline Number&
Number::operator=(unsigned int i) {
  n = i;
  return *this;
}

inline Number&
Number::operator=(long i) {
  n = i;
  return *this;
}

inline Number&
Number::operator=(unsigned long i) {
  n = i;
  return *this;
}

inline
Number::Number(const GiNaC::numeric& gn)
  : n(gn) {
}

inline
Number::~Number() {
}

inline bool
Number::Number::is_positive() const {
  return n.is_positive();
}

inline bool
Number::Number::is_integer() const {
  return n.is_integer();
}

inline bool
Number::Number::is_pos_integer() const {
  return n.is_pos_integer();
}

inline bool
Number::Number::is_nonneg_integer() const {
  return n.is_nonneg_integer();
}

inline bool
Number::Number::is_even() const {
  return n.is_even();
}

inline bool
Number::Number::is_odd() const {
  return n.is_odd();
}

inline bool
Number::Number::is_prime() const {
  return n.is_prime();
}

inline bool
Number::Number::is_rational() const {
  return n.is_rational();
}

inline bool
Number::Number::is_real() const {
  return n.is_real();
}

inline bool
Number::Number::is_cinteger() const {
  return n.is_cinteger();
}

inline bool
Number::Number::is_crational() const {
  return n.is_crational();
}

//  inline Expr
//  numer_denom() const {
//    return n.numer_denom();
//  }

inline int
Number::to_int() const {
  return n.to_int();
}

inline long
Number::to_long() const {
  return n.to_long();
}

inline Number
Number::real() const {
  return n.real();
}

inline Number
Number::imag() const {
  return n.imag();
}

inline Number
Number::numer() const {
  return n.numer();
}

inline Number
Number::denom() const {
  return n.denom();
}

inline Number
irem(const Number& a, const Number& b) {
  return irem(a, b);
}

inline Number
abs(const Number& n) {
  return abs(n);
}

inline Number
power(const Number& x, const Number& y) {
  return GiNaC::pow(x.n, y.n);
}

inline Number
gcd(const Number& x, const Number& y) {
  return GiNaC::gcd(x.n, y.n);
}

inline Number
lcm(const Number& x, const Number& y) {
  return GiNaC::lcm(x.n, y.n);
}

inline Number
factorial(const Number& n) {
  return GiNaC::factorial(n.n);
};

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Number_inlines_hh)
