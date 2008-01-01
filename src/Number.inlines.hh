/* Number class implementation: inline functions.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#ifndef PURRS_Number_inlines_hh
#define PURRS_Number_inlines_hh

#include "Expr.defs.hh"
#include <climits>
#include <cassert>
#include <stdexcept>

namespace Parma_Recurrence_Relation_Solver {


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
Number::Number(const char* s)
  : n(s) {
}

inline
Number::Number(long n, long d)
  : n(n, d) {
}

inline
Number::Number(const Number& y)
  : n(y.n) {
}

inline
Number::Number(const GiNaC::numeric& gn)
  : n(gn) {
}

inline
Number::~Number() {
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

inline Number&
Number::operator=(const char* s) {
  n = s;
  return *this;
}

inline Number&
operator+=(Number& x, const Number& y) {
  x.n += y.n;
  return x;
}

inline Number&
operator-=(Number& x, const Number& y) {
  x.n -= y.n;
  return x;
}

inline Number&
operator*=(Number& x, const Number& y) {
  x.n *= y.n;
  return x;
}

inline Number&
operator/=(Number& x, const Number& y) {
  x.n /= y.n;
  return x;
}

inline Number
operator+(const Number& x) {
  return x;
}

inline Number
operator-(const Number& x) {
  return -x.n;
}

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
operator==(const Number& x, long i) {
  return x.n.operator==(i);
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

inline Number
abs(const Number& x) {
  return GiNaC::abs(x.n);
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
exact_pwr(const Number& x, const Number& y) {
  assert(x.is_complex_rational() && y.is_integer());
  // The following assertion is here for documentation purposes.
  // The current version of GiNaC already checks for these cases,
  // possibly emitting an exception.
  assert((!x.is_zero() || !y.is_zero())
	 && ((!x.is_zero()) || (x.is_zero() && y.is_positive())));
  return GiNaC::pow(x.n, y.n);
}

inline bool
Number::is_zero() const {
  return n.is_zero();
}

inline bool
Number::is_positive() const {
  return n.is_positive();
}

inline bool
Number::is_negative() const {
  return n.is_negative();
}

inline bool
Number::is_integer() const {
  return n.is_integer();
}

inline bool
Number::is_positive_integer() const {
  return n.is_pos_integer();
}

inline bool
Number::is_nonnegative_integer() const {
  return n.is_nonneg_integer();
}

inline bool
Number::is_even() const {
  return n.is_even();
}

inline bool
Number::is_odd() const {
  return n.is_odd();
}

inline bool
Number::is_prime() const {
  // FIXME: this test is due to a problem in CLN 1.1.5 and/or GiNaC 1.0.11.
  // See http://www.cs.unipr.it/pipermail/purrs-devel/2002-October/000411.html.
  if (!is_positive_integer())
    return false;
  return n.is_prime();
}

inline bool
Number::is_rational() const {
  return n.is_rational();
}

inline bool
Number::is_real() const {
  return n.is_real();
}

inline bool
Number::is_complex_integer() const {
  return n.is_cinteger();
}

inline bool
Number::is_complex_rational() const {
  return n.is_crational();
}

inline int
Number::to_int() const {
  // FIXME: this test is necessary to circumvent a misfeature
  // of CLN 1.1.5 and/or GiNaC 1.0.11.
  // See http://www.cs.unipr.it/pipermail/purrs-devel/2002-October/000412.html.
  if (!is_integer() || (n > 0 && n > INT_MAX)  || (n < 0 && n < INT_MIN))
    throw std::domain_error("Cannot convert to an `int' "
			    "in PURRS::Number::to_int()");
  return n.to_int();
}

inline unsigned int
Number::to_unsigned_int() const {
  // FIXME: this test is necessary to circumvent a misfeature
  // of CLN 1.1.5 and/or GiNaC 1.0.11.
  // See http://www.cs.unipr.it/pipermail/purrs-devel/2002-October/000412.html.
  if (!is_integer() || (n > 0 && n > INT_MAX)  || (n < 0 && n < INT_MIN))
    throw std::domain_error("Cannot convert to an `int' "
			    "in PURRS::Number::to_int()");
  return unsigned(n.to_int());
}

inline long
Number::to_long() const {
  // FIXME: this test is necessary to circumvent a misfeature
  // of CLN 1.1.5 and/or GiNaC 1.0.11.
  // See http://www.cs.unipr.it/pipermail/purrs-devel/2002-October/000412.html.
  if (!is_integer() || (n > 0 && n > LONG_MAX)  || (n < 0 && n < LONG_MIN))
    throw std::domain_error("Cannot convert to a `long' "
			    "in PURRS::Number::to_long()");
  return n.to_long();
}

inline double
Number::to_double() const {
  return n.to_double();
}

inline Number
Number::real() const {
  return n.real();
}

inline Number
Number::imaginary() const {
  return n.imag();
}

inline Number
Number::numerator() const {
  return n.numer();
}

inline Number
Number::denominator() const {
  return n.denom();
}

inline Number
mod(const Number& x, const Number& y) {
  // FIXME: this test is due to a problem in CLN 1.1.5 and/or GiNaC 1.0.11.
  // See http://www.cs.unipr.it/pipermail/purrs-devel/2002-October/000417.html.
  if (y.is_zero())
    throw std::runtime_error("Division by zero in PURRS::Number::mod()");
  return GiNaC::mod(x.n, y.n);
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Number_inlines_hh)
