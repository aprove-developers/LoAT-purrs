/* Number class declaration.
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

//! Returns \f$ x \f$.
Number operator+(const Number& x);

//! Returns \f$ -x \f$.
Number operator-(const Number& x);

//! Returns \f$ x + y \f$.
Number operator+(const Number& x, const Number& y);

//! Returns \f$ x - y \f$.
Number operator-(const Number& x, const Number& y);

//! Returns \f$ x \cdot y \f$.
Number operator*(const Number& x, const Number& y);

//! If \f$ y \neq 0 \f$, returns \f$ x / y \f$.
/*!
  \exception std::runtime_error thrown if \f$ y = 0 \f$.
*/
Number operator/(const Number& x, const Number& y);

//! Assigns \f$ x + y \f$ to \f$ x \f$ and returns the result.
Number& operator+=(Number& x, const Number& y);

//! Assigns \f$ x - y \f$ to \f$ x \f$ and returns the result.
Number& operator-=(Number& x, const Number& y);

//! Assigns \f$ x \cdot y \f$ to \f$ x \f$ and returns the result.
Number& operator*=(Number& x, const Number& y);

//! If \f$ y \neq 0 \f$, assigns \f$ x / y \f$ to \f$ x \f$
//! and returns the result.
/*!
  \exception std::runtime_error thrown if \f$ y = 0 \f$.
*/
Number& operator/=(Number& x, const Number& y);

//! Returns <CODE>true</CODE> if and only if \f$ x = y \f$.
bool operator==(const Number& x, const Number& y);

//! Returns <CODE>true</CODE> if and only if \f$ x \neq y \f$.
bool operator!=(const Number& x, const Number& y);

//! Returns <CODE>true</CODE> if and only if \f$ x > y \f$.
bool operator>(const Number& x, const Number& y);

//! Returns <CODE>true</CODE> if and only if \f$ x < y \f$.
bool operator<(const Number& x, const Number& y);

//! Returns <CODE>true</CODE> if and only if \f$ x \geq y \f$.
bool operator>=(const Number& x, const Number& y);

//! Returns <CODE>true</CODE> if and only if \f$ x \leq y \f$.
bool operator<=(const Number& x, const Number& y);

//! Returns the absolute value of \f$ x \f$.
Number abs(const Number& x);

//! Returns the least common multiple of \f$ x \f$ and \f$ y \f$.
Number factorial(const Number& n);

//! Returns the greatest common divisor of \f$ x \f$ and \f$ y \f$.
Number gcd(const Number& x, const Number& y);

//! Returns the least common multiple of \f$ x \f$ and \f$ y \f$.
Number lcm(const Number& x, const Number& y);

//! If ..., returns \f$ x^y \f$.
Number power(const Number& x, const Number& y);


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

  //! Builds the integer corresponding to the decimal integer numeral in \p s.
  Number(const char* s);

  //! If \f$ d \neq 0 \f$, builds the rational number \f$ n/d \f$.
  /*!
    \exception std::runtime_error thrown if \f$ d = 0 \f$.
  */
  Number(long n, long d);

  //! Copy-constructor.
  Number(const Number& s);

  //! Destructor.
  ~Number();

  //! Ordinary assignment operator.
  Number& operator=(const Number& y);

  //! Assigns \p i to \p *this.
  Number& operator=(int i);

  //! Assigns \p i to \p *this.
  Number& operator=(unsigned int i);

  //! Assigns \p i to \p *this.
  Number& operator=(long i);

  //! Assigns \p i to \p *this.
  Number& operator=(unsigned long i);

  //! Assignes to \p *this the integer corresponding
  //! to the decimal integer numeral in \p s.
  Number& operator=(const char* s);

  //! Pre-increment operator.
  Number& operator++();

  //! Pre-decrement operator.
  Number& operator--();

  //! Post-increment operator.
  Number operator++(int);

  //! Post-decrement operator.
  Number operator--(int);

  //! Returns <CODE>true</CODE> if and only if \p *this is positive.
  bool is_positive() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is an integer.
  bool is_integer() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is a positive integer.
  bool is_positive_integer() const;

  //! Returns <CODE>true</CODE> if and only if \p *this
  //! is a nonnegative integer.
  bool is_nonnegative_integer() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is an even integer.
  bool is_even() const;

  bool is_odd() const;
  bool is_prime() const;
  bool is_rational() const;
  bool is_real() const;
  bool is_complex_integer() const;

  int to_int() const;
  long to_long() const;
  Number real() const;
  Number imaginary() const;
  Number numerator() const;
  Number denominator() const;

private:
  friend class Expr;

  friend Number operator+(const Number& x);
  friend Number operator-(const Number& x);

  friend Number operator+(const Number& x, const Number& y);
  friend Number operator-(const Number& x, const Number& y);
  friend Number operator*(const Number& x, const Number& y);
  friend Number operator/(const Number& x, const Number& y);

  friend Number& operator+=(Number& x, const Number& y);
  friend Number& operator-=(Number& x, const Number& y);
  friend Number& operator*=(Number& x, const Number& y);
  friend Number& operator/=(Number& x, const Number& y);

  friend bool operator==(const Number& x, const Number& y);
  friend bool operator!=(const Number& x, const Number& y);
  friend bool operator>(const Number& x, const Number& y);
  friend bool operator<(const Number& x, const Number& y);
  friend bool operator>=(const Number& x, const Number& y);
  friend bool operator<=(const Number& x, const Number& y);

  friend Number abs(const Number& x);
  friend Number factorial(const Number& x);
  friend Number gcd(const Number& x, const Number& y);
  friend Number lcm(const Number& x, const Number& y);
  friend Number power(const Number& x, const Number& y);

private:
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
