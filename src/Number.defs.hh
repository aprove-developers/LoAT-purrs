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

//! \brief
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

//! Returns <CODE>true</CODE> if and only if \f$ x = i \f$.
bool operator==(const Number& x, long i);

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

//! If \f$ x \f$ is a natural number, returns \f$ x! \f$.
/*!
  \exception std::range_error thrown if \f$ x \f$ is not a natural number.
*/
Number factorial(const Number& x);

//! If \f$ x \f$ and \f$ y \f$ are integer, returns the greatest common
//! divisor of \f$ x \f$ and \f$ y \f$.
//! Returns 1 otherwise.
Number gcd(const Number& x, const Number& y);

//! If \f$ x \f$ and \f$ y \f$ are integer, returns the least common
//! multiple of \f$ x \f$ and \f$ y \f$.
//! Returns \f$ x \cdot y \f$ otherwise.
Number lcm(const Number& x, const Number& y);

//! If \f$ x \f$ and \f$ y \f$ are not zero or \f$ x = 0 \f$ and \f$ y \f$
//! is a positive rational number, returns \f$ x^y \f$.
/*!
  \exception std::runtime_error thrown if \f$ x = y = 0 \f$.
  \exception std::logic_error   thrown if \f$ x = 0 \f$ and \f$ y \f$
                                is not a positive rational number.
*/
Number pwr(const Number& x, const Number& y);

//! If \f$ x \f$ and \f$ y \f$ are integer and \f$ y \neq 0 \f$,
//! returns the remainder of the division of \f$ x \f$ by \f$ y \f$.
//! If \f$ y = 0 \f$ an exception is thrown.
//! Returns zero otherwise.
//! If the result is not zero, its sign is the sign of \f$ y \f$.
/*!
  \exception std::runtime_error thrown if \f$ y = 0 \f$.
*/
Number mod(const Number& x, const Number& y);

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
  Number(const Number& y);

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

  //! Returns <CODE>true</CODE> if and only if \p *this is zero.
  bool is_zero() const;

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

  //! Returns <CODE>true</CODE> if and only if \p *this is an odd integer.
  bool is_odd() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is a prime natural.
  bool is_prime() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is a rational number.
  bool is_rational() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is a real number.
  bool is_real() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is a complex number
  //! with integral real and imaginary parts.
  bool is_complex_integer() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is a complex number
  //! with rational real and imaginary parts.
  bool is_complex_rational() const;

  //! Returns the <CODE>int</CODE> corresponding to \p *this, if any.
  /*!
    \exception std::domain_error thrown if \p *this is not convertible
                                 to <CODE>int</CODE>.
  */
  int to_int() const;

  //! Returns the <CODE>long</CODE> corresponding to \p *this, if any.
  /*!
    \exception std::domain_error thrown if \p *this is not convertible
                                 to <CODE>long</CODE>.
  */
  long to_long() const;

  //! Returns the real part of \p *this, seen as a complex number.
  Number real() const;

  //! Returns the imaginary part of \p *this, seen as a complex number.
  Number imaginary() const;

  //! Returns the numerator of \p *this, seen as a fraction.
  Number numerator() const;

  //! Returns the denominator of \p *this, seen as a fraction.
  Number denominator() const;

  static const Number I;

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

  friend bool operator==(const Number& x, long i);

  friend Number abs(const Number& x);
  friend Number factorial(const Number& x);
  friend Number gcd(const Number& x, const Number& y);
  friend Number lcm(const Number& x, const Number& y);
  friend Number pwr(const Number& x, const Number& y);
  friend Number mod(const Number& x, const Number& y);


private:
  friend bool operator==(const Expr& e, const Number& n);

  GiNaC::numeric n;

  //! Builds the number corresponding to \p gn.
  Number(const GiNaC::numeric& gn);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Number.inlines.hh"

#endif // !defined(PURRS_Number_defs_hh)
