/* Boundary class declaration.
   Copyright (C) 2002 Roberto Bagnara <bagnara@cs.unipr.it>

This file is part of the Parma Interval Library (PIL).

The PIL is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The PIL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the Parma Interval Library
site: http://www.cs.unipr.it/pil/ . */

#ifndef _Boundary_defs_hh
#define _Boundary_defs_hh 1

#include <iosfwd>
#include <cmath>
#include <cstddef>

using std::istream;
using std::ostream;

#include "round.h"
#include "IEEE754.h"

#include "Boundary.types.hh"

#if BOUNDARY_TYPE_IS_DOUBLE
#define BOUNDARY_MIN NINF_D
#define BOUNDARY_MAX PINF_D
#define BOUNDARY_MAX_ODD MAX_ODD_D
#define BOUNDARY_MAX_FRACT MAX_FRACT_D
#else
#define BOUNDARY_MIN NINF_F
#define BOUNDARY_MAX PINF_F
#define BOUNDARY_MAX_ODD MAX_ODD_F
#define BOUNDARY_MAX_FRACT MAX_FRACT_F
#endif

class Boundary {
public:
#if BOUNDARY_TYPE_IS_DOUBLE
  typedef double Value;
#else
  typedef float Value;
#endif

protected:
  // Generic flag for generic boundary (useful for compare)
  enum Flag { NEG=-1, ZERO=0, POS=1 };

public:

  Value boundary() const;
  bool is_closed() const;
  void set_min();
  void set_max();
  bool is_min() const;
  bool is_max() const;
  bool is_integer() const;
  static bool is_integer(Value v);

protected:
  Value value;
  Flag flag;

  Boundary() { }
  Boundary(Value value, Flag flag);

  friend int cmp (const Boundary& x, const Boundary& y);

  friend bool operator == (const Boundary& x, const Boundary& y);
  friend bool operator != (const Boundary& x, const Boundary& y);
  friend bool operator <  (const Boundary& x, const Boundary& y);
  friend bool operator >  (const Boundary& x, const Boundary& y);
  friend bool operator <= (const Boundary& x, const Boundary& y);
  friend bool operator >= (const Boundary& x, const Boundary& y);

  friend bool operator == (const Boundary& x, Value y);
  friend bool operator != (const Boundary& x, Value y);
  friend bool operator <  (const Boundary& x, Value y);
  friend bool operator >  (const Boundary& x, Value y);
  friend bool operator <= (const Boundary& x, Value y);
  friend bool operator >= (const Boundary& x, Value y);
  friend bool operator < (const LBoundary& x, const UBoundary& y);
  friend bool operator > (const UBoundary& x, const LBoundary& y);
  friend bool operator <= (const UBoundary& x, const LBoundary& y);
  friend bool operator >= (const LBoundary& x, const UBoundary& y);
  friend bool operator < (const LBoundary& x, Value y);
  friend bool operator > (const UBoundary& x, Value y);
  friend bool operator <= (const UBoundary& x, Value y);
  friend bool operator >= (const LBoundary& x, Value y);
  friend class LBoundary;
  friend class UBoundary;
};

class LBoundary : public Boundary {
public:
  // Closed or open flag
  enum Open_Closed { CLOSED=Boundary::ZERO, OPEN=Boundary::POS };
  void set_openclosed(Open_Closed flag);
  LBoundary() { }
  LBoundary(Boundary::Value value, Open_Closed flag=CLOSED);
#if !BOUNDARY_TYPE_IS_DOUBLE
  LBoundary(double value, Open_Closed flag=CLOSED);
#endif
  LBoundary(const char* value, Open_Closed flag=CLOSED);

  void set(const UBoundary& b);
  void set(Boundary::Value b);
  void set(const LBoundary& b, bool integer);
  void set_gt(Boundary::Value b, bool integer);

  LBoundary& floor_assign();
  LBoundary& ceil_assign();
  LBoundary& integer_assign();
  LBoundary& operator += (Boundary::Value y);
  LBoundary& operator += (const Boundary& y);
  LBoundary& operator -= (Boundary::Value y);
  LBoundary& operator -= (const Boundary& y);
  LBoundary& operator *= (Boundary::Value y);
  LBoundary& operator *= (const Boundary& y);
  LBoundary& operator /= (Boundary::Value y);
  LBoundary& operator /= (const Boundary& y);

  LBoundary& pow_assign (Boundary::Value y);
  LBoundary& pow_assign (const Boundary& y);

  LBoundary& negate();
  LBoundary& exp_assign();
  LBoundary& log_assign();
  LBoundary& sqrt_assign();
  LBoundary& sin_assign();
  LBoundary& cos_assign();
  LBoundary& tan_assign();

private:
  void round_set() const;
  friend ostream& operator << (ostream& s, const LBoundary& b);
  friend istream& operator >> (istream& s, LBoundary& b);
  friend ostream& operator << (ostream& s, const class Interval& x);

  friend class UBoundary;
};

class UBoundary : public Boundary {
public:
  // Closed or open flag
  enum Open_Closed { CLOSED=Boundary::ZERO, OPEN=Boundary::NEG };
  void set_openclosed(Open_Closed flag);
  UBoundary() { }
  UBoundary(Boundary::Value value, Open_Closed flag=CLOSED);
#if !BOUNDARY_TYPE_IS_DOUBLE
  UBoundary(double value, Open_Closed flag=CLOSED);
#endif
  UBoundary(const char* value, Open_Closed flag=CLOSED);

  void set(const LBoundary& b);
  void set(Boundary::Value b);
  void set(const UBoundary& b, bool integer);
  void set_lt(Boundary::Value b, bool integer);

  UBoundary& floor_assign();
  UBoundary& ceil_assign();
  UBoundary& integer_assign();
  UBoundary& operator += (Boundary::Value y);
  UBoundary& operator += (const Boundary& y);
  UBoundary& operator -= (Boundary::Value y);
  UBoundary& operator -= (const Boundary& y);
  UBoundary& operator *= (Boundary::Value y);
  UBoundary& operator *= (const Boundary& y);
  UBoundary& operator /= (Boundary::Value y);
  UBoundary& operator /= (const Boundary& y);

  UBoundary& pow_assign (Boundary::Value y);
  UBoundary& pow_assign (const Boundary& y);

  UBoundary& negate();
  UBoundary& exp_assign();
  UBoundary& log_assign();
  UBoundary& sqrt_assign();
  UBoundary& sin_assign();
  UBoundary& cos_assign();
  UBoundary& tan_assign();

private:
  void round_set() const;
  friend ostream& operator << (ostream& s, const UBoundary& b);
  friend istream& operator >> (istream& s, UBoundary& b);

  friend class LBoundary;
};

#if 0
LBoundary operator - (const LBoundary& x);
LBoundary floor (const LBoundary& x);
LBoundary ceil (const LBoundary& x);
LBoundary integer (const LBoundary& x);
LBoundary exp (const LBoundary& x);
LBoundary log (const LBoundary& x);
LBoundary sqrt (const LBoundary& x);
LBoundary sin (const LBoundary& x);
LBoundary cos (const LBoundary& x);
LBoundary tan (const LBoundary& x);
UBoundary operator - (const UBoundary& x);
UBoundary floor (const UBoundary& x);
UBoundary ceil (const UBoundary& x);
UBoundary integer (const UBoundary& x);
UBoundary exp (const UBoundary& x);
UBoundary log (const UBoundary& x);
UBoundary sqrt (const UBoundary& x);
UBoundary sin (const UBoundary& x);
UBoundary cos (const UBoundary& x);
UBoundary tan (const UBoundary& x);
LBoundary operator + (const LBoundary& x, Boundary::Value y);
LBoundary operator + (const LBoundary& x, const Boundary& y);
LBoundary operator - (const LBoundary& x, Boundary::Value y);
LBoundary operator - (const LBoundary& x, const Boundary& y);
LBoundary operator * (const LBoundary& x, Boundary::Value y);
LBoundary operator * (const LBoundary& x, const Boundary& y);
LBoundary operator / (const LBoundary& x, Boundary::Value y);
LBoundary operator / (const LBoundary& x, const Boundary& y);
LBoundary pow (const LBoundary& x, Boundary::Value y);
LBoundary pow (const LBoundary& x, const Boundary& y);
UBoundary operator + (const UBoundary& x, Boundary::Value y);
UBoundary operator + (const UBoundary& x, const Boundary& y);
UBoundary operator - (const UBoundary& x, Boundary::Value y);
UBoundary operator - (const UBoundary& x, const Boundary& y);
UBoundary operator * (const UBoundary& x, Boundary::Value y);
UBoundary operator * (const UBoundary& x, const Boundary& y);
UBoundary operator / (const UBoundary& x, Boundary::Value y);
UBoundary operator / (const UBoundary& x, const Boundary& y);
UBoundary pow (const UBoundary& x, Boundary::Value y);
UBoundary pow (const UBoundary& x, const Boundary& y);
#endif

#if !OUTLINE
#include "Boundary.inlines.hh"
#endif

#endif
