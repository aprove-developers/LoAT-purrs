/* Boundary class implementation: inline functions.
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
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.

For the most up-to-date information see the Parma Interval Library
site: http://www.cs.unipr.it/pil/ . */

#include <cmath>
#include <cassert>

using std::sqrt;
using std::sin;
using std::cos;

extern double my_strtod(const char *nptr, char **endptr);
#ifndef POW_BUG_FIXED
extern double my_pow(double a, double b);
#define pow(a,b) my_pow(a,b)
#endif

// Simple
#define COMPUTE_SIMPLE(stmt) \
do { \
  stmt; \
} while(0)

// Directed rounding
#define COMPUTE_DR(stmt) \
do { \
  round_set(); \
  stmt; \
} while(0)

// Directed rounding and inexact result detection
#define COMPUTE_DR_IE(stmt) \
do { \
  round_set(); \
  if (is_closed()) { \
    reset_inexact(); \
    stmt; \
    if (is_inexact()) \
      set_openclosed(OPEN); \
  } \
  else { \
    stmt; \
  } \
} while(0)

#if USE_DIRECTED_ROUNDING
#if CHECK_INEXACT_RESULT
#if IE_ASSIGN
#define COMPUTE_ASSIGN(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_ASSIGN(stmt) COMPUTE_DR(stmt)
#endif
#if IE_STRTOD
#define COMPUTE_STRTOD(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_STRTOD(stmt) COMPUTE_DR(stmt)
#endif
#if IE_ADD
#define COMPUTE_ADD(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_ADD(stmt) COMPUTE_DR(stmt)
#endif
#if IE_SUB
#define COMPUTE_SUB(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_SUB(stmt) COMPUTE_DR(stmt)
#endif
#if IE_MUL
#define COMPUTE_MUL(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_MUL(stmt) COMPUTE_DR(stmt)
#endif
#if IE_DIV
#define COMPUTE_DIV(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_DIV(stmt) COMPUTE_DR(stmt)
#endif
#if IE_EXP
#define COMPUTE_EXP(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_EXP(stmt) COMPUTE_DR(stmt)
#endif
#if IE_LOG
#define COMPUTE_LOG(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_LOG(stmt) COMPUTE_DR(stmt)
#endif
#if IE_SQRT
#define COMPUTE_SQRT(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_SQRT(stmt) COMPUTE_DR(stmt)
#endif
#if IE_SIN
#define COMPUTE_SIN(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_SIN(stmt) COMPUTE_DR(stmt)
#endif
#if IE_COS
#define COMPUTE_COS(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_COS(stmt) COMPUTE_DR(stmt)
#endif
#if IE_TAN
#define COMPUTE_TAN(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_TAN(stmt) COMPUTE_DR(stmt)
#endif
#if IE_POW
#define COMPUTE_POW(stmt) COMPUTE_DR_IE(stmt)
#else
#define COMPUTE_POW(stmt) COMPUTE_DR(stmt)
#endif

#else

#define COMPUTE_ASSIGN(stmt) COMPUTE_DR(stmt)
#define COMPUTE_STRTOD(stmt) COMPUTE_DR(stmt)
#define COMPUTE_ADD(stmt) COMPUTE_DR(stmt)
#define COMPUTE_SUB(stmt) COMPUTE_DR(stmt)
#define COMPUTE_MUL(stmt) COMPUTE_DR(stmt)
#define COMPUTE_DIV(stmt) COMPUTE_DR(stmt)
#define COMPUTE_EXP(stmt) COMPUTE_DR(stmt)
#define COMPUTE_LOG(stmt) COMPUTE_DR(stmt)
#define COMPUTE_SQRT(stmt) COMPUTE_DR(stmt)
#define COMPUTE_SIN(stmt) COMPUTE_DR(stmt)
#define COMPUTE_COS(stmt) COMPUTE_DR(stmt)
#define COMPUTE_TAN(stmt) COMPUTE_DR(stmt)
#define COMPUTE_POW(stmt) COMPUTE_DR(stmt)
#endif

#else

#define COMPUTE_ASSIGN(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_STRTOD(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_ADD(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_SUB(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_MUL(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_DIV(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_EXP(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_LOG(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_SQRT(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_SIN(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_COS(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_TAN(stmt) COMPUTE_SIMPLE(stmt)
#define COMPUTE_POW(stmt) COMPUTE_SIMPLE(stmt)
#endif

// Complete constructor (protected)
INLINE
Boundary::Boundary(Value v, Flag f)
: value(v),
  flag(f)
{ }

INLINE Boundary::Value
Boundary::boundary() const {
  return value;
}

// Set boundary to minimum possible
INLINE void
Boundary::set_min() {
  value = BOUNDARY_MIN;
#if IEEE_NUMBERS
  flag = ZERO
#else
  flag = POS;
#endif
}

// Set boundary to maximum possible
INLINE void
Boundary::set_max() {
  value = BOUNDARY_MAX;
#if IEEE_NUMBERS
  flag = ZERO
#else
  flag = NEG;
#endif
}

INLINE bool
Boundary::is_min() const {
  return value == BOUNDARY_MIN;
}

INLINE bool
Boundary::is_max() const {
  return value == BOUNDARY_MAX;
}

INLINE Boundary::Value my_floor(Boundary::Value v) {
  if (v < -BOUNDARY_MAX_FRACT || v > BOUNDARY_MAX_FRACT)
    return v;
  return floor(v);
}

INLINE Boundary::Value my_ceil(Boundary::Value v) {
  if (v < -BOUNDARY_MAX_FRACT || v > BOUNDARY_MAX_FRACT)
    return v;
  return ceil(v);
}

INLINE bool
Boundary::is_integer(Value v) {
  return v != BOUNDARY_MIN && v != BOUNDARY_MAX &&
    (v < -BOUNDARY_MAX_FRACT || v > BOUNDARY_MAX_FRACT || v == floor(v));
}

// Return true if the Boundary have an integer value
INLINE bool
Boundary::is_integer() const {
  return is_integer(value);
}

// Return the open or closed flag appropriate
INLINE bool
Boundary::is_closed() const {
  return flag == ZERO;
}

// Set boundary open or closed
INLINE void
LBoundary::set_openclosed(LBoundary::Open_Closed f) {
  flag = (f == CLOSED?ZERO:POS);
}

// Set the correct round mode for this type of boundary
INLINE void
LBoundary::round_set() const {
  round_down();
}

// Complete constructor
INLINE
LBoundary::LBoundary(Boundary::Value value, LBoundary::Open_Closed f)
: Boundary(value, f == CLOSED ? ZERO : POS) {
#if !IEEE_NUMBERS
  assert(f == OPEN || value != BOUNDARY_MIN);
#endif
}

#if !BOUNDARY_TYPE_IS_DOUBLE
INLINE
LBoundary::LBoundary(double val, LBoundary::Open_Closed f) {
  flag = (f == CLOSED ? ZERO : POS);
  round_save();
  COMPUTE_ASSIGN((value = val));
  round_restore();
#if !IEEE_NUMBERS
  assert(flag == POS || value != BOUNDARY_MIN);
#endif
}
#endif

// Complete constructor with string value
INLINE
LBoundary::LBoundary(const char* val, LBoundary::Open_Closed f) {
  flag = (f == CLOSED?ZERO:POS);
  round_save();
  COMPUTE_STRTOD((value = my_strtod(val, 0)));
  round_restore();
#if !IEEE_NUMBERS
  assert(f == OPEN || value != BOUNDARY_MIN);
#endif
}

// Set boundary open or closed
INLINE void
UBoundary::set_openclosed(UBoundary::Open_Closed f) {
  flag = (f == CLOSED?ZERO:NEG);
}

// Round to smallest integer not less
INLINE LBoundary&
LBoundary::ceil_assign() {
  if (value != BOUNDARY_MIN && 
      value >= -BOUNDARY_MAX_FRACT && value <= BOUNDARY_MAX_FRACT) {
    if (is_closed())
      value = ceil(value);
    else {
      value = floor(value)+1.0f;
      set_openclosed(CLOSED);
    }
  }
  return *this;
}

// Round to largest integer not greater
INLINE LBoundary&
LBoundary::floor_assign() {
  if (value != BOUNDARY_MIN) {
    value = my_floor(value);
    set_openclosed(CLOSED);
  }
  return *this;
}

// Set a LBoundary equal to a UBoundary
INLINE void
LBoundary::set(const UBoundary& b) {
  value = b.value;
  set_openclosed(b.is_closed()?CLOSED:OPEN);
}

INLINE void
LBoundary::set(Boundary::Value b) {
  value = b;
  set_openclosed(CLOSED);
}

INLINE void
LBoundary::set(const LBoundary& b, bool integer) {
  value = b.value;
  flag = b.flag;
  if (integer)
    ceil_assign();
}

// Set this boundary an epsilon greater than b (if possible)
// The boundary is the lower one of an integer interval if specified
INLINE void
LBoundary::set_gt(Boundary::Value b, bool integer) {
  if (integer && b != BOUNDARY_MIN) {
    if (value < -BOUNDARY_MAX_ODD-1 || value > BOUNDARY_MAX_ODD)
      return;
    if (value == BOUNDARY_MAX_FRACT)
      value = ceil(value);
    else
      value = floor(b+1.0f);
    set_openclosed(CLOSED);
  }
  else {
    value = b;
    set_openclosed(OPEN);
  }
}

// Set the correct round mode for this type of boundary
INLINE void
UBoundary::round_set() const {
  round_up();
}

// Complete constructor
INLINE
UBoundary::UBoundary(Boundary::Value value, UBoundary::Open_Closed f)
: Boundary(value, f == CLOSED ? ZERO : NEG) {
#if !IEEE_NUMBERS
  assert(f == OPEN || value != BOUNDARY_MAX);
#endif
}

#if !BOUNDARY_TYPE_IS_DOUBLE
INLINE
UBoundary::UBoundary(double val, UBoundary::Open_Closed f) {
  flag = (f == CLOSED ? ZERO : NEG);
  round_save();
  COMPUTE_ASSIGN((value = val));
  round_restore();
#if !IEEE_NUMBERS
  assert(flag == NEG || value != BOUNDARY_MAX);
#endif
}
#endif

// Complete constructor with string value
INLINE
UBoundary::UBoundary(const char* val, UBoundary::Open_Closed f) {
  flag = (f == CLOSED?ZERO:NEG);
  round_save();
  COMPUTE_STRTOD((value = my_strtod(val, 0)));
  round_restore();
#if !IEEE_NUMBERS
  assert(f == OPEN || value != BOUNDARY_MAX);
#endif
}

// Round to largest integer not greater
INLINE UBoundary&
UBoundary::floor_assign() {
  if (value != BOUNDARY_MAX && 
      value >= -BOUNDARY_MAX_FRACT && value <= BOUNDARY_MAX_FRACT) {
    if (is_closed())
      value = floor(value);
    else {
      value = ceil(value)-1.0f;
      set_openclosed(CLOSED);
    }
  }
  return *this;
}

// Set a UBoundary equal to a LBoundary
INLINE void
UBoundary::set(const LBoundary& b) {
  value = b.value;
  set_openclosed(b.is_closed()?CLOSED:OPEN);
}

INLINE void
UBoundary::set(Boundary::Value b) {
  value = b;
  set_openclosed(CLOSED);
}

INLINE void
UBoundary::set(const UBoundary& b, bool integer) {
  value = b.value;
  flag = b.flag;
  if (integer)
    floor_assign();
}

// Set this boundary an epsilon less than b (if possible)
// The boundary is the upper one of an integer interval if specified
INLINE void
UBoundary::set_lt(Boundary::Value b, bool integer) {
  if (integer && b != BOUNDARY_MAX) {
    if (value < -BOUNDARY_MAX_ODD || value > BOUNDARY_MAX_ODD+1)
      return;
    if (value == -BOUNDARY_MAX_FRACT)
      value = floor(value);
    else
      value = my_ceil(b-1.0f);
    set_openclosed(CLOSED);
  }
  else {
    value = b;
    set_openclosed(OPEN);
  }
}

// Return a negative number if x < y, 0 if x = y, a positive number if x > y
INLINE int
cmp (const Boundary& x, const Boundary& y) {
  Boundary::Value d = x.value-y.value;
  if (d < 0.0)
    return -1;
  else if (d > 0.0)
    return 1;
  else
    return x.flag-y.flag;
}

// Return true if x = y
INLINE bool
operator == (const Boundary& x, const Boundary& y) {
  return x.value == y.value && x.flag == y.flag;
}

// Return true if x != y
INLINE bool
operator != (const Boundary& x, const Boundary& y) {
  return !(x == y);
}

// Return true if x < y
INLINE bool
operator < (const Boundary& x, const Boundary& y) {
  return x.value < y.value ||
    (x.value == y.value && x.flag < y.flag);
}

// Return true if x < y (specialized version)
INLINE bool
operator < (const LBoundary& x, const UBoundary& y) {
  return x.value < y.value;
}

// Return true if x > y
INLINE bool
operator > (const Boundary& x, const Boundary& y) {
  return x.value > y.value ||
    (x.value == y.value && x.flag > y.flag);
}

// Return true if x > y (specialized version)
INLINE bool
operator > (const UBoundary& x, const LBoundary& y) {
  return x.value > y.value;
}

// Return true if x <= y
INLINE bool
operator <= (const Boundary& x, const Boundary& y) {
  return !(x > y);
}

// Return true if x <= y (specialized version)
INLINE bool
operator <= (const UBoundary& x, const LBoundary& y) {
  return !(x > y);
}

// Return true if x >= y
INLINE bool
operator >= (const Boundary& x, const Boundary& y) {
  return !(x < y);
}

// Return true if x <= y (specialized version)
INLINE bool
operator >= (const LBoundary& x, const UBoundary& y) {
  return !(x < y);
}

// Return true if x == y
INLINE bool
operator == (const Boundary& x, Boundary::Value y) {
  return x.value == y && x.flag == Boundary::ZERO;
}

// Return true if x != y
INLINE bool
operator != (const Boundary& x, Boundary::Value y) {
  return !(x == y);
}

// Return true if x < y
INLINE bool
operator < (const Boundary& x, Boundary::Value y) {
  return x.value < y ||
    (x.value == y && x.flag < Boundary::ZERO);
}

// Return true if x < y (specialized version)
INLINE bool
operator < (const LBoundary& x, Boundary::Value y) {
  return x.value < y;
}

// Return true if x > y
INLINE bool
operator > (const Boundary& x, Boundary::Value y) {
  return x.value > y ||
    (x.value == y && x.flag > Boundary::ZERO);
}

// Return true if x > y (specialized version)
INLINE bool
operator > (const UBoundary& x, Boundary::Value y) {
  return x.value > y;
}

// Return true if x <= y
INLINE bool
operator <= (const Boundary& x, Boundary::Value y) {
  return !(x > y);
}

// Return true if x <= y (specialized version)
INLINE bool
operator <= (const UBoundary& x, Boundary::Value y) {
  return !(x > y);
}

// Return true if x >= y
INLINE bool
operator >= (const Boundary& x, Boundary::Value y) {
  return !(x < y);
}

// Return true if x >= y (specialized version)
INLINE bool
operator >= (const LBoundary& x, Boundary::Value y) {
  return !(x < y);
}

// Round to closest integer toward 0
INLINE LBoundary&
LBoundary::integer_assign() {
  if (*this >= 0.0f)
    floor_assign();
  else
    ceil_assign();
  return *this;
}

// Round to smallest integer not less
INLINE UBoundary&
UBoundary::ceil_assign() {
  if (value != BOUNDARY_MAX) {
    value = my_ceil(value);
    set_openclosed(CLOSED);
  }
  return *this;
}

// Round to closest integer toward 0
INLINE UBoundary&
UBoundary::integer_assign() {
  if (*this > 0.0f)
    floor_assign();
  else
    ceil_assign();
  return *this;
}

// Negate (unary minus)
// No need to set rounding direction
INLINE LBoundary&
LBoundary::negate() {
  value = -value;
  return *this;
}

INLINE UBoundary&
UBoundary::negate() {
  value = -value;
  return *this;
}

// Build arithmetic assign operator
INLINE LBoundary&
LBoundary::operator += (Boundary::Value x) {
  COMPUTE_ADD((value = value + x));
  return *this;
}

INLINE LBoundary&
LBoundary::operator += (const Boundary& x) {
  COMPUTE_ADD((value = value + x.value));
  if (is_closed())
    set_openclosed(x.is_closed()?CLOSED:OPEN);
  return *this;
}

INLINE LBoundary&
LBoundary::operator -= (Boundary::Value x) {
  COMPUTE_SUB((value = value - x));
  return *this;
}

INLINE LBoundary&
LBoundary::operator -= (const Boundary& x) {
  COMPUTE_SUB((value = value - x.value));
  if (is_closed())
    set_openclosed(x.is_closed()?CLOSED:OPEN);
  return *this;
}

INLINE LBoundary&
LBoundary::operator *= (Boundary::Value x) {
  if (x == 0.0f) {
    value = 0.0f;
    set_openclosed(CLOSED);
  }
  else if (value != 0.0f)
    COMPUTE_MUL((value = value * x));
  return *this;
}

INLINE LBoundary&
LBoundary::operator *= (const Boundary& x) {
  if (x.value == 0.0f) {
    value = 0.0f;
    if (x.is_closed())
      set_openclosed(CLOSED);
  }
  else if (value != 0.0f) {
    COMPUTE_MUL((value = value * x.value));
    if (is_closed())
      set_openclosed(x.is_closed()?CLOSED:OPEN);
  }
  return *this;
}

INLINE LBoundary&
LBoundary::operator /= (Boundary::Value x) {
  if (value != 0.0f)
    COMPUTE_DIV((value = value / x));
  return *this;
}

INLINE LBoundary&
LBoundary::operator /= (const Boundary& x) {
  if (value != 0.0f) {
    COMPUTE_DIV((value = value / x.value));
    if (is_closed())
      set_openclosed(x.is_closed()?CLOSED:OPEN);
  }
  return *this;
}

INLINE UBoundary&
UBoundary::operator += (Boundary::Value x) {
  COMPUTE_ADD((value = value + x));
  return *this;
}

INLINE UBoundary&
UBoundary::operator += (const Boundary& x) {
  COMPUTE_ADD((value = value + x.value));
  if (is_closed())
    set_openclosed(x.is_closed()?CLOSED:OPEN);
  return *this;
}

INLINE UBoundary&
UBoundary::operator -= (Boundary::Value x) {
  COMPUTE_SUB((value = value - x));
  return *this;
}

INLINE UBoundary&
UBoundary::operator -= (const Boundary& x) {
  COMPUTE_SUB((value = value - x.value));
  if (is_closed())
    set_openclosed(x.is_closed()?CLOSED:OPEN);
  return *this;
}

INLINE UBoundary&
UBoundary::operator *= (Boundary::Value x) {
  if (x == 0.0f) {
    value = 0.0f;
    set_openclosed(CLOSED);
  }
  else if (value != 0.0f)
    COMPUTE_MUL((value = value * x));
  return *this;
}

INLINE UBoundary&
UBoundary::operator *= (const Boundary& x) {
  if (x.value == 0.0f) {
    value = 0.0f;
    if (x.is_closed())
      set_openclosed(CLOSED);
  }
  else if (value != 0.0f) {
    COMPUTE_MUL((value = value * x.value));
    if (is_closed())
      set_openclosed(x.is_closed()?CLOSED:OPEN);
  }
  return *this;
}

INLINE UBoundary&
UBoundary::operator /= (Boundary::Value x) {
  if (value != 0.0f)
    COMPUTE_DIV((value = value / x));
  return *this;
}

INLINE UBoundary&
UBoundary::operator /= (const Boundary& x) {
  if (value != 0.0f) {
    COMPUTE_DIV((value = value / x.value));
    if (is_closed())
      set_openclosed(x.is_closed()?CLOSED:OPEN);
  }
  return *this;
}

// Build some matematic unary function with assign
INLINE LBoundary&
LBoundary::exp_assign() {
  COMPUTE_EXP((value = exp(value)));
  return *this;
}

INLINE LBoundary&
LBoundary::log_assign() {
  COMPUTE_LOG((value = log(value)));
  return *this;
}

INLINE LBoundary&
LBoundary::sqrt_assign() {
  COMPUTE_SQRT((value = sqrt(value)));
  return *this;
}

INLINE LBoundary&
LBoundary::sin_assign() {
  COMPUTE_SIN((value = sin(value)));
  return *this;
}

INLINE LBoundary&
LBoundary::cos_assign() {
  COMPUTE_COS((value = cos(value)));
  return *this;
}

INLINE LBoundary&
LBoundary::tan_assign() {
  COMPUTE_TAN((value = tan(value)));
  return *this;
}

INLINE UBoundary&
UBoundary::exp_assign() {
  COMPUTE_EXP((value = exp(value)));
  return *this;
}

INLINE UBoundary&
UBoundary::log_assign() {
  COMPUTE_LOG((value = log(value)));
  return *this;
}

INLINE UBoundary&
UBoundary::sqrt_assign() {
  COMPUTE_SQRT((value = sqrt(value)));
  return *this;
}

INLINE UBoundary&
UBoundary::sin_assign() {
  COMPUTE_SIN((value = sin(value)));
  return *this;
}

INLINE UBoundary&
UBoundary::cos_assign() {
  COMPUTE_COS((value = cos(value)));
  return *this;
}

INLINE UBoundary&
UBoundary::tan_assign() {
  COMPUTE_TAN((value = tan(value)));
  return *this;
}

// Build some matematic binary function with assign
INLINE LBoundary&
LBoundary::pow_assign(Boundary::Value x) {
  if (x == 0.0f) {
    value = 1.0f;
    set_openclosed(CLOSED);
  }
  else if (value != 0.0f && value != 1.0f) {
    if (x == 0.5f)
      COMPUTE_SQRT((value = sqrt(value)));
    else if (x == 2.0f)
      COMPUTE_MUL((value *= value));
    else
      COMPUTE_POW((value = pow(value, x)));
  }
  return *this;
}

INLINE LBoundary&
LBoundary::pow_assign(const Boundary& x) {
  if (value == 0.0f)
    ;
  else if (x.value == 0.0f) {
    value = 1.0f;
    if (x.is_closed())
      set_openclosed(CLOSED);
  }
  else if (value != 1.0f) {
    if (x == 0.5f)
      COMPUTE_SQRT((value = sqrt(value)));
    else if (x == 2.0f)
      COMPUTE_MUL((value *= value));
    else
      COMPUTE_POW((value = pow(value, x.value)));
    if (is_closed())
      set_openclosed(x.is_closed()?CLOSED:OPEN);
  }
  return *this;
}

INLINE UBoundary&
UBoundary::pow_assign(Boundary::Value x) {
  if (x == 0.0f) {
    value = 1.0f;
    set_openclosed(CLOSED);
  }
  else if (value != 0.0f && value != 1.0f) {
    if (x == 0.5f)
      COMPUTE_SQRT((value = sqrt(value)));
    else if (x == 2.0f)
      COMPUTE_MUL((value *= value));
    else
      COMPUTE_POW((value = pow(value, x)));
  }
  return *this;
}

INLINE UBoundary&
UBoundary::pow_assign(const Boundary& x) {
  if (x.value == 0.0f) {
    value = 1.0f;
    if (x.is_closed())
      set_openclosed(CLOSED);
  }
  else if (value != 0.0f && value != 1.0f) {
    if (x == 0.5f)
      COMPUTE_SQRT((value = sqrt(value)));
    else if (x == 2.0f)
      COMPUTE_MUL((value *= value));
    else
      COMPUTE_POW((value = pow(value, x.value)));
    if (is_closed())
      set_openclosed(x.is_closed()?CLOSED:OPEN);
  }
  return *this;
}

#if 0

// Unary functors without assignment built from ones with it.
INLINE LBoundary
operator - (const LBoundary& x) {
  LBoundary z(x);
  z.negate();
  return z;
}

INLINE LBoundary
floor(const LBoundary& x) {
  LBoundary z(x);
  z.floor_assign();
  return z;
}

INLINE LBoundary
ceil(const LBoundary& x) {
  LBoundary z(x);
  z.ceil_assign();
  return z;
}

INLINE LBoundary
integer(const LBoundary& x) {
  LBoundary z(x);
  z.integer_assign();
  return z;
}

INLINE LBoundary
exp(const LBoundary& x) {
  LBoundary z(x);
  z.exp_assign();
  return z;
}

INLINE LBoundary
log(const LBoundary& x) {
  LBoundary z(x);
  z.log_assign();
  return z;
}

INLINE LBoundary
sqrt(const LBoundary& x) {
  LBoundary z(x);
  z.sqrt_assign();
  return z;
}

INLINE LBoundary
sin(const LBoundary& x) {
  LBoundary z(x);
  z.sin_assign();
  return z;
}

INLINE LBoundary
cos(const LBoundary& x) {
  LBoundary z(x);
  z.cos_assign();
  return z;
}

INLINE LBoundary
tan(const LBoundary& x) {
  LBoundary z(x);
  z.tan_assign();
  return z;
}

INLINE UBoundary
operator - (const UBoundary& x) {
  UBoundary z(x);
  z.negate();
  return z;
}

INLINE UBoundary
floor(const UBoundary& x) {
  UBoundary z(x);
  z.floor_assign();
  return z;
}

INLINE UBoundary
ceil(const UBoundary& x) {
  UBoundary z(x);
  z.ceil_assign();
  return z;
}

INLINE UBoundary
integer(const UBoundary& x) {
  UBoundary z(x);
  z.integer_assign();
  return z;
}

INLINE UBoundary
exp(const UBoundary& x) {
  UBoundary z(x);
  z.exp_assign();
  return z;
}

INLINE UBoundary
log(const UBoundary& x) {
  UBoundary z(x);
  z.log_assign();
  return z;
}

INLINE UBoundary
sqrt(const UBoundary& x) {
  UBoundary z(x);
  z.sqrt_assign();
  return z;
}

INLINE UBoundary
sin(const UBoundary& x) {
  UBoundary z(x);
  z.sin_assign();
  return z;
}

INLINE UBoundary
cos(const UBoundary& x) {
  UBoundary z(x);
  z.cos_assign();
  return z;
}

INLINE UBoundary
tan(const UBoundary& x) {
  UBoundary z(x);
  z.tan_assign();
  return z;
}

// Binary functors without assignment built from ones with it.
INLINE LBoundary
operator + (const LBoundary& x, Boundary::Value y) {
  LBoundary z(x);
  z.operator +=(y);
  return z;
}

INLINE LBoundary
operator + (const LBoundary& x, const Boundary& y) {
  LBoundary z(x);
  z.operator +=(y);
  return z;
}

INLINE LBoundary
operator - (const LBoundary& x, Boundary::Value y) {
  LBoundary z(x);
  z.operator -=(y);
  return z;
}

INLINE LBoundary
operator - (const LBoundary& x, const Boundary& y) {
  LBoundary z(x);
  z.operator -=(y);
  return z;
}

INLINE LBoundary
operator * (const LBoundary& x, Boundary::Value y) {
  LBoundary z(x);
  z.operator *=(y);
  return z;
}

INLINE LBoundary
operator * (const LBoundary& x, const Boundary& y) {
  LBoundary z(x);
  z.operator *=(y);
  return z;
}

INLINE LBoundary
operator / (const LBoundary& x, Boundary::Value y) {
  LBoundary z(x);
  z.operator /=(y);
  return z;
}

INLINE LBoundary
operator / (const LBoundary& x, const Boundary& y) {
  LBoundary z(x);
  z.operator /=(y);
  return z;
}

INLINE LBoundary
pow(const LBoundary& x, Boundary::Value y) {
  LBoundary z(x);
  z.pow_assign(y);
  return z;
}

INLINE LBoundary
pow(const LBoundary& x, const Boundary& y) {
  LBoundary z(x);
  z.pow_assign(y);
  return z;
}

INLINE UBoundary
operator + (const UBoundary& x, Boundary::Value y) {
  UBoundary z(x);
  z.operator +=(y);
  return z;
}

INLINE UBoundary
operator + (const UBoundary& x, const Boundary& y) {
  UBoundary z(x);
  z.operator +=(y);
  return z;
}

INLINE UBoundary
operator - (const UBoundary& x, Boundary::Value y) {
  UBoundary z(x);
  z.operator -=(y);
  return z;
}

INLINE UBoundary
operator - (const UBoundary& x, const Boundary& y) {
  UBoundary z(x);
  z.operator -=(y);
  return z;
}

INLINE UBoundary
operator * (const UBoundary& x, Boundary::Value y) {
  UBoundary z(x);
  z.operator *=(y);
  return z;
}

INLINE UBoundary
operator * (const UBoundary& x, const Boundary& y) {
  UBoundary z(x);
  z.operator *=(y);
  return z;
}

INLINE UBoundary
operator / (const UBoundary& x, Boundary::Value y) {
  UBoundary z(x);
  z.operator /=(y);
  return z;
}

INLINE UBoundary
operator / (const UBoundary& x, const Boundary& y) {
  UBoundary z(x);
  z.operator /=(y);
  return z;
}

INLINE UBoundary
pow(const UBoundary& x, Boundary::Value y) {
  UBoundary z(x);
  z.pow_assign(y);
  return z;
}

INLINE UBoundary
pow(const UBoundary& x, const Boundary& y) {
  UBoundary z(x);
  z.pow_assign(y);
  return z;
}

#endif
