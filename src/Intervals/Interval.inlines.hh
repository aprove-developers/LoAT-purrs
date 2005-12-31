/* Interval class implementation: inline functions.
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

#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdlib>

using std::cerr;
using std::endl;
using std::fabs;

// Set interval to an empty one
INLINE void
Interval::set_empty() {
  lower.set_max();
  upper.set_min();
  flag = INTEGER;
}

// Return true if interval is empty
INLINE bool
Interval::is_empty() const {
  return lower > upper;
}

INLINE bool
Interval::is_integer() const {
  return flag == INTEGER;
}

INLINE bool
Interval::is_lower_min() const {
  return lower.is_min();
}

INLINE bool
Interval::is_upper_max() const {
  return upper.is_max();
}

// Set lower boundary to the minimum possible
INLINE void
Interval::set_lower_min() {
  lower.set_min();
}

// Set upper boundary to the maximum possible
INLINE void
Interval::set_upper_max() {
  upper.set_max();
}

// Set interval to R or Z depending on nf
INLINE void
Interval::set_full(Type nf) {
  set_lower_min();
  set_upper_max();
  flag = nf;
}

// Return true if is non empty or normalized empty (+inf, -inf)
INLINE bool
Interval::ok() const {
  return
    lower <= upper
    || (lower.is_max() && upper.is_min())
    || (cerr << lower << " " << upper << endl, false);
}

// Return true if the interval contain only a value
INLINE bool
Interval::is_one_point() const {
  return lower == upper;
}

// Set interval type to INTEGER, if contain only an integer value,
// REAL otherwise
INLINE void
Interval::set_maybe_real() {
  flag = (is_one_point() && lower.is_integer()) ? INTEGER : REAL;
}

// Set interval type to INTEGER truncating boundary
INLINE void
Interval::set_integer() {
  lower.ceil_assign();
  upper.floor_assign();
  flag = INTEGER;
  if (lower > upper)
    set_empty();
}

INLINE
Interval::Interval() {
  set_full(REAL);
}

// Constructs an interval containing only one value.
// NOTE: This is called implicitly to convert a Boundary::Value to an interval.
INLINE
Interval::Interval(Boundary::Value x)
: lower(LBoundary(x, LBoundary::CLOSED)),
  upper(UBoundary(x, UBoundary::CLOSED)),
  flag(lower.is_integer()?INTEGER:REAL)
{ }

// Constructs an interval containing an integer value.
INLINE
Interval::Interval(int x)
: lower(LBoundary((Boundary::Value) x, LBoundary::CLOSED)),
  upper(UBoundary((Boundary::Value) x, UBoundary::CLOSED)),
  flag(INTEGER)
{ }

#ifdef Integer_INTERFACE
// Constructs an interval containing an integer value.
INLINE
Interval::Interval(const Integer& x)
  // KLUDGE!!!
: lower(LBoundary((Boundary::Value) to_double(x), LBoundary::CLOSED)),
  upper(UBoundary((Boundary::Value) to_double(x), UBoundary::CLOSED)),
  flag(INTEGER)
{ }
#endif

#if !BOUNDARY_TYPE_IS_DOUBLE

extern double id_double(double);

INLINE
Interval::Interval(double x)
: lower(LBoundary(x, LBoundary::CLOSED)),
  upper(UBoundary(id_double(x), UBoundary::CLOSED)) {
  set_maybe_real();
}
#endif

// Constructor for interval containing only one value.
// NOTE: This is called implicitly to convert a string to an interval.
INLINE
Interval::Interval(const char* x)
: lower(LBoundary(x, LBoundary::CLOSED)),
  upper(UBoundary(x, UBoundary::CLOSED)) {
  set_maybe_real();
}

// Complete interval constructor
INLINE
Interval::Interval(const LBoundary& lb, const UBoundary& ub, Type flag)
: lower(lb),
  upper(ub) {
  assert(lb <= ub);
  if (flag == INTEGER)
    set_integer();
  else
    set_maybe_real();
}

INLINE bool
Interval::is_strong() const {
  return is_one_point();
}

INLINE Boundary::Value
Interval::value() const {
  assert(is_strong());
  return lower.boundary();
}

INLINE void
Interval::destructure(Boundary::Value& lb_value, bool& lb_closed, 
		      Boundary::Value& ub_value, bool& ub_closed,
		      bool& integer) const {
  lb_value = lower.boundary();
  lb_closed = lower.is_closed();
  ub_value = upper.boundary();
  ub_closed = upper.is_closed();
  integer = (flag == INTEGER);
}

// Return true if interval is the full one (of type f)
INLINE bool
Interval::is_full(Type f) const {
  return lower.is_min() && upper.is_max() && flag == f;
}

// Return true if x is strictly contained in y
INLINE bool
subset(const Interval& x, const Interval& y) {
  assert(x.ok() && y.ok());
  if (x.flag == Interval::REAL) {
    if (y.flag == Interval::INTEGER)
      return false;
  }
  else if (y.flag == Interval::REAL)
    return (x.upper <= y.upper && x.lower >= y.lower);
  return ((x.upper <  y.upper && x.lower >= y.lower) ||
	  (x.upper <= y.upper && x.lower >  y.lower));
}

// Return true if y is strictly contained in x
INLINE bool
superset(const Interval& x, const Interval& y) {
  return subset(y, x);
}

// Return true if x is contained in or equal to y
INLINE bool
subset_eq(const Interval& x, const Interval& y) {
  assert(x.ok() && y.ok());
  if (x.flag == Interval::REAL) {
    if (y.flag == Interval::INTEGER)
      return false;
  }
  return (x.upper <= y.upper && x.lower >= y.lower);
}

// Return true if y is contained in or equal to x
INLINE bool
superset_eq(const Interval& x, const Interval& y) {
  return subset_eq(y, x);
}

// Return true if x and y contains the same values
INLINE bool
operator == (const Interval& x, const Interval& y) {
  assert(x.ok() && y.ok());
  return (x.lower == y.lower && x.upper == y.upper && x.flag == y.flag);
}

// Return true if x and y not contains the same values
INLINE bool
operator != (const Interval& x, const Interval& y) {
  // RB: was assert(!x.ok() && !y.ok());
  assert(x.ok() && y.ok());
  return !(x == y);
}

// Predicates meaning: forall X in x forall Y in y : p(X, Y)
// Less than
INLINE bool
lt(const Interval& x, const Interval& y) {
  assert(x.ok() && y.ok());
  return x.upper < y.lower;
}

// Greater than
INLINE bool
gt(const Interval& x, const Interval& y) {
  return lt(y, x);
}

// Less or equal
INLINE bool
le(const Interval& x, const Interval& y) {
  assert(x.ok() && y.ok());
  return x.upper <= y.lower;
}

// Greater or equal
INLINE bool
ge(const Interval& x, const Interval& y) {
  return le(y, x);
}

// Equal
INLINE bool
eq(const Interval& x, const Interval& y) {
  return le(x, y) && ge(x, y);
}

// Not equal
INLINE bool
ne(const Interval& x, const Interval& y) {
  return lt(x, y) || gt(x, y);
}

// Predicates meaning: forall X in x : p(X, y)
// Less than
INLINE bool
lt(const Interval& x, Boundary::Value y) {
  assert(x.ok());
  return x.upper < y;
}

// Greater than
INLINE bool
gt(const Interval& x, Boundary::Value y) {
  assert(x.ok());
  return x.lower > y;
}

// Less or equal
INLINE bool
le(const Interval& x, Boundary::Value y) {
  assert(x.ok());
  return x.upper <= y;
}

// Greater or equal
INLINE bool
ge(const Interval& x, Boundary::Value y) {
  assert(x.ok());
  return x.lower >= y;
}

// Equal
INLINE bool
eq(const Interval& x, Boundary::Value y) {
  assert(x.ok());
  return le(x, y) && ge(x, y);
}

// Not equal
INLINE bool
ne(const Interval& x, Boundary::Value y) {
  assert(x.ok());
  return lt(x, y) || gt(x, y);
}

#if 0
// Return the stricter relation between x and y
INLINE Rel
findrel(const Interval& x, const Interval& y) {
  assert(!x.is_empty() && !y.is_empty());
  int c1 = cmp(x.upper, y.lower);
  int c2 = cmp(x.lower, y.upper);
  if (c1 < 0)		// xu < yl
    return LT;
  else {		// xu >= yl
    if (c2 > 0)		// xu >= yl && xl > yu
      return GT;
    else {		// xu >= yl && xl <= yu
      if (c1 == 0) {	// xu = yl && xl <= yu
	if (c2 == 0)	// xu = yl && xl = yu
	  return EQ;
	else		// xu = yl && xl < yu
	  return LE;
      }
      else {		// xu > yl && xl <= yu
	if (c2 == 0)	// xu >= yl && xl = yu
	  return GE;
	else		// xu > yl && xl < yu
	  return UN;
      }
    }
  }
}
#endif

// Lexicographical comparison
INLINE int
lcompare(const Interval& x, const Interval& y) {
  assert(x.ok() && y.ok());
  return (x.lower < y.lower ? -1 :
	  y.lower < x.lower ? 1 :
	  x.upper < y.upper ? -1 :
	  y.upper < x.upper ? 1 :
	  0);
}

// 3 state comparison predicate
// Less than
INLINE Bool3
lt3(const Interval& x, const Interval& y) {
  return lt(x, y) ? ALWAYS : ge(x, y) ? NEVER : MAYBE;
}
// Greater than
INLINE Bool3
gt3(const Interval& x, const Interval& y) {
  return gt(x, y) ? ALWAYS : le(x, y) ? NEVER : MAYBE;
}
// Less or equal
INLINE Bool3
le3(const Interval& x, const Interval& y) {
  return le(x, y) ? ALWAYS : gt(x, y) ? NEVER : MAYBE;
}
// Greater or equal
INLINE Bool3
ge3(const Interval& x, const Interval& y) {
  return ge(x, y) ? ALWAYS : lt(x, y) ? NEVER : MAYBE;
}
// Equal
INLINE Bool3
eq3(const Interval& x, const Interval& y) {
  return eq(x, y) ? ALWAYS : ne(x, y) ? NEVER : MAYBE;
}
// Not equal
INLINE Bool3
ne3(const Interval& x, const Interval& y) {
  return ne(x, y) ? ALWAYS : eq(x, y) ? NEVER : MAYBE;
}

// Less than
INLINE Bool3
lt3(const Interval& x, Boundary::Value y) {
  return lt(x, y) ? ALWAYS : ge(x, y) ? NEVER : MAYBE;
}
// Greater than
INLINE Bool3
gt3(const Interval& x, Boundary::Value y) {
  return gt(x, y) ? ALWAYS : le(x, y) ? NEVER : MAYBE;
}
// Less or equal
INLINE Bool3
le3(const Interval& x, Boundary::Value y) {
  return le(x, y) ? ALWAYS : gt(x, y) ? NEVER : MAYBE;
}
// Greater or equal
INLINE Bool3
ge3(const Interval& x, Boundary::Value y) {
  return ge(x, y) ? ALWAYS : lt(x, y) ? NEVER : MAYBE;
}
// Equal
INLINE Bool3
eq3(const Interval& x, Boundary::Value y) {
  return eq(x, y) ? ALWAYS : ne(x, y) ? NEVER : MAYBE;
}
// Not equal
INLINE Bool3
ne3(const Interval& x, Boundary::Value y) {
  return ne(x, y) ? ALWAYS : eq(x, y) ? NEVER : MAYBE;
}

#if 0
INLINE Bool3
cmp3(const Interval& x, Rel r, const Interval& y) {
  assert(defined(r));
  switch(r) {
  case LT:
    return lt3(x, y);
  case EQ:
    return eq3(x, y);
  case GT:
    return gt3(x, y);
  case LE:
    return le3(x, y);
  case NE:
    return ne3(x, y);
  case GE:
    return ge3(x, y);
  default:
    abort();
    return NEVER;
  }
}

INLINE Bool3
cmp3(const Interval& x, Rel r, Boundary::Value y) {
  assert(defined(r));
  switch(r) {
  case LT:
    return lt3(x, y);
  case EQ:
    return eq3(x, y);
  case GT:
    return gt3(x, y);
  case LE:
    return le3(x, y);
  case NE:
    return ne3(x, y);
  case GE:
    return ge3(x, y);
  default:
    abort();
    return NEVER;
  }
}
#endif

// Intersect this interval with y
INLINE Interval&
Interval::operator &= (const Interval& y) {
  if (y.lower > lower)
    lower = y.lower;
  if (y.upper < upper)
    upper = y.upper;
  if (flag == REAL && y.flag == REAL) {
    if (lower > upper)
      set_empty();
    else
      set_maybe_real();
  }
  else
    set_integer();
  return *this;
}

// Return true if the intersection between x and y is not empty
INLINE bool
operator && (const Interval& x, const Interval& y) {
  Interval i(x);
  i &= y;
  return !i.is_empty();
}

// Return true if y is contained in x
INLINE bool
operator && (const Interval& x, Boundary::Value y) {
  return (x.lower <= y &&
	  x.upper >= y &&
	  (x.flag == Interval::REAL || Boundary::is_integer(y)));
}

// Join this interval with y
// NOTE: it fill the hole
INLINE Interval&
Interval::operator |= (const Interval& y) {
  assert(ok() && y.ok());
  if (y.lower < lower)
    lower = y.lower;
  if (y.upper > upper)
    upper = y.upper;
  if (y.flag == REAL)
    flag = REAL;
  return *this;
}

INLINE Interval&
Interval::operator += (Boundary::Value y) {
  round_save();
  lower += y;
  upper += y;
  round_restore();
  if (flag != Interval::INTEGER || !Boundary::is_integer(y))
    set_maybe_real();
  else
    set_integer();
  return *this;
}

INLINE Interval&
Interval::operator -= (Boundary::Value y) {
  round_save();
  lower -= y;
  upper -= y;
  round_restore();
  if (flag != Interval::INTEGER || !Boundary::is_integer(y))
    set_maybe_real();
  else
    set_integer();
  return *this;
}

// a op b = { A op B | A in a && B in b }
// Addition
INLINE Interval&
Interval::operator += (const Interval& y) {
  assert(ok() && y.ok());
  if (y.is_one_point())
    return *this += y.value();
#if IEEE_NUMBERS
  if ((eq(*this, BOUNDARY_MIN) && eq(y, BOUNDARY_MAX)) ||
      (eq(*this, BOUNDARY_MAX) && eq(y, BOUNDARY_MIN)))
    set_empty();
  else
#endif
    {
      round_save();
      lower += y.lower;
      upper += y.upper;
      round_restore();
      if (flag != INTEGER || y.flag != INTEGER)
	set_maybe_real();
      else
	set_integer();
    }
  return *this;
}

// Subtraction
INLINE Interval&
Interval::operator -= (const Interval& yb) {
  assert(ok() && yb.ok());
  if (yb.is_one_point())
    return *this -= yb.value();
#if IEEE_NUMBERS
  if ((eq(*this, BOUNDARY_MIN) && eq(yb, BOUNDARY_MIN)) ||
      (eq(*this, BOUNDARY_MAX) && eq(yb, BOUNDARY_MAX)))
    set_empty();
  else
#endif
    {
      round_save();
      // To handle correctly a -= a
      const Interval& y = (this == &yb) ? Interval(yb) : yb;
      lower -= y.upper;
      upper -= y.lower;
      round_restore();
      if (flag != INTEGER || y.flag != INTEGER)
	set_maybe_real();
      else
	set_integer();
    }
  return *this;
}

INLINE Interval&
Interval::operator *= (Boundary::Value y) {
  if (y<0.0f) {
    LBoundary lb(lower);
    lower.set(upper);
    upper.set(lb);
  }
  round_save();
  lower *= y;
  upper *= y;
  round_restore();
  if (flag != Interval::INTEGER || !Boundary::is_integer(y))
    set_maybe_real();
  else
    set_integer();
  return *this;
}

INLINE Interval&
Interval::operator /= (Boundary::Value y) {
  if (y == 0.0f)
    set_empty();
  else {    
    if (y<0.0f) {
      LBoundary lb(lower);
      lower.set(upper);
      upper.set(lb);
    }
    round_save();
    lower /= y;
    upper /= y;
    round_restore();
    int exp;
    if (flag == Interval::INTEGER && 
	y >= -1.0f && y <= 1.0f && fabs(frexp(y, &exp)) == 0.5f)
      set_integer();
    else
      set_maybe_real();
  }
  return *this;
}

// f(a, b) = { f(A, B) | A in a && B in b }
// Minimum
INLINE Interval&
Interval::min_assign (const Interval& y) {
  assert(ok() && y.ok());
  if (y.lower < lower) {
    if (y.upper < upper)
      *this = y;
    else {
      lower = y.lower;
      if (y.flag == REAL)
	flag = REAL;
    }
  }
  else if (y.upper < upper) {
    upper = y.upper;
    if (y.flag == REAL)
      flag = REAL;
  }
  return *this;
}

// Maximum
INLINE Interval&
Interval::max_assign (const Interval& y) {
  assert(ok() && y.ok());
  if (y.lower > lower) {
    if (y.upper > upper)
      *this = y;
    else {
      lower = y.lower;
      if (y.flag == REAL)
	flag = REAL;
    }
  }
  else if (y.upper > upper) {
    upper = y.upper;
    if (y.flag == REAL)
      flag = REAL;
  }
  return *this;
}

// op a = { op A | A in a }
// Negation (unary minus)
INLINE Interval&
Interval::negate () {
  LBoundary lb(lower);
  lower.set(upper);
#if 0    // I think this is not useful: IEEE 754 has a sign bit to toggle
  round_save();
#endif
  lower.negate();
  upper.set(lb);
  upper.negate();
#if 0
  round_restore();
#endif
  return *this;
}

// f(a) = { f(A) | A in a }
// Absolute value
INLINE Interval&
Interval::abs_assign() {
  assert(ok());
  if (lt(*this,0.0))
    negate();
  else if (!gt(*this, 0.0)) {
    UBoundary ub;
    ub.set(lower);
    ub.negate();
    if (ub > upper)
      upper = ub;
    lower.set(0.0f);
  }
  return *this;
}

// Square (x*x)
INLINE Interval&
Interval::sqr_assign() {
  assert(ok());
  round_save();
  if (lt(*this, 0.0f)) {
    LBoundary lb(lower);
    lower.set(upper);
    lower *= lower;
    upper.set(lb);
    upper *= lower;
  }
  else if (gt(*this, 0.0)) {
    lower *= lower;
    upper *= upper;
  }
  else {
    UBoundary ub;
    ub.set(lower);
    ub.negate();
    if (ub > upper)
      upper = ub;
    upper *= upper;
    lower.set(0.0f);
  }
  round_restore();
  if (flag != INTEGER)
    set_maybe_real();
  else
    set_integer();
  return *this;
}

// Square root
INLINE Interval&
Interval::sqrt_assign() {
  assert(ok());
  if (lt(*this, 0.0))
    set_empty();
  else {
    round_save();
    if (gt(*this, 0.0))
      lower.sqrt_assign();
    else
      lower.set(0.0f);
    upper.sqrt_assign();
    round_restore();
    set_maybe_real();
  }
  return *this;
}

// e^x
INLINE Interval&
Interval::exp_assign() {
  assert(ok());
  round_save();
  lower.exp_assign();
  upper.exp_assign();
  round_restore();
  set_maybe_real();
  return *this;
}

// Natural logarithm
INLINE Interval&
Interval::log_assign() {
  assert(ok());
#if IEEE_NUMBERS
  if (eq(*this, 0.0)) {
    *this = Interval(BOUNDARY_MIN);
  }
  else if (lt(*this, 0.0))
    set_empty();
#else
  if (le(*this, 0.0))
    set_empty();
#endif
  else {
    round_save();
    if (!gt(*this, 0.0))
      set_lower_min();
    else
      lower.log_assign();
    upper.log_assign();
    round_restore();
    set_maybe_real();
  }
  return *this;
}

// Round to largest integer not greater than
INLINE Interval&
Interval::floor_assign() {
  assert(ok());
  lower.floor_assign();
  upper.floor_assign();
  flag = INTEGER;	
  return *this;
}

// Round to smallest integer not less than
INLINE Interval&
Interval::ceil_assign() {
  assert(ok());
  lower.ceil_assign();
  upper.ceil_assign();
  flag = INTEGER;	
  return *this;
}

// Round to closest integer toward 0
INLINE Interval&
Interval::integer_assign() {
  assert(ok());
  lower.integer_assign();
  upper.integer_assign();
  flag = INTEGER;	
  return *this;
}

// x div y = int(x/y)
INLINE Interval&
Interval::div_assign(const Interval& y) {
  assert(ok());
  *this /= y;
  integer_assign();
  return *this;
}

// x mod y = x-y*div(x,y)
// TODO: better implementation
INLINE Interval&
Interval::mod_assign(const Interval& y) {
  assert(ok());
  Interval tmp(*this);
  tmp.div_assign(y);
  tmp *= y;
  *this -= tmp;
  return *this;
}


#if 1

// Unary functors without assignment built from ones with it.
INLINE Interval
operator - (const Interval& x) {
  Interval z(x);
  z.negate();
  return z;
}

INLINE Interval
abs(const Interval& x) {
  Interval z(x);
  z.abs_assign();
  return z;
}

INLINE Interval
sqr(const Interval& x) {
  Interval z(x);
  z.sqr_assign();
  return z;
}

INLINE Interval
sqrt(const Interval& x) {
  Interval z(x);
  z.sqrt_assign();
  return z;
}

INLINE Interval
exp(const Interval& x) {
  Interval z(x);
  z.exp_assign();
  return z;
}

INLINE Interval
log(const Interval& x) {
  Interval z(x);
  z.log_assign();
  return z;
}

INLINE Interval
floor(const Interval& x) {
  Interval z(x);
  z.floor_assign();
  return z;
}

INLINE Interval
ceil(const Interval& x) {
  Interval z(x);
  z.ceil_assign();
  return z;
}

INLINE Interval
integer(const Interval& x) {
  Interval z(x);
  z.integer_assign();
  return z;
}

INLINE Interval
sin(const Interval& x) {
  Interval z(x);
  z.sin_assign();
  return z;
}

INLINE Interval
cos(const Interval& x) {
  Interval z(x);
  z.cos_assign();
  return z;
}

INLINE Interval
tan(const Interval& x) {
  Interval z(x);
  z.tan_assign();
  return z;
}

// Binary functors without assignment built from ones with it.
INLINE Interval
widen(const Interval& x, const Interval& y) {
  Interval z(x);
  z.widen_assign(y);
  return z;
}

INLINE Interval
narrow(const Interval& x, const Interval& y) {
  Interval z(x);
  z.narrow_assign(y);
  return z;
}

INLINE Interval
operator & (const Interval& x, const Interval& y) {
  Interval z(x);
  z.operator &=(y);
  return z;
}

INLINE Interval
operator | (const Interval& x, const Interval& y) {
  Interval z(x);
  z.operator |=(y);
  return z;
}

INLINE Interval
operator + (const Interval& x, const Interval& y) {
  Interval z(x);
  z.operator +=(y);
  return z;
}

INLINE Interval
operator - (const Interval& x, const Interval& y) {
  Interval z(x);
  z.operator -=(y);
  return z;
}

INLINE Interval
operator * (const Interval& x, const Interval& y) {
  Interval z(x);
  z.operator *=(y);
  return z;
}

INLINE Interval
operator / (const Interval& x, const Interval& y) {
  Interval z(x);
  z.operator /=(y);
  return z;
}

INLINE Interval
min(const Interval& x, const Interval& y) {
  Interval z(x);
  z.min_assign(y);
  return z;
}

INLINE Interval
max(const Interval& x, const Interval& y) {
  Interval z(x);
  z.max_assign(y);
  return z;
}

INLINE Interval
pow(const Interval& x, const Interval& y) {
  Interval z(x);
  z.pow_assign(y);
  return z;
}

INLINE Interval
div(const Interval& x, const Interval& y) {
  Interval z(x);
  z.div_assign(y);
  return z;
}

INLINE Interval
mod(const Interval& x, const Interval& y) {
  Interval z(x);
  z.mod_assign(y);
  return z;
}

#endif
