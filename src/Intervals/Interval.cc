/* Interval class implementation (non-inline functions).
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

#include <config.h>

#include "Interval.defs.hh"

#if OUTLINE
#include "Interval.inlines.hh"
#endif

using std::ws;

// exist X in x | X = offset+k*period with k integer
// WARNING: integer interval flag is ignored
bool
Interval::contains_multiple(Boundary::Value offset,
			    Boundary::Value period)
{
  assert(ok());
  LBoundary lb(lower);
  lb -= offset;
  lb /= period;
  if (lb.is_closed() && lb.is_integer())
    return true;
  lb.floor_assign();
  UBoundary ub(upper);
  ub -= offset;
  ub /= period;
  ub.floor_assign();
  return lb < ub;
}

#if 0
// Stated that between x and y exists a relation rel, update x and y
// deleting all values that cannot respect this relation
void
refine(Interval& x, Rel rel, Interval& y)
{
  assert(!x.is_empty() && !y.is_empty());
  assert(defined(rel));
  switch (rel) {
  case LT:
    if (x.upper >= y.upper) {
      x.upper.set_lt(y.upper.boundary(), x.flag == Interval::INTEGER);
      if (x.lower > x.upper) {
	x.set_empty();
	y.set_empty();
	return;
      }
      if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    if (y.lower <= x.lower) {
      y.lower.set_gt(x.lower.boundary(), y.flag == Interval::INTEGER);
      assert(y.upper >= y.lower);
      if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case GT:
    if (x.lower <= y.lower) {
      x.lower.set_gt(y.lower.boundary(), x.flag == Interval::INTEGER);
      if (x.lower > x.upper) {
	x.set_empty();
	y.set_empty();
	return;
      }
      if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    if (y.upper >= x.upper) {
      y.upper.set_lt(x.upper.boundary(), y.flag == Interval::INTEGER);
      assert(y.upper >= y.lower);
      if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case LE:
    if (x.upper > y.upper) {
      x.upper.set(y.upper, x.flag == Interval::INTEGER);
      if (x.lower > x.upper) {
	x.set_empty();
	y.set_empty();
	return;
      }
      if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    if (y.lower < x.lower) {
      y.lower.set(x.lower, y.flag == Interval::INTEGER);
      assert(y.upper >= y.lower);
      if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case GE:
    if (x.lower < y.lower) {
      x.lower.set(y.lower, x.flag == Interval::INTEGER);
      if (x.lower > x.upper) {
	x.set_empty();
	y.set_empty();
	return;
      }
      if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    if (y.upper > x.upper) {
      y.upper.set(x.upper, y.flag == Interval::INTEGER);
      assert(y.upper >= y.lower);
      if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case EQ:
    if (x.lower < y.lower)
      x.lower = y.lower;
    else
      y.lower = x.lower;
    if (x.upper > y.upper)
      x.upper = y.upper;
    else
      y.upper = x.upper;
    if (x.flag == Interval::REAL && y.flag == Interval::REAL) {
      if (x.lower > x.upper) {
	x.set_empty();
	y.set_empty();
      }
      else {
	x.set_maybe_real();
	y.set_maybe_real();
      }
    }
    else {
      x.set_integer();
      y.set_integer();
    }
    break;

  case NE:
    if (x.is_one_point()) {
      if (y.lower == x.lower)
	y.lower.set_gt(x.lower.boundary(), y.flag == Interval::INTEGER);
      if (y.upper == x.upper)
	y.upper.set_lt(x.upper.boundary(), y.flag == Interval::INTEGER);
      if (y.lower > y.upper) {
	y.set_empty();
	y.set_empty();
	return;
      }
    }
    if (y.is_one_point()) {
      if (x.lower == y.lower)
	x.lower.set_gt(y.lower.boundary(), x.flag == Interval::INTEGER);
      if (x.upper == y.upper)
	x.upper.set_lt(y.upper.boundary(), x.flag == Interval::INTEGER);
      assert(y.upper >= y.lower);
    }
    break;

  default:
    ;
  }
}

void
refine1(Interval& x, Rel rel, const Interval& y)
{
  assert(!x.is_empty() && !y.is_empty());
  assert(defined(rel));
  switch (rel) {
  case LT:
    if (x.upper >= y.upper) {
      x.upper.set_lt(y.upper.boundary(), x.flag == Interval::INTEGER);
      if (x.lower > x.upper)
	x.set_empty();
      else if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    break;

  case GT:
    if (x.lower <= y.lower) {
      x.lower.set_gt(y.lower.boundary(), x.flag == Interval::INTEGER);
      if (x.lower > x.upper)
	x.set_empty();
      else if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    break;

  case LE:
    if (x.upper > y.upper) {
      x.upper.set(y.upper, x.flag == Interval::INTEGER);
      if (x.lower > x.upper)
	x.set_empty();
      else if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    break;

  case GE:
    if (x.lower < y.lower) {
      x.lower.set(y.lower, x.flag == Interval::INTEGER);
      if (x.lower > x.upper)
	x.set_empty();
      else if (x.flag != Interval::INTEGER)
	x.set_maybe_real();
    }
    break;

  case EQ:
    if (x.lower < y.lower)
      x.lower = y.lower;
    if (x.upper > y.upper)
      x.upper = y.upper;
    if (x.flag == Interval::REAL && y.flag == Interval::REAL) {
      if (x.lower > x.upper)
	x.set_empty();
      else
	x.set_maybe_real();
    }
    else
      x.set_integer();
    break;

  case NE:
    if (y.is_one_point()) {
      if (x.lower == y.lower)
	x.lower.set_gt(y.lower.boundary(), x.flag == Interval::INTEGER);
      if (x.upper == y.upper)
	x.upper.set_lt(y.upper.boundary(), x.flag == Interval::INTEGER);
      if (x.lower > x.upper)
	x.set_empty();
    }
    break;

  default:
    ;
  }
}

void
refine2(const Interval& x, Rel rel, Interval& y)
{
  assert(!x.is_empty() && !y.is_empty());
  assert(defined(rel));
  switch (rel) {
  case LT:
    if (y.lower <= x.lower) {
      y.lower.set_gt(x.lower.boundary(), y.flag == Interval::INTEGER);
      if (y.lower > y.upper)
	y.set_empty();
      else if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case GT:
    if (y.upper >= x.upper) {
      y.upper.set_lt(x.upper.boundary(), y.flag == Interval::INTEGER);
      if (y.lower > y.upper)
	y.set_empty();
      else if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case LE:
    if (y.lower < x.lower) {
      y.lower.set(x.lower, y.flag == Interval::INTEGER);
      if (y.lower > y.upper)
	y.set_empty();
      else if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case GE:
    if (y.upper > x.upper) {
      y.upper.set(x.upper, y.flag == Interval::INTEGER);
      if (y.lower > y.upper)
	y.set_empty();
      else if (y.flag != Interval::INTEGER)
	y.set_maybe_real();
    }
    break;

  case EQ:
    if (x.lower > y.lower)
      y.lower = x.lower;
    if (x.upper < y.upper)
      y.upper = x.upper;
    if (x.flag == Interval::REAL && y.flag == Interval::REAL) {
      if (y.lower > y.upper)
	y.set_empty();
      else
	y.set_maybe_real();
    }
    else
      y.set_integer();
    break;

  case NE:
    if (x.is_one_point()) {
      if (y.lower == x.lower)
	y.lower.set_gt(x.lower.boundary(), y.flag == Interval::INTEGER);
      if (y.upper == x.upper)
	y.upper.set_lt(x.upper.boundary(), y.flag == Interval::INTEGER);
      if (y.lower > y.upper)
	y.set_empty();
    }
    break;

  default:
    ;
  }
}
#endif

// Widen this interval (with respect to changes from old)
// If new boundaries are wider than old ones, widen further to next (or same)
// magic value: -inf, -1(o|c), 0(o|c), 1(o|c), +inf

static Boundary::Value widen_stop_points[] = { -1.0f, 0.0f, 1.0f };

Interval&
Interval::widen_assign(const Interval& old)
{
  assert(!is_empty() && !old.is_empty());
  if (lower < old.lower) {
    int n = sizeof(widen_stop_points)/sizeof(widen_stop_points[0]);
    int k;
    for (k = n-1; k >= 0; --k) {
      if (lower >= widen_stop_points[k]) {
 	if (lower > widen_stop_points[k])
 	  lower.set_gt(widen_stop_points[k], flag == INTEGER);
	break;
      }
    }
    if (k < 0)
      set_lower_min();
  }
  if (upper > old.upper) {
    int n = sizeof(widen_stop_points)/sizeof(widen_stop_points[0]);
    int k;
    for (k = 0; k < n; ++k) {
      if (upper <= widen_stop_points[k]) {
 	if (upper < widen_stop_points[k])
 	  upper.set_lt(widen_stop_points[k], flag == INTEGER);
	break;
      }
    }
    if (k == n)
      set_upper_max();
  }
  return *this;
}
// Narrow this interval (with respect to changes from old)
// Use new stricter boundaries only if changes is relevant,
// i.e. it cross magic values: -inf, -1, 0, 1, +inf
Interval&
Interval::narrow_assign(const Interval& old)
{
  assert(!is_empty() && !old.is_empty());
  if (lower > old.lower) {
    if (!(old.lower.is_min() ||
	  (old.lower <= -1.0 && lower >= -1.0) ||
	  (old.lower <= 0.0 && lower >= 0.0) ||
	  (old.lower <= 1.0 && lower >= 1.0)))
      lower = old.lower;
  }
  else
    lower = old.lower;
  if (upper < old.upper) {
    if (!(old.upper.is_max() ||
	  (old.upper >= -1.0 && upper <= -1.0) ||
	  (old.upper >= 0.0 && upper <= 0.0) ||
	  (old.upper >= 1.0 && upper <= 1.0)))
      upper = old.upper;
  }
  else
    upper = old.upper;
  if (lower > upper)
    set_empty();
  return *this;
}

// a op b = { A op B | A in a && B in b }
// Multiplication
/**
+---+-----------+-----------------+-----------+
| * |    -y-    |       -y+       |    +y+    |
+---+-----------+-----------------+-----------+
|-x-|xu*yu,xl*yl|   xl*yu,xl*yl   |xl*yu,xu*yl|
+---+-----------+-----------------+-----------+
|-x+|xu*yl,xl*yl|min(xl*yu,xu*yl),|xl*yu,xu*yu|
|   |           |max(xl*yl,xu*yu) |           |
+---+-----------+-----------------+-----------+
|+x+|xu*yl,xl*yu|   xu*yl,xu*yu   |xl*yl,xu*yu|
+---+-----------+-----------------+-----------+
**/
Interval&
Interval::operator *= (const Interval& yb)
{
  assert(ok() && yb.ok());
#if IEEE_NUMBERS
  if ((eq(*this, 0.0f) && (eq(yb, BOUNDARY_MIN) || eq(yb, BOUNDARY_MAX))) ||
      (eq(yb, 0.0f) && (eq(*this, BOUNDARY_MIN) || eq(*this, BOUNDARY_MAX))))
    set_empty();
  else
#endif
  if (yb.is_one_point())
    return *this *= yb.value();
  {
    // To handle correctly a *= a
    const Interval& y = (this == &yb) ? Interval(yb) : yb;
    round_save();
    if (lt(y,0.0)) {
      if (lt(*this,0.0)) {
	// -x- -y-
	LBoundary lb(lower);
	lower.set(upper);
	upper.set(lb);
	lower *= y.upper;
	upper *= y.lower;
      }
      else if (gt(*this,0.0)) {
	// +x+ -y-
	LBoundary lb(lower);
	lower.set(upper);
	upper.set(lb);
	lower *= y.lower;
	upper *= y.upper;
      }
      else {
	// -x+ -y-
	LBoundary lb(lower);
	lower.set(upper);
	upper.set(lb);
	lower *= y.lower;
	upper *= y.lower;
      }
    }
    else if (gt(y,0.0)) {
      if (lt(*this, 0.0)) {
	// -x- +y+
	lower *= y.upper;
	upper *= y.lower;
      }
      else if (gt(*this, 0.0)) {
	// +x+ +y+
	lower *= y.lower;
	upper *= y.upper;
      }
      else {
	// -x+ +y+
	lower *= y.upper;
	upper *= y.upper;
      }
    }
    else {
      if (lt(*this, 0.0)) {
	// -x- -y+
	upper.set(lower);
	lower *= y.upper;
	upper *= y.lower;
      }
      else if (gt(*this, 0.0)) {
	// +x+ -y+
	lower.set(upper);
	lower *= y.lower;
	upper *= y.upper;
      }
      else {
	// -x+ -y+
	UBoundary ub;
	ub.set(lower);
	LBoundary lb(y.lower);
	lower *= y.upper;
	lb *= upper;
	upper *= y.upper;
	ub *= y.lower;
	if (lb < lower)
	  lower = lb;
	if (ub > upper)
	  upper = ub;
      }
    }
    round_restore();
    if (flag != Interval::INTEGER || y.flag != Interval::INTEGER)
      set_maybe_real();
    else
      set_integer();
  }
  return *this;
}

// Division
// WARNING: no check for cases like that:
// [6] / [1..3] gives [2,6] (not [2..6])
/**

+---+-----------+-----------+
| / |    -y-    |    +y+    |
+---+-----------+-----------+
|-x-|xu/yl,xl/yu|xl/yl,xu/yu|
+---+-----------+-----------+
|-x+|xu/yu,xl/yu|xl/yl,xu/yl|
+---+-----------+-----------+
|+x+|xu/yu,xl/yl|xl/yu,xu/yl|
+---+-----------+-----------+

!ieee && integer(y) && y contains 0
+---+-----------+-----------------------+-----------+
| / |    -y-    |  int     -y+   !ieee  |    +y+    |
+---+-----------+-----------------------+-----------+
|-x-|xu/yl,-xl  |xl         ,xu/yu      |xl   ,xu/yu|
+---+-----------+-----------------------+-----------+
|-x+|-xu  ,-xl  |min(xl,-xu),max(-xl,xu)|xl   ,xu   |
+---+-----------+-----------------------+-----------+
|+x+|-xu  ,xl/yl|xl/yu      ,xu         |xl/yu,xu   |
+---+-----------+-----------------------+-----------+
**/
Interval&
Interval::operator /= (const Interval& yb)
{
  assert(ok() && yb.ok());
  if (yb.is_one_point())
    return *this /= yb.value();
  // To handle correctly a /= a
  const Interval& y = (this == &yb) ? Interval(yb) : yb;
  // WARNING: no handling of signed zeroes
  if (y && 0.0) {
#if IEEE_NUMBERS
    set_full();
#else
    if (yb.flag == INTEGER) {
      round_save();
      if (le(y, 0.0)) {
	if (lt(*this, 0.0)) {
	  // -x- -y-
	  LBoundary lb(lower);
	  lower.set(upper);
	  upper.set(lb);
	  lower /= y.lower;
	  upper.negate();
	}
	else if (gt(*this, 0.0)) {
	  // +x+ -y-
	  LBoundary lb(lower);
	  lower.set(upper);
	  upper.set(lb);
	  lower.negate();
	  upper /= y.lower;
	}
	else {
	  // -x+ -y-
	  LBoundary lb(lower);
	  lower.set(upper);
	  upper.set(lb);
	  lower.negate();
	  upper.negate();
	}
      }
      else if (ge(y, 0.0)) {
	if (lt(*this, 0.0)) {
	  // -x- +y+
	  upper /= y.upper;
	}
	else if (gt(*this, 0.0)) {
	  // +x+ +y+
	  lower /= y.upper;
	}
	else {
	  // -x+ +y+
	}
      }
      else {
	if (lt(*this, 0.0)) {
	  // -x- -y+
	  upper /= y.upper;
	}
	else if (gt(*this, 0.0)) {
	  // +x+ -y+
	  lower /= y.upper;
	}
	else {
	  // -x+ -y+
	  UBoundary ub;
	  ub.set(lower);
	  ub.negate();
	  LBoundary lb;
	  lb.set(upper);
	  lb.negate();
	  if (ub > upper)
	    upper = ub;
	  if (lb < lower)
	    lower = lb;
	}
      }
      round_restore();
      set_maybe_real();
    }
    else
      set_full();
#endif
  }
  else {
    round_save();
    if (lt(y, 0.0)) {
      if (lt(*this, 0.0)) {
	// -x- -y-
	LBoundary lb(lower);
	lower.set(upper);
	upper.set(lb);
	lower /= y.lower;
	upper /= y.upper;
      }
      else if (gt(*this, 0.0)) {
	// +x+ -y-
	LBoundary lb(lower);
	lower.set(upper);
	upper.set(lb);
	lower /= y.upper;
	upper /= y.lower;
      }
      else {
	// -x+ -y-
	LBoundary lb(lower);
	lower.set(upper);
	upper.set(lb);
	lower /= y.upper;
	upper /= y.upper;
      }
    }
    else {
      if (lt(*this, 0.0)) {
	// -x- +y+
	lower /= y.lower;
	upper /= y.upper;
      }
      else if (gt(*this, 0.0)) {
	// +x+ +y+
	lower /= y.upper;
	upper /= y.lower;
      }
      else {
	// -x+ +y+
	lower /= y.lower;
	upper /= y.lower;
      }
    }
    round_restore();
    set_maybe_real();
  }
  return *this;
}


// Power
// WARNING: no support for negative base
/**
            +-------------------------------+
            |               y               |
            +---------------+---------------+
            |-inf(odd),0    |0,+inf(odd)    |
+-+---------+---------------+---------------+
|x|-inf,  -1|x- y- (-1,0)   |x+ y- (-inf,-1)|
| |  -1,   0|x- y+ (-inf,-1)|x+ y+ (-1,0)   |
+-+---------+---------------+---------------+

            +-------------------------------+
            |               y               |
            +---------------+---------------+
            |-inf(even),0   |0,+inf(even)   |
+-+---------+---------------+---------------+
|x|-inf,  -1|x+ y+ (0,1)    |x- y+ (1,+inf) |
| |  -1,   0|x+ y- (1,+inf) |x- y- (0, 1)   |
+-+---------+---------------+---------------+

            +-------------------------------+
            |               y               |
            +---------------+---------------+
            |-inf,0         |0,+inf         |
+-+---------+---------------+---------------+
|x|   0,   1|x- y- (1,+inf) |x+ y- (0,1)    |
| |   1,+inf|x- y+ (0,1)    |x+ y+ (1,+inf) |
+-+---------+---------------+---------------+
**/
Interval&
Interval::pow_assign(const Interval& yb)
{
  assert(ok() && yb.ok());
  // WARNING: no handling of signed zeroes
  // To handle correctly a.pow_assign(a)
  const Interval& y = (this == &yb) ? Interval(yb) : yb;
  bool zero;
  zero = (*this && 0.0);
  if (zero && !gt(y, 0.0))
    set_full();
  else {
#if 1
    if (y.is_one_point()) {
      if (y.is_integer()) {
	if (fmod(y.lower.boundary(), 2.0) == 0.0)
	  abs_assign();
      }
      else {
	if (lower < 0.0f) {
	  lower.set(0.0f);
	  if (is_empty()) {
	    set_empty();
	    return *this;
	  }
	}
	if (y.lower < 0.0f) {
	  LBoundary lb(lower);
	  lower.set(upper);
	  upper.set(lb);
	}
      }
      round_save();
      lower.pow_assign(y.lower);
      upper.pow_assign(y.lower);
      round_restore();
      return *this;
    }
#endif
    warn_check(!ge(*this, 0.0), "pow(x, y) with x < 0 NOT YET IMPLEMENTED",
	       set_full(); return *this);
    warn_check(!ge(y, 0.0), "pow(x, y) with y < 0 NOT YET IMPLEMENTED",
	       set_full(); return *this);
    round_save();
    if (gt(*this, 1.0)) {
      lower.pow_assign(y.lower);
      upper.pow_assign(y.upper);
    }
    else if (lt(*this, 1.0)) {
      lower.pow_assign(y.upper);
      upper.pow_assign(y.lower);
    }
    else {
      lower.pow_assign(y.upper);
      upper.pow_assign(y.upper);
    }
    round_restore();
    if (flag != INTEGER || y.flag != INTEGER || !ge(y,0.0))
      set_maybe_real();
    else
      set_integer();	// Needed: pow(2,3) in some library don't give 8
  }
  return *this;
}

// Cosine
Interval&
Interval::cos_assign()
{
  assert(ok());
  round_save();
  bool tomin = contains_multiple(M_PI, M_PI*2.0);
  bool tomax = contains_multiple(0.0, M_PI*2.0);
  if (tomin && tomax) {
      lower.set(-1.0f);
      upper.set(1.0f);
  }
#if IEEE_NUMBERS
  else if (lower.is_min() || upper.is_max())
    set_empty();
#endif
  else {
#if USE_DIRECTED_ROUNDING
    LBoundary lb(lower);
    UBoundary ub(upper);
#endif
    lower.cos_assign();
    upper.cos_assign();
    if (lower > upper) {
#if USE_DIRECTED_ROUNDING
      if (!tomin) {
	lower.set(ub);
	lower.cos_assign();
      }
      if (!tomax) {
	upper.set(lb);
	upper.cos_assign();
      }
#else
      LBoundary lb(lower);
      if (!tomin)
	lower.set(upper);
      if (!tomax)
	upper.set(lb);
#endif
    }
    if (tomax)
      upper.set(1.0f);
    if (tomin)
      lower.set(-1.0f);
  }
  round_restore();
  set_maybe_real();
  return *this;
}

// Sine
Interval&
Interval::sin_assign()
{
  assert(ok());
  round_save();
  bool tomin = contains_multiple(-M_PI_2, M_PI*2.0);
  bool tomax = contains_multiple(M_PI_2, M_PI*2.0);
  if (tomin && tomax) {
      lower.set(-1.0f);
      upper.set(1.0f);
  }
#if IEEE_NUMBERS
  else if (lower.is_min() || upper.is_max())
    set_empty();
#endif
  else {
#if USE_DIRECTED_ROUNDING
    LBoundary lb(lower);
    UBoundary ub(upper);
#endif
    lower.sin_assign();
    upper.sin_assign();
    if (lower > upper) {
#if USE_DIRECTED_ROUNDING
      if (!tomin) {
	lower.set(ub);
	lower.sin_assign();
      }
      if (!tomax) {
	upper.set(lb);
	upper.sin_assign();
      }
#else
      LBoundary lb(lower);
      if (!tomin)
	lower.set(upper);
      if (!tomax)
	upper.set(lb);
#endif
    }
    if (tomax)
      upper.set(1.0f);
    if (tomin)
      lower.set(-1.0f);
  }
  round_restore();
  set_maybe_real();
  return *this;
}

// Tangent
Interval&
Interval::tan_assign()
{
  assert(ok());
  round_save();
  if (contains_multiple(M_PI_2, M_PI))
    set_full();
#if IEEE_NUMBERS
  else if (lower.is_min() || upper.is_max())
    set_empty();
#endif
  else {
    lower.tan_assign();
    upper.tan_assign();
  }
  round_restore();
  set_maybe_real();
  return *this;
}

// Output operator
ostream&
operator << (ostream& s, const Interval& x)
{
  if (x.is_empty())
    return s << "[]";
  else if (x.is_one_point())
    return s << x.lower.value;
  else
    return (s << x.lower
	    << (x.flag == Interval::REAL ? ", " : " .. ")
	    << x.upper);
}

// Input operator
istream&
operator >> (istream& s, Interval& x)
{
  char l[256], u[256];
  LBoundary::Open_Closed lf;
  UBoundary::Open_Closed uf;
  Interval::Type nf;
  char ch;
  s >> ws;
  s >> ch;
  if (ch == '(')
    lf = LBoundary::OPEN;
  else if (ch == '[')
    lf = LBoundary::CLOSED;
  else {
    s.putback(ch);
    s >> l;
    x = Interval(l);
    return s;
  }
  s >> l;
  s >> ws;
  s >> ch;
  if (ch == '.') {
    s >> ch;
    if (ch != '.') {
      //      s.clear(_bad);
      return s;
    }
    nf = Interval::INTEGER;
    s >> u;
    s >> ws;
    s >> ch;
  }
  else if (ch == ',') {
    nf = Interval::REAL;
    s >> u;
    s >> ws;
    s >> ch;
  }
  else {
#if 0
    // FIXME: ISO C++ forbids assignment of arrays
    u = l;
#endif
    nf = Interval::REAL;
  }
  if (ch == ')')
    uf = UBoundary::OPEN;
  else if (ch == ']')
    uf = UBoundary::CLOSED;
  else {
    //    s.clear(_bad);
    return s;
  }
  x = Interval(LBoundary(l, lf), UBoundary(u, uf), nf);
  return s;
}

#if 0
//  sz-1
// (sum *vars[i]*coeffs[i]) + k (rel) 0
//  i=0
void linear_refine(const Interval& k, Rel rel, int sz,
		   Boundary::Value* coeffs, Interval** vars)
{
  assert(!k.is_empty());
  int i;
  Interval p[sz];
  for (i = 0; i < sz; ++i) {
    assert(coeffs[i] != 0.0);
    assert(!vars[i]->is_empty());
    p[i] = *vars[i];
    p[i] *= coeffs[i];
  }
  for (i = 0; i < sz; ++i) {
    Interval q(k);
    int j;
    for (j = 0; j < sz; ++j) {
      if (i == j)
	continue;
      q += p[j];
    }
    q /= -coeffs[i];
    if (coeffs[i]>0)
      refine1(*vars[i], rel, q);
    else
      refine2(q, rel, *vars[i]);
    if (vars[i]->is_empty()) {
      for (i = 0; i < sz; ++i)
	vars[i]->set_empty();
      break;
    }
  }
  round_restore();
}
#endif
