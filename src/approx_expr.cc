/* To be written.
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

#include <config.h>

#include "globals.hh"
#include "Interval.defs.hh"
#include <cmath>
#include <complex>
#include <cln/rational.h>
#include <ginac/ginac.h>

union usi4_double {
  unsigned short int usi4[4];
  double d;
};

#if __BYTE_ORDER == LITTLE_ENDIAN
static usi4_double PI_lower_bound =  { { 0x2d18, 0x5444, 0x21fb, 0x4009 } };
static usi4_double PI_upper_bound =  { { 0x2d19, 0x5444, 0x21fb, 0x4009 } };

static usi4_double E_lower_bound = { { 0x5769, 0x8b14, 0xbf0a, 0x4005 } };
static usi4_double E_upper_bound = { { 0x576a, 0x8b14, 0xbf0a, 0x4005 } };
#else
static usi4_double PI_lower_bound =  { { 0x4009, 0x21fb, 0x5444, 0x2d18 } };
static usi4_double PI_upper_bound =  { { 0x4009, 0x21fb, 0x5444, 0x2d19 } };

static usi4_double E_lower_bound = { { 0x4005, 0xbf0a, 0x8b14, 0x5769 } };
static usi4_double E_upper_bound = { { 0x4005, 0xbf0a, 0x8b14, 0x576a } };
#endif

typedef std::complex<Interval> CInterval;

using namespace GiNaC;

Interval
approximate_integer(const GNumber& n) {
  // Kludge!!!
  return Interval(to_int(n));
}

Interval
approximate_rational(const GNumber& n) {
  if (n.is_integer())
    return approximate_integer(n);
  else if (n.is_rational())
    return approximate_integer(numer(n)) / approximate_integer(denom(n));
  else
    abort();
}

CInterval
approximate(const GNumber& n) {
  if (n.is_real())
    return CInterval(approximate_rational(real(n)),
		     0);
  else
    return CInterval(approximate_rational(real(n)),
		     approximate_rational(imag(n)));
}

CInterval
approximate(const GExpr& e) {
  CInterval r;
  if (is_exactly_a<numeric>(e))
    return approximate(ex_to<numeric>(e));
  else if (is_exactly_a<add>(e)) {
    r = CInterval(0, 0);
    for (unsigned i = 0, n = e.nops(); i < n; ++i)
      r += approximate(e.op(i));
  }
  else if (is_exactly_a<mul>(e)) {
    r = CInterval(1, 0);
    for (unsigned i = 0, n = e.nops(); i < n; ++i)
      r *= approximate(e.op(i));
  }
#if 0
  else if (is_exactly_a<power>(e)) {
    static GExpr one_half = GNumber(1)/2;
    const GExpr& base = e.op(0);
    const GExpr& exponent = e.op(1);
    if (exponent == one_half)
      return sqrt(approximate(base));
    else
      return pow(approximate(e.op(0)), approximate(e.op(1)));
  }
  else if (is_exactly_a<function>(e)) {
    const GExpr& arg = e.op(0);
    if (is_ex_the_function(e, abs))
      return abs(approximate(arg));
    else if (is_ex_the_function(e, exp))
      return exp(approximate(arg));
    else if (is_ex_the_function(e, log))
      return log(approximate(arg));
    else if (is_ex_the_function(e, sin))
      return sin(approximate(arg));
    else if (is_ex_the_function(e, cos))
      return cos(approximate(arg));
    else if (is_ex_the_function(e, tan))
      return tan(approximate(arg));
    else
      abort();
  }
#endif
  else if (is_exactly_a<constant>(e)) {
    if (e == Pi)
      return CInterval(Interval(PI_lower_bound.d, PI_upper_bound.d), 0);
    else if (e == Euler)
      return CInterval(Interval(E_lower_bound.d, E_upper_bound.d), 0);
    else
      abort();
  }
  else
    abort();

  return r;
}
