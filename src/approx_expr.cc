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

#define FILIB_NAMESPACES 1

#include "globals.hh"
#include "Interval.h"
//#include <cmath>
#include <complex>
#include <cln/rational.h>
#include <ginac/ginac.h>

typedef filib::Interval Interval;
typedef std::complex<Interval> CInterval;

#if 0
template<typename _Tp>
inline complex<_Tp>
tan(const complex<_Tp>& __z) {
  return sin(__z) / cos(__z);
}
#endif

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
  else if (is_exactly_a<power>(e)) {
    static GExpr one_half = GNumber(1)/2;
    const GExpr& base = e.op(0);
    const GExpr& exponent = e.op(1);
#if 0
    if (exponent == one_half)
      return sqrt(approximate(base));
    else
#endif
      return pow(approximate(base), approximate(exponent));
  }
  else if (is_exactly_a<function>(e)) {
    const GExpr& arg = e.op(0);
#if 0
    if (is_ex_the_function(e, abs))
      return abs(approximate(arg));
    else
#endif
    if (is_ex_the_function(e, exp))
      return exp(approximate(arg));
#if 0
    else if (is_ex_the_function(e, log))
      return log(approximate(arg));
#endif
    else if (is_ex_the_function(e, sin))
      return sin(approximate(arg));
    else if (is_ex_the_function(e, cos))
      return cos(approximate(arg));
    else if (is_ex_the_function(e, tan))
      return std::tan(approximate(arg));
    else
      abort();
  }
  else if (is_exactly_a<constant>(e)) {
    if (e == Pi)
      return Interval::PI();
#if 0
    else if (e == Euler)
      return CInterval(Interval(E_lower_bound.d, E_upper_bound.d), 0);
#endif
    else
      abort();
  }
  else
    abort();

  return r;
}
