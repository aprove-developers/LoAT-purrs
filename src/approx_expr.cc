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
#include <cln/rational.h>
#include <ginac/ginac.h>
#include "cimath.h"

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
  if (GiNaC::is_exactly_a<GiNaC::numeric>(e))
    return approximate(GiNaC::ex_to<GiNaC::numeric>(e));
  else if (GiNaC::is_exactly_a<GiNaC::add>(e)) {
    r = CInterval(0, 0);
    for (unsigned i = 0, n = e.nops(); i < n; ++i)
      r += approximate(e.op(i));
  }
  else if (GiNaC::is_exactly_a<GiNaC::mul>(e)) {
    r = CInterval(1, 0);
    for (unsigned i = 0, n = e.nops(); i < n; ++i)
      r *= approximate(e.op(i));
  }
  else if (GiNaC::is_exactly_a<GiNaC::power>(e)) {
    static GExpr one_half = GNumber(1)/2;
    const GExpr& base = e.op(0);
    const GExpr& exponent = e.op(1);
#if 0
    return pow(approximate(base), approximate(exponent));
#else
    std::cout << base << "^" << exponent << std::endl;
    CInterval abase = approximate(base);
    CInterval aexponent = approximate(exponent);
    std::cout << abase << "^" << aexponent << " = ";
    CInterval result = pow(abase, aexponent);
    std::cout << result << std::endl;
    return result;
#endif
  }
  else if (GiNaC::is_exactly_a<GiNaC::function>(e)) {
    const GExpr& arg = e.op(0);
    CInterval aarg = approximate(arg);
#if 0
    if (is_ex_the_function(e, abs))
      return abs(aarg);
    else
#endif
    if (is_ex_the_function(e, exp))
      return exp(aarg);
    else if (is_ex_the_function(e, log))
      return ln(aarg);
    else if (is_ex_the_function(e, sin))
      return sin(aarg);
    else if (is_ex_the_function(e, cos))
      return cos(aarg);
    else if (is_ex_the_function(e, tan))
      return tan(aarg);
    else if (is_ex_the_function(e, acos))
      return acos(aarg);
    else
      abort();
  }
  else if (GiNaC::is_exactly_a<constant>(e)) {
    if (e == Pi)
      return CInterval(Interval::PI(), Interval::ZERO());
#if 0
    else if (e == Euler)
      return CInterval(Interval(E_lower_bound.d, E_upper_bound.d),
		       Interval::ZERO());
#endif
    else
      abort();
  }
  else
    abort();

  return r;
}
