/* Finding the xxx factors of a polynomial.
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
#include "poly_factor.hh"
#include <NTL/ZZXFactoring.h>

using namespace NTL;

static long
ZZ_to_long(const ZZ& zz) {
  assert(zz >= LONG_MIN && zz <= LONG_MAX);
  static const int bs = sizeof(unsigned long)/sizeof(unsigned char);
  static unsigned char b[bs];
  BytesFromZZ(b, zz, bs);
  unsigned long r = 0;
  for (unsigned int i = bs-1; i-- > 0; )
    r += b[i] << (i*8);
  return (zz < 0) ? -r : r;
}

int
poly_factor(const GExpr& p, const GSymbol& x, std::vector<GExpr>& factors) {
  assert(p.info(GiNaC::info_flags::integer_polynomial));
  ZZX ntl_p;
  for (int i = p.ldegree(x), d = p.degree(x); i<= d; ++i) {
    GExpr e_i = p.coeff(x, i);
    assert(GiNaC::is_a<GiNaC::numeric>(e_i));
    GNumber a_i = GiNaC::ex_to<GiNaC::numeric>(e_i);
    if (a_i < LONG_MIN || a_i > LONG_MAX)
      return 1;
    SetCoeff(ntl_p, i, to_long(a_i));
  }
  vec_pair_ZZX_long ntl_factors;
  ZZ c;
  factor(c, ntl_factors, ntl_p);
  int num_factors = ntl_factors.length();
  if (num_factors == 1)
    return 1;
  factors.reserve(num_factors);
  for (int k = 0; k < num_factors; ++k) {
    // We assume the input polynomial has no repeated factors.
    assert(ntl_factors[k].b == 1);
    const ZZX& ntl_factor = ntl_factors[k].a;
    GExpr factor = 0;
    long d = deg(ntl_factor);
    // We assume that zero is not a factor of the input polynomial,
    // thus d != -1.  Moreover, the NTL function factor() confines
    // numeric factors into its first parameter (here named `c');
    // this excludes the case d == 0.
    assert(d > 0);
    for (long i = d; i >= 0; --i)
      if (long a_i = ZZ_to_long(coeff(ntl_factor, i)))
	factor += a_i*pow(x, i);
    factors.push_back(factor);
  }
  return num_factors;
}
