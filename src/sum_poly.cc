/* Definition of some utility functions.
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




/*!
  This file contains the routines that sum formally expressions of the 
  form \f$ p(n) x^n \f$, where \f$ p \f$ is a polynomial (possibly 
  constant) and \f$ x \f$ is a real number (possibly 1) or a 
  <CODE> GExpr </CODE>. 
  The sum is always over the range \f$ [0, N ]\f$. 
*/

#include "sum_poly.hh"

/*!
  This routine computes the falling product 
  \f$ x \cdot (x - 1) \cdots (x - k + 1) \f$. 
  If \f$ k\le 0 \f$, the routine returns 1, which is the usual 
  convention for an empty product. 
  The result is returned in the expression <CODE> q </CODE>. 
*/

static void
falling_product(const GExpr& x, unsigned k, GExpr& q) {
  q = 1;
  for (unsigned i = 0; i < k; ++i) 
    q *= (x-i);
}

/*!
  This routine computes \f$\sum_{j=0}^n j_{(k)} x^j \f$ for a non negative 
  integer \f$ k \f$ by means of a formula explained in purrs.tex, \S4.3.2. 
  The closed formula is returned in the expression <CODE> q </CODE>. 
*/

static void 
sum_falling_prod_times_exp(const GSymbol& n, GNumber k, 
			   const GSymbol& x, GExpr& q) {
  
  q = pow(1-x, -k-1);
  GExpr r;
  for (unsigned i = 0; i < k + 1; ++i) {
    falling_product(n + 1, i, r);
    q -= r / factorial(i) * pow(x, n+1-i) * pow(1-x, i-k-1);
  }
  q *= factorial(k) * pow(x, k);
  q.expand();
}

/*!
  This routine computes the polynomial decomposition 
  \f$ p(x) = \sum_{k=0}^d b_k x_{(k)} \f$ where \f$ d \f$ is the degree 
  of the polynomial \f$ p \f$, as explained in purrs.tex, \S4.3.2. 
  The coefficients \f$ b_0\f$, \f$\dots\f$, \f$ b_d \f$ are stored 
  in the vector <CODE> summands </CODE>. 
*/

static void 
poly_dec(const GExpr& p, const GSymbol& x, vector<GExpr>& summands) {

  unsigned d = p.degree(x);
  GExpr q = p;
  GExpr r = p.coeff(x,0);
  summands[0] = r;
  for (unsigned i = 0; i < d; ) {
    q -= r;
    q /= x - i;
    ++i;
    r = q.subs(x == i);
    summands[i] = r;
  }
}

/*!
  This routine computes \f$\sum_{j=0}^n p(j) \f$ for a non 
  negative integer \f$ n \f$ and a polynomial \$ p \f$ by means 
  of a formula explained in purrs.tex, \S4.3.3. 
  The closed formula is returned in the expression <CODE> q </CODE>. 
*/

static void 
sum_poly(const GExpr& p, const GSymbol& x, const GSymbol& n, GExpr& q) {
  
  unsigned d = p.degree(x);
  vector<GExpr> summands(d+1);
  poly_dec(p, x, summands);
  q = 0;
  for (unsigned i = 0; i <= d; ++i) {
    GExpr r;
    falling_product(n+1, i+1, r);
    q += numeric(1, i+1) * summands[i] * r;
  }
  q.expand();
}

/*!
  At this point we can sum exactly any linear combination of 
  products of polynomial and exponentials (including extreme 
  cases of constant polynomials or exponentials). 
*/

/*!
  This routine computes the closed formula for 
  \f$\sum_{j=0}^n p(j) \alpha^j \f$, 
  where \f$ p \f$ is a polynomial (possibly constant) and 
  \f$ \alpha \f$ is a GNumber (possibly 1). 
  The closed formula is returned in the expression <CODE> q </CODE>. 
*/

void
sum_poly_times_exponentials(const GExpr& p, const GSymbol& x, 
			    const GSymbol& n, const GExpr& alpha, GExpr& q) {
  if (alpha == 1) 
    sum_poly(p, x, n, q);
  // we just have to compute the sum of the values of the polynomial 
  else {
    GExpr r;
    unsigned d = p.degree(x);
    vector<GExpr> summands(d+1);
    poly_dec(p, x, summands);
    q = 0;
    for (unsigned i = 0; i <= d; ++i) {
      sum_falling_prod_times_exp(n, i, x, r);
      r.expand();
      q += r * pow(x, i) * summands[i];
    }
  }
  q = expand(q.subs(x == alpha));
}
