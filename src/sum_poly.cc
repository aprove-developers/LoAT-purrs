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


/*!
  This file contains the routines that sum formally expressions of the 
  form \f$ p(n) x^n \f$, where \f$ p \f$ is a polynomial (possibly 
  constant) and \f$ x \f$ is a real number (possibly 1) or an
  <CODE> Expr </CODE>.
  The sum is always over the range \f$ n \in [0, N] \f$.
*/

#include "sum_poly.hh"
#include "Expr.defs.hh"
#include "Recurrence.defs.hh"
#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

/*!
  This routine computes the falling product 
  \f$ x_{(k)} := x \cdot (x - 1) \cdots (x - k + 1) \f$,
  where \f$ k \f$ is an integer. 
  If \f$ k\le 0 \f$, the routine returns 1, which is the usual 
  convention for an empty product. 
  The result is returned in the expression <CODE>q</CODE>. 
*/

void
falling_product(const Expr& x, const Number& k, Expr& q) {
  q = 1;
  for (Number i = 0; i < k; ++i) 
    q *= (x-i);
}

/*!
  This routine computes \f$\sum_{j=0}^n j_{(k)} x^j \f$ for a non negative 
  integer \f$ k \f$ by means of a formula explained in purrs.tex, \S4.3.2. 
  The closed formula is returned in the expression <CODE>q</CODE>. 
*/

void
sum_falling_prod_times_exp(const Number& k, const Symbol& x, Expr& q) {

  q = pwr(1-x, -k-1);
  Expr r;
  for (Number i = 0; i <= k; ++i) {
    falling_product(Recurrence::n + 1, i, r);
    q -= r / factorial(i) * pwr(x, Recurrence::n + 1 - i) * pwr(1 - x, i - k - 1);
  }
  q *= factorial(k) * pwr(x, k);
}

/*!
  This routine computes the polynomial decomposition 
  \f$ p(x) = \sum_{k=0}^d b_k x_{(k)} \f$ where \f$ d \f$ is the degree 
  of the polynomial \f$ p \f$, as explained in purrs.tex, \S4.3.2. 
  The coefficients \f$ b_0\f$, \f$\dots\f$, \f$ b_d \f$ are stored 
  in the vector <CODE> summands </CODE>. 
*/

void 
poly_dec(const Expr& p, const Symbol& x, std::vector<Expr>& summands) {

  unsigned d = p.degree(x);
  Expr q = p;
  Expr r = p.coeff(x,0);
  summands[0] = r;
  for (unsigned i = 0; i < d; ) {
    q -= r;
    q /= x - i;
    ++i;
    r = q.subs(x, i);
    summands[i] = r;
  }
}

/*!
  This routine computes \f$\sum_{j=0}^n p(j) \f$ for a non 
  negative integer \f$ n \f$ and a polynomial \f$ p \f$ by means 
  of a formula explained in purrs.tex, \S4.3.3. 
  The closed formula is returned in the expression <CODE> q </CODE>. 
*/

void 
sum_poly(const Expr& p, const Symbol& x, Expr& q) {
  
  unsigned d = p.degree(x);
  std::vector<Expr> summands(d+1);
  poly_dec(p, x, summands);
  q = 0;
  for (unsigned i = 0; i <= d; ++i) {
    Expr r;
    falling_product(Recurrence::n + 1, i + 1, r);
    q += Number(1, i + 1) * summands[i] * r;
  }
}

} // anonymous namespace

/*!
  At this point we can sum exactly any linear combination of 
  products of polynomials and exponentials (including extreme 
  cases of constant polynomials or exponentials). 
*/

/*!
  This routine computes the closed formula for 
  \f$\sum_{j=0}^n p(j) \alpha^j \f$, 
  where \f$ p \f$ is a polynomial (possibly constant) and 
  \f$ \alpha \f$ is a Number (possibly 1). 
  The closed formula is returned in the expression <CODE> q </CODE>. 
*/

PURRS::Expr 
PURRS::sum_poly_times_exponentials(const Expr& p, const Symbol& x, 
				   const Expr& alpha) {

  Expr q;
  if (alpha == 1)
    sum_poly(p, x, q);
  // we just have to compute the sum of the values of the polynomial 
  else {
    Expr r;
    unsigned d = p.degree(x);
    std::vector<Expr> summands(d+1);
    poly_dec(p, x, summands);
    q = 0;
    for (unsigned i = 0; i <= d; ++i) {
      sum_falling_prod_times_exp(i, x, r);
      r = r.expand();
      q += r * summands[i];
    }
  }
  q = q.subs(x, alpha).expand();
  return q;
}

/*!
  This routine computes the closed formula for 
  \f$\sum_{j=0}^n p(j) \alpha^j \cos(j\theta) \f$, 
  where \f$ p \f$ is a polynomial (possibly constant),  
  \f$ \alpha \f$ and \f$ \theta \f$ are <CODE> Expr </CODE>.
*/

PURRS::Expr
PURRS::sum_poly_times_exponentials_times_cos(const Expr& p, const Symbol& x, 
					     const Expr& alpha, const Expr& theta) {
  Expr q = 0;
  if (theta.is_zero()) {
    q = sum_poly_times_exponentials(p, x, alpha);
    return q;
  } 
  unsigned d = p.degree(x);
  std::vector<Expr> summands(d+1);
  poly_dec(p, x, summands);
  q = 0;
  for (unsigned i = 0; i <= d; ++i) {
    Expr r = pwr(x, Recurrence::n + 2) * cos(Recurrence::n * theta)
      - pwr(x, Recurrence::n + 1) * cos((Recurrence::n + 1) * theta);
    r += 1 - x * cos(theta);
    r /= pwr(x,2) - 2* x * cos(theta) + 1;
    r = pwr(x,i) * r.diff(x,i);
    r = r.expand();
    q += r * summands[i];
  }
  q = q.subs(x, alpha).expand();
  return q;
}

/*!
  This routine computes the closed formula for 
  \f$\sum_{j=0}^n p(j) \alpha^j \sin(j\theta) \f$, 
  where \f$ p \f$ is a polynomial (possibly constant),  
  \f$ \alpha \f$ and \f$ \theta \f$ are <CODE> Expr </CODE>.
*/

PURRS::Expr
PURRS::sum_poly_times_exponentials_times_sin(const Expr& p, const Symbol& x, 
					     const Expr& alpha, const Expr& theta) {
  Expr q = 0;
  if (theta.is_zero()) {
    return q;
  } 
  unsigned d = p.degree(x);
  std::vector<Expr> summands(d+1);
  poly_dec(p, x, summands);
  q = 0;
  for (unsigned i = 0; i <= d; ++i) {
    Expr r = pwr(x, Recurrence::n + 2) * sin(Recurrence::n * theta)
      - pwr(x, Recurrence::n + 1) * sin((Recurrence::n + 1) * theta);
    r -= x * sin(theta);
    r /= pwr(x,2) - 2* x * cos(theta) + 1;
    r = pwr(x,i) * r.diff(x,i);
    r = r.expand();
    q += r * summands[i];
  }
  q = q.subs(x, alpha).expand();
  return q;
}
