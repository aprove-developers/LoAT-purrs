/* To be written.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#ifndef NOISY
#define NOISY 0
#endif

#include "sum_poly.hh"
#include "Expr.defs.hh"
#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

/*!
  This file contains the routines that sum formally expressions of the
  form \f$ p(n) x^n \f$, where \f$ p \f$ is a polynomial (possibly
  constant) and \f$ x \f$ is a real number (possibly 1) or an
  <CODE> Expr </CODE>.
  The sum is always over the range \f$ n \in [0, N] \f$.
*/

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
    q *= x - i;
}

/*!
  This routine computes \f$\sum_{j=0}^N j_{(k)} x^j \f$ for a non negative
  integer \f$ k \f$ by means of a formula explained in purrs.tex,
  section~4.3.2.
  The closed formula is returned in the expression <CODE>q</CODE>.
*/
void
sum_falling_prod_times_exp(const Number& k, const Symbol& x, const Symbol& N,
			   Expr& q) {
  q = pwr(1 - x, - k - 1);
  Expr r;
  for (Number i = 0; i <= k; ++i) {
    falling_product(N + 1, i, r);
    q -= r / factorial(i) * pwr(x, N + 1 - i)
      * pwr(1 - x, i - k - 1);
  }
  q *= factorial(k) * pwr(x, k);
}

/*!
  This routine computes the polynomial decomposition
  \f$ p(x) = \sum_{k=0}^d b_k x_{(k)} \f$ where \f$ d \f$ is the degree
  of the polynomial \f$ p \f$, as explained in purrs.tex, section~4.3.2.
  The coefficients \f$ b_0\f$, \f$\dots\f$, \f$ b_d \f$ are stored
  in the vector <CODE> summands </CODE>.
*/
void
poly_dec(const Expr& p, const Symbol& x, std::vector<Expr>& summands) {
  unsigned int d = p.degree(x);
  Expr q = p;
  Expr r = p.coeff(x, 0);
  summands[0] = r;
  for (unsigned int i = 0; i < d; ) {
    q -= r;
    q /= x - i;
    ++i;
    r = q.substitute(x, i);
    summands[i] = r;
  }
}

/*!
  This routine computes \f$\sum_{j=0}^N p(j) \f$ for a non
  negative integer \f$ N \f$ and a polynomial \f$ p \f$ by means
  of a formula explained in purrs.tex, section~4.3.3.
  The closed formula is returned in the expression <CODE> q </CODE>.
*/
void
sum_poly(const Expr& p, const Symbol& x, const Symbol& N, Expr& q) {
  unsigned int d = p.degree(x);
  std::vector<Expr> summands(d+1);
  poly_dec(p, x, summands);
  q = 0;
  for (unsigned int i = 0; i <= d; ++i) {
    Expr r;
    falling_product(N + 1, i + 1, r);
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
  This routine computes \f$\sum_{j=0}^N p(j) \f$ for a non
  negative integer \f$ N \f$ and a polynomial \f$ p \f$ using
  Levy's algorithm based on iterated symbolic integration.
*/
PURRS::Expr
PURRS::sum_poly_alt(const Expr& p, const Symbol& x, const Symbol& N) {

  unsigned int deg = p.degree(x);
  Number coefficients[deg + 2];
  Expr symbolic_sum = (N + 1) * p.coeff(x, 0);;
  
  coefficients[0] = 0;
  coefficients[1] = 1;
  for (unsigned int i = 1; i < deg + 1; ++i) {
    Number sum = 0;
    // Loop for symbolic integration
    for (int j = i + 2; --j > 1; ) {
      coefficients[j] = i * coefficients[j - 1] / j;
      sum += coefficients[j];
    }
    // Compute the `constant of integration' by fitting the formula for n = 1
    coefficients[1] = 1 - sum;
    coefficients[0] = 0;
    Expr partial_sum = 0;
    for (int j = i + 2; --j > 0; )
      partial_sum += pwr(N, j) * coefficients[j];
    symbolic_sum += partial_sum * p.coeff(x, i);
  } 
  return(symbolic_sum);
}

/*!
  This routine computes the closed formula for
  \f$\sum_{j=0}^N p(j) \alpha^j \f$,
  where \f$ p \f$ is a polynomial (possibly constant) and
  \f$ \alpha \f$ is a Number (possibly 1).
  The closed formula is returned in the expression <CODE> q </CODE>.
*/
PURRS::Expr
PURRS::sum_poly_times_exponentials(const Expr& p, const Symbol& x,
				   const Symbol& N, const Expr& alpha) {
  Expr q;
  if (alpha == 1)
    sum_poly(p.expand(), x, N, q);
  else {
    // We just have to compute the sum of the values of the polynomial.
    Expr r;
    unsigned int d = p.expand().degree(x);
    std::vector<Expr> summands(d+1);
    poly_dec(p.expand(), x, summands);
    q = 0;
    for (unsigned int i = 0; i <= d; ++i) {
      sum_falling_prod_times_exp(i, x, N, r);
      r = r.expand();
      q += r * summands[i];
    }
  }
  q = q.substitute(x, alpha).expand();
  return q;
}

/*!
  This routine computes the closed formula for
  \f$\sum_{j=0}^N p(j) \alpha^j \cos(j\theta) \f$,
  where \f$ p \f$ is a polynomial (possibly constant),
  \f$ \alpha \f$ and \f$ \theta \f$ are <CODE> Expr </CODE>.
*/
PURRS::Expr
PURRS::sum_poly_times_exponentials_times_cos(const Expr& p, const Symbol& x,
					     const Symbol& N,
					     const Expr& alpha,
					     const Expr& theta) {
  Expr q;
  if (theta.is_zero())
    q = sum_poly_times_exponentials(p, x, N, alpha);
  else {
    unsigned int d = p.expand().degree(x);
    std::vector<Expr> summands(d + 1);
    poly_dec(p.expand(), x, summands);
    q = 0;
    for (unsigned int i = 0; i <= d; ++i) {
      Expr r = pwr(x, N + 2) * cos(N * theta)
	- pwr(x, N + 1) * cos((N + 1) * theta);
      r += 1 - x * cos(theta);
      r /= pwr(x,2) - 2* x * cos(theta) + 1;
      r = pwr(x,i) * r.diff(x,i);
      r = r.expand();
      q += r * summands[i];
    }
    q = q.substitute(x, alpha).expand();
  }
  return q;
}

/*!
  This routine computes the closed formula for
  \f$\sum_{j=0}^N p(j) \alpha^j \sin(j\theta) \f$,
  where \f$ p \f$ is a polynomial (possibly constant),
  \f$ \alpha \f$ and \f$ \theta \f$ are <CODE> Expr </CODE>.
*/
PURRS::Expr
PURRS::sum_poly_times_exponentials_times_sin(const Expr& p, const Symbol& x,
					     const Symbol& N,
					     const Expr& alpha,
					     const Expr& theta) {
  Expr q = 0;
  if (!theta.is_zero()) {
    unsigned int d = p.expand().degree(x);
    std::vector<Expr> summands(d + 1);
    poly_dec(p.expand(), x, summands);
    q = 0;
    for (unsigned int i = 0; i <= d; ++i) {
      Expr r = pwr(x, N + 2) * sin(N * theta)
	- pwr(x, N + 1) * sin((N + 1) * theta);
      r -= x * sin(theta);
      r /= pwr(x, 2) - 2* x * cos(theta) + 1;
      r = pwr(x, i) * r.diff(x, i);
      r = r.expand();
      q += r * summands[i];
    }
    q = q.substitute(x, alpha).expand();
  }
  return q;
}
