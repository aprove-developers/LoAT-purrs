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

#ifndef NOISY
#define NOISY 0
#endif

#include "ep_decomp.hh"

#include "util.hh"
#include "Expr.defs.hh"
#include "Recurrence.defs.hh"

// TEMPORARY
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

void
exp_poly_decomposition_factor(const Expr& base, const Expr& e,
			      std::vector<Expr>& alpha,
			      std::vector<Expr>& p,
			      std::vector<Expr>& q) {
  unsigned alpha_size = alpha.size();
  unsigned position = alpha_size;
  bool found = false;
  for (unsigned i = alpha_size; i-- > 0; )
    if (base == alpha[i]) {
      position = i;
      found = true;
      break;
    }
  if (!found) {
    alpha.push_back(base);
    p.push_back(0);
    q.push_back(0);
  }
  // Here `alpha[position]' contains `base' and the polynomial and
  // possibly not polynomial parts of `e' can be added to
  // `p[position]' and `q[position]', respectively.
  Expr polynomial;
  Expr possibly_non_polynomial;
  isolate_polynomial_part(e, Recurrence::n, polynomial, possibly_non_polynomial);
  p[position] += polynomial;
  q[position] += possibly_non_polynomial;
}

void
exp_poly_decomposition_summand(const Expr& e,
			       std::vector<Expr>& alpha,
			       std::vector<Expr>& p,
			       std::vector<Expr>& q) {
  static Expr exponential = pwr(wild(0), Recurrence::n);
  Expr_List substitution;
  unsigned num_factors = e.is_a_mul() ? e.nops() : 1;
  if (num_factors == 1) {
    if (e.match(exponential, substitution)) {
      // We have found something of the form `power(base, n)'.
      Expr base = get_binding(substitution, 0);
      assert(!base.is_zero());
      if (base.is_scalar_representation(Recurrence::n)) {
	// We have found something of the form `power(base, n)'
	// and `base' is good for the decomposition.
	exp_poly_decomposition_factor(base, 1, alpha, p, q);
	return;
      }
    }
  }
  else
    for (unsigned i = num_factors; i-- > 0; ) {
      if (clear(substitution), e.op(i).match(exponential, substitution)) {
	// We have found something of the form `power(base, n)'.
	Expr base = get_binding(substitution, 0);
	assert(!base.is_zero());
	if (base.is_scalar_representation(Recurrence::n)) {
	  // We have found something of the form `power(base, n)'
	  // and `base' is good for the decomposition: determine
	  // `r = e/power(base, n)'.
	  Expr r = 1;
	  for (unsigned j = num_factors; j-- > 0; )
	    if (i != j)
	      r *= e.op(j);
	  exp_poly_decomposition_factor(base, r, alpha, p, q);
	  return;
	}
      }
    }
  // No proper exponential found: this is treated like `power(1, n)*e'.
  exp_poly_decomposition_factor(1, e, alpha, p, q);
}

} // anonymous namespace

/*!
  Let \f$ e(n) \f$ be the expression in \p n contained in \p e,
  which is assumed to be already expanded.
  This function computes a decomposition
  \f$ e(n) = \sum_{i=0}^k \alpha_i^n \bigl(p_i(n) + q_i(n)\bigr) \f$, where
  - \f$ \alpha_i \f$ is a expression valid for to be an exponential's base.
    (syntactically different from \p 0);
  - \f$ \alpha_i \neq \alpha_j \f$ if \f$ i \neq j \f$;
  - \f$ p_i(n) \f$ is (syntactically) a polynomial in \f$ n \f$.

  The expressions corresponding to \f$ \alpha_i \f$, \f$ p_i \f$ and
  \f$ q_i \f$ are stored in the \f$ i \f$-th position of the vectors
  \p alpha, \p p and \p q, respectively.
*/
void
PURRS::exp_poly_decomposition(const Expr& e,
			      std::vector<Expr>& alpha,
			      std::vector<Expr>& p,
			      std::vector<Expr>& q) {
  unsigned num_summands = e.is_a_add() ? e.nops() : 1;
  // An upper bound to the number of exponentials is the number of
  // summands in `e': reserve space in the output vectors so that
  // no reallocations will be required.
  alpha.reserve(num_summands);
  p.reserve(num_summands);
  q.reserve(num_summands);
  if (num_summands > 1)
    for (unsigned i = num_summands; i-- > 0; )
      exp_poly_decomposition_summand(e.op(i), alpha, p, q);
  else
    exp_poly_decomposition_summand(e, alpha, p, q);
}
