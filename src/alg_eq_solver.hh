/* Declaration of the main function of the algebraic equation solver.
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

#ifndef PURRS_alg_eq_solver_hh
#define PURRS_alg_eq_solver_hh 1

#include <vector>
#include "globals.hh"

class Polynomial_Root {
private:
  GExpr value_;
  GNumber multiplicity_;

public:
  Polynomial_Root(const GExpr& e, GNumber m = 1)
    : value_(e), multiplicity_(m) {
  }

  const GExpr& value() const {
    return value_;
  }
  GExpr& value() {
    return value_;
  }
  GNumber multiplicity() const {
    return multiplicity_;
  }
  void set_multiplicity(GNumber m) {
    multiplicity_ = m;
  }
};

std::ostream& operator<<(std::ostream& s, const Polynomial_Root& r);


//! Let \p p be a polynomial with integer coefficients in \p x and
//! \p roots be a (possibly non-empty) vector.
//! This function appends to \p roots some or all the roots of \p p.
//! Let \f$n\f$ be the degree of \p p and let \f$h\f$ and \f$k\f$
//! be the value of <CODE>roots.size()</CODE> on entry and on exit
//! to and from find_roots(), respectively.
//! If \f$n = k-h\f$, then the positions \f$h, h+1, \ldots, k-1\f$
//! of \p roots contain <EM>all</EM> the (possibly complex) roots of \p p
//! and the function returns <CODE>true</CODE>.
//! If \f$n \neq k-h\f$, then the positions \f$h, h+1, \ldots, k-1\f$
//! of \p roots contain <EM>some</EM> of the roots of \p p and the function
//! returns <CODE>true</CODE> if it was possible to prove that <EM>all</EM>
//! the roots of \p p of maximal modulus are among those inserted
//! into \p roots;
//! the function returns <CODE>false</CODE> otherwise.
//! The parameter \p all_distinct is set to <CODE>true</CODE> if it was
//! possible to prove that all the roots of \p p are distinct;
//! \p all_distinct is set to <CODE>false</CODE> otherwise.
bool
find_roots(const GExpr& p, const GSymbol& x,
	   std::vector<Polynomial_Root>& roots, bool& all_distinct);

//! Finds divisors of positive integer \p n.
void
find_divisors(GNumber n, std::vector<GNumber>& divisors);
#endif
