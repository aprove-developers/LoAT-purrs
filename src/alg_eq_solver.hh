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

#include "Expr.defs.hh"
#include "Number.defs.hh"
#include <vector>

namespace Parma_Recurrence_Relation_Solver {

enum Number_Type {
  /*!
    The expression is a rational number.
  */
  RATIONAL,

  /*!
    The expression is not a rational number.
  */
  NON_RATIONAL
};

class Polynomial_Root {
private:
  //! Value of the root.
  Expr value_;

  //! Type of the root: <CODE>RATIONAL</CODE> or <CODE>NON_RATIONAL</CODE>.
  Number_Type type_;

  //! Multiplicity of the root.
  unsigned multiplicity_;

public:
  Polynomial_Root(const Expr& e, Number_Type t, unsigned m = 1)
    : value_(e), type_(t), multiplicity_(m) {
  }

  const Expr& value() const {
    return value_;
  }

  Expr& value() {
    return value_;
  }

  unsigned multiplicity() const {
    return multiplicity_;
  }

  void set_multiplicity(unsigned m) {
    multiplicity_ = m;
  }

  Number_Type type() const {
    return type_;
  }

  bool is_rational() const {
    return type() == RATIONAL;
  }

  bool is_non_rational() const {
    return type() == NON_RATIONAL;
  }

  void set_type(Number_Type t) {
    type_ = t;
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
find_roots(const Expr& p, const Symbol& x,
	   std::vector<Polynomial_Root>& roots, bool& all_distinct);

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_alg_eq_solver_hh)
