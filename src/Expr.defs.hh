/* *****************
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

#ifndef PURRS_Expr_defs_hh
#define PURRS_Expr_defs_hh 1

#include "Symbol.types.hh"
#include "Expr_List.types.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

class Parma_Recurrence_Relation_Solver::Expr {
public:
  //! Default constructor.
  Expr();

  //! Builds the integer expression \p i.
  Expr(int i);

  //! Copy-constructor.
  Expr(const Expr& x);

  //! Destructor.
  ~Expr();

  //! Assignment operator.
  Expr& operator=(const Expr& x);

  bool is_a_symbol(const Expr& e);
  bool is_a_number(const Expr& e);
  bool is_a_constant(const Expr& e);
  bool is_a_add(const Expr& e);
  bool is_a_mul(const Expr& e);
  bool is_a_power(const Expr& e);
  bool is_a_function(const Expr& e);
  bool is_a_matrix(const Expr& e);
  bool is_a_Expr_List(const Expr& e);
  bool is_a_relational(const Expr& e);

  //bool is_exactly_a< >(const Expr& e);
  ////
  bool info(unsigned flag);

  unsigned nops() const;
  Expr op(unsigned i) const;
  bool is_equal(const Expr& e) const;
  bool is_zero() const;
  Expr subs(const Expr& e) const;
  Expr subs(const Expr_List& symbols, const Expr_List& replacements) const;
  bool match(const Expr& pattern) const;
  bool match(const Expr& pattern, Expr_List& replacements) const;
  bool has(const Expr& pattern) const;
  // lo uso?
  //bool find(const Expr& pattern, Expr& found) const;

  Expr expand();
  //
  Expr collect(const Expr& lst, bool distributed);

  int degree(const Expr& e) const;
  int ldegree(const Expr& e) const;
  Expr coeff(const Expr& e, int k) const;
  Expr lcoeff(const Expr& e) const;
  Expr tcoeff(const Expr& e) const;
  Expr quo(const Expr& a, const Expr& b, const Symbol& x) const;
  Expr rem(const Expr& a, const Expr& b, const Symbol& x) const;
  Expr prem(const Expr& a, const Expr& b, const Symbol& x) const;
  bool divide(const Expr& a, const Expr& b, const Expr& q) const;
  Expr primpart(const Symbol& x);
  Expr gcd(const Expr& a, const Expr& b) const;
  Expr lcm(const Expr& a, const Expr& b) const;
  Expr sqrfree(const Expr& e, const Expr_List& lst) const;
  Expr numer() const;
  Expr denom() const;
  Expr numer_denom() const;

  Expr to_rational(Expr_List& lst);
  //
  int to_int() const;
  long to_long() const;

  // solve(), lsolve()
private:
  GiNaC::ex e;

  Expr(const GiNaC::ex& ge);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Expr.inlines.hh"

#endif // !defined(PURRS_Expr_defs_hh)
