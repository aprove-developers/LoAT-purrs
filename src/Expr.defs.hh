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

#include "Expr.types.hh"
#include "Symbol.types.hh"
#include "Expr_List.types.hh"
#include "Number.types.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

// binary arithmetic operators Expr with Expr
Expr operator+(const Expr& lh, const Expr& rh);
Expr operator-(const Expr& lh, const Expr& rh);
Expr operator*(const Expr& lh, const Expr& rh);
Expr operator/(const Expr& lh, const Expr& rh);

// binary arithmetic assignment operators with Expr
Expr& operator+=(Expr& lh, const Expr& rh);
Expr& operator-=(Expr& lh, const Expr& rh);
Expr& operator*=(Expr& lh, const Expr& rh);
Expr& operator/=(Expr& lh, const Expr& rh);

class Parma_Recurrence_Relation_Solver::Expr {
public:
  //! Default constructor.
  Expr();

  //! Builds the integer expression \p i.
  Expr(int i);

  //! Copy-constructor.
  Expr(const Expr& exp);

  //! Destructor.
  ~Expr();

  //! Assignment operator.
  Expr& operator=(const Expr& exp);

   bool is_a_symbol() const;
  bool is_a_number() const;
  bool is_a_constant() const;
  bool is_a_add() const;
  bool is_a_mul() const;
  bool is_a_power() const;
  bool is_a_function() const;
  bool is_a_matrix() const;
  bool is_a_Expr_List() const;
  bool is_a_relational() const;

  Number ex_to_number();

  // info
  bool is_integer_polynomial() const;
  bool is_rational_polynomial() const;

  unsigned nops() const;
  Expr op(unsigned i) const;
  bool is_equal(const Expr& e) const;
  bool is_zero() const;
  Expr subs(const Expr& e) const;
  Expr subs(const Expr_List& symbols, const Expr_List& replacements) const;
  bool match(const Expr& pattern) const;
  bool match(const Expr& pattern, Expr_List& replacements) const;
  bool has(const Expr& pattern) const;

  Expr expand() const;
  Expr collect(const Expr& lst) const;
  int degree(const Symbol& symb) const;
  int ldegree(const Symbol& symb) const;
  Expr coeff(const Symbol& symb, int k) const;
  Expr lcoeff(const Symbol& symb) const;
  Expr tcoeff(const Symbol& symb) const;
  Expr primpart(const Symbol& symb) const;
  Expr numer() const;
  Expr denom() const;
  Expr numer_denom() const;

  Expr to_rational(Expr_List& lst);

  // solve(), lsolve()
  
private:
  GiNaC::ex e;

  friend class Number;

  Expr(const GiNaC::ex& ge);
};

Expr pow(const Expr& b, const Expr& e);
Expr pow(const Symbol& b, const unsigned i);

Expr quo(const Expr& a, const Expr& b, const Symbol& x);
Expr quo(const Expr& a, const Symbol& x, const Symbol& y);
Expr rem(const Expr& a, const Expr& b, const Symbol& x);
Expr prem(const Expr& a, const Expr& b, const Symbol& x);
Expr gcd(const Expr& a, const Expr& b);
Expr lcm(const Expr& a, const Expr& b);
Expr sqrfree(const Expr& e, const Expr_List& lst);

Expr lsolve(const Expr_List& lst1, const Expr_List& lst2);

} // namespace Parma_Recurrence_Relation_Solver

#include "Expr.inlines.hh"

#endif // !defined(PURRS_Expr_defs_hh)
