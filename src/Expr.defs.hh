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
#include "Constant.types.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

// output stream operators
std::ostream& operator<<(std::ostream& os, const Expr& exp);

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

Expr operator+(const Expr &lh);
Expr operator-(const Expr &lh);

class Parma_Recurrence_Relation_Solver::Expr {
public:
  //! Default constructor.
  Expr();

  //! Builds the integer expression \p i.
  Expr(int i);

  //! Builds the numeric expression \p n.
  Expr(const Number& n);

  //! Builds the symbolic expression \p s.
  Expr(const Symbol& s);

  //! Builds the constant expression \p k.
  Expr(const Constant& k);

  //! Builds the list expression \p lst.
  Expr(const Expr_List& lst);

  //! Builds the list expression \p lst.
  Expr(const std::string& st, const Expr_List& lst);

  //! Copy-constructor.
  Expr(const Expr& exp);

  //! Destructor.
  ~Expr();

  //! Assignment operator.
  Expr& operator=(const Expr& exp);

  Expr operator[](int i) const;

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
  bool is_exactly_a_number() const;
  bool is_exactly_a_constant() const;
  bool is_exactly_a_add() const;
  bool is_exactly_a_mul() const;
  bool is_exactly_a_power() const;
  bool is_exactly_a_function() const;

  // FIXME: const?
  Number ex_to_number() const;

  // info
  bool is_integer_polynomial() const;
  bool is_rational_polynomial() const;
  bool is_relation_equal() const;

  unsigned nops() const;
  Expr op(unsigned i) const;
  bool is_equal(const Expr& e) const;
  bool is_zero() const;
  Expr subs(const Expr& exp1, const Expr& exp2) const;
  Expr subs(const Expr_List& to_replace, const Expr_List& replacements) const;
  bool match(const Expr& pattern) const;
  bool match(const Expr& pattern, Expr_List& replacements) const;
  bool has(const Expr& pattern) const;

  Expr expand() const;
  Expr collect(const Expr_List& lst) const;
  int degree(const Symbol& symb) const;
  int ldegree(const Symbol& symb) const;
  Expr coeff(const Symbol& symb, int k) const;
  Expr lcoeff(const Symbol& symb) const;
  Expr tcoeff(const Symbol& symb) const;
  Expr primpart(const Symbol& symb) const;
  Expr numer() const;
  Expr denom() const;
  Expr numer_denom() const;
  Expr lhs() const;
  Expr rhs() const;

  Expr diff(const Symbol& symb, unsigned nth = 1);

  Expr to_rational(Expr_List& lst);

  // solve(), lsolve()
  
private:
  GiNaC::ex e;

  friend class Number;
  friend class Constant;
  friend class Expr_List;
  friend class Matrix;

  Expr(const GiNaC::ex& ge);
};

// FIXME: meglio con argomento di default `unsigned label = 0'?
Expr wild(unsigned label);

Expr pow(const Expr& b, const Expr& e);
Expr pow(const Symbol& b, const unsigned i);
Expr sqrt(const Expr& e);
Expr sin(const Expr& e);
Expr cos(const Expr& e);
Expr acos(const Expr& e);
Expr tan(const Expr& e);
Expr exp(const Expr& e);
Expr ln(const Expr& e);
Expr factorial(const Expr& e);
// FIXME: ??
Expr factorial(const unsigned i);

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
