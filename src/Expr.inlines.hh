/* *************************: inline functions.
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

#ifndef PURRS_Expr_inlines_hh
#define PURRS_Expr_inlines_hh

#include "Number.defs.hh"
#include "Expr_List.defs.hh"
#include "Symbol.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

Expr
operator+(const Expr& lh, const Expr& rh) {
  return lh + rh;
}

Expr
operator-(const Expr& lh, const Expr& rh) {
  return lh - rh;
}

Expr
operator*(const Expr& lh, const Expr& rh) {
  return lh * rh;
}

Expr
operator/(const Expr& lh, const Expr& rh) {
  return lh / rh;
}

Expr&
operator+=(const Expr& lh, const Expr& rh) {
  return lh += rh;
}

Expr&
operator-=(const Expr& lh, const Expr& rh) {
  return lh -= rh;
}

Expr&
operator*=(const Expr& lh, const Expr& rh) {
  return lh *= rh;
}

Expr&
operator/=(const Expr& lh, const Expr& rh) {
  return lh /= rh;
}

inline
Expr::Expr() {
}

inline
Expr::Expr(int i)
  : e(i) {
}

inline
Expr::Expr(const Expr& exp)
  : e(exp.e) {
};

inline Expr&
Expr::operator=(const Expr& exp) {
  e = exp.e;
  return *this;
};

inline
Expr::Expr(const GiNaC::ex& ge)
  : e(ge) {
}

inline
Expr::~Expr() {
}

inline bool
Expr::is_a_symbol() const {
  return GiNaC::is_a<GiNaC::symbol>(e);
}

inline bool
Expr::is_a_number() const {
  return GiNaC::is_a<GiNaC::numeric>(e);
}

inline bool
Expr::is_a_constant() const {
  return GiNaC::is_a<GiNaC::constant>(e);
}

inline bool
Expr::is_a_add() const {
  return GiNaC::is_a<GiNaC::add>(e);
}

inline bool
Expr::is_a_mul() const {
  return GiNaC::is_a<GiNaC::mul>(e);
}

inline bool
Expr::is_a_power() const {
  return GiNaC::is_a<GiNaC::power>(e);
}

inline bool
Expr::is_a_function() const {
  return GiNaC::is_a<GiNaC::function>(e);
}

inline bool
Expr::is_a_matrix() const {
  return GiNaC::is_a<GiNaC::matrix>(e);
}

inline bool
Expr::is_a_Expr_List() const {
  return GiNaC::is_a<GiNaC::lst>(e);
}

inline bool
Expr::is_a_relational() const {
  return GiNaC::is_a<GiNaC::relational>(e);
}

inline Number
Expr::ex_to_number() {
  return GiNaC::ex_to<GiNaC::numeric>(e);
}

//info
inline bool
Expr::is_integer_polynomial() const {
  return e.info(GiNaC::info_flags::integer_polynomial);
}
//info
inline bool
Expr::is_rational_polynomial() const {
  return e.info(GiNaC::info_flags::rational_polynomial);
}

inline unsigned
Expr::nops() const {
  return e.nops();
}

inline Expr
Expr::op(unsigned i) const {
  return e.op(i);
}

inline bool
Expr::is_equal(const Expr& exp) const {
  return e.is_equal(exp.e);
}

inline bool
Expr::is_zero() const {
  return e.is_zero();
}

inline Expr
Expr::subs(const Expr& exp) const {
  return e.subs(exp.e);
}

inline Expr
Expr::subs(const Expr_List& symbols,
	   const Expr_List& replacements) const {
  return e.subs(symbols.l, replacements.l);
}

inline bool
Expr::match(const Expr& pattern) const {
  return e.match(pattern.e);
}

inline bool
Expr::match(const Expr& pattern, Expr_List& replacements) const {
  return e.match(pattern.e, replacements.l);
}

inline bool
Expr::has(const Expr& pattern) const {
  return e.has(pattern.e);
}

inline Expr
Expr::expand() const {
  return e.expand();
}

inline Expr
Expr::collect(const Expr& exp) const {
  return e.collect(exp.e);
}

inline int
Expr::degree(const Symbol& symb) const {
  return e.degree(symb.s);
}

inline int
Expr::ldegree(const Symbol& symb) const {
  return e.ldegree(symb.s);
}

inline Expr
Expr::coeff(const Symbol& symb, int k) const {
  return e.coeff(symb.s, k);
}

inline Expr
Expr::lcoeff(const Symbol& symb) const {
  return e.lcoeff(symb.s);
}

inline Expr
Expr::tcoeff(const Symbol& symb) const {
  return e.tcoeff(symb.s);
}

inline Expr
Expr::primpart(const Symbol& symb) const {
  return e.primpart(symb.s);
}

inline Expr
Expr::numer() const {
  return e.numer();
}

inline Expr
Expr::denom() const {
  return e.denom();
}

inline Expr
Expr::numer_denom() const {
  return e.numer_denom();
}

inline Expr
Expr::to_rational(Expr_List& lst) {
  return e.to_rational(lst.l);
}

inline Expr
pow(const Expr& b, const Expr& e) {
  return pow(b, e);
}

inline Expr
pow(const Symbol& b, const unsigned i) {
  return pow(b, i);
}

inline Expr
quo(const Expr& a, const Expr& b, const Symbol& x) {
  return quo(a, b, x);
}

inline Expr
quo(const Expr& a, const Symbol& x, const Symbol& y) {
  return quo(a, x, y);
}

inline Expr
rem(const Expr& a, const Expr& b, const Symbol& x) {
  return rem(a, b, x);
}

inline Expr
prem(const Expr& a, const Expr& b, const Symbol& x) {
  return prem(a, b, x);
}

inline Expr
gcd(const Expr& a, const Expr& b) {
  return gcd(a, b);
}

inline Expr
lcm(const Expr& a, const Expr& b) {
  return lcm(a, b);
}

inline Expr
sqrfree(const Expr& exp, const Expr_List& lst) {
  return sqrfree(exp, lst);
}

inline Expr
lsolve(const Expr_List& lst1, const Expr_List& lst2) {
  return lsolve(lst1, lst2);
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Expr_inlines_hh)
