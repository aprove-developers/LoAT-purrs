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
#include "Constant.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

inline std::ostream&
operator<<(std::ostream& os, const Expr& exp) {
  os << exp;
  return os;  
};

inline Expr
operator+(const Expr& lh, const Expr& rh) {
  return lh + rh;
}

inline Expr
operator-(const Expr& lh, const Expr& rh) {
  return lh - rh;
}

inline Expr
operator*(const Expr& lh, const Expr& rh) {
  return lh * rh;
}

inline Expr
operator/(const Expr& lh, const Expr& rh) {
  return lh / rh;
}

inline Expr
operator+(const Expr& lh) {
  return +lh;
}

inline Expr
operator-(const Expr& lh) {
  return -lh;
}

inline Expr&
operator+=(const Expr& lh, const Expr& rh) {
  return lh += rh;
}

inline Expr&
operator-=(const Expr& lh, const Expr& rh) {
  return lh -= rh;
}

inline Expr&
operator*=(const Expr& lh, const Expr& rh) {
  return lh *= rh;
}

inline Expr&
operator/=(const Expr& lh, const Expr& rh) {
  return lh /= rh;
}

inline bool
operator==(const Expr& lh, const Expr& rh) {
  return lh == rh;
}

inline bool
operator!=(const Expr& lh, const Expr& rh) {
  return lh != rh;
}

inline bool
operator<(const Expr& lh, const Expr& rh) {
  return lh < rh;
}

inline bool
operator>(const Expr& lh, const Expr& rh) {
  return lh > rh;
}

inline bool
operator<=(const Expr& lh, const Expr& rh) {
  return lh <= rh;
}

inline bool
operator>=(const Expr& lh, const Expr& rh) {
  return lh >= rh;
}

inline
Expr::Expr() {
}

inline
Expr::Expr(int i)
  : e(i) {
}

inline
Expr::Expr(const Number& n)
  : e(n.n) {
}

inline
Expr::Expr(const Symbol& s)
  : e(s.s) {
}

inline
Expr::Expr(const Constant& k)
  : e(k.c) {
}

#if 0
inline
Expr::Expr(const Expr_List& lst)
  : e(lst.l) {
}
#endif

inline
Expr::Expr(const std::string& st, const Expr_List& lst)
  : e(st, lst.l) {
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

inline Expr
Expr::operator[](int i) const {
  return e[i];
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
Expr::ex_to_number() const {
  assert(GiNaC::is_a<GiNaC::numeric>(e));
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
//info
inline bool
Expr::is_relation_equal() const {
  return e.info(GiNaC::info_flags::relation_equal);
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
Expr::subs(const Expr& exp1, const Expr& exp2) const {
  return e.subs(exp1.e == exp2.e);
}

inline Expr
Expr::subs(const Expr_List& to_replace,
	   const Expr_List& replacements) const {
  return e.subs(to_replace.l, replacements.l);
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
Expr::collect(const Expr_List& lst) const {
  return e.collect(lst.l);
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
Expr::lhs() const {
  return e.lhs();
}

inline Expr
Expr::rhs() const {
  return e.rhs();
}

inline Expr
Expr::to_rational(Expr_List& lst) {
  return e.to_rational(lst.l);
}

inline Expr
Expr::diff(const Symbol& symb, unsigned nth) {
  return e.diff(symb.s, nth);
}

inline Expr
wild(unsigned label) {
  return wild(label);
}

inline Expr
power(const Expr& b, const Expr& e) {
  return GiNaC::pow(b.e, e.e);
}

inline Expr
sqrt(const Expr& e) {
  return GiNaC::sqrt(e.e);
};

inline Expr
sin(const Expr& e) {
  return GiNaC::sin(e.e);
};

inline Expr
cos(const Expr& e) {
  return GiNaC::cos(e.e);
};

inline Expr
acos(const Expr& e) {
  return GiNaC::acos(e.e);
};

inline Expr
tan(const Expr& e) {
  return GiNaC::tan(e.e);
};

inline Expr
exp(const Expr& e) {
  return GiNaC::exp(e.e);
};

inline Expr
log(const Expr& e) {
  return GiNaC::log(e.e);
};

inline Expr
quo(const Expr& a, const Expr& b, const Symbol& x) {
  return quo(a, b, x);
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
