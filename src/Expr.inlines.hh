/* Expr class implementation: inline functions.
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

#include <stdexcept>

namespace GiNaC {

ex x_eval(const ex& e);
ex x_evalf(const ex& e);
ex x_deriv(const ex&, unsigned int);

DECLARE_FUNCTION_1P(x);

} // namespace GiNaC

namespace Parma_Recurrence_Relation_Solver {

inline
Expr::Expr() {
}

inline
Expr::Expr(int i)
  : Base(i) {
}

inline
Expr::Expr(const Number& n)
  : Base(n.n) {
}

inline
Expr::Expr(const Symbol& s)
  : Base(s.s) {
}

inline
Expr::Expr(const Constant& k)
  : Base(k.c) {
}

inline
Expr::Expr(const std::string& st, const Expr_List& lst)
  : Base(st, lst.l) {
}

// FIXME: temporary
inline
Expr::Expr(const Expr& lh, const Expr& rh)
  :  Base(GiNaC::relational(static_cast<const Base&>(lh),
			    static_cast<const Base&>(rh))) {
}

inline
Expr::Expr(const Expr& y)
  : Base(y) {
};

inline Expr&
Expr::operator=(const Expr& y) {
  Base::operator=(y);
  return *this;
};

inline
Expr::Expr(const Base& ge)
  : Base(ge) {
}

inline
Expr::Expr(const GiNaC::function& gf)
  : Base(gf) {
}

inline
Expr::~Expr() {
}

inline std::ostream&
operator<<(std::ostream& s, const Expr& e) {
  s << static_cast<const Expr::Base&>(e);
  return s;
};

inline Expr
operator+(const Expr& x) {
  return x;
}

inline Expr
operator-(const Expr& x) {
  return -static_cast<const Expr::Base&>(x);
}

inline Expr
operator+(const Expr& x, const Expr& y) {
  return static_cast<const Expr::Base&>(x) + static_cast<const Expr::Base&>(y);
}

inline Expr
operator-(const Expr& x, const Expr& y) {
  return static_cast<const Expr::Base&>(x) - static_cast<const Expr::Base&>(y);
}

inline Expr
operator*(const Expr& x, const Expr& y) {
  return static_cast<const Expr::Base&>(x) * static_cast<const Expr::Base&>(y);
}

inline Expr
operator/(const Expr& x, const Expr& y) {
  return static_cast<const Expr::Base&>(x) / static_cast<const Expr::Base&>(y);
}

inline Expr&
operator+=(Expr& x, const Expr& y) {
  static_cast<Expr::Base&>(x) += static_cast<const Expr::Base&>(y);
  return x;
}

inline Expr&
operator-=(Expr& x, const Expr& y) {
  static_cast<Expr::Base&>(x) -= static_cast<const Expr::Base&>(y);
  return x;
}

inline Expr&
operator*=(Expr& x, const Expr& y) {
  static_cast<Expr::Base&>(x) *= static_cast<const Expr::Base&>(y);
  return x;
}

inline Expr&
operator/=(Expr& x, const Expr& y) {
  static_cast<Expr::Base&>(x) /= static_cast<const Expr::Base&>(y);
  return x;
}

#if 0
inline bool
operator==(const Expr& x, const Expr& y) {
  return static_cast<const Expr::Base&>(x) == static_cast<const Expr::Base&>(y);
}

inline bool
operator!=(const Expr& x, const Expr& y) {
  return static_cast<const Expr::Base&>(x) != static_cast<const Expr::Base&>(y);
}
#endif

inline Expr
Expr::operator[](int i) const {
  return Base::operator[](i);
}

inline bool
Expr::is_a_symbol() const {
  return GiNaC::is_a<GiNaC::symbol>(*this);
}

inline bool
Expr::is_a_number() const {
  return GiNaC::is_a<GiNaC::numeric>(*this);
}

inline bool
Expr::is_a_constant() const {
  return GiNaC::is_a<GiNaC::constant>(*this);
}

inline bool
Expr::is_a_add() const {
  return GiNaC::is_a<GiNaC::add>(*this);
}

inline bool
Expr::is_a_mul() const {
  return GiNaC::is_a<GiNaC::mul>(*this);
}

inline bool
Expr::is_a_power() const {
  return GiNaC::is_a<GiNaC::power>(*this);
}

inline bool
Expr::is_a_function() const {
  return GiNaC::is_a<GiNaC::function>(*this);
}

inline bool
Expr::is_a_matrix() const {
  return GiNaC::is_a<GiNaC::matrix>(*this);
}

inline bool
Expr::is_a_Expr_List() const {
  return GiNaC::is_a<GiNaC::lst>(*this);
}

inline bool
Expr::is_a_relational() const {
  return GiNaC::is_a<GiNaC::relational>(*this);
}

inline Number
Expr::ex_to_number() const {
  assert(GiNaC::is_a<GiNaC::numeric>(*this));
  return GiNaC::ex_to<GiNaC::numeric>(*this);
}

// FIXME: info, temporary
inline bool
Expr::is_integer_polynomial() const {
  return info(GiNaC::info_flags::integer_polynomial);
}
// FIXME: info, temporary
inline bool
Expr::is_rational_polynomial() const {
  return info(GiNaC::info_flags::rational_polynomial);
}
// FIXME: info, temporary
inline bool
Expr::is_relation_equal() const {
  return info(GiNaC::info_flags::relation_equal);
}

inline unsigned
Expr::nops() const {
  return Base::nops();
}

inline Expr
Expr::op(unsigned i) const {
  return Base::op(i);
}

inline bool
Expr::is_equal(const Expr& exp) const {
  return Base::is_equal(exp);
}

inline bool
Expr::is_zero() const {
  return Base::is_zero();
}

inline Expr
Expr::subs(const Expr& exp1, const Expr& exp2) const {
  return Base::subs(exp1 == exp2);
}

inline Expr
Expr::subs(const Expr_List& to_replace,
	   const Expr_List& replacements) const {
  return Base::subs(to_replace.l, replacements.l);
}

inline bool
Expr::match(const Expr& pattern) const {
  return Base::match(pattern);
}

inline bool
Expr::match(const Expr& pattern, Expr_List& replacements) const {
  return Base::match(pattern, replacements.l);
}

inline bool
Expr::has(const Expr& pattern) const {
  return Base::has(pattern);
}

inline Expr
Expr::expand() const {
  return Base::expand();
}

inline Expr
Expr::collect(const Expr_List& lst) const {
  return Base::collect(lst.l);
}

inline int
Expr::degree(const Symbol& symb) const {
  return Base::degree(symb.s);
}

inline int
Expr::ldegree(const Symbol& symb) const {
  return Base::ldegree(symb.s);
}

inline Expr
Expr::coeff(const Symbol& symb, int k) const {
  return Base::coeff(symb.s, k);
}

inline Expr
Expr::lcoeff(const Symbol& symb) const {
  return Base::lcoeff(symb.s);
}

inline Expr
Expr::tcoeff(const Symbol& symb) const {
  return Base::tcoeff(symb.s);
}

inline Expr
Expr::primpart(const Symbol& symb) const {
  return Base::primpart(symb.s);
}

inline Expr
Expr::numerator() const {
  return Base::numer();
}

inline Expr
Expr::denominator() const {
  return Base::denom();
}

inline void
Expr::numerator_denominator(Expr& numer, Expr& denom) const {
  numer =  Base::numer_denom().op(0);
  denom =  Base::numer_denom().op(1);
}

inline Expr
Expr::lhs() const {
  return Base::lhs();
}

inline Expr
Expr::rhs() const {
  return Base::rhs();
}

inline Expr
Expr::to_rational(Expr_List& lst) {
  return Base::to_rational(lst.l);
}

inline Expr
Expr::diff(const Symbol& symb, unsigned nth) {
  return Base::diff(symb.s, nth);
}

inline Expr
wild(unsigned label) {
  return wild(label);
}

inline Expr
power(const Expr& x, const Expr& y) {
  return GiNaC::pow(static_cast<const Expr::Base>(x),
		    static_cast<const Expr::Base>(y));
}

inline Expr
sqrt(const Expr& x) {
  return GiNaC::sqrt(static_cast<const Expr::Base>(x));
};

inline Expr
sin(const Expr& x) {
  return GiNaC::sin(static_cast<const Expr::Base>(x));
};

inline Expr
cos(const Expr& x) {
  return GiNaC::cos(static_cast<const Expr::Base>(x));
};

inline Expr
acos(const Expr& x) {
  return GiNaC::acos(static_cast<const Expr::Base>(x));
};

inline Expr
tan(const Expr& x) {
  return GiNaC::tan(static_cast<const Expr::Base>(x));
};

inline Expr
exp(const Expr& x) {
  return GiNaC::exp(static_cast<const Expr::Base>(x));
};

inline Expr
log(const Expr& x) {
  return GiNaC::log(static_cast<const Expr::Base>(x));
};

inline Expr
quo(const Expr& a, const Expr& b, const Symbol& x) {
  return GiNaC::quo(a, b, x.s);
}

inline Expr
rem(const Expr& a, const Expr& b, const Symbol& x) {
  return GiNaC::rem(a, b, x.s);
}

inline Expr
prem(const Expr& a, const Expr& b, const Symbol& x) {
  return GiNaC::prem(a, b, x.s);
}

inline Expr
gcd(const Expr& a, const Expr& b) {
  return GiNaC::gcd(a, b);
}

inline Expr
lcm(const Expr& a, const Expr& b) {
  return GiNaC::lcm(a, b);
}

inline Expr
sqrfree(const Expr& x, const Expr_List& lst) {
  return GiNaC::sqrfree(x, lst.l);
}

inline Expr
lsolve(const Expr_List& lst1, const Expr_List& lst2) {
  return GiNaC::lsolve(lst1.l, lst2.l);
}

inline Expr
x(const Expr& y) {
  return x(static_cast<const Expr::Base>(y));
};

inline bool
Expr::is_the_abs_function() const {
  return is_ex_the_function(*this, abs);
}

inline bool
Expr::is_the_exp_function() const {
  return is_ex_the_function(*this, exp);
}

inline bool
Expr::is_the_log_function() const {
  return is_ex_the_function(*this, log);
}

inline bool
Expr::is_the_sin_function() const {
  return is_ex_the_function(*this, sin);
}

inline bool
Expr::is_the_cos_function() const {
  return is_ex_the_function(*this, cos);
}

inline bool
Expr::is_the_tan_function() const {
  return is_ex_the_function(*this, tan);
}

inline bool
Expr::is_the_acos_function() const {
  return is_ex_the_function(*this, acos);
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Expr_inlines_hh)
