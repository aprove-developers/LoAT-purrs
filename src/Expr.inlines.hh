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

DECLARE_FUNCTION_1P(x);
DECLARE_FUNCTION_4P(sum);
DECLARE_FUNCTION_4P(prod);

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
Expr::Expr(unsigned int i)
  : Base(i) {
}

inline
Expr::Expr(long i)
  : Base(i) {
}

inline
Expr::Expr(unsigned long i)
  : Base(i) {
}

inline
Expr::Expr(const Number& y)
  : Base(y.n) {
}

inline
Expr::Expr(const Symbol& y)
  : Base(y.s) {
}

inline
Expr::Expr(const Constant& y)
  : Base(y.c) {
}

inline
Expr::Expr(const std::string& s, const Expr_List& y)
  : Base(s, y.l) {
}

// FIXME: temporary
inline
Expr::Expr(const Expr& x, const Expr& y)
  :  Base(GiNaC::relational(static_cast<const Base&>(x),
			    static_cast<const Base&>(y))) {
}

inline
Expr::Expr(const Expr& y)
  : Base(y) {
}

inline Expr&
Expr::operator=(const Expr& y) {
  Base::operator=(y);
  return *this;
}

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

inline bool
Expr::is_a_symbol() const {
  return GiNaC::is_a<GiNaC::symbol>(*this);
}

inline bool
Expr::is_a_number() const {
  return GiNaC::is_a<GiNaC::numeric>(*this);
}

inline bool
Expr::is_a_number(Number& n) const {
  if (GiNaC::is_a<GiNaC::numeric>(*this)) {
    n = GiNaC::ex_to<GiNaC::numeric>(*this);
    return true;
  }
  return false;
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

inline Symbol
Expr::ex_to_symbol() const {
  assert(GiNaC::is_a<GiNaC::symbol>(*this));
  return GiNaC::ex_to<GiNaC::symbol>(*this);
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

inline std::ostream&
operator<<(std::ostream& s, const Expr& x) {
  s << static_cast<const Expr::Base&>(x);
  return s;
}

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

inline bool
operator==(const Expr& e, const Symbol& s) {
  return e.is_a_symbol() && GiNaC::ex_to<GiNaC::symbol>(e).is_equal(s.s);
}

inline bool
operator!=(const Expr& e, const Symbol& s) {
  return !(e == s);
}

inline bool
operator==(const Symbol& s, const Expr& e) {
  return e == s;
}

inline bool
operator!=(const Symbol& s, const Expr& e) {
  return !(s == e);
}

inline bool
operator==(const Expr& e, const Constant& c) {
  return e.is_a_constant() && GiNaC::ex_to<GiNaC::constant>(e).is_equal(c.c);
}

inline bool
operator!=(const Expr& e, const Constant& c) {
  return !(e == c);
}

inline bool
operator==(const Constant& c, const Expr& e) {
  return e == c;
}

inline bool
operator!=(const Constant& c, const Expr& e) {
  return !(c == e);
}

inline bool
operator==(const Expr& e, const Number& n) {
  return e.is_a_number() && GiNaC::ex_to<GiNaC::numeric>(e).is_equal(n.n);
}

inline bool
operator==(const Number& n, const Expr& e) {
  return e == n;
}

inline bool
operator==(const Expr& e, long i) {
  return e.is_a_number()
    && GiNaC::ex_to<GiNaC::numeric>(e).operator==(i);
}

inline bool
operator==(const Expr& x, const Expr& y) {
  return x.Base::is_equal(y);
}

inline bool
operator!=(const Expr& x, const Expr& y) {
  return !(x == y);
}

inline Expr
Expr::operator[](int i) const {
  return Base::operator[](i);
}

inline Functor
Expr::functor() const {
  assert(is_a_function());
  return GiNaC::ex_to<GiNaC::function>(*this).get_serial();
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
Expr::is_zero() const {
  return Base::is_zero();
}

inline Expr
Expr::subs(const Expr& x, const Expr& y) const {
  return Base::subs(GiNaC::operator==(x, y));
}

inline Expr
Expr::subs(const Expr_List& x,
	   const Expr_List& y) const {
  return Base::subs(x.l, y.l);
}

inline bool
Expr::match(const Expr& x) const {
  return Base::match(x);
}

inline bool
Expr::match(const Expr& x, Expr_List& y) const {
  return Base::match(x, y.l);
}

inline bool
Expr::has(const Expr& x) const {
  return Base::has(x);
}

inline Expr
Expr::expand() const {
  return Base::expand();
}

inline Expr
Expr::collect(const Expr_List& x) const {
  return Base::collect(x.l);
}

inline unsigned
Expr::degree(const Symbol& x) const {
  int d = Base::degree(x.s);
  assert(d >= 0);
  return unsigned(d);
}

inline unsigned
Expr::ldegree(const Symbol& x) const {
  int d = Base::ldegree(x.s);
  assert(d >= 0);
  return unsigned(d);
}

inline Expr
Expr::coeff(const Symbol& x, int k) const {
  return Base::coeff(x.s, k);
}

inline Expr
Expr::lcoeff(const Symbol& x) const {
  return Base::lcoeff(x.s);
}

inline Expr
Expr::tcoeff(const Symbol& x) const {
  return Base::tcoeff(x.s);
}

inline Expr
Expr::primpart(const Symbol& x) const {
  return Base::primpart(x.s);
}

inline Expr
Expr::content(const Symbol& x) const {
  return Base::content(x.s);
}

inline void
Expr::numerator_denominator(Expr& x, Expr& y) const {
  x =  Base::numer_denom().op(0);
  y =  Base::numer_denom().op(1);
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
Expr::diff(const Symbol& x, unsigned nth) {
  return Base::diff(x.s, nth);
}

inline Expr
wild(unsigned label) {
  return GiNaC::wild(label);
}

inline Expr
apply(Functor f, const Expr& x) {
  return GiNaC::function(f, static_cast<const Expr::Base>(x));
}

inline Expr
pwr(const Expr& x, const Expr& y) {
  return GiNaC::pow(static_cast<const Expr::Base>(x),
		    static_cast<const Expr::Base>(y));
}

inline Expr
sqrt(const Expr& x) {
  return GiNaC::sqrt(static_cast<const Expr::Base>(x));
}

inline Expr
sin(const Expr& x) {
  return GiNaC::sin(static_cast<const Expr::Base>(x));
}

inline Expr
cos(const Expr& x) {
  return GiNaC::cos(static_cast<const Expr::Base>(x));
}

inline Expr
acos(const Expr& x) {
  return GiNaC::acos(static_cast<const Expr::Base>(x));
}

inline Expr
tan(const Expr& x) {
  return GiNaC::tan(static_cast<const Expr::Base>(x));
}

inline Expr
exp(const Expr& x) {
  return GiNaC::exp(static_cast<const Expr::Base>(x));
}

inline Expr
log(const Expr& x) {
  return GiNaC::log(static_cast<const Expr::Base>(x));
}

inline Expr
factorial(const Expr& x) {
  return GiNaC::factorial(static_cast<const Expr::Base>(x));
}

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
sqrfree(const Expr& x, const Expr_List& y) {
  return GiNaC::sqrfree(x, y.l);
}

inline Expr
lsolve(const Expr_List& x, const Expr_List& y) {
  return GiNaC::lsolve(x.l, y.l);
}

inline Expr
x(const Expr& y) {
  return x(static_cast<const Expr::Base>(y));
}

inline Expr
sum(const Expr& index, const Expr& lower, const Expr& upper,
    const Expr& summand) {
  return sum(static_cast<const Expr::Base>(index),
	     static_cast<const Expr::Base>(lower),
	     static_cast<const Expr::Base>(upper),
	     static_cast<const Expr::Base>(summand));
}

inline Expr
prod(const Expr& index, const Expr& lower, const Expr& upper,
     const Expr& factor) {
  return prod(static_cast<const Expr::Base>(index),
	      static_cast<const Expr::Base>(lower),
	      static_cast<const Expr::Base>(upper),
	      static_cast<const Expr::Base>(factor));
}

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

inline bool
Expr::is_the_x_function() const {
  return is_ex_the_function(*this, x);
}

inline bool
Expr::is_the_sum_function() const {
  return is_ex_the_function(*this, sum);
}

inline bool
Expr::is_the_prod_function() const {
  return is_ex_the_function(*this, prod);
}

inline void
Expr::latex_print(std::ostream& s) {
  print(GiNaC::print_latex(s));
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Expr_inlines_hh)
