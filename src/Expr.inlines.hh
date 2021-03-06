/* Expr class implementation: inline functions.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
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
#include <cassert>

namespace GiNaC {

DECLARE_FUNCTION_1P(floor);
DECLARE_FUNCTION_2P(Sc);
DECLARE_FUNCTION_2P(mod);
DECLARE_FUNCTION_2P(binom);
DECLARE_FUNCTION_4P(sum);
DECLARE_FUNCTION_4P(prod);
DECLARE_FUNCTION_2P(max);
DECLARE_FUNCTION_2P(min);

// Use overloading to allow for multiple arguments in x().
class x1_SERIAL { public: static unsigned serial; };
template<typename T1>
inline function x(const T1& p1) {
  return function(x1_SERIAL::serial, ex(p1));
}

class x2_SERIAL { public: static unsigned serial; };
  template<typename T1, typename T2>
  inline function x(const T1& p1, const T2& p2) {
    return function(x2_SERIAL::serial, ex(p1), ex(p2));
  }

template<> inline bool is_the_function<class x_SERIAL>(const ex& e)
  {
    return is_the_function<x1_SERIAL>(e) || is_the_function<x2_SERIAL>(e);
  }
                                                                                                                             



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
Expr::is_a_symbol(Symbol& s) const {
  if (GiNaC::is_a<GiNaC::symbol>(*this)) {
    s = GiNaC::ex_to<GiNaC::symbol>(*this);
    return true;
  }
  return false;
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

inline Functor
Expr::functor() const {
  assert(is_a_function());
  return GiNaC::ex_to<GiNaC::function>(*this).get_serial();
}

inline std::string
Expr::get_function_name() const {
  assert(is_a_function());
  return GiNaC::ex_to<GiNaC::function>(*this).get_name();
}

inline unsigned int
Expr::nops() const {
  return Base::nops();
}

inline Expr
Expr::op(unsigned int i) const {
  assert(!is_a_function() && !is_a_power());
  // If `i' is out of range (0, nops()-1) in GiNaC 1.1.5 returns `0'.
  if (i > nops()-1)
    throw std::out_of_range("PURRS::Expr::op(): the index of `op()' "
			    "must be between 0 and `nops()-1'");
  return Base::op(i);
}

#ifdef UNSAFE_ARG
inline Expr&
Expr::arg(unsigned int i) {
  assert(is_a_function() || is_a_power());
  // FIXME: if `i' is out of range (0, nops()-1) in GiNaC 1.1.5 happens a
  // segmentation fault. 
  if (i > nops()-1)
    throw std::out_of_range("PURRS::Expr::arg(): the index of `arg()' "
			    "must be between 0 and `nops()-1'");
  return static_cast<Expr&>(let_op(i));
}

inline const Expr&
Expr::arg(unsigned int i) const {
  assert(is_a_function() || is_a_power());
  // FIXME: if `i' is out of range (0, nops()-1) in GiNaC 1.1.5 happens a
  // segmentation fault. 
  if (i > nops()-1)
    throw std::out_of_range("PURRS::Expr::arg(): the index of `arg()' "
			    "must be between 0 and `nops()-1'");
  return static_cast<Expr&>(const_cast<Expr&>(*this).let_op(i));
}
#else
inline Expr
Expr::arg(unsigned int i) const {
  assert(is_a_function() || is_a_power());
  // FIXME: if `i' is out of range (0, nops()-1) in GiNaC 1.1.5 happens a
  // segmentation fault. 
  if (i > nops()-1)
    throw std::out_of_range("PURRS::Expr::arg(): the index of `arg()' "
			    "must be between 0 and `nops()-1'");
  return Base::op(i);
}
#endif

inline bool
Expr::is_zero() const {
  return Base::is_zero();
}

inline bool
Expr::has(const Expr& x) const {
  return Base::has(x);
}

inline Expr
Expr::expand() const {
  return Base::expand();
}

inline bool
Expr::is_expanded() const {
  return *this == expand();
}

inline Expr
Expr::collect(const Expr_List& x) const {
  return Base::collect(x.l);
}

inline unsigned int
Expr::degree(const Symbol& x) const {
  assert(is_expanded());
  assert(is_polynomial(x.s));
  int d = Base::degree(x.s);
  assert(d >= 0);
  return unsigned(d);
}

inline unsigned int
Expr::ldegree(const Symbol& x) const {
  assert(is_expanded());
  assert(is_polynomial(x.s));
  int d = Base::ldegree(x.s);
  assert(d >= 0);
  return unsigned(d);
}

inline Expr
Expr::coeff(const Symbol& x, int k) const {
  assert(is_expanded());
  assert(is_polynomial(x.s));
  return Base::coeff(x.s, k);
}

inline Expr
Expr::lcoeff(const Symbol& x) const {
  assert(is_expanded());
  assert(is_polynomial(x.s));
  return Base::lcoeff(x.s);
}

inline Expr
Expr::tcoeff(const Symbol& x) const {
  assert(is_expanded());
  assert(is_polynomial(x.s));
  return Base::tcoeff(x.s);
}

inline Expr
Expr::primpart(const Symbol& x) const {
  assert(is_expanded());
  assert(is_polynomial(x.s));
  return Base::primpart(x.s);
}

inline Expr
Expr::content(const Symbol& x) const {
  assert(is_expanded());
  assert(is_polynomial(x.s));
  return Base::content(x.s);
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
Expr::numerator_denominator(Expr& x, Expr& y) const {
  const Expr& tmp = Base::numer_denom();
  x = tmp.op(0);
  y = tmp.op(1);
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
Expr::diff(const Symbol& x, unsigned int nth) {
  return Base::diff(x.s, nth);
}

inline Expr
Expr::unsafe_fp_approximation() const {
  return Base::evalf();
}

inline Functor
find_functor(const std::string& name, unsigned int num_args) {
  return GiNaC::function::find_function(name, num_args);
}

inline Expr
apply(Functor f, const Expr& x) {
  return GiNaC::function(f, static_cast<const Expr::Base>(x));
}

inline Expr
apply(Functor f, const Expr& x1, const Expr& x2) {
  return GiNaC::function(f,
			 static_cast<const Expr::Base>(x1),
			 static_cast<const Expr::Base>(x2));
}

inline Expr
apply(Functor f, const Expr& x1, const Expr& x2, const Expr& x3) {
  return GiNaC::function(f,
			 static_cast<const Expr::Base>(x1),
			 static_cast<const Expr::Base>(x2),
			 static_cast<const Expr::Base>(x3));
}

inline Expr
apply(Functor f,
      const Expr& x1, const Expr& x2, const Expr& x3, const Expr& x4) {
  return GiNaC::function(f,
			 static_cast<const Expr::Base>(x1),
			 static_cast<const Expr::Base>(x2),
			 static_cast<const Expr::Base>(x3),
			 static_cast<const Expr::Base>(x4));
}

inline Expr
apply(Functor f, const std::vector<Expr>& x) {
  unsigned int x_size = x.size();
  assert(x_size > 1);
  switch (x_size) {
  case 2:
    return apply(f, x[0], x[1]);
    break;
  case 3:
    return apply(f, x[0], x[1], x[2]);
    break;
  case 4:
    return apply(f, x[0], x[1], x[2], x[3]);
    break;
  default:
    {
      std::vector<GiNaC::ex> tmp_x(x_size);
      for (unsigned int i = 0; i < x_size; ++i)
	tmp_x[i] = static_cast<const Expr::Base>(x[i]);
      return GiNaC::function(f, tmp_x);
    }
  }
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
binom(const Expr& n, const Expr& k) {
  return binom(static_cast<const Expr::Base>(n),
	       static_cast<const Expr::Base>(k));
}

inline Expr
gamma(const Expr& x) {
  return GiNaC::tgamma(static_cast<const Expr::Base>(x));
}

inline Expr
quo(const Expr& a, const Expr& b, const Symbol& x) {
  assert(a.is_polynomial(x.s));
  assert(b.is_polynomial(x.s));
  return GiNaC::quo(a, b, x.s);
}

inline Expr
rem(const Expr& a, const Expr& b, const Symbol& x) {
  assert(a.is_polynomial(x.s));
  assert(b.is_polynomial(x.s));
  return GiNaC::rem(a, b, x.s);
}

inline Expr
prem(const Expr& a, const Expr& b, const Symbol& x) {
  assert(a.is_polynomial(x.s));
  assert(b.is_polynomial(x.s));
  return GiNaC::prem(a, b, x.s);
}

inline Expr
gcd(const Expr& a, const Expr& b) {
  assert(a.is_multivariate_polynomial());
  assert(b.is_multivariate_polynomial());
  return GiNaC::gcd(a, b);
}

inline Expr
lcm(const Expr& a, const Expr& b) {
  assert(a.is_multivariate_polynomial());
  assert(b.is_multivariate_polynomial());
  return GiNaC::lcm(a, b);
}

inline Expr
lsolve(const Expr_List& x, const Expr_List& y) {
  return GiNaC::lsolve(x.l, y.l);
}

inline Expr
x(const Expr& y) {
  Number num;
  if (y.is_a_number(num))
    assert(num.is_nonnegative_integer());
  return x(static_cast<const Expr::Base>(y));
}

inline Expr
floor(const Expr& x) {
  return floor(static_cast<const Expr::Base>(x));
}

inline Expr
Sc(const Expr& x, const Expr& y) {
  return Sc(static_cast<const Expr::Base>(x),
	    static_cast<const Expr::Base>(y));
}

inline Expr
mod(const Expr& x, const Expr& y) {
  return mod(static_cast<const Expr::Base>(x),
	     static_cast<const Expr::Base>(y));
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

inline Expr
max(const Expr& x, const Expr& y) {
  return max(static_cast<const Expr::Base>(x),
	     static_cast<const Expr::Base>(y));
}

inline Expr
min(const Expr& x, const Expr& y) {
  return min(static_cast<const Expr::Base>(x),
	     static_cast<const Expr::Base>(y));
}

inline bool
Expr::is_the_abs_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, abs);
}

inline bool
Expr::is_the_exp_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, exp);
}

inline bool
Expr::is_the_log_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, log);
}

inline bool
Expr::is_the_sin_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, sin);
}

inline bool
Expr::is_the_cos_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, cos);
}

inline bool
Expr::is_the_tan_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, tan);
}

inline bool
Expr::is_the_acos_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, acos);
}

inline bool
Expr::is_the_factorial_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, factorial);
}

inline bool
Expr::is_the_gamma_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, tgamma);
}

inline bool
Expr::is_the_binom_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, binom);
}

inline bool
Expr::is_the_x_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, x);
}

// FIXME: Use underlying GiNaC methods if possible.
inline bool
Expr::is_the_x1_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, x) && ((*this).nops() == 1);
}

// FIXME: Use underlying GiNaC methods if possible.
inline bool
Expr::is_the_x2_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, x) && ((*this).nops() == 2);
}

inline bool
Expr::is_the_floor_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, floor);
}

inline bool
Expr::is_the_Sc_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, Sc);
}

inline bool
Expr::is_the_mod_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, mod);
}

inline bool
Expr::is_the_sum_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, sum);
}

inline bool
Expr::is_the_prod_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, prod);
}

inline bool
Expr::is_the_max_function() const {
  using namespace GiNaC;
  return is_ex_the_function(*this, max);
}

inline void
Expr::latex_print(std::ostream& s) {
  print(GiNaC::print_latex(s));
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Expr_inlines_hh)
