/* Recurrence class implementation (non-inline functions).
   Copyright (C) 2001, 2002 Roberto Bagnara <bagnara@cs.unipr.it>

This file is part of the Parma Polyhedra Library (PPL).

The PPL is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The PPL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.

For the most up-to-date information see the Parma Polyhedra Library
site: http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#include "Expr.defs.hh"

#include "Symbol.types.hh"
#include "Expr_List.types.hh"

#include <ginac/ginac.h>

namespace PURRS = Parma_Recurrence_Relation_Solver;

bool
PURRS::Expr::is_a_symbol(const Expr& exp) {
  return GiNaC::is_a<GiNaC::symbol>(exp.e);
}

bool
PURRS::Expr::is_a_number(const Expr& exp) {
  return GiNaC::is_a<GiNaC::numeric>(exp.e);
}

bool
PURRS::Expr::is_a_constant(const Expr& exp) {
  return GiNaC::is_a<GiNaC::constant>(exp.e);
}

bool
PURRS::Expr::is_a_add(const Expr& exp) {
  return GiNaC::is_a<GiNaC::add>(exp.e);
}

bool
PURRS::Expr::is_a_mul(const Expr& exp) {
  return GiNaC::is_a<GiNaC::mul>(exp.e);
}

bool
PURRS::Expr::is_a_power(const Expr& exp) {
  return GiNaC::is_a<GiNaC::power>(exp.e);
}

bool
PURRS::Expr::is_a_function(const Expr& exp) {
  return GiNaC::is_a<GiNaC::function>(exp.e);
}

bool
PURRS::Expr::is_a_matrix(const Expr& exp) {
  return GiNaC::is_a<GiNaC::matrix>(exp.e);
}

bool
PURRS::Expr::is_a_Expr_List(const Expr& exp) {
  return GiNaC::is_a<GiNaC::lst>(exp.e);
}

bool
PURRS::Expr::is_a_relational(const Expr& exp) {
  return GiNaC::is_a<GiNaC::relational>(exp.e);
}

//bool
//PURRS::Expr::info(unsigned flag) {
//}

unsigned
PURRS::Expr::nops() const {
  return e.nops();
}

PURRS::Expr
PURRS::Expr::op(unsigned i) const {
  Expr tmp;
  tmp.e = e.op(i);
  return tmp;
}

bool
PURRS::Expr::is_equal(const Expr& exp) const {
  return e.is_equal(exp.e);
}

bool
PURRS::Expr::is_zero() const {
  return e.is_zero();
}

PURRS::Expr
PURRS::Expr::subs(const Expr& exp) const {
  
}
#if 0
PURRS::Expr
PURRS::Expr::subs(const Expr_List& symbols,
		  const Expr_List& replacements) const {
}

bool
PURRS::Expr::match(const Expr& pattern) const {
}

bool
PURRS::Expr::match(const Expr& pattern, Expr_List& replacements) const {
}

bool
PURRS::Expr::has(const Expr& pattern) const {
}

// lo uso?
//bool
//PURRS::Expr::find(const Expr& pattern, Expr& found) const {
//}

PURRS::Expr
PURRS::Expr::expand() {
}

//
PURRS::Expr
PURRS::Expr::collect(const Expr& lst, bool distributed) {
}

int
PURRS::Expr::degree(const Expr& exp) const {
}

int
PURRS::Expr::ldegree(const Expr& exp) const {
}

PURRS::Expr
PURRS::Expr::coeff(const Expr& exp, int k) const {
}

PURRS::Expr
PURRS::Expr::lcoeff(const Expr& exp) const {
}

PURRS::Expr
PURRS::Expr::tcoeff(const Expr& exp) const {
}

PURRS::Expr
PURRS::Expr::quo(const Expr& a, const Expr& b, const Symbol& x) const {
}

PURRS::Expr
PURRS::Expr::rem(const Expr& a, const Expr& b, const Symbol& x) const {
}

PURRS::Expr
PURRS::Expr::prem(const Expr& a, const Expr& b, const Symbol& x) const {
}

bool
PURRS::Expr::divide(const Expr& a, const Expr& b, const Expr& q) const {
}

PURRS::Expr
PURRS::Expr::primpart(const Symbol& x) {
}

PURRS::Expr
PURRS::Expr::gcd(const Expr& a, const Expr& b) const {
}

PURRS::Expr
PURRS::Expr::lcm(const Expr& a, const Expr& b) const {
}

PURRS::Expr
PURRS::Expr::sqrfree(const Expr& exp, const Expr_List& lst) const {
}

PURRS::Expr
PURRS::Expr::numer() const {
}

PURRS::Expr
PURRS::Expr::denom() const {
}

PURRS::Expr
PURRS::Expr::numer_denom() const {
}

PURRS::Expr
PURRS::Expr::to_rational(Expr_List& lst) {
}

//
int
PURRS::Expr::to_int() const {
}

long
PURRS::Expr::to_long() const {
}
#endif
// solve(), lsolve()
