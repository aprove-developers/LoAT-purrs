/* Definition of some utility functions.
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

#include "util.hh"

using namespace GiNaC;

/*!
  Computes the gcd between the integers \p n and \p m.
*/
int
gcd(int n, int m) {
  int r = m;
  while (r != 0){
    r = n % m;
    n = m; 
    m = r;
  }
  return n;  
}

/*!
  Computes the lcm among the numbers in the vector \p v.
*/
GNumber
lcm(const std::vector<GNumber>& v) {
  GNumber n = 1;
  for (unsigned i = v.size(); i-- > 0; )
    n = lcm(n, v[i]);
  return n;
}

GExpr
cubic_root(const GExpr& e) {
  static GExpr one_third = GExpr(1)/3;
  return pow(e, one_third);
}

void
clear(GList& l) {
  for (unsigned n = l.nops(); n-- > 0; )
    l.remove_first();
  assert(l.nops() == 0);
}

/*!
  We assume that \p substitution has been produced by
  <CODE>GiNaC::match()</CODE> and that the binding for the wildcard of
  index \p wild_index is in the position \p wild_index of \p substitution.
*/
GExpr
get_binding(const GList& substitution, unsigned wild_index) {
  assert(wild_index < substitution.nops());
  assert(substitution.op(wild_index).info(GiNaC::info_flags::relation_equal));
  assert(substitution.op(wild_index).lhs() == GiNaC::wild(wild_index));
  return substitution.op(wild_index).rhs();
}

/*!
  Return <CODE>true</CODE> if the <CODE>GiNaC::GExpr</CODE> \p p is a
  scalar for poly in var, <CODE>false</CODE> otherwise.
*/
static bool
is_scalar_for_poly(const GExpr& e, const GSymbol& var) {
  bool found;
  if (is_a<numeric>(e))
    found = true;
  else if (is_a<constant>(e))
    found = true;
  else if (is_a<symbol>(e) && !e.is_equal(var))
    found = true;
  else if (is_a<power>(e))
    found = is_scalar_for_poly(e.op(0), var)
      && is_scalar_for_poly(e.op(1), var);
  else if (is_a<function>(e))
    found = is_scalar_for_poly(e.op(0), var);
  else if (is_a<add>(e) || is_a<mul>(e)) {
    found = true;
    for (unsigned i = e.nops(); i-- > 0; )
      found = found && is_scalar_for_poly(e.op(i), var);
  }
  else
    found = false;
  
  return found;
}


/*!
  Return <CODE>true</CODE> if the <CODE>GiNaC::GExpr</CODE> \p p is a
  polynomial in var, <CODE>false</CODE> otherwise.
*/
static bool
is_polynomial(const GExpr& e, const GSymbol& var) {
  bool found;
  if (is_scalar_for_poly(e, var))
    found = true;
  else if (e.is_equal(var))
    found = true;
  else if (is_a<power>(e)) {
    bool exp_ok = false;
    if (is_a<numeric>(e.op(1))) {
      GNumber exp = GiNaC::ex_to<GiNaC::numeric>(e.op(1));
      if (exp.is_pos_integer())
	exp_ok = true;
    }
    if (is_polynomial(e.op(0), var) && exp_ok)
      found = true;
    else
      found = false;
  }
  else if (is_a<add>(e) || is_a<mul>(e)) {
    found = true;
    for (unsigned i = e.nops(); i-- > 0; )
      found = found && is_polynomial(e.op(i), var);
  }
  else
    found = false;

  return found;
}

/*!
  This function realized the definition of <EM>polynomial_in_var</EM>.
  Given an expression \p p and a symbol \p var, this function builds two
  other expression \p poly and \p no_poly that contain the polynomial part
  and the non-polynomial part of \p p regarding the variable \p var.
  The polynomial part of an expression is the sum of those terms that are
  polynomials in a variable in according to the following definition
  in two steps (when the polynomial part lacks \p poly is zero and so
  also for non polynomial part).
  Step 1
  We consider the variable \p var.
  Definition of the object <EM>scalar for poly in var</EM> inductively as
  follows:
  - every number is a scalar for poly in var;
  - every symbolic constant is a scalar for poly in var;
  - every parameter different from \f$ var \f$ is a scalar for poly in var;
  - if \f$ f \f$ is any unary function and \f$ x \f$ is scalar for poly in var,
    then \f$ f(x) \f$ is a scalar for poly in var;
  - if \f$ a \f$ and \f$ b \f$ are scalars for poly in var then
    \f$ a + b \f$, \f$ a * b \f$ and \f$ a^b \f$ are scalars for poly in var.
  Step 2
  Definition of <EM>polynomial_in_var</EM> inductively as follows:
  - every scalar for poly in var is a polynomial in var;
  - \f$ var \f$ is a polynomial in var;
  - if \f$ a \f$ is a polynomial in var and \f$ b \f$ is a positive integer,
    then \f$ a^b \f$ is a polynomial in var;
  - if \f$ a \f$ and \f$ b \f$ are polynomials in var then
    \f$ a + b \f$ and \f$ a * b \f$ are polynomials in var.
*/
void
assign_polynomial_part(const GExpr& p, const GSymbol& var,
		       GExpr& poly, GExpr& no_poly) {
  if (is_a<add>(p)) {
    poly = 0;
    no_poly = 0;
    for (unsigned i = p.nops(); i-- > 0; ) {
      if (is_polynomial(p.op(i), var))
	poly += p.op(i);
      else
	no_poly += p.op(i);
    }
  }
  else {
    poly = 1;
    no_poly = 1;
    if (is_polynomial(p, var)) {
      poly *= p;
      no_poly = 0;
    }
    else {
      no_poly *= p;
      poly = 0;
    }
  }
}
