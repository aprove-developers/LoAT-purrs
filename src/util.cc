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
  We assume that \p substitution has been produced by GiNaC::match()
  and that the binding for the wildcard of index \p wild_index
  is in the position \p wild_index of \p substitution.
*/
GExpr
get_binding(const GList& substitution, unsigned wild_index) {
  assert(wild_index < substitution.nops());
  assert(substitution.op(wild_index).info(GiNaC::info_flags::relation_equal));
  assert(substitution.op(wild_index).lhs() == GiNaC::wild(wild_index));
  return substitution.op(wild_index).rhs();
}

/*!
  Return <CODE>true</CODE> if the <CODE>GiNaC::GExpr</CODE> \p p is an
  object <CODE>scalar_for_poly_in_var</CODE>, <CODE>false</CODE> otherwise.
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
  Return <CODE>true</CODE> if the <CODE>GiNaC::GExpr</CODE> \p p is an
  object <CODE>polynomial_in_var</CODE>, <CODE>false</CODE> otherwise.
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
  Give a <CODE>GiNaC::GExpr</CODE> \p p and a <CODE>GiNaC::GSymbol</CODE>
  \p var, builds two other <CODE>GiNaC::GExpr</CODE> \p p_poly and \p p_no_poly
  that contain the polynomial part and the non-polynomial part of \p p
  regarding the variable \p var.
  The polynomial part of an expression is the sum of those terms that are
  polynomials in a variable in according to the following definition
  in two steps (when the polynomial part lacks \p p_poly is zero and so
  also for non polynomial part).
  Step 1
  We consider the variable \p var.
  Definition the object <CODE>scalar_for_poly_in_var</CODE> in inductive way:
  - every <CODE>GiNaC::numeric</CODE> is a
    <CODE>scalar_for_poly_in_var</CODE>;
  - every <CODE>GiNaC::constant</CODE> is a
    <CODE>scalar_for_poly_in_var</CODE>;
  - every <CODE>GiNaC::symbol</CODE> different from \p var is a
    <CODE>scalar_for_poly_in_var</CODE>;
  - given \f$ e \f$ a <CODE>GiNaC::power</CODE>, if \f$ e.op(0) \f$ and
    \f$ e.op(1) \f$ are <CODE>scalar_for_poly_in_var</CODE>
    then \f$ e \f$ is a <CODE>scalar_for_poly_in_var</CODE>;
  - given \f$ f \f$ a <CODE>GiNaC::function</CODE>,
    if \f$ f.op(0) \f$ is <CODE>scalar_for_poly_in_var</CODE>
    then \f$ f \f$ is a <CODE>scalar_for_poly_in_var</CODE>;
  - given the binary operations sum (\f$ + \f$) and multiplication (\f$ * \f$),
    if \f$ a \f$ and \f$ b \f$ are <CODE>scalar_for_poly_in_var</CODE> then
    \f$ a + b \f$ and \f$ a * b \f$ are <CODE>scalar_for_poly_in_var</CODE>.
  Step 2
  Definition of <CODE>polynomial_in_var</CODE> in inductive way
  - every <CODE>scalar_for_poly_in_var</CODE> is a
    <CODE>polynomial_in_var</CODE>;
  - \p var is a <CODE>polynomial_in_var</CODE>;
  - given \f$ e \f$ a <CODE>GiNaC::power</CODE>,
    if \f$ e.op(0) \f$ is <CODE>polynomial_in_var</CODE> and
    \f$ e.op(1) \f$ is <CODE>GiNaC::numeric</CODE> such that
    <CODE>e.op(1).GiNaC::is_pos_integer()</CODE>
    then \f$ e \f$ is a <CODE>valid_base</CODE>;
  - given the binary operations sum (\f$ + \f$) and multiplication (\f$ * \f$),
    if \f$ a \f$ and \f$ b \f$ are <CODE>polynomial_in_var</CODE> then
    \f$ a + b \f$ and \f$ a * b \f$ are <CODE>polynomial_in_var</CODE>.
*/
void
assign_poly_part_and_no_poly_part(const GExpr& p, const GSymbol& var,
				  GExpr& p_poly, GExpr& p_no_poly) {
  if (is_a<add>(p)) {
    p_poly = 0;
    p_no_poly = 0;
    for (unsigned i = p.nops(); i-- > 0; ) {
      if (is_polynomial(p.op(i), var))
	p_poly += p.op(i);
      else
	p_no_poly += p.op(i);
    }
  }
  else {
    p_poly = 1;
    p_no_poly = 1;
    if (is_polynomial(p, var)) {
      p_poly *= p;
      p_no_poly = 0;
    }
    else {
      p_no_poly *= p;
      p_poly = 0;
    }
  }
}
