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
  Computes the lcm among the integers in the vector \p numbers.
*/
GNumber
lcm(const std::vector<GNumber>& numbers) {
  for (unsigned i = numbers.size() - 1; i-- > 0; )
    assert(numbers[i].is_integer());

  GNumber n = numbers[numbers.size() - 1];
  for (unsigned i = numbers.size() - 1; i-- > 0; )
    n = lcm(n, numbers[i]);
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
  Return <CODE>true</CODE> if the <CODE>GiNaC::GExpr</CODE> \p p is a
  polynomial in a variable 'var', <CODE>false</CODE> otherwise.
*/
static bool
find_polynomial_part(const GExpr& p, const GSymbol& var, bool poly) {
  assert(poly);
  // 'is_a<add>' is necessary even if the 'add' are considerated in
  // function 'assign_poly_part_and_no_poly_part' because could be,
  // for example, (sqrt(6)+2)^a
  if (is_a<mul>(p) || is_a<add>(p))
    for (unsigned i = p.nops(); i-- > 0; )
      poly = find_polynomial_part(p.op(i), var, poly);
  else if (is_a<function>(p)) {
    // Checks the argument of the function and if it contains 'var' then
    // the function is not a polynomial..
    if (p.op(0).has(var))
      poly = false; 
  }
  // NOTE: I suppose that there are not nested powers.
  else if (is_a<power>(p)) {
    poly = find_polynomial_part(p.op(0), var, poly);
    poly = find_polynomial_part(p.op(1), var, poly);
    // If the exponent contains 'var', then 'p' is not a polynomial in 'var'.
    if (p.op(1).has(var))
      poly = false;
    // The exponent is a GiNaC::numeric.
    // The only case in which the power is not a polynomial in 'var' is
    // when the exponent is not a positive integer and the base contains
    // 'var'.
    if (is_a<numeric>(p.op(1))) {
      GNumber exp = GiNaC::ex_to<GiNaC::numeric>(p.op(1));
      if (!exp.is_pos_integer() && p.op(0).has(var))
	poly = false;
    }
  }
  return poly;
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
  1. every <CODE>GiNaC::numeric</CODE> is a
     <CODE>scalar_for_poly_in_var</CODE>;
  2. every <CODE>GiNaC::symbol</CODE> different from \p var is a
     <CODE>scalar_for_poly_in_var</CODE>;
  3. give a function <CODE>f</CODE>, <CODE>f(scalar_for_poly_in_var)</CODE>
     is a <CODE>scalar_for_poly_in_var</CODE>.
     <CODE>f<CODE> is one of those listed in 5.10 of GiNaC 1.0.8.
  Step 2
  Definition of <CODE>polynomial in var</CODE> in inductive way
  1. every <CODE>scalar_for_poly_in_var</CODE> is a
     <CODE>polynomial_in_var</CODE>;
  2. \p var is a <CODE>polynomial_in_var</CODE>;
  3. give the binary operations sum (\f$ + \f$) and multiplication (\f$ * \f$),
     <CODE>scalar_for_poly_in_var + scalar_for_poly_in_var</CODE>
     is a </CODE>polynomial_in_var</CODE> and
     <CODE>scalar_for_poly_in_var * scalar_for_poly_in_var</CODE>
     is a <CODE>polynomial_in_var</CODE>;
  4. give i of the type <CODE>GiNaC::numeric</CODE> with
     <CODE>i.info(info_flags::posint) == true</CODE>
     <CODE>power(polynomial_in_var, i)</CODE> is a
     <CODE>polynomial_in_var</CODE>.
*/
void
assign_poly_part_and_no_poly_part(const GExpr& p, const GSymbol& var,
				  GExpr& p_poly, GExpr& p_no_poly) {
  if (is_a<add>(p)) {
    p_poly = 0;
    p_no_poly = 0;
    for (unsigned i = p.nops(); i-- > 0; ) {
      bool poly = true;
      poly = find_polynomial_part(p.op(i), var, poly);
      if (poly)
	p_poly += p.op(i);
      else
	p_no_poly += p.op(i);
    }
  }
  else {
    p_poly = 1;
    p_no_poly = 1;
    bool poly = true;
    poly = find_polynomial_part(p, var, poly);
    if (poly) {
      p_poly *= p;
      p_no_poly = 0;
    }
    else {
      p_no_poly *= p;
      p_poly = 0;
    }
  }
}
