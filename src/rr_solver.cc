/* Definition of the main recurrence relation solver.
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

#include <config.h>

#include "globals.hh"
#include "util.hh"
#include "alg_eq_solver.hh"
#include <climits>

// TEMPORARY
#include <iostream>

using namespace GiNaC;

static GExpr
get_binding(const GList& l, unsigned wild_index) {
  assert(wild_index < l.nops());
  assert(l.op(wild_index).info(GiNaC::info_flags::relation_equal));
  assert(l.op(wild_index).lhs() == GiNaC::wild(wild_index));
  return l.op(wild_index).rhs();
}

static bool
get_linear_decrement(const GExpr& e, const GSymbol& n, GNumber& decrement) {
  static GExpr n_plus_d = n + GiNaC::wild(0);
  GList substitution;
  if (match(e, n_plus_d, substitution)) {
    GExpr d = get_binding(substitution, 0);
    if (GiNaC::is_a<GiNaC::numeric>(d)) {
      decrement = -GiNaC::ex_to<GiNaC::numeric>(d);
      return true;
    }
  }
  return false;
}

/*!
  Transforms experssions of the form \f$a^(b*x)\f$ in the expression \p p
  into \f$(a^b)^x\f$.
*/
GExpr
split_exp(GExpr& p, const GSymbol& x, const GList powers) {
  GExpr q = 0;

  for (size_t i = powers.nops(); i-- > 0; ) {
    // Finds coefficient of exponential 'powers.op(i)'.
    GList l;
    GExpr coeff;
    if (p.find(powers.op(i) * wild(), l)) {
      coeff = l.op(0);
      coeff = coeff.subs(wild(0)*pow(wild(1), wild(2)*x) == wild(0));
    }
    GExpr tmp = powers.op(i);
    tmp = tmp.subs(pow(wild(0),x*wild(1)) == pow(pow(wild(0),wild(1)),x));
    q += tmp*coeff;
    p -= coeff*powers.op(i);
  }
  // Now p does not contain other exponential of the form a^(b*x). 
  q += p;

  return q;
}

/*!
  Transforms expressions of the form \f$a^x*b^x\f$ in the expression \p p
  into \f$(a^b)^x\f$.
*/
GExpr 
union_exp(GExpr& p, const GSymbol& x, const GList powers, const unsigned n) {
  GExpr q = 0;

  for (size_t i = powers.nops(); i-- > 0; ) {
    GExpr tmp = powers.op(i);
    if (n == 1) {
      tmp = tmp.subs(pow(wild(0),x)*pow(wild(1),x) == pow(wild(0)*wild(1),x));
      q += tmp;
    } 
    else {
      tmp = tmp.subs(pow(wild(0),x)*pow(wild(1),x)*wild(2) == 
		     pow(wild(0)*wild(1),x));
      // Finds coefficient of exponential 'powers.op(i)'.  
      GList l;
      GExpr coeff;
      if (p.find(powers.op(i) * wild(), l)) {
	coeff = l.op(0);
	coeff = coeff.subs(wild(0)*pow(wild(1), x)*pow(wild(2), x) == wild(0));
      }
      q += tmp*coeff;
    }
    p -= powers.op(i);
  }
  // Now p does not contain other exponential of the form a^x*b^x.
  q += p;

  return q;
}

void
check_exp_inhomogeneous_term(GExpr& e, const GSymbol& n) {
  GList powers;

  // Transforms a^(bn+c) into (a^b)^n*a^c.
  static GExpr a_bn = pow(wild(0), n*wild(1));
  static GExpr a_n_c = pow(wild(0), n+wild(1));
  static GExpr a_bn_c = pow(wild(0), n*wild(1)+wild(2));
  if (e.find(a_bn_c, powers) || e.find(a_n_c, powers)) 
    e = e.expand();
  if (e.find(a_bn, powers))
    e = split_exp(e, n, powers);

  // Transforms a^n*b^n into (a*b)^n.
  static GExpr a_n_b_n = pow(wild(0), n)*pow(wild(1), n);
  static GExpr c_a_n_b_n = pow(wild(0), n)*pow(wild(1), n)*wild(2);
  while (e.find(a_n_b_n, powers) || e.find(c_a_n_b_n, powers)) { 
    if (e.find(a_n_b_n, powers))
	e = union_exp(e, n, powers, 1);
    else
      e = union_exp(e, n, powers, 2);
    clear(powers);
  }
}

/*!
  This function makes a matrix with two rows and a number of columns
  does not exceed the number of exponentials in the inhomogeneous term
  plus one.
  The function gives the decomposition 
  \f$ e(n) = sum_{i=0}^k \alpha_i^j \cdot p(n)_i \f$ with
  - \f$ \alpha_i \ne \alpha_j \f$ if \f$ i \ne j \f$
  - \p p does not contains exponentials.
  It returns the matrix whose \f$ i\f$-th column contains 
  \f$ \alpha_i^n \f$ and \f$ p(n)_i \f$
  for \f$ i = 1, \ldots, k \f$.
*/
GMatrix
decomposition_inhomogeneous_term(const GExpr& e, const GSymbol& n) {
  GExpr p,q;
  GList(powers);

  p = e;
  GExpr pattern = pow(wild(), n);
  p.find(pattern, powers);
  p = p.collect(powers);

  if (powers.nops() == 0) {
    // In this case there are no axponentials.
    return GMatrix(2, 1, lst(1, p));
  }
  else {
    GMatrix terms_divided(2, powers.nops()+1, powers);
    // Impose that the last element of the first row is the constant
    // exponential 1.
    terms_divided(0, powers.nops()) = 1;
    GList part;
    q = p;
    for (size_t i = powers.nops(); i-- > 0; ) {
      clear(part);
      // 'part' has always only one elements because we have collected
      // 'p' with respect to the exponentials.
      p.find(powers.op(i) * wild(), part);
      if (part.nops() != 0) {
	GExpr coeff = part.op(0);
	coeff = coeff.subs(wild(1)*pow(wild(2), n) == wild(1));
	terms_divided(1, i) = coeff;
	q -= part.op(0);
      }
      else {
	// In this case 'p' does not contain the constant exponential.
	terms_divided(1, i) = 1;
	q -= powers.op(i);
      }
    }
    // Now 'q' does not contains any exponentials or product of exponential 
    // times other expressions.
    if (!q.is_zero())
      terms_divided(1, powers.nops()) = q;
  
    return terms_divided;
  }
}

bool
solve(const GExpr& rhs, const GSymbol& n) {
  static GExpr x_i = x(GiNaC::wild(0));
  static GExpr x_i_plus_r = x_i + GiNaC::wild(1);
  static GExpr a_times_x_i = GiNaC::wild(1)*x_i;
  static GExpr a_times_x_i_plus_r = a_times_x_i + GiNaC::wild(2);

  static GList substitution;

  int degree = -1;
  std::vector<GNumber> coefficients;
  GExpr e = rhs;
  bool failed = false;
  do {
    GExpr i;
    GExpr a;
    // The following matches are attempted starting from the most common,
    // then the second most common and so forth.
    if (clear(substitution), match(e, x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = get_binding(substitution, 1);
    }
    else if (clear(substitution), match(e, a_times_x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = get_binding(substitution, 2);
    }
    else if (clear(substitution), match(e, a_times_x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = 0;
    }
    else if (clear(substitution), match(e, x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = 0;
    }
    else
      break;

    GNumber decrement;
    if (!get_linear_decrement(i, n, decrement)) {
      failed = true;
      break;
    }
    std::cout << "decrement = " << decrement << std::endl;
    if (!decrement.is_integer()
	|| decrement < 0
	|| decrement >= coefficients.max_size()) {
      failed = true;
      break;
    }
    if (!GiNaC::is_a<GiNaC::numeric>(a)) {
      failed = true;
      break;
    }
    GNumber coefficient = GiNaC::ex_to<GiNaC::numeric>(a);
    // FIXME: turn this assertion into something more appropriate.
    assert(decrement >= LONG_MIN && decrement <= LONG_MAX);
    unsigned long index = decrement.to_long();
    if (degree < 0 || index > unsigned(degree))
      degree = index;
    if (index > coefficients.size())
      coefficients.insert(coefficients.end(),
			  index - coefficients.size(),
			  GNumber(0));
    if (index == coefficients.size())
      coefficients.insert(coefficients.end(), coefficient);
    else
      coefficients[index] += coefficient;
  } while (e != 0);

  // See if what is left is the inhomogeneous term,
  // i.e., if all the occurrences of `x(e)' are such that
  // `e' is numeric.
  GList occurrences;
  if (e.find(x_i, occurrences)) {
    for (unsigned i = 0, n = occurrences.nops(); i < n; ++i)
      if (!is_a<numeric>(occurrences.op(i))) {
	failed = true;
	break;
      }
  }

  if (failed)
    return false;

  std::cout << "Degree = " << degree << std::endl;
  std::cout << "Coefficients = ";
  for (int i = 1; i <= degree; ++i)
    std::cout << coefficients[i] << " ";
  std::cout << std::endl;
  std::cout << "Inhomogeneous term = " << e << std::endl;

  // The factors of the form a^(bx+c) (a,b,c numeric) must be transformed
  // into (a^b)^x*a^c. GiNaC tranforms only a^(bx+c) in a^c*a^(bx) but not
  // a^(bx) into (a^b)^x.
  check_exp_inhomogeneous_term(e, n);

  // Now certainly the inhomogeneous term contains only exponential
  // with exponent x.
  // 'decomposition' is a matrix with two rows and a number of columns
  // which is at most the number of exponentials in the inhomogeneous term
  // plus one.
  // In every column there is a exponential in the first row and its
  // coefficient in the second row. In the last column there is the
  // constant exponential with its coefficients. 
  GMatrix decomposition = decomposition_inhomogeneous_term(e, n);
  std::cout << "Inhomogeneous term's decomposition"
	    << decomposition << std::endl; 

  /*
  // Build the expression here.
  static GSymbol x("x");
  GExpr p = 0;
  for (int i = 1; i <= degree; ++i)
    p += coefficients[i]*pow(x, i);
  
  std::vector<Polynomial_Root> roots;
  bool all_distinct;
  if (!find_roots(p, x, roots, all_distinct))
    return false;

  std::cout << degree << " " << roots.size() << " " << all_distinct << std::endl;
  */
  return true;
}
