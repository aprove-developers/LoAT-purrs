/* To be written
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

#ifndef NOISY
#define NOISY 0
#endif

#include "zeilberger.hh"
#include "globals.hh"
#include "simplify.hh"
#include "numerator_denominator.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"

#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

/*!
  By definition, an expression \f$ F(m, k) \f$ is an hypergeometric term
  in both argument, i.e., \f$ F(m+1, k) / F(m, k) \f$ and
  \f$ K(m, k+1) / F(m, k) \f$ are both rational functions of \f$ m \f$
  and \f$ k \f$.
  This function returns <CODE>true</CODE> if \p F_m_k is a hypergeometric
  term and in this case, called
  \f[
    t(k) = a_0 F(m, k) + a_1 F(m+1, k) + \dots a_j F(m+j, k),
  \f]
  computes \f$ p_0(k) \f$, \f$ r(k) \f$ and \f$ s(k) \f$ so that
  \f[
    \frac{t(k+1)}{t(k)} = \frac{p(k+1)}{p_0(k)} \frac{r(k)}{s(k)}.
  \f]
  and store them in \p p_0_k, \p r_k and \p s_k.
  Returns <CODE>false</CODE> if \p F_m_k is not a hypergeometric term.
*/
bool
zeilberger_step_one(const Expr& F_m_k,
		    const Symbol& m, const Symbol& k,
		    const std::vector<Symbol>& coefficients,
		    Expr& p_0_k, Expr& r_k, Expr& s_k) {
  DD_MSGVAR("ZEIL step one -> ", F_m_k);
  Expr first
    = simplify_binomials_factorials_exponentials(F_m_k.substitute(m, m+1))
    * pwr(simplify_binomials_factorials_exponentials(F_m_k), -1);
  Expr second
    = simplify_binomials_factorials_exponentials(F_m_k.substitute(k, k+1))
    * pwr(simplify_binomials_factorials_exponentials(F_m_k), -1);
  // FIXME: we must understand the better simplification to use in this case.
  first = simplify_numer_denom(first);
  second = simplify_numer_denom(second);
  if (!first.is_rational_function(m) || !second.is_rational_function(k))
    return false;

  // F_m_k is an hypergeometric term in both argument and then the
  // Zeilberger's algorithm is applicable.
  Expr tmp
    = simplify_binomials_factorials_exponentials(F_m_k)
    * pwr(simplify_binomials_factorials_exponentials(F_m_k.substitute(m, m-1)),
	  -1);
  tmp = simplify_numer_denom(tmp);

  Expr r_1;
  Expr r_2;
  second.numerator_denominator(r_1, r_2);
  if (r_1.is_a_mul())
    for (unsigned int i = r_1.nops(); i-- > 0; )
      if (r_1.op(i) == -1) {
	r_1 *= -1;
	r_2 *= -1;
	break;
      }
  D_VAR(r_1);
  D_VAR(r_2);
  Expr s_1;
  Expr s_2;
  tmp.numerator_denominator(s_1, s_2);
  if (s_1.is_a_mul())
    for (unsigned int i = s_1.nops(); i-- > 0; )
      if (s_1.op(i) == -1) {
	s_1 *= -1;
	s_2 *= -1;
	break;
      }
  D_VAR(s_1);
  D_VAR(s_2);

  // Computes `p_0_k'.
  index_type order = coefficients.size() - 1;
  for (index_type j = 0; j <= order; ++j) {
    Expr prod_s_1 = 1;
    for (index_type i = 0; i < j; ++i)
      prod_s_1 *= s_1.substitute(m, m+j-i);
    Expr prod_s_2 = 1;
    for (index_type i = j+1; i <= order; ++i)
      prod_s_2 *= s_2.substitute(m, m+i);
    p_0_k += coefficients[j] * prod_s_1 * prod_s_2;
  }

  // Computes `r_k' and `s_k'.
  Expr prod_for_r = 1;
  Expr prod_for_s = 1;
  for (index_type i = 1; i <= order; ++i) {
    const Expr& tmp = s_2.substitute(m, m+i);
    prod_for_r *= tmp;
    prod_for_s *= tmp.substitute(k, k+1);
  }
  r_k = r_1 * prod_for_r;
  s_k = r_2 * prod_for_s;
  
  return true;
}

bool
parametric_gosper_step_two(const Symbol& /*m*/, const Expr& /*r_m*/,
			   Expr& /*a_m*/, Expr& /*b_m*/, Expr& /*c_m*/) {

  return true;
}


bool
parametric_gosper_step_three(const Symbol& /*m*/,
			     const Expr& /*a_m*/, const Expr& /*b_m*/,
			     const Expr& /*c_m*/, Expr& /*x_m*/) {

  return true;
}


} // anonymous namespace

/*!
  \param F_m_k    proper hypergeometric sequence.
  \param m        first variable of \p F_m_k.
  \param k        second variable of \p F_m_k.

  \return         ...
*/
bool
PURRS::zeilberger_algorithm(const Expr& F_m_k,
			    const Symbol& m, const Symbol& k) {
  Expr p_0_k = 0;
  Expr r_k = 0;
  Expr s_k = 0;
  // FIXME: temporary.
  // We must consider the maximum order for the
  // recurrence ... and, starting from the lower, if the algorithm fails,
  // to increase the order until `order' and to repeat the algorithm.
  index_type order = 1; // TEMPORARY
  std::vector<Symbol> coefficients(order + 1);
  for (index_type i = 0; i < order+1; ++i)
    coefficients[i] = Symbol();
  if (!zeilberger_step_one(F_m_k, m, k, coefficients, p_0_k, r_k, s_k))
    return false;
  DD_VAR(p_0_k);
  DD_VAR(r_k);
  DD_VAR(s_k);

  Expr p_1_k;
  Expr p_2_k;
  Expr p_3_k;
#if 0
  if (!parametric_gosper_step_two(k, r_k / s_k, p_1_k, p_2_k, p_3_k))
    // Problem in the computation of the resultant and its roots.
    return false;
#else
  const Expr& tmp = simplify_all(r_k / s_k);
  p_1_k = 1;
  p_2_k = numerator(tmp);
  p_3_k = denominator(tmp);
#endif
  DD_VAR(p_1_k);
  DD_VAR(p_2_k);
  DD_VAR(p_3_k);

  Expr b_k;
  if (parametric_gosper_step_three(k, p_2_k.expand(), p_3_k.expand(),
				   (p_0_k * p_1_k).expand(), b_k)) {
  }

  return true;
}
