/* To be written.
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

#include "factorize.hh"
#include "numerator_denominator.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

void
assign_common_factor_and_rem(const Expr& base_factor,
			     const Expr& exponent_factor,
			     const Expr& exp_common_factor,
			     Expr& common_factor, Expr& rem_summand) {
  std::vector<Expr> bases;
  std::vector<Expr> exponents;
  // Find bases and exponents of each factor of `common_factor'.
  split_bases_exponents(common_factor, bases, exponents);
  bool new_common_factor = true;
  for (unsigned i = bases.size(); i-- > 0; ) {
    if (base_factor == bases[i]) {
      new_common_factor = false;
      const Expr& exponents_i = exponents[i];
      Number num_exp_f;
      Number num_i;
      if (exponent_factor.is_a_number(num_exp_f)
	  && exponents_i.is_a_number(num_i)
	  && num_exp_f.is_positive_integer()
	  && num_i.is_positive_integer())
	if (num_exp_f > num_i)
	  rem_summand *= pwr(base_factor, num_exp_f - num_i);
	else
	  rem_summand *= pwr(base_factor, num_i - num_exp_f);
      else if (exponent_factor != exponents_i)
	rem_summand *= pwr(base_factor, exponent_factor - exponents_i);
      break;
    }
  }
  if (new_common_factor) {
    common_factor *= pwr(base_factor, exp_common_factor);
    if (exp_common_factor != exponent_factor)
      rem_summand *= pwr(base_factor, exponent_factor - exp_common_factor);
  }
}

/*!
  Let \f$ e = e_1 + e_2 + \dots + e_k \f$ be the expression contained
  in \p e and let \p factor be a factor of one subexpression of \p e,
  i. e., \f$ \exists i \in \{ 1, \cdots, k \} \itc e_i = factor * r \f$,
  where \f$ r \f$ represents the possibly other factors of \f$ e_i \f$.
  This function looks if \p factor is a common factor to every term
  of \p e and, in this case, if it is not already in \p common_factor,
  add it. If \p factor is not a common_factor to every term then
  it is stored in \p rem_summand.
*/
void
in_all_factors(const Expr& e,
	       const Expr& base_factor, const Expr& exponent_factor,
	       Expr& common_factor, Expr& rem_summand) {
  D_MSG("***I");
  D_VAR(base_factor);
  D_VAR(exponent_factor);
  Expr exp_common_factor = exponent_factor;
  bool is_in_every_term = true;
  // Run over all terms of `e' in order to look if `factor' is in
  // every term.
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& e_i = e.op(i);
    D_VAR(e_i);
    std::vector<Expr> bases;
    std::vector<Expr> exponents;
    // Find bases and exponents of each factor of `e_i'.
    split_bases_exponents(e_i, bases, exponents);
    D_VEC(bases, 0, bases.size()-1);
    D_VEC(exponents, 0, exponents.size()-1);
    bool found = false;
    for (unsigned j = bases.size(); j-- > 0; ) {
      if (bases[j] == base_factor) {
	Number num_exp_common_f;
	Number num_exp_j;
	if (exp_common_factor.is_a_number(num_exp_common_f)
	    && exponents[j].is_a_number(num_exp_j)
	    && num_exp_common_f.is_positive_integer()
	    && num_exp_j.is_positive_integer()) {
	  exp_common_factor = num_exp_j > num_exp_common_f
	    ? num_exp_common_f : num_exp_j;
	  found = true;
	  break;
	}
	else if (exp_common_factor == exponents[j]) {
	  found = true;
	  break;
	}
      }
    }
    D_VAR(exp_common_factor);
    if (!found) {
      is_in_every_term = false;
      break;
    }
  }
  if (is_in_every_term)
    assign_common_factor_and_rem(base_factor, exponent_factor,
				 exp_common_factor,
				 common_factor, rem_summand);
  else
    rem_summand *= pwr(base_factor, exponent_factor);
  D_VAR(common_factor);
  D_VAR(rem_summand);
  D_MSG("***F");
}

/*!
  Let \f$ e = e_1 + e_2 + \dots + e_k \f$ be the expression contained
  in \p e and let \p term be one subexpression of \p e, i. e.,
  \f$ \exists i \in \{ 1, \cdots, k \} \itc e_i = term \f$.
  This function returns in \p common_factor the factors of \p term
  that are contained also in all other terms of \p e, i. e. the common
  factors, and leaves in \p rem_summand the other factors of \p term.
*/
void
collect_common_factor(const Expr& e, const Expr& term,
		      Expr& common_factor, Expr& rem_summand) {
  std::vector<Expr> bases;
  std::vector<Expr> exponents;
  // Find bases and exponents of each factor of `term'.
  split_bases_exponents(term, bases, exponents);
  D_MSGVAR("COLLECT I ", term);
//    D_VEC(bases, 0, bases.size()-1);
//    D_VEC(exponents, 0, exponents.size()-1);
  for (unsigned i = bases.size(); i-- > 0; )
    in_all_factors(e, bases[i], exponents[i], common_factor, rem_summand);
  D_MSGVAR("COLLECT F ", common_factor);
  D_VAR(rem_summand);
} 

} // anonymous namespace

/*!
  This function factorizes as much possible the expression \p e,
  where \p e is not a ratios of expression.
  Returns the two values \p common_factor and \p remainder
  such that the product of them is equal to \p e.
*/
void
PURRS::factorize_no_ratio_ex(const Expr& e,
			     Expr& common_factor, Expr& remainder) {
  assert(denominator(e) == 1);
  D_MSG("");
  D_MSGVAR("INPUT ", e);
  Expr e_factorized = e.expand();
  if (e.is_rational_polynomial()) {
    e_factorized = sqrfree(e_factorized);
    D_MSG("sqrfree");
  }
  D_VAR(e_factorized);

  common_factor = 1;
  remainder = 0;
  if (e_factorized.is_a_add()) {
    for (unsigned i = e_factorized.nops(); i-- > 0; ) {
      Expr rem_summand = 1;
      collect_common_factor(e_factorized, e_factorized.op(i),
			    common_factor, rem_summand);
      remainder += rem_summand;
    }
  }
  else if (e_factorized.is_a_power()) {
    Expr base;
    Number exponent;
    if (e_factorized.arg(1).is_a_number(exponent)
	&& exponent.is_positive_integer())
      base = e_factorized.arg(0);
    else {
      base = e_factorized;
      exponent = 1;
    }
    Expr rem_summand = 1;
    if (base.is_a_add()) {
      factorize_no_ratio_ex(base, common_factor, remainder);
      if (exponent != 1) {
	common_factor = pwr(common_factor, exponent);
	remainder = pwr(remainder, exponent);
      }
    }
    else
      remainder += e_factorized;
  }
  else if (e_factorized.is_a_mul()) {
    Expr rem_all_factors = 1;
    for (unsigned i = e_factorized.nops(); i-- > 0; ) {
      const Expr& factor = e_factorized.op(i);
      Expr common = 1;
      Expr rem = 0;
      factorize_no_ratio_ex(factor, common, rem);
      common_factor *= common;
      rem_all_factors *= rem; 
    }
    remainder = rem_all_factors;
  }
  else
    remainder = e;

  D_MSGVAR("OUTPUT ", common_factor);
  D_VAR(remainder);
  D_MSG("");
}

/*!
  This function factorizes as much possible the expression \p e
  and returns the two values \p common_factor and \p remainder
  such that the product of them is equal to \p e.
*/
void
PURRS::factorize(const Expr& e,
		 Expr& common_factor, Expr& remainder) {
  D_MSGVAR("factorize ", e);
  Expr numerator;
  Expr denominator;
  numerator_denominator_purrs(e, numerator, denominator);
  Expr common_num;
  Expr rem_num;
  factorize_no_ratio_ex(numerator, common_num, rem_num);
  Expr common_den;
  Expr rem_den;
  factorize_no_ratio_ex(denominator, common_den, rem_den);
  if (common_den.is_a_number() && rem_den.is_a_number()) {
    common_factor = common_num / (common_den * rem_den);
    remainder = rem_num;
  }
  else {
    common_factor = common_num / common_den;
    remainder = rem_num / rem_den;
  }
}
