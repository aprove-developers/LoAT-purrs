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

#include "numerator_denominator.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include <vector>

// TEMPORARY
#include <iostream>

namespace Parma_Recurrence_Relation_Solver {

/*!
  Given \f$ e \f$ a factor of the \f$ i \f$-th term of a more general
  expression \f$ f \f$, this function splits its numerator and denominator
  putting them in the \f$ i \f$-th position of the vectors \p numerators or
  \p denominators.
  Returns the denominator of \f$ e \f$ while in the \f$ i \f$-th position
  of the vectors \p numerators and \p denominators are built the numerator
  and the denominator of \f$ f \f$.
*/
static Expr
find_denominator_single_factor(const Expr& e, unsigned position,
			       std::vector<Expr>& numerators,
			       std::vector<Expr>& denominators) {
  Number exponent;
  // `e' is a denominator.
  if (e.is_a_power() && e.arg(1).is_a_number(exponent)
      && !exponent.is_nonnegative_integer()) {
    const Expr& den = pwr(e.arg(0), -exponent);
    denominators[position] *= den;
    return den;
  }
  // `e' is not a denominator.
  else {
    numerators[position] *= e;
    return 1;
  }
}

/*!
  Given \f$ e \f$ the \f$ i \f$-th term of a more general expression \f$ f \f$,
  this function splits its numerator and denominator
  putting them in the \f$ i \f$-th position of the vectors \p numerators or
  \p denominators.
  Returns the denominator of \f$ e \f$.
*/
static Expr
find_denominator_single_term(const Expr& e, unsigned position,
			     std::vector<Expr>& numerators,
			     std::vector<Expr>& denominators) {
  Expr denominator = 1;
  if (e.is_a_mul())
    for (unsigned i = e.nops(); i-- > 0; )
      denominator *= find_denominator_single_factor(e.op(i), position,
						    numerators, denominators);
  else
    denominator = find_denominator_single_factor(e, position,
						 numerators, denominators);
  return denominator;
}

/*!
  We consider \f$ d \f$ and \f$ f \f$ splitted in its base and exponent.
  Returns <CODE>false</CODE> if \f$ f \f$ is already a factor of \f$ d \f$.
  Returns <CODE>true</CODE> if \f$ f \f$ is not a factor of \f$ d \f$:
  in this case, if in \f$ d \f$ there is a factor with the same base \f$ b \f$
  of \f$ f \f$ and the exponents \f$ e_d \f$ and \f$ e_f \f$ are both positive
  integers with \f$ e_f > e_d \f$, then in \p factor_to_add is stored the
  factor \f$ b^{e_f - e_d} \f$.
*/
static bool
check_factor_is_to_add(const Expr& base_f, const Expr& exponent_f,
		       const Number& numeric_exponent_f,
		       const Expr& d, Expr& factor_to_add) {
  bool add_factor = true;
  Expr base_d;
  Expr exponent_d;
  // `d' is not a `power'.
  if (d.is_a_power()) {
    base_d = d.arg(0);
    exponent_d = d.arg(1);
    // The two factor have same bases.
    if (base_d == base_f) {
      Number numeric_exponent_d;
      // Both exponents are numeric.
      if (numeric_exponent_f != 0
	  && exponent_d.is_a_number(numeric_exponent_d)) {
	// Both exponents are positive integers.
	if (numeric_exponent_f.is_positive_integer()
	    && numeric_exponent_d.is_positive_integer()) {
	  if (numeric_exponent_f != numeric_exponent_d) {
	    if (numeric_exponent_f > numeric_exponent_d)
	      factor_to_add = pwr(base_f,
				  numeric_exponent_f - numeric_exponent_d);
	    else
	      add_factor = false;
	  }
	  else
	    add_factor = false;
	}
      }
      // Both exponents are not numeric.
      else if (numeric_exponent_f == 0
	       && !exponent_d.is_a_number())
	if (exponent_f == exponent_d)
	  add_factor = false;
    }
  }
  // `d' is not a `power'.
  else
    if (d == base_f)
      if (numeric_exponent_f.is_positive_integer()) {
	if (numeric_exponent_f != 1)
	  factor_to_add = pwr(base_f, numeric_exponent_f - 1);
	// Also `f' is not a power and `d = base_f^numeric_exponent_f'.
	else
	  add_factor = false;
      }
  return add_factor;
}

static Expr
take_common_single_factor(const Expr& base_f, const Expr& exponent_f,
			  const Number& numeric_exponent_f,
			  const Expr& d) {
  Expr common_factor = d;
  Expr factor_to_add = 0;
  // `d' is a `mul'.
  if (d.is_a_mul()) {
    bool add_factor = true;
    for (unsigned i = d.nops(); i-- > 0; )
      if (check_factor_is_to_add(base_f, exponent_f, numeric_exponent_f,
				 d.op(i), factor_to_add)) {
	if (!factor_to_add.is_zero()) {
	  common_factor *= factor_to_add;
	  add_factor = false;
	  break;
	}
      }
      else {
	add_factor = false;
	// The factor does not added to `common_factor'.
	break;
      }
    if (add_factor)
      if (numeric_exponent_f == 0)
	common_factor *= pwr(base_f, exponent_f);
      else
	common_factor *= pwr(base_f, numeric_exponent_f);
  }
  // `d' is not a `mul'.
  else
    if (check_factor_is_to_add(base_f, exponent_f, numeric_exponent_f,
			       d, factor_to_add)) {
      if (factor_to_add == 0)
	if (numeric_exponent_f == 0)
	  common_factor *= pwr(base_f, exponent_f);
	else
	  common_factor *= pwr(base_f, numeric_exponent_f);
      else
	common_factor *= factor_to_add;
    }
  return common_factor;
}

/*!
  Individuates base and exponent (differing if it is numeric or not)
  of the expression \p factor.
  If exponent is numeric \p exponent is setted to \f$ 0 \f$, otherwise
  \p numeric_exponent is setted to \f$ 0 \f$.
*/
static void
split_factor(const Expr& factor, Expr& base, Expr& exponent,
	     Number& numeric_exponent) {
  if (factor.is_a_power()) {
    base = factor.arg(0);
    exponent = factor.arg(1);
    if (exponent.is_a_number(numeric_exponent))
      exponent = 0;
    else
      numeric_exponent = 0;
  }
  else {
    base = factor;
    exponent = 0;
    numeric_exponent = 1;
  }
}

static void
take_common_single_factors(const Expr& d, const Expr& f,
			   Expr& common_factors) {
  Expr temp_d = d;
  Expr base_f;
  Expr exponent_f;
  Number numeric_exponent_f;
  split_factor(f, base_f, exponent_f, numeric_exponent_f);
  if (base_f.is_a_mul())
    for (unsigned i = base_f.nops(); i-- > 0; ) {
      common_factors = take_common_single_factor(base_f.op(i), exponent_f,
						  numeric_exponent_f, temp_d);
      temp_d = common_factors;
    }
  else
    common_factors = take_common_single_factor(base_f, exponent_f,
					       numeric_exponent_f, d);
}

/*!
  We consider \f$ d = {d_1}^{e_1} \dots {d_k}^{e_k} \f$ and \f$ f = b^e \f$.
  This function builds a new expression taking all factors, common and
  not common, of the two expressions with the maximum exponent.
  More exactly:
  - if \f$ \exists i \in \{1, \dotsc , k\} \st {d_i}^{e_i} = f \f$,
    then returns \f$ d \f$;
  - if \f$ \exists i \in \{1, \dotsc , k\} \st d_i = b \f$, \f$ e_i \f$ and
    \f$ e \f$ are different integer positive numbers,
    if \f$ e_i > e \f$ then returns \f$ d \f$, otherwise
    returns \f$ {d_1}^{e_1} \dots {d_i}^e \dots {d_k}^{e_k}\f$;
  - returns \f$ d f \f$ in all other cases.
*/
static Expr
take_common_factors(const Expr& d, const Expr& f) {
  Expr common_factors = 1;
  Expr temp_d = d;
  if (f.is_a_mul())
    for (unsigned i = f.nops(); i-- > 0; ) {
      const Expr& factor = f.op(i);
      take_common_single_factors(temp_d, factor, common_factors);
      temp_d = common_factors;
    }
  else
    take_common_single_factors(d, f, common_factors);
  return common_factors;
}

/*!
  We consider \f$ d = {b_d}^{e_d} \f$ and \f$ f = {b_f}^{e_f} \f$.
  This function builds a new expression \f$ r \f$ in the following way:
  - if \f$ b_d != b_f \f$ then \f$ r = d \f$;
  - if \f$ b_d = b_f = b \f$ and the two exponent are positive integers,
    then \f$ r = b^{e_d - e_f} \f$;
  - in all other cases \f$ r = 1 \f$.
*/
static Expr
find_factor_for_numerator(const Expr& d, const Expr& f) {
  Expr base_d;
  Expr exponent_d;
  Number numeric_exponent_d;
  split_factor(d, base_d,
	       exponent_d, numeric_exponent_d);
  Expr rem = 1;
  if (f.is_a_mul()) {
    bool d_factor_of_f = false;
    for (unsigned i = f.nops(); i-- > 0; )
      if (d != f.op(i)) {
	Expr base_f;
	Expr exponent_f;
	Number numeric_exponent_f;
	split_factor(f.op(i), base_f,
		     exponent_f, numeric_exponent_f);
	if (base_f == base_d && numeric_exponent_f.is_positive_integer()
	    && numeric_exponent_d.is_positive_integer()) {
	  rem *= pwr(base_f, numeric_exponent_d - numeric_exponent_f);
	  d_factor_of_f = true;
	}
      }
      else
	d_factor_of_f = true;
    if (!d_factor_of_f)
      rem *= d;
  }
  else
    if (d != f) {
      Expr base_f;
      Expr exponent_f;
      Number numeric_exponent_f;
      split_factor(f, base_f,
		   exponent_f, numeric_exponent_f);
      if (base_f == base_d && numeric_exponent_f.is_positive_integer()
	  && numeric_exponent_d.is_positive_integer())
	rem *= pwr(base_f, numeric_exponent_d - numeric_exponent_f);
      else
	rem *= d;
    }
  return rem;
}


/*!
  Let \f$ e(n) = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.
  The vectors \p numerators and \p denominators contain
  \f$ n_1, \dots, n_k \f$ and \f$ d_1, \dots, d_k \f$, respectively;
  \p denominator contain the common denominator of \f$ d_1, \dots, d_k \f$
  in according with the explanation of the function
  <CODE>take_common_factors()</CODE>.
  This function computes \f$ \sum_{i = 1}^k \frac{d}{d_i} n_i \f$.
*/
static Expr
find_numerator(const std::vector<Expr>& numerators,
	       const std::vector<Expr>& denominators,
	       const Expr& denominator) {
  assert(numerators.size() == denominators.size());
  Expr numerator = 0;
  for (unsigned i = numerators.size(); i-- > 0; ) {
    const Expr& i_th_denominator = denominators[i];
    Expr multiply_to_numerator = 1;
    // If `denominator == denominators[i]' we do not have factors to multiply
    // to numerator.
    if (denominator != i_th_denominator)
      if (denominator.is_a_mul()) {
	for (unsigned j = denominator.nops(); j-- > 0; )
	  multiply_to_numerator *= find_factor_for_numerator(denominator.op(j),
							     i_th_denominator);
      }
      else
	multiply_to_numerator *= find_factor_for_numerator(denominator,
							   i_th_denominator);
    numerator += numerators[i] * multiply_to_numerator;
  }
  return numerator;
}

static void
numerator_denominator_factor(const Expr& e,
			     Expr& numerator, Expr& denominator) {
  // The dimension of `numerators' and `denominators' is terms'number of `e'
  // and will contain the numerator and the denominator of each term of `e'.  
  std::vector<Expr> numerators;
  std::vector<Expr> denominators;
  numerator = 1;
  denominator = 1;
  if (e.is_a_add()) {
    numerators.insert(numerators.begin(), e.nops(), Number(1));
    denominators.insert(denominators.begin(), e.nops(), Number(1));
    for (unsigned i = e.nops(); i-- > 0; ) {
      Expr factor_denominator = find_denominator_single_term(e.op(i), i,
 							     numerators,
							     denominators);
      Expr old_denominator = denominator;
      denominator = take_common_factors(old_denominator, factor_denominator);
    }
    numerator = find_numerator(numerators, denominators, denominator);
  }
  else if (e.is_a_power()) {
    numerators.push_back(1);
    denominators.push_back(1);
    Expr base_num;
    Expr base_den;
    numerator_denominator_factor(e.arg(0), base_num, base_den);
    Expr exp_num;
    Expr exp_den;
    numerator_denominator_factor(e.arg(1), exp_num, exp_den);
    denominator
      = find_denominator_single_term(pwr(base_num/base_den, exp_num/exp_den),
				     0, numerators, denominators);
    numerator = numerators[0];
  }
  else if (e.is_a_mul())
    for (unsigned i = e.nops(); i-- > 0; ) {
      if (e.op(i).is_a_power()) {
	numerators.push_back(1);
	denominators.push_back(1);
	Expr num;
	Expr den;
	numerator_denominator_factor(e.op(i), num, den);
	denominator *= find_denominator_single_term(num/den, 0,
						    numerators, denominators);
	numerator = numerators[0];
      }
      else {
	numerators.push_back(1);
	denominators.push_back(1);
	denominator *= find_denominator_single_term(e.op(i), 0,
						    numerators, denominators);
	numerator = numerators[0];
      }
    }
  else {
    numerators.push_back(1);
    denominators.push_back(1);
    denominator = find_denominator_single_term(e, 0, numerators, denominators);
    numerator = numerators[0];
  }
}

void
numerator_denominator_purrs(const Expr& e, Expr& numerator, Expr& denominator) {
  // The dimension of `numerators' and `denominators' is terms'number of `e'
  // and will contain the numerator and the denominator of each term of `e'.  
  std::vector<Expr> numerators;
  std::vector<Expr> denominators;
  numerator = 1;
  denominator = 1;
  if (e.is_a_add()) {
    numerators.insert(numerators.begin(), e.nops(), Number(1));
    denominators.insert(denominators.begin(), e.nops(), Number(1));
    for (unsigned i = e.nops(); i-- > 0; ) {
      Expr tmp_denominator;
      numerator_denominator_factor(e.op(i), numerator, tmp_denominator);
      numerators[i] = numerator;
      denominators[i] = tmp_denominator;
      Expr& old_denominator = denominator;
      denominator = take_common_factors(old_denominator, tmp_denominator);
    }
    numerator = find_numerator(numerators, denominators, denominator);
  }
  else
    numerator_denominator_factor(e, numerator, denominator);
}

/*!
  Let \f$ e(n) = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.
  This function returns \f$ n = \sum_{i = 1}^k \frac{d}{d_i} n_i \f$
  where \f$ d \f$ is the product of common and not common factors of
  \f$ d_1, \cdots, d_k \f$ with maximum exponent.
*/
Expr
numerator(const Expr& e) {
  Expr numerator;
  Expr denominator;
  numerator_denominator_purrs(e.distribute_mul_over_add(),
			      numerator, denominator);
  return numerator;
}

/*!
  Let \f$ e(n) = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.
  This function returns the product of common and not common factors of
  \f$ d_1, \cdots, d_k \f$ with maximum exponent.
*/
Expr
denominator(const Expr& e) {
  Expr numerator;
  Expr denominator;
  numerator_denominator_purrs(e.distribute_mul_over_add(),
			      numerator, denominator);
  return denominator;
}

/*!
  Let \f$ e = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.  This
  function returns an expression equivalent to \f$ e \f$ of the form
  \f$ \frac{n}{d} \f$.  Moreover, \f$ d \f$ is the product of common
  and not common factors of \f$ d_1, \cdots, d_k \f$ with maximum
  exponent and \f$ n = \sum_{i = 1}^k \frac{d}{d_i} n_i \f$.
*/
Expr
transform_in_single_fraction(const Expr& e) {
  Expr numerator;
  Expr denominator;
  numerator_denominator_purrs(e.distribute_mul_over_add(),
			      numerator, denominator);
  return numerator / denominator;
}

} // namespace Parma_Recurrence_Relation_Solver
