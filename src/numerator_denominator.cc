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
// See FIXME below.
//#include <algorithm>

// TEMPORARY
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

/*!
  Given \f$ e \f$ a factor of the \f$ i \f$-th term of a more general
  expression \f$ f \f$, this function splits its numerator and denominator
  putting them in the \f$ i \f$-th position of the vectors \p numerators or
  \p denominators, respectively.
  Returns the denominator of \f$ e \f$ while in the \f$ i \f$-th position
  of the vectors \p numerators and \p denominators are built the numerator
  and the denominator of \f$ f \f$.
*/
Expr
find_denominator_single_factor(const Expr& e, unsigned position,
			       std::vector<Expr>& numerators,
			       std::vector<Expr>& denominators) {
  Number num;
  // `e' is a denominator, i. e., a power with exponent a negative
  // rational number.
  if (e.is_a_power() && e.arg(1).is_a_number(num) && num.is_rational()
      && !num.is_positive()) {
    Expr den = pwr(e.arg(0), -num);
    denominators[position] *= den;
    return den;
  }
  // `e' is a rational number.
  else if (e.is_a_number(num) && num.is_rational()) {
    numerators[position] *= num.numerator();
    Expr den = num.denominator();
    denominators[position] *= den;
    return den;
  }
  // `e' is not a denominator and is not a rational number.
  else {
    numerators[position] *= e;
    return 1;
  }
}

/*!
  Given \f$ e \f$ the \f$ i \f$-th term of a more general expression \f$ f \f$,
  this function splits its numerator and denominator
  putting them in the \f$ i \f$-th position of the vectors \p numerators or
  \p denominators, respectively.
  Returns the denominator of \f$ e \f$.
*/
Expr
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
  We consider two factorized expressions \f$ e_1 \f$ and \f$ e_2 \f$
  and two vectors for each expression containing bases and exponents
  of their factors.
  This function builds a new expression with all factors, common and
  not common, of \f$ e_1 \f$ and \f$ e_2 \f$, with the maximum exponent.
*/
Expr
take_common_and_not_factors(std::vector<Expr>& bases_1,
			    std::vector<Expr>& exponents_1,
			    std::vector<Expr>& bases_2,
			    std::vector<Expr>& exponents_2) {
  Expr e = 1;
  // We consider `base_1' the i-th element of `bases_1'.
  for (unsigned i = bases_1.size(); i-- > 0; ) {
    const Expr& base_1 = bases_1[i];
    Number exponent_1;
    // The exponent of `base' is a positive integer number.  
    if (exponents_1[i].is_a_number(exponent_1)
	&& exponent_1.is_positive_integer()) {
      bool added = false;
      for (unsigned j = bases_2.size(); j-- > 0; ) {
	const Expr& base_2 = bases_2[j];
	// We have found an element of `bases_2' sintactically equal
	// to `base_1'.
	if (base_2 == base_1) {
	  Number exponent_2;
	  // Both exponents are positive integer numbers.
	  if (exponents_2[j].is_a_number(exponent_2)
	      && exponent_2.is_positive_integer()) {
	    if (exponent_1 >= exponent_2)
	      e *= pwr(base_1, exponent_1);
	    else
	      e *= pwr(base_1, exponent_2);
	  }
	  // At least one of the two exponents is not a positive
	  // integer number.
	  else {
	    e *= pwr(base_1, exponent_1);
	    e *= pwr(base_2, exponents_2[j]);
	  }
	  // A factor with exponent equal to `0' means that it was
	  // already considered.
	  exponents_2[j] = 0;
	  added = true;
	  break;
	}
      }
      // We have not found an element of `bases_2' equal to `base_1'.
      if (!added)
	e *= pwr(base_1, exponent_1);
    }
    // The exponent of `base_1' is not a positive integer number.
    else
      e *= pwr(base_1, exponents_1[i]);
  }
  // All the factors of the first vectors (`bases_1' and `exponents_1') were already
  // considered. Now we must to consider the factors of the seconds vectors
  // (`bases_2' and `exponents_2') with the exponents not equal to `0', i. e.,
  // not considered in the previous loop because they did not have common factors
  // with the first vectors. 
  for (unsigned i = bases_2.size(); i-- > 0; )
    if (!exponents_2[i].is_zero())
      e *= pwr(bases_2[i], exponents_2[i]);
  return e;
}

void
split_bases_exponents_factor(const Expr& e,
			     std::vector<Expr>& bases,
			     std::vector<Expr>& exponents) {
  Number e_num;
  if (e.is_a_number(e_num) && e_num.is_integer()) {
    std::vector<Number> e_num_bases;
    std::vector<int> e_num_exponents;
    partial_factor(e_num, e_num_bases, e_num_exponents);
#if 1
    for (unsigned i = e_num_bases.size(); i -- > 0; ) {
      bases.push_back(e_num_bases[i]);
      exponents.push_back(e_num_exponents[i]);
    }
#else
    // FIXME: this two rows are not equivalent to those in the #if 1?
    copy(e_num_bases.begin(), e_num_bases.end(), bases.begin());
    copy(e_num_exponents.begin(), e_num_exponents.end(), exponents.begin());
#endif
  }
  else
    if (e.is_a_power()) {
      bases.push_back(e.arg(0));
      exponents.push_back(e.arg(1));
    }
    else {
      bases.push_back(e);
      exponents.push_back(1);
    }
}

/*!
  Given an expression \f$ e \f$, this function returns bases and exponents of
  each factor of \f$ e \f$ in a pair of vectors.
*/
void
split_bases_exponents(const Expr& e,
		      std::vector<Expr>& bases, std::vector<Expr>& exponents) {
  if (e.is_a_mul())
    for (unsigned i = e.nops(); i-- > 0; )
      split_bases_exponents_factor(e.op(i), bases, exponents);
  else
    split_bases_exponents_factor(e, bases, exponents);
}

/*!
  We consider \f$ d = {d_1}^{e_1} \dots {d_k}^{e_k} \f$ and
  \f$ f = {f_1}^{g_1} \dots {f_h}^{g_h} \f$.
  Let \f$ f^g \f$ be the generic factor of \f$ f \f$.
  This function builds a new expression taking all factors, common and
  not common, of the two expressions with the maximum exponent.
  More exactly:
  - if \f$ \exists i \in \{1, \dotsc , k\} \st {d_i}^{e_i} = f^g \f$,
    then returns \f$ d \f$;
  - if \f$ \exists i \in \{1, \dotsc , k\} \st d_i = f \f$, \f$ e_i \f$ and
    \f$ g \f$ are different integer positive numbers,
    if \f$ e_i > g \f$ then returns \f$ d \f$, otherwise
    returns \f$ {d_1}^{e_1} \dots {d_i}^g \dots {d_k}^{e_k}\f$;
  - returns \f$ d f \f$ in all other cases.
*/
Expr
take_common_and_not_factors(const Expr& d, const Expr& f) {
  std::vector<Expr> d_bases;
  std::vector<Expr> d_exponents;
  std::vector<Expr> f_bases;
  std::vector<Expr> f_exponents;
  split_bases_exponents(d, d_bases, d_exponents);
  split_bases_exponents(f, f_bases, f_exponents);
  return take_common_and_not_factors(d_bases, d_exponents,
				     f_bases, f_exponents);
}

/*!
  We consider \f$ d = {b_d}^{e_d} \f$ and \f$ f = {b_f}^{e_f} \f$.
  This function builds a new expression \f$ r \f$ in the following way:
  - if \f$ b_d != b_f \f$ then \f$ r = d \f$;
  - if \f$ b_d = b_f = b \f$ and the two exponent are positive integers,
    then \f$ r = b^{e_d - e_f} \f$;
  - in all other cases \f$ r = 1 \f$.
*/
Expr
find_factor_for_numerator(const Expr& d, const Expr& f) {
  assert(d != f);
  Expr rem = 1;
  std::vector<Expr> d_bases;
  std::vector<Expr> d_exponents;
  std::vector<Expr> f_bases;
  std::vector<Expr> f_exponents;
  split_bases_exponents(d, d_bases, d_exponents);
  split_bases_exponents(f, f_bases, f_exponents);
  for (unsigned i = d_bases.size(); i -- > 0; ) {
    const Expr& d_base = d_bases[i];
    bool added = false;
    for (unsigned j = f_bases.size(); j -- > 0; ) {
      const Expr& f_base = f_bases[j];
      if (d_base == f_base)
	if (d_exponents[i] == f_exponents[j]) {
	  added = true;
	  break;
	}
	else {
	  Number d_exponent;
	  Number f_exponent;
	  if (d_exponents[i].is_a_number(d_exponent)
	      && f_exponents[j].is_a_number(f_exponent)
	      && d_exponent.is_positive_integer()
	      && f_exponent.is_positive_integer()) {
	    rem *= pwr(d_base, d_exponent - f_exponent);
	    added = true;
	    break;
	  }
	  else {
	    rem *= pwr(d_base, d_exponent);
	    added = true;
	    break;
	  }
	}
    }
    if (!added)
      rem *= pwr(d_base, d_exponents[i]);
  }
  return rem;
}

/*!
  Let \f$ e = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.
  The vectors \p numerators and \p denominators contain
  \f$ n_1, \dots, n_k \f$ and \f$ d_1, \dots, d_k \f$, respectively;
  \p denominator contain the common denominator of \f$ d_1, \dots, d_k \f$
  in according with the explanation of the function
  <CODE>take_common_and_not_factors()</CODE>.
  This function returns \f$ \sum_{i = 1}^k \frac{d}{d_i} n_i \f$.
*/
Expr
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
      multiply_to_numerator *= find_factor_for_numerator(denominator,
							 i_th_denominator);
    numerator += numerators[i] * multiply_to_numerator;
  }
  return numerator;
}

/*!
  Let \f$ e \f$ be an expression, this function finds recursively its
  numerator and denominator.
*/
void
numerator_denominator_term(const Expr& e,
			   Expr& numerator, Expr& denominator) {
  // The dimension of `numerators' and `denominators' is terms'number of `e'
  // and will contain the numerator and the denominator of each term of `e'.  
  std::vector<Expr> numerators;
  std::vector<Expr> denominators;
  denominator = 1;
  // Since this function is recursive, even if the first time `e' can not
  // be an `add', for the recursive calls `e' can be an `add'.
  if (e.is_a_add()) {
    numerators.insert(numerators.begin(), e.nops(), Number(1));
    denominators.insert(denominators.begin(), e.nops(), Number(1));
    for (unsigned i = e.nops(); i-- > 0; ) {
      numerator_denominator_term(e.op(i), numerators[i], denominators[i]);
      denominator = take_common_and_not_factors(denominator,
						denominators[i]);
    }
    numerator = find_numerator(numerators, denominators, denominator);
  }
  else if (e.is_a_power()) {
    numerators.push_back(1);
    denominators.push_back(1);
    Expr base_num;
    Expr base_den;
    numerator_denominator_term(e.arg(0), base_num, base_den);
    Expr exp_num;
    Expr exp_den;
    numerator_denominator_term(e.arg(1), exp_num, exp_den);
    denominator
      = find_denominator_single_term(pwr(base_num/base_den, exp_num/exp_den),
				     0, numerators, denominators);
    numerator = numerators[0];
  }
  else if (e.is_a_mul())
    for (unsigned i = e.nops(); i-- > 0; ) {
      const Expr& factor = e.op(i);
      if (factor.is_a_power()) {
	numerators.push_back(1);
	denominators.push_back(1);
	Expr num;
	Expr den;
	numerator_denominator_term(factor, num, den);
	denominator *= find_denominator_single_term(num/den, 0,
						    numerators, denominators);
	numerator = numerators[0];
      }
      else {
	numerators.push_back(1);
	denominators.push_back(1);
	denominator *= find_denominator_single_term(factor, 0,
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

} // anonymous namespace

/*!
  Let \f$ e = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.
  This function returns two expressions:
  in \p denominator puts the product \f$ d \f$ of common and not common
  factors of \f$ d_1, \cdots, d_k \f$ with maximum exponent;
  in \p numerator puts the \f$ \sum_{i = 1}^k \frac{d}{d_i} n_i \f$.
  The expression \f$ n / d \f$ is equivalent to \f$ e \f$.
*/
void
PURRS::numerator_denominator_purrs(const Expr& e,
				   Expr& numerator, Expr& denominator) {
  // The dimension of `numerators' and `denominators' is terms'number of `e'
  // and will contain the numerator and the denominator of each term of `e'.  
  std::vector<Expr> numerators;
  std::vector<Expr> denominators;
  denominator = 1;
  if (e.is_a_add()) {
    numerators.insert(numerators.begin(), e.nops(), Number(1));
    denominators.insert(denominators.begin(), e.nops(), Number(1));
    for (unsigned i = e.nops(); i-- > 0; ) {
      // Find numerator and denominator of i-th term of `e'.
      numerator_denominator_term(e.op(i), numerators[i], denominators[i]);
      // Find common denominator, i.e., the product of common and not common
      // factors with maximum exponent, between `denominator' and the
      // denominator of the i-th term.
      denominator = take_common_and_not_factors(denominator, denominators[i]);
    }
    // Now we have numerator and denominator of each term and also the common
    // denominator, so we can find the numerator.
    numerator = find_numerator(numerators, denominators, denominator);
  }
  else
    // Find numerator and denominator of `e'.
    numerator_denominator_term(e, numerator, denominator);
}

/*!
  Let \f$ e = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.
  This function returns \f$ n = \sum_{i = 1}^k \frac{d}{d_i} n_i \f$
  where \f$ d \f$ is the product of common and not common factors of
  \f$ d_1, \cdots, d_k \f$ with maximum exponent.
*/
PURRS::Expr
PURRS::numerator(const Expr& e) {
  Expr numerator;
  Expr denominator;
  numerator_denominator_purrs(e.distribute_mul_over_add(),
			      numerator, denominator);
  return numerator;
}

/*!
  Let \f$ e = \frac{n_1}{d_1} + \dots + \frac{n_k}{d_k} \f$.
  This function returns the product of common and not common factors of
  \f$ d_1, \cdots, d_k \f$ with maximum exponent.
*/
PURRS::Expr
PURRS::denominator(const Expr& e) {
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
PURRS::Expr
PURRS::transform_in_single_fraction(const Expr& e) {
  Expr numerator;
  Expr denominator;
  numerator_denominator_purrs(e.distribute_mul_over_add(),
			      numerator, denominator);
  return numerator / denominator;
}
