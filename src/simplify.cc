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

#include <config.h>

#include "simplify.hh"

#include "numerator_denominator.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include "Recurrence.defs.hh"
#include <vector>

// TEMPORARY
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

Expr
simplify_expanded_ex_for_input(const Expr& e, bool input);

Expr
simplify_expanded_ex_for_output(const Expr& e, bool input);

Expr
simplify_logarithm_in_expanded_ex(const Expr& e);

//! \brief
//! Given an expression that is a power this function consider its exponent.. 
//! Returns <CODE>true</CODE> if the rule `E2' of the rules's set <EM>Expand</EM>
//! is applicable and in this case \p e will contain the new exponent for
//! the expression; returns <CODE>false</CODE> otherwise.
/*!
  If \f$ e = e_1 \cdots e_k \f$ and
  \f$ \exists i \in \{1, \dotsc , k\} \st e_i = n \f$,
  then this function returns <CODE>true</CODE> and \f$ e \f$ becomes
  \f$ e_1 \cdots e_{i-1} \cdot e_{i+1} \cdots e_k \f$.
  If \f$ e = n \f$ then returns <CODE>true</CODE> and \f$ e \f$
  becomes \f$ 1 \f$.
  For all the other cases returns <CODE>false</CODE> and \f$ e \f$
  does not change.
*/
bool
erase_factor(Expr& e) {
  if (e.is_a_mul()) {
    unsigned num_factors = e.nops();
    unsigned i;
    for (i = 0; i < num_factors; ++i)
      if (e[i] == Recurrence::n)
	break;
    if (i < num_factors) {
      // Found an occurrence of the symbol `n'.
      Expr r = 1;
      for (unsigned j = 0; j < num_factors; ++j)
	if (i != j)
	  r *= e[j];
      e = r;
      return true;
    }
  }
  else if (e == Recurrence::n) {
    e = 1;
    return true;
  }
  return false;
}

//! Application of a part of the rule `E1' of the rules's set <EM>Expand</EM>.
/*!
  Returns <CODE>true</CODE> if \p base and \p exponent are rational numbers
  and \f$ base^exponent \f$ is again a rational number.
  Returns <CODE>false</CODE> otherwise (for example \f$ base = 3 \f$ and
  \f$ exponent = 1/2 \f$).
*/
bool
perfect_root(const Expr& base, const Number& exponent) {
  if (exponent.is_rational()) {
    Number num;
    if (pwr(base, exponent).is_a_number(num) && num.is_rational())
      return true;
  }
  return false;
}

/*!
  Given the base \p base of a power, the numeric and non-numeric part of the
  exponent, \p num_exp and \p not_num_exp respectively, this function
  returns a semantically equivalent form of
  \f$ base^{num\_exp * not\_num\_exp} \f$, according to some conditions 
  that depend on the values of the boolean \p input and \p is_numeric_base.
  If \p input is <CODE>true</CODE> \p n is considered like a special symbol,
  i. e., it is always collected with respect to the others parameters.
  If \p is_numeric_base is <CODE>true</CODE>, the function computes
  explicitly \f$ base^{num\_exp} \f$, and then raises the result to
  the power \f$ not\_num\_exp \f$.
*/
Expr
return_power(bool is_numeric_base, const Expr& base,
	     const Expr& num_exp, const Expr& not_num_exp, bool input) {
  Expr not_num_exp_minus_n = not_num_exp;
  bool n_removed = false;
  if (input)
    n_removed = erase_factor(not_num_exp_minus_n);
  // We do not want to collect the special symbol `n' or it is
  // not in `not_num_exp'.
  if (!input || !n_removed)
    if (is_numeric_base)
      return pwr(pwr(base, num_exp), not_num_exp);
    else
      return pwr(base, num_exp * not_num_exp);
  // We collect the special symbol `n'. 
  else
    if (is_numeric_base)
      return pwr(pwr(pwr(base, num_exp), not_num_exp_minus_n), Recurrence::n);
    else
      return pwr(pwr(base, num_exp * not_num_exp_minus_n), Recurrence::n);
}

/*!
  Separates numeric factors and non-numeric factors of an expression \p e and
  it stores them in \p num and \p not_num, respectively.
*/
void
split_exponent(const Expr& e, Expr& num, Expr& not_num) {
  if (e.is_a_number())
    num *= e;
  else
    not_num *= e;
}

/*!
  Let \p base be a power.
  If \p base is a nested power then we find the "real" base and, at the
  same time, build the exponent that will be obtained from the product of
  \p numeric_exponent and \p not_numeric_exponent.
*/
void
find_real_base_and_build_exponent(Expr& base, Expr& numeric_exponent,
				  Expr& not_numeric_exponent) {
  assert(base.is_a_power());
  while (base.is_a_power()) {
    // Split the exponent in two parts: numeric and non-numeric.
    const Expr& exponent = base.arg(1);
    if (exponent.is_a_mul())
      for (unsigned i = exponent.nops(); i-- > 0; )
        split_exponent(exponent.op(i), numeric_exponent, not_numeric_exponent);
    else
      split_exponent(exponent, numeric_exponent, not_numeric_exponent);
    base = base.arg(0);
  }
}

/*!
  \p base is a <CODE>mul</CODE>.
  We consider three vectors with a number of elements like the factor's
  number of \p base: \p vect_base, \p vect_num_exp and \p vect_not_num_exp.
  Initially \p vect_num_exp and \p vect_not_num_exp will have each element
  equal to \p num_exponent or \p not_num_exponents, respectively.
  The function looks at each factor of \p base and there are two cases:
  -  if it is a power, puts its base in the respective position of the vector
     \p vect_base and updates the respective values of \p vect_num_exp and
     \p vect_not_num_exp with the exponent of \p base's factor;
  -  if it is not a power, puts it in the respective position of the vector
     \p vect_base and leaves unchanged the other vectors.
  If \p input is <CODE>true</CODE> then we use the simplifications
  which collect the special symbol \p n; otherwise, i. e. if \p input
  is <CODE>false</CODE>, \p n is like the other parameters.
*/
Expr
simplify_powers_with_bases_mul(const Expr& base, const Expr& num_exponent,
			       const Expr& not_num_exponent, bool input) {
  assert(base.is_a_mul());
  std::vector<Expr> vect_base(base.nops());
  std::vector<Expr> vect_num_exp(base.nops());
  std::vector<Expr> vect_not_num_exp(base.nops());
  for (unsigned i = base.nops(); i-- > 0; ) {
    vect_num_exp[i] = num_exponent;
    vect_not_num_exp[i] = not_num_exponent;
  }
  for (unsigned i = base.nops(); i-- > 0; ) {
    const Expr& base_factor = base.op(i);
    if (base_factor.is_a_power()) {
      Expr tmp_base = base_factor;
      find_real_base_and_build_exponent(tmp_base, vect_num_exp[i],
					vect_not_num_exp[i]);
      vect_base[i] = tmp_base;
    }
    else
      vect_base[i] = base_factor;
  }
  // The numeric and non-numeric part of the exponent of each factor
  // of the base is determined.
  Expr result = 1;
  for (unsigned i = base.nops(); i-- > 0; )
    if (!vect_base[i].is_a_number())
      result *= return_power(false, vect_base[i],
			     vect_num_exp[i], vect_not_num_exp[i], input);
    else
      result *= return_power(true, vect_base[i],
			     vect_num_exp[i], vect_not_num_exp[i], input);
  return result;
}

/*!
  This function applies the rules \f$ \textbf{E1}, \textbf{E2},
  \textbf{E4} \f$ and \f$ \textbf{E5} \f$ of the rules' set
  <EM>Expand</EM>.  The expression \p e is a power:
  this function finds the base and the exponent of the power (\p e could
  be a series of nested powers). The exponents may contain products but
  not sums, since the expression \p e is expanded.  The function
  separates the exponents into two parts: it puts numeric factors in
  \p num_exponent, and non-numeric factors in \p not_num_exponent.
  Afterwards it tests the base: if it is not a multiplication the
  function returns; otherwise it must raise every factor of the base to
  the exponents.
  If \p input is <CODE>true</CODE> then we use the simplifications
  which collect the special symbol \p n; otherwise, i. e. \p input
  is <CODE>false</CODE>, \p n is like the other parameters.
*/
Expr
simplify_powers(const Expr& e, bool input) {
  assert(e.is_a_power());
  // Accumulate here the numerical part of the exponent.
  Expr num_exponent = 1;
  // Accumulate here the non-numerical part of the exponent.
  Expr not_num_exponent = 1;
  // Since `e' can be a nested power, find the "real" base of `e' and
  // build the exponent for the new power, which is not nested.
  Expr base = e;
  find_real_base_and_build_exponent(base, num_exponent, not_num_exponent);
  D_VAR(base);
  D_VAR(num_exponent);
  D_VAR(not_num_exponent);

  // The base is a multiplication.
  if (base.is_a_mul())
    return simplify_powers_with_bases_mul(base, num_exponent, not_num_exponent,
					  input);
  // The base is not a multiplication: is not necessary to call the function
  // `simplify_powers_with_bases_mul' that uses vectors.
  else
    if (base.is_a_number()) {
      Number num_exp = num_exponent.ex_to_number();
      // The function `perfect_root' allows to apply the rule `E1'.
      if (num_exp.is_integer() || perfect_root(base, num_exp))
	return return_power(true, base, num_exp, not_num_exponent, input);
      else
	return return_power(false, base, num_exp, not_num_exponent, input);
    }
    else
      return return_power(false, base, num_exponent, not_num_exponent, input);
}

/*!
  Returns in the vectors \p bases and \p exponents the bases and the
  exponents of the eventual multiplication's factors which are powers
  contained in \p e.
*/
void
find_bases_and_exponents(const Expr& e, std::vector<Expr>& bases,
			 std::vector<Expr>& exponents) {
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& factor = e.op(i);
    if (factor.is_a_power()) {
      bases.push_back(factor.arg(0));
      exponents.push_back(factor.arg(1));
    }
  }
}

//! Applies the rules \f$ \textbf{C3} \f$ of the set of rules <EM>Collect</EM>.
/*!
  Applies the rules \f$ \textbf{C3} \f$ of the set of rules
  <EM>Collect</EM> to the expression \p e under the condition
  that the common exponent of the powers is not an integer
  because, in this case, the power is automatically decomposed,
  i. e., \f$ (a*b)^4 \f$ is automatically transformed in
  \f$ a^4*b^4 \f$.
  Returns a new expression \p e_rewritten containing the
  modified expression \p e.
*/
Expr
collect_same_exponents(const Expr& e) {
  assert(e.is_a_mul());
  // Builds two vectors containing the bases and the exponents of
  // the eventual multiplication's factors which are powers. 
  std::vector<Expr> bases;
  std::vector<Expr> exponents;
  find_bases_and_exponents(e, bases, exponents);
  Expr e_rewritten = 1;
  // At the end of the following cycle, `e_rewritten' will contain the powers
  // of `e' with the same exponents simplified in only one power with the
  // base equal to the previous bases multiplicated among themselves 
  // (rule `C3').
  for (unsigned i = exponents.size(); i-- > 0; ) {
    for (unsigned j = i; j-- > 0; ) {
      if (exponents[j] == exponents[i]) {
	// We have found two equal exponents. We store in the i-th
	// position of the vectors `bases' and `exponents' the new power
	// with base equal to the product of the base in the i-th and in
	// the j-th position and exponent unchanged.
	bases[i] *= bases[j];
	exponents[j] = 0;
      }
      else if (exponents[j] == -exponents[i]) {
	// We have found two opposite exponents. We store in the i-th
	// position of the vectors `bases' and `exponents' the new power.
	// If the i-th exponent is more brief than the j-th then the new
	// power has the base equal to the ratio of the i-th and the j-th
	// base and exponent equal to the i-th.
	// If the j-th exponent is more brief than the i-th then the new
	// power has the base equal to the ratio of the j-th base and the
	// i-th and exponent equal to the j-th.
	if (exponents[j].nops() > exponents[i].nops())
	  bases[i] /= bases[j];
	else {
	  bases[i] = bases[j] / bases[i];
	  exponents[i] = exponents[j];
	}
	exponents[j] = 0;
      }
      // We have marked with `0' the j-th position of the vector `exponents' in
      // order to remember to us that it contains garbage.
    }
    bases[i] = simplify_numer_denom(bases[i]);
    e_rewritten *= pwr(bases[i], exponents[i]);
  }
  // Now adds to `e_rewritten' the factors of `e' not considered in the
  // previous cycle, i.e., the factors which are not powers.
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& factor = e.op(i);
    if (!factor.is_a_power())
      e_rewritten *= factor;
  }
  return e_rewritten;
}

//! \brief
//! Applies the rules \f$ \textbf{C1} \f$ and \f$ \textbf{C2} \f$ of the
//! set of rules <EM>Collect</EM>.
/*!
  Applies the rules \f$ \textbf{C1} \f$ and \f$ \textbf{C2} \f$ of the set
  of rules <EM>Collect</EM> to the <CODE>Expr</CODE> \p e that is
  certainly a <CODE>mul</CODE>.
  It returns a new expression \p e_rewritten containing the modified
  expression \p e.
*/
Expr
collect_same_base(const Expr& e) {
  assert(e.is_a_mul());
  // Builds two vectors containing the bases and the exponents of
  // the eventual multiplication's factors which are powers. 
  std::vector<Expr> bases;
  std::vector<Expr> exponents;
  find_bases_and_exponents(e, bases, exponents);
  Expr e_rewritten = 1;
  // At the end of the following cycle, `e_rewritten' will contain all 
  // the powers of `e' with the same bases simplified in only one power
  // with exponent equal to the sum of the previous powers' exponents
  // (rule `C2').
  for (unsigned i = exponents.size(); i-- > 0; ) {
    for (unsigned j = i; j-- > 0; )
      if (bases[j] == bases[i]) {
	// We have found two equal bases. We store in the i-th
	// position of the vectors `bases' and `exponents' the new power
	// with exponent equal to the sum of the exponent in the i-th and in
	// the j-th position and base unchanged.
	// We mark with `0' the j-th position of the vector `exponents' in
	// order to remember to us that it contains garbage.
	exponents[i] += exponents[j];
	exponents[j] = 0;
      }
    e_rewritten *= pwr(bases[i], exponents[i]);
  }
  // Now adds to `e_rewritten' the factors of `e' not considered in the
  // previous cycle, i.e., the factors which are not powers applying,
  // when possible, the rule `C1'.
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& factor_e = e.op(i);
    if (!factor_e.is_a_power()) {
      // We must consider those factors that are not powers but are equal
      // to some base of the vector `bases', in order to apply the rule `C1',
      // i.e., `a * a^e = a^(e + 1)'. In this case, we add `1' to the
      // correspondenting exponent.
      bool to_sum = false;
      Expr old_factor;
      Expr new_factor;
      for (unsigned j = bases.size(); j-- > 0; ) {
	const Expr& bases_j = bases[j];
	const Expr& exponents_j = exponents[j];
	if (bases_j == factor_e && !exponents_j.is_zero()) {
	  to_sum = true;
  	  old_factor = pwr(bases_j, exponents_j);
	  Number base;
	  // If the base is an integer number then automatically the expression
	  // is transformed, for instance, `2^(3/2)' is transformed in
	  // `2*sqrt(2)': in this case we do not add 1 to the exponent.
	  if (bases_j.is_a_number(base)
	      && (!base.is_integer() || base == -1))
	    new_factor = pwr(bases_j, exponents_j + 1);
	  else if (!bases_j.is_a_number() || !exponents_j.is_a_number())
	    new_factor = pwr(bases_j, exponents_j + 1);
	  else
	    new_factor = pwr(bases_j, exponents_j) * factor_e;
	}
      }
      // Applies rule `C1'.
      if (to_sum)
	e_rewritten = e_rewritten.substitute(old_factor, new_factor);
      else
	e_rewritten *= factor_e;
    }
  }
  return e_rewritten;
}

//! Applies the rules of the set of rules <EM>Collect</EM>.
/*!
  Applies the rules of the set <EM>Collect</EM> to an expression
  \p e that is certainly a <CODE>mul</CODE>.
  Returns a new expression containing the modified expression \p e. 
*/
Expr
collect_base_exponent(const Expr& e) {
  assert(e.is_a_mul());
  Expr e_rewritten = e;

  // We have a better simplification if we apply `collect_same_exponents()'
  // before than `collect_same_base()' (ex. `2*2^n*(1/2)^n').
  // Applies rules `C3'.
  e_rewritten = collect_same_exponents(e_rewritten);
  D_MSGVAR("e_rewritten dopo same exponents: ", e_rewritten);

  // After the simplifications by the function `collect_same_exponents()'
  // `e_rewritten' could not be a `mul'. 
  if (e_rewritten.is_a_mul())
    // Applies rules `C1' and `C2'.    
    e_rewritten = collect_same_base(e_rewritten);
  D_MSGVAR("e_rewritten dopo same base: ", e_rewritten);

  return e_rewritten;
}

/*!
  Compute the standard form for \f$ n^{1/k} \f$, where \f$ k \f$ is a non-zero 
  integer (possibly negative), and \f$ n \f$ is an integer whose partial 
  factorization \f$ n = b_1^{e_1} \cdots b_r^{e_r} \f$ computed by 
  <CODE>partial_factor()</CODE> is given. 
  If \f$ k > 0 \f$, the exponents in the vector \p exponents are 
  reduced modulo \f$ k \f$, and the value 
  \f$ m = b_1^{[e_1/k]} \cdots b_r^{[e_r/k]} \f$ is returned. 
  If \f$ k < 0 \f$ (that is, if the number \f$ n \f$ appears in the 
  denominator of the fraction whose \f$ k \f$-th root we are trying to 
  reduce in normal form) and \f$ n \f$ factors as above, exponents and 
  remainders are adjusted so that the possible irrational part will only 
  appear in the numerator of the final result. 
  In particular, if for example \f$ k \f$ does not divide \f$  e_1 \f$, 
  instead of returning the factor \f$ b_1^{[e_1/|k|]} \f$ and replacing 
  the exponent \f$ e_1 \f$ by \f$ e_1 \mod |k| \f$, we return 
  \f$ b_1^{1 + [e_1/|k|]} \f$ and replace the exponent by 
  \f$ |k| - (e_1 \mod |k|) \f$. 
*/
Expr
to_std_form(const Number& k, const std::vector<Number>& bases, 
	    std::vector<int>& exponents) {
  
  assert(k != 0);
  assert(k.is_integer());
  int abs_k = ::abs(k.to_int());
  Expr m = 1;
  for (unsigned i = bases.size(); i-- > 0; ) {
    int remainder = exponents[i] % abs_k;
    int quotient  = exponents[i] / abs_k;
    if (k < 0 && remainder != 0) { // adjust quotient and remainder
      ++quotient;
      remainder = abs_k - remainder;
    }
    exponents[i] = remainder;
    m *= pwr(bases[i], quotient);
  }
  return m;
}

/*!
  Computes a `simple' form for \f$ r^{1/k} \f$, where \p r is a
  rational number, and \p k is a non-zero (possibly negative) integer.
  We need to define the <EM>standard form</EM> of \f$ r^{1/k} \f$,
  which we denote by \f$ \stdform\bigl( r^{1/k} \bigr) \f$, so that
  \f$ r^{1/k} = \stdform\bigl( r^{1/k} \bigr) \f$ and the expression on
  the right is uniquely defined and as simple as possible.
  The standard form will have the shape of the product of a complex \p sign
  such that \f$ sign \in \{1, -1, i, -i\} \f$, a positive rational
  number \f$ r_1 \f$ and \f$ m^{1/j} \f$, where \f$ m \f$ is
  a positive integer which is \f$ j \f$-free (i. e., it is not divisible
  by the \f$ j \f$-th power of any integer larger than 1),
  and \f$ j \f$ is the smallest factor of \p k for which such an integer
  \f$ m \f$ exists.

  The definition is then as follows: we assume that
  \f$ r = num / den \f$ where \p num and \p den are non-zero, coprime integers.
  - If \p num and \p den are both positive, and \f$ k > 0 \f$, then
    \f[
      \stdform \left( \Bigl( \frac{num}{den} \Bigr)^{1/k} \right)
      \defrel{=}
      \frac{n_1}{d_1} \cdot m^{1/j},
    \f]
    where \f$ n_1, d_1 \f$ and \f$ m \f$ are positive integers, \f$ j \f$ is a
    divisor of \p k, \f$ m \f$ is \f$ j \f$-free, and \f$ m \f$ is not a
    \f$ j_1 \f$-th power, for any positive integer \f$ j_1 \f$ with
    \f$ \gcd(j, j_1) > 1 \f$.
  - If \p num and \p den are both positive, and \f$ k < 0 \f$, then
    \f[
      \stdform \left( \Bigl( \frac{num}{den} \Bigr)^{1/k} \right)
      \defrel{=}
      \stdform \left( \Bigl( \frac{den}{num} \Bigr)^{1/(-k)} \right).
    \f]
  - If \f$ num / den < 0 \f$ and \p k is odd, then
    \f[
    \stdform \left( \Bigl( \frac{num}{den} \Bigr)^{1/k} \right)
    \defrel{=}
    - \stdform \left( \Bigl| \frac{num}{den} \Bigr|^{1/k} \right).
    \f]
  - If \f$ num / den < 0 \f$ and \p k is even, then
    \f[
    \stdform \left( \Bigl( \frac{num}{den} \Bigr)^{1/k} \right)
    \defrel{=}
    \pm i \cdot \stdform \left( \Bigl| \frac{num}{den} \Bigr|^{1/k} \right),
    \f]
    where the sign is chosen according to the sign of \p k.
*/
Expr 
reduce_to_standard_form(const Number& root_index, const Number& r) {
  assert(root_index.is_integer());
  assert(root_index != 0);
  int k = root_index.to_int();
  if (k < 0)
    assert(r != 0);
  if (k > 0)
    if (r == 0)
      return 0;
  // FIXME: deal with complex numbers
  if (!r.is_real()) {
    Expr index = 1 / root_index;
    return pwr(r, index);
  }

  Number sign = r > 0 ? 1 : -1;
  Number num = r.numerator();
  Number den = r.denominator();
  Number g = gcd(num, den);
  num /= g;
  den /= g; // clear any common factors from num, den.
  if (k % 2 == 0 && sign == -1) // complex sign if k is even and r < 0.
    sign = k > 0 ? Number::I : -Number::I;
  if (k < 0) { // swap numerator and denominator, and change sign of k.
    Number i = num;
    num = den;
    den = i;
    k *= -1;
  } // now, num, den and k are all positive.
  if (k == 1)
    return sign * num * pwr(den, -1);
  
  std::vector<Number> num_bases;
  std::vector<int> num_exponents;
  std::vector<Number> den_bases;
  std::vector<int> den_exponents;  
  // Partial factor and reduce numerator and denominator.
  partial_factor(num, num_bases, num_exponents);
  unsigned num_size = num_bases.size();
  Expr reduced_num = to_std_form(k, num_bases, num_exponents);
  // Here <CODE>to_std_form</CODE> is called with a negative value of k 
  // because we are dealing with the denominator of r.
  partial_factor(den, den_bases, den_exponents);
  unsigned den_size = den_bases.size();
  Expr reduced_den = to_std_form(-k, den_bases, den_exponents);
  
  // Try one last simplification: if all exponents have a common factor 
  // with the root index, remove it.
  int gc = k;
  for (unsigned i = 0; i < num_size && gc > 1; ++i)
    gc = gcd(gc, num_exponents[i]);
  for (unsigned i = 0; i < den_size && gc > 1; ++i)
    gc = gcd(gc, den_exponents[i]);
  
  if (gc > 1) {
    k /= gc;
    for (unsigned i = 0; i < num_size; ++i)
      num_exponents[i] /= gc;
    for (unsigned i = 0; i < den_size; ++i)
      den_exponents[i] /= gc;
  }
  
  // The object `irr_part' is surely numeric but we have to use an expression
  // because otherwise the number, since it is irrational, is rounded.
  Expr irr_part = 1;
  for (unsigned i = 0; i < num_size; ++i)
    irr_part *= pwr(num_bases[i], num_exponents[i]);
  for (unsigned i = 0; i < den_size; ++i)
    irr_part *= pwr(den_bases[i], den_exponents[i]);
  Expr q = sign * reduced_num * pwr(reduced_den, -1);
  if (irr_part.ex_to_number() > 1)
    q *= pwr(irr_part, Number(1, k));
  return q;
}

/*!
  The purpose of this function is the computation of the standard form
  of the product of two irrational numbers \f$ a_1 = b_1^{n_1 / d_1} \f$
  and \f$ a_2 = b_2^{n_2 / d_2} \f$, where \f$ b_1 \f$ and \f$ b_2 \f$
  are numerics (not necessarily rational numbers), and \f$ n_1 \f$,
  \f$ d_1 \f$, \f$ n_2 \f$ and \f$ d_2 \f$ are (possibly negative)
  integers.
  The arguments contain \f$ b_1 \f$, \f$ n_1 / d_1 \f$, \f$ b_2 \f$,
  \f$ n_2 / d_2 \f$, respectively.
  \f$ a_1 \f$ and \f$ a_2 \f$ are transformed in two equivalent
  irrational numbers with the same denominator in the exponents
  (reduction to same index of two roots), and then passed to 
  <CODE>reduce_to_standard_form()</CODE>.
*/
Expr
red_prod(const Number& base1, const Number& exp1, 
	 const Number& base2, const Number& exp2) {
  assert(exp1 != 0);
  assert(exp2 != 0);
  assert(exp1.is_rational());
  assert(exp2.is_rational());
  Number base_1 = base1;  
  Number base_2 = base2;
  Number k1_num = exp1.numerator();
  Number k1_den = exp1.denominator();
  Number k2_num = exp2.numerator();
  Number k2_den = exp2.denominator();
  // We want that the sign of the exponent is stored in the denominator
  // (while the methods `numerator()' and `denominator()' store the sign
  // in the numerator).
  // This is needed in the function `reduce_to_standard_form()'.
  if (!k1_num.is_positive()) {
    k1_num *= -1;
    k1_den *= -1;
  }
  if (!k2_num.is_positive()) {
    k2_num *= -1;
    k2_den *= -1;
  }
  base_1 = pwr(base_1, k1_num);
  base_2 = pwr(base_2, k2_num);
  Number g = gcd(k1_den, k2_den);
  Number k = k1_den * k2_den / g;
  Number b1 = pwr(base_1, k2_den / g);
  Number b2 = pwr(base_2, k1_den / g);
  Number b = b1 * b2;
  return reduce_to_standard_form(k, b);
}

/*!
  Applies the rules of the set <EM>Irrationals</EM> to each \p e's factor
  which is a <CODE>power</CODE>.
  This function is called to simplify products: if any factor within
  the product is a power with numeric base and exponent, it calls 
  <CODE>red_prod()</CODE>, and keeps accumulating partial results
  in the pair of expressions <CODE>base_1</CODE> and <CODE>exp_1</CODE>.
  It is guaranteed that <CODE>base_1^exp_1</CODE> is always in its
  standard form, as defined in the comment to the function
  <CODE>reduce_to_standard_form()</CODE>.
*/
Expr
reduce_product(const Expr& e) {
  assert(e.is_a_mul());
  Expr factor_to_reduce = 1;
  Expr factor_no_to_reduce = 1;
  Number base_1 = 1;
  Number exp_1 = 1;
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& factor_e = e.op(i);
    if (factor_e.is_a_power()) {
      // Base and exponent of `factor_e' are both numerics.
      Number base_2;
      Number exp_2;
      if (factor_e.arg(0).is_a_number(base_2) &&
	  factor_e.arg(1).is_a_number(exp_2)) {
	Expr to_reduce = red_prod(base_1, exp_1, base_2, exp_2);
	// red_prod returns `numeric' or `numeric^numeric' or
	// `numeric * numeric^numeric'.
	if (to_reduce.is_a_mul()) {
	  assert(to_reduce.nops() == 2); 
	  for (unsigned j = 2; j-- > 0; ) {
	    const Expr& factor = to_reduce.op(j);
	    if (factor.is_a_power()) {
	      const Expr& base_factor = factor.arg(0);
	      const Expr& exponent_factor = factor.arg(1);
	      base_1 = base_factor.ex_to_number();
	      exp_1  = exponent_factor.ex_to_number();
	      factor_to_reduce = pwr(base_factor, exponent_factor);
	    }
	    else
	      factor_no_to_reduce *= factor;
	  }
	}
	else if (to_reduce.is_a_power()) {
	  const Expr& base_to_reduce = to_reduce.arg(0);
	  const Expr& exponent_to_reduce = to_reduce.arg(1);
	  base_1 = base_to_reduce.ex_to_number();
	  exp_1 = exponent_to_reduce.ex_to_number();
	  factor_to_reduce = pwr(base_to_reduce, exponent_to_reduce);
	}
	else {
	  base_1 = to_reduce.ex_to_number();
	  exp_1 = 1;
	  factor_to_reduce = pwr(to_reduce, 1);
	}
      }
      else
	// Base and exponent of `factor_e' are not both numerics.
	factor_no_to_reduce *= factor_e;
    }
    // `factor_e' is not a `power'.
    else
      factor_no_to_reduce *= factor_e;
  }
  return factor_to_reduce * factor_no_to_reduce;
}

/*!
  Applies all rules of term rewriting system \f$ \mathfrak{R}_o \f$ which
  are applicable on factors.
  Returns an expression containing the modified expression \p e.
*/
Expr
manip_factor(const Expr& e, bool input) {
  assert(e.is_a_mul());
  Expr e_rewritten = 1;
  // Simplifies each factor that is a `power'.
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& factor_e = e.op(i);
    if (factor_e.is_a_power()) {
      Expr base = simplify_expanded_ex_for_output(factor_e.arg(0), input);
      Expr exp = simplify_expanded_ex_for_output(factor_e.arg(1), input);
      Number num_base;
      Number num_exp;
      if (base.is_a_number(num_base) && exp.is_a_number(num_exp))
	e_rewritten *= red_prod(num_base, num_exp, 1, 1);
      else {
	Expr tmp = pwr(base, exp);
	if (tmp.is_a_power())
	  e_rewritten *= simplify_powers(tmp, input);
	else
	  e_rewritten *= tmp;
      }
    }
    else
      e_rewritten *= factor_e;
  }
  D_MSGVAR("e_rewritten dopo nested: ", e_rewritten);
  // From this time forward we do not know if `e_rewritten' is a again `mul'.
  // Simplifies recursively the factors which are functions simplifying
  // their arguments.
  if (e_rewritten.is_a_mul()) {
    Expr factor_function = 1;
    Expr factor_no_function = 1;
    for (unsigned i = e_rewritten.nops(); i-- > 0; ) {
      const Expr& factor_e_rewritten = e_rewritten.op(i);
      if (factor_e_rewritten.is_a_function())
	if (factor_e_rewritten.nops() == 1)
	  factor_function
	    *= apply(factor_e_rewritten.functor(),
		     simplify_expanded_ex_for_output(factor_e_rewritten.arg(0),
						     input));
	else {
	  unsigned num_argument = factor_e_rewritten.nops();
	  std::vector<Expr> argument(num_argument);
	  for (unsigned i = 0; i < num_argument; ++i)
	    argument[i]
	      = simplify_expanded_ex_for_output(factor_e_rewritten.arg(i),
						input);
	  factor_function *= apply(factor_e_rewritten.functor(), argument);
	}
      else
	factor_no_function *= factor_e_rewritten;
    }
    e_rewritten = factor_function * factor_no_function;
  }
  else if (e_rewritten.is_a_function()) {
    if (e_rewritten.nops() == 1)
      e_rewritten = apply(e_rewritten.functor(),
			  simplify_expanded_ex_for_output(e_rewritten.arg(0),
							  input));
    else {
      unsigned num_argument = e_rewritten.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned i = 0; i < num_argument; ++i)
	argument[i] = simplify_expanded_ex_for_output(e_rewritten.arg(i), input);
      e_rewritten = apply(e_rewritten.functor(), argument);
    }
  }
  D_MSGVAR("e_rewritten dopo function: ", e_rewritten);
  // Special case: the exponential `exp' is a `function' but it has
  // the same properties of the powers.
  if (e_rewritten.is_a_mul()) {
    Expr argument = 0;
    Expr rem = 1;
    for (unsigned i = e_rewritten.nops(); i-- > 0; ) {
      const Expr& factor_e_rewritten = e_rewritten.op(i);
      if (factor_e_rewritten.is_the_exp_function())
	argument += factor_e_rewritten.arg(0);
      else
	rem *= factor_e_rewritten;
    }
    e_rewritten = exp(argument) * rem;
    D_MSGVAR("e_rewritten dopo `exp': ", e_rewritten);
  }
  // Simplifies eventual product of irrational numbers.
  if (e_rewritten.is_a_add()) {
    Expr terms = 0;
    for (unsigned i = e_rewritten.nops(); i-- > 0; ) { 
      const Expr& term_e_rewritten = e_rewritten.op(i);
      if (term_e_rewritten.is_a_mul())
	terms += reduce_product(term_e_rewritten);
      else
	terms += term_e_rewritten;
    }
    e_rewritten = terms;
  }
  else
    if (e_rewritten.is_a_mul())
      e_rewritten = reduce_product(e_rewritten);
  // Simplifies eventual powers with same base or same exponents applying
  // the rule of the term rewriting system \f$ \mathfrak{R}_o \f$.
  if (e_rewritten.is_a_add()) {
    Expr terms = 0;
    for (unsigned i = e_rewritten.nops(); i-- > 0; ) { 
      const Expr& term_e_rewritten = e_rewritten.op(i);
      if (term_e_rewritten.is_a_mul())
	terms += collect_base_exponent(term_e_rewritten);
      else
	terms += term_e_rewritten;
    }
    e_rewritten = terms;
  }
  else
    if (e_rewritten.is_a_mul())
      e_rewritten = collect_base_exponent(e_rewritten);
  return e_rewritten;
}

/*!
  Crosses the tree of the expanded expression \p e recursively to find
  subexpressions to which we apply the rules of the terms rewriting
  system \f$ \mathfrak{R}_i \f$. More exactly here the rules of the set
  <EM>Expand</EM> are implemented because the rules of the set
  <EM>Automatic</EM> are automatically executed.  We remark that the
  rule \f$ \textbf{E4} \f$ is automatically executed when the exponent is
  integer, while the rules \f$ \textbf{E3} \f$ and \f$ \textbf{E6} \f$
  are executed by the method <CODE>expand()</CODE> (\f$ \textbf{E3} \f$
  only partially because for instance \f$ expand(3^(4*x+2*a)) =
  3^(2*a)*3^(4*x) \f$): hence we here consider only powers.
  \p input is always <CODE>true</CODE> and this means that \p n is a
  special symbol, i. e., in the simplifications it is always collected
  with respect to the others parameters.  The function returns a
  new expression containing the modified expression \p e.
*/
Expr
simplify_expanded_ex_for_input(const Expr& e, bool input) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten += simplify_expanded_ex_for_input(e.op(i), input);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten *= simplify_expanded_ex_for_input(e.op(i), input);
  }
  else if (e.is_a_power())
    return simplify_powers(e, input);
  else if (e.is_a_function()) {
    if (e.nops() == 1)
      return apply(e.functor(), simplify_expanded_ex_for_input(e.arg(0), input));
    else {
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned i = 0; i < num_argument; ++i)
	argument[i] = simplify_expanded_ex_for_input(e.arg(i), input);
      return apply(e.functor(), argument);
    }
  }
  else
    return e;
  return e_rewritten;
}

/*!
  Crosses the tree of the expanded expression \p e recursively to find
  subexpressions to which we want to apply the rules of the terms
  rewriting system \f$ \mathfrak{R}_o \f$. The remarks about the
  function <CODE>simplify_expanded_ex_for_input()</CODE> are valid here too,
  because all the rules of the term rewriting system \f$ \mathfrak{R}_i
  \f$, except for \f$ \textbf{E2} \f$, are also in \f$ \mathfrak{R}_o
  \f$.  \p input is always <CODE>false</CODE> and this means that \p n
  is not considerated like a special symbol but like the others
  parameters.  The function returns a new expression containing
  the modified expression \p e.
*/
Expr
simplify_expanded_ex_for_output(const Expr& e, bool input) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten += simplify_expanded_ex_for_output(e.op(i), input);
  }
  else if (e.is_a_mul())
    // We can not call `simplify_expanded_ex_for_output()' on every factor
    // because otherwise it is not possible to transform products.
    return manip_factor(e, input);
  else if (e.is_a_power()) {
    // Is necessary to call `red_prod()' in order to rewrite irrational
    // numbers, i.e., powers with base and exponent both numerics. 
    Expr base = simplify_expanded_ex_for_output(e.arg(0), input);
    Expr exp = simplify_expanded_ex_for_output(e.arg(1), input);
    Number num_base;
    Number num_exp;
    if (base.is_a_number(num_base) && exp.is_a_number(num_exp))
      e_rewritten = red_prod(num_base, num_exp, 1, 1);
    else {
      Expr tmp = pwr(base, exp);
      if (tmp.is_a_power())
	e_rewritten = simplify_powers(tmp, input);
      else
	e_rewritten = tmp;
    }
    // Necessary for the output: for example if `e = sqrt(18)^a' then
    // with the following rows `e_rewritten = (3*sqrt(2))^a'
    // (without the following rows `e_rewritten = sqrt(2)^a*3^a').
    if (e_rewritten.is_a_mul())
      return collect_base_exponent(e_rewritten);
  }
  else if (e.is_a_function()) {
    if (e.nops() == 1)
      return apply(e.functor(), simplify_expanded_ex_for_output(e.arg(0), input));
    else {
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned i = 0; i < num_argument; ++i)
	argument[i] = simplify_expanded_ex_for_output(e.arg(i), input);
      return apply(e.functor(), argument);
    }
  }
  else
    return e;
  return e_rewritten;
}

Expr
rewrite_factorial(const Number& a, const Number& b) {
  Expr prod = factorial(a * Recurrence::n);
  if (b > 0)
    for (Number j = 1; j <= b; ++j)
      prod *= a * Recurrence::n + j; 
  else
    for (Number j = 0; j < abs(b); ++j)
      prod *= pwr(a * Recurrence::n - j, -1);
  return prod;
}

bool
check_form_of_mul(const Expr& e, Number& a) {
  assert(e.is_a_mul() && e.nops() == 2);
  const Expr& first = e.op(0);
  const Expr& second = e.op(1);
  if ((first == Recurrence::n && second.is_a_number(a)
       && a.is_positive_integer())
      || (second == Recurrence::n && first.is_a_number(a)
	  && a.is_positive_integer()))
    return true;
  else
    return false;
}

/*!
  Given the factorial expression \f$ (a n + b)! \f$, this function 
  returns a new expression that contains explicitly \f$ (a n)! \f$.
  We use the rewrite rule explained in the comment for
  <CODE>rewrite_factorials()</CODE>
*/
Expr
decompose_factorial(const Expr& e) {
  assert(e.is_the_factorial_function());
  const Expr& argument = e.arg(0);
  if (argument.is_a_add() && argument.nops() == 2) {
    const Expr& first = argument.op(0);
    const Expr& second = argument.op(1);
    Number b;
    if (first.is_a_number(b) && b.is_integer())
      if (second.is_a_mul() && second.nops() == 2) {
	Number a;
	// Checks if `second' has the form `a*n with `a' a positive
	// integer number'. 
	if (check_form_of_mul(second, a))
	  return rewrite_factorial(a, b);
      }
      else if (second == Recurrence::n)
	return rewrite_factorial(1, b);
    if (second.is_a_number(b) && b.is_integer())
      if (first.is_a_mul() && first.nops() == 2) {
	Number a;
	// Checks if `second' has the form `a*n with `a' a positive
	// integer number'. 
	if (check_form_of_mul(first, a))
	  return rewrite_factorial(a, b);
      }
      else if (first == Recurrence::n)
	return rewrite_factorial(1, b);
  }
  return e;
}

/*!
  Applies the following rewrite rules:
  - \f[
      a^{b n + c} = a^{b n} \cdot a^c;
    \f]
  - let \f$ a \in \Nset \setminus \{0\} \f$ and let \f$ b \in \Zset \f$:
    \f[
      \begin{cases}
        (a n + b)!
        =
        (a n)! \cdot (a n + 1) \cdots (a n + b),
          \quad \text{if } b \in \Nset \setminus \{0\}; \\
        (a n)!,
          \quad \text{if } b = 0; \\
        \dfrac{(a n)!}{(a n) \cdot (a n - 1) \cdots (a n + b + 1))},
          \quad \text{if } b \in \Zset \setminus \Nset.
      \end{cases}
    \f]
*/
Expr
rewrite_factorials_and_exponentials(const Expr& e) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten += rewrite_factorials_and_exponentials(e.op(i));
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten *= rewrite_factorials_and_exponentials(e.op(i));
  }
  else if (e.is_a_power()) {
    e_rewritten = pwr(rewrite_factorials_and_exponentials(e.arg(0)),
		      rewrite_factorials_and_exponentials(e.arg(1)));
    if (e_rewritten.is_a_power() && e_rewritten.arg(1).is_a_add())
      return e_rewritten.expand();
  }
  else if (e.is_a_function())
    if (e.is_the_factorial_function())
      return decompose_factorial(e);
    else
      if (e.nops() == 1)
	return apply(e.functor(),
		     rewrite_factorials_and_exponentials(e.arg(0)));
      else {
	unsigned num_argument = e.nops();
	std::vector<Expr> argument(num_argument);
	for (unsigned i = 0; i < num_argument; ++i)
	  argument[i] = rewrite_factorials_and_exponentials(e.arg(i));
	return apply(e.functor(), argument);
      }
  else
    return e;
  return e_rewritten;
}

Expr
change_base_logarithm(const Expr& base, const Expr exponent,
		      const Number& new_base, const Number& new_exponent = 1) {
  if (exponent.is_a_mul() && exponent.nops() == 2) {
    if (exponent.op(0).is_the_log_function()
	&& exponent.op(1) == 1 / log(new_base))
      return pwr(exponent.op(0).arg(0), new_exponent);
    if (exponent.op(1).is_the_log_function()
	&& exponent.op(0) == 1 / log(new_base))
      return pwr(exponent.op(1).arg(0), new_exponent);
  }
  return pwr(simplify_logarithm_in_expanded_ex(base),
	     simplify_logarithm_in_expanded_ex(exponent));
}

void
find_log_factors(const Expr& e, Expr& log_factors, Expr& rem_factors) {
  assert(e.is_a_mul());
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& factor = e.op(i);
    if ((factor.is_a_power() && factor.arg(0).is_the_log_function()
	 && factor.arg(1) == -1)
	|| factor.is_the_log_function())
      log_factors *= factor;
    else
      rem_factors *= factor;
  }
}

Expr
prepare_change_base_logarithm(const Expr& base, const Expr& exponent) {
  Number base_num;
  // The base is a positive number.
  if (base.is_a_number(base_num) && base_num.is_positive()) {
    std::vector<Number> bases;
    std::vector<int> exponents;
    // The base is integer.
    if (base_num.is_integer() && exponent.is_a_mul()) {
      Expr log_factors_exp = 1;
      Expr rem_factors_exp = 1;
      find_log_factors(exponent, log_factors_exp, rem_factors_exp);
      partial_factor(base_num, bases, exponents);
      if (bases.size() == 1) {
	Number new_base = bases[0];
	return pwr(change_base_logarithm(base, log_factors_exp, new_base,
					 exponents[0]), rem_factors_exp);
      }
      else
	return pwr(change_base_logarithm(base, log_factors_exp, base_num),
		   rem_factors_exp);
    }
    // The base is rational with the numerator equal to `1'.
    if (base_num.is_rational() && base_num.numerator() == 1
	&& exponent.is_a_mul()) {
      Number denom_base = base_num.denominator();
      Expr log_factors_exp = 1;
      Expr rem_factors_exp = 1;
      find_log_factors(exponent, log_factors_exp, rem_factors_exp);
      partial_factor(denom_base, bases, exponents);
      if (bases.size() == 1) {
	Number new_base = bases[0];
	return pwr(change_base_logarithm(base, log_factors_exp, new_base,
					 exponents[0]),
		   -rem_factors_exp);
      }
      else
	return pwr(change_base_logarithm(base, log_factors_exp, base_num),
		   -rem_factors_exp);
    }
  }
  return pwr(simplify_logarithm_in_expanded_ex(base),
	     simplify_logarithm_in_expanded_ex(exponent));
}

/*!
  Applies the following logarithm's property:
  \f[
    \begin{cases}
      log(exp(1)^a) = a, \\
      log(a^b) = b log(a), \\
      (a^b)^{log c / log a} = c^b, \quad \text{where } a \in \Rset, a > 0
        \text{and } b \in \Nset, \\
      (a^{-1})^{log c / log a} = c^{-1}, \quad \text{where } a \in \Rset, a > 0.
    \end{cases}
  \f]
*/
Expr
simplify_logarithm_in_expanded_ex(const Expr& e) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten += simplify_logarithm_in_expanded_ex(e.op(i));
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten *= simplify_logarithm_in_expanded_ex(e.op(i));
  }
  else if (e.is_a_power())
    // Apply the third and fourth properties.
    return prepare_change_base_logarithm(e.arg(0), e.arg(1));
  else if (e.is_a_function()) {
    if (e.is_the_log_function()) {
      const Expr& arg_log = e.arg(0);
      // Apply the second property.
      if (arg_log.is_a_power()) {
	const Expr& base = arg_log.arg(0);
	const Expr& exponent = arg_log.arg(1);
	Number num_base;
	if (base.is_a_number(num_base)) {
	  // Factorize the base of the argument of the logarithm.
	  std::vector<Number> bases;
	  std::vector<int> exponents;
	  partial_factor(num_base, bases, exponents);
	  if (exponents.size() == 1)
	    return exponent * exponents[0] * log(bases[0]);
	  else {
	    Number new_base = bases[0];
	    for (unsigned i = exponents.size(); i-- > 1; )
	      if (exponents[i] != exponents[0])
		return e;
	      else
		new_base *= bases[i];
	    return exponent * exponents[0] * log(new_base);
	  }
	}
	else
	  return exponent * log(base);
	// Apply the first property.
	if (base.is_the_exp_function() && base.arg(0) == 1)
	  return exponent;
      }
      Number arg_log_num;
      if (arg_log.is_a_number(arg_log_num)
	  && arg_log_num.is_positive_integer()) {
	// Factorize the base of the argument of the logarithm.
	std::vector<Number> bases;
	std::vector<int> exponents;
	partial_factor(arg_log_num, bases, exponents);
	if (exponents.size() == 1)
	  return exponents[0] * log(bases[0]);
	else {
	  Number new_base = bases[0];
	  for (unsigned i = exponents.size(); i-- > 1; )
	    if (exponents[i] != exponents[0])
	      return e;
	    else
	      new_base *= bases[i];
	  return exponents[0] * log(new_base);
	}
      }
      else
	return e;
    }
    else if (e.nops() == 1)
      return apply(e.functor(), simplify_logarithm_in_expanded_ex(e.arg(0)));
    else {
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned i = 0; i < num_argument; ++i)
	argument[i] = simplify_logarithm_in_expanded_ex(e.arg(i));
      return apply(e.functor(), argument);
    }
  }
  else 
    e_rewritten = e;
  return e_rewritten;
}

} // anonymous namespace



PURRS::Expr
PURRS::simplify_ex_for_input(const Expr& e, bool input) {
  return simplify_expanded_ex_for_input(e.expand(), input);
}

PURRS::Expr
PURRS::simplify_ex_for_output(const Expr& e, bool input) {
  return simplify_expanded_ex_for_output(e.expand(), input);
}

/*!
  Using the function <CODE>numerator_denomominator_purrs()</CODE> we
  are able to simplify better rational expression.
*/
PURRS::Expr
PURRS::simplify_numer_denom(const Expr& e) {
#if 1
  // Since we need both numerator and denominator, to call 'numer_denom'
  // is faster than to use 'numer()' and 'denom()' separately.
  Expr numer_e;
  Expr denom_e;
  e.numerator_denominator(numer_e, denom_e);
  Expr num = numer_e.expand();
  Expr den = denom_e.expand();
  return num * pwr(den, -1);
#else
  Expr numer_e;
  Expr denom_e;
  numerator_denominator_purrs(e, numer_e, denom_e);
  Expr num = numer_e.expand();
  Expr den = denom_e.expand();
  return num * pwr(den, -1);
  //  return transform_in_single_fraction(e);
#endif
}

/*!
  Let \p e be an expression containing ratios of factorials or ratios
  of exponentials.
  This function looks for possible instances of factorials of type
  \f$ (a n + b)! \f$, with \f$ a \in \Nset \setminus \{0\} \f$ and
  \f$ b \in \Zset \f$, and for possible instances of exponentials
  both in the numerator and denominator of \p e, in order to collect
  common factors, which are then erased.
  We remark that, for this type of simplifications, we do not call
  the function <CODE>simplify_expanded_ex_for_output()</CODE> because it
  would expand the expression while we want to maintain the product of factors.
  Returns a <CODE>Expr</CODE> that contains the modified expression \p e.
*/
PURRS::Expr
PURRS::simplify_factorials_and_exponentials(const Expr& e) {
  Expr e_numerator;
  Expr e_denominator;
  numerator_denominator_purrs(e, e_numerator, e_denominator);
  e_numerator = rewrite_factorials_and_exponentials(e_numerator);
  e_denominator = rewrite_factorials_and_exponentials(e_denominator);
  return e_numerator / e_denominator;
}

PURRS::Expr
PURRS::simplify_logarithm(const Expr& e) {
  return simplify_logarithm_in_expanded_ex(e.expand());
}

/*!
  Executes consecutively all simplifications described in the comment
  of the functions <CODE>simplify_numer_denom()</CODE>,
  <CODE>simplify_factorials_and_exponentials()</CODE> and
  <CODE>simplify_expanded_on_output_ex()</CODE>.
*/
PURRS::Expr
PURRS::simplify_all(const Expr& e) {
  Expr e_rewritten = e;
  e_rewritten = simplify_numer_denom(e_rewritten);
  e_rewritten = simplify_factorials_and_exponentials(e_rewritten);
  e_rewritten = simplify_ex_for_output(e_rewritten, false);
  e_rewritten = simplify_numer_denom(e_rewritten).expand();
  return e_rewritten;
}
