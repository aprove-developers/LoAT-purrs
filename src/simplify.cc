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

#include "sum_poly.hh"
#include "util.hh"
#include "ep_decomp.hh"
#include "gosper.hh"
#include "Expr.defs.hh"
#include "Recurrence.defs.hh"
#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

#define Napier exp(Expr(1))

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
//! Returns <CODE>true</CODE> if the rule `E2' of the rules's set
//! <EM>Expand</EM> is applicable and in this case \p e will contain the
//! new exponent for the expression; returns <CODE>false</CODE> otherwise.
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
    unsigned int num_factors = e.nops();
    unsigned int i;
    for (i = 0; i < num_factors; ++i)
      if (e.op(i) == Recurrence::n)
	break;
    if (i < num_factors) {
      // Found an occurrence of the symbol `n'.
      Expr r = 1;
      for (unsigned int j = 0; j < num_factors; ++j)
	if (i != j)
	  r *= e.op(j);
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
      for (unsigned int i = exponent.nops(); i-- > 0; )
        split_exponent(exponent.op(i), numeric_exponent, not_numeric_exponent);
    else
      split_exponent(exponent, numeric_exponent, not_numeric_exponent);
    // Note that doing `base = base.arg(0)' here would not work.
    // In fact, base could have only one reference and, in this case,
    // it would be deallocated after taking a reference to its `arg(0)'
    // subexpression but before completing the assignment.
    Expr arg_0 = base.arg(0);
    base = arg_0;
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
  for (unsigned int i = base.nops(); i-- > 0; ) {
    vect_num_exp[i] = num_exponent;
    vect_not_num_exp[i] = not_num_exponent;
  }
  for (unsigned int i = base.nops(); i-- > 0; ) {
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
  for (unsigned int i = base.nops(); i-- > 0; )
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
  assert(e.is_expanded());
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
  // The second condition of the `if' is necessary in order to do not
  // consider the case `base = -a' with `a' not number because it is
  // considered like a product: `base = - 1 * a').
  if (base.is_a_mul()
      && !(base.nops() == 2 && (base.op(0) == -1 || base.op(1) == -1)))
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
  for (unsigned int i = e.nops(); i-- > 0; ) {
    const Expr& factor = e.op(i);
    // Special case: e = -(...).
    if (factor == -1) {
      bases.push_back(factor);
      exponents.push_back(1);
    }
    else if (factor.is_a_power()) {
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
  When \p collect_only_n is <CODE>true</CODE> then the function
  collect powers with the same exponent only when it is equal
  to \p Recurrence::n.
*/
Expr
collect_same_exponents(const Expr& e, bool collect_only_n = false) {
  assert(e.is_a_mul());
  D_MSGVAR("*** inizio collect_same_exponents: ", e);
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
  for (unsigned int i = exponents.size(); i-- > 0; ) {
    for (unsigned int j = i; j-- > 0; ) {
      if (exponents[j] == exponents[i]) {
	if (!collect_only_n
	    || (collect_only_n && exponents[j] == Recurrence::n)) {
	  // We have found two equal exponents. We store in the i-th
	  // position of the vectors `bases' and `exponents' the new power
	  // with base equal to the product of the base in the i-th and in
	  // the j-th position and exponent unchanged.
	  bases[i] *= bases[j];
	  exponents[j] = 0;
	}
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
  D_VAR(e_rewritten);

  // Now adds to `e_rewritten' the factors of `e' not considered in the
  // previous cycle, i.e., the factors which are not powers and the factor
  // `-1' (if any).
  for (unsigned int i = e.nops(); i-- > 0; ) {
    const Expr& factor = e.op(i);
    if (!factor.is_a_power() && factor != -1)
      e_rewritten *= factor;
  }
  D_MSGVAR("*** fine collect_same_exponents: ", e_rewritten);
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
  D_MSGVAR("*** inizio collect_same_base: ", e);
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
  for (unsigned int i = exponents.size(); i-- > 0; ) {
    for (unsigned int j = i; j-- > 0; )
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
  D_VAR(e_rewritten);

  // Now adds to `e_rewritten' the factors of `e' not considered in the
  // previous cycle, i.e., the factors which are not powers applying,
  // when possible, the rule `C1'.
  for (unsigned int i = e.nops(); i-- > 0; ) {
    const Expr& factor_e = e.op(i);
    ///////
    if (!factor_e.is_a_power() && factor_e != -1) {
      // We must consider those factors that are not powers but are equal
      // to some base of the vector `bases', in order to apply the rule `C1',
      // i.e., `a * a^e = a^(e + 1)'. In this case, we add `1' to the
      // correspondenting exponent.
      bool to_sum = false;
      Expr old_factor;
      Expr new_factor;
      for (unsigned int j = bases.size(); j-- > 0; ) {
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
  D_MSGVAR("*** fine collect_same_base: ", e_rewritten);
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
  for (unsigned int i = bases.size(); i-- > 0; ) {
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
  unsigned int num_size = num_bases.size();
  Expr reduced_num = to_std_form(k, num_bases, num_exponents);
  // Here <CODE>to_std_form</CODE> is called with a negative value of k 
  // because we are dealing with the denominator of r.
  partial_factor(den, den_bases, den_exponents);
  unsigned int den_size = den_bases.size();
  Expr reduced_den = to_std_form(-k, den_bases, den_exponents);
  
  // Try one last simplification: if all exponents have a common factor 
  // with the root index, remove it.
  int gc = k;
  for (unsigned int i = 0; i < num_size && gc > 1; ++i)
    gc = gcd(gc, num_exponents[i]);
  for (unsigned int i = 0; i < den_size && gc > 1; ++i)
    gc = gcd(gc, den_exponents[i]);
  
  if (gc > 1) {
    k /= gc;
    for (unsigned int i = 0; i < num_size; ++i)
      num_exponents[i] /= gc;
    for (unsigned int i = 0; i < den_size; ++i)
      den_exponents[i] /= gc;
  }
  
  // The object `irr_part' is surely numeric but we have to use an expression
  // because otherwise the number, since it is irrational, is rounded.
  Expr irr_part = 1;
  for (unsigned int i = 0; i < num_size; ++i)
    irr_part *= pwr(num_bases[i], num_exponents[i]);
  for (unsigned int i = 0; i < den_size; ++i)
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
  base_1 = exact_pwr(base_1, k1_num);
  base_2 = exact_pwr(base_2, k2_num);
  Number g = gcd(k1_den, k2_den);
  Number k = k1_den * k2_den / g;
  Number b1 = exact_pwr(base_1, k2_den / g);
  Number b2 = exact_pwr(base_2, k1_den / g);
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
  for (unsigned int i = e.nops(); i-- > 0; ) {
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
	  for (unsigned int j = 2; j-- > 0; ) {
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
  D_MSGVAR("inizio manip ", e);
  Expr e_rewritten = 1;
  // Simplifies each factor that is a `power'.
  for (unsigned int i = e.nops(); i-- > 0; ) {
    const Expr& factor_e = e.op(i);
    if (factor_e.is_a_power()) {
      const Expr& base
	= simplify_expanded_ex_for_output(factor_e.arg(0).expand(), input);
      const Expr& exp
	= simplify_expanded_ex_for_output(factor_e.arg(1).expand(), input);
      Number num_base;
      Number num_exp;
      if (base.is_a_number(num_base) && exp.is_a_number(num_exp))
	e_rewritten *= red_prod(num_base, num_exp, 1, 1);
      else {
	const Expr& tmp = pwr(base, exp).expand();
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
    for (unsigned int i = e_rewritten.nops(); i-- > 0; ) {
      const Expr& factor_e_rewritten = e_rewritten.op(i);
      if (factor_e_rewritten.is_a_function())
	if (factor_e_rewritten.nops() == 1)
	  factor_function
	    *= apply(factor_e_rewritten.functor(),
		     simplify_expanded_ex_for_output(factor_e_rewritten.arg(0)
						     .expand(), input));
	else {
	  unsigned int num_argument = factor_e_rewritten.nops();
	  std::vector<Expr> argument(num_argument);
	  for (unsigned int i = 0; i < num_argument; ++i)
	    argument[i]
	      = simplify_expanded_ex_for_output(factor_e_rewritten.arg(i)
						.expand(), input);
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
			  simplify_expanded_ex_for_output(e_rewritten.arg(0)
							  .expand(), input));
    else {
      unsigned int num_argument = e_rewritten.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
	argument[i] = simplify_expanded_ex_for_output(e_rewritten.arg(i)
						      .expand(), input);
      e_rewritten = apply(e_rewritten.functor(), argument);
    }
  }
  D_MSGVAR("e_rewritten dopo function: ", e_rewritten);
  // Special case: the exponential `exp' is a `function' but it has
  // the same properties of the powers.
  if (e_rewritten.is_a_mul()) {
    Expr argument = 0;
    Expr rem = 1;
    for (unsigned int i = e_rewritten.nops(); i-- > 0; ) {
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
    for (unsigned int i = e_rewritten.nops(); i-- > 0; ) { 
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
  D_MSGVAR("e_rewritten dopo `reduce_product': ", e_rewritten);
  // Simplifies eventual powers with same base or same exponents applying
  // the rule of the term rewriting system \f$ \mathfrak{R}_o \f$.
  if (e_rewritten.is_a_add()) {
    Expr terms = 0;
    for (unsigned int i = e_rewritten.nops(); i-- > 0; ) { 
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
  D_MSGVAR("fine manip ", e_rewritten);
  return e_rewritten;
}

/*!
  Crosses the tree of the expanded expression \p e recursively to find
  subexpressions to which we apply the rules of the terms rewriting
  system \f$ \mathfrak{R}_i \f$ or to which we can compute symbolic sums.
  More exactly here the rules of the set <EM>Expand</EM> are implemented
  because the rules of the set <EM>Automatic</EM> are automatically executed.
  We remark that the rule \f$ \textbf{E4} \f$ is automatically executed when
  the exponent is integer, while the rules \f$ \textbf{E3} \f$ and
  \f$ \textbf{E6} \f$ are executed by the method <CODE>expand()</CODE>
  (\f$ \textbf{E3} \f$ only partially because for instance
  \f$ expand(3^(4*x+2*a)) = 3^(2*a)*3^(4*x) \f$): hence we here consider
  only powers.
  \p input is always <CODE>true</CODE> and this means that \p n is a
  special symbol, i. e., in the simplifications it is always collected
  with respect to the others parameters.  The function returns a
  new expression containing the modified expression \p e.
*/
Expr
simplify_expanded_ex_for_input(const Expr& e, bool input) {
  assert(e.is_expanded());
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += simplify_expanded_ex_for_input(e.op(i).expand(), input);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten *= simplify_expanded_ex_for_input(e.op(i).expand(), input);
    // In the case of expressions for input we are interesting to collect
    // powers with the same exponents when this is equal to `Recurrence::n'.
    if (e_rewritten.is_a_mul())
      e_rewritten = collect_same_exponents(e_rewritten, true);
  }
  else if (e.is_a_power())
    return simplify_powers(e, input);
  else if (e.is_a_function()) {
    if (e.nops() == 1)
      return apply(e.functor(),
		   simplify_expanded_ex_for_input(e.arg(0).expand(), input));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
	argument[i] = simplify_expanded_ex_for_input(e.arg(i).expand(), input);
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
  assert(e.is_expanded());
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += simplify_expanded_ex_for_output(e.op(i).expand(), input);
  }
  else if (e.is_a_mul())
    // We can not call `simplify_expanded_ex_for_output()' on every factor
    // because otherwise it is not possible to transform products.
    return manip_factor(e, input);
  else if (e.is_a_power()) {
    // Is necessary to call `red_prod()' in order to rewrite irrational
    // numbers, i.e., powers with base and exponent both numerics.
    const Expr& base = simplify_expanded_ex_for_output(e.arg(0).expand(),
						       input);
    const Expr& exp = simplify_expanded_ex_for_output(e.arg(1).expand(),
						      input);
    Number num_base;
    Number num_exp;
    if (base.is_a_number(num_base) && exp.is_a_number(num_exp))
      e_rewritten = red_prod(num_base, num_exp, 1, 1);
    else {
      const Expr& tmp = pwr(base, exp).expand();
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
      return apply(e.functor(),
		   simplify_expanded_ex_for_output(e.arg(0).expand(), input));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
	argument[i] = simplify_expanded_ex_for_output(e.arg(i).expand(),
						      input);
      return apply(e.functor(), argument);
    }
  }
  else
    return e;
  return e_rewritten;
}

/*!
  Given the factorial expression \f$ (a + b)! \f$, where \f$ a \f$ is
  a non-numeric expression and \f$ b \f$ an integer, this function 
  returns a new expression that contains explicitly \f$ a! \f$.
  We use the rewrite rule explained in the comment for
  <CODE>rewrite_factorials()</CODE>
*/
Expr
decompose_factorial(const Expr& e) {
  assert(e.is_the_factorial_function());
  const Expr& argument = e.arg(0);
  if (argument.is_a_add()) {
    Number num = 0;
    Expr new_arg_fact = 0;
    for (unsigned int i = argument.nops(); i-- > 0; ) {
      Number tmp_num;
      if (argument.op(i).is_a_number(tmp_num) && tmp_num.is_integer())
	// We are sure that `num' is again `0' because automatically
	// two integer numbers would have been added.
	num = tmp_num;
      else
	new_arg_fact += argument.op(i);
    }
    if (num.is_positive()) {
      Expr new_factors = 1;
      for (Number i = 1; i <= num; ++i)
	new_factors *= new_arg_fact + i;
      return factorial(new_arg_fact) * new_factors;
    }
    else {
      Expr new_factors = 1;
      for (Number i = 0; i < -num; ++i)
	new_factors *= pwr(new_arg_fact - i, -1);
      return factorial(new_arg_fact) * new_factors;
    }
  }
  return e;
}

/*!
  Applies the following rewrite rules:
  - \f[
      a^{b n + c} = a^{b n} \cdot a^c;
    \f]
  - let \f$ a \f$ be a non-numeric expression and let \f$ b \in \Zset \f$:
    \f[
      \begin{cases}
        (a + b)!
        =
        a! \cdot (a + 1) \cdots (a + b),
          \quad \text{if } b \in \Nset \setminus \{0\}; \\
        a!,
          \quad \text{if } b = 0; \\
        \dfrac{a!}{a \cdot (a - 1) \cdots (a + b + 1))},
          \quad \text{if } b \in \Zset \setminus \Nset.
      \end{cases}
    \f]
*/
Expr
rewrite_factorials_and_exponentials(const Expr& e) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += rewrite_factorials_and_exponentials(e.op(i));
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
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
	unsigned int num_argument = e.nops();
	std::vector<Expr> argument(num_argument);
	for (unsigned int i = 0; i < num_argument; ++i)
	  argument[i] = rewrite_factorials_and_exponentials(e.arg(i));
	return apply(e.functor(), argument);
      }
  else
    return e;
  return e_rewritten;
}

/*!
  Applies the following rewrite rules:
  \f[
    \binom{n}{k}
    \begin{cases}
      \frac{n!}{(n-k)! k!},
        \quad \text{if } n \text{ and } k \text{ are not both numerics or }
	n \in \Zset, n \geq 0 \text{ or } k \in \Zset, k \geq 0; \\
      (-1)^k \frac{(k-n-1)!}{(-n-1)! k!},
        \quad \text{if } n \in \Zset, n < 0; \\
      \binom{n}{k},
      \quad \text{otherwise}.
    \end{cases}
  \f]
  Note that it is not possible that both the arguments of the function
  \f$ \binom() \f$ are numerics.
*/
Expr
rewrite_binomials(const Expr& e) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += rewrite_binomials(e.op(i));
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten *= rewrite_binomials(e.op(i));
  }
  else if (e.is_a_power())
    return pwr(rewrite_binomials(e.arg(0)), rewrite_binomials(e.arg(1)));
  else if (e.is_a_function())
    if (e.is_the_binom_function()) {
      const Expr& m = e.arg(0);
      const Expr& k = e.arg(1);
      assert(!(m.is_a_number() && k.is_a_number()));
      Number num_m;
      Number num_k;
      if (!m.is_a_number() && !k.is_a_number())
	return factorial(m) / (factorial(m-k)*factorial(k));
      else if (m.is_a_number(num_m) && num_m.is_integer())
	if (num_m >= 0)
	  return factorial(num_m) / (factorial(m-k)*factorial(k));
	else
	  return pwr(-1, k)*factorial(k-num_m-1)
	    / (factorial(-num_m-1)*factorial(k));
      else if (k.is_a_number(num_k) && num_k.is_nonnegative_integer())
	return factorial(m) / (factorial(m-num_k)*factorial(num_k));
      else
	return e;
    }
    else
      if (e.nops() == 1)
	return apply(e.functor(),
		     rewrite_binomials(e.arg(0)));
      else {
	unsigned int num_argument = e.nops();
	std::vector<Expr> argument(num_argument);
	for (unsigned int i = 0; i < num_argument; ++i)
	  argument[i] = rewrite_binomials(e.arg(i));
	return apply(e.functor(), argument);
      }
  else
    return e;
  return e_rewritten;
}

Expr
factorize_base_arg_log(const Number& base, const Expr& exponent = 1) {
  assert(base.is_integer());
  std::vector<Number> bases;
  std::vector<int> exponents;
  partial_factor(base, bases, exponents);
  if (exponents.size() == 1)
    return exponent * exponents[0] * log(bases[0]);
  else {
    Number new_base = bases[0];
    for (unsigned int i = exponents.size(); i-- > 1; )
      if (exponents[i] != exponents[0])
	return log(pwr(base, exponent));
      else
	new_base *= bases[i];
    return exponent * exponents[0] * log(new_base);
  }
}

Expr
logarithm_of_product(const Expr& arg_log) {
  assert(arg_log.is_a_mul());
  Expr sum_log = 0;
  for (unsigned int i = arg_log.nops(); i-- > 0; ) {
    const Expr& factor = arg_log.op(i);
    if (factor.is_a_power() || factor.is_a_number())
      sum_log += simplify_logarithm_in_expanded_ex(log(factor));
    else
      sum_log += log(factor);
  }
  return sum_log;
}

/*!
  Applies the following logarithm's property:
  \f[
    \begin{cases}
      log(exp(1)^a) = a, \\
      log(a^b) = b log(a), \\
      log(a * b) = log(a) + log(b), \\
      log(a / b) = log(a) - log(b).
    \end{cases}
  \f]
  The last two properties are applied also if \p all_properties is
  <CODE>true</CODE>.
*/
Expr
apply_elementary_prop(const Expr& e, bool all_properties = true) {
  assert(e.is_the_log_function());
  Expr arg_log;
  if (e.arg(0).is_a_power())
    arg_log = simplify_powers(e.arg(0), true);
  else
    arg_log = e.arg(0);
  // Apply the second property.
  Number arg_log_num;
  if (arg_log.is_a_power()) {
    const Expr& base = arg_log.arg(0);
    const Expr& exponent = arg_log.arg(1);
    Number num_base;
    if (base.is_a_number(num_base)) {
      if (num_base.is_integer())
	// Factorize the base of the argument of the logarithm.
	return factorize_base_arg_log(num_base, exponent);
      else {
	assert(num_base.is_rational());
	return factorize_base_arg_log(num_base.numerator(), exponent)
	  - factorize_base_arg_log(num_base.denominator(), exponent);
      }
    }
    else
      return exponent * log(base);
    // Apply the first property.
    if (base.is_the_exp_function() && base.arg(0) == 1)
      return exponent;
  }
  else if (arg_log.is_a_number(arg_log_num)) {
    if (arg_log_num.is_positive_integer())
      // Factorize the base of the argument of the logarithm and apply
      // the second property.
      return factorize_base_arg_log(arg_log_num);
    if (all_properties && arg_log_num.denominator() != 1) {
      Number numer = arg_log_num.numerator();
      Number denom = arg_log_num.denominator();
      assert(denom.is_positive_integer());
      // Apply the third and the fourth properties.
      if (numer.is_positive_integer())
	return factorize_base_arg_log(numer)
	  - factorize_base_arg_log(denom);
      else
	return numer - factorize_base_arg_log(denom);
    }
  }
  // Apply the third and the fourth properties.
  if (all_properties && arg_log.is_a_mul())
    return logarithm_of_product(arg_log);
  return log(arg_log);
}

Expr
simpl_power_in_logarithm(const Expr& base, const Expr exponent,
			 const Number& new_base,
			 const Number& new_exponent = 1) {
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

/*!
  Returns <CODE>true</CODE> if exists an integer \p b such that
  \f$ a^b = c \f$; returns <CODE>false</CODE> otherwise.
*/
bool
is_perfect_power(const Number& c, const Number& a, Number& b) {
  assert(a.is_positive() && c.is_positive());
  // Special case.
  if (a == 1 && c != 1)
    return false;
  for (int i = 1; ; ++i)
    if (exact_pwr(a, i) == c) {
      b = i;
      return true;
    }
    else if (a > 1 && exact_pwr(a, i) > c)
      break;
    else if (a < 1 && exact_pwr(a, i) < c)
      break;
  for (int i = -1; ; --i)
    if (exact_pwr(a, i) == c) {
      b = i;
      return true;
    }
    else if (a > 1 && exact_pwr(a, i) < c)
      break;
    else if (a < 1 && exact_pwr(a, i) > c)
      break;
  return false;
}

/*!
  Applies the following logarithm's property:
  \f[
    \begin{cases}
      a^{c log b} = b^{c log a},
        \quad \text{if } b { contains the special symbol Recurrence::n or if }
	a \in \Rset_+ \wedge b \text{ not a number}; \\
      (a^b)^{log c / log a} = c^b, \quad \text{where } a \in \Rset, a > 0
        \quad \text{and } b \in \Nset; \\
      (a^{-1})^{log c / log a} = c^{-1},
        \quad \text{where } a \in \Rset, a > 0.
    \end{cases}
  \f]
*/
Expr
prepare_simpl_power_in_logarithm(const Expr& base, const Expr& exponent) {
  // Apply the properties `log(exp(1)^a) = a' and
  // `log(a^b) = b log(a)' to the exponent.
  Expr exponent_simpl = 1;
  if (exponent.is_a_mul())
    for (unsigned int i = exponent.nops(); i-- > 0; ) {
      const Expr& factor = exponent.op(i);
      if (factor.is_the_log_function())
	exponent_simpl *= apply_elementary_prop(factor, false);
      else if (factor.is_a_power() && factor.arg(0).is_the_log_function())
	exponent_simpl *= pwr(apply_elementary_prop(factor.arg(0), false),
			      factor.arg(1));
      else
	exponent_simpl *= factor;
    }
  else
    if (exponent.is_the_log_function())
      exponent_simpl = apply_elementary_prop(exponent, false);
    else
      exponent_simpl = exponent;

  // Apply the first property.
  Expr base_simpl = base;
  if (exponent_simpl.is_a_mul()) {
    Expr arg_factor_log = 1;
    Expr rem = 1;
    for (unsigned int i = exponent_simpl.nops(); i-- > 0; ) {
      const Expr& factor = exponent_simpl.op(i);
      if (factor.is_the_log_function()
	  && (factor.arg(0).has(Recurrence::n)
	      || (base.is_a_number() && !factor.arg(0).is_a_number())))
	if (arg_factor_log == 1)
	  arg_factor_log *= factor.arg(0);
	else
	  // In the exponent thare is more than one factor of the shape
	  // `log(b)', where `b' contains the special symbol `Recurrence::n'.
	  break;
      else
	rem *= factor;
    }
    if (arg_factor_log != 1) {
      base_simpl = arg_factor_log;
      exponent_simpl = rem * log(base);
    }
  }

  // Apply the second and the third properties.
  Number base_num;
  if (base_simpl.is_a_number(base_num) && base_num.is_positive()
      && exponent_simpl.is_a_mul()) {
    // `arg_log_num' will contain the argument of the eventual
    // factor `log(.)';
    // `arg_log_den' will contain the argument of the eventual
    // factor `1/log(.)';
    // `rem' will contain all factors different from `log(.)' and `1/log(.)'.
    // If in the exponent there is more than one factor of the shape
    // `log(.)' or more than one factor of the shape `1/log(.)', the
    // simplification's process is stopped.
    Expr arg_log_num = 1;
    Number arg_log_den = 1;
    Expr rem = 1;
    bool stop_simplification = false;
    for (unsigned int i = exponent_simpl.nops(); i-- > 0; ) {
      const Expr& factor = exponent_simpl.op(i);
      if (factor.is_the_log_function())
	if (arg_log_num == 1)
	  arg_log_num = factor.arg(0);
	else {
	  // In the exponent thare is more than one factor of the shape
	  // `log(c)'.
	  stop_simplification = true;
	  break;
	}
      else if (factor.is_a_power()
	       && factor.arg(0).is_the_log_function()
	       && factor.arg(0).arg(0).is_a_number()
	       && factor.arg(1) == -1)
	if (arg_log_den == 1)
	  arg_log_den = factor.arg(0).arg(0).ex_to_number();
	else {
	  // In the exponent thare is more than one factor of the shape
	  // `1/log(a)'.
	  stop_simplification = true;
	  break;
	}
      else
	rem *= factor;
    }
    Number exp;
    if (!stop_simplification && arg_log_num != 1 && arg_log_den != 1
	&& is_perfect_power(base_num, arg_log_den, exp))
      return pwr(arg_log_num, exp * rem);
  }
  else if (base_simpl == Napier) {
    if (exponent_simpl.is_the_log_function())
      return exponent_simpl.arg(0);
    if (exponent_simpl.is_a_mul()) {
      Expr log_factors_exp = 1;
      Expr rem_factors_exp = 1;
      for (unsigned int i = exponent_simpl.nops(); i-- > 0; ) {
	const Expr& factor = exponent_simpl.op(i);
	if ((factor.is_a_power() && factor.arg(0).is_the_log_function()
	     && factor.arg(1) == -1)
	    || factor.is_the_log_function())
	  log_factors_exp *= factor;
	else
	  rem_factors_exp *= factor;
      }
      if (log_factors_exp.nops() == 1)
	return pwr(log_factors_exp.arg(0), rem_factors_exp);
    }
  }
  return pwr(simplify_logarithm_in_expanded_ex(base_simpl),
	     simplify_logarithm_in_expanded_ex(exponent_simpl));
}

/*!
  This function applies the following rewrite rule:
  \f[
    \frac{\log(a/b)}{\log(b/a)} = -1.
  \f]
*/
Expr
simpl_quotient_of_logs(const Expr& e) {
  assert(e.is_a_mul());
  // Divide the factors of \p e in 3 different expressions:
  // `log_factors' will contain the factors of the form `log(a)';
  // `inv_log_factors' will contain the factors of the form `1 / log(a)';
  // all the other factors will be contained in `e_simplified'.
  Expr e_simplified = 1;
  Expr log_factors = 1;
  Expr inv_log_factors = 1;
  for (unsigned int i = e.nops(); i-- > 0; ) {
    const Expr& factor = e.op(i);
    if (factor.is_the_log_function())
      log_factors *= factor;
    else if (factor.is_a_power()
	     && factor.arg(0).is_the_log_function() && factor.arg(1) == -1)
      inv_log_factors *= factor;
    else
      e_simplified *= factor;
  }
  // No simplification to apply.
  if (log_factors == 1 || inv_log_factors == 1)
    return e_simplified * log_factors * inv_log_factors;

  // There is more than one factor of the form `log(a)'.
  if (log_factors.is_a_mul())
    for (unsigned int i = log_factors.nops(); i-- > 0; ) {
      const Expr& arg_log_1 = log_factors.op(i).arg(0);
      // There is more than one factor of the form `1 / log(a)'.
      if (inv_log_factors.is_a_mul())
	for (unsigned int j = inv_log_factors.nops(); j-- > 0; ) {
	  const Expr& arg_log_2 = inv_log_factors.op(j).arg(0).arg(0);
	  if (arg_log_1.numerator() == arg_log_2.denominator()
	      && arg_log_2.numerator() == arg_log_1.denominator()) {
	    e_simplified *= -1;
	    inv_log_factors /= inv_log_factors.op(j);
	    break;
	  }
	  else {
	    e_simplified *= log_factors.op(i) * inv_log_factors.op(j);
	    inv_log_factors /= inv_log_factors.op(j);
	    break;
	  }
	}
      // There is only 1 factor of the form `1 / log(a)'.
      else {
	if (inv_log_factors == 1)
	  e_simplified *= log_factors.op(i);
	else {
	  const Expr& arg_log_2 = inv_log_factors.arg(0).arg(0);
	  if (arg_log_1.numerator() == arg_log_2.denominator()
	      && arg_log_2.numerator() == arg_log_1.denominator()) {
	    e_simplified *= -1;
	    inv_log_factors = 1;
	  }
	  else {
	    e_simplified *= log_factors.op(i) * inv_log_factors;
	    inv_log_factors = 1;
	  }
	}
      }
    }
  // There is only 1 factor of the form `log(a)'.
  else {
    const Expr& arg_log_1 = log_factors.arg(0);
    // There is more than one factor of the form `1 / log(a)'.
    if (inv_log_factors.is_a_mul())
      for (unsigned int j = inv_log_factors.nops(); j-- > 0; ) {
	const Expr& arg_log_2 = inv_log_factors.op(j).arg(0).arg(0);
	if (arg_log_1.numerator() == arg_log_2.denominator()
	    && arg_log_2.numerator() == arg_log_1.denominator())
	  return e_simplified * -1 * inv_log_factors / inv_log_factors.op(j);
	else
	  return e_simplified * log_factors * inv_log_factors
	    / inv_log_factors.op(j);
      }
    // There is only 1 factor of the form `1 / log(a)'.
    else {
      const Expr& arg_log_2 = inv_log_factors.arg(0).arg(0);
      if (arg_log_1.numerator() == arg_log_2.denominator()
	  && arg_log_2.numerator() == arg_log_1.denominator())
	return e_simplified * -1;
      else
	return e_simplified * log_factors * inv_log_factors;
    }
  }
  return e_simplified * inv_log_factors;
}

/*!
  Applies the following logarithm's property:
  \f[
    \begin{cases}
      log(exp(1)^a) = a, \\
      log(a^b) = b log(a), \\
      log(a * b) = log(a) + log(b), \\
      log(a / b) = log(a) - log(b), \\
      a^{c log b} = b^{c log a},
        \quad \text{if } b { contains the special symbol Recurrence::n}, \\
      (a^b)^{log c / log a} = c^b, \quad \text{where } a \in \Rset, a > 0
        \quad \text{and } b \in \Nset, \\
      (a^{-1})^{log c / log a} = c^{-1},
        \quad \text{where } a \in \Rset, a > 0.
    \end{cases}
  \f]
*/
Expr
simplify_logarithm_in_expanded_ex(const Expr& e) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_rewritten += simplify_logarithm_in_expanded_ex(e.op(i));
  }
  else if (e.is_a_mul()) {
    const Expr& tmp = simpl_quotient_of_logs(e);
    if (tmp.is_a_mul()) {
      e_rewritten = 1;
      for (unsigned int i = tmp.nops(); i-- > 0; )
	e_rewritten *= simplify_logarithm_in_expanded_ex(tmp.op(i));
    }
    else
      return simplify_logarithm_in_expanded_ex(tmp);
  }
  else if (e.is_a_power())
    // Apply the last three properties: note that is important that
    // these properties are applied before than the third and the fourth
    // properties.
    return prepare_simpl_power_in_logarithm(e.arg(0), e.arg(1));
  else if (e.is_a_function()) {
    // Apply the firsts four properties.
    if (e.is_the_log_function())
      return apply_elementary_prop(e);
    else if (e.nops() == 1)
      return apply(e.functor(), simplify_logarithm_in_expanded_ex(e.arg(0)));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
	argument[i] = simplify_logarithm_in_expanded_ex(e.arg(i));
      return apply(e.functor(), argument);
    }
  }
  else 
    e_rewritten = e;
  return e_rewritten;
}

/*!
  Let \f$ e(x) \f$ be the expression in \p x contained in \p e.
  This function computes two expressions \f$ e_1 \f$ and \f$ e_2 \f$
  such that \f$ e = e_1 \cdot e_2 \f$: \f$ e_1 \f$ contains all factors of
  \f$ e \f$ that do not depend from the symbol \p x; \f$ e_2 \f$
  contains all factors of \f$ e \f$ that depend from the symbol \p x.
*/
void
get_out_factors_from_argument(const Expr& e, const Expr& x,
			      Expr& in, Expr& out) {
  if (e.is_a_mul())
    for (unsigned int i = e.nops(); i-- > 0; ) {
      const Expr& factor = e.op(i);
      if (factor.has(x))
	in *= factor;
      else
	out *= factor;
    }
  else
    if (e.has(x))
      in *= e;
    else
      out *= e;
}

Expr
split_sum(const Expr& e) {
  Expr e_rewritten;
  const Expr& first_term = e.arg(2).op(0);
  const Expr& second_term = e.arg(2).op(1);
  Number numeric_term;
  Symbol symbolic_term;
  if (first_term.is_a_number() && second_term.is_a_symbol()) {
    numeric_term = first_term.ex_to_number();
    symbolic_term = second_term.ex_to_symbol();
  }
  else if (first_term.is_a_symbol() && second_term.is_a_number()) {
    numeric_term = second_term.ex_to_number();
    symbolic_term = first_term.ex_to_symbol();
  }
  else {
    Expr factors_in = 1;
    Expr factors_out = 1;
    get_out_factors_from_argument(e.arg(3), e.arg(0),
				  factors_in, factors_out);
    if (factors_in == 1)
      return factors_out * (e.arg(2) - e.arg(1) + 1);
    else
      return factors_out * sum(e.arg(0), e.arg(1), e.arg(2), factors_in);
  }
  if (numeric_term.is_integer()) {
    Expr factors_in = 1;
    Expr factors_out = 1;
    get_out_factors_from_argument(e.arg(3), e.arg(0),
				  factors_in, factors_out);
    if (factors_in == 1)
      e_rewritten += factors_out * (symbolic_term - e.arg(1) + 1);
    else
      e_rewritten += factors_out
	* sum(e.arg(0), e.arg(1), Expr(symbolic_term), factors_in);
    if (numeric_term.is_positive_integer())
      for (Number j = 1; j <= numeric_term; ++j)
	e_rewritten += e.arg(3).substitute(e.arg(0), symbolic_term + j);
    else
      for (Number j = numeric_term + 1; j <= 0 ; ++j)
	e_rewritten -= e.arg(3).substitute(e.arg(0), symbolic_term + j);
  }
  return e_rewritten;
}

Expr
compute_sum(const Expr& e) {
  // If `upper' is of the form `m + h' or `m - h' with `m' a symbol
  // and `h' an integer, isolate `m' and `h' respectively in the expressions
  // `symb_part_of_upper' and `numeric_part_of_upper'.
  // If the form of `upper' is different from that aforesaid one, then
  // the system can not to compute the sum.
  Expr upper = e.arg(2);
  Symbol symb_part_of_upper;
  Number numeric_part_of_upper = 0;
  if (upper.is_a_symbol())
    symb_part_of_upper = upper.ex_to_symbol();
  else if (upper.is_a_add() && upper.nops() == 2
	   && upper.op(0).is_a_symbol()
	   && upper.op(1).is_a_number(numeric_part_of_upper)
	   && numeric_part_of_upper.is_integer())
    symb_part_of_upper = upper.op(0).ex_to_symbol();
  else if (upper.is_a_add() && upper.nops() == 2
	   && upper.op(1).is_a_symbol()
	   && upper.op(0).is_a_number(numeric_part_of_upper)
	   && numeric_part_of_upper.is_integer())
    symb_part_of_upper = upper.op(1).ex_to_symbol();
  else
    return e;

  std::vector<Expr> base_of_exps;
  std::vector<Expr> exp_poly_coeff;
  std::vector<Expr> exp_no_poly_coeff;
  exp_poly_decomposition(e.arg(3).expand(), e.arg(0).ex_to_symbol(),
			 base_of_exps, exp_poly_coeff, exp_no_poly_coeff);

  Expr e_rewritten = 0;
  if (vector_not_all_zero(exp_poly_coeff)) {
    for (unsigned int i = exp_poly_coeff.size(); i-- > 0; ) {
      Symbol k;
      Expr coeff_k = exp_poly_coeff[i].substitute(e.arg(0), k);
      Expr solution = sum_poly_times_exponentials(coeff_k, k,
						  symb_part_of_upper,
						  base_of_exps[i]);
      // `sum_poly_times_exponentials' computes the sum until
      // `symb_part_of_upper', whereas we want that the sum start from
      // `symb_part_of_upper + numeric_part_of_upper'. 
      solution
	= solution.substitute(symb_part_of_upper,
			      symb_part_of_upper + numeric_part_of_upper);
      // `sum_poly_times_exponentials' computes the sum from 0, whereas
      // we want that the sum start from the lower bound.
      for (Number h = 0; h < e.arg(1).ex_to_number(); ++h)
	solution -= coeff_k.substitute(k, h) * pwr(base_of_exps[i], h);
      e_rewritten += solution;
    }
  }
  if (vector_not_all_zero(exp_no_poly_coeff)) {
    for (unsigned int i = exp_poly_coeff.size(); i-- > 0; ) {
      Expr gosper_solution;
      if (!gosper_algorithm(symb_part_of_upper,
			    e.arg(3).substitute(e.arg(0), symb_part_of_upper),
			    e.arg(1).ex_to_number(), upper, gosper_solution))
	return e;
      e_rewritten += gosper_solution;
    }
  }
  return e_rewritten;
}

//! \brief
//! If \p simplification is equal to <CODE>REWRITE_UPPER_LIMIT</CODE>
//! rewrite the sum with the upper limit of the form \f$ m + j \f$
//! (\f$ m \f$ is a symbol and \f$ j \in \Zset \f$), so that the
//! upper limit is \f$ m \f$; if \p simplification is equal to
//! <CODE>COMPUTE_SUM</CODE> split the sum in as many sums as
//! the addends of the summand and compute, when possible, symbolic sums.
/*!
 If \p simplification is equal to <CODE>REWRITE_UPPER_LIMIT</CODE> this
 function applies the following property for the function representing
 finite sums:
  \f[
    \begin{cases}
      \sum_{k = a}^b f(k) =
      \sum_{k = a}^m f(k) - f(m) - f(m-1) - \cdots - f(m-j+1),
        \quad \text{if } b = m + j
	\text{with m a symbol and j is a positive integer}; \\
      \sum_{k = a}^b f(k) = \sum_{k = a}^m f(k) + f(m+1) + \cdots + f(m+j),
        \quad \text{if } b = m + j
	\text{with m a symbol and j is a negative integer}.
    \end{cases}
  \f]

  If \p simplification is equal to <CODE>COMPUTE_SUM</CODE> this function
  computes, when possible, symbolic sums using two different techniques:
  if the summand is polynomial, exponential or product of them uses the
  method exposed in <CODE>sum_poly.{hh, cc}</CODE>; otherwise the Gosper's
  algorithm exposed in <CODE>gosper.{hh, cc}</CODE>.
*/
Expr
simplify_sum_in_expanded_ex(const Expr& e,
			    const Sum_Simplification_Kind simplification) {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
    e_rewritten += simplify_sum_in_expanded_ex(e.op(i), simplification);
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
    e_rewritten *= simplify_sum_in_expanded_ex(e.op(i), simplification);
  }
  else if (e.is_a_power())
  return pwr(simplify_sum_in_expanded_ex(e.arg(0), simplification),
	       simplify_sum_in_expanded_ex(e.arg(1), simplification));
  else if (e.is_a_function()) {
    if (e.nops() == 1)
      return apply(e.functor(),
		   simplify_sum_in_expanded_ex(e.arg(0), simplification));
    else {
      if (e.is_the_sum_function()) {
	// Splits the sum in as many sums as the addends of the
	// summand.
	if (simplification == COMPUTE_SUM && e.arg(3).is_a_add()) {
	  e_rewritten = 0;
	  for (unsigned int i = e.arg(3).nops(); i-- > 0; )
	    e_rewritten
	      += simplify_sum_in_expanded_ex(sum(e.arg(0), e.arg(1), e.arg(2),
						 e.arg(3).op(i)),
					     simplification);
	    return e_rewritten;
	}
	// If the summand does not contain functions `x()' with the
	// index of the sum in the argument (weighted-average recurrence)
	// and we are at the first visit to the expression, then the
	// sum has been inserted by the user and we try to compute it.
	if (!e.arg(3).has_x_function(e.arg(0))
	    && simplification == COMPUTE_SUM)
	  return compute_sum(e);
	// `upper' is a sum of two addends: if we are at the point of
	// recurrence's verification rewrite this sum so that the upper bound
	// is only a symbol.
	if (e.arg(2).is_a_add() && e.arg(2).nops() == 2
	    && simplification == REWRITE_UPPER_LIMIT)
	  return split_sum(e);
      }
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i)
      argument[i] = simplify_sum_in_expanded_ex(e.arg(i), simplification);
      return apply(e.functor(), argument);
    }
  }
  else 
    e_rewritten = e;
  return e_rewritten;
}

//! \brief
//! Look for symbolic sums in the expression and, when possible, pack
//! them all into a single symbolic sum. When verifying solutions we
//! can thus simplify sums we cannot compute one by one.
Expr
simplify_collect_sums_in_expanded_ex(const Expr& e) {
  // FIXME: Bare-bone version tailored to verification needs. Add recursion
  //        (to simplify sums nested into other expressions) and variable
  //        indices.
  if (e.is_a_add()) {
    // Preliminary check: see whether we have sums.
    Expr addend;
    bool work_on_this = false;
    for (unsigned int i = e.nops(); i-- > 0; ) {
      addend = e.op(i);
      // FIXME: is_the_sum_function() does not recognize -sum(i,1,n,i).
      if (addend.is_the_sum_function()) {
	if (addend.arg(1) == 1 && addend.arg(2) == Recurrence::n) {
	  work_on_this = true;
	  break;
	}
      }
    }
    if (!work_on_this)
      return e;
    // Find and merge symbolic sums going from 1 to n.
    Expr e_without_sums = 0;
    Expr e_sum = 0;
    Symbol sum_variable;
    unsigned int merged_sums=0;
    for (unsigned int i = e.nops(); i-- > 0; ) {
      addend = e.op(i);
      if (addend.is_the_sum_function()) {
	if (addend.arg(1) == 1 && addend.arg(2) == Recurrence::n) {
	  merged_sums++;
	  // FIXME: add a coefficient when adding support for a*sum(...).
	  e_sum += addend.arg(3).substitute(addend.arg(0), sum_variable);
	}
	else {
	  e_without_sums += addend;
	};
      }
      else e_without_sums += addend;
    }
    // FIXME: Casting to Expr seems necessary.
    e_sum = sum((Expr) sum_variable, (Expr) 1, (Expr) Recurrence::n, e_sum);
    // Return the original expression unless we actually managed to merge
    // some sums.
    if (merged_sums > 1)
      return e_without_sums + e_sum;
    else
      return e;
  }
  else return e;
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
  // Since we need both numerator and denominator, calling 'numer_denom'
  // is faster than using 'numer()' and 'denom()' separately.
  Expr numer_e;
  Expr denom_e;
  e.numerator_denominator(numer_e, denom_e);
  Expr num = numer_e.expand();
  Expr den = denom_e.expand();
  return num * pwr(den, -1);
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
PURRS::simplify_binomials_factorials_exponentials(const Expr& e) {
#if 0
  Expr e_numerator;
  Expr e_denominator;
  e.numerator_denominator(e_numerator, e_denominator);
  e_numerator = rewrite_factorials_and_exponentials(e_numerator);
  e_denominator = rewrite_factorials_and_exponentials(e_denominator);
  return e_numerator / e_denominator;
#else
  Expr e_rewritten = rewrite_binomials(e);
  return rewrite_factorials_and_exponentials(e_rewritten);
#endif
}

PURRS::Expr
PURRS::simplify_logarithm(const Expr& e) { 
  return simplify_logarithm_in_expanded_ex(e.expand()).expand();
}

PURRS::Expr
PURRS::simplify_sum(const Expr& e,
		    const Sum_Simplification_Kind simplification) {
  return simplify_sum_in_expanded_ex(e.expand(), simplification).expand();
}

PURRS::Expr
PURRS::simplify_collect_sums(const Expr& e) {
  return simplify_collect_sums_in_expanded_ex(e.expand()).expand();
}

/*!
  Executes consecutively all simplifications described in the comment
  of the functions <CODE>simplify_numer_denom()</CODE>,
  <CODE>simplify_binomials_factorials_exponentials()</CODE>,
  <CODE>simplify_ex_for_output()</CODE> and
  <CODE>simplify_collect_sums()</CODE> and
  <CODE>simplify_logarithm()</CODE>.
*/
PURRS::Expr
PURRS::simplify_all(const Expr& e) {
  Expr e_rewritten = e;
  e_rewritten = simplify_numer_denom(e_rewritten);
  e_rewritten = simplify_binomials_factorials_exponentials(e_rewritten);
  e_rewritten = simplify_collect_sums(e_rewritten);
  e_rewritten = simplify_ex_for_output(e_rewritten, false);
  e_rewritten = simplify_logarithm(e_rewritten);
  return e_rewritten;
}
