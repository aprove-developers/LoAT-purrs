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

#include "simplify.hh"

#include "util.hh"
#include "Expr.defs.hh"
#include <vector>

// TEMPORARY
#include <iostream>

namespace Parma_Recurrence_Relation_Solver {

#ifndef NOISY
#define NOISY 0
#endif

static const unsigned
FACTOR_THRESHOLD = 100;

Expr
simplify_on_input_ex(const Expr& e, const Symbol& n, const bool& input);

Expr
simplify_on_output_ex(const Expr& e, const Symbol& n, const bool& input);


/*!
  Separates numeric factors and not numeric factors of an expression \p e. 
*/
static void
split_exponent(Expr& num, Expr& not_num, const Expr& e) {
  if (!e.is_a_number())
    not_num *= e;
  else
    num *= e;
}

/*!
  There are three cases:
  1. the <CODE>Expr</CODE> \p not_num_exponent is a <CODE>GiNaC::mul</CODE>
     and has like a factor \p n:
     returns the <CODE>Expr</CODE> \p not_num_exp_minus_n with
     all the factors of \p not_num_exponent minus \p n;
  2. the <CODE>Expr</CODE> \p not_num_exponent is not a
     <CODE>GiNaC::mul</CODE> and it is equal to \p n:
     returns the <CODE>Expr</CODE> \p not_num_exp_minus_n equal to \f$ 1 \f$; 
  3. the <CODE>Expr</CODE> \p not_num_exponent not contains \p n:
     returns the <CODE>Expr</CODE> \p not_num_exp_minus_n equal to
     \p not_num_exponent.
*/
static Expr
found_and_erase_n(const Expr& not_num_exponent, const Symbol& n) {
 bool found = false;
 if (not_num_exponent.is_a_mul()) {
   for (unsigned i = not_num_exponent.nops(); i-- > 0; )
     if (not_num_exponent[i].is_equal(n))
       found = true;
 }
 else
   if (not_num_exponent.is_equal(n))
     found = true;

 Expr not_num_exp_minus_n = 1;
 if (found) {
   if (not_num_exponent.is_a_mul()) {
     for (unsigned i = not_num_exponent.nops(); i-- > 0; )
       if (!not_num_exponent[i].is_equal(n))
	 not_num_exp_minus_n *= not_num_exponent[i];
   }
 }
 else
   not_num_exp_minus_n = not_num_exponent;
 
 return not_num_exp_minus_n;
}

/*!
  Returns <CODE>true</CODE> if \p base and \p exp_num are rational numbers
  and \f$ base^exp_num \f$ is again a rational number; <CODE>false</CODE>
  otherwise (for example \f$ base = 3 \f$ and \f$ exp_num = 1/2 \f$).
*/
static bool
perfect_root(const Expr& base, const Number& exp_num) {
  if (exp_num.is_rational()) {
    Expr pow_base_num_exp = Parma_Recurrence_Relation_Solver::power(base, exp_num);
    if (pow_base_num_exp.is_a_number()) {
      Number tmp = pow_base_num_exp.ex_to_number();
      if (tmp.is_rational())
	return true;
    }
  }
  return false;
}

/*!
  Given the base \p base, the numeric and not numeric part of the
  exponent, \p num_exp and \p not_num_exp respectively, this function
  returns the right power in according with some conditions checked
  by the boolean \p input and \p is_numeric_base.
*/
static Expr
return_power(const bool& is_numeric_base, const bool& input,
	     const Expr& num_exp, const Expr& not_num_exp,
	     const Expr& base, const Symbol& n) {
  Expr not_num_exp_minus_n;
  if (input)	
    not_num_exp_minus_n = found_and_erase_n(not_num_exp, n);
  // We do not want put in evidence the special symbol `n' or it is
  // not in `vect_not_num_exp[i]'.
  if (!input || not_num_exp_minus_n.is_equal(not_num_exp))
    if (is_numeric_base)
      return Parma_Recurrence_Relation_Solver::power(Parma_Recurrence_Relation_Solver::power(base, num_exp), not_num_exp);
    else
      return Parma_Recurrence_Relation_Solver::power(base, num_exp * not_num_exp);
  // We put in evidence the special symbol `n'. 
  else
    if (is_numeric_base)
      return Parma_Recurrence_Relation_Solver::power(Parma_Recurrence_Relation_Solver::power(Parma_Recurrence_Relation_Solver::power(base, num_exp), not_num_exp_minus_n), n);
    else
      return Parma_Recurrence_Relation_Solver::power(Parma_Recurrence_Relation_Solver::power(base, num_exp * not_num_exp_minus_n), n);
}

/*!
  \p base is a <CODE>GiNaC::mul</CODE>.
  We consider three vectors with a number of elements like the factor's
  number of \p base: \p vect_base, \p vect_num_exp and \p vect_not_num_exp.
  Initially \p vect_num_exp and \p vect_not_num_exp will have each element
  equal to \p num_exponent or \p not_num_exponents, respectively.
  Looks each factor of \p e and there are two cases:
  -  if it is a power, puts its base in the respective position of the vector
     \p vect_base and upgrades the rispective values of \p vect_num_exp and
     \p vect_not_num_exp with the exponent of \p base's factor;
  -  if it is not a power, puts it in the respective position of the vector
     \p vect_base and left unchanged the others vector.
  If \p input is <CODE>true</CODE> then we use the simplifications
  which put in evidence the special symbol \p n; otherwise, i. e. \p input
  is <CODE>false</CODE>, \p n is like the other parameters.
 */
static Expr
simpl_powers_base(const Expr& base, const Expr& num_exponent,
                  const Expr& not_num_exponent, const Symbol& n,
		  const bool& input) {
  std::vector<Expr> vect_base(base.nops());
  std::vector<Expr> vect_num_exp(base.nops());
  std::vector<Expr> vect_not_num_exp(base.nops());
  for (unsigned i = base.nops(); i-- > 0; ) {
    vect_num_exp[i] = num_exponent;
    vect_not_num_exp[i] = not_num_exponent;
  }
  for (unsigned i = base.nops(); i-- > 0; )
    if (base.op(i).is_a_power()) {
      Expr tmp = base.op(i);
      while (tmp.is_a_power()) {
        // The exponent of the factor `base.op(i)' is a multiplication.
        if (tmp.op(1).is_a_mul())
          for (unsigned j = tmp.op(1).nops(); j-- > 0; )
            split_exponent(vect_num_exp[i], vect_not_num_exp[i],
                           tmp.op(1).op(j));
        // The exponent of the factor `base.op(i)' is not a multiplication.
        else
          split_exponent(vect_num_exp[i], vect_not_num_exp[i], tmp.op(1));
        tmp = tmp.op(0);
      }
      vect_base[i] = tmp;
    }
    else
      vect_base[i] = base.op(i);

  // Now, for each factor of the base, is individualized numeric
  // and not numeric part of its exponent.
  Expr tot = 1;
  for (unsigned i = base.nops(); i-- > 0; )
    if (!vect_base[i].is_a_number())
      tot *= return_power(false, input, vect_num_exp[i], vect_not_num_exp[i],
			  vect_base[i], n);
    else
      tot *= return_power(true, input, vect_num_exp[i], vect_not_num_exp[i],
			  vect_base[i], n);
  return tot;
}

/*!
  Applies the rules \f$ \textbf{E1}, \textbf{E2}, \textbf{E4} \f$ and
  \f$ \textbf{E5} \f$ of the rules'set <EM>Expand</EM>.
  The <CODE>Expr</CODE> \p e is a <CODE>GiNaC::power</CODE>:
  it finds the base and the exponent of the power (\p e could be a serie
  of nested powers). While it does this operation divides the exponents
  (that can be multiplications but not addictions because the expression
  \p e is expanded) in two parts: in \p num_exponent put
  numeric factors and in \p not_num_exponent put not numeric factors.
  Therefore tests the base: if it is not a multiplication the checks and the
  simplifications are finished, otherwise we must elevate every factor of the
  base to the exponents.
  If \p input is <CODE>true</CODE> then we use the simplifications
  which put in evidence the special symbol \p n; otherwise, i. e. \p input
  is <CODE>false</CODE>, \p n is like the other parameters.
*/
static Expr
pow_simpl(const Expr& e, const Symbol& n, const bool& input) {
  Expr num_exponent = 1;
  Expr not_num_exponent = 1;
  Expr base = e;
  // At the end of the following cycle `base' will contains the `real' base
  // of the power.
  while (base.is_a_power()) {
    // Checks and divides the exponents in two part: those numeric and  
    // those not numeric.
    if (base.op(1).is_a_mul())
      for (unsigned i = base.op(1).nops(); i-- > 0; )
        split_exponent(num_exponent, not_num_exponent, base.op(1).op(i));
    else
      split_exponent(num_exponent, not_num_exponent, base.op(1));
    base = base.op(0);
  };
#if NOISY
  std::cout << "base " << base << std::endl;
  std::cout << "num_exp " << num_exponent << std::endl;
  std::cout << "not_num_exp " << not_num_exponent << std::endl;
#endif

  // The base is a multiplication.
  if (base.is_a_mul())
    return simpl_powers_base(base, num_exponent, not_num_exponent, n, input);
  // The base is not a multiplication: is not necessary to use the vectors,
  // i.e., call the function `simpl_powers_base'.
  else {
    if (base.is_a_number()) {
      Number num_exp = num_exponent.ex_to_number();
      // The function `perfect_root' allows to apply the rule `E1'.
      if (num_exp.is_integer() || perfect_root(base, num_exp))
	return return_power(true, input, num_exp, not_num_exponent,
			    base, n);
      else
	return return_power(false, input, num_exp, not_num_exponent,
			    base, n);
    }
    else
      return return_power(false, input, num_exponent, not_num_exponent,
			  base, n);
  }
}

/*!
  Applies the rule \f$ \textbf{C3} \f$ of the set of rules
  <EM>Collect</EM> to the <CODE>Expr</CODE> \p e  under condition
  that the common exponent to the powers is not integer
  because, in this case, <CODE>GiNaC</CODE> automatically decomposes the
  power, i. e., \f$ (a*b)^4 \f$ is automatically transformed in
  \f$ a^4*b^4 \f$.
  Returns a new <CODE>Expr</CODE> \p ris containing the modified
  expression \p e and the modified vectors \p bases and \p exponents will
  be use by the function <CODE>collect_same_exponent()</CODE> called soon
  afterwards this.
*/
static Expr
collect_same_exponents(const Expr& e, std::vector<Expr>& bases,
		       std::vector<Expr>& exponents) {
  assert(e.is_a_mul());
  Expr ris = 1;
  // At the end of the cycle `ris' will contain the powers of `e', with
  // the same exponents, simplified in only one power with the base equal
  // to the previous bases multiplicated among themselves (rule `C3').
  unsigned i = exponents.size();
  while (i > 0) {
    --i;
    if (!exponents[i].is_zero()) {
      for (unsigned j = i; j-- > 0; )	
	// In the vector `bases' in position `i' we put the new base.
	// In the vactor `exponents' in position `j' we put `0'.
	if (exponents[j].is_equal(exponents[i])) {
	  bases[i] *= bases[j];
	  exponents[j] = 0;
	}
      ris *= Parma_Recurrence_Relation_Solver::power(bases[i], exponents[i]);
    }
  }
  // FIXME: si potrebbe migliorare togliendo da `exponents'
  // gli elementi nulli e quelli nelle stesse posizioni da `bases'?
  // Ne vale la pena?

  // Now adds to `ris' the factor of `e' not considered in the
  // previous simplifications, i.e., the factor which are not powers.
  for (i = e.nops(); i-- > 0; )
    if (!e.op(i).is_a_power())
      ris *= e.op(i);
  return ris;
}

/*!
  Applies the rules \f$ \textbf{C1} \f$ and \f$ \textbf{C2} \f$ of the set
  of rules <EM>Collect</EM> to the <CODE>Expr</CODE> \p e that is
  certainly a <CODE>GiNaC::mul</CODE>.
  The vectors \p bases and \p exponents contain rispectively all bases and
  exponents of the powers that are in \p e and, at the end, will contain
  the new bases and exponents of the powers in \p e after the simplification.
  This function is called after <CODE>collect_same_exponents()</CODE>.
  Returns a new <CODE>Expr</CODE> \p ris containing the modified
  expression \p e.
*/
static Expr
collect_same_base(const Expr& e, std::vector<Expr>& bases,
		  std::vector<Expr>& exponents) {
  assert(e.is_a_mul());
  // At the end of the cycle `while', `ris' will contain all the powers of `e'.
  // The powers with the same bases will simplified in only one power with
  // the summed exponents (rule `C2').
  Expr ris = 1;
  unsigned i = bases.size();
  while (i > 0) {
    --i;
    if (!exponents[i].is_zero()) {
      for (unsigned j = i; j-- > 0; )
	// In the vector `exponents' in position `i' we put the new
	// exponent while in opsition `j' we put `0'.	
	if (bases[j].is_equal(bases[i])) {
	  exponents[i] += exponents[j];
	  exponents[j] = 0;
	}
      ris *= Parma_Recurrence_Relation_Solver::power(bases[i], exponents[i]);
    }
  }
  // FIXME: si potrebbe migliorare togliendo da `exponents'
  // gli elementi nulli e quelli nelle stesse posizioni da `bases'?
  // Ne vale la pena?

  // Now adds to `ris' the factor of `e' not considered in the
  // previous simplification, i.e., the factor which are not powers .
  for (i = e.nops(); i-- > 0; ) {
    if (!e.op(i).is_a_power()) {
      // We must consider those factors that are not powers but are equal
      // to some base of the vector `bases', in order to apply the rule `C1',
      // i.e., `a * a^e = a^(e + 1)'. In this case, we add `1' to the
      // correspondenting exponent.
      bool to_sum = false;
      for (unsigned j = bases.size(); j-- > 0; )
	if (bases[j].is_equal(e.op(i)) && !exponents[j].is_zero()) {
	  to_sum = true;
	  // If the base is `integer' GiNaC automatically transforms for
	  // instance `2^(3/2)' in `2*sqrt(2)': in this case we do not add
	  // 1 to the exponent. 
	  if (!bases[j].is_a_number() || !exponents[j].is_a_number())
	    exponents[j] = exponents[j] + 1;
	  else {
	    Number base = bases[j].ex_to_number();
	    if (!base.is_integer() || base == -1)
	      exponents[j] = exponents[j] + 1;
	  }
	}
      // Applies rule `C1'.
      if (to_sum)
	ris = ris.subs(Parma_Recurrence_Relation_Solver::power(e.op(i), wild(0)), Parma_Recurrence_Relation_Solver::power(e.op(i), wild(0) + 1));
      else
	ris *= e.op(i);
    }
  }
  return ris;
}

/*!
  Applies the rules of the set <EM>Collect</EM> to <CODE>Expr</CODE>
  \p e, that is certainly a <CODE>GiNaC::mul</CODE>.
  Returns a new <CODE>Expr</CODE> containing the modified expression \p e. 
*/
static Expr
collect_base_exponent(const Expr& e) {
  assert(e.is_a_mul());
  Expr tmp = e;
  // Builds two vectors containing the bases and the exponents of
  // the eventual multiplication's factors which are powers. 
  std::vector<Expr> bases;
  std::vector<Expr> exponents;
  for (unsigned i = tmp.nops(); i-- > 0; )
    if (tmp.op(i).is_a_power()) {
      bases.push_back(tmp.op(i).op(0));
      exponents.push_back(tmp.op(i).op(1));
    }
  // We have a better simplification if we apply `collect_same_exponents()'
  // before than `collect_same_base()' (ex. `2*2^n*(1/2)^n').
  // Applies rule `C3'.
  tmp = collect_same_exponents(tmp, bases, exponents);
#if NOISY
  std::cout << "tmp dopo same exponents... " << tmp << std::endl;
#endif
  // After the simplifications by the function `collect_same_exponents()'
  // `tmp' could not be a `mul'. 
  if (tmp.is_a_mul())
    // Applies rules `C1' and `C2'.    
    tmp = collect_same_base(tmp, bases, exponents);
#if NOISY
  std::cout << "tmp dopo same base... " << tmp << std::endl;
#endif
  return tmp;
}

/*!
  Construct a partial factorization of the integer \p n.
  \p n is tested for divisibility by 2 and by odd integers between 3
  and <CODE>FACTOR_THRESHOLD</CODE>.
  The partially factored form is returned in the pair of vectors 
  \p bases and \p exponents, of <CODE>Number</CODE>s
  and <CODE>int</CODE>s respectively.
*/
static void 
partial_factor(const Number n, std::vector<Number>& bases,
	       std::vector<int>& exponents) {
  
  assert(n.is_integer());
  Number m = abs(n);
  assert(m != 0);
  int k = 0;
  while (mod(m, 2) == 0) { // the case 2 is handled separately 
    m /= 2;
    ++k;
  }
  if (k>0) {
    bases.push_back(2);
    exponents.push_back(k);
  }
  for (unsigned i=3; (i < FACTOR_THRESHOLD) && (i*i <= m); i += 2) {
    k = 0;
    while (mod(m, i) == 0) { // test for divisibility by the odd integer i
      m /= i;
      ++k;
    }
    if (k > 0) {
      bases.push_back(i);
      exponents.push_back(k);
    }
  }
  if (m > 1) { // here n has not necessarily been factored completely 
    bases.push_back(m);
    exponents.push_back(1);
  }
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
static Expr
to_std_form(const Number k, const std::vector<Number>& bases, 
	    std::vector<int>& exponents) {
  
  assert(k != 0);
  assert(k.is_integer());
  int abs_k = ::abs(k.to_int());
  Expr m = 1;
  for (unsigned i = 0; i < bases.size(); ++i) {
    int remainder = exponents[i] % abs_k;
    int quotient  = exponents[i] / abs_k;
    if (( k < 0 ) && ( remainder != 0 )) { // adjust quotient and remainder
      ++quotient;
      remainder = abs_k - remainder;
    }
    exponents[i] = remainder;
    m *= Parma_Recurrence_Relation_Solver::power(bases[i], quotient);
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
static Expr 
reduce_to_standard_form(const Number root_index, const Number r) {
  assert(root_index.is_integer());
  assert(root_index != 0);
  int k = root_index.to_int();
  if (k < 0)
    assert(r != 0);
  if (k > 0)
    if (r == 0)
      return 0;
  Number num = r.numerator();
  Number den = r.denominator();

  // FIXME: deal with complex numbers
  if (!r.is_real()) {
    Expr index = 1 / root_index;
    return Parma_Recurrence_Relation_Solver::power(r, index);
  }
  Number sign = r > 0 ? 1 : -1;
  Number g = gcd(num, den);
  num /= g;
  den /= g; // clear any common factors from num, den.
  if ((k % 2 == 0) && (sign == -1)) // complex sign if k is even and r < 0.
    sign = k > 0 ? Number::I : -Number::I;
  if (k < 0) { // swap numerator and denominator, and change sign of k.
    Number i = num;
    num = den;
     den = i;
     k *= -1;
  } // now, num, den and k are all positive.
  
  if (k == 1)
    return sign * num * Parma_Recurrence_Relation_Solver::power(den, -1);
  
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
  for (unsigned i=0; (i < num_size) && (gc > 1); ++i)
    gc = gcd(gc, num_exponents[i]);
  for (unsigned i=0; (i < den_size) && (gc > 1); ++i)
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
    irr_part *= Parma_Recurrence_Relation_Solver::power(num_bases[i], num_exponents[i]);
  for (unsigned i = 0; i < den_size; ++i)
    irr_part *= Parma_Recurrence_Relation_Solver::power(den_bases[i], den_exponents[i]);
  Expr q = sign * reduced_num * Parma_Recurrence_Relation_Solver::power(reduced_den, -1);
  if (irr_part.ex_to_number() > 1)
    q *= Parma_Recurrence_Relation_Solver::power(irr_part, Number(1, k));
  return q;
}

/*!
  ...
*/
static Expr
red_prod(const Number& base1, const Number& exp1, 
	 const Number& base2, const Number& exp2) {
  Number base_1 = base1;  
  Number base_2 = base2;
  assert(exp1 != 0);
  assert(exp2 != 0);
  
  Number k1_num = exp1.numerator();
  Number k1_den = exp1.denominator();
  Number k2_num = exp2.numerator();
  Number k2_den = exp2.denominator();
  
  base_1 = Parma_Recurrence_Relation_Solver::power(base_1, k1_num);
  base_2 = Parma_Recurrence_Relation_Solver::power(base_2, k2_num);
  
  Number g = gcd(k1_den, k2_den);
  Number k = k1_den * k2_den / g;
  Number b1 = Parma_Recurrence_Relation_Solver::power(base_1, k2_den / g);
  Number b2 = Parma_Recurrence_Relation_Solver::power(base_2, k1_den / g);
  Number b = b1 * b2;

  return reduce_to_standard_form(k, b);
}

/*!
  Applies the rules of the set <EM>Irrationals</EM> to
  <CODE>Expr</CODE> \p e if \p e is a <CODE>GiNaC::power</CODE> or to
  each \p e's factor which is a <CODE>GiNaC::power</CODE> if \p e is a
  <CODE>GiNaC::mul</CODE>.
*/
static Expr
reduce_product(const Expr& e) {
  if (e.is_a_mul()) {
    Expr tmp = e;
    Expr factor_to_reduce = 1;
    Expr factor_no_to_reduce = 1;
    Number base_1 = 1;
    Number exp_1 = 1;
    for (unsigned i = tmp.nops(); i-- > 0; )
      if (tmp.op(i).is_a_power()) {
	// Base and exponent of `tmp.op(i)' are both numerics.
	if (tmp.op(i).op(0).is_a_number() &&
	    tmp.op(i).op(1).is_a_number()) {
	  Number base_2 = tmp.op(i).op(0).ex_to_number();
	  Number exp_2  = tmp.op(i).op(1).ex_to_number();
	  Expr to_reduce = red_prod(base_1, exp_1, base_2, exp_2);
	  // red_prod returns `numeric' or `numeric^numeric' or
	  // `numeric * numeric^numeric'.
	  if (to_reduce.is_a_mul()) {
	    assert(to_reduce.nops() == 2);
	    for (unsigned j = 2; j-- > 0; )
	      if (to_reduce.op(j).is_a_power()) {
		base_1 = to_reduce.op(j).op(0).ex_to_number();
		exp_1  = to_reduce.op(j).op(1).ex_to_number();
		factor_to_reduce = Parma_Recurrence_Relation_Solver::power(to_reduce.op(j).op(0),
				       to_reduce.op(j).op(1));
	      }
	      else
		factor_no_to_reduce *= to_reduce.op(j);
	  }
	  else if (to_reduce.is_a_power()) {
	    base_1 = to_reduce.op(0).ex_to_number();
	    exp_1 = to_reduce.op(1).ex_to_number();
	    factor_to_reduce = Parma_Recurrence_Relation_Solver::power(to_reduce.op(0), to_reduce.op(1));
	  }
	  else {
	    base_1 = to_reduce.ex_to_number();
	    exp_1 = 1;
	    factor_to_reduce = Parma_Recurrence_Relation_Solver::power(to_reduce, 1);
	  }
	}
	// Base and exponent of `tmp.op(i)' are not both numerics.
	else
	  factor_no_to_reduce *= tmp.op(i);
      }
    // `tmp.op(i)' is not a `GiNaC::power'.
      else
	factor_no_to_reduce *= tmp.op(i);
    return factor_to_reduce * factor_no_to_reduce;
  }
  else if (e.is_a_power())
    if (e.op(0).is_a_number() && e.op(1).is_a_number()) {
      Number base = e.op(0).ex_to_number();
      Number exp  = e.op(1).ex_to_number();
      return red_prod(base, exp, 1, 1);
    }
  return e;
}

/*!
  Applies all rules of term rewriting system \f$ \mathfrak{R}_o \f$ which
  are applicable on factors.
  Returns a <CODE>Expr</CODE> that contains the modified expression \p e.
*/
static Expr
manip_factor(const Expr& e, const Symbol& n, const bool& input) {
  assert(e.is_a_mul());
  Expr tmp = 1;

  // Simplifies each factor that is a `GiNaC::power'.
  for (unsigned i = e.nops(); i-- > 0; )
    if (e.op(i).is_a_power()) {
      Expr base = simplify_on_output_ex(e.op(i).op(0), n, input);
      Expr exp = simplify_on_output_ex(e.op(i).op(1), n, input);
      if (base.is_a_number() && exp.is_a_number())
	tmp *= reduce_product(Parma_Recurrence_Relation_Solver::power(base, exp));
      else
	tmp *= pow_simpl(Parma_Recurrence_Relation_Solver::power(base, exp), n, input);
    }
    else
      tmp *= e.op(i);
#if NOISY
  std::cout << "tmp dopo nested... " << tmp << std::endl;
#endif
  // From this time forward we do not know if `tmp' is a again `mul'.
  
  // Simplifies recursively the factors which are functions simplifying
  // their arguments.
  if (tmp.is_a_mul()) {
    Expr factor_function = 1;
    Expr factor_no_function = 1;
    for (unsigned i = tmp.nops(); i-- > 0; )
      if (tmp.op(i).is_a_function()) {
	Expr argument = simplify_on_output_ex(tmp.op(i).op(0), n, input);
	factor_function *= tmp.op(i).subs(tmp.op(i).op(0), argument);
      }
      else
	factor_no_function *= tmp.op(i);
    tmp = factor_function * factor_no_function;
  }
  else if (tmp.is_a_function()) {
    Expr argument = simplify_on_output_ex(tmp.op(0), n, input);
    tmp = tmp.subs(tmp.op(0), argument);
  }
#if NOISY
  std::cout << "tmp dopo function... " << tmp << std::endl << std::endl;
#endif
  // Special case: the exponential `exp' is a `GiNaC::function' but it has
  // the same properties of the powers.
  if (tmp.is_a_mul()) {
    Expr argument = 0;
    Expr rem = 1;
    for (unsigned i = tmp.nops(); i-- > 0; ) {
      Expr_List l;
      if (tmp.op(i).match(exp(wild(0)), l))
	argument += l.op(0).rhs();
      else
	rem *= tmp.op(i);
    }
    tmp = exp(argument) * rem;
#if NOISY
    std::cout << "tmp dopo `exp'... " << tmp << std::endl << std::endl;
#endif
  }
  
  // Simplifies eventual powers with same base or same exponents.
  if (tmp.is_a_add()) {
    Expr terms = 0;
    for (unsigned i = tmp.nops(); i-- > 0; ) 
      if (tmp.op(i).is_a_mul())
	terms += collect_base_exponent(tmp.op(i));
      else
	terms += tmp.op(i);
    tmp = terms;
  }
  else
    if (tmp.is_a_mul())
      tmp = collect_base_exponent(tmp);
  
  return tmp;
}

/*!
  Crosses the tree of the expanded expression \p e recursevely to find
  subexpressions to which we apply the rules of the terms rewriting system
  \f$ \mathfrak{R}_i \f$. More exactly here the rules of the set
  \emph{Expand} are implemented because the rules of the set <EM>Automatic</EM>
  are automatically executed by <CODE>GiNaC</CODE>.
  We observe that the rule \f$ \textbf{E4} \f$ is automatically executed
  by <CODE>GiNaC</CODE> if the exponent is integer while the rules
  \f$ \textbf{E3} \f$ and \f$ \textbf{E6} \f$ are executed by the method
  <CODE>expand()</CODE> (\f$ \textbf{E3} \f$ only partially because for
  instance \f$ expand(3^(4*x+2*a)) = 3^(2*a)*3^(4*x) \f$):
  hence we here consider only <CODE>GiNaC::power</CODE>.
  \p input is always <CODE>true</CODE> and this means that \p n is a special
  symbol, i. e., in the simplifications is always put in evidence in respect
  of the others parameters.
  Returns a <CODE>Expr</CODE> that contains the modified expression \p e.
*/
Expr
simplify_on_input_ex(const Expr& e, const Symbol& n, const bool& input) {
  Expr ris;
  if (e.is_a_add()) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      ris += simplify_on_input_ex(e.op(i), n, input);
  }
  else if (e.is_a_mul()) {
    ris = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      ris *= simplify_on_input_ex(e.op(i), n, input);
  }
  else if (e.is_a_power())
      return pow_simpl(e, n, input);
  else if (e.is_a_function()) {
    Expr f = e;
    Expr tmp = simplify_on_input_ex(e.op(0), n, input);
    ris = f.subs(f.op(0), tmp);
  }
  else
    ris += e;

  return ris;
}

/*!
  Crosses the tree of the expanded expression \p e recursevely to find
  subexpressions which we want to apply the rules of the terms rewriting system
  \f$ \mathfrak{R}_o \f$. The observations about the function
  <CODE>simplify_on_input_ex()</CODE> are correct here too, because all the
  rules of the term rewriting system \f$ \mathfrak{R}_i \f$, minus
  \f$ \textbf{E2} \f$, are also in \f$ \mathfrak{R}_o \f$.
  \p input is always <CODE>false</CODE> and this means that \p n is not
  considerated like a special symbol but like the others parameters.
  Returns a <CODE>Expr</CODE> that contains the modified expression \p e.
*/
Expr
simplify_on_output_ex(const Expr& e, const Symbol& n, const bool& input) {
  Expr ris;
  if (e.is_a_add()) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      ris += simplify_on_output_ex(e.op(i), n, input);
  }
  else if (e.is_a_mul())
    // We can not call `simplify_on_output_ex' on every factor because
    // otherwise it is not possible to transform products.
    ris = manip_factor(e, n, input);
  else if (e.is_a_power()) {
    Expr base = simplify_on_output_ex(e.op(0), n, input);
    Expr exp = simplify_on_output_ex(e.op(1), n, input);
    if (base.is_a_number() && exp.is_a_number())
      ris = reduce_product(Parma_Recurrence_Relation_Solver::power(base, exp));
    else
      ris = pow_simpl(Parma_Recurrence_Relation_Solver::power(base, exp), n, input);
    // Necessary for l'output: for example if `e = sqrt(18)^a' then
    // `ris = sqrt(2)^a*3^a'.
    if (ris.is_a_mul())
      ris = collect_base_exponent(ris);
  }
  else if (e.is_a_function()) {
    Expr f = e;
    Expr tmp = simplify_on_output_ex(e.op(0), n, input);
    ris = f.subs(f.op(0), tmp);
  }
  else
    ris += e;
  return ris;
}


/*!
  Using the function <CODE>GiNaC::numer_denom()</CODE> we are able to
  simplify better rational expression.
*/
Expr
simplify_numer_denom(const Expr& e) {
  // Since we need both numerator and denominator, to call 'numer_denom'
  // is faster than to use 'numer()' and 'denom()' separately.
  Expr numer_e;
  Expr denom_e;
  e.numerator_denominator(numer_e, denom_e);
  Expr num = numer_e.expand();
  Expr den = denom_e.expand();
  Expr ris = num * Parma_Recurrence_Relation_Solver::power(den, -1);
  return ris;
}


/*!
  Given \p a, \p b and \p c elements of the factorial expression
  \f$ (a n + b)! \f$, this function returns a new expression that put
  in evidence \f$ (a n)! \f$.
*/
static Expr
decompose_factorial(const Symbol& n,
		    const Expr& a, const Expr& b, const Expr& c) {
  Expr prod = c * factorial(a*n);
  if (b.is_a_number()) {
    Number tmp_b = b.ex_to_number();
    if (tmp_b > 0)
      for (Number j = 1; j <= tmp_b; ++j)
	prod *= a*n + j; 
    else
      for (Number j = 0; j < abs(tmp_b); ++j)
	prod *= Parma_Recurrence_Relation_Solver::power(a*n - j, -1);
  }
  return prod;
}

/*!
  Applies the following rewrite rule where \f$ a \in \Nset \setminus \{0\} \f$
  and \f$ b \in \Zset \f$:
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
static Expr
rewrite_factorials(const Expr& e, const Symbol& n) {
  Expr new_e;
  Expr fact_of_sum = factorial(n + wild(0));
  Expr a_times_fact_of_sum = wild(1) * fact_of_sum;
  Expr fact_of_sum_coeff = factorial(wild(0) * n + wild(1));
  Expr a_times_fact_of_sum_coeff = wild(2) * fact_of_sum_coeff;
  Expr_List substitution;
  if (e.is_a_add()) {
    new_e = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      if (clear(substitution), e.op(i).match(fact_of_sum, substitution))
	new_e += decompose_factorial(n, 1, get_binding(substitution, 0), 1);
      else if (clear(substitution),
	       e.op(i).match(a_times_fact_of_sum, substitution))
	new_e += decompose_factorial(n, 1,
				     get_binding(substitution, 0),
				     get_binding(substitution, 1));
      else if (clear(substitution),
	       e.op(i).match(fact_of_sum_coeff, substitution))
	new_e += decompose_factorial(n, get_binding(substitution, 0),
				     get_binding(substitution, 1), 1);
      else if (clear(substitution),
	       e.op(i).match(a_times_fact_of_sum_coeff, substitution))
	new_e += decompose_factorial(n, get_binding(substitution, 0),
				     get_binding(substitution, 1),
				     get_binding(substitution, 2));
      else
	new_e += e.op(i);
  }
  else if (e.is_a_mul()) {
    new_e = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      if (clear(substitution), e.op(i).match(fact_of_sum, substitution))
	new_e *= decompose_factorial(n, 1, get_binding(substitution, 0), 1);
      else if (clear(substitution),
	       e.op(i).match(fact_of_sum_coeff, substitution))
	new_e *= decompose_factorial(n, get_binding(substitution, 0),
				     get_binding(substitution, 1), 1);
      else
	new_e *= e.op(i);
  }
  else {
    new_e = 0;
    if (clear(substitution), e.match(fact_of_sum, substitution))
      new_e += decompose_factorial(n, 1, get_binding(substitution, 0), 1);
    else if (clear(substitution), e.match(fact_of_sum_coeff, substitution))
      new_e += decompose_factorial(n, get_binding(substitution, 0),
				   get_binding(substitution, 1), 1);
    else
      new_e += e;
  }
  return new_e;
}

/*!
  Applies the following rewrite rule
  \f[
    a^{b n + c} = a^{b n} \cdot a^c
  \f]
  without to expand all expression \p e.
*/
static Expr
simpl_exponentials(const Expr& e, const Symbol& n) {
  Expr new_e;
  if (e.is_a_add()) {
    new_e = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      new_e += simpl_exponentials(e.op(i), n);
  }
  else if (e.is_a_mul()) {
    new_e = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      if (e.op(i)
	  .match(Parma_Recurrence_Relation_Solver::power(wild(0), n + wild(1)))
	  || e.op(i)
	  .match(Parma_Recurrence_Relation_Solver::power(wild(0),
							 wild(2) * n + wild(1))))
	new_e *= e.op(i).expand();
      else
	new_e *= e.op(i);
  }
  else
    if (e.match(Parma_Recurrence_Relation_Solver::power(wild(0), n + wild(1)))
	|| e.match(Parma_Recurrence_Relation_Solver::power(wild(0),
							   wild(2) * n + wild(1))))
      new_e += e.expand();
    else
      new_e += e;
  return new_e;
}

/*!
  Let \p e be an expression containing ratios of factorials or ratios
  of exponentials.
  This function tries in the numerator and denominator of \p e eventual
  factorials of the type \f$ (a n + b)! \f$, whit
  \f$ a \in \Nset \setminus \{0\} \f$ and \f$ b \in \Zset \f$, and
  eventual exponentials in order to put in evidence common factors to
  numerator and denominator to erase.
  We remark that, for this type of simplifications, is not good to call
  the function <CODE>simplify_on_output_ex()</CODE> because it expanded
  the expression while we want to mantain the product of factors. 
  Returns a <CODE>Expr</CODE> that contains the modified expression \p e.
*/
Expr
simplify_factorials_and_exponentials(const Expr& e, const Symbol& n) {
  Expr e_numerator;
  Expr e_denominator;
  e.numerator_denominator(e_numerator, e_denominator);

  e_numerator = rewrite_factorials(e_numerator, n);
  e_numerator = simpl_exponentials(e_numerator, n);
  e_denominator = rewrite_factorials(e_denominator, n);
  e_denominator = simpl_exponentials(e_denominator, n);

  Expr new_e = e_numerator / e_denominator;

  return new_e;
}

} // namespace Parma_Recurrence_Relation_Solver
