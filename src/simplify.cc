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

#include "simplify.hh"
#include <vector>

// TEMPORARY
#include <iostream>

using namespace GiNaC;

#define NOISY 0

static const unsigned
FACTOR_THRESHOLD = 100;

GExpr
simplify_on_output_ex(const GExpr& e);

GExpr
simplify_on_input_ex(const GExpr& e);



static void
split_exponent(GExpr& num, GExpr& not_num, const GExpr& e) {
  if (!is_a<numeric>(e))
    not_num *= e;
  else
    num *= e;
}

/*!
  Applies the rules 9 and 10 of the terms rewriting system
  \f$ \mathfrak{R}_o \f$.
  Returns a new <CODE>GExpr</CODE> <CODE>ris</CODE> with the eventual
  simplifications.
 */
static GExpr
collect_same_base(const GExpr& e, const std::vector<GExpr>& bases,
		  std::vector<GExpr>& exponents) {
  assert(is_a<mul>(e));
  // At the end of the cycle 'ris' will contain the powers of 'e', with
  // the same bases, simplified in only one power with the exponents summed
  // (rule 10).
  GExpr ris = 1;
  std::vector<GExpr> tmp_exp(exponents);
  unsigned i = bases.size();
  while (i > 0) {
    --i;
    GExpr exp = tmp_exp[i];
    GExpr base = bases[i];
    for (unsigned j = i; j-- > 0; )
      if (bases[j].is_equal(base)) {
	exp += tmp_exp[j];
	tmp_exp[j] = 0;
	exponents[j] = 0;
	exponents[i] = 0;
      }
    ris *= pow(base, exp);
    tmp_exp[i] = 0;
  }
  // Now adds to 'ris' the factor of 'e' not considered in the
  // previous simplifications.
  for (i = e.nops(); i-- > 0; ) {
    if (!is_a<power>(e.op(i))) {
      bool to_sum = false;
      for (unsigned j = bases.size(); j-- > 0; )
	if (bases[j].is_equal(e.op(i))) {
	  to_sum = true;
	  if (!is_a<numeric>(bases[j]))
	    exponents[j] = exponents[j] + 1;
	  else {
	    GNumber base = GiNaC::ex_to<GiNaC::numeric>(bases[j]);
	    if (!base.is_integer())
	      exponents[j] = exponents[j] + 1;
	  }
	}
      // (rule 9)
      if (to_sum)
	ris = ris.subs(pow(e.op(i), wild()) == pow(e.op(i), (wild() + 1)));
      else
	ris *= e.op(i);
    }
    else {
      bool insert = true;
      for (unsigned j = bases.size(); j-- > 0; )
	if (bases[j].is_equal(e.op(i).op(0)))
	  insert = false;
      if (insert)
	ris *= e.op(i);
    }
  }
  return ris;
}

/*!
  Applies the rule 12 of the terms rewriting system \f$ \mathfrak{R}_o \f$
  under condition that the common exponent of the powers is not integer
  because, in this case, <CODE>GiNaC</CODE> automatically decomposes the
  power.
  Returns a new <CODE>GExpr</CODE> <CODE>ris</CODE> with the eventual
  simplifications.
 */
static GExpr
collect_same_exponents(const GExpr& e, std::vector<GExpr>& bases,
		       std::vector<GExpr>& exponents) {
  assert(is_a<mul>(e));
  GExpr ris = 1;
  // At the end of the cycle 'ris' will contain the powers of 'e', with
  // the same exponents, simplified in only one power with the base equal
  // to the two previous base multiplicated (rule 12).
  unsigned i = exponents.size();
  while (i > 0) {
    --i;
    if (exponents[i] != 0) {
      GExpr exp = exponents[i];
      GExpr base = bases[i];
      bool found = false;
      for (unsigned j = i; j-- > 0; )
	if (exponents[j].is_equal(exp)) {
	  found = true;
	  base *= bases[j];
	  exponents[j] = 0;
	}
      if (found)
	ris *= pow(base, exp);
      else
	exponents[i] = 0;
    }
  }
  for (i = e.nops(); i-- > 0; )
    if (!is_a<power>(e.op(i)))
      ris *= e.op(i);
    else {
      bool insert = true;
      for (unsigned j = exponents.size(); j-- > 0; )
	if (exponents[j].is_equal(e.op(i).op(1)))
	  insert = false;
      if (insert)
	ris *= e.op(i);
    }
  return ris;
}

static GExpr
collect_base_exponent(const GExpr& e) {
  GExpr tmp = e;
  assert(is_a<mul>(tmp));
  // Builds two vector that contain the bases and the exponents of
  // the eventual multiplication's factors that are powers. 
  std::vector<GExpr> bases;
  std::vector<GExpr> exponents;
  for (unsigned i = tmp.nops(); i-- > 0; )
    if (is_a<power>(tmp.op(i))) {
      bases.push_back(tmp.op(i).op(0));
      exponents.push_back(tmp.op(i).op(1));
    }
  // Rules 9 and 10.    
  tmp = collect_same_base(tmp, bases, exponents);
#if NOISY
  std::cout << "tmp dopo same base... " << tmp << std::endl;
#endif
  // Rule 12.
  if (is_a<mul>(tmp))
    tmp = collect_same_exponents(tmp, bases, exponents);
#if NOISY
  std::cout << "tmp dopo same exponents... " << tmp << std::endl;
#endif
  return tmp;
}

/*!
  \p e is a multiplication and at least one of the factors is a
  <CODE>GiNaC::power</CODE>.
  We consider three vectors with a number of elements like the factor's
  number of \p e: \p vect_base, \p vect_num_exp and \p vect_not_num_exp.
  Initially \p vect_num_exp and \p vect_not_num_exp will have each element
  equal to \p num_exponent or \p not_num_exponents, respectively.
  Looks each factor of \p e and, if it is a power, puts its base in
  the respective position of the vector \p vect_base and upgrades the
  rispective values of \p vect_num_exp and \p vect_not_num_exp with the
  exponent of \p e's factor.
  If \p e's factor is not a power puts it in the respective position
  of the vector \p vect_base and left unchanged the others vector.
 */
static GExpr
simpl_powers_base(const GExpr& e, const GExpr& num_exponent,
                  const GExpr& not_num_exponent) {
  std::vector<GExpr> vect_base(e.nops());
  std::vector<GExpr> vect_num_exp(e.nops());
  std::vector<GExpr> vect_not_num_exp(e.nops());
  for (unsigned i = e.nops(); i-- > 0; ) {
    vect_num_exp[i] = num_exponent;
    vect_not_num_exp[i] = not_num_exponent;
  }
  for (unsigned i = e.nops(); i-- > 0; )
    if (is_a<power>(e.op(i))) {
      GExpr tmp = e.op(i);
      while (is_a<power>(tmp)) {  
        // The exponent of the factor 'e.op(i)' is a multiplication.
        if (is_a<mul>(tmp.op(1)))
          for (unsigned j = tmp.op(1).nops(); j-- > 0; )
            split_exponent(vect_num_exp[i], vect_not_num_exp[i],
                           tmp.op(1).op(j));
        // The exponent of the factor 'e.op(i)' is not a multiplication.
        else
          split_exponent(vect_num_exp[i], vect_not_num_exp[i], tmp.op(1));
        tmp = tmp.op(0);
      }
      vect_base[i] = tmp;
    }
    else
      vect_base[i] = e.op(i);
  GExpr tot = 1;
  for (unsigned i = e.nops(); i-- > 0; )
    if (!is_a<numeric>(vect_base[i]))
      tot *= pow(vect_base[i], vect_num_exp[i] * vect_not_num_exp[i]);
    else
      tot *= pow(pow(vect_base[i], vect_num_exp[i]), vect_not_num_exp[i]);
  return tot;
}

/*!
  Applies the rule \f$ 6, 7 \f$ and \f$ 8 \f$ of the terms rewriting system
  \f$ \mathfrak{R}_o \f$.
  The <CODE>GExpr</CODE> \p e is a <CODE>GiNaC::power</CODE>:
  finds the base and the exponent of the power (\p e could be a serie
  of nested powers). While does this operation divides the exponents
  (that can be multiplications but not addictions because the expression
  \p e is expanded) in two parts: in \p num_exponent put
  numeric factors and in \p not_num_exponent put not numeric factors.
  Then checks the base: if it is not a multiplication checks and
  simplifications are finished, otherwise we must raise every factor of the
  base to the exponents. In the case there is some base's factor that is a
  <CODE>GiNaC::power</CODE> we must do further simplifications.
*/
static GExpr
pow_simpl(const GExpr& e) {
  GExpr num_exponent = 1;
  GExpr not_num_exponent = 1;
  GExpr b = e;
  // At the end of the following cycle 'b' will contains the base of
  // the power.
  while (is_a<power>(b)) {
    // Checks and divides the exponents in two part: those numeric and  
    // those not numeric.
    if (is_a<mul>(b.op(1)))
      for (unsigned i = b.op(1).nops(); i-- > 0; )
        split_exponent(num_exponent, not_num_exponent, b.op(1).op(i));
    else
      split_exponent(num_exponent, not_num_exponent, b.op(1));
    b = b.op(0);
  };
#if NOISY
  std::cout << "base " << b << std::endl;
  std::cout << "num_exp " << num_exponent << std::endl;
  std::cout << "not_num_exp " << not_num_exponent << std::endl;
#endif
  // Simplifies evental powers with same base or with same exponent
  // in 'not_num_exponent'.
  if (is_a<mul>(not_num_exponent))
    not_num_exponent = collect_base_exponent(not_num_exponent);
  // The base is a multiplication.
  if (is_a<mul>(b)) {
    GExpr tot = 1;
    bool base_power = false;
    for (unsigned i = b.nops(); i-- > 0; )
      if (is_a<power>(b.op(i)))
	base_power = true;
    // Some factor of the base is a power.
    if (base_power)
      tot = simpl_powers_base(b, num_exponent, not_num_exponent);
    // No factor of the base is a power.
    else {
      for (unsigned i = b.nops(); i-- > 0; )
	if (!is_a<numeric>(b.op(i)))
	  tot *= pow(b.op(i), num_exponent * not_num_exponent);
	else
	  tot *= pow(pow(b.op(i), num_exponent), not_num_exponent);
    }
    return tot;
  }
  // The base is not a multiplication.
  else
    if (is_a<numeric>(b)) {
      GNumber exp_num = GiNaC::ex_to<GiNaC::numeric>(num_exponent);
      if (exp_num.is_integer()) 
	return pow(pow(b, exp_num), not_num_exponent);
      else
	return pow(b, exp_num * not_num_exponent);
    }
    else
      return pow(b, num_exponent * not_num_exponent);
}

/*!
  Compute the gcd between the integers \f$n\f$ and \f$m\f$.
*/
static int
gcd(int n, int m) {
  
  int r = abs(m);
  while (r != 0) {
    r = n % m;
    n = m; 
    m = r;
  }
  return n;
}

/*!
  Construct a partial factorization of the integer \f$n\f$.
  \f$n\f$ is tested for divisibility by 2 and by odd integers between 3
  and <CODE>FACTOR_THRESHOLD</CODE>.
  The partially factored form is returned in the pair of vectors 
  <CODE>bases</CODE> and <CODE>exponents</CODE>, of <CODE>GNumber</CODE>s
  and <CODE>int</CODE>s respectively.
*/

static void 
partial_factor(const GNumber n, std::vector<GNumber>& bases,
	       std::vector<int>& exponents) {
  
  assert(n.is_integer());
  GNumber m = abs(n);
  assert(m != 0);
  int k = 0;
  while (irem(m, 2) == 0) { // the case 2 is handled separately 
    m /= 2;
    ++k;
  }
  if (k>0) {
    bases.push_back(2);
    exponents.push_back(k);
  }
  for (unsigned i=3; (i < FACTOR_THRESHOLD) && (i*i <= m); i += 2) {
    k = 0;
    while (irem(m, i) == 0) { // test for divisibility by the odd integer i
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
  <CODE>partial_factor</CODE> is given. 
  If \f$ k > 0 \f$, the exponents in the vector <CODE>exponents</CODE> are 
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

static GExpr
to_std_form(const GNumber k, const std::vector<GNumber>& bases, 
	    std::vector<int>& exponents) {
  
  assert(k != 0);
  assert(k.is_integer());
  int abs_k = abs(k.to_int());
  GExpr m = 1;
  for (unsigned i = 0; i < bases.size(); ++i) {
    int remainder = exponents[i] % abs_k;
    int quotient  = exponents[i] / abs_k;
    if (( k < 0 ) && ( remainder != 0 )) { // adjust quotient and remainder
      ++quotient;
      remainder = abs_k - remainder;
    }
    exponents[i] = remainder;
    m *= pow(bases[i], quotient);
  }
  return m;
}

/*!
  Transform \f$ r^{1/k} \f$ into its standard form. 
  If \f$ r = n / d \f$ where \f$ n \f$ and \f$ d \f$ are non-zero, coprime 
  integers, the standard form of \f$ r^{1/k} \f$ is defined as follows: 
  1. if \f$ n / d < 0 \f$ and \f$ k \f$ is odd, it is the opposite of the 
  standard form of \f$ |n / d| ^ (1/k) \f$. 
  2. if \f$ n / d < 0 \f$ and \f$ k \f$ is even, it is \f$ I \f$ times the 
  standard form of \f$ |n / d| ^ (1/k) \f$. 
  3. if \f$ k < 0 \f$ it is the standard form of \f$ (d/n)^{1/(-k)} \f$. 
  4. if \f$ n \f$ and \f$ d \f$ are both positive, the standard form is 
  the expression \f$ (n_1 / d_1) \cdot m^{1/k} \f$ where \f$ n_1 \f$, 
  \f$ d_1 \f$ and \f$ m \f$ are positive integers and \f$ m \f$ is 
  \f$ k\f$-free, that is, \f$ m \f$ is not divisible by the 
  \f$ k \f$-th power of any integer larger than 1. 
  (Note that if \f$ k = 1 \f$ then necessarily \f$ m = 1 \f$). 
*/

static GExpr 
reduce_to_standard_form(const GNumber root_index, const GNumber r) {
  
  assert(root_index.is_integer());
  assert(root_index != 0);
  int k = root_index.to_int();
  if (k < 0)
    assert(r != 0);
  if (k > 0)
    if (r == 0)
      return 0;
  GExpr   n_d = numer_denom(r);
  GNumber num = GiNaC::ex_to<GiNaC::numeric>(n_d.op(0));
  GNumber den = GiNaC::ex_to<GiNaC::numeric>(n_d.op(1));
  // FIXME: deal with complex numbers
  if (!r.is_real()) {
    GExpr index = 1 / root_index;
    return pow(r, index);
  }
  GNumber sign = r > 0 ? 1 : -1;
  GNumber g    = gcd(num, den);
  num /= g;
  den /= g; // clear any common factors from num, den
  if ((k % 2 == 0) && (sign == -1)) // complex sign if k is even and r < 0
    sign = k > 0 ? I : -I;
  if (k < 0) { // swap numerator and denominator, and change sign of k 
    GNumber i = num;
    num = den;
     den = i;
     k *= -1;
  } // now, num, den and k are all positive
  
  if (k == 1)
    return sign * num * pow(den, -1);
  
  std::vector<GNumber> num_bases;
  std::vector<int> num_exponents;
  std::vector<GNumber> den_bases;
  std::vector<int> den_exponents;
  
  // partial factor and reduce numerator and denominator
  partial_factor(num, num_bases, num_exponents);
  unsigned num_size = num_bases.size();
  GExpr reduced_num = to_std_form(k, num_bases, num_exponents);
  
  // here <CODE>to_std_form</CODE> is called with a negative value of k 
  // because we are dealing with the denominator of r 
  partial_factor(den, den_bases, den_exponents);
  unsigned den_size = den_bases.size();
  GExpr reduced_den = to_std_form(-k, den_bases, den_exponents);
  
  // Try one last simplification: if all exponents have a common factor 
  // with the root index, remove it
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
  
  GExpr irr_part = 1;
  for (unsigned i=0; i < num_size; ++i)
    irr_part *= pow(num_bases[i], num_exponents[i]);
  for (unsigned i=0; i < den_size; ++i)
    irr_part *= pow(den_bases[i], den_exponents[i]);
  
  GExpr q = sign * reduced_num * pow(reduced_den, -1);
  if (irr_part > 1)
    q *= pow(irr_part, numeric(1) / k);
  
  return q;
}

static GExpr
red_prod(const GNumber& base1, const GNumber& exp1, 
	 const GNumber& base2, const GNumber& exp2) {
  GNumber base_1 = base1;  
  GNumber base_2 = base2;
  assert(exp1 != 0);
  assert(exp2 != 0);
  
  GExpr   n_d_1  = numer_denom(exp1);
  GNumber k1_num = GiNaC::ex_to<GiNaC::numeric>(n_d_1.op(0));
  GNumber k1_den = GiNaC::ex_to<GiNaC::numeric>(n_d_1.op(1));
  GExpr   n_d_2  = numer_denom(exp2);
  GNumber k2_num = GiNaC::ex_to<GiNaC::numeric>(n_d_2.op(0));
  GNumber k2_den = GiNaC::ex_to<GiNaC::numeric>(n_d_2.op(1));
  
  base_1 = pow(base_1, k1_num);
  base_2 = pow(base_2, k2_num);
  
  GNumber g  = gcd(k1_den, k2_den);
  GNumber k  = k1_den * k2_den / g;
  GNumber b1 = GiNaC::ex_to<GiNaC::numeric>(pow(base_1, k2_den / g));
  GNumber b2 = GiNaC::ex_to<GiNaC::numeric>(pow(base_2, k1_den / g));
  GNumber b = b1 * b2;
  return reduce_to_standard_form(k, b);
}

static GExpr
reduce_product(const GExpr& a) {
  if (is_a<mul>(a)) {
    GExpr tmp = a;
    GExpr factor_to_reduce = 1;
    GExpr factor_no_to_reduce = 1;
    GNumber base_1 = 1;
    GNumber exp_1 = 1;
    for (unsigned i = tmp.nops(); i-- > 0; )
      if (is_a<power>(tmp.op(i))) {
	// Base and exponent of 'tmp.op(i)' are both numeric.
	if (GiNaC::is_a<GiNaC::numeric>(tmp.op(i).op(0)) &&
	    GiNaC::is_a<GiNaC::numeric>(tmp.op(i).op(1))) {
	  GNumber base_2 = GiNaC::ex_to<GiNaC::numeric>(tmp.op(i).op(0));
	  GNumber exp_2  = GiNaC::ex_to<GiNaC::numeric>(tmp.op(i).op(1));
	  GExpr to_reduce = red_prod(base_1, exp_1, base_2, exp_2);
	  // red_prod restituisce 
	  // numerico            oppure
	  // numerico^numerico   oppure
	  // numerico * numerico^numerico
	  // es. 3^(1/4)*6^(3/4) restituisce 3*8^(1/4)
	  // FIXME: giusto?
	  if (is_a<mul>(to_reduce)) {
	    assert(to_reduce.nops() == 2);
	    for (unsigned j = 2; j-- > 0; )
	      if (is_a<power>(to_reduce.op(j))) {
		assert(is_a<numeric>(to_reduce.op(j).op(0)));
		assert(is_a<numeric>(to_reduce.op(j).op(1)));
		base_1 = GiNaC::ex_to<GiNaC::numeric>(to_reduce.op(j).op(0));
		exp_1  = GiNaC::ex_to<GiNaC::numeric>(to_reduce.op(j).op(1));
		factor_to_reduce = pow(to_reduce.op(j).op(0),
				       to_reduce.op(j).op(1));
	      }
	      else
		factor_no_to_reduce *= to_reduce.op(j);
	  }
	  else if (is_a<power>(to_reduce)) {
	    assert(is_a<numeric>(to_reduce.op(0)));
	    assert(is_a<numeric>(to_reduce.op(1)));
	    base_1 = GiNaC::ex_to<GiNaC::numeric>(to_reduce.op(0));
	    exp_1 = GiNaC::ex_to<GiNaC::numeric>(to_reduce.op(1));
	    factor_to_reduce = pow(to_reduce.op(0), to_reduce.op(1));
	  }
	  else {
	    assert(is_a<numeric>(to_reduce));
	    base_1 = GiNaC::ex_to<GiNaC::numeric>(to_reduce);
	    exp_1 = 1;
 	  factor_to_reduce = pow(to_reduce, 1);
	  }
	}
	// Base and exponent of 'tmp.op(i)' are not both numeric.
	else
	  factor_no_to_reduce *= tmp.op(i);
      }
    // 'tmp.op(i)' is not a 'GiNaC::power'.
      else
	factor_no_to_reduce *= tmp.op(i);
    return factor_to_reduce * factor_no_to_reduce;
  }
  else if (is_a<power>(a))
    if (GiNaC::is_a<GiNaC::numeric>(a.op(0)) &&
	GiNaC::is_a<GiNaC::numeric>(a.op(1))) {
      GNumber base = GiNaC::ex_to<GiNaC::numeric>(a.op(0));
      GNumber exp  = GiNaC::ex_to<GiNaC::numeric>(a.op(1));
      return red_prod(base, exp, 1, 1);
    }
  return a;
}

/*!
  Give an expression \p e builds two other expressions \p numerica and
  \p symbolic that contain the numeric and the symbolic part of \p e,
  rispectively.  
*/
static void
split(const GExpr& e, GExpr& numerica, GExpr& symbolic) {
  if (is_a<add>(e)) {
    GExpr tmp_num_term = 1;
    GExpr tmp_symb_term = 1;
    GExpr tmp_num = 0;
    GExpr tmp_symb = 0;
    for (unsigned i = e.nops(); i-- > 0; ) {
      split(e.op(i), tmp_num_term, tmp_symb_term);
      if (tmp_symb_term.is_equal(1))
	tmp_num += e.op(i);
      else
	tmp_symb += e.op(i);
    }
    if (tmp_symb.is_zero())
      numerica *= e;
    else
      symbolic *= e;
  }
  else if (is_a<mul>(e))
    for (unsigned i = e.nops(); i-- > 0; )
      split(e.op(i), numerica, symbolic);
  else if (is_a<function>(e)) {
    GExpr tmp_num = 1;
    GExpr tmp_symb = 1;
    split(e.op(0), tmp_num, tmp_symb);
    if (tmp_symb.is_equal(1))
      numerica *= e;
    else
      symbolic *= e;
  }
  else if (is_a<power>(e)) {
    GExpr tmp_num_b = 1;
    GExpr tmp_symb_b = 1;
    GExpr tmp_num_e = 1;
    GExpr tmp_symb_e = 1;
    split(e.op(0), tmp_num_b, tmp_symb_b);
    split(e.op(1), tmp_num_e, tmp_symb_e);
    if (tmp_symb_b.is_equal(1) && tmp_symb_e.is_equal(1))
      numerica *= e;
    else
      symbolic *= e;
  }
  else {
    if (is_a<numeric>(e))
      numerica *= e;
    else
      symbolic *= e;
  }
}

static GExpr
manip_factor(const GExpr& e) {
  assert(is_a<mul>(e));
  GExpr tmp = 1;
  // Simplifies each factor that is a 'GiNaC::power'.
  for (unsigned i = e.nops(); i-- > 0; )
    if (is_a<power>(e.op(i))) {
      GExpr base = simplify_on_output_ex(e.op(i).op(0));
      GExpr exp = simplify_on_output_ex(e.op(i).op(1));
      tmp *= pow_simpl(pow(base, exp));
    }
    else
      tmp *= e.op(i);
#if NOISY
  std::cout << "tmp dopo nested... " << tmp << std::endl;
#endif
  // Goes recursively to simplify factors that are functions simplifying
  // their arguments.
  if (is_a<mul>(tmp)) {
    GExpr factor_function = 1;
    GExpr factor_no_function = 1;
    for (unsigned i = tmp.nops(); i-- > 0; )
      if (is_a<function>(tmp.op(i))) {
	GExpr argument = simplify_on_output_ex(tmp.op(i).op(0));
	factor_function *= tmp.op(i).subs(tmp.op(i).op(0) == argument);
      }
      else
	factor_no_function *= tmp.op(i);
    tmp = factor_function * factor_no_function;
  }
  else if (is_a<function>(tmp)) {
    GExpr argument = simplify_on_output_ex(tmp.op(0));
    tmp = tmp.subs(tmp.op(0) == argument);
  }
#if NOISY
  std::cout << "tmp dopo function... " << tmp << std::endl << std::endl;
#endif
  // Special case: the exponential 'exp' is a GiNaC::function but have
  // the same properties of the other powers.
  if (is_a<mul>(tmp)) {
    GExpr argument = 0;
    GExpr rem = 1;
    for (unsigned i = tmp.nops(); i-- > 0; ) {
      GList l;
      if (tmp.op(i).match(exp(wild()), l))
	argument += l.op(0).rhs();
      else
	rem *= tmp.op(i);
    }
    tmp = exp(argument) * rem;
#if NOISY
    std::cout << "tmp dopo 'exp'... " << tmp << std::endl << std::endl;
#endif
  }
  // Divides numeric and symbolic part of the expression 'e'.
  GExpr numerica = 1;
  GExpr symbolic = 1;
  split(tmp, numerica, symbolic);
#if NOISY
  std::cout << std::endl << "symbolic " << symbolic << std::endl;
  std::cout << "numerica " << numerica << std::endl << std::endl;
#endif
  // Simplifies eventual powers with same base or same exponents in
  // a symbolic part.
  if (is_a<mul>(symbolic))
    symbolic = collect_base_exponent(symbolic);
#if NOISY
  std::cout << std::endl << "symbolic dopo simpl " << symbolic << std::endl;
#endif
  // If there are in 'symbolic' factors numeric or powers with base
  // and exponent numeric, then multiplies this factors to 'numerica'.
  if (is_a<mul>(symbolic))
    for (unsigned i = symbolic.nops(); i-- > 0; ) {
      if (is_a<power>(symbolic.op(i))) {
	if (is_a<numeric>(symbolic.op(i).op(0))
	    && is_a<numeric>(symbolic.op(i).op(1))) {
	  numerica *= symbolic.op(i);
	  symbolic = symbolic.subs(symbolic.op(i) == 1);
	}
      }
      else
	if (is_a<numeric>(symbolic.op(i))) {
	  numerica *= symbolic.op(i);
	  symbolic = symbolic.subs(symbolic.op(i) == 1);
	}
    }
#if NOISY
  std::cout << std::endl << "symbolic dopo simpl " << symbolic << std::endl;
  std::cout << "numerica " << numerica << std::endl << std::endl;
#endif
  // Simplifies eventual powers with same base or same exponents in
  // a numeric part.
  if (is_a<mul>(numerica))
    numerica = collect_base_exponent(numerica);
#if NOISY
  std::cout << std::endl << "numerica dopo simpl " << numerica << std::endl;
#endif
  // Simplifies roots.
  numerica = reduce_product(numerica);
#if NOISY
  std::cout << std::endl << "numerica dopo roots " << numerica << std::endl;
#endif

  return symbolic * numerica;
}

/*!
  Crosses the tree of the expanded expression \p e recursevely to find
  subexpressions to which we apply the rules of the terms rewriting system
  \f$ /mathfrak{R}_i \f$. More exactly here the rules \f$ 6, 7 \f$ and
  \f$ 8 \f$ of the system are implemented because the first four roules are
  automatically executed by <CODE>GiNaC</CODE> and the rule \f$ 5 \f$ is
  partially executed by the method <CODE>expand()</CODE> (partially because
  \f$ expand(3^(4*x+2*a) = 3^(2*a)*3^(4*x) \f$).
  Returns a <CODE>GExpr</CODE> that contains the modified expression \p.
*/
GExpr
simplify_on_input_ex(const GExpr& e) {
  GExpr ris;
  if (is_a<add>(e)) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      ris += simplify_on_input_ex(e.op(i));
  }
  else if (is_a<mul>(e)) {
    ris = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      ris *= simplify_on_input_ex(e.op(i));
  }
  else if (is_a<power>(e))
      return pow_simpl(e);
  else if (is_a<function>(e)) {
    GExpr f = e;
    GExpr tmp = simplify_on_input_ex(e.op(0));
    ris = f.subs(f.op(0) == tmp);
  }
  else
    ris += e;
  return ris;
}

/*!
  Crosses the tree of the expanded expression \p e recursevely to find
  subexpressions to which we apply the rules of the terms rewriting system
  \f$ /mathfrak{R}_o \f$. Beyond to the made observations about the function
  \p simplify_on_input_ex that are valid also here we observe that also
  the rule \f$ 11 \f$ is automatically executed by <CODE>GiNaC</CODE> while
  here are implemented the rules \f$ 9, 10, 12, 13, 14, 15, 16 \f$.
  Returns a <CODE>GExpr</CODE> that contains the modified expression \p.
*/
GExpr
simplify_on_output_ex(const GExpr& e) {
  GExpr ris;
  if (is_a<add>(e)) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      ris += simplify_on_output_ex(e.op(i));
  }
  else if (is_a<function>(e)) {
    GExpr f = e;
    GExpr tmp = simplify_on_output_ex(e.op(0));
    ris = f.subs(f.op(0) == tmp);
  }
  else if (is_a<power>(e)) {
    GExpr base = simplify_on_output_ex(e.op(0));
    GExpr exp = simplify_on_output_ex(e.op(1));
    if (is_a<numeric>(base) && is_a<numeric>(exp)) {
      GNumber base_1 = GiNaC::ex_to<GiNaC::numeric>(base);
      GNumber exp_1 = GiNaC::ex_to<GiNaC::numeric>(exp);
      ris = red_prod(1, 1, base_1, exp_1);
    }
    else
      ris = pow_simpl(pow(base, exp));
  }
  else if (is_a<mul>(e)) {
    ris = manip_factor(e);
  }
  else
    ris += e;
  return ris;
}
