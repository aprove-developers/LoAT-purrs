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
#include "util.hh"
#include <vector>

// TEMPORARY
#include <iostream>

using namespace GiNaC;

#define NOISY 0

static const unsigned
FACTOR_THRESHOLD = 100;

static void
split_exponent(GExpr& num, GExpr& not_num, const GExpr& e) {
  if (!is_a<numeric>(e))
    not_num *= e;
  else
    num *= e;
}

static GExpr
get_binding(const GList& l, unsigned wild_index) {
  assert(wild_index < l.nops());
  assert(l.op(wild_index).info(GiNaC::info_flags::relation_equal));
  return l.op(wild_index).rhs();
}

/*!
  \f$ f \f$ is a multiplication and at least one of the factors is a
  <CODE>GiNaC::power</CODE>.
  We consider three vectors with a number of elements like the factor's
  number of \f$ f \f$: <CODE>vect_base</CODE>, <CODE>vect_num_exp</CODE>
  and <CODE>vect_not_num_exp</CODE>.
  Initially <CODE>vect_num_exp</CODE> and <CODE>vect_not_num_exp</CODE>
  will have each element equal to <CODE>num_exponent</CODE> or
  <CODE>not_num_exponents</CODE>, respectively.
  Looks each factor of \f$ f \f$ and, if it is a power, puts its base in
  the respective position of the vector <CODE>vect_base</CODE> and
  upgrades the rispective values of <CODE>vect_num_exp</CODE> and
  <CODE>vect_not_num_exp</CODE> with the exponent of \f$ f \f$'s factor.
  If \f$ f \f$'s factor is not a power puts it in the respective position
  of the vector <CODE>vect_base</CODE> and left unchanged the others vector.
 */
static GExpr
simpl_powers_base(const GExpr& f, const GExpr& num_exponent,
                  const GExpr& not_num_exponent) {
  std::vector<GExpr> vect_base(f.nops());
  std::vector<GExpr> vect_num_exp(f.nops());
  std::vector<GExpr> vect_not_num_exp(f.nops());
  for (unsigned i = f.nops(); i-- > 0; ) {
    vect_num_exp[i] = num_exponent;
    vect_not_num_exp[i] = not_num_exponent;
  }
  for (unsigned i = f.nops(); i-- > 0; )
    if (is_a<power>(f.op(i))) {
      GExpr tmp = f.op(i);
      while (is_a<power>(tmp)) {  
        // The exponent of the factor 'f.op(i)' is a multiplication.
        if (is_a<mul>(tmp.op(1)))
          for (unsigned j = tmp.op(1).nops(); j-- > 0; )
            split_exponent(vect_num_exp[i], vect_not_num_exp[i],
                           tmp.op(1).op(j));
        // The exponent of the factor 'f.op(i)' is not a multiplication.
        else
          split_exponent(vect_num_exp[i], vect_not_num_exp[i], tmp.op(1));
        tmp = tmp.op(0);
      }
      vect_base[i] = tmp;
    }
    else
      vect_base[i] = f.op(i);
  GExpr tot = 1;
  for (unsigned i = f.nops(); i-- > 0; )
    if (!is_a<numeric>(vect_base[i]))
      tot *= pow(vect_base[i], vect_num_exp[i] * vect_not_num_exp[i]);
    else
      tot *= pow(pow(vect_base[i], vect_num_exp[i]), vect_not_num_exp[i]);
  return tot;
}

/*!
  Applies the rule 8, 9 and 10 of the terms rewriting system.
  The <CODE>GExpr</CODE> \f$ e \f$ is a <CODE>GiNaC::power</CODE>:
  finds the base and the exponent of the power (\f$ e \f$ could be a serie
  of nested powers). While does this operation divides the exponents
  (that can be multiplications but not addictions because the expression
  \f$ e \f$ is expanded) in two parts: in <CODE>num_exponent</CODE> put
  numeric factors and in <CODE>not_num_exponent</CODE> put not numeric
  factors.
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
    // No factors of the base is a power.
    else
      for (unsigned i = b.nops(); i-- > 0; )
	if (!is_a<numeric>(b.op(i)))
	  tot *= pow(b.op(i), num_exponent * not_num_exponent);
	else
	  tot *= pow(pow(b.op(i), num_exponent), not_num_exponent);
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
  Crosses the tree of the expression expanded \f$ e \f$ recursevely to find
  subexpressions which to apply the rule 8, 9 and 10 of the terms rewriting
  system \f$ /mathfrak{R} \f$.
  Returns a <CODE>GExpr</CODE> that contains the input expression \f$ e \f$
  modified.  
*/
//FIXME: da sistemare il caso per es. (cos(2^(3*x))^3 oppure
//       3^(cos(2^(3*x)))
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

//   else if (is_a<power>(e)) {
//     ris = 1;
//     if (!is_a<function>(e.op(0)))
//       return pow_simpl(e);
//     else    
//  }
  else if (is_a<function>(e)) {
    ris = 1;
    GExpr f = e;
    GExpr tmp = simplify_on_input_ex(e.op(0));
    ris *= f.subs(f.op(0) == tmp);
  }
  else
    ris += e;
  return ris;
}

/*!
  Applies the rules 5 and 6 of the terms rewriting system \f$ \mathfrak{R} \f$.
  Returns a new <CODE>GExpr</CODE> <CODE>ris</CODE> with the eventual
  simplifications.
 */
static GExpr
collect_same_base(const GExpr& e, const std::vector<GExpr>& bases,
		  std::vector<GExpr>& exponents) {
  assert(is_a<mul>(e));
  // At the end of the cycle 'ris' will contain the powers of 'e', with
  // the same bases, simplified in only one power with the exponents summed
  // (rule 6).
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
  for (i = e.nops(); i-- > 0; )
    if (!is_a<power>(e.op(i))) {
      bool to_sum = false;
      for (unsigned j = bases.size(); j-- > 0; )
	if (bases[j].is_equal(e.op(i)))
	  to_sum = true;
      // (rule 5.)
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
  return ris;
}

/*!
  Applies the rule 12 of the terms rewriting system \f$ \mathfrak{R} \f$
  under condition that the common exponent of the powers is not integer
  because, in this case, <CODE>GiNaC</CODE> automatically decomposes the
  power.
  Returns a new <CODE>GExpr</CODE> <CODE>ris</CODE> with the eventual
  simplifications.
 */
static GExpr
collect_same_exponents(const GExpr& e, std::vector<GExpr>& bases,
		       std::vector<GExpr>& exponents) {
  GExpr ris = 1;
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
  if (is_a<mul>(e))
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
  else
    ris *= e;
  return ris;
}

static GExpr
collect_base_exponent(const GExpr& e) {
  GExpr tmp = e;
  // Simplifies nested powers.
  for (unsigned i = tmp.nops(); i-- > 0; )
    if (is_a<power>(tmp.op(i)))
      if (is_a<power>(tmp.op(i).op(0))
	  || has(tmp.op(i).op(0), sqrt(wild())))
	tmp = tmp.subs(tmp.op(i) == pow_simpl(tmp.op(i)));
      else if (is_a<mul>(tmp.op(i).op(0)))
	tmp = tmp.subs(tmp.op(i).op(0)
		       == collect_base_exponent(tmp.op(i).op(0)));
#if NOISY
  std::cout << "tmp dopo nested... " << tmp << std::endl;
#endif
  GExpr factor_function = 1;
  for (unsigned i = tmp.nops(); i-- > 0; ) {
    // Goes recursively to simplify factors that are functions simplifying
    // the argument of the functions.
    // Put this factor simplified in a new GExpr 'factor_function' and
    // remove this factor from 'tmp'.
    if (is_a<function>(tmp.op(i))) {
      GExpr t = simplify_on_output_ex(tmp.op(i).op(0));
      factor_function *= tmp.op(i).subs(tmp.op(i).op(0) == t);
      tmp = tmp.subs(tmp.op(i) == 1);
    }
  }
#if NOISY
  std::cout << "factor_function... " << factor_function << std::endl;
  std::cout << "tmp dopo ... " << tmp << std::endl;
#endif
  // Now 'tmp' contain only the factor that are not functions
  // and are not powers with exponents not integer.
  GList lst;
  if (find(tmp, wild(0) * pow(wild(1), wild(2)), lst)) {
    std::vector<GExpr> bases;
    std::vector<GExpr> exponents;
    for (unsigned i = lst.op(0).nops(); i-- > 0; )
      if (is_a<power>(lst.op(0).op(i))) {
	bases.push_back(lst.op(0).op(i).op(0));
	exponents.push_back(lst.op(0).op(i).op(1));
      }
    // Rules 5 and 6.    
    tmp = collect_same_base(tmp, bases, exponents);
#if NOISY
    std::cout << "tmp dopo same base... " << tmp << std::endl;
#endif
    // Rule 12. 
    tmp = collect_same_exponents(tmp, bases, exponents);
#if NOISY
    std::cout << "tmp dopo same exponents... " << tmp << std::endl;
#endif
  }
  tmp *= factor_function;

  return tmp;
}

/*!
  Crosses the tree of the expression expanded \f$ e \f$ recursively to find
  subexpressions which to apply the rule 5, 6 and 12 of the terms rewriting
  system \f$ \mathfrak{R} \f$ (the rule 11 is automatically applied by GiNaC).
  Note that the function <CODE>collect_base_exponents()</CODE> calls the
  function that applies the rules 8, 9 and 10.
  Returns a <CODE>GExpr</CODE> that contains the input expression \f$ e \f$
  modified.  
*/ 
GExpr
simplify_on_output_ex(const GExpr& e) {
  GExpr ris;
  if (is_a<add>(e)) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      ris += simplify_on_output_ex(e.op(i));
  }
  else if (is_a<power>(e))
    // FIXME: la condizione qui dovrebbe essere
    // 'simplify_on_output_ex(pow_simpl(e))' ma cosi' ciclerebbero semplici
    // potenze come x^2. Se metto condizioni di arresto come
    // 'if (is_a<power>(e.op(0)))' non fa potenze come '2^(3*x)' che pero'
    // non posso gestire dopo ma solo contemporaneamente alle potenze annidate.
    // Cosi' pero' non mette a posto potenze come 'sqrt(x^a*5^a)' per la quale
    // occorrerebbero due giri. 
    return pow_simpl(e);
  else if (is_a<function>(e)) {
    ris = 1;
    GExpr f = e;
    GExpr tmp = simplify_on_output_ex(e.op(0));
    ris *= f.subs(f.op(0) == tmp);
  }
  else if (is_a<mul>(e))
    return collect_base_exponent(e);
  else
    ris += e;
  return ris;
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

// FIXME: rifare ricorsivamente...
static GExpr
reduce_product(const GExpr& a) {
  GExpr tmp = a;
  GList substitution;
  clear(substitution);
  if (match(tmp, pow(wild(0), wild(1)) * pow(wild(2), wild(3)) * wild(4),
	    substitution)) {
    GExpr base_1 = get_binding(substitution, 0);
    GExpr exp_1  = get_binding(substitution, 1);
    GExpr base_2 = get_binding(substitution, 2);
    GExpr exp_2  = get_binding(substitution, 3);
    GExpr factor = get_binding(substitution, 4);
    if (GiNaC::is_a<GiNaC::numeric>(base_1) && 
	GiNaC::is_a<GiNaC::numeric>(base_2) &&
	GiNaC::is_a<GiNaC::numeric>(exp_1)  &&
	GiNaC::is_a<GiNaC::numeric>(exp_2)) {
      GNumber base1 = GiNaC::ex_to<GiNaC::numeric>(base_1);
      GNumber exp1  = GiNaC::ex_to<GiNaC::numeric>(exp_1);
      GNumber base2 = GiNaC::ex_to<GiNaC::numeric>(base_2);
      GNumber exp2  = GiNaC::ex_to<GiNaC::numeric>(exp_2);
      tmp = factor * red_prod(base1, exp1, base2, exp2);
    }
  }
  else if (clear(substitution), match(tmp, pow(wild(0), wild(1)), 
				      substitution)) {
    GExpr base   = get_binding(substitution, 0);
    GExpr exp    = get_binding(substitution, 1);
    if (GiNaC::is_a<GiNaC::numeric>(base) && 
	GiNaC::is_a<GiNaC::numeric>(exp)) {
      assert(exp != 0);
      GNumber b = GiNaC::ex_to<GiNaC::numeric>(base);
      GNumber k = GiNaC::ex_to<GiNaC::numeric>(exp);
      GExpr   n_d = numer_denom(k);
      GNumber k_num = GiNaC::ex_to<GiNaC::numeric>(n_d.op(0));
      GNumber k_den = GiNaC::ex_to<GiNaC::numeric>(n_d.op(1));
      b = pow(b, k_num);
      tmp = reduce_to_standard_form(k_den, b);
     }
  }
  else if (clear(substitution), match(tmp, pow(wild(0), wild(1)) * wild(2), 
				      substitution)) {
    GExpr base   = get_binding(substitution, 0);
    GExpr exp    = get_binding(substitution, 1);
    GExpr factor = get_binding(substitution, 2);
    if (GiNaC::is_a<GiNaC::numeric>(base) && 
	GiNaC::is_a<GiNaC::numeric>(exp)) {
      assert(exp != 0);
      GNumber b = GiNaC::ex_to<GiNaC::numeric>(base);
      GNumber k = GiNaC::ex_to<GiNaC::numeric>(exp);
      GExpr   n_d = numer_denom(k);
      GNumber k_num = GiNaC::ex_to<GiNaC::numeric>(n_d.op(0));
      GNumber k_den = GiNaC::ex_to<GiNaC::numeric>(n_d.op(1));
      b = pow(b, k_num);
      tmp = factor * reduce_to_standard_form(k_den, b);
    }
  }
  return tmp;
}

//    // DIFFERENT VERSION OF reduce_product()
//    //
//    GExpr tmp = a;
//    GExpr factor_to_reduce = 1;
//    GExpr factor_no_to_reduce = 1;
//    GNumber base_1 = 1;
//    GNumber exp_1 = 1;
//    for (unsigned i = tmp.nops(); i-- > 0; )
//      if (is_a<power>(tmp.op(i)))
//        // Base and exponent of 'tmp.op(i)' are both numeric.
//        if (GiNaC::is_a<GiNaC::numeric>(tmp.op(i).op(0)) &&
//  	  GiNaC::is_a<GiNaC::numeric>(tmp.op(i).op(1))) {
//  	GNumber base_2 = GiNaC::ex_to<GiNaC::numeric>(tmp.op(i).op(0));
//  	GNumber exp_2 = GiNaC::ex_to<GiNaC::numeric>(tmp.op(i).op(1));
//  	factor_to_reduce = red_prod(base_1, exp_1, base_2, exp_2);
//  	// red_prod restituisce 
//  	// numerico            oppure
//  	// numerico^numerico   oppure
//  	// numerico * numerico^numerico
//  	// es. 3^(1/4)*6^(3/4) restituisce 3*8^(1/4)
//  	// FIXME: giusto?
//  	if (is_a<mul>(factor_to_reduce)) {
//  	  assert(factor_to_reduce.nops() == 2);
//  	  for (unsigned j = 2; j-- > 0; )
//  	    if (is_a<power>(factor_to_reduce.op(j))) {
//  	      assert(is_a<numeric>(factor_to_reduce.op(j).op(0)));
//  	      assert(is_a<numeric>(factor_to_reduce.op(j).op(1)));
//  	      base_1 = GiNaC::ex_to<GiNaC::numeric>(factor_to_reduce.op(j).op(0));
//  	      exp_1 = GiNaC::ex_to<GiNaC::numeric>(factor_to_reduce.op(j).op(1));
//  	    }
//  	}
//  	else if (is_a<power>(factor_to_reduce)) {
//  	  assert(is_a<numeric>(factor_to_reduce.op(0)));
//  	  assert(is_a<numeric>(factor_to_reduce.op(1)));
//  	  base_1 = GiNaC::ex_to<GiNaC::numeric>(factor_to_reduce.op(0));
//  	  exp_1 = GiNaC::ex_to<GiNaC::numeric>(factor_to_reduce.op(1));
//  	}
//  	else {
//  	  assert(is_a<numeric>(factor_to_reduce));
//  	  base_1 = GiNaC::ex_to<GiNaC::numeric>(factor_to_reduce);
//  	  exp_1 = 1;
//  	}
//        }
//    // Base and exponent of 'tmp.op(i)' are not both numeric.
//        else
//  	factor_no_to_reduce *= tmp.op(i);
//    // 'tmp.op(i)' is not a 'GiNaC::power'.
//      else
//        factor_no_to_reduce *= tmp.op(i);
//    return factor_to_reduce * factor_no_to_reduce;
//  }


GExpr
simplify_roots(const GExpr& e) {
  GExpr ris;
  if (is_a<add>(e)) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      ris += simplify_roots(e.op(i));
  }
  else if (is_a<function>(e)) {
    ris = 1;
    GExpr f = e;
    GExpr tmp = simplify_roots(e.op(0));
    ris *= f.subs(f.op(0) == tmp);
  }
  else if (is_a<mul>(e)) {
    // Se si usa la nuova reduce_product().
    //return reduce_product(e);

//      // Se si usa la vecchia reduce_product().
    ris = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      ris *= simplify_roots(reduce_product(e.op(i)));
    // FIXME: vedere il modo migliore e piu' efficiente...
    // Questo non lo e'!
    //ris = reduce_product(ris);
  }
  else if (is_a<power>(e)) {
    ris = 1;
    if (is_a<numeric>(e.op(0)) && is_a<numeric>(e.op(1))) {
      GNumber base = GiNaC::ex_to<GiNaC::numeric>(e.op(0));
      GNumber exp = GiNaC::ex_to<GiNaC::numeric>(e.op(1));
      ris *= red_prod(1, 1, base, exp);
    }
    else
      ris *= pow(simplify_roots(e.op(0)), simplify_roots(e.op(1)));
  }
  else
    ris += e;
  return ris;
}
