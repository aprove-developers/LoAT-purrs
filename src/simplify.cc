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

static const unsigned
FACTOR_THRESHOLD = 100;

static void
split_exponent(GExpr& num, GExpr& not_num, const GExpr& e) {
  if (!is_a<numeric>(e))
    not_num *= e;
  else
    num *= e;
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
  GExpr f = e;
  // At the end of the following cycle 'f' will contains the base of
  // the power.
  while (is_a<power>(f)) {
    // Checks and divides the exponents in two part: those numeric and  
    // those not numeric.
    if (is_a<mul>(f.op(1)))
      for (unsigned i = f.op(1).nops(); i-- > 0; )
        split_exponent(num_exponent, not_num_exponent, f.op(1).op(i));
    else
      split_exponent(num_exponent, not_num_exponent, f.op(1));
    f = f.op(0);
  };
#if NOISY
  std::cout << "base " << f << endl;
  std::cout << "num_exp " << num_exponent << endl;
  std::cout << "not_num_exp " << not_num_exponent << endl;
#endif
  // The base is a multiplication.
  if (is_a<mul>(f)) {
    GExpr tot = 1;
    bool base_power = false;
    for (unsigned i = f.nops(); i-- > 0; )
      if (is_a<power>(f.op(i)))
	base_power = true;
    // Some factor of the base is a power.
    if (base_power)
      tot = simpl_powers_base(f, num_exponent, not_num_exponent);
    // No factors of the base is a power.
    else
      for (unsigned i = f.nops(); i-- > 0; )
	if (!is_a<numeric>(f.op(i)))
	  tot *= pow(f.op(i), num_exponent * not_num_exponent);
	else
          tot *= pow(pow(f.op(i), num_exponent), not_num_exponent);
    return tot;
  }
  // The base is not a multiplication.
  else
    if (!is_a<numeric>(f))
      return pow(f, num_exponent * not_num_exponent);
    else
      return pow(pow(f, num_exponent), not_num_exponent);
  return 0;
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
    for (unsigned i = e.nops(); i-- > 0; ) {
#if NOISY
      std::cout << "e " << e.op(i) << endl;
#endif
      ris += simplify_on_input_ex(e.op(i));
    }
  }
  else if (is_a<mul>(e)) {
    ris = 1;
    for (unsigned i = e.nops(); i-- > 0; ) {
#if NOISY
      std::cout << "e " << e.op(i) << endl;
#endif
      ris *= simplify_on_input_ex(e.op(i));
    }
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

// Rule 5 and 6.
static GExpr
collect_same_base(const GExpr& e) {
  GExpr tmp = e;
  std::cout << "Rule 5 o 6" << std::endl;
  // Rule 5.
  while (tmp.has(wild(1) * pow(wild(1), wild(2))))
    tmp = tmp.subs(wild(1) * pow(wild(1), wild(2))
		   == pow(wild(1), 1 + wild(2)));
  while (tmp.has(wild(0) * wild(1) * pow(wild(1), wild(2)))) {
    //std::cout << "5.2" << endl;
    tmp = tmp.subs(wild(0) * wild(1) * pow(wild(1), wild(2))
                   == wild(0) * pow(wild(1), 1 + wild(2)));
  }
  // Rule 6.
  while (tmp.has(wild(0) * pow(wild(1), wild(2)) * pow(wild(1), wild(3))))
    tmp = tmp.subs(wild(0) * pow(wild(1), wild(2)) * pow(wild(1), wild(3))
                   == wild(0) * pow(wild(1), wild(2) + wild(3)));
  return tmp;
}

// Rule 12.
// Applied under condition that the exponent is not an integer.
static GExpr
collect_same_exponent(const GExpr& e) {
  GExpr tmp = e;
  GList lst;
  bool again = true;
  while (again)
    if (tmp.find(pow(wild(1), wild(0)) * pow(wild(2), wild(0)) * wild(3),
		 lst)) {
      if (is_a<numeric>(lst.op(0).op(0).op(1))) {
	GNumber num = ex_to<numeric>(lst.op(0).op(0).op(1));
	if (!num.is_integer())
	  tmp = tmp.subs(pow(wild(1), wild(0))
			 * pow(wild(2), wild(0)) * wild(3)
			 == pow(wild(1)*wild(2), wild(0)) * wild(3));
	else
	  again = false;
      }
      else
	tmp = tmp.subs(pow(wild(1), wild(0)) * pow(wild(2), wild(0)) * wild(3)
		       == pow(wild(1)*wild(2), wild(0)) * wild(3));
    }
    else
      again = false;
//  while (tmp.has(pow(wild(1), wild(0)) * pow(wild(2), wild(0)) * wild(3))) {
//     std::cout << "tmp prima: " << tmp << endl;
//     tmp = tmp.subs(pow(wild(1), wild(0)) * pow(wild(2), wild(0)) * wild(3)
//                    == pow(wild(1)*wild(2), wild(0)) * wild(3));
//     std::cout << "tmp dopo: " << tmp << endl;
//   } 
//   std::cout << "ris_12: " << tmp << endl;
  return tmp;
}

/*!
  Crosses the tree of the expression expanded \f$ e \f$ recursevely to find
  subexpressions which to apply the rule 5, 6 and 12 of the terms rewriting
  system \f$ /mathfrak{R} \f$ (the rule 11 is automatically applied by GiNaC).
  Returns a <CODE>GExpr</CODE> that contains the input expression \f$ e \f$
  modified.  
*/ 
GExpr
simplify_on_output_ex(const GExpr& e) {
  GExpr ris;
  if (is_a<add>(e)) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; ) {
#if NOISY
      std::cout << "e " << e.op(i) << endl;
#endif
      ris += simplify_on_output_ex(e.op(i));
    }
  }
  else if (is_a<mul>(e)) {
    GList lst;
    // Rule 5 and 6.
    if (e.has(wild(1) * pow(wild(1), wild(2)))
	// FIXME: questa seconda condizione non prende casi come
	// 5*x*x^a mentre prende b*x*x^a.
	|| e.has(wild(0) * wild(1) * pow(wild(1), wild(2)))
	|| e.has(wild(0) * pow(wild(1), wild(2)) * pow(wild(1), wild(3))))
      ris += collect_same_base(e);

    // Rule 12.
    else if (e.has(pow(wild(1), wild(0)) * pow(wild(2), wild(0)) * wild(3)))
      ris += collect_same_exponent(e);
    else
      ris += e;
  }
  else
    ris += e;
  return ris;
}      

static GExpr
get_binding(const GList& l, unsigned wild_index) {
  assert(wild_index < l.nops());
  assert(l.op(wild_index).info(GiNaC::info_flags::relation_equal));
  assert(l.op(wild_index).lhs() == GiNaC::wild(wild_index));
  return l.op(wild_index).rhs();
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
  while(irem(m, 2) == 0) { // the case 2 is handled separately 
    m /= 2;
    ++k;
  }
  if (k>0) {
    bases.push_back(2);
    exponents.push_back(k);
  }
  for (unsigned i=3; (i < FACTOR_THRESHOLD) && (i*i <= m); i += 2) {
    k = 0;
    while(irem(m, i) == 0) { // test for divisibility by the odd integer i
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
red_prod(const GExpr& base_1, const GExpr& exp_1, 
	 const GExpr& base_2, const GExpr& exp_2) {
  
  GExpr base1 = base_1;
  GExpr base2 = base_2;
  GExpr exp1 = exp_1;
  GExpr exp2 = exp_2;
  
  GExpr q = 0;
  if (GiNaC::is_a<GiNaC::numeric>(base1) && 
      GiNaC::is_a<GiNaC::numeric>(base2) &&
      GiNaC::is_a<GiNaC::numeric>(exp1)  &&
      GiNaC::is_a<GiNaC::numeric>(exp2)) {
    assert(exp1 != 0);
    assert(exp2 != 0);
    
    GExpr   n_d_1  = numer_denom(exp1);
    GNumber k1_num = GiNaC::ex_to<GiNaC::numeric>(n_d_1.op(0));
    GNumber k1_den = GiNaC::ex_to<GiNaC::numeric>(n_d_1.op(1));
    GExpr   n_d_2  = numer_denom(exp2);
    GNumber k2_num = GiNaC::ex_to<GiNaC::numeric>(n_d_2.op(0));
    GNumber k2_den = GiNaC::ex_to<GiNaC::numeric>(n_d_2.op(1));
    
    base1 = pow(base1, k1_num);
    base2 = pow(base2, k2_num);
    
    GNumber g  = gcd(k1_den, k2_den);
    GNumber k  = k1_den * k2_den / g;
    GNumber b1 = GiNaC::ex_to<GiNaC::numeric>(pow(base1, k2_den / g));
    GNumber b2 = GiNaC::ex_to<GiNaC::numeric>(pow(base2, k1_den / g));
    GNumber b = b1 * b2;
    q = reduce_to_standard_form(k, b);
  }
  else q = pow(base1, exp1) * pow(base2, exp2);
  return q;
}

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
    tmp = factor * red_prod(base_1, exp_1, base_2, exp_2);
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

 // These routines simplify expressions that GiNaC leaves alone: 
 // for instance, they simplify 9^(1/4)*sqrt(3) to 3, which is correct, 
 // while GiNaC leaves it unchanged, and evalf(9^(1/4)*sqrt(3)) produces 
 // an approximate answer which is slightly too large. 
 // Actually, evalf(9^(1/4)*sqrt(3) - 3) yields the answer 
 // 2.168404344971008868E-19


GExpr
simplify_roots(const GExpr& e){
  GExpr ris;
  if (is_a<add>(e)) {
    ris = 0;
    for (unsigned i = e.nops(); i-- > 0; ) {
#if NOISY
      cout << "e(add) " << e.op(i) << endl;
#endif
      ris += simplify_roots(e.op(i));
    }
  }
  else if (is_a<function>(e)) {
    ris = 1;
    GExpr f = e;
    GExpr tmp = simplify_roots(e.op(0));
    ris *= f.subs(f.op(0) == tmp);
  }
  else if (is_a<mul>(e))
    return reduce_product(e);
  else if (is_a<power>(e)) {
    ris = 1;
    ris *= pow(simplify_roots(e.op(0)), simplify_roots(e.op(1)));
  }
  else
    ris += e;
  return ris;
}
