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
#define NOISY 1
#endif

#include "zeilberger.hh"
#include "factorize_giac.hh"
#include "globals.hh"
#include "simplify.hh"
#include "util.hh"
#include "Recurrence.defs.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"

#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

/*!
  By definition, an expression \f$ F(m, k) \f$ is an hypergeometric term
  in both arguments if \f$ F(m+1, k) / F(m, k) \f$ and
  \f$ F(m, k+1) / F(m, k) \f$ are both rational functions of \f$ m \f$
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
  DD_MSGVAR("Zeilberger step one -> ", F_m_k);
  Expr first
    = simplify_binomials_factorials_exponentials(F_m_k.substitute(m, m+1))
    * pwr(simplify_binomials_factorials_exponentials(F_m_k), -1);
  Expr second
    = simplify_binomials_factorials_exponentials(F_m_k.substitute(k, k+1))
    * pwr(simplify_binomials_factorials_exponentials(F_m_k), -1);
  // FIXME: we must decide on the best simplification to use in this case.
  first = simplify_numer_denom(first);
  second = simplify_numer_denom(second);
  if (!first.is_rational_function(m) || !second.is_rational_function(k))
    return false;

  // F_m_k is an hypergeometric term in both arguments and we can thus
  // apply Zeilberger's algorithm.
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

  // FIXME: Assert that r_1,r_2, s_1, s_2 are polynomials.

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

/*!
  Reduce the rational function \f$ q(k) \f$ to canonical form as
  follows.
  Every rational function \f$ q(k) \f$ can be written in the form
  \f$ q(k) = p_1(k+1)/p_1(k) * p_2(k)/p_3(k) \f$ where 
  \f$ p_1, p_2, p_3 \f$ are polynomials in \f$ k \f$ and
  \f$ gdc(p_2(n), p_3(n+h)) = 1 \f$ for all nonnegative integers \f$ k \f$.
*/
bool
parametric_gosper_step_two(const Symbol& k, const Expr& q_k,
			   Expr& p1_k, Expr& p2_k, Expr& p3_k) {
  // FIXME: Consider all special cases as well.
  Expr num;
  Expr den;
  q_k.numerator_denominator(num, den);

  Expr a_k = num;
  Expr b_k = den;
  Expr c_k = 1;

  for (Number h = 1; h <= 100; ++h) {
    // FIXME: Zeilberger trick.
    if (gcd(num, den).degree(k) > 0) {
      Expr temp_b_k = (b_k.substitute(k, k + h)).expand();
      Expr s = gcd(a_k, temp_b_k);
      a_k = quo(a_k, s, k);
      Expr temp_s = s.substitute(k, k - h);
      b_k = quo(b_k, temp_s, k);
      // Using this trick integer_roots[i] is a Number and not an Expr.
      for (Number j = 1; j <= h; ++j)
	c_k *= s.substitute(k, k - j);  
   }
  }

#if 0
  Symbol h("h");
  Expr shifted_den = den.substitute(k, k+h);
  Expr res;
  res = sylvester_matrix_resultant(num, shifted_den, k);
  D_VAR(res);

  // Factorize resultant in order to get integer roots.
  // Not needed for the time being.
  Expr numeric;
  Expr non_numeric;
  factorize_giac(res, numeric, non_numeric);
  //  numeric = 1;
  //  non_numeric = h+2;
  D_VAR(numeric);
  D_VAR(non_numeric);
#endif

#if 0
  integer_roots.push_back(-Recurrence::n-2);
  Expr a_k;
  Expr b_k;

  // normalize f and g, and store conversion factor in Z
  Expr lead_num = num.lcoeff(k);
  Expr lead_den = den.lcoeff(k);
  Expr Z = lead_num * pwr(lead_den, -1);
  a_k = num * pwr(lead_num, -1);
  b_k = den * pwr(lead_den, -1);
  // Computation of the output polynomials.
  c_k = 1;
  unsigned int integer_roots_size = integer_roots.size();
  for (unsigned int i = 0; i < integer_roots_size; ++i) {
    Expr temp_b_k = (b_k.substitute(k, k + integer_roots[i])).expand();
    Expr s = general_gcd(a_k, temp_b_k, k);
    a_k = quo(a_k, s, k);
    Expr temp_s = s.substitute(k, k - integer_roots[i]);
    b_k = quo(b_k, temp_s, k);
    // FIXME: integer_roots[i] is an Expr, not a Number!
    //    for (Number j = 1; j <= integer_roots[i]; ++j)
    //     c_k *= s.substitute(k, k - j);
  }
  a_k *= Z;
  // The polynomials `a_m' and `b_m' may have rational coefficients.
  // Multiply numerator and denominator of the fraction `a_m/b_m'
  // by suitable integers, so that the output values of `a_m' and `b_m'
  // have integer coefficients. 

  if (!a_m.is_integer_polynomial(m)) {
    Expr a_m_factor;
    a_m = convert_to_integer_polynomial(a_m, m, a_m_factor);
    a_m *= a_m_factor.numerator();
    b_m *= a_m_factor.denominator();
  } 

  if (!b_m.is_integer_polynomial(m)) {
    Expr b_m_factor;
    b_m = convert_to_integer_polynomial(b_m, m, b_m_factor);
    a_m *= b_m_factor.denominator();
    b_m *= b_m_factor.numerator();
  }
#endif

  D_VAR(a_k);
  D_VAR(b_k);
  D_VAR(c_k);

  p2_k = a_k;
  p3_k = b_k;
  p1_k = c_k;

  return true;
}


// Coefficients represent the originary unknown related to the order J of
// sought recurrence.
bool
parametric_gosper_step_three(const Symbol& m, const std::vector<Symbol>& coefficients,
			     const Expr& p2_m, const Expr& p3_m,
			     const Expr& p_m, Expr& b_m, std::vector<Expr>& coefficients_values) {
  std::cout << "Zeilberger step 3\n";
  D_VAR(m);
  D_VAR(p2_m);
  D_VAR(p3_m);
  D_VAR(p_m);


  // Gosper's algorithm, step 3.1.
  // Finds the degree of `x(n)'.
  unsigned int deg_a = p2_m.degree(m);
  unsigned int deg_b = p3_m.degree(m);
  unsigned int deg_c = p_m.degree(m);
  Expr lead_a = p2_m.lcoeff(m);
  Expr lead_b = p3_m.lcoeff(m);
  Number deg_x = -1;
  // On output, if a possible degree exists,
  // `deg_x' will contain a non-negative number.
  // FIXME: lead_a and lead_b might contain parameters.
  if (deg_a != deg_b || lead_a != lead_b) {
    if (deg_a >= deg_b && deg_c >= deg_a)
      deg_x = deg_c - deg_a;
    if (deg_b > deg_a && deg_c >= deg_b)
      deg_x = deg_c - deg_b;
  }
  else {
    // `deg_a = deg_b' and `lead_a = lead_b'. [FIXME: consider parameters]
    D_VAR(deg_a);
    Expr A = p2_m.coeff(m, deg_a - 1);
    Expr B = p3_m.substitute(m, m - 1).expand().coeff(m, deg_a - 1);
    // FIXME; B_A is not a Number. Only nonnegative integers are to be
    // considered (A=B, 5.4).
   //    Number B_A = ((B - A) * pwr(lead_a, -1)).ex_to_number();
    Number possible_deg = Number(deg_c) - Number(deg_a) + 1;
    //    if (B_A.is_nonnegative_integer())
    //      if (B_A > possible_deg)
    //	possible_deg = B_A;
    deg_x = possible_deg >= 0 ? possible_deg : -1;
  }
  if (deg_x == -1)
    return false;

  DD_MSGVAR("Looking for a polynomial solution of degree (at most) ", deg_x);

  //Number deg_x = 2;
  //unsigned int deg_a = p2_m.degree(m);
  //unsigned int deg_b = p3_m.degree(m);
  //unsigned int deg_c = p_m.degree(m);
  // `number_of_coeffs' is the number of coefficients of 
  // the polynomial `p', that is, 1 + deg_x.
  unsigned int number_of_coeffs = (1 + deg_x).to_unsigned_int();

  // FIXME: What?
  // Compute the real number of unknowns of the polynomial equation
  // `a(m) * p(m+1) - b(m-1) * p(m) - c(m) = 0'.
  // In general, this is larger than `number_of_coeffs' because 
  // the polynomial equation above may have large degree.
  unsigned int number_of_unknowns = number_of_coeffs;
  number_of_unknowns += deg_a > deg_b ? deg_a : deg_b;
  number_of_unknowns = number_of_unknowns > deg_c ? number_of_unknowns : deg_c;

  // FIXME: Calculate upper bound for degree.
  //  unsigned int upper_bound_degree_b = 10;

  //  unsigned int number_of_unknowns = 10;
  //  unsigned int number_of_coeffs = 9;
  //  unsigned int number_of_unknowns = number_of_coeffs;
  //  number_of_unknowns += deg_a > deg_b ? deg_a : deg_b;
  //  number_of_unknowns = number_of_unknowns > deg_c ? number_of_unknowns : deg_c;

  //   Expr_List unknowns;
  Expr_List unknowns;
  for (unsigned int i = 0; i < number_of_unknowns; ++i) {
    static char buf[2];
    buf[0] = 'A'+i;
    buf[1] = '\0';
    unknowns.append(Symbol(buf));
    D_VAR(unknowns.op(unknowns.nops()-1));
  }
  // Add our originary unknowns to the coefficients. These MUST be added last.
  for (unsigned int i = 0; i < coefficients.size(); ++i) {
    unknowns.append(coefficients[i]);
    number_of_unknowns++;
    D_VAR(unknowns.op(unknowns.nops()-1));
  }

  // Builds the generic polynomial `p' of degree `deg_x'.
  b_m = 0;
  for (unsigned int i = 0; i < number_of_coeffs; ++i)
    b_m += pwr(m, i) * unknowns.op(i);

  D_VAR(p2_m);
  D_VAR(b_m);
  D_VAR(p_m);
  Expr b_m_forward = b_m.substitute(m, m+1).expand();
  D_VAR(b_m_forward);
  Expr p3_back = p3_m.substitute(m, m-1).expand();
  D_VAR(p3_back);

  // Considers the recurrence relation to solve.
  Expr rr = p2_m * b_m_forward - p3_back * b_m - p_m;
  rr = rr.expand();
  D_VAR(rr);
  D_VAR(rr.coeff(m,0));
  D_VAR(rr.coeff(m,1));
  D_VAR(rr.coeff(m,2));
  D_VAR(rr.coeff(m,3));

  // Builds the lists to put in the matrix `rr_coefficients' and `rhs'.
  Expr_List equations;
  for (unsigned int i = 0; i < number_of_unknowns; ++i) {
    Expr lhs = rr.coeff(m, i);
    // D_VAR(lhs);
    equations.prepend(Expr(lhs, 0));
    D_VAR(equations.op(0));
  }
    
  //Expr solution;
  Expr solution = lsolve(equations, unknowns);  
  if (solution.nops() == 0) {
    DD_MSGVAR("Solution not found: ", solution);
    return false;
  }
  else {
    DD_MSGVAR("Polynomial solution found: ", solution);
  }

  std::vector<unsigned int> dummy_vars(number_of_unknowns);
  // Replace solutions in the form "A==A" by putting A=1.
  for (unsigned int j = 0; j < solution.nops(); ++j)
    if (solution.op(j).op(0) == solution.op(j).op(1))
      for (unsigned int i = 0; i < number_of_unknowns; ++i)
	if (solution.op(j).op(1) == unknowns.op(i))
	  dummy_vars.push_back(i);

  // Builds the solution `x(n)'.
  for (unsigned int i = 0; i < number_of_coeffs; ++i)
    b_m = b_m.substitute(unknowns.op(i), solution.op(i).op(1));

  for (unsigned int i = 0; i < dummy_vars.size(); ++i) {
    b_m = b_m.substitute(unknowns.op(dummy_vars[i]), 1);
  }

  D_VAR(b_m);

  // If the polynomial solution is zero, it is invalid and we must try 
  // a recurrence of higher order.
  if (b_m.is_zero())
    return false;

  // Give unknown coefficients the values just found.
  //  for (unsigned int i = 0; i < number_of_coeffs; ++i)
  for (unsigned int j = 0; j < coefficients.size(); ++j) {
    Expr coeff_expr = solution.op(number_of_unknowns - coefficients.size() + j).op(1);
    D_VAR(coeff_expr);
    for (unsigned int i = 0; i < dummy_vars.size(); ++i)
      coeff_expr = coeff_expr.substitute(unknowns.op(dummy_vars[i]), 1);
    coefficients_values[j]=coeff_expr;
  }

  return true;
}


} // anonymous namespace

bool zeilberger_for_fixed_order(const Expr& F_m_k,
				const Symbol& m, const Symbol& k,
				const index_type order,
				Expr& solution) {
  Expr p_0_k = 0;
  Expr r_k = 0;
  Expr s_k = 0;

  std::vector<Symbol> coefficients(order + 1);
  std::vector<Expr> coefficients_values(order + 1);
#if 1
  // FIXME: Temporary, to use familiar names when order==1.
  if (order == 1) {
    coefficients[0] = Symbol("F");
    coefficients[1] = Symbol("G");
  }
  else
#endif
    for (index_type i = 0; i < order + 1; ++i)
      coefficients[i] = Symbol();
  if (!zeilberger_step_one(F_m_k, m, k, coefficients, p_0_k, r_k, s_k))
    return false;
  DD_VAR(p_0_k);
  DD_VAR(r_k);
  DD_VAR(s_k);

  Expr p_1_k;
  Expr p_2_k;
  Expr p_3_k;
#if 1
  if (!parametric_gosper_step_two(k, r_k / s_k, p_1_k, p_2_k, p_3_k))
    // Problem in the computation of the resultant and its roots.
    return false;
#else
  const Expr& tmp = simplify_all(r_k / s_k);
  p_1_k = 1;
  p_2_k = tmp.numerator();
  p_3_k = tmp.denominator();
#endif
  DD_VAR(p_1_k);
  DD_VAR(p_2_k);
  DD_VAR(p_3_k);

  Expr p_k = p_0_k * p_1_k;
  Expr b_k;
  if (parametric_gosper_step_three(k, coefficients, p_2_k.expand(), p_3_k.expand(),
				   p_k.expand(), b_k, coefficients_values)) {
    DD_MSGVAR("Polynomial solution is: ", b_k);
    for (unsigned int i = 0; i < coefficients_values.size(); ++i) {
      D_VAR(coefficients_values[i]);
    }
  }
  else 
    // Could not solve the linear system. Maybe we should just try an higher order.
    return false;

  for (unsigned int i = 0; i < coefficients.size(); ++i) {
    p_2_k = p_2_k.substitute(coefficients[i], coefficients_values[i]);
    p_3_k = p_2_k.substitute(coefficients[i], coefficients_values[i]);
    p_k = p_2_k.substitute(coefficients[i], coefficients_values[i]);
  }

  // Construct the telescoped recurrence.
  const Symbol& n = Recurrence::n;
  Expr rec_expr = 0;
  for (unsigned int i = 1; i < order + 1; ++i) {
    rec_expr += coefficients_values[i] * x(n+i);
    D_VAR(rec_expr);
  }
  assert(!coefficients_values[0].is_zero());
  rec_expr = - rec_expr / coefficients_values[0];

  // Expr rec_expr =  - (x(n+1) * coefficients_values[1]) / coefficients_values[0];
  Recurrence rec(rec_expr);
  Expr exact_solution;
  std::map<index_type, Expr> initial_conditions;
  rec.compute_exact_solution();
  // FIXME: Are exactly N initial conditions needed for a recurrence of order N?
  // Revise this when number_initial_conditions() is written.
  for (unsigned int i = 0; i < order; ++i) {
    unsigned int evaluation_point = i + rec.first_valid_initial_condition();
    D_VAR(evaluation_point);
    // FIXME: The support of F_m_k for m==evaluation_point must be in [0, max_support].
    // Calculate it explicitly if possible.
    Expr max_support = 10;
    initial_conditions[evaluation_point] = sum( (Expr) k, (Expr) 0, max_support, 
						F_m_k.substitute(m, (Number) evaluation_point));
  }
  rec.set_initial_conditions(initial_conditions);
  rec.exact_solution(exact_solution);
  exact_solution = simplify_all(exact_solution);
  D_VAR(exact_solution);

  solution = exact_solution;

  return true;
}


/*!
  \param F_m_k    proper hypergeometric sequence.
  \param m        first variable of \p F_m_k.
  \param k        second variable of \p F_m_k.
  \param solution Variable where the solution will be stored.
  \return         ...
*/
bool
PURRS::zeilberger_algorithm(const Expr& F_m_k,
			    const Symbol& m, const Symbol& k, Expr& solution) {
  // FIXME: Calculate the highest possible order for the target recurrence.
  index_type highest_possible_order = 10;

  unsigned int i = 1;

  for (i = 1; i < highest_possible_order; ++i) {
    DD_MSGVAR("Zeilberger: trying order ", i);
    if (zeilberger_for_fixed_order(F_m_k, m, k, i, solution)) {
      break;
    }
    // No further attempts can be made.
    if (i == highest_possible_order - 1)
      return false;
  }

  DD_MSGVAR("Zeilberger: found solution for order ", i);
  D_VAR(solution);
  return true;
}
       
