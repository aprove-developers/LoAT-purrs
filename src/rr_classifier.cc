/* Definition of the main recurrence relation solver.
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

#include <config.h>

#ifndef NOISY
#define NOISY 0
#endif

#include "util.hh"
#include "simplify.hh"
#include "numerator_denominator.hh"
#include "finite_order.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include "Finite_Order_Info.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Non_Linear_Info.defs.hh"
#include "Recurrence.defs.hh"
#include "Recurrence.inlines.hh"

#include <climits>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>

// TEMPORARY
#include <iostream>
#include <fstream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

#define Napier exp(Expr(1))

namespace {
using namespace PURRS;

/*!
  Returns <CODE>true</CODE> if \p e is of the form \f$ n - d \f$ with
  \f$ d \f$ an integer different from \f$ 0 \f$: in this case assign
  the opposite of \f$ d \f$ to \p decrement.
  Returns <CODE>false</CODE> otherwise.
*/
bool
get_constant_decrement(const Expr& e, Number& decrement) {
  if (e.is_a_add() && e.nops() == 2) {
    // `e' is of the form a+b.
    const Expr& a = e.op(0);
    const Expr& b = e.op(1);
    Expr d;
    if (a == Recurrence::n)
      d = b;
    else if (b == Recurrence::n)
      d = a;
    else
      return false;
    Number i;
    if (d.is_a_number(i) && i.is_integer()) {
      decrement = -i;
      return true;
    }
  }
  return false;
}

bool
has_constant_decrement(const Expr& e) {
  if (e.is_a_add() && e.nops() == 2) {
    // `e' is of the form a+b.
    const Expr& a = e.op(0);
    const Expr& b = e.op(1); 
    Expr d;
    if (a == Recurrence::n)
      d = b;
    else if (b == Recurrence::n)
      d = a;
    else
      return false;
    Number i;
    if (d.is_a_number(i) && i.is_integer())
      return true;
  }
  return false;
}

/*!
  Returns <CODE>true</CODE> if \p e is of the form \f$ n / d \f$ with
  \f$ d \f$ a positive rational: in this case assign the inverse of
  \f$ d \f$ to \p divisor.
  Returns <CODE>false</CODE> otherwise.
*/
bool
get_constant_divisor(const Expr& e, Number& divisor) {
  assert(e.is_a_mul() && e.nops() == 2);
  // `e' is of the form a*b.
  const Expr& a = e.op(0);
  const Expr& b = e.op(1);
  Number d;
  if (a == Recurrence::n && b.is_a_number(d) && d.is_rational()
      && d.is_positive()) {
    divisor = 1 / d;
    return true;
  }
  else if (b == Recurrence::n && a.is_a_number(d) && d.is_rational()
	   && d.is_positive()) {
    divisor = 1 / d;
    return true;
  }
  else
    return false;
}

#if 0
void
impose_condition(const std::string&) {
}
#endif

/*!
  Performes a change of variable for the `x' functions in the
  expression \p e substituting every occurrence of the expression
  \p s in the arguments of the `x' functions with the expression \p r.
  Returns an expression containing the result of the substitution.
*/
Expr
change_variable_function_x(const Expr& e, const Expr& s, const Expr& r) {
  Expr e_substituted;
    if (e.is_a_add()) {
    e_substituted = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_substituted += change_variable_function_x(e.op(i), s, r);
  }
  else if (e.is_a_mul()) {
    e_substituted = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_substituted *= change_variable_function_x(e.op(i), s, r);
  }
  else if (e.is_a_power())
    return pwr(change_variable_function_x(e.arg(0), s, r),
	       change_variable_function_x(e.arg(1), s, r));
  else if (e.is_a_function())
    if (e.nops() == 1)
      if (e.is_the_x_function())
	return apply(e.functor(), e.arg(0).substitute(s, r));
      else
	return apply(e.functor(),
		     change_variable_function_x(e.arg(0), s, r));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int j = 0; j < num_argument; ++j)
	argument[j] = change_variable_function_x(e.arg(j), s, r);
      return apply(e.functor(), argument);
    }
    else
      return e;
    return e_substituted;
}

/*!
  If \p possibly_dec is greater than \p max_decrement then \p possibly_dec
  becomes the new maximum decrement and it is assigned to \p max_decrement,
  in this case \p possibly_coeff becomes the new \p coefficient.
*/
void
assign_max_decrement_and_coeff(const Expr& possibly_dec,
			       const Expr& possibly_coeff,
			       int& max_decrement, Expr& coefficient) {
  Number decrement;
  get_constant_decrement(possibly_dec, decrement);
  int dec = -decrement.to_int();
  if (dec > max_decrement) {
    max_decrement = dec;
    coefficient = possibly_coeff;
  }
}

void
find_max_decrement_and_coeff_factor(const Expr& e,
				    int& max_decrement, Expr& coefficient) {
  assert(!e.is_a_add());
  if (e.is_a_mul()) {
    Expr possibly_coeff = 1;
    Expr possibly_argument = 0;
    for (unsigned int i = e.nops(); i-- > 0; ) {
      const Expr& factor = e.op(i);
      if (factor.is_the_x_function())
	possibly_argument = factor.arg(0);
      else
	possibly_coeff *= factor;
      if (!possibly_argument.is_zero())
	assign_max_decrement_and_coeff(possibly_argument, possibly_coeff,
				       max_decrement, coefficient);
    }
  }
  else
    if (e.is_the_x_function())
      assign_max_decrement_and_coeff(e.arg(0), 1,
				     max_decrement, coefficient);
}

/*!
  Let \p e be the right hand side of a linear recurrence.
  This functions seeks the largest positive integer \f$ j \f$ such that
  \f$ x(n+j) \f$ occurs in \p e and its coefficient.
  These two values are stored in \p max_decrement and \p coefficient,
  respectively.
*/
void
find_max_decrement_and_coeff(const Expr& e,
			     int& max_decrement, Expr& coefficient) {
  if (e.is_a_add())
    for (unsigned int i = e.nops(); i-- > 0; )
      find_max_decrement_and_coeff_factor(e.op(i),
					  max_decrement, coefficient);
  else
    find_max_decrement_and_coeff_factor(e, max_decrement, coefficient);
}

/*!
  Assuming that \p rhs contains occurrences of \f$ x(n-k) \f$
  where \f$ k \f$ is a negative integer, this function
  performs suitable changes of variables that preserve the meaning of
  the recurrence relation, but transforms it into its <EM>standard
  form</EM> \f$ x(n) = new_rhs \f$, where \f$ new_rhs \f$
  does not contain any instance of \f$ x(n-k) \f$, with a
  negative integer \f$ k \f$.
*/
void
eliminate_negative_decrements(const Expr& rhs, Expr& new_rhs) {
  assert(rhs.is_expanded());
  // Seeks `max_decrement', i.e., the largest positive integer `j' such that
  // `x(n+j)' occurs in `rhs' with a coefficient `coefficient' which is not
  // syntactically 0.
  int max_decrement = INT_MIN;
  Expr coefficient;
  find_max_decrement_and_coeff(rhs, max_decrement, coefficient);
  // The changes of variables includes replacing `n' by `n-max_decrement',
  // changing sign, and division by `coefficient'.
  new_rhs = change_variable_function_x(rhs, Recurrence::n,
				       Recurrence::n - max_decrement);
  new_rhs *= -1;
  new_rhs = new_rhs.substitute(x(Recurrence::n),
			       - x(Recurrence::n-max_decrement)
			 * pwr(coefficient, -1));
  new_rhs /= coefficient;
}

/*!
  Let \f$ e = a x(n) + b \f$ be the expression contained in \p e.
  Initially, \p coeff_x_n is equal to \f$ 0 \f$ and \p remainder is
  equal to \p e.
  This function finds \f$ a \f$ and \f$ b \f$ and stores them in
  \p coeff_x_n and \p remainder, respectively.
*/
void
find_coeff_x_n_and_remainder(const Expr& e,
			     Expr& coeff_x_n, Expr& remainder) {
  for (unsigned int i = e.nops(); i-- > 0; ) {
    const Expr& term = e.op(i);
    if (term.is_a_mul()) {
      Expr tmp = 1;
      bool found_x_n = false;
      for (unsigned int j = term.nops(); j-- > 0; ) {
	const Expr& factor = term.op(j);
	if (factor == x(Recurrence::n))
	  found_x_n = true;
	else
	  tmp *= factor;
      }
      if (found_x_n) {
	coeff_x_n = tmp;
	remainder -= term;
      }
    }
    else
      if (term == x(Recurrence::n)) {
	coeff_x_n = 1;
	remainder -= term;
      }
  }
}

/*!
  Here we assume that \p rhs contains occurrences of \f$ x(n) \f$ itself.
  Therefore the recurrence may be impossible or indeterminated.
  This function decides if this is the case and, if so, it returns
  \f$ 1 \f$ or \f$ 2 \f$, respectively.
  If the recurrence is solvable, it is rewritten into its normal form,
  which is then written in \f$ new_rhs \f$, and the function returns
  \f$ 0 \f$.
  Returns \f$ 3 \f$ in all the other cases for which the recurrence is
  considered too complex.
*/
unsigned int
eliminate_null_decrements(const Expr& rhs, Expr& new_rhs) {
  assert(rhs.is_expanded());
  // Collect the terms `x(n)' so that the right hand side of the recurrence
  // `rhs' is in the form `rhs = a*x(n) + b' and that `b' does different to
  // zero and does not contain `x(n)'.
  new_rhs = rhs.collect(Expr_List(x(Recurrence::n)));

  // The following cases are possible:
  // 1. If `a = 1' and `b = 0' then the recurrence is indeterminate. 
  // 2. If `a = 1' and `b != 0' does not contain any occurrence of `x(n-k)'
  //    where `k' is a positive integer, the recurrence is impossible.
  // 3. If `a = 1' and `b' contains `x(n-k)' for some positive integer `k'
  //    and with a coefficient that is not syntactically 0, we remove
  //    `x(n)' from both sides of `x(n) = rhs', and then rewrite the
  //    recurrence into its standard form.
  // 4. If `a != 1' we move `a*x(n)' to the left-hand side, and divide
  //    through by `1 - a', obtaining the standard form, which is 
  //    `(rhs - a*x(n)) / (1-a)'.
  if (new_rhs.is_a_add()) {
    // Finds `a' and `b'.
    Expr a = 0;
    Expr b = new_rhs;
    find_coeff_x_n_and_remainder(new_rhs, a, b);
    D_VAR(a);
    D_VAR(b);
    if (a == 1) {
      // Case 2. and Case 3.
      bool found_x = false;
      if (b.is_a_add())
	for (unsigned int i = b.nops(); i-- > 0; ) {
	  const Expr& term = b.op(i);
	  if (term.is_a_mul()) {
	    for (unsigned int j = term.nops(); j-- > 0; ) {
	      const Expr& t_j = term.op(j);
	      if (t_j.is_the_x_function()) {
		const Expr& t_j_arg = t_j.arg(0);
		if (t_j_arg == Recurrence::n 
		    || has_constant_decrement(t_j_arg)) {
		  found_x = true;
		  break;
		}
		else
		  return 3;
	      }
	    }
	  }
	  else
	    if (term.is_the_x_function()) {
	      const Expr& arg = term.arg(0);
	      if (arg == Recurrence::n  || has_constant_decrement(arg)) {
		found_x = true;
		break;
	      }
	      else
		return 3;
	    }
	}
      else if (b.is_a_mul())
	for (unsigned int i = b.nops(); i-- > 0; ) {
	  const Expr& b_i = b.op(i);
	  if (b_i.is_the_x_function()) {
	    const Expr& b_i_arg = b_i.arg(0);
	    if (b_i_arg == Recurrence::n 
		|| has_constant_decrement(b_i_arg)) {
	      found_x = true;
	      break;
	    }
	    else
	      return 3;
	  }
	}
      else if (b.is_the_x_function()) {
	const Expr& arg = b.arg(0);
	if (arg == Recurrence::n  || has_constant_decrement(arg))
	  found_x = true;
	else
	  return 3;
      }
      // Case 2.
      if (!found_x)
	return 1;
      // Case 3.
      else {
	// Seeks `max_decrement', i.e., the largest integer `j' (it may be
	// non positive) such that `x(n+j)' occurs in `b' with a coefficient
	// `coefficient' which is not syntactically 0.
	int max_decrement = INT_MIN;
	Expr coefficient;
	find_max_decrement_and_coeff(b, max_decrement, coefficient);
	// Rewrites the recurrence into its standard form:
	// removes from `b' the term that will be the right hand side of the
	// recurrence, i.e. `x(n+max_decrement)'; changes variable replacing
	// `n+max_decrement' by `n', changes sign and divides for the
	// coefficient of `x(n+max_decrement)'.
	new_rhs = b - coefficient * x(Recurrence::n+max_decrement);
	new_rhs = new_rhs.substitute(Recurrence::n,
				     Recurrence::n-max_decrement);
	new_rhs *= -1;
	new_rhs /= coefficient;
      }
    }
    // Case 4.
    else
      new_rhs = b * pwr(1 - a, -1);
  }
  // Let `rhs = a*x(n)'.
  else if (new_rhs.is_a_mul())
    for (unsigned int i = new_rhs.nops(); i-- > 0; ) {
      const Expr& factor = new_rhs.op(i);
      if (factor == x(Recurrence::n)) {
	new_rhs = 0;
	break;
      }
    }
  // Let `rhs = x(n)'.
  else if (new_rhs == x(Recurrence::n))
    return 2;
  return 0;
}
  
void
insert_coefficients(const Expr& coeff, unsigned long index,
		    std::vector<Expr>& coefficients) {
  // The vector `coefficients' contains in the `i'-th position the
  // coefficient of `x(n-i)'.  The first position always contains 0.
  if (index > coefficients.size())
    coefficients.insert(coefficients.end(),
			index - coefficients.size(),
			Number(0));
  if (index == coefficients.size())
    coefficients.push_back(coeff);
  else
    coefficients[index] += coeff;
}

/*!
  - If there is in \p e an \f$ x \f$ function (with the argument containing
    \f$ n \f$) inside an other function different from \f$ x \f$ function
    or if there is in \p e a power with an \f$ x \f$ function in the base
    or in the exponent of it then returns \f$ 0 \f$;
  - if there is in \p e an \f$ x \f$ function (with the argument containing
    \f$ n \f$) inside an other \f$ x \f$ function return \f$ 1 \f$;
  - in all the other cases returns \f$ 2 \f$.
*/
unsigned int
x_function_in_powers_or_functions(const Expr& e) {
  // There is an `x' function (with the argument containing `n') inside an
  // other function.
  if (e.is_a_function()) {
    for (unsigned int i = e.nops(); i-- > 0; ) {
      // `operand' is the argument (if unary function) or the i-th
      // argument (otherwise) of the function `e'.
      const Expr& operand = e.arg(i);
      if (operand.has_x_function(Recurrence::n))
	if (e.is_the_x_function())
	  return 1;
	else
	  return 0;
    }
  }
  // There is a power with an `x' function in the base or in the exponent.
  else if (e.is_a_power()) {
    const Expr& base = e.arg(0);
    const Expr& exponent = e.arg(1);
    if (base.has_x_function(Recurrence::n)
	|| exponent.has_x_function(Recurrence::n))
      return 0;
  }
  return 2;
}


//! \brief
//! Returns \f$ 1 \f$ if finds the non linear term of the form
//! \f$ x(x(a)) \f$ with \f$ a \f$ containing the special symbol \f$ n \f$;   
//! returns \f$ 0 \f$ if finds all the other type of non linear term;
//! returns \f$ 2 \f$ otherwise.
unsigned int
find_non_linear_term(const Expr& e) {
  assert(!e.is_a_add());
  // Even if we find a `legal' (i.e. not containing two nested `x' function)
  // non-linear term we do not exit from this function because we must be
  // sure that there are not non-linear term not legal.
  // If `non_linear_term' at the end of this function is again `2' means
  // that we have not find non-linear terms.
  unsigned int non_linear_term = 2;
  unsigned int num_factors = e.is_a_mul() ? e.nops() : 1;
  if (num_factors == 1) {
    unsigned int tmp = x_function_in_powers_or_functions(e);
    // Nested `x' function: not legal non-linear term.
    if (tmp == 1)
      return 1;
    // We have found legal non-linear term.
    else if (tmp == 0)
      non_linear_term = 0;
  }
  else {
    bool found_function_x = false;
    for (unsigned int j = num_factors; j-- > 0; ) {
      const Expr& factor = e.op(j);
      unsigned int tmp = x_function_in_powers_or_functions(factor);
      // Nested `x' function: not legal non-linear term.
      if (tmp == 1)
	return 1;
      // We have found legal non-linear term.
      else if (tmp == 0)
	non_linear_term = 0;
      if (factor.is_the_x_function())
	// We have found legal non-linear term product of two ore more
	// `x' functions.
	if (found_function_x)
	  non_linear_term = 0;
	else
	  if (factor.arg(0).has(Recurrence::n))
	    found_function_x = true;
    }
  }
  return non_linear_term;
}

//! \brief
//! Returns <CODE>true</CODE> if the non linear recurrence \f$ x(n) = rhs \f$
//! is rewritable as linear recurrence \f$ x(n) = new_rhs \f$.
//! Returns <CODE>false</CODE> otherwise.
/*!
  Cases of rewritable non linear recurrences:
  -  \f$ x(n) = a x(n-k_1)^b_1 \cdots x(n-k_h)^b_h \f$,
     where \f$ k_1, \dots, k_h, b_1, \dots, b_h \f$ are positive integers
     and \f$ h > 1 \f$;
  -  \f$ x(n) = x(n-k)^b \f$,
     where \f$ b \f$ is a rational number while \f$ k \f$ is a
     positive integers.
  The two cases above hold also if instead of terms like `x(n-k)' there are
  term like `x(n/k)'.
  \p coeff_and_base_ is used in two different ways:
  In the case of simple non-linear recurrence of the form
  \f$ x(n) = c x(n-1)^{\alpha} \f$ it contains the pair \f$ c, \alpha \f$;
  in all the other cases the numeric value of the pair holds \f$ 0 \f$
  while the second element contains the value that will be the
  logarithm's base or the exponential's base used in the rewriting
  of the non-linear recurrence in the correspondent linear recurrence.
*/
bool
rewrite_non_linear_recurrence(const Recurrence& rec, const Expr& rhs,
			      Expr& new_rhs,
			      std::pair<Number, Expr>& coeff_and_base,
			      std::vector<Symbol>& auxiliary_symbols) {
  // First case.
  if (rhs.is_a_mul()) {
    bool simple_cases = false;
    Number common_exponent = 1;
    for (unsigned int i = rhs.nops(); i-- > 0; ) {
      const Expr& factor = rhs.op(i);
      Number num_exp;
      if (factor.is_a_power() && factor.arg(0).is_the_x_function()
	  && factor.arg(1).is_a_number(num_exp)) {
	assert(num_exp.is_rational());
	simple_cases = true;
	// If one of the two number is not integer, then `lcm()'
	// returns their product.
	common_exponent = lcm(num_exp, common_exponent);
      }
      else if (factor.is_the_x_function())
	simple_cases = true;
      // Recurrence that the system is not able to transform in linear. 
      else if (factor.has_x_function(Recurrence::n)) {
	simple_cases = false;
	break;
      }
    }
    // Consider the special case `x(n) = c x(n-1)^a' with `a' and `c'
    // constants (`a != 1').
    // In this case is not necessary the `linearization' of the recurrence
    // because we already know the solution in function of `c' and `a'.
    if (simple_cases && rhs.nops() == 2)
      if (rhs.op(0).is_a_number() && rhs.op(1).is_a_power()
	  && rhs.op(1).arg(0).is_the_x_function()) {
	coeff_and_base.first = rhs.op(0).ex_to_number();
	assert(rhs.op(1).arg(1).is_a_number());
	if (rhs.op(1).arg(1).ex_to_number().is_negative())
	  common_exponent *= -1;
	coeff_and_base.second = common_exponent;
	return true;
      }
      else if (rhs.op(1).is_a_number() && rhs.op(0).is_a_power()
	       && rhs.op(0).arg(0).is_the_x_function()) {
	coeff_and_base.first = rhs.op(1).ex_to_number();
	assert(rhs.op(0).arg(1).is_a_number());
	if (rhs.op(0).arg(1).ex_to_number().is_negative())
	  common_exponent *= -1;
	coeff_and_base.second = common_exponent;
	return true;
      }
    new_rhs = 0;
    if (simple_cases)
      if (common_exponent == 1) {
	// Substitute any function `x()' with `exp(1)^{x()}'.
	coeff_and_base.second = Napier;
	Expr tmp = substitute_x_function(rhs, Napier, true);
	tmp = simplify_ex_for_input(tmp, true);
	for (unsigned int i = tmp.nops(); i-- > 0; ) {
	  Number num;
	  if (tmp.op(i).is_a_number(num) && num.is_negative()) {
	    Symbol s = rec.insert_auxiliary_definition(num);
	    auxiliary_symbols.push_back(s);
	    new_rhs += log(s);
	  }
	  else
	    new_rhs += log(tmp.op(i));
	}
	new_rhs = simplify_logarithm(new_rhs);
	
	return true;
      }
      else {
	// Substitute any function `x()' with `common_exponent^{x()}'.
	coeff_and_base.second = common_exponent;
	Expr tmp = substitute_x_function(rhs, abs(common_exponent), true);
	tmp = simplify_ex_for_input(tmp, true);
	for (unsigned int i = tmp.nops(); i-- > 0; ) {
	  Number num;
	  if (tmp.op(i).is_a_number(num) && num.is_negative()) {
	    Symbol s = rec.insert_auxiliary_definition(num);
	    auxiliary_symbols.push_back(s);
	    new_rhs += log(s) / log(abs(common_exponent));
	  }
	  else
	    new_rhs += log(tmp.op(i)) / log(abs(common_exponent));
	}
	new_rhs = simplify_logarithm(new_rhs);
	return true;
      }
  }
  // Second case.
  else if (rhs.is_a_power()) {
    Number num_exp;
    if (rhs.arg(0).is_the_x_function() && rhs.arg(1).is_a_number(num_exp)) {
      coeff_and_base.second = num_exp;
      coeff_and_base.first = 1;
      return true;
    }
    else
      return false;
  }
  return false;
}

//! \brief
//! Returns <CODE>true</CODE> if the linear infinite order recurrence,
//! whose right-hand-side is stored in \p rhs,
//! belongs to the class the system is able to compute, i.e. is in
//! \ref normal_form "normal form" or is possible to rewrite it in
//! normal form \f[ x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n) \f];
//! returns <CODE>false</CODE> otherwise.
/*!
  \param rhs                  the right hand side of a infinite order
                              recurrence.
  \param term_sum             the term of \p rhs containing the sum.  
  \param weight               the expression multiplied for the sum
                              (\f$ f(n) \f$).
  \param inhomogeneous        the inhomogeneous term of the infinite
                              order recurrence (\f$ g(n) \f$).
  \param rhs_first_order      the right hand side of the first order
                              recurrence associated to the infinite
			      order recurrence.
  \param first_valid_index    the least non-negative integer \f$ j \f$
                              such that the infinite order recurrence
			      is well-defined for \f$ n \geq j \f$.

  \return                     <CODE>true</CODE> if the infinite order
                              recurrence is in <EM>normal form</EM> or
			      is possible to rewrite it in normal form
			      <EM>normal form</EM>; returns
			      <CODE>false</CODE> otherwise.

  The system is able to compute linear infinite order recurrence
  in normal form
  \f[
    x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n).
  \f]
  This is possible using that, for \f$ n > 1 \f$, the sequence
  \f$ x \f$ also satisfies the linear recurrence of first order
  \f[
    x(n) = \frac{f(n)}{f(n-1)} (1+f(n-1)) x(n-1)
      + f(n) \left( \frac{g(n)}{f(n)} - \frac{g(n-1)}{f(n-1)} \right).
  \f]
  For \f$ n = 1 \f$ it must consider \f$ x(1) = f(1) x(0) + g(1) \f$.

  Moreover, this function transform, when possible, infinite order
  recurrence in normal form.
  This transformation is performed in the following way:
  - If the upper limit of the sum is \f$ n \f$ and \f$ f(n) \neq 1 \f$
    then, independently from the lower limit of the sum
    \f$ n_0 \in \Nset \f$, the recurrence
    \f[
      x(n) = f(n) \sum_{k=n_0}^n x(k) + g(n),
    \f]
    where \f$ f(n) \neq 1 \f$, is transformed in
    \f[
      x(n) = \frac{f(n)}{1-f(n)} \sum_{k=n_0}^{n-1} x(k) + \frac{g(n)}{1-f(n)}.
    \f]
  - If the upper limit of the sum is \f$ n \f$ and \f$ f(n) = 1 \f$
    then the previous method in not applicable.
    A method is pull outside from the sum the terms valued in \f$ n \f$
    and in \f$ n -1 \f$:
    \f[
      x(n) = \sum_{k=n_0}^{n-2} x(k) + x(n-1) + x(n) + g(n).
    \f]
    Hence, perform the shift \f$ n = n+1 \f$ so that
    \f[
      x(n) = -\sum_{k=n_0+1}^{n-1} x(k) - g(n+1).
    \f]
  - If the lower bound of the sum \f$ n_0 \f$ is greater than \f$ 0 \f$
    the recurrence
    \f[
      x(n) = f(n) \sum_{k=n_0}^{n-1} x(k) + g(n),
    \f]
    is transformed in
    \f[
      x(n) = f(n+n_0) \sum_{k=0}^{n-1} x(k) + g(n+n_0).
    \f]

  In conclusion the algorithm for classifying and for finding the solution
  of the infinite order recurrence previews the following steps:
  - eventual rewriting in the normal form of the recurrence;
  - computation of the right hand side of the associated first order
    recurrence;
  - shift forward of the first order recurence: \f$ n = n + 1 \f$;
  - computation of the solution of the first order recurrence;
  - shift backward of the solution: \f$ n = n - 1 \f$;
  - substitution of the initail condition \f$ x(1) = f(1) x(0) + g(1) \f$.
*/
bool
rewrite_infinite_order_recurrence(const Expr& rhs, const Expr& term_sum,
				  Expr& weight, Expr& inhomogeneous,
				  Expr& rhs_first_order,
				  index_type& first_valid_index,
				  bool& rewritten) {
  Expr upper = term_sum.arg(2);
  int lower = term_sum.arg(1).ex_to_number().to_int();
  inhomogeneous = rhs - weight * term_sum;

  // The recurrence is too complex for the system in the following cases:
  // - if the weight `f(n)' or the inhomogeneous term `g(n)' contain
  //   other functions `x()' with `n' in the argument;
  // - if the weight `f(n)' or the inhomogeneous term `g(n)' contain
  //   the parameters;
  // - the upper bound of the sum is different from `n' and `n-1'.
  if (weight.has_x_function(Recurrence::n)
      || inhomogeneous.has_x_function(Recurrence::n)
      || has_parameters(weight)
      || (upper != Recurrence::n && upper != Recurrence::n-1))
    return false;

  // Find the weight `f(n)' and the inhomogeneous term of the
  // recurrence transformed so that to have the upper limit of the
  // sum equal to `n'.
  if (upper == Recurrence::n)
    if (weight != 1) {
      rewritten = true;
      const Expr& tmp = 1 - weight;
      weight /= tmp;
      inhomogeneous /= tmp;
    }
    else {
      rewritten = true;
      weight *= -1;
      lower += 1;
      inhomogeneous = -inhomogeneous.substitute(Recurrence::n,
						Recurrence::n + 1) ;
    }
  
  // Find the weight `f(n)' and the inhomogeneous term of the
  // recurrence transformed so that to have the lower limit of the
  // sum equal to `0'.
  if (lower > 0) {
    rewritten = true;
    weight = weight.substitute(Recurrence::n, Recurrence::n + lower);
    inhomogeneous = inhomogeneous.substitute(Recurrence::n,
					     Recurrence::n + lower);
  }
  
  Number z = 0;
  // Find the largest positive or null integer that cancel the numerator of
  // `f(n)' and store it in `z' if it is bigger than the current `z'.
  if (!largest_positive_int_zero(numerator(weight), Recurrence::n, z))
    return false;
  // `z' will contain the largest positive or null integer, if it exists,
  // that cancel the denominator of `f(n)'.
  // If this integer does not exist then `z' is left to 0.
  if (!largest_positive_int_zero(denominator(weight), Recurrence::n, z))
    return false;

  if (!largest_positive_int_zero(denominator(inhomogeneous), Recurrence::n, z))
    return false;
  first_valid_index = z.to_unsigned_int();

  // Find the element will form the first order recurrence.
  const Expr& weight_shifted = weight.substitute(Recurrence::n,
						 Recurrence::n-1);
  Expr coeff_first_order
    = simplify_all(weight / weight_shifted * (1 + weight_shifted));
  // Shift forward: `n -> n + 1'.
  coeff_first_order =
    coeff_first_order.substitute(Recurrence::n, Recurrence::n+1);

  const Expr& inhomogeneous_shifted
    = inhomogeneous.substitute(Recurrence::n, Recurrence::n-1);
  Expr inhomog_first_order
    = simplify_all(weight * (inhomogeneous / weight
			     - inhomogeneous_shifted / weight_shifted));
  // Shift forward: `n -> n + 1'.
  inhomog_first_order
    = inhomog_first_order.substitute(Recurrence::n, Recurrence::n+1);
  
  rhs_first_order
    = coeff_first_order * x(Recurrence::n-1) + inhomog_first_order;

  return true;
}

} // anonymous namespace

PURRS::Recurrence::Classifier_Status
PURRS::Recurrence::compute_order(const Number& decrement, index_type& order,
				 unsigned long& index,
				 unsigned long max_size) {
  if (decrement < 0)
    return HAS_NEGATIVE_DECREMENT;
  // Make sure that (1) we can represent `decrement' as a long, and
  // (2) we will be able to store the coefficient into the
  // appropriate position of the `coefficients' vector.
  if (decrement >= LONG_MAX || decrement >= max_size)
    return CL_HAS_HUGE_DECREMENT;
  
  // The `order' is defined as the maximum value of `index'.
  index = decrement.to_long();
  if (order == 0 || index > unsigned(order))
    order = index;
  return CL_SUCCESS;
}

/*!
  Analyzes the \f$ i \f$-th addend of the right hand side \p rhs of the
  recurrence \p *this collecting and updating all necessary informations
  of which the system it will have need during the computations of the
  solution or of the bounds, or during the verification of the obtained
  results.
*/
PURRS::Recurrence::Classifier_Status
PURRS::Recurrence::classification_summand(const Expr& addend, Expr& rhs,
					  Expr& inhomogeneous,
					  index_type& order,
					  std::vector<Expr>& coefficients,
					  int& gcd_among_decrements,
					  int num_term,
					  std::map<Number, Expr>&
					  homogeneous_terms) const {
  // `non_linear_term == 0' or `non_linear_term == 1' indicate
  // two different cases of non-linearity.
  unsigned int non_linear_term = find_non_linear_term(addend);
  if (non_linear_term == 0) {
    // We will store here the right hand side of the linear recurrence
    // obtained transforming that one non-linear.
    Expr new_rhs;
    // We will store here the symbols associated to the eventual negative
    // numbers that will be the arguments of the logarithms.
    std::vector<Symbol> auxiliary_symbols;
    // In the case of simple non-linear recurrence of the form
    // `x(n) = c x(n-1)^a' `coeff_and_base' will contain the pair `c' and a';
    // in all the other cases the numeric value of the pair will hold `0'
    // while the second element will contain the value that will be the
    // logarithm's base or the exponential's base used in the rewriting
    // of the non-linear recurrence in the correspondent linear recurrence.
    std::pair<Number, Expr> coeff_and_base;
    if (rewrite_non_linear_recurrence(*this, rhs, new_rhs, coeff_and_base,
				      auxiliary_symbols)) {
      set_non_linear_finite_order();
      non_linear_p = new Non_Linear_Info(Recurrence(new_rhs), coeff_and_base,
					 auxiliary_symbols);
      return CL_SUCCESS;
    }
    else
      return CL_TOO_COMPLEX;
  }
  // This is the case of nested `x' function with argument containing `n'.
  else if (non_linear_term == 1)
    return CL_MALFORMED_RECURRENCE;

  unsigned int num_factors = addend.is_a_mul() ? addend.nops() : 1;
  Number num;
  if (num_factors == 1)
    if (addend.has_non_rational_numbers())
      return CL_MALFORMED_RECURRENCE;
    else if (addend.is_the_x_function()) {
      const Expr& argument = addend.arg(0);
      if (argument == n)
	return HAS_NULL_DECREMENT;
      else if (has_parameters(argument))
	return CL_TOO_COMPLEX;
      // Check if this term has the form `x(n + d)'.
      else if (argument.is_a_add() && argument.nops() == 2) {
	Number decrement;
	if (get_constant_decrement(argument, decrement)) {
	  unsigned long index;
	  Classifier_Status status
	    = compute_order(decrement, order, index,
			    coefficients.max_size());
	  if (status != CL_SUCCESS)
	    return status;
	  if (classifier_status_ == NOT_CLASSIFIED_YET || type_ == ORDER_ZERO)
	    set_linear_finite_order_const_coeff();
	  else if (is_functional_equation())
	    return CL_TOO_COMPLEX;
	  // `num_term == 0' if `r' is the unique term of `rhs'
	  // or if it is the first term of `rhs' (i.e. is the
	  // first time that the system entry in this function).
	  if (num_term == 0)
	    gcd_among_decrements = index;
	  else
	    gcd_among_decrements = gcd(gcd_among_decrements, index);
	  insert_coefficients(1, index, coefficients);
	}
	else
	  return CL_HAS_NON_INTEGER_DECREMENT;
      }
      // Check if this term has the form `x(n / d)'.
      else if (argument.is_a_mul() && argument.nops() == 2) {
	Number divisor;
	if (get_constant_divisor(argument, divisor)) {
	  if (classifier_status_ == NOT_CLASSIFIED_YET || type_ == ORDER_ZERO)
	    set_functional_equation();
	  else if (is_linear_finite_order())
	    return CL_TOO_COMPLEX;
	  homogeneous_terms
	    .insert(std::map<Number, Expr>::value_type(divisor, 1));
	}
	else
	  return CL_TOO_COMPLEX;
      }
      else if (argument.has(n))
	return CL_TOO_COMPLEX;
      else {
	inhomogeneous += addend;
	if (classifier_status_ == NOT_CLASSIFIED_YET)
	  set_order_zero();
      }
    } // ended case of `addend' `x' function.
  // Check if the summand has the `x' function with the argument
  // dependently from the index of the sum.
    else if (addend.is_the_sum_function() && addend.arg(2).has(n)
	     && addend.arg(3).has_x_function(addend.arg(0))) {
      if (classifier_status_ == NOT_CLASSIFIED_YET || type_ == ORDER_ZERO) {
	// If there are many terms equal to the sum stored in `addend'
	// we must collect them in order to find the weight `f(n)' of
	// the recurrence of infinite order
	// `x(n) = f(n) sum(k, n_0, u(n), x(k)) + g(n)'.
	Expr weight;
	Expr rhs_rewritten = rhs.collect_term(addend, weight);
	Expr rhs_first_order;
	bool rewritten;
	index_type first_valid_index;
	if (rewrite_infinite_order_recurrence(rhs_rewritten, addend, weight,
					      inhomogeneous, rhs_first_order,
					      first_valid_index, rewritten)) { 
	  if (first_valid_index > 0)
	    return CL_DOMAIN_ERROR;
	  infinite_order_p
	    = new Infinite_Order_Info(Recurrence(rhs_first_order), weight);
	  set_linear_infinite_order();
	  if (rewritten) {
	    bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	    rec_rewritten = true;
	    set_original_rhs(rhs);
	    Symbol h;
	    rhs = weight * PURRS::sum(h, 0, n-1, x(h)) + inhomogeneous;
	  }
	  return CL_SUCCESS;
	}
	else
	  return CL_TOO_COMPLEX;
      }
      else
	return CL_TOO_COMPLEX;
    }
    else {
      inhomogeneous += addend;
      if (classifier_status_ == NOT_CLASSIFIED_YET)
	set_order_zero();
    }
  else {
    Expr no_x_factor = 1;
    bool has_x = false;
    bool has_n = false;
    unsigned long index;
    Number divisor;
    for (unsigned int i = num_factors; i-- > 0; ) {
      const Expr& factor = addend.op(i);
      if (factor.has_non_rational_numbers())
	return CL_MALFORMED_RECURRENCE;
      else if (factor.is_the_x_function()) {
	const Expr& argument = factor.arg(0);
	if (argument == n)
	  return HAS_NULL_DECREMENT;
	else if (has_parameters(argument))
	  return CL_TOO_COMPLEX;
	else if (argument.is_a_add() && argument.nops() == 2) {
	  Number decrement;
	  if (get_constant_decrement(argument, decrement)) {
	    // The non linear terms have already been considered before.
	    assert(!has_x);
	    Classifier_Status status
	      = compute_order(decrement, order, index,
			      coefficients.max_size());
	    if (status != CL_SUCCESS)
	      return status;
	    if (classifier_status_ != NOT_CLASSIFIED_YET
		&& is_functional_equation())
	      return CL_TOO_COMPLEX;
	    if (num_term == 0)
	      gcd_among_decrements = index;
	    else
	      gcd_among_decrements = gcd(gcd_among_decrements, index);
	    has_x = true;
	  }
	  else
	    return CL_HAS_NON_INTEGER_DECREMENT;
	}
	else if (argument.is_a_mul() && argument.nops() == 2) {
	  if (get_constant_divisor(argument, divisor)) {
	    // The non linear terms have already been considered before.
	    assert(!has_x);
	    if (classifier_status_ == NOT_CLASSIFIED_YET
		|| type_ == ORDER_ZERO)
	      set_functional_equation();
	    else if (is_linear_finite_order())
	      return CL_TOO_COMPLEX;
	    has_x = true;
	  }
	  else
	    return CL_TOO_COMPLEX;
	}
	else if (argument.has(n))
	  return CL_TOO_COMPLEX;
	else
	  no_x_factor *= factor;
      } // ended case of `factor' `x' function.
      // Check if the summand has the `x' function with the argument
      // dependently from the index of the sum.
      else if (factor.is_the_sum_function() && factor.arg(2).has(n)
	       && factor.arg(3).has_x_function(factor.arg(0))) {
	if (classifier_status_ == NOT_CLASSIFIED_YET || type_ == ORDER_ZERO) {
	  // If there are many terms equal to the sum in `factor' stored
	  // in `rhs', we must collect them in order to find the weight
	  // `f(n)' of the recurrence of infinite order
	  // `x(n) = f(n) sum(k, n_0, u(n), x(k)) + g(n)'.
	  Expr weight;
	  Expr rhs_rewritten = rhs.collect_term(factor, weight);
	  // There are not other sums equal to `factor'.
	  if (weight == 1)
	    for (unsigned int j = num_factors; j-- > 0; )
	      if (addend.op(j) != factor)
		weight *= addend.op(j);
	  Expr rhs_first_order;
	  bool rewritten;
	  index_type first_valid_index;
	  if (rewrite_infinite_order_recurrence(rhs_rewritten, factor, weight,
						inhomogeneous, rhs_first_order,
						first_valid_index,
						rewritten)) { 
	    if (first_valid_index > 0)
	      return CL_DOMAIN_ERROR;
	    infinite_order_p
	      = new Infinite_Order_Info(Recurrence(rhs_first_order), weight);
	    set_linear_infinite_order();
	    if (rewritten) {
	      bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	      rec_rewritten = true;
	      set_original_rhs(rhs);
	      Symbol h;
	      rhs = weight * PURRS::sum(h, 0, n-1, x(h)) + inhomogeneous;
	    }
	    return CL_SUCCESS;
	  }
	  else
	    return CL_TOO_COMPLEX;
	}
	else
	  return CL_TOO_COMPLEX;
      }
      else {
	if (factor.has(n))
	  has_n = true;
	no_x_factor *= factor;
      }
    }
    if (has_x) {
      if (classifier_status_ != NOT_CLASSIFIED_YET
	  && is_functional_equation())
	homogeneous_terms
	  .insert(std::map<Number, Expr>::value_type(divisor, no_x_factor));
      else {
	insert_coefficients(no_x_factor, index, coefficients);
	if (classifier_status_ == NOT_CLASSIFIED_YET
	    || (classifier_status_ != NOT_CLASSIFIED_YET
		&& !is_linear_finite_order_var_coeff()))
	  if (has_n)
	    set_linear_finite_order_var_coeff();
	  else
	    set_linear_finite_order_const_coeff();
      }
    }
    else {
      inhomogeneous += no_x_factor;
      if (classifier_status_ == NOT_CLASSIFIED_YET)
	set_order_zero();
    }
  }
  return CL_SUCCESS;
}

/*!
  Classifies the recurrence \p *this.
  Returns:
  - <CODE>CL_SUCCESS</CODE>         if the recurrence is linear of
                                           finite or infinite order;
					   non-linear of finite order;
					   in the case of functional equation;
  - <CODE>CL_HAS_NON_INTEGER_DECREMENT</CODE> if the right-hand side of the
                                           recurrence contains at least an
					   occurrence of <CODE>x(n-k)</CODE>,
					   where <CODE>k</CODE> is not an
					   integer;
  - <CODE>HAS_NEGATIVE_DECREMENT</CODE>    if the right-hand side of the
                                           recurrence contains at least an
					   occurrence of <CODE>x(n-k)</CODE>,
					   where <CODE>k</CODE> is a negative
					   integer;
  - <CODE>HAS_NULL_DECREMENT</CODE>        if the right-hand side of the
                                           recurrence contains at least an
					   occurrence of <CODE>x(n)</CODE>;
  - <CODE>CL_HAS_HUGE_DECREMENT</CODE>        if the right-hand side of the
                                           recurrence contains at least an
					   occurrence of <CODE>x(n-k)</CODE>,
					   where <CODE>k</CODE> is too big to
					   be handled by the standard
					   solution techniques;
  - <CODE>CL_MALFORMED_RECURRENCE</CODE>      if the recurrence does not have
                                           any sense;
  - <CODE>CL_DOMAIN_ERROR</CODE>              if the recurrence is not
                                           well-defined;
  - <CODE>CL_UNSOLVABLE_RECURRENCE</CODE>     if the recurrence is not solvable;
  - <CODE>CL_INDETERMINATE_RECURRENCE</CODE>  if the recurrence is indeterminate,
                                           hence it has infinite solutions.
  - <CODE>CL_TOO_COMPLEX</CODE>    in all the other cases.

  For each class of recurrences for which the system returns
  <CODE>CL_SUCCESS</CODE>, it is initialized a pointer to an opportune
  class containing all necessary informations about the recurrence
  of which the system it will have need during the computations of the
  solution or of the bounds, or during the verification of the obtained
  results.
  
  If the recurrence is <EM>linear of finite order</EM> in the structure
  <CODE>Finite_Order_Info</CODE> are stored the order, the coefficients,
  the greates common divisor among the positive integer \f$ k \f$
  of the terms of the form \f$ x(n-k) \f$ contained in the right hand side
  and the first valid index (i. e., the least non-negative integer \f$ j \f$
  such that the recurrence is well-defined for \f$ n \geq j \f$.)
  setted to \f$ 0 \f$ by default.

  If the recurrence is <EM>linear of infinite order</EM> there is not
  a general method to solve it. At the moment we solve only recurrences
  of infinite order of the shape
  \f[
    T(n) = f(n) \sum_{k=0}^{n-1} T(k) + g(n)
  \f]
  transforming it in the linear recurrence of first order
  \f[
    T(n) = \frac{f(n)}{f(n-1)} (1+f(n-1)) T(n-1)
      + f(n) \left( \frac{g(n)}{f(n)} - \frac{g(n-1)}{f(n-1)} \right).
  \f].
  In the structure <CODE>Infinite_Order_Info</CODE> are stored the right
  hand side, the coefficient and the inhomogeneous term of the finite order
  recurrence in which \p *this is transformed, the weight \f$ f(n) \f$.

  If the recurrence is <EM>non-linear of finite order</EM> there is not
  a general method to solve it. At the moment we solve only a small class
  of them transforming them in linear recurrences.
  In the structure <CODE>Non_Linear_Info</CODE> are stored the right hand
  side of the linear recurrence, the logarithm's base or the exponential's
  base used in the transformation, a vector of symbols associated to the
  eventual negative numbers that will be the arguments of the logarithms.

  In the case of <EM>functional equation</EM> in the structure
  <CODE>Functional_Equation_Info</CODE> are stored the values
  \f$ a \f$ and \f$ b \f$ of the terms \f$ a x_{n/b} \f$ contained
  in the right hand side, the positive integer starting from which the
  inhomogeneous term is a non negative, non decreasing function (setted to
  \f$ 1 \f$ by default).

  In all the previous cases is besides computed the non-homogeneous part
  of the recurrence.
*/
PURRS::Recurrence::Classifier_Status
PURRS::Recurrence::classify() const {
  Expr& rhs = const_cast<Expr&>(recurrence_rhs);
  // Simplifies expanded expressions, in particular rewrites nested powers.
  rhs = simplify_ex_for_input(recurrence_rhs, true);
  // Splits the sum in many sums how many are the addends of the summand
  // and computes, when possible, symbolic sums.
  rhs = simplify_sum(rhs, COMPUTE_SUM);

  // Date for linear finite order recurrences.

  // Initialize the computation of the order of the linear part of the
  // recurrence.  This works like the computation of a maximum: it is
  // the maximum `k', if it exists, such that `rhs = a*x(n-k) + b' where `a'
  // is not syntactically 0; if not exists `k' such that `rhs = a*x(n-k) + b',
  // then `order' is left to `0'.
  index_type order = 0;
  // We will store here the coefficients of linear part of the recurrence.
  std::vector<Expr> coefficients;

  // We will store here the greatest common denominator among the decrements
  // `d' of the terms `x(n-d)' contained in the linear part of the
  // recurrence.
  int gcd_among_decrements = 0;

  // Date for functional equations.

  std::map<Number, Expr> homogeneous_terms;

  // Date for all types of recurrences
  Expr inhomogeneous = 0;

  Classifier_Status status;

  unsigned int num_summands = rhs.is_a_add() ? rhs.nops() : 1;
  if (num_summands > 1)
    // It is necessary that the following loop starts from `0'.
    for (unsigned int i = 0; i < num_summands; ++i) {
      if ((status = classification_summand(rhs.op(i), rhs, inhomogeneous,
					   order, coefficients,
					   gcd_among_decrements, i,
					   homogeneous_terms))
	  != CL_SUCCESS)
	return status;
      // As soon as the system notices the this is recurrences
      // of infinite order, stops the classification because already
      // have all necessary information in order to continue.
      if (is_linear_infinite_order())
	break;
    }
  else
    if ((status = classification_summand(rhs, rhs, inhomogeneous,
					 order, coefficients,
					 gcd_among_decrements, 0,
					 homogeneous_terms))
	!= CL_SUCCESS)
      return status;

  set_inhomogeneous_term(inhomogeneous);
  D_MSGVAR("Inhomogeneous term: ", inhomogeneous_term);

  // In the case of non linear recurrence or infinite order recurrence
  // we have already done the operation `new ...'.
  if (is_functional_equation())
    functional_eq_p = new Functional_Equation_Info(homogeneous_terms);
  else if (is_linear_finite_order())
    finite_order_p = new Finite_Order_Info(order, coefficients,
					   gcd_among_decrements);

  assert(is_linear_finite_order() || is_functional_equation()
	 || is_non_linear_finite_order() || is_linear_infinite_order());

  return CL_SUCCESS;
}

/*!
  Classifies the recurrence \p *this calling the method
  <CODE>classify()</CODE>.
  If the function <CODE>classify()</CODE> returns the value
  <CODE>HAS_NEGATIVE_DECREMENT</CODE> or the value
  <CODE>HAS_NULL_DECREMENT</CODE>, this function tries to rewrite the
  recurrence in the normal form \f$ x(n) = r \f$, where \f$ r \f$
  does not contain terms <CODE>x(n-k)</CODE>, where <CODE>k</CODE>
  is not a positive integer.

  Returns:
  - <CODE>CL_SUCCESS</CODE>         if the recurrence is linear of
                                           finite or infinite order;
					   non-linear of finite order;
					   in the case of functional equation;
  - <CODE>CL_HAS_NON_INTEGER_DECREMENT</CODE> if the right-hand side of the
                                           recurrence contains at least an
					   occurrence of <CODE>x(n-k)</CODE>,
					   where <CODE>k</CODE> is not an
					   integer;
  - <CODE>CL_HAS_HUGE_DECREMENT</CODE>        if the right-hand side of the
                                           recurrence contains at least an
					   occurrence of <CODE>x(n-k)</CODE>,
					   where <CODE>k</CODE> is too big to
					   be handled by the standard
					   solution techniques;
  - <CODE>CL_MALFORMED_RECURRENCE</CODE>      if the recurrence does not have
                                           any sense;
  - <CODE>CL_DOMAIN_ERROR</CODE>              if the recurrence is not
                                           well_defined;
  - <CODE>CL_UNSOLVABLE_RECURRENCE</CODE>     if the recurrence is not solvable;
  - <CODE>CL_INDETERMINATE_RECURRENCE</CODE>  if the recurrence is indeterminate,
                                           hence it has infinite solutions.
  - <CODE>CL_TOO_COMPLEX</CODE>    in all the other cases. 
  
  For each class of recurrences for which the system returns
  <CODE>CL_SUCCESS</CODE>, it is initialized a pointer to an opportune
  class containing all necessary informations about the recurrence
  of which the system it will have need during the computations of the
  solution or of the bounds, or during the verification of the obtained
  results.
*/
PURRS::Recurrence::Classifier_Status
PURRS::Recurrence::classify_and_catch_special_cases() const {
  // If `*this' is already classified returns as soon as its status.
  if (classifier_status_ != NOT_CLASSIFIED_YET)
    return classifier_status_;

  bool exit_anyway = false;
  Classifier_Status status;
  do {
    status = classify();
    switch (status) {
    case CL_SUCCESS:
      break;
    case CL_HAS_NON_INTEGER_DECREMENT:
    case CL_HAS_HUGE_DECREMENT:
    case CL_MALFORMED_RECURRENCE:
    case CL_DOMAIN_ERROR:
    case CL_TOO_COMPLEX:
      exit_anyway = true;
      break;
    case HAS_NEGATIVE_DECREMENT:
      {
	Expr new_rhs;
	eliminate_negative_decrements(recurrence_rhs, new_rhs);
	bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	rec_rewritten = true;
	Expr& rhs = const_cast<Expr&>(recurrence_rhs);
	rhs = new_rhs;
	classifier_status_ = NOT_CLASSIFIED_YET;
	type_ = UNKNOWN;
	status = classify_and_catch_special_cases();
      }
      break;
    case HAS_NULL_DECREMENT:
      {
	Expr new_rhs;
	unsigned int result
	  = eliminate_null_decrements(recurrence_rhs, new_rhs);
	if (result == 0) {
	  bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	  rec_rewritten = true;
	  Expr& rhs = const_cast<Expr&>(recurrence_rhs);
	  rhs = new_rhs;
	  classifier_status_ = NOT_CLASSIFIED_YET;
	  type_ = UNKNOWN;
	  status = classify_and_catch_special_cases();
	}
	else if (result == 1)
	  status = CL_UNSOLVABLE_RECURRENCE;
	else if (result == 2)
	  status = CL_INDETERMINATE_RECURRENCE;
	else
	  status = CL_TOO_COMPLEX;
	exit_anyway = true;
      }
      break;

    default:
      throw std::runtime_error("PURRS internal error: "
			       "catch_special_cases().");
      break;
    }
  } while (!exit_anyway && status != CL_SUCCESS);

  classifier_status_ = status;
  return status;
}
