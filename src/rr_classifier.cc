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
#include "Cached_Expr.defs.hh"
#include "Order_Reduction_Info.defs.hh"
#include "Recurrence.defs.hh"

#include <climits>
#include <string>
#include <vector>

// TEMPORARY
#include <iostream>
#include <fstream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

/*!
  Returns <CODE>true</CODE> if \p e is of the form \f$ n - d \f$ with
  \f$ d \f$ an integer: in this case assign the opposite of \f$ d \f$ to
  \p decrement.
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

/*!
  Returns <CODE>true</CODE> if \p e is of the form \f$ n / d \f$ with
  \f$ d \f$ a positive integer: in this case assign \f$ d \f$ to \p divisor.
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
      && numerator(d) == 1) {
    divisor = d.denominator();
    assert(divisor.is_positive_integer());
    return true;
  }
  else if (b == Recurrence::n && a.is_a_number(d) && d.is_rational()
	   && numerator(d) == 1) {
    divisor = d.denominator();
    assert(divisor.is_positive_integer());
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
    for (unsigned i = e.nops(); i-- > 0; )
      e_substituted += change_variable_function_x(e.op(i), s, r);
  }
  else if (e.is_a_mul()) {
    e_substituted = 1;
    for (unsigned i = e.nops(); i-- > 0; )
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
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned j = 0; j < num_argument; ++j)
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
    for (unsigned i = e.nops(); i-- > 0; ) {
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
    for (unsigned i = e.nops(); i-- > 0; )
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
  for (unsigned i = e.nops(); i-- > 0; ) {
    const Expr& term = e.op(i);
    if (term.is_a_mul()) {
      Expr tmp = 1;
      bool found_x_n = false;
      for (unsigned j = term.nops(); j-- > 0; ) {
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
*/
unsigned
eliminate_null_decrements(const Expr& rhs, Expr& new_rhs) {
  // Let `rhs = a*x(n) + b' and that `b' does different to zero
  // and does not contain `x(n)'.  The following cases are possible:
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
  if (rhs.is_a_add()) {
    // Finds `a' and `b'.
    Expr a = 0;
    Expr b = rhs;
    find_coeff_x_n_and_remainder(rhs, a, b);
    if (a == 1) {
      // Case 2. and Case 3.
      bool found_x = false;
      if (b.is_a_add())
	for (unsigned i = b.nops(); i-- > 0; ) {
	  const Expr& term = b.op(i);
	  if (term.is_a_mul()) {
	    for (unsigned j = term.nops(); j-- > 0; )
	      if (term.op(j).is_the_x_function()) {
		found_x = true;
		break;
	      }
	  }
	  else
	    if (term.is_the_x_function()) {
	      found_x = true;
	      break;
	    }
	}
      else if (b.is_a_mul())
	for (unsigned i = b.nops(); i-- > 0; ) {
	  if (b.op(i).is_the_x_function()) {
	    found_x = true;
	    break;
	  }
	}
      else if (b.is_the_x_function())
	found_x = true;
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
  else if (rhs.is_a_mul())
    for (unsigned i = rhs.nops(); i-- > 0; ) {
      const Expr& factor = rhs.op(i);
      if (factor == x(Recurrence::n))
	new_rhs = 0;
    }
  // Let `rhs = x(n)'.
  else if (rhs == x(Recurrence::n))
    return 2;
  return 0;
}

/*!
  Returns <CODE>true</CODE> if \p e contains the \f$ x \f$ function
  with its argument containing the symbol \p s.
  Returns <CODE>false</CODE> otherwise. 
*/
bool
has_x_function(const Expr& e, const Symbol& s) {
  if (e.is_a_add() || e.is_a_mul()) {
    for (unsigned i = e.nops(); i-- > 0; )
      if (has_x_function(e.op(i), s))
	return true;
    return false;
  }
  else if (e.is_a_power()) {
    if (has_x_function(e.arg(0), s) || has_x_function(e.arg(1), s))
      return true;
    return false;
  }
  else if (e.is_a_function()) {
    if (e.is_the_x_function() && e.arg(0).has(s))
      return true;
    return false;
  }
  else
    return false;
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
  Returns <CODE>true</CODE> in two cases:
  - there is in \p e an \f$ x \f$ function (with the argument containing
    \f$ n \f$) inside an other function;
  - there is in \p e a power with an \f$ x \f$ function in the base.
  Returns <CODE>false</CODE> in all the other cases. 
*/
bool
x_function_in_powers_or_functions(const Expr& e) {
  // First case.
  if (e.is_a_function()) {
    for (unsigned i = e.nops(); i-- > 0; ) {
      const Expr& operand = e.arg(i);
      if (operand.is_the_x_function())
	if (operand.arg(0).has(Recurrence::n))
	  return true;
    }
  }
  // Second case.
  else if (e.is_a_power()) {
    const Expr& base = e.arg(0);
    if (base.is_the_x_function())
      if (base.arg(0).has(Recurrence::n))
	return true;
  }
  return false;
}

/*!
  Returns <CODE>true</CODE> if finds non linear term;
  returns <CODE>false</CODE> otherwise.
*/
bool
find_non_linear_recurrence(const Expr& e) {
  unsigned num_summands = e.is_a_add() ? e.nops() : 1;
  if (num_summands > 1)
    for (unsigned i = num_summands; i-- > 0; ) {
      const Expr& term = e.op(i);
      unsigned num_factors = term.is_a_mul() ? term.nops() : 1;
      if (num_factors == 1) {
	if (x_function_in_powers_or_functions(term))
	  return true;
      }
      else {
	bool found_function_x = false;
	for (unsigned j = num_factors; j-- > 0; ) {
	  const Expr& factor = term.op(j);
	  if (x_function_in_powers_or_functions(factor))
	    return true;
	  if (factor.is_the_x_function())
	    if (found_function_x)
	      return true;
	    else
	      if (factor.arg(0).has(Recurrence::n))
		found_function_x = true;
	}
      }
    }
  else {
    unsigned num_factors = e.is_a_mul() ? e.nops() : 1;
    if (num_factors == 1) {
      if (x_function_in_powers_or_functions(e))
	return true;
    }
    else {
      bool found_function_x = false;
      for (unsigned j = num_factors; j-- > 0; ) {
	const Expr& factor = e.op(j);
	if (x_function_in_powers_or_functions(factor))
	  return true;
	if (factor.is_the_x_function())
	  if (found_function_x)
	    return true;
	  else
	    if (factor.arg(0).has(Recurrence::n))
	      found_function_x = true;
      }
    }
  }
  return false;
}

//! \brief
//! Returns <CODE>true</CODE> if the non linear recurrence \f$ x(n) = rhs \f$
//! is rewritable as linear recurrence \f$ x(n) = new_rhs \f$.
//! Returns <CODE>false</CODE> otherwise.
/*!
  Cases of rewritable non linear recurrences:
  -  \f$ x(n) = a * x(n-k_1)^b_1 * ... * x(n-k_h)^b_h \f$,
     where \f$ k_1, \dots, k_h, b_1, \dots, b_h \f$ are positive integers
     and \f$ h > 1 \f$ or \f$ b_1 > 1 \f$;
  -  x(n) = x(n-k)^b, where \f$ b > 1 \f$ and \f$ k \f$ are positive integers.
*/
bool
rewrite_non_linear_recurrence(const Expr& rhs, Expr& new_rhs, Expr& base) {
  D_MSGVAR("*** ", rhs);
  // First case.
  if (rhs.is_a_mul()) {
    bool simple_cases = false;
    Number common_exponent = 1;
    for (unsigned i = rhs.nops(); i-- > 0; ) {
      const Expr& factor = rhs.op(i);
      Number num_exp;
      if (factor.is_a_power() && factor.arg(0).is_the_x_function()
	  && factor.arg(1).is_a_number(num_exp)
	  && num_exp.is_positive_integer()) {
	simple_cases = true;
	common_exponent = lcm(num_exp, common_exponent);
      }
      if (factor.is_the_x_function()) {
	simple_cases = true;
	common_exponent = lcm(1, common_exponent);
      }
    }
    D_VAR(common_exponent);
    new_rhs = 0;
    if (simple_cases)
      if (common_exponent == 1) {
	// Substitute any `x' function `f' with `exp(1)^f'.
	base = exp(Expr(1));
	Expr tmp = substitute_x_function(rhs, exp(Expr(1)), true);
	tmp = simplify_ex_for_input(tmp, true);
	for (unsigned i = tmp.nops(); i-- > 0; )
	  new_rhs += log(tmp.op(i));
	new_rhs = simplify_logarithm(new_rhs);
	D_VAR(new_rhs);
	return true;
      }
      else {
	D_VAR(common_exponent);
	// Substitute any `x' function `f' with `common_exponent^f'.
	base = common_exponent;
	Expr tmp = substitute_x_function(rhs, common_exponent, true);
	tmp = simplify_ex_for_input(tmp, true);
	for (unsigned i = tmp.nops(); i-- > 0; )
	  new_rhs += log(tmp.op(i)) / log(Expr(common_exponent));
	new_rhs = simplify_logarithm(new_rhs);
	D_VAR(new_rhs);
	return true;
      }
  }
  // Second case.
  else if (rhs.is_a_power()) {
    const Expr& base_rhs = rhs.arg(0);
    const Expr& exponent_rhs = rhs.arg(1);
    Number num_exp;
    if (base_rhs.is_the_x_function() && exponent_rhs.is_a_number(num_exp)
	&& num_exp.is_positive_integer()) {
      unsigned exponent_x_function = num_exp.to_unsigned();
      // Substitute any `x' function `f' with `exponent_x_function^f'.
      base = exponent_x_function;
      Expr tmp = substitute_x_function(rhs, exponent_x_function, true);
      tmp = simplify_ex_for_input(tmp, true);
      new_rhs += log(tmp) / log(Expr(exponent_x_function));
      new_rhs = simplify_logarithm(new_rhs);
      D_VAR(new_rhs);
      return true;
    }
    else
      return false;
  }
  return false;
}

} // anonymous namespace

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::compute_order(const Number& decrement, unsigned int& order,
				 unsigned long& index,
				 unsigned long max_size) {
  if (decrement < 0)
    return HAS_NEGATIVE_DECREMENT;
  // Make sure that (1) we can represent `decrement' as a long, and
  // (2) we will be able to store the coefficient into the
  // appropriate position of the `coefficients' vector.
  if (decrement >= LONG_MAX || decrement >= max_size)
    return HAS_HUGE_DECREMENT;
  
  // The `order' is defined as the maximum value of `index'.
  index = decrement.to_long();
  if (order == 0 || index > unsigned(order))
    order = index;
  return SUCCESS;
}

PURRS::Recurrence::Solver_Status
PURRS::Recurrence::classification_summand(const Expr& r, Expr& e,
					  std::vector<Expr>& coefficients,
					  unsigned int& order,
					  int& gcd_among_decrements,
					  int num_term,
					  Expr& coefficient,
					  unsigned& divisor_arg) const {
  unsigned num_factors = r.is_a_mul() ? r.nops() : 1;
  if (num_factors == 1)
    if (r.is_the_x_function()) {
      const Expr& argument = r.arg(0);
      if (argument == n)
	return HAS_NULL_DECREMENT;
      else if (has_parameters(argument))
	return TOO_COMPLEX;
      // Check if this term has the form `x(n + d)'.
      else if (argument.is_a_add() && argument.nops() == 2) {
	Number decrement;
	if (get_constant_decrement(argument, decrement)) {
	  unsigned long index;
	  Solver_Status status
	    = compute_order(decrement, order, index, coefficients.max_size());
	  if (status != SUCCESS)
	    return status;
	  if (num_term == 0)
	    gcd_among_decrements = index;
	  else
	    gcd_among_decrements = gcd(gcd_among_decrements, index);
	  insert_coefficients(1, index, coefficients);
	  if (is_order_zero() || is_unknown())
	    set_linear_finite_order_const_coeff();
	  else if (is_functional_equation())
	    return TOO_COMPLEX;
	}
	else
	  return HAS_NON_INTEGER_DECREMENT;
      }
      // Check if this term has the form `x(n / d)'.
      else if (argument.is_a_mul() && argument.nops() == 2) {
	Number divisor;
	if (get_constant_divisor(argument, divisor)) {
	  coefficient = 1;
	  divisor_arg = divisor.to_unsigned(); 
	}
	else
	  return TOO_COMPLEX;
	if (is_order_zero() || is_unknown())
	  set_functional_equation();
	// `else if (!is_functional_equation())' allows also more than one
	// term of the form `x(n/b)'.
	else
	  return TOO_COMPLEX;
      }
      else if (argument.has(n))
	return TOO_COMPLEX;
      else
	e += r;
    } // end case of r `x' function.
    else if (r.is_the_sum_function() && r.arg(2) == n) {
      // Check if the summand has the `x' function with the argument
      // dependently from the index of the sum.
      if (has_x_function(r.arg(3), r.arg(0).ex_to_symbol())) {
	D_MSG("infinite order");
	return TOO_COMPLEX;
      }
    }
    else
      e += r;
  else {
    Expr possibly_coeff = 1;
    bool found_function_x = false;
    bool found_n = false;
    unsigned long index;
    for (unsigned i = num_factors; i-- > 0; ) {
      const Expr& factor = r.op(i);
      if (factor.is_the_x_function()) {
	const Expr& argument = factor.arg(0);
	if (argument == n)
	  return HAS_NULL_DECREMENT;
	else if (has_parameters(argument))
	  return TOO_COMPLEX;
	else if (argument.is_a_add() && argument.nops() == 2) {
	  Number decrement;
	  if (get_constant_decrement(argument, decrement)) {
#if 0
	    if (found_function_x) {
	      set_non_linear_finite_order();
	      return SUCCESS;
	    }
#else
	    assert(!found_function_x);
#endif
	    Solver_Status status
	      = compute_order(decrement, order, index,
			      coefficients.max_size());
	    if (status != SUCCESS)
	      return status;
	    if (num_term == 0)
	      gcd_among_decrements = index;
	    else
	      gcd_among_decrements = gcd(gcd_among_decrements, index);
	    found_function_x = true;
	  }
	  else
	    return HAS_NON_INTEGER_DECREMENT;
	}
	else if (argument.is_a_mul() && argument.nops() == 2) {
	  Number divisor;
	  if (get_constant_divisor(argument, divisor)) {
	    divisor_arg = divisor.to_unsigned();
	    found_function_x = true;
	    if (is_order_zero() || is_unknown())
	      set_functional_equation();
	    // `else if (!is_functional_equation())' allows also more than one
	    // term of the form `x(n/b)'.
	    else
	      return TOO_COMPLEX;
	  }
	  else
	    return TOO_COMPLEX;
	}
	else if (argument.has(n))
	  return TOO_COMPLEX;
	else
	  possibly_coeff *= factor;
      }
      else if (factor.is_the_sum_function() && factor.arg(2) == n) {
	// Check if the summand has the `x' function with the argument
	// dependently from the index of the sum.
	if (has_x_function(factor.arg(3), factor.arg(0).ex_to_symbol())) {
	  D_MSG("infinite order");
	  return TOO_COMPLEX;
	}
      }
      else {
	if (factor.has(n))
	  found_n = true;
	possibly_coeff *= factor;
      }
    }
    if (found_function_x)
      if (is_functional_equation())
	coefficient = possibly_coeff;
      else {
	insert_coefficients(possibly_coeff, index, coefficients);
	if (!is_linear_finite_order_var_coeff())
	  if (found_n)
	    set_linear_finite_order_var_coeff();
	  else
	    set_linear_finite_order_const_coeff();
      }
    else
      e += possibly_coeff;
  }
  return SUCCESS;
}

/*!
  Returns <CODE>SUCCESS</CODE> if the recurrence is linear and of finite
  order or if is the case of functional equation.
  In the first case are stored in the structure <CODE>Finite_Order_Info</CODE>
  the order of the recurrence, the first initial condition
  (i. e., the smallest positive integer for which the recurrence is
  well-defined) and the coefficients.
  In the second case are stored
  in the structure <CODE>Functional_Equation_Info</CODE> the values
  \f$ a \f$ and \f$ b \f$ of the functional equation
  \f$ x_n = a x_{n/b} + d n^e \f$.
  In both the case is besides computed the non-homogeneous part
  of the recurrence.
  If the function not returns <CODE>SUCCESS</CODE> then is one of the
  following cases: <CODE>HAS_NEGATIVE_DECREMENT</CODE>,
  <CODE>HAS_HUGE_DECREMENT</CODE>, <CODE>HAS_NULL_DECREMENT</CODE>,
  <CODE>HAS_NON_INTEGER_DECREMENT</CODE>, <CODE>TOO_COMPLEX</CODE>.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::classification_recurrence(const Expr& rhs,
					     int& gcd_among_decrements) const {
  // Check if the inhomogeneous term contains floating point numbers.
  if (rhs.has_floating_point_numbers())
    return MALFORMED_RECURRENCE;

  if (find_non_linear_recurrence(rhs)) {
    set_non_linear_finite_order();
    D_MSG("NON LINEAR");
    return SUCCESS;
  }

  // Initialize the computation of the order of the linear part of the
  // recurrence.  This works like the computation of a maximum: it is
  // the maximum `k', if it exists, such that `rhs = a*x(n-k) + b' where `a'
  // is not syntactically 0; if not exists `k' such that `rhs = a*x(n-k) + b',
  // then `order' is left to `0'.
  unsigned int order = 0;
  // We will store here the coefficients of linear part of the recurrence.
  std::vector<Expr> coefficients;

  // We will store here the coefficient of the functional equation
  // `x(n) = a x(n/b) + d n^e'.
  Expr coefficient;
  // We will store here the positive integer divisor of the argument of the
  // function `x' in the functional equation `x(n) = a x(n/b) + d n^e'.
  unsigned divisor_arg;

  Expr inhomogeneous = 0;

  Solver_Status status;

  unsigned num_summands = rhs.is_a_add() ? rhs.nops() : 1;
  if (num_summands > 1)
    // It is necessary that the following loop starts from `0'.
    for (unsigned i = 0; i < num_summands; ++i) {
      if ((status = classification_summand(rhs.op(i), inhomogeneous,
					   coefficients, order,
					   gcd_among_decrements, i,
					   coefficient, divisor_arg))
	  != SUCCESS)
	return status;
    }
  else
    if ((status = classification_summand(rhs, inhomogeneous,
					 coefficients, order,
					 gcd_among_decrements, 0,
					 coefficient, divisor_arg))
	!= SUCCESS)
      return status;

  set_inhomogeneous_term(inhomogeneous);
  D_MSGVAR("Inhomogeneous term: ", inhomogeneous_term);
  
  if (!is_functional_equation())
    // `inhomogeneous_term' is a function of `n', the parameters and of
    // `x(k_1)', ..., `x(k_m)' where `m >= 0' and `k_1', ..., `k_m' are
    //  non-negative integers.
    if (order == 0)
      set_order_zero();

  if (is_linear_finite_order())
    finite_order_p = new Finite_Order_Info(order, 0, coefficients);
  else if (is_functional_equation())
    functional_eq_p = new Functional_Equation_Info(coefficient, divisor_arg);
  
  return SUCCESS;
}

/*!
  This function solves recurrences of SOME TYPE provided they are
  supplied in SOME FORM. (Explain.)
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::classify() const {
  D_VAR(recurrence_rhs);
  // Simplifies expanded expressions, in particular rewrites nested powers.
  Expr rhs = simplify_ex_for_input(recurrence_rhs, true);
  D_VAR(rhs);

  // We will store here the greatest common denominator among the decrements
  // `d' of the terms `x(n-d)' contained in the linear part of the
  // recurrence.
  int gcd_among_decrements = 0;

  Solver_Status status;

  // The function `classification_recurrence()' returns `SUCCESS' if
  // is the order is finite or if is the case of functional equation.
  if ((status = classification_recurrence(rhs, gcd_among_decrements))
      != SUCCESS)
    return status;
  assert(is_linear_finite_order() || is_functional_equation()
	 || is_non_linear_finite_order());

  D_VAR(gcd_among_decrements);
  if (finite_order_p != 0)
    // If the greatest common divisor among the decrements is greater than one,
    // the order reduction is applicable.
    // FIXME: the order reduction is for the moment applied only to
    // recurrences with constant coefficients because the recurrences
    // with variable coefficients are not allowed with parameters.
    if (gcd_among_decrements > 1 && is_linear_finite_order_const_coeff()) {
      recurrence_rhs_rewritten = true;
      Symbol r = insert_auxiliary_definition(mod(n, gcd_among_decrements));
      order_reduction_p
	= new Order_Reduction_Info(rhs, gcd_among_decrements, r);
      // Build the new recurrence substituting `n' not contained in the
      // `x' functions with `gcd_among_decrements * n + r' and `x(n-k)' with
      // `x(n - k / gcd_among_decrements)'.
      recurrence_rhs = rewrite_reduced_order_recurrence(rhs, r,
							gcd_among_decrements);
      status = classify();
    }

  if (is_non_linear_finite_order()) {
    Expr new_rhs;
    Expr base;
    if (rewrite_non_linear_recurrence(rhs, new_rhs, base)) {
      recurrence_rhs_rewritten = true;
      non_linear_p = new Non_Linear_Info(recurrence_rhs, base);
      recurrence_rhs = new_rhs;
      status = classify();
    }
    else
      return TOO_COMPLEX;
  }
  
  return status;
}

/*!
  This function solves recurrences of SOME TYPE provided they
  are supplied in SOME FORM. (Explain.)
  It does that by repeatedly calling solve() and handling
  the errors that may arise.
*/
PURRS::Recurrence::Solver_Status
PURRS::Recurrence::classify_and_catch_special_cases() const {
  bool exit_anyway = false;
  Solver_Status status;
  do {
    status = classify();
    switch (status) {
    case SUCCESS:
      break;
    case HAS_NON_INTEGER_DECREMENT:
    case HAS_HUGE_DECREMENT:
    case TOO_COMPLEX:
      exit_anyway = true;
      break;
    case HAS_NEGATIVE_DECREMENT:
      {
	Expr new_rhs;
	eliminate_negative_decrements(recurrence_rhs, new_rhs);
	recurrence_rhs_rewritten = true;
	recurrence_rhs = new_rhs;
	status = classify_and_catch_special_cases();
      }
      break;
    case HAS_NULL_DECREMENT:
      {
	Expr new_rhs;
	unsigned result = eliminate_null_decrements(recurrence_rhs, new_rhs);
	if (result == 0) {
	  recurrence_rhs_rewritten = true;
	  recurrence_rhs = new_rhs;
	  status = classify_and_catch_special_cases();
	}
	else if (result == 1)
	  status = UNSOLVABLE_RECURRENCE;
	else
	  status = INDETERMINATE_RECURRENCE;
	exit_anyway = true;
      }
      break;
    case MALFORMED_RECURRENCE:
      exit_anyway = true;
      break;

    default:
      throw std::runtime_error("PURRS internal error: "
			       "catch_special_cases().");
      break;
    }
  } while (!exit_anyway && status != SUCCESS);
  return status;
}
