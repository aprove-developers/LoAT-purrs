/* Definition of the main recurrence relation solver.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#include "util.hh"
#include "simplify.hh"
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
  if (a == Recurrence::n && b.is_a_number(d)
      && d.is_rational() && !d.is_integer() && d.is_positive()) {
    divisor = 1 / d;
    return true;
  }
  else if (b == Recurrence::n && a.is_a_number(d)
	   && d.is_rational() && !d.is_integer()&& d.is_positive()) {
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
	return PURRS::apply(e.functor(), e.arg(0).substitute(s, r));
      else
	return PURRS::apply(e.functor(),
		     change_variable_function_x(e.arg(0), s, r));
    else {
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int j = 0; j < num_argument; ++j)
	argument[j] = change_variable_function_x(e.arg(j), s, r);
      return PURRS::apply(e.functor(), argument);
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
  if (dec >= max_decrement) {
    coefficient = possibly_coeff;
    if (dec > max_decrement) {
      max_decrement = dec;
    }
  }
}

void
find_max_decrement_and_coeff_factor(const Expr& e,
				    int& max_decrement, Expr& coefficient) {
  Expr coeff = coefficient;
  int max_dec = max_decrement;
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
				       max_dec, coeff);
    }
    if (max_dec > max_decrement)
      coefficient = coeff;
    else
      coefficient += coeff;
  }
  else {
    if (e.is_the_x_function())
      assign_max_decrement_and_coeff(e.arg(0), 1, max_dec, coeff);
    coefficient = coeff;
  }

  max_decrement = max_dec;
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
  the recurrence relation, but transforms it into its
  <EM>standard form</EM> \f$ x(n) = new_rhs \f$, where \f$ new_rhs \f$
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
#if DEBUG
  DD_MSGVAR("Eliminate negative decrements in ", rhs);
  DD_MSGVAR("Max decrement: ", max_decrement);
  DD_MSGVAR("Coefficient: ", coefficient);
#endif

  // We wish to transform `x(n) = a_1 * x(n+1) + ... + a_m * x(n+m)'
  // into `x(t) = (-x(t-m) + a_1 * x (t-m+1) + ... ) / (-a_m)' where
  // `t=n+m' and the coefficients a_m may depend on n.
  // If `x(n)' appears in the rhs, if should have already been removed
  // by now.

  // 1: Let `m = max_decrement'. Put `t = n + m'. (For efficiency,
  //    actually replace `n' with `n+m' everywhere).
  const Symbol& n = Recurrence::n;
  new_rhs = rhs.substitute(n, n - max_decrement);
  coefficient = coefficient.substitute(n, n - max_decrement);

  // 2: Move the old `x(n)' (now `x(n-m)') to the right hand side.
  new_rhs -= x(n - max_decrement);

  // 3: Divide everything by the leading coefficient.
  new_rhs /= -coefficient;

  // 4: The coefficient of `x(n)' (the old `x(n+m)') is now -1. Move it to
  //    the left hand side (i.e., remove it from the right hand side). For
  //    safety, explicitly add `x(n)' to the rhs to accomplish this.
  // new_rhs = new_rhs.substitute(x(n), 0);
  new_rhs = new_rhs.expand().collect(Expr_List(x(n))) + x(n);

  // FIXME: Check that `x(n)' no longer appears in the rhs.
  new_rhs = new_rhs.expand();
#if DEBUG
  DD_MSGVAR("Rewritten recurrence: ", new_rhs);
#endif

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

//! Kinds of non-linear term.
enum Kinds_Of_Term {
  //! \brief
  //! It is a not legal non-linear term of the form \f$ x(x(a)) \f$ with
  //! \f$ a \f$ containing the special symbol \f$ n \f$; 
  NESTED,

  //! \brief
  //! It is a legal non-linear term:
  //! -  \f$ x(n) = a x(n-k_1)^b_1 \cdots x(n-k_h)^b_h \f$,
  //!    where \f$ k_1, \dots, k_h, b_1, \dots, b_h \f$ are positive integers
  //!    and \f$ h > 1 \f$;
  //! -  \f$ x(n) = x(n-k)^b \f$,
  //!    where \f$ b \f$ is a rational number while \f$ k \f$ is a
  //!    positive integers.
  //! The two cases above hold also if instead of terms like `x(n-k)'
  //! there are term like `x(n/k)'.
  LEGAL_NON_LINEAR_TERM,

  //! It is a non-linear term of the form \f$ x(n)^b \f$ (\f$ b \neq 1 \f$).
  INDETERMINATE,

  //! It is a non-linear term of the form \f$ b^x(n) \f$ (\f$ b \neq 1 \f$).
  UNSOLVABLE,

  //! It is a linear term.
  LINEAR_TERM
};

/*!
  - If there is in \p e an \f$ x \f$ function (with the argument containing
    \f$ n - d \f$, \f$ d > 1 \f$) inside an other function different from
    \f$ x \f$ function or if there is in \p e a power with an \f$ x \f$
    function in the base or in the exponent of it then returns
    <CODE>LEGAL_NON_LINEAR_TERM</CODE>;
  - if there is in \p e an \f$ x \f$ function (with the argument containing
    \f$ n \f$) inside an other \f$ x \f$ function return <CODE>NESTED</CODE>;
  - in all the other cases returns <CODE>LINEAR_TERM</CODE>.
*/
Kinds_Of_Term
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
	  return NESTED;
	else
	  return LEGAL_NON_LINEAR_TERM;
    }
  }
  // There is a power with an `x' function in the base or in the exponent.
  else if (e.is_a_power()) {
    const Expr& base = e.arg(0);
    const Expr& exponent = e.arg(1);
    if (base.has_x_function(Recurrence::n))
      if (base == x(Recurrence::n))
	return INDETERMINATE;
      else
	return LEGAL_NON_LINEAR_TERM;
    if (exponent.has_x_function(Recurrence::n))
      if (exponent == x(Recurrence::n))
	return UNSOLVABLE;
      else
	return LEGAL_NON_LINEAR_TERM;
  }
  return LINEAR_TERM;
}

Kinds_Of_Term
find_non_linear_term(const Expr& e) {
  assert(!e.is_a_add());
  // Even if we find a `legal' (i.e. not containing two nested `x' function)
  // non-linear term we do not exit from this function because we must be
  // sure that there are not non-linear term not legal.
  // If `term' at the end of this function is again `LINEAR_TERM' means
  // that we have not find non-linear terms.
  Kinds_Of_Term term = LINEAR_TERM;
  unsigned int num_factors = e.is_a_mul() ? e.nops() : 1;
  if (num_factors == 1) {
    Kinds_Of_Term tmp = x_function_in_powers_or_functions(e);
    // We have found a nested `x' function (not legal non-linear term) or
    // a legal non-linear term.
    if (tmp != LINEAR_TERM)
      return tmp;
  }
  else {
    bool found_function_x = false;
    for (unsigned int j = num_factors; j-- > 0; ) {
      const Expr& factor = e.op(j);
      Kinds_Of_Term tmp = x_function_in_powers_or_functions(factor);
      // Nested `x' function: not legal non-linear term.
      if (tmp == NESTED)
	return tmp;
      // We have found legal non-linear term.
      else if (tmp == LEGAL_NON_LINEAR_TERM)
	term = tmp;
      if (factor.is_the_x_function())
	// We have found legal non-linear term product of two ore more
	// `x' functions.
	if (found_function_x && factor.arg(0).has(Recurrence::n))
	  term = LEGAL_NON_LINEAR_TERM;
	else
	  if (factor.arg(0).has(Recurrence::n))
	    found_function_x = true;
    }
  }
  return term;
}

//! Kinds of non-linear recurrence.
enum Kinds_Of_Non_Linear_Rec {
  //! Recurrence that the system is able to solve.
  OK_NON_LINEAR_REC,

  //! \brief
  //! Wrong recurrence (e.g. \f$ x(n) = c x(n-1)^a \f$, with
  //! \f$ a \in \Qset \setminus \Zset \f$, \f$ c < 0 \f$).
  DOMAIN_ERROR_NON_LINEAR_REC,

  //! Recurrence that the system is not able to solve.
  TOO_COMPLEX_NON_LINEAR_REC
};

//! \brief
//! Returns <CODE>true</CODE> if the non linear recurrence \f$ x(n) = rhs \f$
//! is rewritable as linear recurrence. Returns <CODE>false</CODE> otherwise.
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
  in the case of simple non-linear recurrence of the form
  \f$ x(n) = c x(n-1)^{\alpha} \f$ it contains the pair \f$ c, \alpha \f$
  (in this case \p new_rhs will be \f$ 0 \f$);
  in all the other cases the numeric value of the pair holds \f$ 0 \f$
  while the second element contains the value that will be the
  logarithm's base or the exponential's base used in the rewriting
  of the non-linear recurrence in the correspondent linear recurrence
  (which will be stored in \p new_rhs).
*/
Kinds_Of_Non_Linear_Rec
rewrite_non_linear_recurrence(const Recurrence& rec, const Expr& rhs,
			      Expr& new_rhs,
			      std::pair<Number, Expr>& coeff_and_base,
			      std::vector<Symbol>& auxiliary_symbols,
			      index_type& first_valid_index) {
  Number z;
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
      else
	if (!find_domain_in_N(factor, Recurrence::n, z))
	  return TOO_COMPLEX_NON_LINEAR_REC;
    }
    first_valid_index = z.to_unsigned_int();

    // Consider the special case `x(n) = c x(n-1)^a' with `a' and `c'
    // constants (`a != 1').
    // In this case is not necessary the `linearization' of the recurrence
    // because we already know the solution in function of `c' and `a'.
    // The linear recurrence is necessary if the solution must be evaluated
    // on initial conditions with indeces bigger than the indeces of the
    // symbolic initial conditions occurring in the solution.
    if (simple_cases && rhs.nops() == 2) {
      Number coeff;
      if (rhs.op(0).is_a_number(coeff) && rhs.op(1).is_a_power()
	  && rhs.op(1).arg(0).is_the_x_function()
	  && rhs.op(1).arg(0).arg(0) == Recurrence::n-1) {
	// Store `c'.
	coeff_and_base.first = coeff;
	// Store `a'.
	assert(rhs.op(1).arg(1).is_a_number());
	coeff_and_base.second = rhs.op(1).arg(1);
	if (coeff.is_negative()
	    && !rhs.op(1).arg(1).ex_to_number().is_integer())
	  return DOMAIN_ERROR_NON_LINEAR_REC;
	// Store the right hand side of the associated linear recurrence
	// (when possible).
	if (coeff_and_base.second != -1) {
	  Expr tmp = coeff_and_base.first;
	  if (coeff_and_base.first.is_negative()) {
	    Symbol s = rec.insert_auxiliary_definition(coeff_and_base.first);
	    auxiliary_symbols.push_back(s);
	    tmp = s;
	  }
	  new_rhs = coeff_and_base.second * x(Recurrence::n-1)
	    + log(tmp) / log(abs(rhs.op(1).arg(1).ex_to_number()));
	}
	return OK_NON_LINEAR_REC;
      }
      else if (rhs.op(1).is_a_number(coeff) && rhs.op(0).is_a_power()
	       && rhs.op(0).arg(0).is_the_x_function()
	       && rhs.op(0).arg(0).arg(0) == Recurrence::n-1) {
	// Store `c'.
	coeff_and_base.first = coeff;
	// Store `a'.
	assert(rhs.op(0).arg(1).is_a_number());
	coeff_and_base.second = rhs.op(0).arg(1);
	if (coeff.is_negative()
	    && !rhs.op(0).arg(1).ex_to_number().is_integer())
	  return DOMAIN_ERROR_NON_LINEAR_REC;
	// Store the right hand side of the associated linear recurrence
	// (when possible).
	if (coeff_and_base.second != -1) {
	  Expr tmp = coeff_and_base.first;
	  if (coeff_and_base.first.is_negative()) {
	    Symbol s = rec.insert_auxiliary_definition(coeff_and_base.first);
	    auxiliary_symbols.push_back(s);
	    tmp = s;
	  }
	  new_rhs = coeff_and_base.second * x(Recurrence::n-1)
	    + log(tmp) / log(abs(rhs.op(0).arg(1).ex_to_number()));
	}
	D_VAR(new_rhs);
	return OK_NON_LINEAR_REC;
      }
    }
    new_rhs = 0;
    if (simple_cases)
      if (common_exponent == 1) {
	// Substitute any function `x()' with `exp(1)^{x()}'.
	coeff_and_base.second = Napier;
	Expr tmp = substitute_x_function(rhs, Napier, EXPONENT);
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
	
	return OK_NON_LINEAR_REC;
      }
      else {
	// Substitute any function `x()' with `common_exponent^{x()}'.
	coeff_and_base.second = common_exponent;
	Expr tmp = substitute_x_function(rhs, abs(common_exponent), EXPONENT);
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
	return OK_NON_LINEAR_REC;
      }
  }
  // Second case.
  else if (rhs.is_a_power()) {
    Number num_exp;
    if (rhs.arg(0).is_the_x_function() && rhs.arg(1).is_a_number(num_exp)) {
      coeff_and_base.second = num_exp;
      const Expr& argument = rhs.arg(0).arg(0);
      Number d;
      if (get_constant_decrement(argument, d)) {
	if (d == 1)
	  coeff_and_base.first = 1;
	new_rhs = num_exp * x(Recurrence::n-d);
      }
      else if (get_constant_divisor(argument, d))
	new_rhs = num_exp * x(Recurrence::n/d);   
      else
	return TOO_COMPLEX_NON_LINEAR_REC;
      return OK_NON_LINEAR_REC;
    }
    else
      return TOO_COMPLEX_NON_LINEAR_REC;
  }
  return TOO_COMPLEX_NON_LINEAR_REC;
}

//! \brief
//! Returns <CODE>true</CODE> if \p rhs is the right-hand side of a
//! weighted-average recurrence or it is possible to rewrite it in
//! the form of weighted-average recurrence
//! \f[
//!   x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n).
//! \f]
//! Returns <CODE>false</CODE> otherwise.
/*!
  \param rhs                     the expression that we want to rewrite,
                                 if possible, in the right-hand side of
				 a weighted-average recurrence.
  \param term_sum                the term of \p rhs containing the sum.  
  \param weight                  the weight \f$ f(n) \f$ of the
                                 weighted-average recurrence obtained
				 rewriting, if necessary, \p rhs.
  \param original_weight         the weight of the right-hand side of the
                                 recurrence, initially stored in \p rhs,
				 before the rewriting so that the lower
				 limit of the sum is \f$ 0 \f$, but after
				 the rewriting so that the upper limit of
				 the sum is \f$ n-1 \f$.
  \param inhomogeneous           the inhomogeneous term \f$ g(n) \f$ of the
                                 weighted-average recurrence obtained
				 rewriting, if necessary, \p rhs.
  \param original_inhomogeneous  the inhomogeneous term  of the
                                 recurrence, initially stored in \p rhs,
				 before the rewriting so that the lower
				 limit of the sum is \f$ 0 \f$, but after
				 the rewriting so that the upper limit of
				 the sum is \f$ n-1 \f$.
  \param rhs_first_order         the right hand side of the first order
                                 recurrence associated to the weighted-average
				 recurrence.
  \param first_valid_index       see
                                 \ref first_valid_index "first_valid_index".
  \param domain_problem          <CODE>true</CODE> if
                                 <CODE>first_valid_index</CODE> is largest
				 than the lower limit of the sum;
				 <CODE>false</CODE> otherwise.
  \param rewritten               <CODE>true</CODE> if \p rhs is not in the
                                 form of weighted-average recurrence
				 and this function rewrite it in a
				 weighted-average recurrence.

  \return                        <CODE>true</CODE> if \p rhs is a
                                 weighted-average recurrence or
				 is possible to rewrite it in a
				 weighted-average recurrence;
				 returns <CODE>false</CODE> otherwise.

  The system is able to compute weighted-average recurrence
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

  Moreover, this function rewrite a recurrence, when possible, in the form
  of weighted-average.
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
  Note that there is an important difference between the two rewritings:
  the rewriting that works on the upper limit of the sum transform the
  recurrence in normal form so that the right-hand side does not
  contain terms \f$ x(n) \f$; the rewriting which works on the lower limit
  of the sum is applied so that in the following computations will be
  possible to work at the same mode on each weighted-average recurrence.
  In the last case, once we have found the solution of the rewritten
  recurrence, we have to come back to the solution of the original
  recurrence.

  In conclusion the algorithm for classifying and for finding the solution
  of the weighted-average recurrence previews the following steps:
  - possible rewriting in the form of weighted-average recurrence;
  - computation of the right hand side of the associated first order
    recurrence;
  - shift forward of the first order recurence: \f$ n = n + 1 \f$;
  - computation of the solution of the first order recurrence;
  - shift backward of the solution: \f$ n = n - 1 \f$;
  - substitution of the initail condition \f$ x(1) = f(1) x(0) + g(1) \f$.
  - Another possible shift if the original recurrence has had
    the lower limit of the sum different from \f$ 0 \f$ in order
    to come back to the solution of the original recurrence.
*/
bool
rewrite_weighted_average_recurrence(const Expr& rhs, const Expr& term_sum,
				    Expr& weight, Expr& original_weight,
				    Expr& inhomogeneous,
				    Expr& original_inhomogeneous,
				    Expr& rhs_first_order,
				    index_type& first_valid_index,
				    bool& domain_problem,
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

  // Assert that the argument of the sum is linear in x(k): it cannot be
  // a weighted-average recurrence if it is not linear.
  Symbol x_k;
  Expr summand = term_sum.arg(3).substitute(x(term_sum.arg(0)), x_k).expand();
  if (!summand.is_polynomial(x_k) || summand.degree(x_k) > 1)
    return false;

  // Find the weight `f(n)' and the inhomogeneous term of the
  // recurrence transformed so that to have the upper limit of the
  // sum equal to `n'.
  if (upper == Recurrence::n) {
    if (weight != 1) {
      const Expr& tmp = 1 - weight;
      weight /= tmp;
      inhomogeneous /= tmp;
    }
    else {
      weight *= -1;
      lower += 1;
      inhomogeneous = -inhomogeneous.substitute(Recurrence::n,
						Recurrence::n + 1) ;
    }
    rewritten = true;
  }
  original_weight = weight;
  original_inhomogeneous = inhomogeneous;

  Number z = 0;
  // Find the largest positive or null integer that cancel the numerator of
  // `f(n)' and store it in `z' if it is bigger than the current `z'.
  if (!find_domain_in_N(weight.numerator(), Recurrence::n, z))
    return false;
  // `z' will contain the largest positive or null integer, if it exists,
  // that cancel the denominator of `f(n)'.
  // If this integer does not exist then `z' is left to 0.
  if (!find_domain_in_N(weight.denominator(), Recurrence::n, z))
    return false;

  if (!find_domain_in_N(inhomogeneous.denominator(), Recurrence::n, z))
    return false;
  first_valid_index = z.to_unsigned_int();
  
  // Find the weight `f(n)' and the inhomogeneous term of the
  // recurrence transformed so that to have the lower limit of the
  // sum equal to `0'.
  if (lower > 0) {
    rewritten = true;
    weight = weight.substitute(Recurrence::n, Recurrence::n + lower);
    inhomogeneous = inhomogeneous.substitute(Recurrence::n,
					     Recurrence::n + lower);
  }
  if (int(first_valid_index) > lower)
    domain_problem = true;
  
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
    return CL_HAS_NEGATIVE_DECREMENT;
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
  // `term == 0' or `term == 1' indicate
  // two different cases of non-linearity.
  unsigned int term = find_non_linear_term(addend);
  if (term == LEGAL_NON_LINEAR_TERM) {
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
    index_type first_valid_index = 0;
    switch (rewrite_non_linear_recurrence(*this, rhs, new_rhs, coeff_and_base,
					  auxiliary_symbols,
					  first_valid_index)) {
    case OK_NON_LINEAR_REC:
      set_non_linear_finite_order();
      D_VAR(new_rhs);
      non_linear_p = new Non_Linear_Info(Recurrence(new_rhs), coeff_and_base,
					 auxiliary_symbols);
      set_first_valid_index(first_valid_index);
      return CL_SUCCESS;
    case DOMAIN_ERROR_NON_LINEAR_REC:
      return CL_DOMAIN_ERROR;
    case TOO_COMPLEX_NON_LINEAR_REC:
      return CL_TOO_COMPLEX;
    default:
      throw std::runtime_error("PURRS internal error: "
			       "classification_summand().");
    }
  }
  // This is the case of nested `x' function with argument containing `n'.
  else if (term == NESTED)
    return CL_MALFORMED_RECURRENCE;
  else if (term == INDETERMINATE)
    return CL_INDETERMINATE_RECURRENCE;
  else if (term == UNSOLVABLE)
    return CL_UNSOLVABLE_RECURRENCE;
  

  unsigned int num_factors = addend.is_a_mul() ? addend.nops() : 1;
  Number num;
  if (num_factors == 1)
    if (addend.has_non_rational_numbers())
      return CL_MALFORMED_RECURRENCE;
    else if (addend.is_the_x_function()) {
      const Expr& argument = addend.arg(0);
      if (argument == n)
	return CL_HAS_NULL_DECREMENT;
      else if (!argument.has(n))
	return CL_MALFORMED_RECURRENCE;
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
	// the weighted-average recurrence
	// `x(n) = f(n) sum(k, n_0, u(n), x(k)) + g(n)'.
	Expr orig_weight;
	const Expr& rhs_rewritten = rhs.collect_term(addend, orig_weight);
	Expr rhs_first_order;
	bool rewritten = false;
	index_type first_valid_index;
	bool domain_problem = false;
	Expr weight = orig_weight;
	Expr orig_inhomogeneous;
	// If it is possible, rewrites the recurrence in
	// the form of weighted-average recurrence.
	if (rewrite_weighted_average_recurrence(rhs_rewritten, addend,
						weight, orig_weight,
						inhomogeneous,
						orig_inhomogeneous,
						rhs_first_order,
						first_valid_index,
						domain_problem, rewritten)) { 
	  if (domain_problem)
	    return CL_DOMAIN_ERROR;
	  weighted_average_p
	    = new Weighted_Average_Info(Recurrence(rhs_first_order), weight);
	  set_weighted_average();
	  set_first_valid_index(std::max(first_valid_index,
					 addend.arg(1).ex_to_number()
					 .to_unsigned_int()));
	  if (rewritten) {
	    bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	    rec_rewritten = true;
	    // To save these informations is important in order to verify
	    // the solution (or a bound) of the original recurrence.
	    set_original_rhs(orig_weight, orig_inhomogeneous,
			     addend.arg(1).ex_to_number().to_unsigned_int(),n);
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
	  return CL_HAS_NULL_DECREMENT;
	else if (!argument.has(n))
	  return CL_MALFORMED_RECURRENCE;
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
	  // `f(n)' of the weighted-average recurrence
	  // `x(n) = f(n) sum(k, n_0, u(n), x(k)) + g(n)'.
	  Expr orig_weight;
	  const Expr& rhs_rewritten = rhs.collect_term(factor, orig_weight);
	  // There are not other sums equal to `factor'.
	  if (orig_weight == 1)
	    for (unsigned int j = num_factors; j-- > 0; )
	      if (addend.op(j) != factor)
		orig_weight *= addend.op(j);
	  Expr rhs_first_order;
	  bool rewritten = false;
	  index_type first_valid_index;
	  bool domain_problem = false;
	  Expr weight = orig_weight;
	  Expr orig_inhomogeneous;
	  // If it is possible, rewrites the recurrence in
	  // the form of weighted-average recurrence.
	  if (rewrite_weighted_average_recurrence(rhs_rewritten, factor,
						  weight, orig_weight,
						  inhomogeneous,
						  orig_inhomogeneous,
						  rhs_first_order,
						  first_valid_index,
						  domain_problem,
						  rewritten)) { 
	    if (domain_problem)
	      return CL_DOMAIN_ERROR;
	    weighted_average_p
	      = new Weighted_Average_Info(Recurrence(rhs_first_order), weight);
	    set_weighted_average();
	    set_first_valid_index(std::max(first_valid_index,
					   factor.arg(1).ex_to_number()
					   .to_unsigned_int()));
	    if (rewritten) {
	      bool& rec_rewritten = const_cast<bool&>(recurrence_rewritten);
	      rec_rewritten = true;
	      // To save these informations is important in order to verify
	      // the solution (or a bound) of the original recurrence.
	      set_original_rhs(orig_weight, orig_inhomogeneous,
			       factor.arg(1).ex_to_number().to_unsigned_int(),n);
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
  - <CODE>CL_SUCCESS</CODE>                   if the recurrence is linear of
                                              finite order; weighted-average;
					      non-linear of finite order; in
					      the case of functional equation;
  - <CODE>CL_HAS_NON_INTEGER_DECREMENT</CODE> if the right-hand side of the
                                              recurrence contains at least an
					      occurrence of \f$ x(n-k) \f$,
					      where \f$ k \f$ is not an
					      integer;
  - <CODE>CL_HAS_NEGATIVE_DECREMENT</CODE>    if the right-hand side of the
                                              recurrence contains at least an
					      occurrence of \f$ x(n-k) \f$,
					      where \f$ k \f$ is a negative
					      integer;
  - <CODE>CL_HAS_NULL_DECREMENT</CODE>        if the right-hand side of the
                                              recurrence contains at least an
					      occurrence of \f$ x(n) \f$;
  - <CODE>CL_CL_HAS_HUGE_DECREMENT</CODE>     if the right-hand side of the
                                              recurrence contains at least an
					      occurrence of \f$ x(n-k) \f$,
					      where \f$ k \f$ is too big to
					      be handled by the standard
					      solution techniques;
  - <CODE>CL_MALFORMED_RECURRENCE</CODE>      if the recurrence is not
                                              \ref syntactically_correct "syntactically correct";
  - <CODE>CL_DOMAIN_ERROR</CODE>              if the recurrence is not
                                              well-defined;
  - <CODE>CL_UNSOLVABLE_RECURRENCE</CODE>     if the recurrence is not
                                              solvable;
  - <CODE>CL_INDETERMINATE_RECURRENCE</CODE>  if the recurrence is
                                              indeterminate, hence it has
					      infinite solutions.
  - <CODE>CL_TOO_COMPLEX</CODE>               in all the other cases.

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

  For the <EM>weighted average</EM> recurrences there is not
  a general method to solve it. At the moment we solve weighted average
  recurrences of the shape
  \f[
    T(n) = f(n) \sum_{k=0}^{n-1} T(k) + g(n)
  \f]
  transforming it in the linear recurrence of first order
  \f[
    T(n) = \frac{f(n)}{f(n-1)} (1+f(n-1)) T(n-1)
      + f(n) \left( \frac{g(n)}{f(n)} - \frac{g(n-1)}{f(n-1)} \right).
  \f].
  In the structure <CODE>Weighted_Average_Info</CODE> are stored the right
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
  // Simplifies expanded expressions, in particular simplifies
  // expression containing logarithms and rewrites nested powers.
  rhs = simplify_logarithm(recurrence_rhs);
  rhs = simplify_ex_for_input(rhs, true);
  // Splits the sum in as many sums as the addends of the summand
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

  // We also accept recurrences defined as 
  // x(n)=max(f(x(0),...,x(n-1)), g(x(0),...,x(n-1))). As this classifier
  // will we invoked recursively in this case, we can return immediately.
  if (rhs.is_the_max_function()) {
    type_ = MAX_FUNCTION;
    return CL_SUCCESS;
  }

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
      // As soon as the system notices the this is a weighted-average
      // recurrences, stops the classification because already
      // have all necessary information in order to continue.
      if (is_weighted_average())
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

  // Find the least non-negative integer `z' such that the
  // recurrence is well-defined for `n >= z'.
  // The function used in order to find `z' works like the
  // computation of a maximum: if the system is able to find
  // `z' and its value is bigger than the previous, updates `z';
  // otherwise `z' is left unchanged.
  Number z = 0;
  if (is_linear_finite_order()) {
    // In the case of order equal to `0' we do not perform any
    // check because the solution is simply the right-hand side
    // of the recurrence.
    if (order > 0) {
      // FIXME: check also the numerator (ex. log(n-2))!!!
      // Check the denominator of the inhomogeneous term.
      if (!inhomogeneous.is_zero()) {
	const Expr& denom_inhomogeneous = inhomogeneous.denominator();
	if (has_parameters(denom_inhomogeneous))
	  return CL_TOO_COMPLEX;
	if (!find_domain_in_N(denom_inhomogeneous, n, z))
	  return CL_TOO_COMPLEX;
      }
      if (is_linear_finite_order_var_coeff() && order == 1) {
	// Check the coefficient.
	const Expr& coefficient = coefficients[1];
	if (has_parameters(coefficient))
	  return CL_TOO_COMPLEX;
	if (!find_domain_in_N(coefficient, n, z))
	  return CL_TOO_COMPLEX;
      }
    }
    finite_order_p = new Finite_Order_Info(order, coefficients,
					   gcd_among_decrements);
    set_first_valid_index(z.to_unsigned_int());
  }
  else if (is_functional_equation()) {
    for (std::map<Number, Expr>::const_iterator i = homogeneous_terms.begin(),
	   homogeneous_terms_end = homogeneous_terms.end();
	 i != homogeneous_terms_end; ++i)
      // Check the coefficients `a' of the homogeneous terms `a x(n/b)'
      // in the recurrence. 
      if (!find_domain_in_N(i->second, n, z))
	return CL_TOO_COMPLEX;
    functional_eq_p = new Functional_Equation_Info(homogeneous_terms);
    set_first_valid_index(std::max(z.to_int(), 1));
  }
  // In the case of non linear recurrence or weighted-average recurrence
  // we have already done the operation `new ...' and we have already
  // computed `first_valid_index'.

  assert(is_linear_finite_order() || is_functional_equation()
	 || is_non_linear_finite_order() || is_weighted_average());

  return CL_SUCCESS;
}

/*!
  Classifies the recurrence \p *this calling the method
  <CODE>classify()</CODE>.
  If the function <CODE>classify()</CODE> returns the value
  <CODE>CL_HAS_NEGATIVE_DECREMENT</CODE> or the value
  <CODE>CL_HAS_NULL_DECREMENT</CODE>, this function tries to rewrite the
  recurrence in the normal form \f$ x(n) = r \f$, where \f$ r \f$
  does not contain terms \f$ x(n-k) \f$, where \f$ k \f$
  is not a positive integer.

  Returns:
  - <CODE>CL_SUCCESS</CODE>                   if the recurrence is linear of
                                              finite; weighted-average;
					      non-linear of finite order; in
					      the case of functional equation;
  - <CODE>CL_HAS_NON_INTEGER_DECREMENT</CODE> if the right-hand side of the
                                              recurrence contains at least an
					      occurrence of \f$ x(n-k) \f$,
					      where \f$ k \f$ is not an
					      integer;
  - <CODE>CL_CL_HAS_HUGE_DECREMENT</CODE>     if the right-hand side of the
                                              recurrence contains at least an
					      occurrence of \f$ x(n-k) \f$,
					      where \f$ k \f$ is too big
					      to be handled by the standard
					      solution techniques;
  - <CODE>CL_MALFORMED_RECURRENCE</CODE>      if the recurrence is not
                                              \ref syntactically_correct "syntactically correct";
  - <CODE>CL_DOMAIN_ERROR</CODE>              if the recurrence is not
                                              well_defined;
  - <CODE>CL_UNSOLVABLE_RECURRENCE</CODE>     if the recurrence is not
                                              solvable;
  - <CODE>CL_INDETERMINATE_RECURRENCE</CODE>  if the recurrence is
                                              indeterminate, hence it has
					      infinite solutions.
  - <CODE>CL_TOO_COMPLEX</CODE>               in all the other cases. 
  
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
    case CL_UNSOLVABLE_RECURRENCE:
    case CL_INDETERMINATE_RECURRENCE:
    case CL_DOMAIN_ERROR:
    case CL_TOO_COMPLEX:
      exit_anyway = true;
      break;
    case CL_HAS_NEGATIVE_DECREMENT:
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
    case CL_HAS_NULL_DECREMENT:
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
