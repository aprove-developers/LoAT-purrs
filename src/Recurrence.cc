/* Recurrence class implementation (non-inline functions).
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

#include "Recurrence.defs.hh"
#include "util.hh"
#include "simplify.hh"
#include <algorithm>
#include <iostream>
#include <fstream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

const PURRS::Symbol&
PURRS::Recurrence::n = Symbol("n");

namespace {
using namespace PURRS;

/*!
  Substitutes every occurrence of the symbol \p s in the expression \p e
  with the expression \p r.
  Returns an expression containing the result of the substitution.
*/
Expr
substitute_symbol_with_expression(const Expr& e, const Symbol& s,
				  const Expr& r) {
  Expr e_substituted;
  if (e.is_a_add()) {
    e_substituted = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_substituted += substitute_symbol_with_expression(e.op(i), s, r);
  }
  else if (e.is_a_mul()) {
    e_substituted = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_substituted *= substitute_symbol_with_expression(e.op(i), s, r);
  }
  else if (e.is_a_power())
    return pwr(substitute_symbol_with_expression(e.arg(0), s, r),
	       substitute_symbol_with_expression(e.arg(1), s, r));
  else if (e.is_a_function())
    if (e.nops() == 1)
      return apply(e.functor(),
		   substitute_symbol_with_expression(e.arg(0), s, r));
    else {
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned j = 0; j < num_argument; ++j)
	argument[j] = substitute_symbol_with_expression(e.arg(j), s, r);
      return apply(e.functor(), argument);
    }
  else if (e == s)
    return r;
  else
    return e;
  return e_substituted;
}

/*!
  Returns the same expression \p term if it does not contain
  the function `x'; returns 0 otherwise.
*/
Expr
find_term_without_function_x(const Expr& term) {
  if (term.is_a_mul()) {
    bool found_initial_condition = false;
    for (unsigned j = term.nops(); j-- > 0; )
      if (term.op(j).is_the_x_function()) {
	found_initial_condition = true;
	break;
      }
    if (!found_initial_condition)
      return term;
  }
  else
    if (!term.is_the_x_function())
      return term;
  return 0;
}

#if 0
/*!
  Given an expression \p e, this function substituted each occurence
  of the function \f$ x(n - d) \f$, with \f$ d \f$ a positive integers,
  with an expression contained in the vector \p repls.
  The position of the \p repls' expression to substitute is given by
  the position of the positive integer \f$ d \f$ in the vector \p positions.
*/
Expr
substitute_function_x(const Expr& e, const std::vector<Expr>& repls,
		      const std::vector<unsigned>& positions) {
  Expr e_substituted;
  if (e.is_a_add()) {
    e_substituted = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_substituted += substitute_function_x(e.op(i), repls, positions);
  }
  else if (e.is_a_mul()) {
    e_substituted = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_substituted *= substitute_function_x(e.op(i), repls, positions);
  }
  else if (e.is_a_power())
    return pwr(substitute_function_x(e.arg(0), repls, positions),
	       substitute_function_x(e.arg(1), repls, positions));
  else if (e.is_a_function())
    if (e.is_the_x_function()) {
      const Expr& argument = e.arg(0);
      // Colud be the case of an initial condition.
      if (argument.nops() == 2) {
	// Finds the positive integer `d' of `x(n-d)'.
	unsigned d;
	if (argument.op(0) == Recurrence:: n)
	  d = argument.op(1).ex_to_number().to_int();
	else
	  d = argument.op(0).ex_to_number().to_int();
	// Finds the position of `d' in the vector `positions'.
	unsigned position = *find(positions.begin(), positions.end(), d);
	return repls[position];
      }
      else
	return e;
    }
    else if (e.nops() == 1)
      return apply(e.functor(),
		   substitute_function_x(e.arg(0), repls, positions));
    else {
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned j = 0; j < num_argument; ++j)
	argument[j] = substitute_function_x(e.arg(j), repls, positions);
      return apply(e.functor(), argument);
    }
  else
    return e;
  return e_substituted;
}
#endif

} // anonymous namespace

/*!
  Consider the right hand side \p rhs of the order \f$ k \f$ recurrence
  relation
  \f$ a_1 * x_{n-1} + a_2 * x_{n-2} + \dotsb + a_k * x_{n-k} + p(n) \f$.
  We try to check that the solution is correct.
  - Validation of initial conditions.
    If <CODE>recurrence_rhs</CODE> is equal to \f$ x(0), \cdots, x(k) \f$ for
    \f$ n = 0, \cdots, k-1 \f$ respectively then
    the initial conditions are verified and we continue to check; otherwise
    return <CODE>false</CODE> because the solution can be wrong or it is not
    simplified enough.
  - Since the initial conditions are verified, we erase from
    <CODE>solution</CODE> all terms containing an initial condition.
    In other words, we check that ther remainder of the solution
    satisfies the same recurrence relation, but with the initial conditions
    all equal to \f$ 0 \f$.
    Starting from the partial solution just computed, we substitute
    \f$ x(n-1) \f$, \f$ x(n-2) \f$, \f$ \dots \f$, \f$ x_{n-k} \f$ into
    <CODE>recurrence_rhs</CODE>.
    We next consider the difference between the partial solution
    and the new right hand side:
    - if it is equal to zero     -> returns <CODE>true</CODE>:
                                    the solution is certainly right.
    - if it is not equal to zero (in a syntactical sense)
                                 -> returns <CODE>false</CODE>:   
				    the solution can be wrong or
				    we failed to simplify it.

  FIXME: In the latter case, we will need more powerful tools to
  decide whether the solution is right or it is really wrong.
*/
Recurrence::VERIFY_STATUS
PURRS::Recurrence::verify_solution() const {
  if (solved || solve()) {
    D_VAR(old_recurrence_rhs);
    D_VAR(recurrence_rhs);
    D_VAR(gcd_decrements_old_rhs);
    D_VAR(order());
    D_VAR(first_initial_condition());
    if (order() == 0)
      return CORRECT;
    else {
      // Step 1: validation of initial conditions.
      for (unsigned i = order(); i-- > 0; ) {
	Expr solution_valuated
	  = substitute_symbol_with_expression(solution, n,
					      first_initial_condition() + i);
	solution_valuated = simplify_numer_denom(solution_valuated);
	D_VAR(solution_valuated);
	if (solution_valuated != x(first_initial_condition() + i))
	  return DONT_KNOW;
      }
      // Step 2: find `partial_solution'.
      // The initial conditions are verified. Build the expression
      // `partial_solution' that has all terms of `solution' minus those
      // containing an initial condition.
      Expr partial_solution = 0;
      if (solution.is_a_add())
	for (unsigned i = solution.nops(); i-- > 0; )
	  partial_solution
	    += find_term_without_function_x(solution.op(i));
      else
	partial_solution = find_term_without_function_x(solution);
      partial_solution = simplify_on_output_ex(partial_solution.expand(),
					       false);
      D_VAR(partial_solution);
      // The recurrence is homogeneous.
      if (partial_solution == 0)
	return CORRECT;
      // Step 3: construct the vector `terms_to_sub': each element of it
      // contains `partial_solution' with `n' substituted by `n - d'
      // (the `d' are the decrements of the terms `x(n - d)').
      // These new expressions contained in the vector `terms_to_sub' are
      // substituted to the correspondenting values in `recurrence_rhs'.
#if 0
      std::vector<unsigned> decrements = tdip->get_decrements();
      std::vector<Expr> terms_to_sub(decrements.size());
      for (unsigned i = decrements.size(); i-- > 0; )
	terms_to_sub[i]
	  = substitute_symbol_with_expression(partial_solution,
					      n, n-decrements[i]);
      Expr substituted_rhs = simplify_on_output_ex(recurrence_rhs.expand(),
						   false);
      substituted_rhs = substitute_function_x(recurrence_rhs, terms_to_sub,
					      decrements);
#else
      Expr substituted_rhs;
      std::vector<Expr> terms_to_sub(order());
      // We have not applied the order reduction in order to solve the
      // recurrence. 
      if (gcd_decrements_old_rhs == 0) {
	substituted_rhs = recurrence_rhs;
	gcd_decrements_old_rhs = 1;
      }
      // We have applied the order reduction in order to solve the recurrence.
      else
	substituted_rhs = old_recurrence_rhs;
      for (unsigned i = order(); i-- > 0; ) {
	terms_to_sub[i]
	  = simplify_all(partial_solution 
			 .subs(n, n - (i + 1) * gcd_decrements_old_rhs));
	substituted_rhs
	  = substituted_rhs.subs(x(n - (i + 1) * gcd_decrements_old_rhs),
				 terms_to_sub[i]);
      }
      D_VEC(terms_to_sub, 0, terms_to_sub.size()-1);
#endif
      D_VAR(substituted_rhs);
      Expr diff = simplify_all(partial_solution - substituted_rhs);
      D_VAR(diff);
      if (!diff.is_zero())
	return DONT_KNOW;
      return CORRECT;
    }
  }
  // We failed to solve the recurrence.
  // If the client still insists in asking for the verification...
  return DONT_KNOW;
}

bool
PURRS::Recurrence::OK() const {
#ifndef NDEBUG
  using std::endl;
  using std::cerr;
#endif

  switch(type) {
  case UNKNOWN:
  case ORDER_ZERO:
    if (tdip != 0) {
#ifndef NDEBUG
      cerr << "Recurrence with type unknown or of order zero!" << endl;
#endif
      return false;
    }
  case LINEAR_FINITE_ORDER_CONST_COEFF:
  case LINEAR_FINITE_ORDER_VAR_COEFF:
    if (tdip == 0) {
#ifndef NDEBUG
      cerr << "Recurrence of finite order!" << endl;
#endif
      return false;
    }
    else {
      //      if (! || !)
    }
  default:
    return true;
  }

  return true;
}

void
PURRS::Recurrence::dump(std::ostream& s) const {
  s << "solved = " << (solved ? "true" : "false") << std::endl;
  s << "recurrence_rhs = " << recurrence_rhs << std::endl;
  s << "auxiliary_definitions:" << std::endl;
  blackboard.dump(s);
  //s << "solution = " << solution << std::endl;
  s << std::endl;
}
