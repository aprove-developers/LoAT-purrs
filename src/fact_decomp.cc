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

#include "fact_decomp.hh"
#include "simplify.hh"
#include "util.hh"
#include "Expr.defs.hh"
#include "Symbol.types.hh"
#include "Number.types.hh"
#include <vector>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;

void
factorial_decomposition_argument(const Expr& e, const Symbol& x,
				 const Expr& coeff_factorial,
				 std::vector<Number>& arg,
				 std::vector<Expr>& coeff) {
  // `e' is the argument of the factorial function.
  assert(!e.is_a_add());
  unsigned int arg_size = arg.size();
  unsigned int position = arg_size;
  unsigned int num_factors = e.is_a_mul() ? e.nops() : 1;
  if (num_factors == 1) {
    bool found = false;
    for (unsigned int i = arg_size; i-- > 0; ) {
      Expr e_minus_x = e;
      if (e == x)
	e_minus_x = 1;
      if (arg[i] == e_minus_x) {
	position = i;
	found = true;
	break;
      }
    }
    if (!found) {
      if (e == x)
	arg.push_back(1);
      else
	arg.push_back(0);
      coeff.push_back(0);
    }
    if (e == x)
      coeff[position] += coeff_factorial;
    else
      coeff[position] += factorial(coeff_factorial);
  }
  else {
    bool has_x = false;
    Expr arg_not_x = 1;
    for (unsigned int i = e.nops(); i-- > 0; ) {
      const Expr& factor = e.op(i);
      if (factor == x)
	has_x = true;
      else
	arg_not_x *= factor;
    }
    bool found = false;
    for (unsigned int i = arg_size; i-- > 0; )
      if (has_x && arg_not_x.is_a_number()) {
	if (arg[i] == arg_not_x) {
	  position = i;
	  found = true;
	break;
	}
      }
      else
	if (arg[i] == 0) {
	  position = i;
	  found = true;
	  break;
	}
    if (!found)
      if (has_x && arg_not_x.is_a_number()) {
	arg.push_back(arg_not_x.ex_to_number());
	coeff.push_back(0);
	coeff[position] += coeff_factorial;	
      }
      else {
	arg.push_back(0);
	coeff.push_back(0);
	coeff[position] += coeff_factorial * factorial(arg_not_x);
      }
  }
}

void
factorial_decomposition_summand(const Expr& e, const Symbol& x,
				std::vector<Number>& arg,
				std::vector<Expr>& coeff) {
  unsigned int arg_size = arg.size();
  unsigned int position = arg_size;
  unsigned int num_factors = e.is_a_mul() ? e.nops() : 1;
  if (num_factors == 1) {
    if (e.is_the_factorial_function() && !e.arg(0).is_a_add())
      factorial_decomposition_argument(e.arg(0), x, 1, arg, coeff);
    else {
      bool found = false;
      for (unsigned int i = arg_size; i-- > 0; )
	if (arg[i] == 0) {
	  position = i;
	  found = true;
	  break;
	}
      if (!found) {
	arg.push_back(0);
	coeff.push_back(0);
      }
      coeff[position] += e;
    }
  }
  else {
    Expr argument = 0;
    Expr coefficient = 1;
    for (unsigned int i = num_factors; i-- > 0; ) {
      const Expr& factor = e.op(i);
      if (factor.is_the_factorial_function()) {
	if (!argument.is_zero())
	  // FIXME: Product of factorials.
	  throw "PURRS internal error: product of factorials\n"
	    "in `factorial_decomposition()'.";
	argument = factor.arg(0);
      }
      else
	coefficient *= factor;
    }
    if (!argument.is_zero())
      factorial_decomposition_argument(argument, x, coefficient, arg, coeff);
    else {
      bool found = false;
      for (unsigned int i = arg_size; i-- > 0; )
	if (arg[i] == 0) {
	  position = i;
	  found = true;
	  break;
	}
      if (!found) {
	arg.push_back(0);
	coeff.push_back(0);
      }
      coeff[position] += coefficient;
    }
  }
}

} // anonymous namespace

/*!
  Let \f$ e(x) \f$ be the expression in \p x contained in \p e,
  which is assumed to be already expanded.
  This function first applies the rewrite rule
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
  and then computes a decomposition
  \f[
    e(x) = \sum_{j=1}^r h_j(x) (a_j x)!
  \f]
  where
  \f[
    h_j \in
    \{
      \fund{f}{\Nset}{\Cset}
      \itc \exists r \in \Nset \st \forall j = 1, \dots r
      \itc \exists q_j \in \Cset[x]
      \st \exists \alpha_j \in \Cset \setminus \{0\}
      \st e(x) = \sum_{j=1}^r q_j(x) \alpha_j^x.
    \}
  \f]
  
  The expressions corresponding to \f$ a_j \f$ and \f$ h_j(x) \f$
  are stored in the \f$ i \f$-th position of the vectors \p arg and 
  \p coeff, respectively.
*/
void
PURRS::factorial_decomposition(const Expr& e, const Symbol& x,
			       std::vector<Number>& arg,
			       std::vector<Expr>& coeff) {
  assert(e.is_expanded());
  // This simplification is necessary because rewrite eventual
  // nested powers.
  Expr e_simpl = simplify_ex_for_input(e, true);
  // With the logarithms' simplification non-polynomial expressions
  // can become polynomial expressions for the rule `log(a^b) = b log(a)'.
  e_simpl = simplify_logarithm(e_simpl);
  e_simpl = simplify_binomials_factorials_exponentials(e_simpl);
  
  unsigned int num_summands = e_simpl.is_a_add() ? e_simpl.nops() : 1;
  // An upper bound to the number of exponentials is the number of
  // summands in `e': reserve space in the output vectors so that
  // no reallocations will be required.
  arg.reserve(num_summands);
  coeff.reserve(num_summands);
  if (num_summands > 1)
    for (unsigned int i = num_summands; i-- > 0; )
      factorial_decomposition_summand(e_simpl.op(i), x, arg, coeff);
  else
    factorial_decomposition_summand(e_simpl, x, arg, coeff);
  D_VEC(arg, 0, arg.size()-1);
  D_VEC(coeff, 0, coeff.size()-1);
}
