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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#include <config.h>

#include "factorize_giac.hh"
#include "util.hh"
#include "Blackboard.defs.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Number.defs.hh"
#include <giac/sym2poly.h>
#include <giac/usual.h>
#include <giac/symbolic.h>
#include <utility>
#include <string>

using namespace giac;
namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {
using namespace PURRS;
using namespace giac;
using namespace std;

void
factorize_giac_recursive(const Expr& e,
			 Expr& num_factors, Expr& not_num_factors);

void
print_indent(unsigned int i) {
  for ( ; i-- > 0; )
    cout << " ";
}

void
myprint(const gen& e_giac, unsigned int indent = 0) {
  //  cout << "Subtype = " << e_giac.subtype << endl;
  print_indent(indent);
  if (e_giac.is_integer())
    cout << "integer: " << e_giac << endl;
  else if (e_giac.type == _FRAC) {
    cout << "frac: " << endl;
    print_indent(indent);
    cout << "numerator: ";
    myprint(e_giac._FRACptr->num, indent+2);
    print_indent(indent);
    cout << "denominator: ";
    myprint(e_giac._FRACptr->den, indent+2);
  }
  else if (e_giac.type == _IDNT)
    cout << "identif: " << *(e_giac._IDNTptr->name) << endl;
  else if (e_giac.type == _SYMB) {
    cout << "symbolic: ";
    cout << e_giac._SYMBptr->sommet.ptr->s << endl;
    myprint(e_giac._SYMBptr->feuille, indent+2);
  }
  else if (e_giac.type == _VECT) {
    vecteur &v = *e_giac._VECTptr;
    for (unsigned int i = v.size(); i-- > 0; ) {
      myprint(v[i], indent+2);
    }
  }
}

unary_function_ptr
find_functor_ptr(const std::string& str, unary_function_ptr tab[]) {
  for (int i = 1; *tab != 0; ++tab, ++i) {
    unary_function_ptr tmp = *tab;
    if (tmp.ptr->s == str)
      return tmp;
  }
  return 0;
}

/*!
  Translate the GiNaC expression stored in \p e in a giac expression.
*/
gen
translate(const Expr& e, Expr& num_factors, Expr& not_num_factors,
	  Blackboard& blackboard) {
  gen e_giac;
  if (e.is_a_add()) {
    e_giac = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_giac = e_giac + translate(e.op(i), num_factors, not_num_factors,
				  blackboard);
  }
  if (e.is_a_mul()) {
    e_giac = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_giac = e_giac * translate(e.op(i), num_factors, not_num_factors,
				  blackboard);
  }
  else if (e.is_a_power())
    return pow(translate(e.arg(0), num_factors, not_num_factors, blackboard),
	       translate(e.arg(1), num_factors, not_num_factors, blackboard));
  else if (e.is_a_function()) {
    if (e.nops() == 1) {
      // Recursevely factorize the argument of the function `e'.
      Expr num_factors_arg;
      Expr not_num_factors_arg;
      factorize_giac_recursive(e.arg(0), num_factors_arg, not_num_factors_arg);
      unary_function_ptr functor = find_functor_ptr(e.get_function_name(),
						    archive_function_tab);
      // FIXME: Not to use the constructor of the class `symbolic' but 
      // to use the identificateur (as for the functions with more than one
      // arguments).
      // Using the identificateur is better than the symbolic object
      // because in this way is applicable the factorization's process
      // (on symbolic object the factorization's process does not work).
      if (functor == 0)
#if 1
	return symbolic(unary_function_abstract(e.get_function_name()),
			translate(num_factors_arg*not_num_factors_arg,
				  num_factors, not_num_factors,
				  blackboard));
#else
      {
// 	const Expr& tmp = apply(e.functor(), argument);
// 	return identificateur(blackboard.find_symbol(tmp).get_name());
	
      }
#endif
      else
#if 1
	return symbolic(functor, translate(num_factors_arg*not_num_factors_arg,
					   num_factors, not_num_factors,
					   blackboard));
#else
      {
	
      }
#endif
    }
    else {
      // Recursevely factorize the arguments of the function `e'.
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i) {
	Expr num_factors_arg;
	Expr not_num_factors_arg;
	factorize_giac_recursive(e.arg(i),
				 num_factors_arg, not_num_factors_arg);
	argument[i] = num_factors_arg * not_num_factors_arg;
      }
      const Expr& tmp = apply(e.functor(), argument);
      return identificateur(blackboard.find_symbol(tmp).get_name());
    }
  }
  else if (e.is_a_symbol())
    return identificateur(blackboard.find_symbol(e).get_name());
  else if (e.is_a_number()) {
    Number e_num = e.ex_to_number();
    if (e_num.is_integer())
      return e_num.to_int();
    else if (e_num.is_rational())
      return fraction(e_num.numerator().to_int(),
		      e_num.denominator().to_int());
    else {
      assert(e_num.is_complex_rational());
      return identificateur(blackboard.find_symbol(e).get_name());
    }
  }
  else if (e.is_a_constant())
    return identificateur(blackboard.find_symbol(e).get_name());
  
  return e_giac;
}

#if 0
/*!
  For `factors'. 
*/
gen
translate(const Expr& e, Blackboard& blackboard) {
  gen e_giac;
  if (e.is_a_add()) {
    e_giac = 0;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_giac = e_giac + translate(e.op(i), blackboard);
  }
  if (e.is_a_mul()) {
    e_giac = 1;
    for (unsigned int i = e.nops(); i-- > 0; )
      e_giac = e_giac * translate(e.op(i), blackboard);
  }
  else if (e.is_a_power())
    return pow(translate(e.arg(0), blackboard),
	       translate(e.arg(1), blackboard));
  else if (e.is_a_function()) {
    if (e.nops() == 1) {
      // Recursevely factorize the argument of the function `e'.
      Expr num_factors_arg;
      Expr not_num_factors_arg;
      factorize_giac_recursive(e.arg(0), num_factors_arg, not_num_factors_arg);
      unary_function_ptr functor = find_functor_ptr(e.get_function_name(),
						    archive_function_tab);
      // FIXME: Not to use the constructor of the class `symbolic' but 
      // to use the identificateur (as for the functions with more than one
      // arguments).
      // Using the identificateur is better than the symbolic object
      // because in this way is applicable the factorization's process
      // (on symbolic object the factorization's process does not work).
      if (functor == 0)
	return symbolic(unary_function_abstract(e.get_function_name()),
			translate(num_factors_arg*not_num_factors_arg,
				  blackboard));
      else
	return symbolic(functor, translate(num_factors_arg*not_num_factors_arg,
					   blackboard));
    }
    else {
      // Recursevely factorize the arguments of the function `e'.
      unsigned int num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned int i = 0; i < num_argument; ++i) {
	Expr num_factors_arg;
	Expr not_num_factors_arg;
	factorize_giac_recursive(e.arg(i),
				 num_factors_arg, not_num_factors_arg);
	argument[i] = num_factors_arg * not_num_factors_arg;
      }
      const Expr& tmp = apply(e.functor(), argument);
      return identificateur(blackboard.find_symbol(tmp).get_name());
    }
  }
  else if (e.is_a_symbol())
    return identificateur(blackboard.find_symbol(e).get_name());
  else if (e.is_a_number()) {
    Number e_num = e.ex_to_number();
    if (e_num.is_integer())
      return e_num.to_int();
    else if (e_num.is_rational())
      return fraction(e_num.numerator().to_int(),
		      e_num.denominator().to_int());
    else {
      assert(e_num.is_complex_rational());
      return identificateur(blackboard.find_symbol(e).get_name());
    }
  }
  else if (e.is_a_constant())
    return identificateur(blackboard.find_symbol(e).get_name());
  
  return e_giac;
}
#endif

/*!
  Translate the giac expression stored in \p e_giac in a GiNaC
  expression.
*/
Expr
translate(const gen& e_giac) {
  Expr e;
  if (e_giac.is_integer())
    return e_giac.val;
  // At the moment the only purrs'expressions transformed in
  // giac's fraction are rational numbers.
  else if (e_giac.type == _FRAC)
    return Number(e_giac._FRACptr->num.to_int(),
		  e_giac._FRACptr->den.to_int());
  else if (e_giac.type == _IDNT)
    return Symbol(*(e_giac._IDNTptr->name));
  else if (e_giac.type == _SYMB) {
    if (e_giac.is_symb_of_sommet(at_prod)) {
      gen& arg = e_giac._SYMBptr->feuille;
      e = 1;
      if (arg.type == _VECT) {
	vecteur &v = *arg._VECTptr;
	for (unsigned int i = v.size(); i-- > 0; )
	  e *= translate(v[i]);
      }
    }
    else if (e_giac.is_symb_of_sommet(at_plus)) {
      gen& arg = e_giac._SYMBptr->feuille;
      if (arg.type == _VECT) {
	vecteur &v = *arg._VECTptr;
	for (unsigned int i = v.size(); i-- > 0; )
	  e += translate(v[i]);
      }
    }
    else if (e_giac.is_symb_of_sommet(at_inv))
      return pwr(translate(e_giac._SYMBptr->feuille), -1);
    else if (e_giac.is_symb_of_sommet(at_pow))
      return pwr(translate(e_giac._SYMBptr->feuille._VECTptr->front()),
		 translate(e_giac._SYMBptr->feuille._VECTptr->back()));
    else {
      Functor f = find_functor(e_giac._SYMBptr->sommet.ptr->s);
      Expr arg = translate(e_giac._SYMBptr->feuille);
      return apply(f, arg); 
    }
  }

  return e;
}

void
factorize_giac_recursive(const Expr& e,
			 Expr& num_factors, Expr& not_num_factors) {
  Blackboard blackboard;
  // Translate GiNaC expressions in giac expressions.
  gen giac_e = translate(e, num_factors, not_num_factors, blackboard);
  // Call the giac function in order to factorize expression.
  gen giac_e_factorized = factor(giac_e);
  // The decomposition executed by `giac' considers the numeric factors,
  // but when we come back to the GiNaC's expressions, the numeric
  // factors are automatically distributed. In order to avoid this,
  // we divide the numeric factors from the non-numeric factors.
  gen numeric_factors = 1;
  gen not_numeric_factors = 1;
  if (giac_e_factorized.is_symb_of_sommet(at_prod)) {
    gen& arg = giac_e_factorized._SYMBptr->feuille;
    if (arg.type == _VECT) {
      vecteur &v = *arg._VECTptr;
      for (unsigned int i = v.size(); i-- > 0; )
	if (v[i].is_integer() || v[i].type == _FRAC
	    || v[i].is_symb_of_sommet(at_inv))
	  numeric_factors = numeric_factors * v[i];
	else
	  not_numeric_factors = not_numeric_factors * v[i];
    }
  }
  else
    if (giac_e_factorized.is_integer() || giac_e_factorized.type == _FRAC
	|| giac_e_factorized.is_symb_of_sommet(at_inv))
      numeric_factors = numeric_factors * giac_e_factorized;
    else
      not_numeric_factors = not_numeric_factors * giac_e_factorized;
  num_factors = translate(numeric_factors);
  not_num_factors = translate(not_numeric_factors);
  num_factors = blackboard.rewrite(num_factors);
  not_num_factors = blackboard.rewrite(not_num_factors);
}

#if 0
void
factors_giac_recursive(const Expr& e,
		       std::vector<std::pair<Expr, int> >& factors_e) {
  Blackboard blackboard;
  gen giac_e = translate(e, blackboard);
  gen giac_e_factorized = factors(giac_e);
  translate(giac_e_factorized, factors_e);
  for (unsigned int i = factors_e.size(); i-- > 0; )
    blackboard.rewrite(factors_e[i].first);
}
#endif

} // anonymous namespace


/*!
  Factorizes as much as possible the expression \p e dividing
  numeric factors and not numeric factors, respectively stored in
  \p num_factors and \p not_num_factors.
*/
void
PURRS::factorize_giac(const Expr& e,
		      Expr& num_factors, Expr& not_num_factors) {
  factorize_giac_recursive(e, num_factors, not_num_factors);
  // FIXME: to uncomment the following rows when will be solved
  // the problem about the method `expand()' on functions.
  // assert(e.expand() == (num_factors * not_num_factors).expand());
}

#if 0
void
PURRS::factors_giac(const Expr& e,
		    std::vector<std::pair<Expr, int> >& factors_e) {
  factors_giac_recursive(e, factors_e);
}
#endif
