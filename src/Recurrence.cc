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
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

const PURRS::Symbol&
PURRS::Recurrence::n = Symbol("n");

PURRS::Expr
PURRS::Recurrence::substitute_auxiliary_definition(const Expr& e) const {
  Expr e_after_subs;
  if (e.is_a_add()) {
    e_after_subs = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_after_subs += substitute_auxiliary_definition(e.op(i));
  }
  else if (e.is_a_mul()) {
    e_after_subs = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_after_subs *= substitute_auxiliary_definition(e.op(i));
  }
  else if (e.is_a_power())
    e_after_subs = pwr(substitute_auxiliary_definition(e.op(0)),
		       substitute_auxiliary_definition(e.op(1)));
  else if (e.is_a_function()) {
    // FIXME: evitare la copia e trovare come accedere al funtore.
    // e_after_subs = functor(e)(simplify_on_input_ex(e.op(0), n, input));
    Expr tmp = substitute_auxiliary_definition(e.op(0));
    Expr f = e;
    e_after_subs = f.subs(f.op(0), tmp);
  }
  else if (e.is_a_symbol()) {
    Symbol s = e.ex_to_symbol();
    e_after_subs = get_auxiliary_definition(s);
  }
  else
    e_after_subs = e;
  return e_after_subs;
}

void
PURRS::Recurrence::dump(std::ostream& s) const {
  s << "solved = " << (solved ? "true" : "false") << std::endl;
  s << "recurrence_rhs = " << recurrence_rhs << std::endl;
  s << "auxiliary_definitions = ";
  if (auxiliary_definitions.empty())
    s << "empty" << std::endl;
  else {
    s << std::endl;
    typedef std::map<Symbol, Expr> Map;
    for (Map::const_iterator i = auxiliary_definitions.begin(),
	   ad_end = auxiliary_definitions.end();
	 i != ad_end;
	 ++i)
      s << "  " << i->first << " = " << i->second << std::endl;
  }
  //s << "solution = " << solution << std::endl;
}

