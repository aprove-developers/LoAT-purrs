/* Blackboard class implementation (non-inline functions).
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

#include "Blackboard.defs.hh"
#include "util.hh"
#include "approximate_impl.hh"
#include <iostream>

namespace PURRS = Parma_Recurrence_Relation_Solver;

PURRS::Expr
PURRS::Blackboard::rewrite(Cached<Expr>& ce) const {
  if (timestamp > ce.timestamp) {
    ce.value = rewrite(ce.value);
    ce.timestamp = timestamp;
  }
  return ce.value;
}

PURRS::Expr
PURRS::Blackboard::rewrite(const Expr& e) const {
  Expr e_rewritten;
  if (e.is_a_add()) {
    e_rewritten = 0;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten += rewrite(e.op(i));
  }
  else if (e.is_a_mul()) {
    e_rewritten = 1;
    for (unsigned i = e.nops(); i-- > 0; )
      e_rewritten *= rewrite(e.op(i));
  }
  else if (e.is_a_power())
    e_rewritten = pwr(rewrite(e.arg(0)),
		       rewrite(e.arg(1)));
  else if (e.is_a_function()) {
    if (e.nops() == 1)
      e_rewritten = apply(e.functor(),
			   rewrite(e.arg(0)));
    else {
      unsigned num_argument = e.nops();
      std::vector<Expr> argument(num_argument);
      for (unsigned i = 0; i < num_argument; ++i)
	argument[i] = rewrite(e.arg(i));
      e_rewritten = apply(e.functor(), argument);
    }
  }
  else if (e.is_a_symbol()) {
    Symbol z = e.ex_to_symbol();
    std::map<Symbol, unsigned>::const_iterator i = index.find(z);
    if (i != index.end())
      e_rewritten = rewrite(definitions[i->second].expansion);
    else
      e_rewritten = e;
  }
  else
    e_rewritten = e;
  return e_rewritten;
}

unsigned
PURRS::Blackboard::size_norm(Cached<unsigned>& ce) const {
  if (timestamp > ce.timestamp) {
    ce.value = generic_size_norm(ce.value, *this);
    ce.timestamp = timestamp;
  }
  return ce.value;
}

unsigned
PURRS::Blackboard::size_norm(const Expr& e) const {
  return generic_size_norm(e, *this);
}

void
PURRS::Blackboard::dump(std::ostream& s) const {
  if (definitions.empty())
    s << "Blackboard empty." << std::endl;
  else {
    s << "Blackboard contents:" << std::endl;
    for (std::map<Symbol, unsigned>::const_iterator i = index.begin(),
	   index_end = index.end(); i != index_end; ++i)
      s << "  " << i->first
	<< " = " << definitions[i->second].rhs << std::endl;
  }
  s << std::endl;
}