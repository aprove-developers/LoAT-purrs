/* Cached_Expr class implementation (non-inline functions).
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

#include "Cached_Expr.defs.hh"

#include "Recurrence.defs.hh"

#include <set>

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {

/*!
  Builds a string joining a single lower-case letter with a number
  (if different from `0'), both created considering the value of `i'.
*/
std::string
new_symbol_name(int i) {
  const char symbol_names[] = "abcdfghijklmopqrstuvwyz";
  const size_t num_names
    = sizeof(symbol_names) / sizeof(symbol_names[0]) - 1;

  // Chooses the letter among those defined in `symbol_names'.
  char letter = symbol_names[i % num_names];
  // Chooses the number that will follow the letter.
  int index = i / num_names;
  std::stringstream s;
  if (index == 0)
    s << letter;
  else
    s << letter << index;
  std::string good_symbol_name;
  s >> good_symbol_name;
  return good_symbol_name;
}

} // anonymous namespace

PURRS::Expr
PURRS::Cached_Expr::remove_bad_symbols(const Recurrence& rec) const {
  Symbol::SymbolSet bad_symbols;
  Symbol::SymbolSet good_symbols;
  Expr new_e = expression();
  // Stores all symbols of `new_e' in two sets, dividing `bad' symbols
  // from `good' symbols.
  new_e.collect_symbols(bad_symbols, good_symbols);

  Symbol::SymbolSet::const_iterator bsi = bad_symbols.begin();
  for (int i = 0, iend = bad_symbols.size(); i < iend; ++i, ++bsi) {
    int try_index = i;
    // Builds the new name `good_symbol_name' for the symbol `*bsi'.
    std::string good_symbol_name = new_symbol_name(try_index);
    // Checks if `good_symbol_name' already exists, in this case
    // builds another name.
    while (good_symbols.find(Symbol(good_symbol_name.c_str()))
	   != good_symbols.end()) {
      ++try_index;
      good_symbol_name = new_symbol_name(try_index);
    }
    Symbol new_symb(good_symbol_name.c_str());
    // The new symbol `new_symb' is added to the set `good_symbols'.
    good_symbols.insert(new_symb);
    // Substitutes all occurrences in `new_e' of the bad symbol with
    // the good symbol built.
    new_e = new_e.substitute(*bsi, new_symb);
    // Substitutes the bad symbol with the good symbol built in
    // the blackboard.
    rec.blackboard.substitute(*bsi, new_symb);
  }
  return new_e;
}
