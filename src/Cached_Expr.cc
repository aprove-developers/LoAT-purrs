/* Cached_Expr class implementation (non-inline functions).
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
  std::string new_symbol_name;
  s >> new_symbol_name;
  return new_symbol_name;
}

} // anonymous namespace

PURRS::Expr
PURRS::Cached_Expr::
replace_system_generated_symbols(const Recurrence& rec) const {
  Symbol::SymbolSet system_generated_symbols;
  Symbol::SymbolSet new_symbols;
  Expr new_e = expression();
  // Stores all symbols of `new_e' in two sets, dividing `bad' symbols
  // from `good' symbols.
  new_e.collect_symbols(system_generated_symbols, new_symbols);

  Symbol::SymbolSet::const_iterator bsi = system_generated_symbols.begin();
  for (int i = 0, iend = system_generated_symbols.size();
       i < iend; ++i, ++bsi) {
    int try_index = i;
    // Builds the new name `new_name' for the symbol `*bsi'.
    std::string new_name = new_symbol_name(try_index);
    // Checks if `new_name' already exists, in this case builds another name.
    while (new_symbols.find(Symbol(new_name.c_str())) != new_symbols.end()) {
      ++try_index;
      new_name = new_symbol_name(try_index);
    }
    Symbol new_symb(new_name.c_str());
    // The new symbol `new_symb' is added to the set `new_symbols'.
    new_symbols.insert(new_symb);
    // Substitutes all occurrences in `new_e' of the symbol generated
    // by the system with the new symbol built.
    new_e = new_e.substitute(*bsi, new_symb);
    // Substitutes the symbol generated by the system with the good
    // symbol built in the blackboard.
    rec.blackboard.substitute(*bsi, new_symb);
  }
  return new_e;
}
