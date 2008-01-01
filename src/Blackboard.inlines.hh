/* Blackboard class implementation: inline functions.
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

#ifndef PURRS_Blackboard_inlines_hh
#define PURRS_Blackboard_inlines_hh

#include "util.hh"

#include <iostream>
#include <utility>
#include <deque>
#include <map>

namespace Parma_Recurrence_Relation_Solver {

template <typename T>
inline
Blackboard::Cached<T>::Cached()
  // Please note: the following line reads "timestamp-of-zero".
  // Since the origin of time for the blackboard is 1, this means
  // that cached elements are constructed out-of-date.
  : timestamp(0) {
}

inline
Blackboard::Definition::Definition(const Expr& e)
  : rhs(e) {
}

inline
Blackboard::Blackboard()
  // Please note: the following line reads "timestamp-of-one".
  : timestamp(1) {
}

inline
Blackboard::Blackboard(const Blackboard& y)
  : index(y.index), definitions(y.definitions), timestamp(y.timestamp) {
}

inline
Blackboard::~Blackboard() {
}

inline Blackboard&
Blackboard::operator=(const Blackboard& y) {
  index = y.index;
  definitions = y.definitions;
  timestamp = y.timestamp;
  return *this;
}

inline Symbol
Blackboard::insert_definition(const Expr& e) {
  Symbol new_symbol;
  index.insert(std::map<Symbol, unsigned int>::value_type(new_symbol,
							  definitions.size()));
  definitions.push_back(Definition(e));
  ++timestamp;
  return new_symbol;
}

inline Expr
Blackboard::get_definition(const Symbol& z) const {
  std::map<Symbol, unsigned int>::const_iterator i = index.find(z);
  if (i != index.end())
    return definitions[i->second].rhs;
  else
    return z;
}

inline Symbol
Blackboard::find_symbol(const Expr& e) {
  std::map<Symbol, unsigned int>::const_iterator j = index.begin();
  for (std::deque<Definition>::const_iterator i = definitions.begin(),
	 definitions_end = definitions.end(); i != definitions_end; ++i, ++j) {
    if ((*i).rhs == e)
      return j->first;
  }
  return insert_definition(e);
}

inline void
Blackboard::substitute(const Symbol& system_generated_symbol,
		       const Symbol& new_symbol) {
  std::map<Symbol, unsigned int>::const_iterator i
    = index.find(system_generated_symbol);
  // Could exist a `system_generated_symbol' symbol that are
  // not stored in the blackboard (e.g. the index of the sum).
  if (i != index.end()) {
    // Insert the equation `new_symbol = e', where `e' is the expression
    // that was associated to `system_generated_symbol'.
    index.insert(std::map<Symbol, unsigned int>::value_type(new_symbol,
							    definitions.size()));
    definitions.push_back(Definition(definitions[i->second].rhs));
    ++timestamp;
    // Remove from the blackboard the `system_generated_symbol' element.
    index.erase(system_generated_symbol);
  }
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Blackboard_inlines_hh)
