/* Symbol class implementation (non-inline functions).
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

#include "Symbol.defs.hh"

#include <string>

namespace PURRS = Parma_Recurrence_Relation_Solver;

bool
PURRS::Symbol::is_a_bad_symbol() const {
  const std::string str = get_name();
  // `str' has less than 7 characters: it can not be a bad symbol.
  if (str.size() < 7)
    return false;
  else {
    std::string::size_type pos = str.find("symbol");
    // `str' does not have the substring `symbol' or has the substring
    // `symbol' but not in the first positions of the string:
    // it can not be a bad symbol. 
    if (pos != 0)
      return false;
    const std::string numbers("0123456789");
    // Finds the first non-numeric character in `str' occurring after `symbol'.
    pos = str.find_first_not_of(numbers, 6);
    return (pos == std::string::npos);
  }
}
