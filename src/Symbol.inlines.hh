/* Symbol class implementation: inline functions.
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

#ifndef PURRS_Symbol_inlines_hh
#define PURRS_Symbol_inlines_hh

namespace Parma_Recurrence_Relation_Solver {

inline
Symbol::Symbol() {
}

inline
Symbol::Symbol(const char* n)
  : s(n) {
}

inline
Symbol::Symbol(const std::string& str)
  : s(str) {
}

inline
Symbol::Symbol(const Symbol& y)
  : s(y.s) {
};

inline Symbol&
Symbol::operator=(const Symbol& y) {
  s = y.s;
  return *this;
};

inline
Symbol::Symbol(const GiNaC::symbol& gs)
  : s(gs) {
}

inline
Symbol::~Symbol() {
}

inline std::string
Symbol::get_name() const {
  return s.get_name();
}

inline bool
Symbol::NameCompare::operator()(const Symbol& x, const Symbol& y) const {
  return x.get_name() < y.get_name();
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Symbol_inlines_hh)
