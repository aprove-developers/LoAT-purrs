/* Functional_Equation_Info class implementation (non-inline functions).
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

#include "Functional_Equation_Info.defs.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

void
PURRS::Functional_Equation_Info::
dump_homogeneous_terms(std::ostream& s) const {
  if (!homogeneous_terms.empty()) {
    s << "Homogeneous term:" << std::endl;
    for (std::map<Number, Expr>::const_iterator i = homogeneous_terms.begin(),
	   homogeneous_terms_end = homogeneous_terms.end();
	 i != homogeneous_terms_end; ++i)
      s << "div = " << i->first
	<< ", coeff = " << i->second << std::endl;
  }
}
