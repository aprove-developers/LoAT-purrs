/* Recurrence class implementation: inline functions.
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

#ifndef PURRS_Recurrence_inlines_hh
#define PURRS_Recurrence_inlines_hh

#include "rr_solver.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Recurrence::Recurrence()
  : recurrence_rhs(0),
    solved(false) {
}

inline
Recurrence::Recurrence(const Expr& e)
  : recurrence_rhs(e),
    solved(false) {
}

inline
Recurrence::Recurrence(const Recurrence& y)
  : recurrence_rhs(y.recurrence_rhs),
    system_rhs(y.system_rhs),
    solved(y.solved),
    solution(y.solution) {
}

inline
Recurrence::~Recurrence() {
}

inline Recurrence&
Recurrence::operator=(const Recurrence& y) {
  recurrence_rhs = y.recurrence_rhs;
  system_rhs = y.system_rhs;
  solved = y.solved;
  solution = y.solution;
  return *this;
}

inline void
Recurrence::replace_recurrence(const Expr& e) {
  recurrence_rhs = e;
}

inline void
Recurrence::replace_recurrence(unsigned k, const Expr& e) {
  std::pair<std::map<unsigned, Expr>::iterator, bool> stat
    = system_rhs.insert(std::map<unsigned, Expr>::value_type(k, e));
  if (!stat.second)
    // There was already something associated to `k': overwrite it.
    stat.first->second = e;
}

inline bool
Recurrence::solve() const {
  Symbol n("n");
  if (!solved && solve_try_hard(recurrence_rhs, n, solution) == OK)
    solved = true;
  return solved;
}

inline Expr
Recurrence::exact_solution() const {
  if (solved || solve())
    return solution;
  else
    // Well, if the client insists...
    return recurrence_rhs;
}

inline bool
Recurrence::verify_solution() const {
  if (solved || solve()) {
    // Verify the solution.
    return false;
  }
  // Well, if the client insists...
  return true;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Recurrence_inlines_hh)
