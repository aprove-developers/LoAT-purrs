/* Declaration of the main function of the recurrence relation solver.
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

#ifndef PURRS_rr_solver_hh
#define PURR_rr_solver_hh 1

#include "Expr.types.hh"
#include "Symbol.types.hh"

namespace Parma_Recurrence_Relation_Solver {

enum Solver_Status {
  /*!
    Solution was successful.
  */
  OK,
  /*!
    The right-hand side of the recurrence contains at least an occurrence
    of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is not an integer.
  */
  HAS_NON_INTEGER_DECREMENT,
  /*!
    The right-hand side of the recurrence contains at least an occurrence
    of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is a negative integer.
  */
  HAS_NEGATIVE_DECREMENT,
  /*!
    The right-hand side of the recurrence contains at least an occurrence
    of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is too big to be handled
    by the standard solution techniques.
  */
  HAS_HUGE_DECREMENT,
  /*!
    The right-hand side of the recurrence contains at least an occurrence
    of <CODE>x(n)</CODE>.
  */
  HAS_NULL_DECREMENT,
  /*!
    The recurrence is not linear.
  */
  NON_LINEAR_RECURRENCE,
  /*!
    Catchall: the recurrence is generically too complex.
  */
  TOO_COMPLEX
};

Solver_Status
solve(const Expr& rhs, const Symbol& n, Expr& solution);

Solver_Status
solve_try_hard(const Expr& rhs, const Symbol& n, Expr& solution);

} // namespace Parma_Recurrence_Relation_Solver

#endif
