/* This is the main object PURRS operates upon: a recurrence.
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

#ifndef PURRS_Recurrence_defs_hh
#define PURRS_Recurrence_defs_hh 1

#include "Recurrence.types.hh"
#include "Expr.defs.hh"
#include <map>

namespace Parma_Recurrence_Relation_Solver {

/*!
  Explain...
*/
class Recurrence {
public:
  //! Default constructor: it builds the recurrence \f$ x(n) = 0 \f$.
  Recurrence();

  //! Builds the recurrence \f$ x(n) = e \f$.
  explicit Recurrence(const Expr& e);

  //! Copy-constructor.
  Recurrence(const Recurrence& y);

  //! Destructor.
  ~Recurrence();

  //! Assignment operator.
  Recurrence& operator=(const Recurrence& y);

  //! WRITEME
  void replace_recurrence(const Expr& e);

  //! \brief
  //! Sets to \f$ e \f$ the right-hand side of the recurrence
  //! of index \f$ k \f$.  The system of recurrences will then
  //! include \f$ x_k(n) = e \f$.
  void replace_recurrence(unsigned k, const Expr& e);

  //! 
  Symbol insert_auxiliary_equation(const Expr& e);

  bool solve() const;
  Expr exact_solution() const;
  bool verify_solution() const;

private:
  //! Holds the right-hand side of the global recurrence to be solved.
  //! This may have been set directly by the constructor or it may be the
  //! result of transforming a system into a single recurrence.
  //! The global recurrence is thus of the form
  //! <CODE>x(n) = recurrence_rhs</CODE>.
  Expr recurrence_rhs;

  //! Holds the right-hand sides of a system of  recurrence equations.
  //! If <CODE>i == system_rhs.find(k)</CODE> then
  //! <CODE>x(k,n) = (*i).second()</CODE>
  //! is one of the equations of the system.
  std::map<unsigned, Expr> system_rhs;

  mutable bool solved;

  mutable Expr solution;

#if 0
    Expr poly_char;
    Expr symb_solution;

public:
  insert_exact_solution(const Expr& e);
  insert_side_condition(int n, const Expr& e);
  // Restituisce le condizioni iniziali non assegnate.
  std::set<unsigned int> undefined_side_conditions();
  complex_interval approximate(const std::vector<Expr> side_conditions,
			       int n);
//    Expr big_o(...);
//    Expr big_omega(...);
#endif
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Recurrence.inlines.hh"

#endif // !defined(PURRS_Recurrence_defs_hh)

