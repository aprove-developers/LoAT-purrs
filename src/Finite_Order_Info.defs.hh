/* To be written.
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

#ifndef PURRS_Finite_Order_Info_defs_hh
#define PURRS_Finite_Order_Info_defs_hh 1

#include "Finite_Order_Info.types.hh"
#include <vector>

namespace Parma_Recurrence_Relation_Solver {

class Finite_Order_Info {
public:
  //! Constructor: set \f$ order = k \f$.
  Finite_Order_Info(int k);

  //! Copy-constructor.
  Finite_Order_Info(const Finite_Order_Info& y);

  //! Destructor.
  ~Finite_Order_Info();

  //! Assignment operator.
  Finite_Order_Info& operator=(const Finite_Order_Info& y);

  // Returns the order of the recurrence.
  int get_order();

  //! Sets the vector <CODE>decrements</CODE> from \p dec.
  void set_decrements(const std::vector<unsigned> dec);

  //! Sets the vector <CODE>initial_conditions</CODE> from \p i_c.
  void set_initial_conditions(const std::vector<unsigned> i_c);

  //! Adds \p d to the vector <CODE>decrements</CODE>.
  void add_decrements(unsigned d);

  //! Adds \p i to the vector <CODE>initial_conditions</CODE>.
  void add_initial_conditions(unsigned i);

private:
  friend class Recurrence;

  //! The order of the recurrence. 
  int order;

  //! Stores the positive integers \f$ d \f$ of all the \f$ x(n-d) \f$
  //! contained in the right hand side of the recurrence.
  std::vector<unsigned> decrements;

  //! Stores the positive integers that represent the initial conditions
  //! for the recurrences. For these integers the recurrence is well defined.
  std::vector<unsigned> initial_conditions;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Finite_Order_Info.inlines.hh"

#endif // !defined(PURRS_Finite_Order_Info_defs_hh)
