/* A class for containing the necessary informations about finite
   order recurrences (and that we do not want to compute them
   again).
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
#include "Expr.defs.hh"
#include <vector>

namespace Parma_Recurrence_Relation_Solver {

class Finite_Order_Info {
public:
  //! \brief
  //! Constructor: sets \f$ order_ = k \f$,
  //! \f$ first_initial_condition_ = i_c \f$, and
  //! \f$ coefficients_ = coeffs \f$.
  Finite_Order_Info(int k, unsigned i_c, const std::vector<Expr>& coeffs);

  //! Copy-constructor.
  Finite_Order_Info(const Finite_Order_Info& y);

  //! Destructor.
  ~Finite_Order_Info();

  //! Assignment operator.
  Finite_Order_Info& operator=(const Finite_Order_Info& y);

  //! Returns <CODE>order_</CODE>.
  unsigned int order() const;

  //! Returns <CODE>order_</CODE>.
  unsigned int& order();

  //! Returns <CODE>first_initial_condition_</CODE>.
  unsigned first_initial_condition() const;

  //! Returns <CODE>first_initial_condition_</CODE>.
  unsigned& first_initial_condition();

  //! Sets <CODE>first_initial_condition_</CODE> with \p i_c
  void set_first_initial_condition(unsigned i_c);

  //! Returns <CODE>coefficients_</CODE>.
  const std::vector<Expr>& coefficients() const;

  //! Returns <CODE>coefficients_</CODE>.
  std::vector<Expr>& coefficients();

private:
  //! The order of the recurrence. 
  unsigned int order_;

  //! \brief
  //! The smallest positive integer for which the recurrence is
  //! well-defined: the initial conditions will start from it.
  unsigned first_initial_condition_;

  //! Stores the coefficients of the recurrence.
  std::vector<Expr> coefficients_;

};

} // namespace Parma_Recurrence_Relation_Solver

#include "Finite_Order_Info.inlines.hh"

#endif // !defined(PURRS_Finite_Order_Info_defs_hh)
