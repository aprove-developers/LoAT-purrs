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
  //! Constructor: sets \f$ order{\_} = k \f$,
  //! \f$ first{\_}i{\_}c{\_}for{\_}linear{\_} = i{\_}c \f$,
  //! \f$ coefficients{\_}lfo{\_} = coeffs \f$ and
  //! \f$ gcd{\_}among{\_}decrements{\_} = gcd \f$.
  Finite_Order_Info(int k, unsigned i_c, const std::vector<Expr>& coeffs,
		    unsigned gcd);

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

  //! Returns <CODE>first_i_c_for_linear_</CODE>.
  unsigned first_i_c_for_linear() const;

  //! Returns <CODE>first_i_c_for_linear_</CODE>.
  unsigned& first_i_c_for_linear();

  //! Sets <CODE>first_i_c_for_linear_</CODE> with \p i_c
  void set_first_i_c_for_linear(unsigned i_c);

  //! Returns <CODE>coefficients_</CODE>.
  const std::vector<Expr>& coefficients() const;

  //! Returns <CODE>coefficients_</CODE>.
  std::vector<Expr>& coefficients();

  //! Returns <CODE>gcd_among_decrements_</CODE>.
  unsigned gcd_among_decrements() const;

  //! Returns <CODE>gcd_among_decrements_</CODE>.
  unsigned& gcd_among_decrements();

private:
  //! The order of the recurrence. 
  unsigned int order_;

  //! \brief
  //! The smallest positive integer for which the recurrence is
  //! well-defined: the initial conditions will start from it.
  unsigned first_i_c_for_linear_;

  //! Stores the coefficients of the recurrence.
  std::vector<Expr> coefficients_;

  //! \brief
  //! Stores the greatest common divisor among the positive integer \f$ k \f$
  //! of the terms of the form \f$ x(n-k) \f$ contained in the right
  //! hand side of the recurrence. If the order is zero then
  //! \p gcd_among_decrements_ stores \f$ 0 \f$.
  unsigned gcd_among_decrements_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Finite_Order_Info.inlines.hh"

#endif // !defined(PURRS_Finite_Order_Info_defs_hh)
