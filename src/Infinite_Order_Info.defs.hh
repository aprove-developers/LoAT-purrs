/* A class for containing all informations necessary for to solve
   infinite order recurrences.
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

#ifndef PURRS_Infinite_Order_Info_defs_hh
#define PURRS_Infinite_Order_Info_defs_hh 1

#include "Infinite_Order_Info.types.hh"
#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Infinite_Order_Info {
public:
  //! \brief
  //! Constructor: sets
  //! \f$ rhs_transformed_in_first_order_var_coeffs_ = new_rhs \f$,
  //! \f$ weight_ = weight \f$, \f$ infinite_order_fwdr_ = i_c \f$.
  Infinite_Order_Info(const Expr& new_rhs, const Expr& weight,
		      unsigned i_c);

  //! Copy-constructor.
  Infinite_Order_Info(const Infinite_Order_Info& y);

  //! Destructor.
  ~Infinite_Order_Info();

  //! Assignment operator.
  Infinite_Order_Info& operator=(const Infinite_Order_Info& y);

  //! Returns <CODE>rhs_transformed_in_first_order_var_coeffs_</CODE>.
  Expr rhs_transformed_in_first_order_var_coeffs() const;

  //! Returns <CODE>rhs_transformed_in_first_order_var_coeffs_</CODE>.
  Expr& rhs_transformed_in_first_order_var_coeffs();

  //! Returns <CODE>weight_</CODE>.
  Expr weight() const;

  //! Returns <CODE>weight_</CODE>.
  Expr& weight();

  //! Returns <CODE>infinite_order_fwdr_</CODE>.
  unsigned infinite_order_fwdr() const;

  //! Returns <CODE>infinite_order_fwdr_</CODE>.
  unsigned& infinite_order_fwdr();

  //! Sets <CODE>infinite_order_fwdr_</CODE> with \p i_c
  void set_infinite_order_fwdr(unsigned i_c);

private:
  //! \brief
  //! Contains the right hand side of the recurrence obtained
  //! transforming a infinite order recurrence of the form
  //! \f[
  //!   T(n) = f(n) \sum_{k=0}^{n-1} T(k) + g(n).
  //! \f]
  Expr rhs_transformed_in_first_order_var_coeffs_;

  //! \brief
  //! Contains the common coefficient \f$ f(n) \f$ to all terms \f$ x(i) \f$,
  //! for \f$ i = 0, dots, n-1 \f$, of the infinite order recurrence
  //! of the form
  //! \f[
  //!   T(n) = f(n) \sum_{k=0}^{n-1} T(k) + g(n).
  //! \f]
  Expr weight_;

  //! \brief
  //! Stores the smallest positive integer for which the recurrence is
  //! well-defined: the initial conditions will start from it.
  unsigned infinite_order_fwdr_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Infinite_Order_Info.inlines.hh"

#endif // !defined(PURRS_Infinite_Order_Info_defs_hh)