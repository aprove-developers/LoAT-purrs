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
#include "globals.hh"
#include "Expr.defs.hh"
#include "Recurrence.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Weighted_Average_Info {
public:
  //! \brief
  //! Constructor: sets
  //! \f$ associated_first_order_rec_ = associated_first_order_rec \f$,
  //! \f$ infinite_order_weight_ = infinite_order_weight \f$.
  Weighted_Average_Info(const Recurrence& associated_first_order_rec,
			const Expr& infinite_order_weight);

  //! Copy-constructor.
  Weighted_Average_Info(const Weighted_Average_Info& y);

  //! Destructor.
  ~Weighted_Average_Info();

  //! Assignment operator.
  Weighted_Average_Info& operator=(const Weighted_Average_Info& y);

  //! Returns <CODE>associated_first_order_rec_</CODE>.
  const Recurrence& associated_first_order_rec() const;

  //! Returns <CODE>associated_first_order_rec_</CODE>.
  Recurrence& associated_first_order_rec();

  //! Sets <CODE>original_rhs__</CODE> with \p original_rhs
  void set_original_rhs(const Expr& original_rhs);

  //! Returns <CODE>infinite_order_weight_</CODE>.
  const Expr& infinite_order_weight() const;

  //! Returns <CODE>infinite_order_weight_</CODE>.
  Expr& infinite_order_weight();

  //! Returns <CODE>first_valid_index_inf_order_</CODE>.
  index_type first_valid_index_inf_order() const;

  //! Sets <CODE>first_valid_index_inf_order_</CODE> with \p i_c
  void set_first_valid_index_inf_order(index_type i_c);

private:
  //! \brief
  //! In the case which the system is able to rewrite the infinite order
  //! recurrence \p *this in a first order recurrence, this method stores
  //! the first order recurrence computed (in order to know the cases of
  //! rewritable infinite order recurrences see the function
  //! <CODE>rewrite_infinite_order_recurrence()</CODE>).
  Recurrence associated_first_order_rec_;

  //! \brief
  //! When the recurrence is not in normal form, this data contains
  //! its right hand side before the transformation in normal form.
  //! If the recurrence is already in normal form it is undefined.
  Expr original_rhs_;

  //! \brief
  //! Contains the factor \f$ f(n) \f$ of the infinite order recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n).
  //! \f]
  Expr infinite_order_weight_;

  //! \brief
  //! Stores the least non-negative integer \f$ j \f$ such that
  //! the recurrence is well-defined for \f$ n \geq j \f$.
  index_type first_valid_index_inf_order_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Infinite_Order_Info.inlines.hh"

#endif // !defined(PURRS_Infinite_Order_Info_defs_hh)
