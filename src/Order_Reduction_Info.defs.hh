/* A class for containing all informations necessary in the
   order reduction's method.
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

#ifndef PURRS_Order_Reduction_Info_defs_hh
#define PURRS_Order_Reduction_Info_defs_hh 1

#include "Order_Reduction_Info.types.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Order_Reduction_Info {
public:
  //! \brief
  //! Constructor: sets \f$ old_recurrence_rhs_ = rhs \f$,
  //! \f$ gcd_decrements_old_rhs_ = gcd \f$, and
  //! \f$ symbol_for_mod_ = r \f$.
  Order_Reduction_Info(const Expr& rhs, unsigned gcd, const Symbol& r);

  //! Copy-constructor.
  Order_Reduction_Info(const Order_Reduction_Info& y);

  //! Destructor.
  ~Order_Reduction_Info();

  //! Assignment operator.
  Order_Reduction_Info& operator=(const Order_Reduction_Info& y);

  //! Returns <CODE>old_recurrence_rhs_</CODE>.
  Expr old_recurrence_rhs() const;

  //! Returns <CODE>old_recurrence_rhs_</CODE>.
  Expr& old_recurrence_rhs();

  //! Returns <CODE>gcd_decrements_old_rhs_</CODE>.
  unsigned gcd_decrements_old_rhs() const;

  //! Returns <CODE>gcd_decrements_old_rhs_</CODE>.
  unsigned& gcd_decrements_old_rhs();

  //! Sets <CODE>gcd_decrements_old_rhs_</CODE> with \p g.
  void set_gcd_decrements_old_rhs(unsigned g);

  //! Returns <CODE>symbol_for_mod_</CODE>.
  Symbol symbol_for_mod() const;

  //! Returns <CODE>symbol_for_mod_</CODE>.
  Symbol& symbol_for_mod();

  //! Returns <CODE>solution_order_reduced_</CODE>.
  Expr solution_order_reduced() const;

  //! Returns <CODE>solution_order_reduced_</CODE>.
  Expr& solution_order_reduced();

  //! Sets <CODE>solution_order_reduced_</CODE> with \p e.
  void set_solution_order_reduced(const Expr& e);

  //! Returns <CODE>verified_one_time_</CODE>.
  bool verified_one_time() const;

  //! Returns <CODE>verified_one_time_</CODE>.
  void not_verified_one_time(); 

private:
  //! \brief
  //! Stores the value of the <CODE>recurrence_rhs</CODE> before
  //! to apply the order reduction.
  Expr old_recurrence_rhs_;

  //! \brief
  //! Stores the greatest common divisor among the decrements
  //! <CODE>d</CODE> of the terms <CODE>x(n-d)</CODE> present in the
  //! right hand side of the recurrence before to apply the order reduction.
  unsigned gcd_decrements_old_rhs_;

  //! \brief
  //! Stores the auxiliary symbol used in place of the function
  //! \f$ mod() \f$.
  Symbol symbol_for_mod_;

  //! Stores the solution of the reduced order recurrence.
  Expr solution_order_reduced_;

  //! \brief
  //! If \p verified_one_time_ is true means that the system will
  //! try to verify the non-reduced recurrence; if \p verified_one_time_
  //! is false then the system will verify only the reduced recurrence.
  bool verified_one_time_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Order_Reduction_Info.inlines.hh"

#endif // !defined(PURRS_Order_Reduction_Info_defs_hh)
