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

  Expr old_recurrence_rhs() const;
  Expr& old_recurrence_rhs();

  unsigned gcd_decrements_old_rhs() const;
  unsigned& gcd_decrements_old_rhs();
  void set_gcd_decrements_old_rhs(unsigned g);

  Symbol symbol_for_mod() const;
  Symbol& symbol_for_mod();

  Expr solution_order_reduced() const;
  Expr& solution_order_reduced();
  void set_solution_order_reduced(const Expr& e);

  bool verified_one_time() const;
  void not_verified_one_time(); 

private:
  //! \brief
  //! When is applied the order reduction stores the value of the
  //! <CODE>recurrence_rhs</CODE> before to solve the reduced recurrence;
  //! stores 0 otherwise.
  Expr old_recurrence_rhs_;

  //! \brief
  //! When is applied the order reduction stores the greatest common
  //! divisor among the decrements <CODE>d</CODE> of the terms
  //! <CODE>x(n-d)</CODE> present in the right hand side of the recurrence
  //! before to solve the reduced recurrence;
  //! stores 0 otherwise.
  unsigned gcd_decrements_old_rhs_;

  Symbol symbol_for_mod_;

  //! \brief
  //! When is applied the order reduction stores the solution of the
  //! reduced order recurrence; stores 0 otherwise.
  Expr solution_order_reduced_;

  bool verified_one_time_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Order_Reduction_Info.inlines.hh"

#endif // !defined(PURRS_Order_Reduction_Info_defs_hh)
