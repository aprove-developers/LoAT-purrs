/* A class for containing all informations necessary for to solve
   non linear recurrences.
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

#ifndef PURRS_Non_Linear_Info_defs_hh
#define PURRS_Non_Linear_Info_defs_hh 1

#include "Non_Linear_Info.types.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Non_Linear_Info {
public:
  //! \brief
  //! Constructor: sets \f$ rhs_transformed_in_linear_ = new_rhs \f$,
  //! \f$ base_exp_log_ = base_exp_log \f$,
  //! \f$ auxiliary_symbols_ = auxiliary_symbols \f$.
  Non_Linear_Info(const Expr& new_rhs, const Expr& base_exp_log,
		  const std::vector<Symbol> auxiliary_symbols);

  //! Copy-constructor.
  Non_Linear_Info(const Non_Linear_Info& y);

  //! Destructor.
  ~Non_Linear_Info();

  //! Assignment operator.
  Non_Linear_Info& operator=(const Non_Linear_Info& y);

  //! Returns <CODE>rhs_transformed_in_linear_</CODE>.
  Expr rhs_transformed_in_linear() const;

  //! Returns <CODE>rhs_transformed_in_linear_</CODE>.
  Expr& rhs_transformed_in_linear();

  //! Returns <CODE>base_exp_log_</CODE>.
  Expr base_exp_log() const;

  //! Returns <CODE>base_exp_log_</CODE>.
  Expr& base_exp_log();

  //! Returns <CODE>auxiliary_symbols_</CODE>.
  const std::vector<Symbol>& auxiliary_symbols() const;
  
  //! Returns <CODE>auxiliary_symbols_</CODE>.
  std::vector<Symbol>& auxiliary_symbols();

  //! Returns <CODE>order_if_linear_</CODE>.
  unsigned int order_if_linear() const;

  //! Sets <CODE>order_if_non_linear_</CODE> with \p x.
  void set_order_if_linear(unsigned int x);

  //! Returns <CODE>non_linear_to_linear_fwdr_</CODE>.
  unsigned non_linear_to_linear_fwdr() const;

  //! Sets <CODE>non_linear_to_linear_fwdr_</CODE> with \p i_c
  void set_non_linear_to_linear_fwdr(unsigned i_c);

private:
  Expr rhs_transformed_in_linear_;

  //! \brief
  //! If the rewriting of the non-linear recurrence in a linear
  //! recurrence has success then this data contain the value that
  //! will be the logarithm's base or the exponential's base
  //! used in the rewriting.
  Expr base_exp_log_;

  //! \brief
  //! Stores the symbols associated to the eventual negative numbers
  //! that will be the arguments of the logarithms.
  std::vector<Symbol> auxiliary_symbols_;

  //! \brief
  //! When the non-linear recurrence is rewritable in a linear recurrence
  //! of finite order this data stores the order of the linear recurrence.
  unsigned int order_if_linear_;

  //! \brief
  //! When the non-linear recurrence is rewritable in a linear recurrence
  //! of finite order this data stores the smallest positive integer for
  //! which the recurrence is well-defined: the initial conditions will
  //! start from it.
  unsigned non_linear_to_linear_fwdr_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Non_Linear_Info.inlines.hh"

#endif // !defined(PURRS_Non_Linear_Info_defs_hh)
