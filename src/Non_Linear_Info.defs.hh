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
#include "globals.hh"
#include "Expr.defs.hh"
#include "Symbol.defs.hh"
#include "Recurrence.defs.hh"

#include<utility>

namespace Parma_Recurrence_Relation_Solver {

class Non_Linear_Info {
public:
  //! \brief
  //! Constructor: sets \f$ associated_linear_rec_ = associated_linear_rec \f$,
  //! \f$ coeff_and_base_ = coeff_and_base \f$,
  //! \f$ auxiliary_symbols_ = auxiliary_symbols \f$.
  Non_Linear_Info(const Recurrence& associated_linear_rec,
		  const std::pair<Number, Expr>& coeff_and_base,
		  const std::vector<Symbol> auxiliary_symbols);

  //! Copy-constructor.
  Non_Linear_Info(const Non_Linear_Info& y);

  //! Destructor.
  ~Non_Linear_Info();

  //! Assignment operator.
  Non_Linear_Info& operator=(const Non_Linear_Info& y);

  //! Returns <CODE>associated_linear_rec_</CODE>.
  const Recurrence& associated_linear_rec() const;

  //! Returns <CODE>associated_linear_rec_</CODE>.
  Recurrence& associated_linear_rec();

  //! Returns <CODE>coeff_and_base_.first</CODE>.
  const Number& coeff_simple_non_linear_rec() const;

  //! Returns <CODE>coeff_and_base_.first</CODE>.
  Number& coeff_simple_non_linear_rec();

  //! Returns <CODE>coeff_and_base_.second</CODE>.
  const Expr& base_exp_log() const;

  //! Returns <CODE>coeff_and_base_.second</CODE>.
  Expr& base_exp_log();

  //! Returns <CODE>auxiliary_symbols_</CODE>.
  const std::vector<Symbol>& auxiliary_symbols() const;
  
  //! Returns <CODE>auxiliary_symbols_</CODE>.
  std::vector<Symbol>& auxiliary_symbols();

private:
  //! \brief
  //! In the case which the system is able to rewrite the non-linear
  //! recurrence \p *this in linear, this method stores the linear
  //! recurrence computed (in order to know the cases of rewritable
  //! non-linear recurrences see the function
  //! <CODE>rewrite_non_linear_recurrence()</CODE>).
  Recurrence associated_linear_rec_;

  //! \brief
  //! \p coeff_and_base_ is used in two different ways:
  //! In the case of simple non-linear recurrence of the form
  //! \f$ x(n) = c x(n-1)^{\alpha} \f$ it contains the pair \f$ c, \alpha \f$;
  //! in all the other cases the numeric value of the pair holds \f$ 0 \f$
  //! while the second element contains the value that will be the
  //! logarithm's base or the exponential's base used in the rewriting
  //! of the non-linear recurrence in the correspondent linear recurrence.
  std::pair<Number, Expr> coeff_and_base_;

  //! \brief
  //! Stores the symbols associated to the eventual negative numbers
  //! that will be the arguments of the logarithms.
  std::vector<Symbol> auxiliary_symbols_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Non_Linear_Info.inlines.hh"

#endif // !defined(PURRS_Non_Linear_Info_defs_hh)
