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

#ifndef PURRS_Functional_Equation_Info_defs_hh
#define PURRS_Functional_Equation_Info_defs_hh 1

#include "Functional_Equation_Info.types.hh"
#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

/*!
  x_n = a x_{n/b} + d n^e
*/
class Functional_Equation_Info {
public:
  //! \brief
  //! Constructor: sets \f$ coefficient_ = a \f$, and
  //! \f$ divisor_arg_ = b \f$.
  Functional_Equation_Info(const Expr& a, unsigned b);

  //! Copy-constructor.
  Functional_Equation_Info(const Functional_Equation_Info& y);

  //! Destructor.
  ~Functional_Equation_Info();

  //! Assignment operator.
  Functional_Equation_Info& operator=(const Functional_Equation_Info& y);

  //! Returns <CODE>coefficient_</CODE>.
  Expr coefficient() const;

  //! Returns <CODE>coefficient_</CODE>.
  Expr& coefficient();

  //! Returns <CODE>divisor_arg_</CODE>.
  unsigned divisor_arg() const;

  //! Returns <CODE>divisor_arg_</CODE>.
  unsigned& divisor_arg();

private:
  //! \brief
  //! Stores the coefficient \f$ a \f$ of the equation
  //! \f$ x_n = a x_{n/b} + d n^e \f$. 
  Expr coefficient_;

  //! \brief
  //! Stores the divisor\f$ b \f$ of the argument of the function
  //! \f$ x \f$ in the equation \f$ x_n = a x_{n/b} + d n^e \f$.
  unsigned divisor_arg_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Functional_Equation_Info.inlines.hh"

#endif // !defined(PURRS_Functional_Equation_Info_defs_hh)
