/* A class for containing the necessary informations about functional
   equations (and that we do not want to compute them again).
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
  //! Constructor: sets \f$ rank{\_} = k \f$,
  //! \f$ coefficients{\_}fe{\_} = coeffs \f$,
  //! \f$ divisors{\_}arg{\_} = divisors \f$
  //! and \f$ applicability{\_}condition{\_} = c \f$
  Functional_Equation_Info(unsigned int k, const std::vector<Expr>& coeffs,
			   const std::vector<Number>& divisors,
			   unsigned c = 1);

  //! Copy-constructor.
  Functional_Equation_Info(const Functional_Equation_Info& y);

  //! Destructor.
  ~Functional_Equation_Info();

  //! Assignment operator.
  Functional_Equation_Info& operator=(const Functional_Equation_Info& y);

  //! Returns <CODE>rank_</CODE>.
  unsigned int rank() const;

  //! Returns <CODE>rank_</CODE>.
  unsigned int& rank();

  //! Returns <CODE>coefficients_fe_</CODE>.
  const std::vector<Expr>& coefficients_fe() const;

  //! Returns <CODE>coefficients_fe_</CODE>.
  std::vector<Expr>& coefficients_fe();

  //! Returns <CODE>divisors_arg_</CODE>.
  const std::vector<Number>& divisors_arg() const;

  //! Returns <CODE>divisors_arg_</CODE>.
  std::vector<Number>& divisors_arg();

  //! Returns <CODE>applicability_condition_</CODE>.
  unsigned applicability_condition() const;

  //! Returns <CODE>applicability_condition_</CODE>.
  unsigned& applicability_condition();

  //! Sets <CODE>applicability_condition_</CODE> with \p c.
  void set_applicability_condition(unsigned c);

private:
  //! \brief
  //! The rank of the functional equation, i. e., the number of terms
  //! of the form \f$ a x(n/b) \f$ where \f$ b \f$ is a rational number
  //! larger than one.
  unsigned int rank_;

  //! Stores the coefficients of the functional equation.
  std::vector<Expr> coefficients_fe_;

  //! \brief
  //! Stores the divisors \f$ b \f$ of the argument of the function
  //! \f$ x \f$ in the terms of the form \f$ a x_{n/b} \f$ contained in
  //! the functional equation.
  std::vector<Number> divisors_arg_;

  //! \brief
  //! The positive integer starting from which the inhomogeneous term
  //! is a non negative, non decreasing function. 
  unsigned applicability_condition_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Functional_Equation_Info.inlines.hh"

#endif // !defined(PURRS_Functional_Equation_Info_defs_hh)
