/* A class for containing all informations necessary for to solve
   weighted-average recurrences.
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

#ifndef PURRS_Weighted_Average_Info_defs_hh
#define PURRS_Weighted_Average_Info_defs_hh 1

#include "Weighted_Average_Info.types.hh"
#include "globals.hh"
#include "Expr.defs.hh"
#include "Recurrence.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Weighted_Average_Info {
public:
  //! \brief
  //! Constructor: sets
  //! \f$ associated_first_order_rec_ = associated_first_order_rec \f$,
  //! \f$ weight_ = weight \f$.
  Weighted_Average_Info(const Recurrence& associated_first_order_rec,
			const Expr& weight);

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

  //! Sets <CODE>original_rhs_</CODE> with \p original_rhs.
  void set_original_rhs(const Expr& original_rhs);

  //! Returns <CODE>lower_limit_</CODE>.
  unsigned int lower_limit() const;

  //! Sets <CODE>lower_limit_</CODE> with \p lower.
  void set_lower_limit(unsigned int lower);
  
  //! Returns <CODE>weight_</CODE>.
  const Expr& weight() const;

  //! Returns <CODE>weight_</CODE>.
  Expr& weight();

private:
  //! \brief
  //! Contains the first order recurrence associated to \p *this,
  //! if it is possible to obtain it (in order to know the cases of
  //! rewritable weighted-average recurrences see the function
  //! <CODE>rewrite_weighted_average_recurrence()</CODE>).
  Recurrence associated_first_order_rec_;

  //! \brief
  //! Contains the right hand side of the recurrence before the rewriting
  //! of the system in a weighted-average recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n).
  //! \f]
  //! If the recurrence is already in form of weighted-average recurrence
  //! then this data is undefined.
  Expr original_rhs_;

  //! \brief
  //! Contains the factor \f$ f(n) \f$ of the weighted-average recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n).
  //! \f]
  Expr weight_;

  //! \brief
  //! Contains the lower limit of the sum \f$ n_0 \f$ of a recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n).
  //! \f]
  //! rewritable in a weighted-average recurrence.
  unsigned int lower_limit_;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Weighted_Average_Info.inlines.hh"

#endif // !defined(PURRS_Weighted_Average_Info_defs_hh)
