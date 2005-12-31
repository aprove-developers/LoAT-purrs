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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
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

//! \brief
//! Contains all the components necessary to rebuilding the
//! original recurrence before the transformation of it in
//! weighted-average recurrence.
struct Components_Original_Rec {
  //! \brief
  //! Contains the factor \f$ f(n) \f$ of the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ n_0 \in \Nset \cup \{ 0 \} \f$ and
  //! \f$ u(n) \in \{ n-1, n \} \f$.
  Expr weight_;

  //! \brief
  //! Contains the inhomogeneous term \f$ g(n) \f$ of the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ n_0 \in \Nset \cup \{ 0 \} \f$ and
  //! \f$ u(n) \in \{ n-1, n \} \f$.
  Expr inhomogeneous_;

  //! \brief
  //! Contains the lower limit \f$ n_0 \in \Nset \cup \{ 0 \} \f$
  //! of the sum in the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ u(n) \in \{ n-1, n \} \f$.
  unsigned int lower_limit_;  

  //! \brief
  //! Contains the upper limit \f$ u(n) \in \{ n-1, n \} \f$
  //! of the sum in the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ n_0 \in \Nset \cup \{ 0 \} \f$.
  Expr upper_limit_;
};

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

  //! Returns <CODE>weight_</CODE>.
  const Expr& weight() const;

  //! Returns <CODE>weight_</CODE>.
  Expr& weight();

  //! Returns <CODE>original_weight_</CODE>.
  const Expr& original_weight() const;

  //! Returns <CODE>Components_Original_Rec.original_weight_</CODE>.
  Expr& original_weight();

  //! Returns <CODE>Components_Original_Rec.original_inhomogeneous_</CODE>.
  const Expr& original_inhomogeneous() const;

  //! Returns <CODE>Components_Original_Rec.original_inhomogeneous_</CODE>.
  Expr& original_inhomogeneous();

  //! Returns <CODE>Components_Original_Rec.lower_limit_</CODE>.
  unsigned int lower_limit() const;

  //! \brief
  //! Sets with \p weight, \p inhomogeneous, \p lower and \p upper
  //! the elements of <CODE>Components_Original_Rec</CODE>.
  void set_original_rhs(const Expr& weight, const Expr& inhomogeneous,
			unsigned int lower, const Expr& upper);

private:
  //! \brief
  //! Contains the first order recurrence associated to \p *this,
  //! if it is possible to obtain it (in order to know the cases of
  //! rewritable weighted-average recurrences see the function
  //! <CODE>rewrite_weighted_average_recurrence()</CODE>).
  Recurrence associated_first_order_rec_;

  //! \brief
  //! Contains the factor \f$ f(n) \f$ of the weighted-average recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=0}^{n-1} x(k) + g(n).
  //! \f]
  Expr weight_;

  //! \brief
  //! Pointer to a structure which contains all the components necessary
  //! to rebuilding the original recurrence before the transformation of
  //! it in weighted-average recurrence.
  Components_Original_Rec* components_original_rec_p;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Weighted_Average_Info.inlines.hh"

#endif // !defined(PURRS_Weighted_Average_Info_defs_hh)
