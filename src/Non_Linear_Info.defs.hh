/* A class for containing all informations necessary for to solve
   non linear recurrence.
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
  //! Constructor: sets \f$ original_recurrence_rhs_ = rhs \f$,
  //! \f$ base_exp_log_ = base_exp_log \f$.
  Non_Linear_Info(const Expr& rhs, const Expr& base_exp_log);

  //! Copy-constructor.
  Non_Linear_Info(const Non_Linear_Info& y);

  //! Destructor.
  ~Non_Linear_Info();

  //! Assignment operator.
  Non_Linear_Info& operator=(const Non_Linear_Info& y);

  Expr original_recurrence_rhs() const;
  Expr& original_recurrence_rhs();

  Expr base_exp_log() const;
  Expr& base_exp_log();
 
private:
  //!
  Expr original_recurrence_rhs_;

  //!
  Expr base_exp_log_;

};

} // namespace Parma_Recurrence_Relation_Solver

#include "Non_Linear_Info.inlines.hh"

#endif // !defined(PURRS_Non_Linear_Info_defs_hh)
