/* *****************
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

#ifndef PURRS_Number_defs_hh
#define PURRS_Number_defs_hh 1

#include "Number.types.hh"
#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Parma_Recurrence_Relation_Solver::Number : public Expr {
public:
  //! Ordinary copy-constructor.
  Number();

  //! Copy-constructor.
  Number(const Number& x);

  //! Destructor.
  ~Number();

  //! Assignment operator.
  Number& operator=(const Number& x);

  bool is_positive(const Number& n) const;
  bool is_integer(const Number& n) const;
  bool is_pos_integer(const Number& n) const;
  bool is_nonnes_integer(const Number& n) const;
  bool is_even(const Number& n) const;
  bool is_odd(const Number& n) const;
  bool is_prime(const Number& n) const;
  bool is_rational(const Number& n) const;
  bool is_real(const Number& n) const;
  bool is_cinteger(const Number& n) const;
  bool is_crational(const Number& n) const;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Number.inlines.hh"

#endif // !defined(PURRS_Number_defs_hh)
