/* Sum class declaration.
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

#ifndef PURRS_Sum_defs_hh
#define PURRS_Sum_defs_hh 1

#include "Sum.types.hh"
#include "Expr.defs.hh"
#include "Number.defs.hh"
#include "Symbol.defs.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

//! Returns \f$ x + y \f$.
Sum operator+(const Sum& x, const Sum& y);

//! Returns \f$ x - y \f$.
Sum operator-(const Sum& x, const Sum& y);

//! Returns \f$ x * y \f$.
Sum operator*(const Expr& x, const Sum& y);

//! Simplify sums when possible.
Expr simplify_sum(const Sum& x);

/*!
  We want to define the general summation symbol
  \f[
    \sum_{k = lower\_limit}^{\upper\_limit} f(k).
  \f]
  An object of the class Sum represents the sum above, and is specified by 
  its four arguments:
  - the first argument is a <CODE>Number</CODE>;
  - the second argument is an <CODE>Expr</CODE>;
  - the third argument is the summand \f$ f \f$, and is an <CODE>Expr</CODE>;
  - the fourth argument is a <CODE>Symbol</CODE>.
  We explicitly allow parameters in the third argument, and therefore we
  have to specify the symbol we are summing over.
*/
class Sum {
public:
  //! Default constructor.
  Sum();

  //! Construct empty sum.
  Sum(const Number& lower, const Expr& upper,
      const Expr& summ, const Symbol& ind);

  //! Copy-constructor.
  Sum(const Sum& y);

  //! Destructor.
  ~Sum();

  //! Ordinary assignment operator.
  Sum& operator=(const Sum& y);

private:
  friend Sum operator+(const Sum& x, const Sum& y);
  friend Sum operator-(const Sum& x, const Sum& y);
  friend Sum operator*(const Expr& x, const Sum& y);
  friend Expr simplify_sum(const Sum& x);

  //! Lower limit in the summation.
  Number lower_limit;

  //! Upper limit in the summation.
  Expr upper_limit;

  //! The summand.
  Expr summand;

  //! The symbol we are summing over.
  Symbol index;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Sum.inlines.hh"

#endif // !defined(PURRS_Sum_defs_hh)
