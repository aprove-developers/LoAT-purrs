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
#include <map>

namespace Parma_Recurrence_Relation_Solver {

/*!
  ...
*/
class Functional_Equation_Info {
public:
  Functional_Equation_Info(const std::map<Number, Expr>& hom_terms,
			   unsigned c = 1);

  //! Copy-constructor.
  Functional_Equation_Info(const Functional_Equation_Info& y);

  //! Destructor.
  ~Functional_Equation_Info();

  //! Assignment operator.
  Functional_Equation_Info& operator=(const Functional_Equation_Info& y);

  //! Returns <CODE>applicability_condition_</CODE>.
  unsigned applicability_condition() const;

  //! Returns <CODE>applicability_condition_</CODE>.
  unsigned& applicability_condition();

  //! Sets <CODE>applicability_condition_</CODE> with \p c.
  void set_applicability_condition(unsigned c);

private:
  typedef std::map<Number, Expr> Homogeneous_Terms; 
  
public:
  typedef Homogeneous_Terms::iterator ht_iterator;
  typedef Homogeneous_Terms::const_iterator ht_const_iterator;

  ht_iterator ht_begin();
  ht_iterator ht_end();

  ht_const_iterator ht_begin() const;
  ht_const_iterator ht_end() const;

  //! \brief
  //! Returns the rank of the functional equation, i. e., the number of terms
  //! of the form \f$ a x(n/b) \f$ where \f$ b \f$ is a rational number
  //! larger than one.
  size_t rank() const;
    
private:
  //! \brief
  //! Stores the divisor \f$ b \f$ and the coefficient \f$ a \f$ of the
  //! homogeneous terms \f$ a x(n/b) \f$ in \p *this. 
  std::map<Number, Expr> homogeneous_terms;

  //! \brief
  //! The positive integer starting from which the inhomogeneous term
  //! is a non negative, non decreasing function. 
  unsigned applicability_condition_;

public:
  void dump_homogeneous_terms(std::ostream& s) const;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Functional_Equation_Info.inlines.hh"

#endif // !defined(PURRS_Functional_Equation_Info_defs_hh)
