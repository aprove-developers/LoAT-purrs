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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#ifndef PURRS_Functional_Equation_Info_defs_hh
#define PURRS_Functional_Equation_Info_defs_hh 1

#include "Functional_Equation_Info.types.hh"
#include "globals.hh"
#include "Expr.defs.hh"
#include <map>
#include <string>

namespace Parma_Recurrence_Relation_Solver {

/*!
  ...
*/
class Functional_Equation_Info {
public:
  Functional_Equation_Info(const std::map<Number, Expr>& hom_terms);

  //! Copy-constructor.
  Functional_Equation_Info(const Functional_Equation_Info& y);

  //! Destructor.
  ~Functional_Equation_Info();

  //! Assignment operator.
  Functional_Equation_Info& operator=(const Functional_Equation_Info& y);

  //! Returns <CODE>applicability_condition_</CODE>.
  index_type applicability_condition() const;

  //! Sets <CODE>applicability_condition_</CODE> with \p c.
  void set_applicability_condition(index_type c);

private:
  //! \brief
  //! Contains the two elements of the terms in the form \f$ a x(n/b) \f$:
  //! in the first position there is \f$ b \f$; in the second position there
  //! is \f$ a \f$.
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
  index_type rank() const;

  //! Returns <CODE>definition_Sc_</CODE>.
  std::string definition_Sc() const;
    
  //! \brief
  //! If the divisor \f$ b \f$ of the homogeneous term \f$ a x(n/b) \f$
  //! in \p *this is different from \f$ 2 \f$ then set
  //! <CODE>definition_Sc_</CODE> with the definition of the function
  //! \f$ Sc() \f$.
  void set_definition_Sc();

private:
  //! \brief
  //! Stores the divisor \f$ b \f$ and the coefficient \f$ a \f$ of the
  //! homogeneous terms \f$ a x(n/b) \f$ in \p *this. 
  std::map<Number, Expr> homogeneous_terms;

  //! \brief
  //! The positive integer starting from which the inhomogeneous term
  //! is a non negative, non decreasing function. 
  index_type applicability_condition_;

  //! \brief
  //! If the divisor \f$ b \f$ of the homogeneous term \f$ a x(n/b) \f$
  //! in \p *this is \f$ 2 \f$ then is the empty string; otherwise
  //! the string conatins the definition of the function \f$ Sc() \f$:
  //! \f[
  //!   Sc(n, b) = \lfloor
  //!                \frac{n}{b^{\lfloor \frac{\log n}{\log b} \rfloor}}
  //!              \rfloor.
  //! \f]
  std::string definition_Sc_;
public:
  void dump_homogeneous_terms(std::ostream& s) const;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Functional_Equation_Info.inlines.hh"

#endif // !defined(PURRS_Functional_Equation_Info_defs_hh)
