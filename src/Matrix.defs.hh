/* Matrix class declaration.
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

#ifndef PURRS_Matrix_defs_hh
#define PURRS_Matrix_defs_hh 1

#include "Matrix.types.hh"
#include "Expr_List.types.hh"
#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

class Matrix {
public:
  //! Default constructor.
  Matrix();

  //! Builds the matrix with \p i rows and \p j columns.
  Matrix(unsigned i, unsigned j);

  //! \brief
  //! Builds the matrix with \p i rows and \p j columns containing the
  //! expression in \p lst.
  /*!
    If the list has fewer elements than the matrix, the remaining
    matrix elements are set to zero.
    If the list has more elements than the matrix, the excessive
    elements are thrown away.
  */
  Matrix(unsigned i, unsigned j, const Expr_List& y);

  //! Copy-constructor.
  Matrix(const Matrix& y);

  //! Destructor.
  ~Matrix();

  //! Assignment operator.
  Matrix& operator=(const Matrix& y);

  unsigned num_rows() const;

  unsigned num_columns() const;

  //! Accesses to element in \p r row and \p c column of \p *this for reading.
  /*!
    \exception std::range_error thrown if \p r is bigger or equal to
                                rows'number or \p c is bigger or equal
				columns'number.
  */
  const Expr& operator()(unsigned r, unsigned c) const;

  //! Accesses to element in \p r row and \p c column of \p *this for writing.
  /*!
    \exception std::range_error thrown if \p r is bigger or equal to
                                rows'number or \p c is bigger or equal
				columns'number.
  */
  Expr& operator()(unsigned r, unsigned c);

  //! \brief
  //! Solves a linear system consisting of a m x n matrix \p *this and a
  //! m x p right hand side \p rhs.
  /*!
    \param vars    n x p matrix, all elements must be symbols.
    \param rhs     m x p matrix.
    Returns a n x p solution matrix.

    \exception std::logic_error      thrown if the matrices'dimensions are
                                     incompatible.
    \exception std::invalid_argument thrown if \p vars is not a matrix of
                                     symbols.
    \exception std::runtime_error    thrown if the linear system is
                                     inconsistent. 
  */
  Matrix solve(const Matrix& vars, const Matrix& rhs) const;

private:
  friend class Expr;

  GiNaC::matrix m;

  //! Builds the matrix corresponding to \p gm.
  Matrix(const GiNaC::matrix& gm);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Matrix.inlines.hh"

#endif // !defined(PURRS_Matrix_defs_hh)
