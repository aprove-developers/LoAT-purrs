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

  //! Builds the matrix with \p i rows and \j columns.
  Matrix(unsigned i, unsigned j);

  //! \brief
  //! Builds the matrix with \p i rows and \j columns containing the
  //! expression in \p lst.
  Matrix(unsigned i, unsigned j, const Expr_List& y);

  //! Copy-constructor.
  Matrix(const Matrix& y);

  //! Destructor.
  ~Matrix();

  //! Assignment operator.
  Matrix& operator=(const Matrix& y);

  //! Accesses to element in \p r row and \p c column of \p *this for reading.
  const Expr& operator()(unsigned r, unsigned c) const;

  //! Accesses to element in \p r row and \p c column of \p *this for writing.
  Expr& operator()(unsigned r, unsigned c);

  //! \brief
  //! Solves a linear system consisting of a m x n matrix \p *this and a
  //! m x p right hand side \p rhs.
  /*!
    \param vars    n x p matrix, all elements must be symbols.
    \param rhs     m x p matrix.

    Returns a n x p solution matrix.
  */
  Matrix solve(const Matrix& vars, const Matrix& rhs) const;

private:
  GiNaC::matrix m;

  friend class Expr;

  //! Builds the matrix corresponding to \p gm.
  Matrix(const GiNaC::matrix& gm);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Matrix.inlines.hh"

#endif // !defined(PURRS_Matrix_defs_hh)
