/* Matrix class implementation: inline functions.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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

#ifndef PURRS_Matrix_inlines_hh
#define PURRS_Matrix_inlines_hh

#include "Expr_List.defs.hh"
#include "Expr.defs.hh"

namespace Parma_Recurrence_Relation_Solver {

inline
Matrix::Matrix() {
}

inline
Matrix::Matrix(unsigned int i, unsigned int j)
  : m(i, j) {
};

inline
Matrix::Matrix(unsigned int i, unsigned int j, const Expr_List& y)
  : m(i, j, y.l) {
};

inline
Matrix::Matrix(const Matrix& y)
  : m(y.m) {
};

inline Matrix&
Matrix::operator=(const Matrix& y) {
  m = y.m;
  return *this;
};

inline
Matrix::Matrix(const GiNaC::matrix& gm)
  : m(gm) {
}

inline
Matrix::~Matrix() {
}

inline unsigned int
Matrix::num_rows() const {
  return m.rows();
};

inline unsigned int
Matrix::num_columns() const {
  return m.cols();
};


inline const Expr&
Matrix::operator()(unsigned int r, unsigned int c) const {
  return static_cast<const Expr&>(m(r, c));
};

inline Expr&
Matrix::operator()(unsigned int r, unsigned int c) {
  return static_cast<Expr&>(m(r, c));
};

inline Expr
Matrix::determinant() const {
  return m.determinant();
}

inline Matrix
Matrix::solve(const Matrix& vars, const Matrix& rhs) const {
  return m.solve(vars.m, rhs.m, GiNaC::solve_algo::gauss);
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Matrix_inlines_hh)
