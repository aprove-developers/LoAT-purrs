/* Cached_Expr class implementation: inline functions.
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

#ifndef PURRS_Cached_Expr_inlines_hh
#define PURRS_Cached_Expr_inlines_hh

namespace Parma_Recurrence_Relation_Solver {

inline
Cached_Expr::Cached_Expr()
  : has_expression_(false),
    expression_(),
    has_size_(false),
    size_(0) {
}

inline
Cached_Expr::Cached_Expr(const Cached_Expr& y)
  : has_expression_(y.has_expression_),
    expression_(y.expression_),
    has_size_(y.has_size_),
    size_(y.size_) {
}

inline
Cached_Expr::~Cached_Expr() {
}

inline Cached_Expr&
Cached_Expr::operator=(const Cached_Expr& y) {
  has_expression_ = y.has_expression_;
  expression_ = y.expression_;
  has_size_ = y.has_size_;
  size_ = y.size_;
}

inline bool
Cached_Expr::has_expression() const {
  return has_expression_;
}

inline const Expr&
Cached_Expr::expression() {
  return expression_;
}

inline void
Cached_Expr::set_expression(const Expr& e) {
  expression_ = e;
  has_expression_ = true;
}

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_Cached_Expr_inlines_hh)

