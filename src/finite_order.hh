/* To be written.
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

#ifndef PURRS_finite_order_hh
#define PURRS_finite_order_hh 1

#include "Expr.types.hh"
#include "Number.types.hh"
#include "Symbol.types.hh"
#include "Matrix.types.hh"
#include "Recurrence.types.hh"
#include "alg_eq_solver.hh"
#include <vector>

namespace Parma_Recurrence_Relation_Solver {

bool
characteristic_equation_and_its_roots(unsigned int order,
				      const std::vector<Expr>& coefficients,
				      std::vector<Number>& num_coefficients,
				      Expr& characteristic_eq,
				      std::vector<Polynomial_Root>& roots,
				      bool& all_distinct);

void
compute_symbolic_sum(const Symbol& alpha, const Symbol& lambda,
		     const std::vector<Polynomial_Root>& roots,
		     const std::vector<Expr>& base_of_exps,
		     const std::vector<Expr>& exp_poly_coeff,
		     std::vector<Expr>& symbolic_sum_distinct,
		     std::vector<Expr>& symbolic_sum_no_distinct);

Expr
subs_to_sum_roots_and_bases(const Symbol& alpha, const Symbol& lambda,
			    const std::vector<Polynomial_Root>& roots,
			    const std::vector<Expr>& base_of_exps,
			    std::vector<Expr>& symbolic_sum_distinct,
			    std::vector<Expr>& symbolic_sum_no_distinct);

Matrix
solve_system(bool all_distinct,
	     const std::vector<Number>& coefficients,
	     const std::vector<Polynomial_Root>& roots);

Expr
find_g_n(bool all_distinct, const Matrix& sol,
	 const std::vector<Polynomial_Root>& roots);

void
prepare_for_symbolic_sum(const Expr& g_n,
			 const std::vector<Polynomial_Root>& roots,
			 const std::vector<Expr>& exp_poly_coeff,
			 std::vector<Expr>& poly_coeff_tot);

Expr
compute_non_homogeneous_part(const Expr& g_n, unsigned int order,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff);

Expr
rewrite_reduced_order_recurrence(const Expr& e, const Symbol& r,
				 int gcd_among_decrements);

Expr 
come_back_to_original_variable(const Expr& e, const Symbol& r, const Expr& m,
			       int gcd_among_decrements);

void
substitute_non_rational_roots(const Recurrence& rec,
			      std::vector<Polynomial_Root>& roots);

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_finite_order_hh)
