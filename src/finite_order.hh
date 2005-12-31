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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
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
#include "globals.hh"
#include "alg_eq_solver.hh"
#include <vector>

namespace Parma_Recurrence_Relation_Solver {

bool
characteristic_equation_and_its_roots(index_type order,
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
		     unsigned lower_bound_sum,
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
compute_non_homogeneous_part(const Expr& g_n, index_type order,
			     const std::vector<Expr>& base_of_exps,
			     const std::vector<Expr>& exp_poly_coeff);

Expr
compute_sum_with_transcendental_method(const Number& lower, const Expr& upper,
				       const std::vector<Expr>& base_of_exps,
				       const std::vector<Expr>&
				       exp_no_poly_coeff,
				       const std::vector<Polynomial_Root>&
				       roots);

Expr
write_reduced_order_recurrence(const Expr& e, const Symbol& r,
			       unsigned gcd_among_decrements,
			       const std::vector<Expr>& coefficients,
			       std::vector<Expr>& new_coefficients,
			       Expr& inhomogeneous);

Expr 
come_back_to_original_variable(const Expr& e, const Symbol& r, const Expr& m,
			       unsigned gcd_among_decrements);

void
substitute_non_rational_roots(const Recurrence& rec,
			      std::vector<Polynomial_Root>& roots);

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_finite_order_hh)
