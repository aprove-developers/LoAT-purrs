/* Definition of the main recurrence relation solver.
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

#include "globals.hh"
#include "util.hh"
#include "alg_eq_solver.hh"

static GExpr
get_binding(const GList& l, unsigned wild_index) {
  assert(wild_index < l.nops());
  assert(l.op(wild_index).info(GiNaC::info_flags::relation_equal));
  assert(l.op(wild_index).lhs() == GiNaC::wild(wild_index));
  return l.op(wild_index).rhs();
}

static bool
get_linear_decrement(const GExpr& e, const GSymbol& n, GNumber& decrement) {
  static GExpr n_plus_d = n + GiNaC::wild(0);
  GList substitution;
  if (match(e, n_plus_d, substitution)) {
    GExpr d = get_binding(substitution, 0);
    if (GiNaC::is_a<GiNaC::numeric>(d)) {
      decrement = -GiNaC::ex_to<GiNaC::numeric>(d);
      return true;
    }
  }
  return false;
}

bool
solve(const GExpr& rhs, const GSymbol& n) {
  static GExpr x_i = x(GiNaC::wild(0));
  static GExpr x_i_plus_r = x_i + GiNaC::wild(1);
  static GExpr a_times_x_i = GiNaC::wild(1)*x_i;
  static GExpr a_times_x_i_plus_r = a_times_x_i + GiNaC::wild(2);

  static GList substitution;

  int degree = -1;
  vector<GNumber> coefficients;
  GExpr e = rhs;
  bool failed = false;
  do {
    GExpr i;
    GExpr a;
    // The following matches are attempted starting from the most common,
    // then the second most common and so forth.
    if (clear(substitution), match(e, x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = get_binding(substitution, 1);
    }
    else if (clear(substitution), match(e, a_times_x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = get_binding(substitution, 2);
    }
    else if (clear(substitution), match(e, a_times_x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = 0;
    }
    else if (clear(substitution), match(e, x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = 0;
    }
    else
      break;

    GNumber decrement;
    if (!get_linear_decrement(i, n, decrement)) {
      failed = true;
      break;
    }
    cout << "decrement = " << decrement << endl;
    if (!decrement.is_integer()
	|| decrement < 0
	|| decrement >= coefficients.max_size()) {
      failed = true;
      break;
    }
    if (!GiNaC::is_a<GiNaC::numeric>(a)) {
      failed = true;
      break;
    }
    GNumber coefficient = GiNaC::ex_to<GiNaC::numeric>(a);
    unsigned long index = decrement.to_long();
    if (degree < 0 || index > unsigned(degree))
      degree = index;
    if (index > coefficients.size())
      coefficients.insert(coefficients.end(),
			  index - coefficients.size(),
			  GNumber(0));
    if (index == coefficients.size())
      coefficients.insert(coefficients.end(), coefficient);
    else
      coefficients[index] += coefficient;
  } while (e != 0);

  if (e.has(x_i))
    failed = true;

  if (failed)
    return false;

  cout << "Degree = " << degree << endl;
  cout << "Coefficients = ";
  for (int i = 1; i <= degree; ++i)
    cout << coefficients[i] << " ";
  cout << endl;
  cout << "Inhomogeneous term = " << e << endl;

  // Build the expression here.
  static GSymbol x("x");
  GExpr p = 0;
  for (int i = 1; i <= degree; ++i)
    p += coefficients[i]*pow(x, i);
  
  std::vector<Polynomial_Root> roots;
  bool all_distinct;
  if (!find_roots(p, x, roots, all_distinct))
    return false;

  cout << degree << " " << roots.size() << " " << all_distinct << endl;
  return true;
}
