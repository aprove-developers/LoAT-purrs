/* Simple program to test the algebraic equation solver.
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

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cassert>

#include "purrs_install.hh"

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

#define NOISY 1

#ifdef USE_READLINE
#define INPUT_STREAM rdl
#else
#define INPUT_STREAM cin
#endif

int
main() try {
#ifdef USE_READLINE
  auto_ptr<readlinebuf> rdlb(new readlinebuf());
  istream rdl(rdlb.get());
#endif

  Symbol x("x");
  Expr_List symbols(x);
  Expr p;
  while (INPUT_STREAM) {
    string s;
    getline(INPUT_STREAM, s);

    if (!INPUT_STREAM)
      return 0;

    // Skip comments.
    if (s.find("%") == 0)
      continue;

    if (s.find("symbols ") == 0) {
      for (int i = 7, l = s.length(); i < l; ++i) {
	char c = s[i];
	if (c >= 'a' && c <= 'z') {
	  char name[2];
	  name[0] = c;
	  name[1] = '\0';
	  Symbol new_symbol(name);
	  symbols.append(new_symbol);
	}
      }
      continue;
    }

    p = Expr(s, symbols);
    if (p.is_zero())
      continue;

#if NOISY
    cout << std::endl << "Trying to solve " << p << " = 0" << endl;
#endif

    std::vector<Polynomial_Root> roots;
    bool all_distinct;
    if (!find_roots(p, x, roots, all_distinct)) {
#if NOISY
      cout << "Sorry, this is too difficult." << endl;
#endif
    }
#if NOISY
    else {
      size_t n = roots.size();
      for (size_t i = 0; i < n; ++i) {
	Expr value = roots[i].value();
	Number multiplicity = roots[i].multiplicity();
	cout << "x_" << i+1 << " = " << value;
	if (multiplicity > 1)
	  cout << " (multiplicity " << multiplicity << ")";
	cout << endl;
	if (!value.is_a_number())
#if 0
	  cout << "****  x_" << i+1 << " ~= " << value.evalf() << endl;
#else
	  cout << "****  x_" << i+1 << " ~= " << approximate(value) << endl;
#endif
      }
    }
#endif
  }
  return 0;
} catch (exception &p) {
#if NOISY
  cerr << "Exception caught: " << p.what() << endl;
#endif
  return 1;
}
