/* Simple program to test the recurrence relation solver.
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
#include <cassert>

#include "purrs_install.hh"

using namespace std;
using namespace GiNaC;

#define NOISY 1

int
main() try {
  GSymbol n("n");
  GList symbols(n);
  GExpr rhs;
  while (cin) {
    string s;
    getline(cin, s);

    if (!cin)
      return 0;

    // Skip comments.
    if (s.find("%") == 0)
      continue;

    rhs = GExpr(s, symbols);

#if NOISY
    cout << "Trying to solve x(n) = " << rhs << endl;
#endif

    if (!solve(rhs, n)) {
#if NOISY
      cout << "Sorry, this is too difficult." << endl;
#endif
    }
#if NOISY
    else {
      cout << "Print solutions here." << endl;
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