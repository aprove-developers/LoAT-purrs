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

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

using namespace std;
using namespace GiNaC;

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

  GSymbol n("n");
  GSymbol a("a");
  GSymbol b("b");
  GSymbol c("c");
  GSymbol d("d");
  GList symbols(n, a, b, c, d);
  GExpr rhs;
  while (INPUT_STREAM) {
    string s;
    getline(INPUT_STREAM, s);

    if (!INPUT_STREAM)
      return 0;

    // Skip comments.
    if (s.find("%") == 0)
      continue;

    rhs = GExpr(s, symbols);

    // The recurrence relations are passed  to the function solve
    // only expanded, i.e., in the form
    // x(n-1)*a+x(n-2)*b+...+x(n-k)*h+p(n) with p(n) also expanded.
    rhs = rhs.expand();
#if NOISY
    cout << "Trying to solve x(n) = " << rhs << endl;
#endif

    GExpr solution;
    if (!solve(rhs, n, solution)) {
#if NOISY
      cout << "Sorry, this is too difficult." << endl;
#endif
    }
#if NOISY
    else {
      cout << "*** SOLUTION ***"
	   << endl
	   << solution
	   << endl
	   << "****************"
	   << endl << endl;
    }
#endif
  }
  return 0;
} 
catch (exception &p) {
#if NOISY
  cerr << "Exception caught: " << p.what() << endl;
#endif
}
catch (const char* s) {
  cerr << s << endl;
}
