/* Simple program to test the simplifier.
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
#include <stdexcept>
#include <string>

#include "purrs_install.hh"

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

#ifndef NOISY
#define NOISY 1
#endif

int main() try {
  Expr p;
  Symbol x("x");
  Symbol n("n");
  Symbol a("a");
  Symbol b("b");
  Symbol c("c");
  Symbol d("d");
  string s;

  for (int i = 0; i < 10; ++i) {
    cout << "---------------------------------------------" << endl;
    cout << endl << "Insert an expression: ";
    getline(cin,s);
    Expr_List l(x, n, a, b, c, d);
    p = Expr(s, l);
    // `expand' does the simplification of the rule E6.
    p = p.expand();
#if NOISY
    cout << "p after expand: " << p << endl;
#endif
    // Da fare sulle espressioni in input. 
    //Expr q = simplify_on_input_ex(p, n, true);
    // Da fare sulle espressioni in output.
    Expr q = simplify_on_output_ex(p, n, false);
#if NOISY
    cout << endl << "SOLUZIONE    " << q << endl;
#endif
  }
  return 0;
}
catch (exception &p) {
#if NOISY
  cerr << "Exception caught: " << p.what() << endl;
#endif
}
