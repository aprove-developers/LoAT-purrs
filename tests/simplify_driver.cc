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
#include <cctype>

#include "purrs_install.hh"
#include "ehandlers.hh"

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

#ifndef NOISY
#define NOISY 1
#endif

bool
all_space(const string& s) {
  for (unsigned i = s.length(); i-- > 0; )
    if (!isspace(s[i]))
      return false;
  return true;
}

int
main() try {
  readlinebuf* prdlb = 0;
  istream* pinput_stream;

#ifdef USE_READLINE
  prdlb = new readlinebuf();
  pinput_stream = new istream(prdlb);
#else
  pinput_stream = &cin;
#endif
  istream& input_stream = *pinput_stream;

  Expr e;
  Symbol x("x");
  Symbol n("n");
  Symbol a("a");
  Symbol b("b");
  Symbol c("c");
  Symbol d("d");

  while (input_stream) {
    cout << endl;
    cout << "What kind of expression do you want to simplify?" << endl;
    cout << "1. For the input (to put in evidence the symbol `n')" << endl;
    cout << "2. For the output" << endl;
    cout << "3. Containing factorials" << endl;
    int choice = 0;
    do {
      if (!input_stream)
	return 0;
      input_stream >> choice;
    } while (choice < 1 || choice > 3); 

    cout << endl << "Insert an expression: ";
    Expr_List l(x, n, a, b, c, d);

    string s;
    do {
      getline(input_stream, s);
    } while (all_space(s));

    
    if (!input_stream)
	return 0;
    
    e = Expr(s, l);
    // `expand' does the simplification of the rule E6.
    e = e.expand();
#if NOISY
    cout << "Expanded expression = " << e << endl;
#endif
    Expr solution;
    switch (choice) {
    case 1:
      solution = simplify_on_input_ex(e, n, true);
      break;
    case 2:
      solution = simplify_on_output_ex(e, n, false);
      break;
    case 3:
      solution = simplify_factorials_and_exponentials(e, n);
      break;
    }
#if NOISY
    cout << endl << "SOLUTION = " << solution << endl;
    cout << endl << "---------------------------------------------" << endl;
#endif
  }
  return 0;
}
catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
  return 1;
}
