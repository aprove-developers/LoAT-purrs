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

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

#ifndef NOISY
#define NOISY 1
#endif

int main() try {
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
    Expr_List l(x, n, a, b, c, d);
    string s;
    Expr tmp;
    do {
      getline(input_stream, s);
      if (!input_stream)
	return 0;
      tmp = Expr(s, l);
    } while (!tmp.is_equal(1) && !tmp.is_equal(2) && !tmp.is_equal(3)); 
    Number choice;
    choice = tmp.ex_to_number();
    cout << endl << "Insert an expression: ";
    getline(input_stream, s);
    
    if (!input_stream)
	return 0;
    
    // Skip comments.
    if (s.find("%") == 0)
      continue;
    
    e = Expr(s, l);
    // `expand' does the simplification of the rule E6.
    e = e.expand();
#if NOISY
    cout << "Expanded expression = " << e << endl;
#endif
    Expr solution;
    switch (choice.to_int()) {
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
