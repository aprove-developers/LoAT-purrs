/* Simple program to test the Gosper's algorithm.
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

#include <giac/giac.h>
#include "factorize_giac.hh"

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

using namespace std;
using namespace giac;
using namespace Parma_Recurrence_Relation_Solver;

namespace PURRS = Parma_Recurrence_Relation_Solver;

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

  Symbol x("x");
  Symbol a("a");
  Symbol b("b");
  Symbol c("c");
  Symbol d("d");
  Symbol k("k");

  while (input_stream) {
    Expr solution = 0;
    //    string default_term = "binom(n,k)^2";
    string default_term = "(-1)^(k)*(binom(2*n,k))^3";
    string s;
#if NOISY
    cout << endl << "Insert the hypergeometric term t(n): " << endl;
    cout << "(Blank line to use the hard-coded hypergeometric term " << default_term << " )" << endl;
#endif
    getline(input_stream,s);

    // We may be at end of file.
    if (!input_stream)
      break;
      
    // Skip comments.
    if (s.find("%") == 0)
      continue;

    if (s == "") {
      s = default_term;
    }
    Expr_List l(Recurrence::n, k, a, b, c, d);
    Expr t_n;
    try {
      t_n = Expr(s,l);
    cout << t_n << endl;

    if (zeilberger_algorithm(t_n, Recurrence::n, k, solution))
#if NOISY
      std::cout << endl << "The sum is: " << solution << std::endl;
#endif
    else {
      Symbol h;
#if NOISY
      std::cout << endl << "Error in Zeilberger Algorithm for term " << t_n
		<< std::endl;
    }
#endif
    }
    catch (exception& e) {
      std::cerr << "std::exception caught: " << e.what();
    }
    
    cout << endl << "---------------------------------------------" << endl;
      
  }
  return 0;
}
catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
  return 1;
}
