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
#include <giac/giac.h>

#include "factorize_giac.hh"

#include "purrs_install.hh"

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
    string s;
#if NOISY
    cout << endl << "Insert the hypergeometric term t(n): " << endl;
    cout << "(Blank line to use the hard-coded hypergeometric term binom(n,k)^2  )" << endl;
#endif
    getline(input_stream,s);
    if (s == "") {
      s = "((n!) / ((k!) * (n-k)!))^2";
    }
    Expr_List l(Recurrence::n, k, a, b, c, d);
    Expr t_n = Expr(s,l);
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
#endif
    }
    cout << endl << "---------------------------------------------" << endl;
      
  }
  return 0;
}
catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
  return 1;
}
