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
#include <getopt.h>

#include "purrs_install.hh"
#include "ehandlers.hh"

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

static struct option long_options[] = {
  {"help",           no_argument,       0, 'h'},
  {"interactive",    no_argument,       0, 'i'},
  {0, 0, 0, 0}
};

static const char* usage_string
= "Usage: %s [OPTION]...\n\n"
"  -h, --help              prints this help text\n"
"  -i, --interactive       sets interactive mode on\n";

#define OPTION_LETTERS "hi"

// Interactive mode is on when true.
static bool interactive = false;

static void
my_exit(int status) {
  //(void) purrs_finalize();
  exit(status);
}

static void
process_options(int argc, char* argv[]) {
  int option_index;
  int c;

  while (1) {
    option_index = 0;
    c = getopt_long(argc, argv, OPTION_LETTERS, long_options, &option_index);
    if (c == EOF)
      break;

    switch (c) {
    case 0:
      break;

    case '?':
    case 'h':
      fprintf(stderr, usage_string, argv[0]);
      my_exit(0);
      break;

    case 'i':
      interactive = 1;
      break;

    default:
      abort();
    }
  }

  if (optind < argc) {
    fprintf(stderr, usage_string, argv[0]);
    my_exit(1);
  }
}

int
main(int argc, char *argv[]) try {
  set_handlers();

  process_options(argc, argv);

  readlinebuf* prdlb = 0;
  istream* pinput_stream;

#ifdef USE_READLINE
  if (interactive) {
    prdlb = new readlinebuf();
    pinput_stream = new istream(prdlb);
  }
  else
    pinput_stream = &cin;
#else
  pinput_stream = &cin;
#endif
  istream& input_stream = *pinput_stream;

  Symbol n("n");
  Symbol a("a");
  Symbol b("b");
  Symbol c("c");
  Symbol d("d");
  Expr_List symbols(n, a, b, c, d);
  Expr rhs;
  while (input_stream) {
    string s;
    getline(input_stream, s);

    if (!input_stream)
      return 0;

    // Skip comments.
    if (s.find("%") == 0)
      continue;

    rhs = Expr(s, symbols);

    if (interactive)
      cout << "Trying to solve x(n) = " << rhs << endl;

    Expr solution;
    if (solve_try_hard(rhs, n, solution) != OK) {
      if (interactive)
	cout << "Sorry, this is too difficult." << endl;
    }
    else if (interactive) {
      cout << "*** SOLUTION ***"
	   << endl
	   << solution
	   << endl
	   << "****************"
	   << endl << endl;
    }
  }
  return 0;
} 
catch (exception &p) {
  cerr << "std::exception caught: " << p.what() << endl;
  return 1;
}
catch (const char* s) {
  cerr << s << endl;
  return 1;
}
