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
#include <vector>
#include <cassert>
#include <getopt.h>
#include <cctype>

#include "purrs_install.hh"
#include "ehandlers.hh"

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

#ifdef USE_READLINE
#define INPUT_STREAM rdl
#else
#define INPUT_STREAM cin
#endif

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
      interactive = true;
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

bool
all_space(const string& s) {
  for (unsigned i = s.length(); i-- > 0; )
    if (!isspace(s[i]))
      return false;
  return true;
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

  Symbol x("x");
  Symbol a("a");
  Symbol b("b");
  Symbol c("c");
  Symbol d("d");
  Expr_List symbols(Recurrence::n, a, b, c, d, x);

  while (input_stream) {
    string s;
    getline(input_stream, s);

    // We may be at end of file.
    if (!input_stream)
      return 0;

    // The string may be constituted by white space only.
    if (all_space(s))
      continue;

    // Skip comments.
    if (s.find("%") == 0)
      continue;

    Expr e = Expr(s, symbols);

    // `expand' does the simplification of the rule E6.
    Expr e_expanded = e.expand();
    if (interactive)
      cout << "Expanded expression = " << e_expanded << endl << endl;
    Expr solution_1 = simplify_on_input_ex(e_expanded, Recurrence::n, true);
    Expr solution_2 = simplify_on_output_ex(e_expanded, Recurrence::n, false);
    Expr solution_3 = simplify_factorials_and_exponentials(e, Recurrence::n);

    if (interactive) {
      cout << "Simplifications for the input (to collect the symbol `n')" << endl;
      cout << solution_1 << endl << endl;
      cout << "Simplifications for the output" << endl;
      cout << solution_2 << endl << endl;
      cout << "Factorials and exponentials' simplifications " << endl;
      cout << solution_3 << endl;
      cout << endl << "---------------------------------------------"
	   << endl << endl;
    }
  }
  return 0;
}
catch (exception &p) {
  cerr << "std::exception caught: " << p.what() << endl;
  my_exit(1);
}
catch (const char* s) {
  cerr << s << endl;
  my_exit(1);
}

