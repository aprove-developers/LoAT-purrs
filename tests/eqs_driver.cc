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
#include <getopt.h>

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
  Expr_List symbols(x);
  Expr p;
  while (input_stream) {
    string s;
    getline(input_stream, s);

    if (!input_stream)
      return 0;

    // The string may be constituted by white space only.
    if (all_space(s))
      continue;

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

    if (interactive)
      cout << std::endl << "Trying to solve " << p << " = 0" << endl;

    std::vector<Polynomial_Root> roots;
    bool all_distinct;
    if (!find_roots(p, x, roots, all_distinct)) {
      if (interactive)
	cout << "Sorry, this is too difficult." << endl;
    }
    else {
      size_t n = roots.size();
      for (size_t i = 0; i < n; ++i) {
	Expr value = roots[i].value();
	Number multiplicity = roots[i].multiplicity();
	if (interactive) {
	  cout << "x_" << i+1 << " = " << value;
	  if (multiplicity > 1)
	    cout << " (multiplicity " << multiplicity << ")";
	  cout << endl;
	}
	if (!value.is_a_number())
	  if (interactive) {
	    cout << "****  x_" << i+1 << " ~= " << approximate(value) << endl;
	  }
      }
    }
  }
  return 0;
} catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
  return 1;
}
