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
#include <sstream>
#include <cassert>
#include <cctype>
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
  {"help",              no_argument,       0, 'h'},
  {"interactive",       no_argument,       0, 'i'},
  {"regress-test",      no_argument,       0, 'r'},
  {"verbose",           no_argument,       0, 'v'},
  {0, 0, 0, 0}
};

const char* program_name = 0;

void
print_usage() {
  cerr << "Usage: " << program_name << " [OPTION]...\n\n"
    "  -h, --help              prints this help text\n"
    "  -i, --interactive       sets interactive mode on\n"
    "  -r, --regress-test      sets regression-testing mode on\n"
    "  -v, --verbose           be verbose"
       << endl;
}

#define OPTION_LETTERS "hirv"

// Interactive mode is on when true.
static bool interactive = false;

// Regression-testing mode is on when true.
static bool regress_test = false;

// Verbose mode is on when true.
static bool verbose = false;

static void
my_exit(int status) {
  //(void) purrs_finalize();
  exit(status);
}

static void
process_options(int argc, char* argv[]) {
  int option_index;
  int c;

  while (true) {
    option_index = 0;
    c = getopt_long(argc, argv, OPTION_LETTERS, long_options, &option_index);
    if (c == EOF)
      break;

    switch (c) {
    case 0:
      break;

    case '?':
    case 'h':
      print_usage();
      my_exit(0);
      break;

    case 'i':
      interactive = true;
      break;

    case 'r':
      regress_test = true;
      break;

    case 'v':
      verbose = true;
      break;

    default:
      abort();
    }
  }

  if (optind < argc) {
    print_usage();
    my_exit(1);
  }
}

static unsigned line_number = 0;

static void
message(const string& s) {
  cerr << program_name << ": on line " << line_number << ": " << s << endl;
}

static void
error(const string& s) {
  message(s);
  my_exit(1);
}

static bool expect_solved;
static bool expect_not_solved;

void
set_expectations(const string& s) {
  // No expectations by default.
  expect_solved
    = expect_not_solved
    = false;

  const char* p = s.c_str();
  while (char c = *p++) {
    if (isspace(c))
      return;
    switch (c) {
    case 'y':
      expect_solved = true;
      break;
    case 'n':
      expect_not_solved = true;
      break;
    case '*':
      break;
    default:
      {
	ostringstream m;
	m << "unsupported expectation `" << c << "' in `" << s << "'";
	error(m.str());
      }
    }
  }
}

bool
solve_wrapper(const Recurrence& rec, const Symbol& n) {
  try {
    return rec.solve(n);
  }
  catch (exception& e) {
    if (verbose) {
      ostringstream m;
      m << "std::exception caught: " << e.what();
      message(m.str());
    }
  }
  catch (const char* s) {
    if (verbose)
      message(s);
  }
  return false;
}

int
main(int argc, char *argv[]) try {
  program_name = argv[0];

  set_handlers();

  //purrs_initialize();

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

  while (input_stream) {
    ++line_number;

    string the_line;
    getline(input_stream, the_line);

    // We may be at end of file.
    if (!input_stream)
      my_exit(0);

    // Skip comments.
    if (the_line.find("%") == 0)
      continue;

    istringstream line(the_line);
    string rhs_string;

    if (regress_test) {
      // Read the expectations' string and use it.
      string expectations;
      line >> expectations;

      // Skip empty lines.
      if (!line)
	continue;

      set_expectations(expectations);

      getline(line, rhs_string);
      // Premature end of file?
      if (!line)
	error("premature end of file after expectations' string");
    }
    else
      getline(line, rhs_string);

    Expr rhs;
    try {
      rhs = Expr(rhs_string, symbols);
    }
    catch (exception& e) {
      ostringstream m;
      m << "parse error: " << e.what();
      message(m.str());
      continue;
    }

    if (verbose) {
      if (!interactive)
	cerr << line_number << ": ";
      cout << "trying to solve x(n) = " << rhs << endl;
    }

    Recurrence rec(rhs);

    Expr solution;
    if (!solve_wrapper(rec, n)) {
      if (interactive)
	cout << "Sorry, this is too difficult." << endl;
    }
    else if (interactive) {
      cout << "*** SOLUTION ***"
	   << endl
	   << rec.exact_solution(n)
	   << endl
	   << "****************"
	   << endl << endl;
    }
  }
  my_exit(0);
} 
catch (exception &p) {
  cerr << "std::exception caught: " << p.what() << endl;
  my_exit(1);
}
catch (const char* s) {
  cerr << s << endl;
  my_exit(1);
}
