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
  {"latex",             no_argument,       0, 'l'},
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
    "  -l, --latex             output LaTeX code\n"
    "  -r, --regress-test      sets regression-testing mode on\n"
    "  -v, --verbose           be verbose"
       << endl;
}

#define OPTION_LETTERS "hilrv"

// Interactive mode is on when true.
static bool interactive = false;

// LaTeX mode is on when true.
static bool latex = false;

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

    case 'l':
      latex = true;
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

static bool explodes;
static bool expect_exactly_solved;
static bool expect_not_exactly_solved;
static bool expect_diagnose_unsolvable;
static bool expect_not_diagnose_unsolvable;

void
set_expectations(const string& s) {
  // No expectations by default.
  explodes
    = expect_exactly_solved
    = expect_not_exactly_solved
    = expect_diagnose_unsolvable
    = expect_not_diagnose_unsolvable
    = false;

  const char* p = s.c_str();
  while (char c = *p++) {
    if (isspace(c))
      return;
    switch (c) {
    case 'E':
      explodes = true;
      break;
    case 'y':
      expect_exactly_solved = true;
      break;
    case 'n':
      expect_not_exactly_solved = true;
      break;
    case 'U':
      expect_diagnose_unsolvable = true;
      break;
    case 'u':
      expect_not_diagnose_unsolvable = true;
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

Recurrence::Solver_Status
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
  return Recurrence::TOO_COMPLEX;
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
  program_name = argv[0];

  set_handlers();

  //purrs_initialize();

  process_options(argc, argv);

  unsigned unexpected_exact_solutions = 0;
  unsigned unexpected_exact_failures = 0;
  unsigned unexpected_unsolvability_diagnoses = 0;
  unsigned unexpected_failures_do_diagnose_unsolvability = 0;

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

  if (latex)
    cout << "\\documentclass{article}\n\\begin{document}" << endl;

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
      break;

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

      // Avoid exploding.
      if (explodes)
	continue;

      getline(line, rhs_string);
      // Premature end of file?
      if (!line)
	error("premature end of file after expectations' string");
    }
    else
      getline(line, rhs_string);

    // The string may be constituted by white space only.
    if (all_space(rhs_string))
      continue;

    Expr rhs;
    try {
      // Dirty trick here: GiNaC's parser sucks!
      Expr dummy_expr("0", symbols);
      // We sandwich the string to be parsed within
      // "x(" and ")" to force its complete parsing.
      Expr tmp("x(" + rhs_string + ")", symbols);
      if (tmp == dummy_expr)
	throw std::runtime_error("not detected by GiNaC's parser");
      // Get the mortadella.
      rhs = tmp.op(0);
    }
    catch (exception& e) {
      ostringstream m;
      m << "parse error: " << e.what();
      message(m.str());
      continue;
    }

    if (verbose) {
      if (latex)
	cout << "\\bigskip\n\n\\noindent\n\\textbf{Line";
      if (!interactive)
	cout << line_number << ": ";
      if (latex)
	cout << "} $";
      cout << "x(n) = ";
      if (latex) {
	rhs.latex_print(cout);
	cout << "$";
      }
      else
	cout << rhs;
      cout << endl;
      if (latex)
	cout << endl;
    }

    Recurrence rec(rhs);

    Expr solution;
    switch (solve_wrapper(rec, n)) {
    case Recurrence::OK:
      if (regress_test) {
	if (expect_not_exactly_solved) {
	  if (verbose)
	    cerr << "*** unexpected exact solution" << endl;
	  ++unexpected_exact_solutions;
	}
      }
      if (interactive) {
	cout << "*** SOLUTION ***"
	     << endl
	     << rec.exact_solution(n)
	     << endl
	     << "****************"
	     << endl << endl;
      }
      if (latex) {
	cout << "\\medskip\\indent\n\\(\n x(n) = ";
	rec.exact_solution(n).latex_print(cout);
	cout << "\n\\)\n" << endl;
      }
      break;

    case Recurrence::UNSOLVABLE_RECURRENCE:
      if (expect_not_diagnose_unsolvable) {
	if (verbose)
	  cerr << "*** unexpected unsolvability diagnosis" << endl;
	++unexpected_unsolvability_diagnoses;
      }
      if (interactive)
	cout << "Unsolvable." << endl;
      goto failed;
      break;

 
    default:
      if (expect_diagnose_unsolvable) {
	if (verbose)
	  cerr << "*** unexpected failure to diagnose unsolvability" << endl;
	++unexpected_failures_do_diagnose_unsolvability;
      }
      if (interactive)
	cout << "Sorry, this is too difficult." << endl;

    failed:
      if (regress_test) {
	if (expect_exactly_solved) {
	  if (verbose)
	    cerr << "*** unexpected failure to find exact solution" << endl;
	  ++unexpected_exact_failures;
	}
      }
    }
  } // while (input_stream)

  if (latex)
    cout << "\\end{document}" << endl;

  if (regress_test) {
    bool failed = false;
    if (unexpected_exact_solutions > 0) {
      failed = true;
      cerr << unexpected_exact_solutions
	   << " unexpected exact solutions"
	   << endl;
    }
    if (unexpected_exact_failures > 0) {
      failed = true;
      cerr << unexpected_exact_failures
	   << " unexpected failures to find exact solutions"
	   << endl;
    }
    if (unexpected_unsolvability_diagnoses > 0) {
      failed = true;
      cerr << unexpected_unsolvability_diagnoses
	   << " unexpected unsolvability diagnoses"
	   << endl;
    }
    if (unexpected_failures_do_diagnose_unsolvability > 0) {
      failed = true;
      cerr << unexpected_failures_do_diagnose_unsolvability
	   << " unexpected failures to diagnose unsolvability"
	   << endl;
    }

    if (failed)
      my_exit(1);
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
