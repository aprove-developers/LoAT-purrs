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

#ifdef USE_READLINE
#define INPUT_STREAM rdl
#else
#define INPUT_STREAM cin
#endif

static istream* pinput_stream;

static struct option long_options[] = {
  {"help",           no_argument,          0, 'h'},
  {"interactive",    no_argument,          0, 'i'},
  {"regress-test",   no_argument,          0, 'r'},
  {"verbose",        no_argument,          0, 'v'},
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
static unsigned num_errors = 0;

static void
message(const string& s) {
  cerr << program_name << ": on line " << line_number << ": " << s << endl;
}

static void
error(const string& s) {
  message(s);
  my_exit(1);
}

static bool expect_right_simplification;
static bool simplification_wrong;

void
set_expectations(const string& s) {
  expect_right_simplification
    = simplification_wrong
    = false;

  const char* p = s.c_str();
  while (char c = *p++) {
    if (isspace(c))
      return;
    switch (c) {
    case 'n':
      simplification_wrong = true;
      break;
    case 'y':
      expect_right_simplification = true;
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
all_space(const string& s) {
  for (unsigned i = s.length(); i-- > 0; )
    if (!isspace(s[i]))
      return false;
  return true;
}

Expr
trick_ginac(const Expr_List& symbols, const string& s) {
  // Dirty trick here: GiNaC's parser sucks!
  Expr dummy_expr("0", symbols);
  // We sandwich the string to be parsed within
  // "x(" and ")" to force its complete parsing.
  Expr tmp("x(" + s + ")", symbols);
  if (tmp == dummy_expr)
    throw std::runtime_error("not detected by GiNaC's parser");
  // Get the mortadella.
  return tmp.arg(0);
}

void
get_expression(Expr& e, const Expr_List& symbols) {
  istream& input_stream = *pinput_stream;
  for ( ; ; ) {
    ++line_number;

    string the_line;
    getline(input_stream, the_line);
 
    // We may be at end of file.
    if (!input_stream)
      error("premature end of file: simplified expression is missing");

    istringstream line(the_line);

    // Skip comments.
    if (the_line.find("%") == 0)
      continue;
    
    // The line may be constituted by white space only.
    if (all_space(the_line))
      continue;

    try {
      e = trick_ginac(symbols, the_line);
    }
    catch (exception& except) {
      ostringstream m;
      m << "parse error: " << except.what();
      error(m.str());
    }

    // Everything was OK.
    break;
  }
}

void
check_and_notify(const char* what,
		 const Expr& source, const Expr& expected, const Expr& computed) {
  if (expected != computed) {
    ++num_errors;
    if (verbose)
      cerr << "simplification for " << what << " of " << source
	   << endl
	   << "was expected to be " << expected
	   << endl
	   << "but resulted in " << computed
	   << endl
	   << "(before line " << line_number << ")"
	   <<endl;
  }
}

int
main(int argc, char *argv[]) try {
  program_name = argv[0];

  set_handlers();

  process_options(argc, argv);

  readlinebuf* prdlb = 0;

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
    string s;

    Expr e_simpl_input;
    Expr e_simpl_output;
    Expr e_simpl_factorials;      

    if (regress_test) {
      // Read the expectations' string and use it.
      string expectations;
      line >> expectations;

      // Skip empty lines.
      if (!line)
	continue;

      set_expectations(expectations);
      
      getline(line, s);
      // Premature end of file?
      if (!line)
	error("premature end of file after expectations' string");

      get_expression(e_simpl_input, symbols);
      get_expression(e_simpl_output, symbols);
      get_expression(e_simpl_factorials, symbols);
    }
    else
      getline(line, s);

    // The string may be constituted by white space only.
    if (all_space(s))
      continue;

    Expr ex;
    try {
      ex = trick_ginac(symbols, s);
    }
    catch (exception& e) {
      ostringstream m;
      m << "parse error: " << e.what();
      message(m.str());
      continue;
    }

    if (verbose) {
      if (!interactive)
	cout << line_number - 3 << ": ";
      cout << "e = " << ex << endl;
    }

    Expr ex_expanded = ex.expand();
    if (interactive)
      cout << "Expanded expression = " << ex_expanded << endl;
    
    Expr sol_input = simplify_on_input_ex(ex_expanded, Recurrence::n, true);
    Expr sol_output = simplify_on_output_ex(ex_expanded, Recurrence::n, false);
    Expr sol_factorials = simplify_factorials_and_exponentials(ex,
							       Recurrence::n);
    if (expect_right_simplification) {
      check_and_notify("input", ex, e_simpl_input, sol_input);
      check_and_notify("output", ex, e_simpl_output, sol_output);
      check_and_notify("factorials", ex, e_simpl_factorials, sol_factorials);
    }    
    if (interactive) {
      cout << endl;
      cout << "Simplifications for the input (to collect the symbol `n')"
	   << endl;
      cout << sol_input << endl << endl;
      cout << "Simplifications for the output" << endl;
      cout << sol_output << endl << endl;
      cout << "Factorials and exponentials' simplifications " << endl;
      cout << sol_factorials << endl;
      cout << endl << "---------------------------------------------"
	   << endl << endl;
    }
  }
  if (regress_test)
    if (num_errors > 0) {
      cerr << num_errors
	   << " differences from expected simplifications"
	   << endl;
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

