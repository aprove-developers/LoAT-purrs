/* Simple program to filter some recurrence relations from
   "heap" satisfying determinate properties.
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

// FIXME: capire
static struct option long_options[] = {
  {"help",                            no_argument,       0, 'h'},
  {"interactive",                     no_argument,       0, 'i'},
  {"regress-test",                    no_argument,       0, 'r'},
  {"verbose",                         no_argument,       0, 'v'},
  {"linear_finite_order_const_coeff", no_argument,       0, 'A'},
  {"linear_finite_order_var_coeff",   no_argument,       0, 'B'},
  {"non_linear_finite_order",         no_argument,       0, 'C'},
  {"weighted_average",                no_argument,       0, 'D'},
  {"functional_equation",             no_argument,       0, 'E'},
  {0, 0, 0, 0}
};

const char* program_name = 0;

void
print_usage() {
  cerr << "Usage: " << program_name << " [OPTION]...\n\n"
    "  -h, --help               print this help text\n"
    "  -i, --interactive        set interactive mode on\n"
    "  -r, --regress-test       set regression-testing mode on\n"
    "  -v, --verbose            be verbose\n"
    "  -A, --select A           select only linear finite order recurrences\n"
                                "\t\t\t\twith constant coefficients\n"
    "  -B, --select B           select only linear finite order recurrences\n"
                                "\t\t\t\twith variable coefficients\n"
    "  -C, --select C           select only non linear finite order recurrences\n"
    "  -D, --select D           select only weighted-average recurrences\n"
    "  -E, --select E           select only functional equations\n"
       << endl;
}

#define OPTION_LETTERS "hirvABCDE"

// Interactive mode is on when true.
static bool interactive = false;

// Linear finite order constant coefficients filter mode is on when true.
static bool linear_finite_order_const_coeff = false;

// Linear finite order variable coefficients filter mode is on when true.
static bool linear_finite_order_var_coeff = false;

// Non-linear finite order filter mode is on when true.
static bool non_linear_finite_order = false;

// Weighted-average filter mode is on when true.
static bool weighted_average = false;

// Functional equation filter mode is on when true.
static bool functional_equation = false;

// Regression-testing mode is on when true.
static bool regress_test = false;

// Verbose mode is on when true.
static bool verbose = false;

static void
my_exit(int status) {
  //(void) purrs_finalize();
  exit(status);
}

Expr_List symbols;

static void
init_symbols() {
  symbols.append(Recurrence::n);
  char a[2];
  a[1] = '\0';
  for (char c = 'a'; c <= 'z'; ++c)
    if (c != 'e' && c != 'n' && c != 'x') {
      a[0] = c;
      symbols.append(Symbol(a));
    }
}

static string parse_error_message;

static bool
parse_expression(const string& s, Expr& expr) {
  try {
    //  for (unsigned int i = 0; i < symbols.nops(); ++i)
    //    cout << symbols.op(i) << endl;
    expr = Expr(s, symbols);
    return true;
  }
  catch (exception& e) {
    ostringstream m;
    m << "parse error: " << e.what();
    parse_error_message = m.str();
    return false;
  }
}

static Recurrence* precp = 0;

static void
process_options(int argc, char* argv[]) {
  int option_index;
  int c;
  std::map<index_type, Expr> initial_conditions;

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

    case 'A':
      linear_finite_order_const_coeff = true;
      break;

    case 'B':
      linear_finite_order_var_coeff = true;
      break;

    case 'C':
      non_linear_finite_order = true;
      break;

    case 'D':
      weighted_average = true;
      break;

    case 'E':
      functional_equation = true;
      break;

    default:
      abort();
    }
  }

  if (optind < argc) {
    print_usage();
    my_exit(1);
  }

  if (!initial_conditions.empty())
    precp->set_initial_conditions(initial_conditions);
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
  
  init_symbols();
  
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

  while (input_stream) {
    ++line_number;
    
    string the_line;
    getline(input_stream, the_line);
    
    // We may be at end of file.
    if (!input_stream)
      break;
    
    // Skip comments.
    if (the_line.find("%") == 0) {
      cout << the_line << endl;
      continue;
    }

    istringstream line(the_line);
    string rhs_string;

    // Read the expectations' string and use it.
    string expectations;
    line >> expectations;
    if (regress_test) {
      // Skip empty lines.
      if (!line)
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
    if (!parse_expression(rhs_string, rhs)) {
      message(parse_error_message);
      continue;
    }

    Recurrence rec(rhs);
    Recurrence::Classifier_Status classifier_status
      = rec.classify_and_catch_special_cases();
    if (classifier_status == Recurrence::CL_SUCCESS)
      // *** regress test
      if (regress_test) {
	switch (rec.type_) {
	case Recurrence::ORDER_ZERO:
	case Recurrence::LINEAR_FINITE_ORDER_CONST_COEFF:
	  if (linear_finite_order_const_coeff)
 	    if (verbose)
 	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::LINEAR_FINITE_ORDER_VAR_COEFF:
	  if (linear_finite_order_var_coeff)
	    if (verbose)
	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::NON_LINEAR_FINITE_ORDER:
	  if (non_linear_finite_order)
	    if (verbose)
	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::WEIGHTED_AVERAGE:
	  if (weighted_average)
	    if (verbose)
 	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::FUNCTIONAL_EQUATION:
	  if (functional_equation)
	    if (verbose)
	      cout << expectations << "\t" << rhs << endl;
	  break;
	default:
	  throw std::runtime_error("PURRS internal error: rr_filter");
	}
      }
      else {
	switch (rec.type_) {
	case Recurrence::ORDER_ZERO:
	case Recurrence::LINEAR_FINITE_ORDER_CONST_COEFF:
	  if (linear_finite_order_const_coeff)
	    if (verbose)
	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::LINEAR_FINITE_ORDER_VAR_COEFF:
	  if (linear_finite_order_var_coeff)
	    if (verbose)
	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::NON_LINEAR_FINITE_ORDER:
	  if (non_linear_finite_order)
	    if (verbose)
	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::WEIGHTED_AVERAGE:
	  if (weighted_average)
	    if (verbose)
	      cout << expectations << "\t" << rhs << endl;
	  break;
	case Recurrence::FUNCTIONAL_EQUATION:
	  if (functional_equation)
	    if (verbose)

	  break;
	default:
	  throw std::runtime_error("PURRS internal error: rr_filter");
	}
      }
  }
}
catch (exception &p) {
  cerr << "std::exception caught: " << p.what() << endl;
  my_exit(1);
}
catch (const char* s) {
  cerr << s << endl;
  my_exit(1);
}
