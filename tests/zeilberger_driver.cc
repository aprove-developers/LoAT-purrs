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
#include <getopt.h>

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

static struct option long_options[] = {
  {"interactive",       no_argument,       0, 'i'},
  {"regress-test",      optional_argument, 0, 'r'},
  {"verbose",           no_argument,       0, 'v'},
  {"version",           no_argument,       0, 'V'},
  {0, 0, 0, 0}
};
                                                                                                                             
const char* program_name = 0;
                                                                                                                             
void
print_usage() {
  cerr << "Usage: " << program_name << " [OPTION]...\n\n"
    "  -r, --regress-test	set regression-testing mode on\n"
    "  -v, --verbose            be verbose\n"
    "  -V, --version            show version number and exit"
       << endl;
}
                                                                                                                             
#define OPTION_LETTERS "irvV"
                                                                                                                             
static bool test_mode = false;
                                                                                                                             
// Verbose mode is on when true.
static bool verbose = false;

static void
my_exit(int status) {
  exit(status);
}

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
                                                                                                                             
    case 'r':
      test_mode = true;
      break;
                                                                                                                             
    case 'v':
      verbose = true;
      break;
                                                                                                                             
    case 'V':
      cerr << "zeilberger_driver version " << PACKAGE_VERSION << ".\n"
        "This is the Parma Recurrence Relation Solver simple driver.\n"
        "Copyright (C) 2002-2004 Roberto Bagnara <bagnara@cs.unipr.it>\n"
        "This is free software; see the source for copying conditions.  There is NO\n"
        "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
        "See http://www.cs.unipr.it/purrs/ for more information.\n"
           << endl;
      my_exit(0);
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


// FIXME: Save strings. Put in util.cc.
string itos(const int num)
{
    ostringstream buffer;
    buffer << num;
    return buffer.str();
}
                                                                                                                             
string dump_vector(const vector<unsigned>& vec) {
  const string separator = ", ";
  string tmp="";
  unsigned int vec_size = vec.size();
  if (vec_size == 0)
    return "";
  else
    tmp = itos(vec[0]);
  unsigned int i = 1;
  while (i < vec_size) {
      tmp += separator + itos(vec[i]);
      ++i;
  }
  return tmp;
}

void
errors_summary(const vector<unsigned>& errors_vec, const string& message, bool& failed) {
  if (errors_vec.size() > 0) {
    failed = true;
      cerr << errors_vec.size() << " " << message
           << " at lines " << dump_vector(errors_vec)
           << endl;
  }
}

int main(int argc, char *argv[]) try {

  process_options(argc, argv);

  unsigned int line_number = 0;
  vector<unsigned int> verification_failures;

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
    string s2;
#if NOISY
    cout << endl << "Enter the hypergeometric term F(n, k) to be summed over k: " << endl;
    cout << "(Blank line to use the hard-coded hypergeometric term " << default_term << " )" << endl;
#endif
    getline(input_stream,s);
    ++line_number;

    // We may be at end of file.
    if (!input_stream)
      break;
      
    // Skip comments.
    if (s.find("%") == 0)
      continue;

    // If we are in test mode, the next line contains the expected solution.
    if (test_mode) {
#if NOISY
      cout << "Enter expected solution: " << endl;
#endif
      getline(input_stream, s2);
      ++line_number;
      if (!input_stream)
	break;
    }

    if (!test_mode && s == "") {
      s = default_term;
    }
    Expr_List l(Recurrence::n, k, a, b, c, d);
    Expr t_n;
    Expr expected_solution;
    try {
      t_n = Expr(s,l);
      cout << t_n << endl;
      if (test_mode) {
      expected_solution = Expr(s2,l);
      cout << expected_solution << endl;
      }

      if (zeilberger_algorithm(t_n, Recurrence::n, k, solution)) {
#if NOISY
	cout << endl << "The sum is: " << solution << endl;
	if (test_mode) {
	  if (solution == expected_solution)
	    cout << "Expected solution verified." << endl;
	  else {
	    cout << "*** Expected solution not verified: expected " 
		 << expected_solution 
		 << ", found "  << solution  << "." 
		 << endl;
	    verification_failures.push_back(line_number);
	  }
	}
      }
      else {
	cout << endl << "Error in Zeilberger Algorithm for term " << t_n
	     << endl;
      }
#endif
    }
    catch (exception& e) {
      std::cerr << "std::exception caught: " << e.what();
    }
    
    cout << endl << "---------------------------------------------" << endl;
      
  }

  if (test_mode) {
    bool failed = false;
    // FIXME: Keep track of the algorithm failures, not only in test mode.
    // errors_summary(zeilberger_failures, "unexpected failures in Zeilberger Algorithm", failed);
    errors_summary(verification_failures, "unexpected failures to verify given solution", failed);
    if (failed)
      cerr << "*** Errors found during regression testing for Zeilberger Algorithm." << endl;
}

  return 0;
}
catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
  return 1;
}
