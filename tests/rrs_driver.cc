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
static bool expect_upper_bound;
static bool expect_lower_bound;
static bool expect_provably_correct_result;
static bool expect_provably_wrong_result;
static bool expect_inconclusive_verification;
static bool expect_not_to_be_solved;
static bool expect_diagnose_unsolvable;
static bool expect_not_diagnose_unsolvable;

void
set_expectations(const string& s) {
  // No expectations by default.
  explodes
    = expect_exactly_solved
    = expect_upper_bound
    = expect_lower_bound
    = expect_provably_correct_result
    = expect_provably_wrong_result
    = expect_inconclusive_verification
    = expect_not_to_be_solved
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
    case 'v':
      expect_provably_correct_result = true;
      break;
    case 'w':
      expect_provably_wrong_result = true;
      break;
    case 'd':
      expect_inconclusive_verification = true;
      break;
    case 'u':
      expect_upper_bound = true;
      break;
    case 'l':
      expect_lower_bound = true;
      break;
    case 'n':
      expect_not_to_be_solved = true;
      break;
    case 'K':
      expect_diagnose_unsolvable = true;
      break;
    case 'k':
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
compute_exact_solution_wrapper(const Recurrence& rec) {
  try {
    return rec.compute_exact_solution();
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

ostream&
operator<<(ostream& s, Recurrence::Verify_Status v) {
  switch (v) {
  case Recurrence::PROVABLY_CORRECT:
    s << "PROVABLY_CORRECT";
    break;
  case Recurrence::PROVABLY_INCORRECT:
    s << "PROVABLY_INCORRECT";
    break;
  case Recurrence::INCONCLUSIVE_VERIFICATION:
    s << "INCONCLUSIVE_VERIFICATION";
    break;
  }
  return s;
}

int
main(int argc, char *argv[]) try {
  program_name = argv[0];

  set_handlers();

  //purrs_initialize();

  process_options(argc, argv);

  unsigned unexpected_solution_or_bounds_for_it = 0;
  unsigned unexpected_exact_failures = 0;
  unsigned unexpected_lower_failures = 0;
  unsigned unexpected_upper_failures = 0;
  unsigned unexpected_unsolvability_diagnoses = 0;
  unsigned unexpected_failures_to_diagnose_unsolvability = 0;
  unsigned unexpected_failures_to_verify = 0;
  unsigned unexpected_failures_to_disprove = 0;
  unsigned unexpected_conclusive_verifications = 0;

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

  Symbol a("a");
  Symbol b("b");
  Symbol c("c");
  Symbol d("d");
  Expr_List symbols(Recurrence::n, a, b, c, d);

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
      rhs = Expr(rhs_string, symbols);
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

    Expr lower;
    Expr upper;
    Expr exact_solution;
    // *** regress test
    if (regress_test) {
      if (expect_exactly_solved) {
	if (compute_exact_solution_wrapper(rec) == Recurrence::SUCCESS) {
	  // Get the exact solution.
	  rec.exact_solution(exact_solution);
	  if (expect_provably_correct_result
	      || expect_provably_wrong_result
	      || expect_inconclusive_verification) {
	    Recurrence::Verify_Status status = rec.verify_solution();
	    if (expect_provably_correct_result
		&& status != Recurrence::PROVABLY_CORRECT) {
	      if (verbose)
		cerr << "*** unexpected failure to verify solution: gave "
		     << status << endl;
	      ++unexpected_failures_to_verify;
	    }
	    if (expect_provably_wrong_result
		&& status != Recurrence::PROVABLY_INCORRECT) {
	      if (verbose)
		cerr << "*** unexpected failure to disprove solution: gave "
		     << status << endl;
	      ++unexpected_failures_to_disprove;
	    }
	    if (expect_inconclusive_verification
		&& status != Recurrence::INCONCLUSIVE_VERIFICATION) {
	      if (verbose)
		cerr << "*** unexpected conclusive verification: gave "
		     << status << endl;
	      ++unexpected_conclusive_verifications;
	    }
	  }
	  if (interactive) {
	    cout << "*** EXACT SOLUTION ***"
		 << endl
		 << exact_solution
		 << endl
		 << "****************"
		 << endl << endl;
#if 0
	    cout << "*** APPROXIMATED SOLUTION ***"
		 << endl
		 << rec.approximated_solution()
		 << endl
		 << "*****************************"
		 << endl << endl;
#endif
	  }
	  if (latex) {
	    cout << "\\medskip\\indent\n\\(\n x(n) = ";
	    exact_solution.latex_print(cout);
	  }
	}
	else
	  // There is not the exact solution.
	  if (verbose) {
	    cerr << "*** unexpected failure to find exact solution" << endl;
	    ++unexpected_exact_failures;
	  }
      }
      if (expect_lower_bound)
	if (rec.compute_lower_bound() == Recurrence::SUCCESS) {
	  // Get the lower bound.
	  rec.lower_bound(lower);
	  if (interactive)
	    cout << "*** LOWER BOUND ***"
		 << endl
		 << lower
		 << endl
		 << "*******************"
		 << endl << endl;
	  if (latex) {
	    cout << "\\medskip\\indent\n\\(\n x(n) = ";
	    lower.latex_print(cout);
	  }
	}
	else {
	  // There is not the lower bound.
	  if (verbose)
	    cerr << "*** unexpected failure to find lower bound for "
		 << "the solution" << endl;
	  ++unexpected_lower_failures;
	}
      if (expect_upper_bound)
	if (rec.compute_upper_bound() == Recurrence::SUCCESS) {
	  // Get the upper bound.
	  rec.upper_bound(upper);
	  if (interactive)
	    cout << "*** UPPER BOUND ***"
		 << endl
		 << upper
		 << endl
		 << "*******************"
		 << endl << endl;
	  if (latex) {
	    cout << "\\medskip\\indent\n\\(\n x(n) = ";
	    upper.latex_print(cout);
	  }
	}
	else {
	  // There is not the upper bound.
	  if (verbose)
	    cerr << "*** unexpected failure to find upper bound for "
		 << "the solution" << endl;
	  ++unexpected_upper_failures;
	}
      if (expect_not_to_be_solved)
	if (compute_exact_solution_wrapper(rec) == Recurrence::SUCCESS
	    || rec.compute_lower_bound() == Recurrence::SUCCESS
	    || rec.compute_upper_bound() == Recurrence::SUCCESS) {
	  if (verbose)
	    cerr << "*** unexpected solution or bounds for it" << endl;
	  ++unexpected_solution_or_bounds_for_it;
	}
      if (expect_not_diagnose_unsolvable)
	if (compute_exact_solution_wrapper(rec)
	    == Recurrence::UNSOLVABLE_RECURRENCE
	    || rec.compute_lower_bound() == Recurrence::UNSOLVABLE_RECURRENCE
	    || rec.compute_upper_bound()
	    == Recurrence::UNSOLVABLE_RECURRENCE) {
	  if (verbose)
	    cerr << "*** unexpected unsolvability diagnosis" << endl;
	  ++unexpected_unsolvability_diagnoses;
	  if (interactive)
	    cout << "Unsolvable." << endl;
	}
      if (expect_diagnose_unsolvable)
	if (compute_exact_solution_wrapper(rec) == Recurrence::TOO_COMPLEX
	    || rec.compute_lower_bound() == Recurrence::TOO_COMPLEX
	    || rec.compute_upper_bound() == Recurrence::TOO_COMPLEX) {
	  if (verbose)
	    cerr << "*** unexpected failure to diagnose unsolvability" << endl;
	  ++unexpected_failures_to_diagnose_unsolvability;
	  if (interactive)
	    cout << "Sorry, this is too difficult." << endl;
	}
    } // *** regress test
    else {
      bool solved_or_approximated = false;
      // Try to solve the recurrence exactly.
      if (compute_exact_solution_wrapper(rec) == Recurrence::SUCCESS) {
	// OK: get the exact solution and print it.
	rec.exact_solution(exact_solution);
	solved_or_approximated = true;
	if (interactive) {
	  cout << "*** EXACT SOLUTION ***"
	       << endl
	       << exact_solution
	       << endl
	       << "****************"
	       << endl << endl;
#if 0
	  cout << "*** APPROXIMATED SOLUTION ***"
	       << endl
	       << rec.approximated_solution()
	       << endl
	       << "*****************************"
	       << endl << endl;
#endif
	}
      }
      else {
	// Cannot solve it exacly: try to get lower and upper bound.
	if (rec.compute_lower_bound() == Recurrence::SUCCESS) {
	  // OK: get the lower bound and print it.
	  rec.lower_bound(lower);
	  solved_or_approximated = true;
	  if (interactive)
	    cout << "*** LOWER BOUND ***"
		 << endl
		 << lower
		 << endl
		 << "*******************"
		 << endl << endl;
	}
	if (rec.compute_upper_bound() == Recurrence::SUCCESS) {
	  // OK: get the upper bound and print it.
	  rec.upper_bound(upper);
	  solved_or_approximated = true;
	  if (interactive)
	    cout << "*** UPPER BOUND ***"
		 << endl
		 << upper
		 << endl
		 << "*******************"
		 << endl << endl;
	}
      }
      bool unsolvable = false;
      if (!solved_or_approximated)
	if (compute_exact_solution_wrapper(rec)
	    == Recurrence::UNSOLVABLE_RECURRENCE
	    || rec.compute_lower_bound()
	    == Recurrence::UNSOLVABLE_RECURRENCE
	    || rec.compute_upper_bound()
	    == Recurrence::UNSOLVABLE_RECURRENCE) {
	  unsolvable = true;
	  if (interactive)
	    cout << "Unsolvable." << endl << endl;
	}
      if (!solved_or_approximated && !unsolvable) {
	assert(compute_exact_solution_wrapper(rec)
	       == Recurrence::TOO_COMPLEX
	       || rec.compute_lower_bound()
	       == Recurrence::TOO_COMPLEX
	       || rec.compute_upper_bound()
	       == Recurrence::TOO_COMPLEX);
	if (interactive)
	  cout << "Sorry, this is too difficult." << endl << endl;
      }
    }
    if (interactive)
      rec.dump(cout);
  } // while (input_stream)

  if (latex)
    cout << "\\end{document}" << endl;

  if (regress_test) {
    bool failed = false;
    if (unexpected_failures_to_verify > 0) {
      failed = true;
      cerr << unexpected_failures_to_verify
	   << " unexpected failures to verify solutions"
	   << endl;
    }
    if (unexpected_failures_to_disprove > 0) {
      failed = true;
      cerr << unexpected_failures_to_disprove
	   << "unexpected failures to disprove "
	   << endl;
    }
    if (unexpected_conclusive_verifications > 0) {
      failed = true;
      cerr << unexpected_conclusive_verifications
	   << "unexpected conclusive verifications "
	   << endl;
    }
    if (unexpected_solution_or_bounds_for_it > 0) {
      failed = true;
      cerr << unexpected_solution_or_bounds_for_it
	   << " unexpected solution or bounds for it"
	   << endl;
    }
    if (unexpected_exact_failures > 0) {
      failed = true;
      cerr << unexpected_exact_failures
	   << " unexpected failures to find exact solutions"
	   << endl;
    }
    if (unexpected_lower_failures > 0) {
      failed = true;
      cerr << unexpected_lower_failures
	   << " unexpected failures to find lower bound for solutions"
	   << endl;
    }
    if (unexpected_upper_failures > 0) {
      failed = true;
      cerr << unexpected_upper_failures
	   << " unexpected failures to find upper bound for solutions"
	   << endl;
    }
    if (unexpected_unsolvability_diagnoses > 0) {
      failed = true;
      cerr << unexpected_unsolvability_diagnoses
	   << " unexpected unsolvability diagnoses"
	   << endl;
    }
    if (unexpected_failures_to_diagnose_unsolvability > 0) {
      failed = true;
      cerr << unexpected_failures_to_diagnose_unsolvability
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
