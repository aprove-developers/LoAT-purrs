/* Simple program to test the recurrence relation solver.
   Copyright (C) 2005 Roberto Bagnara <bagnara@cs.unipr.it>

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
#include "timings.hh"

#include <ginac/ginac.h>

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

#define DEBUG 1

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

static struct option long_options[] = {
  {"recurrence",        required_argument, 0, 'R'},
  {"initial-condition", required_argument, 0, 'I'},
  {"help",              no_argument,       0, 'h'},
  {"regress-test",      no_argument,       0, 'r'},
  {"format",            required_argument, 0, 'f'},
  {"version",           no_argument,       0, 'V'},
  {0, 0, 0, 0}
};

enum Formats {TEXT, PROLOG};

const char* program_name = 0;

void
print_usage() {
  cerr << "Usage: " << program_name << " [OPTION]...\n\n"
    "  -R, --recurrence  \"<rec>\"        set the right-hand side of the recurrence\n"
    "                           that has to be solved/approximated\n" 
        "  -h, --help               print this help text\n"
        "  -f, --format=FORMAT      give output in the specified FORMAT:\n"
       "                           text (default), Prolog\n"
    "  -r, --regress-test       set regression-testing mode on\n"
    "  -V, --version            show version number and exit"
       << endl;
}

#define OPTION_LETTERS "I:R:f:hr::vV"

// To avoid mixing incompatible options.
static bool production_mode = false;
static bool test_mode = false;

// When true, the recurrence for the production mode has been specified.
static bool have_recurrence = false;

// When true, the output has the form of prolog term.
static bool prolog_term_required = false;

// Output format, defaults to text.
static int output = TEXT;

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

// If 'n' happens to be among the arguments, it will be temporary 
// replaced by a new symbol.
// Global scope is needed because these symbols will be used troughout
// the program.
Symbol n_replacement;
Expr real_var_symbol;

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

static void
do_not_mix_modes() {
  if (production_mode && test_mode) {
    cerr << program_name
         << ": production mode options (-R, -I, -E, -L, -U, -C, -P, -T) and\n"
         << "test mode options (-r, -v) are mutually exclusive"
         << endl;
    my_exit(1);
  }
}

static void
error_message(const char* msg) {
  cerr << program_name << ": " << msg << endl;
  my_exit(1);
}

static void
invalid_recurrence(const char* culprit) {
  cerr << program_name << ": invalid recurrence `" << culprit << "'" << endl;
  my_exit(1);
}

static bool
invalid_initial_condition(const Expr& e) {
  Number value;
  if (e.is_a_number(value) && !value.is_rational())
    return true;
  if (e == Recurrence::n)
    return true;
  else if (e.is_a_power())
    return invalid_initial_condition(e.arg(0))
      || invalid_initial_condition(e.arg(1));
  else if (e.is_a_function()) {
    if (e.is_the_x1_function() || e.is_the_x2_function())
      return true;
    for (unsigned i = e.nops(); i-- > 0; )
      if (invalid_initial_condition(e.arg(i)))
        return true;
  }
  else if (e.is_a_add() || e.is_a_mul())
    for (unsigned i = e.nops(); i-- > 0; )
      if (invalid_initial_condition(e.op(i)))
        return true;
  return false;
}

static void
invalid_initial_condition(const char* culprit) {
  cerr << program_name << ": invalid initial condition `" << culprit << "';\n"
       << "must be in the form `x(i, {par_1, par_2, ..., par_m}) == k'" << endl
       << "where `i' is an integer index, "
       << "`{par_1, par_2, ..., par_m}' a parameter list" << endl
       << "and `k' a rational number "
       << "or a symbolic expression" << endl
       << "containing the parameters a, b,..., z"
       << " different from `n', `x' and `e'." << endl;
  my_exit(1);
}

static Recurrence* precp = 0;
static Expr prhs = 0;
std::vector<Expr> recs;
std::vector<Expr> conds;


static void
init_production_recurrence() {
  if (precp == 0)
    precp = new Recurrence();
}

  std::map<index_type, Expr> function_args;


static void
process_options(int argc, char* argv[]) {
  int option_index;
  std::map<index_type, Expr> initial_conditions;

  while (true) {
    option_index = 0;
    int c = getopt_long(argc, argv, OPTION_LETTERS, long_options,
			&option_index);
    if (c == EOF)
      break;

    switch (c) {
    case 0:
      break;

    case 'R':
      {
        production_mode = true;
        do_not_mix_modes();
        if (!optarg) {
          cerr << program_name << ": a recurrence must follow `-R'"
               << endl;
          my_exit(1);
        }
        // Here `optarg' points to the beginning of the rhs.
        assert(optarg);
        Expr rec;
	// FIXME: Pre-parse as string to allow using `=' instead of `=='.
        if (!parse_expression(optarg, rec))
          invalid_recurrence(optarg);
#if 0
        have_recurrence = true;
        init_production_recurrence();
        precp->replace_recurrence(rec);
	// FIXME: Save lhs as well.
#endif
	recs.push_back(rec);
      }
      break;
      
    case 'I':
      {
        string cond = optarg;
	/*
        string::size_type eq_pos = cond.find('=');
        if (eq_pos == string::npos)
          invalid_initial_condition(optarg);
        string lhs(cond, 0, eq_pos);
        string rhs(cond, eq_pos+1);
        Expr l;
        if (!parse_expression(lhs, l))
          invalid_initial_condition(optarg);
        if (!l.is_the_x_function())
          invalid_initial_condition(optarg);
        // FIXME: lhs == rhs returns a boolean. Do not evaluate it.
	conds.push_back(lhs == rhs);
	*/
	Expr c;
        if (!parse_expression(optarg, c) || !c.is_a_relational())
          invalid_initial_condition(optarg);
	Expr lhs = c.op(0);
	if (!lhs.is_the_x2_function())
          invalid_initial_condition(optarg);
	Expr rhs = c.op(1);
	// FIXME: Make sure lhs and rhs are well-formed.
	conds.push_back(c);
	//	cerr << conds[conds.size()-1] << " - " << lhs << " - " << rhs << endl;
        /*
        production_mode = true;
        do_not_mix_modes();
        // Here `optarg' points to the beginning of the initial condition.
        assert(optarg);
        Number index;
        if (!l.arg(0).is_a_number(index)
            || !index.is_nonnegative_integer()
            || index > LONG_MAX)
          invalid_initial_condition(optarg);
        Expr r;
        if (!parse_expression(rhs, r))
          invalid_initial_condition(optarg);
        if (invalid_initial_condition(r))
          invalid_initial_condition(optarg);
        init_production_recurrence();
	// Insert the pair (index, r), which represents the initial
	// condition `x(index) = r', in the map `initial_conditions'.
	initial_conditions.insert(std::map<index_type, Expr>
				  ::value_type(index.to_unsigned_int(), r));
	*/
      }
      break;

    case '?':
    case 'h':
      print_usage();
      my_exit(0);
      break;

    case 'f':
      {
	string output_format = optarg;
	if (output_format=="text" || output_format=="t")
	  output = TEXT;
	else if (output_format=="prolog" || output_format=="p") {
	  output = PROLOG;
	  prolog_term_required = true;
	}
	else {
          cerr << program_name << ": Invalid output format '" << output_format << "'."
               << endl;
          my_exit(1);
        }
      }
      break;

    case 'r':
      {
	cerr << program_name << ": Regression testing is not implemented yet." << endl;
	my_exit(1);
	/*
	if (optarg) {
	  char* endptr;
	  tries = strtol(optarg, &endptr, 10);
	};
	if (tries <= 0) {
          cerr << program_name << ": invalid optional tries number in regression testing mode."
               << endl;
          my_exit(1);
	};
      regress_test = true;
      test_mode = true;
      do_not_mix_modes();
	*/
      }
      break;

    case 'v':
      verbose = true;
      // Verbose output can be enabled in both operation modes.
      // test_mode = true;
      // do_not_mix_modes();
      break;

    case 'V':
      cerr << "multivar_driver version " << PACKAGE_VERSION << ".\n"
        "This is the Parma Recurrence Relation Solver driver.\n"
        "Copyright (C) 2002-2005 Roberto Bagnara <bagnara@cs.unipr.it>\n"
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

static bool explodes;
static bool expect_exactly_solved;
static bool expect_upper_bound;
static bool expect_lower_bound;
static bool expect_provably_correct_result;
static bool expect_partial_provably_correct_result;
static bool expect_provably_wrong_result;
static bool expect_inconclusive_verification;
static bool expect_not_to_be_solved;
static bool expect_diagnose_unsolvable;
static bool expect_not_diagnose_unsolvable;
static bool expect_diagnose_indeterminate;
static bool expect_diagnose_malformed;
static bool expect_diagnose_domain_error;

void
set_expectations(const string& s) {
  // No expectations by default.
  explodes
    = expect_exactly_solved
    = expect_upper_bound
    = expect_lower_bound
    = expect_provably_correct_result
    = expect_partial_provably_correct_result
    = expect_provably_wrong_result
    = expect_inconclusive_verification
    = expect_not_to_be_solved
    = expect_diagnose_unsolvable
    = expect_not_diagnose_unsolvable
    = expect_diagnose_indeterminate
    = expect_diagnose_malformed
    = expect_diagnose_domain_error
    = false;

  const char* p = s.c_str();
  while (char c = *p++) {
    if (isspace(c))
      return;
    switch (c) {
    case 'E':
    case 'b':
      explodes = true;
      break;
    case 'y':
      expect_exactly_solved = true;
      break;
    case 'v':
      expect_provably_correct_result = true;
      break;
    case 'p':
      expect_partial_provably_correct_result = true;
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
    case 'B':
    case 'X':
    case 'n':
      expect_not_to_be_solved = true;
      break;
    case 'K':
      expect_diagnose_unsolvable = true;
      break;
    case 'k':
      expect_not_diagnose_unsolvable = true;
      break;
    case 'I':
      expect_diagnose_indeterminate = true;
      break;
    case 'M':
      expect_diagnose_malformed = true;
      break;
    case 'D':
      expect_diagnose_domain_error = true;
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
operator<<(ostream& s, const Recurrence::Verify_Status& v) {
  switch (v) {
  case Recurrence::PROVABLY_CORRECT:
    s << "PROVABLY_CORRECT";
    break;
  case Recurrence::PARTIAL_PROVABLY_CORRECT:
    s << "PARTIAL_PROVABLY_CORRECT";
    break;
  case Recurrence::PROVABLY_INCORRECT:
    s << "PROVABLY_INCORRECT";
    break;
  case Recurrence::INCONCLUSIVE_VERIFICATION:
    s << "INCONCLUSIVE_VERIFICATION";
    break;
  default:
    throw std::runtime_error("PURRS internal error: "
			     "operator<<(ostream&, Verify_Status&).");
  }
  return s;
}

string
set_string_validation(const Recurrence::Verify_Status& v) {
  switch (v) {
  case Recurrence::PROVABLY_CORRECT:
    return "provably_correct";
  case Recurrence::PARTIAL_PROVABLY_CORRECT:
    return "partial_provably_correct";
  case Recurrence::PROVABLY_INCORRECT:
    return "provably_incorrect";
  case Recurrence::INCONCLUSIVE_VERIFICATION:
    return "inconclusive_verification";
  default:
    throw std::runtime_error("PURRS internal error: "
			     "set_string_validation().");
  }
}

//! Kinds of solution or approximation required.
enum Kind { EXACT, LOWER, UPPER };

// Print the solution and all informations about it.
void
output_solution(const Kind& kind, const Expr& solution_or_bound,
		const std::vector<string>& conditions,
		const std::map<index_type, Expr>& initial_conditions,
		const std::vector<string>& blackboard,
		const string& validation) {
  std::ostringstream s;
  switch (kind) {
  case EXACT:
    s << "x(n) = ";
    break;
  case LOWER:
    s << "x(n) >= ";
    break;
  case UPPER:
    s << "x(n) <= ";
    break;
  }
  // Solution or bound.
  s << solution_or_bound << ", ";
  cout << s.str() << endl;
  // List of conditions.
  if (!conditions.empty()) {
    cout << "for ";
    unsigned int conditions_size = conditions.size();
    for (unsigned int i = 0; i < conditions_size; ++i) {
      cout << conditions[i];
      if (i != conditions_size - 1)
	cout << ", ";
    }
    cout << endl;
  }
  // Blackboard.
  if (!blackboard.empty()) {
    cout << "where ";
    unsigned int blackboard_size = blackboard.size();
    for (unsigned int i = 0; i < blackboard_size; ++i) {
      cout << blackboard[i];
      if (i != blackboard_size - 1)
	cout << ", ";
    }
    cout << endl;
  }
  // List of initial conditions.
  if ((kind == LOWER || kind == UPPER) && initial_conditions.size() != 1) {
    for (std::map<index_type, Expr>::const_iterator i
	   = initial_conditions.begin(),
	   initial_conditions_end = initial_conditions.end(), j = i;
	 i != initial_conditions_end; ++i) {
      cout << "x(" << i->first << ")=" << i->second;
      if (++j != initial_conditions_end)
	cout << ", ";
    }
    if (!initial_conditions.empty())
      cout << endl;
  }
  // Result of the validation's process.
  if (!validation.empty())
    cout << validation << endl;
}

// Print the solution and all informations about it in the form of
// prolog term.
void
output_solution_prolog_term(const Kind& kind, const Expr& solution_or_bound,
			    const std::vector<string>& conditions,
			    const std::map<index_type, Expr>&
			    initial_conditions,
			    const std::vector<string>& blackboard,
			    const string& validation) {
  std::ostringstream s;
  // Beginning of the Prolog term.
  switch (kind) {
  case EXACT:
    s << "exact_solution(";
    break;
  case LOWER:
    s << "lower_bound(";
    break;
  case UPPER:
    s << "upper_bound(";
    break;
  }
  // Solution or bound.
  s << solution_or_bound << ", ";
  // List of conditions.
  s << "[";
  unsigned int conditions_size = conditions.size();
  for (unsigned int i = 0; i < conditions_size; ++i) {
    s << conditions[i];
    if (i != conditions_size - 1)
      s << ", ";
  }
  s << "], ";
  // List of initial conditions.
  s << "[";
  for (std::map<index_type, Expr>::const_iterator i
	 = initial_conditions.begin(),
	 initial_conditions_end = initial_conditions.end(), j = i;
       i != initial_conditions_end; ++i) {
    s << "x(" << i->first << ")=" << i->second;
    if (++j != initial_conditions_end)
      s << ", ";
  }
  s << "], ";
  // Blackboard;
  s << "[";
  unsigned int blackboard_size = blackboard.size();
  for (unsigned int i = 0; i < blackboard_size; ++i) {
    s << blackboard[i];
    if (i != blackboard_size - 1)
      s << ", ";
  }
  s << "]";
  // Result of the validation's process.
  if (!validation.empty())
    s << ", " << validation;
  // End of the Prolog term.
  s << ").";
  
  cout << s.str() << endl;
}

void
prepare_for_the_output(const Kind& kind,
		       std::vector<string>& conditions,
		       std::map<index_type, Expr>& initial_conditions,
		       std::vector<string>& blackboard) {
  std::ostringstream s;
  s << "n>=" << precp->first_valid_index_for_solution();
  conditions.push_back(s.str());
  
  // If there is also the exact solution then this is the case
  // of "trivial" bound: the necessary informations on the bound are
  // the informations on the exact solution.
  if (!compute_exact_solution_wrapper(*precp) == Recurrence::SUCCESS
      && (kind == LOWER || kind == UPPER)) {
    string Sc_function = precp->definition_Sc();
    std::ostringstream cond_i_c;
    if (initial_conditions.empty()) {
      if (Sc_function.empty())
	cond_i_c << "x(1)>=0";
      else {
	cond_i_c << "x(" << Sc_function.substr(0, 8) << ")>=0";
      }
      conditions.push_back(cond_i_c.str());
    }
    if (!Sc_function.empty())
      blackboard.push_back(Sc_function);
  }
}

void
result_of_the_verification(unsigned type,
                           const Recurrence::Verify_Status& status,
                           vector<unsigned>& unexpected_failures_to_verify,
                           vector<unsigned>& unexpected_failures_to_partially_verify,
                           vector<unsigned>& unexpected_failures_to_disprove,
                           vector<unsigned>& unexpected_conclusive_verifications) {
  if (expect_provably_correct_result
      && status != Recurrence::PROVABLY_CORRECT) {
    if (verbose)
      if (type == 0)
        cerr << "*** unexpected failure to verify solution: gave "
             << status << endl;
      else if (type == 1)
        cerr << "*** unexpected failure to verify lower bound: gave "
             << status << endl;
      else
        cerr << "*** unexpected failure to verify upper bound: gave "
             << status << endl;
    unexpected_failures_to_verify.push_back(line_number);
  }
  if (expect_partial_provably_correct_result
      && status != Recurrence::PARTIAL_PROVABLY_CORRECT) {
    if (verbose)
      if (type == 0)
        cerr << "*** unexpected failure to partially verify solution: gave "
             << status << endl;
      else if (type == 1)
        cerr << "*** unexpected failure to partially verify lower bound: gave "
             << status << endl;
      else
        cerr << "*** unexpected failure to partially verify upper bound: gave "
             << status << endl;
    unexpected_failures_to_partially_verify.push_back(line_number);
  }
  if (expect_provably_wrong_result
      && status != Recurrence::PROVABLY_INCORRECT) {
    if (verbose)
      if (type == 0)
        cerr << "*** unexpected failure to disprove solution: gave "
             << status << endl;
      else if (type == 1)
        cerr << "*** unexpected failure to disprove lower bound: gave "
             << status << endl;
      else
        cerr << "*** unexpected failure to disprove upper bound: gave "
             << status << endl;
    unexpected_failures_to_disprove.push_back(line_number);
  }
  if (expect_inconclusive_verification
      && status != Recurrence::INCONCLUSIVE_VERIFICATION) {
    if (verbose)
      if (type == 0)
        cerr << "*** unexpected conclusive solution's verification: gave "
             << status << endl;
      else if (type == 1)
        cerr << "*** unexpected conclusive lower bound's verification: gave "
             << status << endl;
      else
        cerr << "*** unexpected conclusive upper bound's verification: gave "
             << status << endl;
    unexpected_conclusive_verifications.push_back(line_number);
  }  
}


bool find_terms_with_x(const Expr& this_term, vector<Expr>& terms_with_x) {
  if (this_term.is_the_x2_function()) {
    // FIXME: Check that it is well-formed.
    // FIXME: Check arity.
    terms_with_x.push_back(this_term);
  }
  else if (this_term.is_the_x1_function()) {
    error_message("The x() function must have two arguments: a unique number and a list of arguments.");
  }
  else if (this_term.nops() == 0)
    return true;
  else { 	   
    for (int i = this_term.nops()-1; i >= 0; --i) {
      // FIXME: why do we use different names for the op() and arg() functions?
      // Write a generalized op().
      Expr inner_term=(this_term.is_a_function() || this_term.is_a_power()?this_term.arg(i):this_term.op(i));
      if (!find_terms_with_x(inner_term, terms_with_x))
	return false;
    }
    return true;
  }
  // Useless. Just to avoid GCC warnings.
  return true;
}

// Insert initial conditions.
bool insert_initial_conditions(Expr& solution) {
  bool replaced = false;
  for (unsigned int i = 0; i < solution.nops(); ++i) {
    const Expr& this_term = solution.op(i);
    if (this_term.is_the_x2_function()) {
      for (unsigned int j = 0; j < conds.size(); ++j) {
	bool matches = true;
	const Expr& term_list = this_term.arg(1);
	const Expr& this_lhs = conds[j].op(0);
	const Expr& this_list = this_lhs.arg(1);
	const Expr& this_rhs = conds[j].op(1);
	if (this_term.arg(0) == this_lhs.arg(0)) {
	  for (unsigned int k = 0; k < this_list.nops(); ++k) {
	    if (this_list.op(k).is_a_number()) {
	      if (!term_list.op(k).is_a_number())
		matches = false;
	      else
		if (term_list.op(k) != this_list.op(k))
		  matches = false;
	    }
	    if (!matches)
	      break;
	  }
	}
	else
	  matches = false;
	if (matches) {
	  Expr new_rhs = this_rhs;
	  for (unsigned int h = 0; h < term_list.nops(); ++h)
	    new_rhs = new_rhs.substitute(this_list.op(h), term_list.op(h));
	  solution = solution.substitute(this_term, new_rhs);
	  replaced = true;
	}
      }
    }
  }
  return replaced;
}

      
Recurrence::Solver_Status multivar_solve(Expr& lhs, Expr& rhs, vector<Expr>& terms_with_x, 
		    const Symbol& n_replacement, Expr& real_var_symbol, Expr& solution) {
  const int num_param = lhs.arg(1).nops();
  vector<int> dummy(num_param);

  // We need a symbol different from any used symbol.
  Symbol different_symbol;

  lhs = lhs.substitute(Recurrence::n, n_replacement);
  rhs = rhs.substitute(Recurrence::n, n_replacement);
  for (unsigned i = terms_with_x.size(); i-- > 0; )
    terms_with_x[i] = terms_with_x[i].substitute(Recurrence::n, n_replacement);
  for (unsigned i = conds.size(); i-- > 0; )
    conds[i] = conds[i].substitute(Recurrence::n, n_replacement);
  
  for (unsigned i = num_param; i-- > 0; ) {
    dummy[i] = true;
    for (vector<Expr>::const_iterator j = terms_with_x.begin(); j != terms_with_x.end(); ++j) 
      if (lhs.arg(0) == j->arg(0)) {
	if (!dummy[i])
	  break;
	for (unsigned k = j->arg(1).nops(); k-- > 0; ) {
	  const Expr& examined_arg = lhs.arg(1).op(i);
	  const Expr& this_arg = j->arg(1).op(k);
	  // An argument can be neglected if it always occurs as itself or as a number...
	  if (k==i) {
	    // FIXME: Can we afford to have greater flexibility here?
	    if (!this_arg.is_a_number() && this_arg != examined_arg) {
	      dummy[i] = false;
	      break;
	    }
	  }
	  else {
	    // ... but it mustn't appear anywhere else as well.
	    if (this_arg != this_arg.substitute(examined_arg, different_symbol)) {
	      dummy[i] = false;
	      break;
	    }
	  }
	}
      }
  }
  std::vector<int> real_var_index;
  
  for (int i = num_param - 1; i >= 0; --i) {
    if (!dummy[i]) {
      if (real_var_index.size() >= 2) {
	error_message("Too many true parameters. Resolution can succeed only if\n"
		      "all the arguments of `x()' but one - or two if they are related - are useless.");
      }
      else
	real_var_index.push_back(i);
    }
  }
  
  // One of the x2() functions will be converted into the x1() function.
  index_type converted_x_index;

  Recurrence rec;
  Recurrence::Solver_Status outcome;
  bool constant_difference = false;
  bool constant_sum = false;
    
  if (real_var_index.size() == 1) {
    for (vector<Expr>::const_iterator i = terms_with_x.begin(); i != terms_with_x.end(); ++i) {
      const Expr& this_term = (*i);
      Expr real_var_expr = this_term.arg(1).op(real_var_index[0]);
      bool symbol_found = false;
      for (int j = real_var_expr.nops() - 1; j >=0; --j) {
	if (real_var_expr.op(j).is_a_symbol()) {
	  if (symbol_found) {
	    assert (real_var_symbol == real_var_expr.op(j));
	  }
	  else
	    real_var_symbol = real_var_expr.op(j);
	}
      }
      rhs = rhs.substitute(this_term, x(real_var_expr.substitute(real_var_symbol, Recurrence::n)));
      // FIXME: Behave properly when multiple x2() functions appear.
      converted_x_index = this_term.arg(0).ex_to_number().to_unsigned_int();
    }
    
    rec.replace_recurrence(rhs);
    
    outcome = compute_exact_solution_wrapper(rec);
    
    if (outcome == Recurrence::SUCCESS) {
      rec.exact_solution(solution);
      
      if (verbose)
	std::cerr << solution << endl;
      
      
      // FIXME: The recurrence must not have been rewritten for this to succeed.
      
      // Restore original arity and symbol names.
      for (unsigned int i = 0; i < solution.nops(); ++i) {
	const Expr& this_term = solution.op(i);
	if (this_term.is_the_x1_function()) {
	  solution = solution.substitute(this_term, 
					 lhs.substitute(real_var_symbol, this_term.arg(0)));
	}
      }

      insert_initial_conditions(solution);
      
      // Replace the substituted symbols back to their place.
      solution = solution.substitute(Recurrence::n, real_var_symbol);
      solution = solution.substitute(n_replacement, Recurrence::n);
      lhs = lhs.substitute(n_replacement, Recurrence::n);
      
    }
  }
  else if (real_var_index.size() == 2) {
    // We can solve the recurrence even if we have two true parameters
    // provided that their difference is constant.
    
    // Check whether the difference or sum is constant.
    constant_difference = true;
    bool invalid_constant_difference = false;
    constant_sum = true;
    bool invalid_constant_sum = false;
    size_t decreasing_variable;
    size_t increasing_variable;
    
    Expr real_var_expr_0;
    Expr real_var_expr_1;
    Expr real_var_symbol_0;
    Expr real_var_symbol_1;
    
    for (vector<Expr>::const_iterator i = terms_with_x.begin(); i != terms_with_x.end(); ++i) {
      const Expr& this_term = (*i);
      real_var_expr_0 = this_term.arg(1).op(real_var_index[0]);
      real_var_expr_1 = this_term.arg(1).op(real_var_index[1]);
      real_var_symbol_0 = lhs.arg(1).op(real_var_index[0]);
      real_var_symbol_1 = lhs.arg(1).op(real_var_index[1]);
      if (real_var_expr_0 - real_var_expr_1 != real_var_symbol_0 - real_var_symbol_1)
	constant_difference = false;
      if (real_var_expr_0 + real_var_expr_1 != real_var_symbol_0 + real_var_symbol_1)
	constant_sum = false;
      if (real_var_expr_0.substitute(real_var_symbol_0, different_symbol) !=  
	  real_var_expr_1.substitute(real_var_symbol_1, different_symbol))
	invalid_constant_difference = true;
      // FIXME: Write a similar check for invalid_constant_sum.
      // FIXME: Do the following check for each term.
      if (compare(real_var_expr_0, real_var_symbol_0) == -1) {
	decreasing_variable = real_var_index[0];
	increasing_variable = real_var_index[1];
      }
      else {
	decreasing_variable = real_var_index[1];
	increasing_variable = real_var_index[0];
      }
    }

    if (!constant_difference && !constant_sum)
      error_message("Too complex. Neither difference nor sum is constant");
    if ((constant_difference && invalid_constant_difference) ||
	(constant_sum && invalid_constant_sum)) {
	error_message("Too complex. Terms on the right hand side depend on too many terms on the left hand side.");
    }
    if (constant_difference) {
#if DEBUG
      cout << "Constant difference." << endl;
#endif
    
      // Pick one of the two variables and solve the recurrence with respect to it.
      for (vector<Expr>::const_iterator i = terms_with_x.begin(); i != terms_with_x.end(); ++i) {
	const Expr& this_term = (*i);
	const Expr& real_var_expr_0 = this_term.arg(1).op(real_var_index[0]);
	real_var_symbol_0 = lhs.arg(1).op(real_var_index[0]);
	rhs = rhs.substitute(this_term, x(real_var_expr_0.substitute(real_var_symbol_0, Recurrence::n)));
      }
      
      rec.replace_recurrence(rhs);
      
      outcome = compute_exact_solution_wrapper(rec);
      
      if (outcome == Recurrence::SUCCESS) {
	rec.exact_solution(solution);
	
	if (verbose)
	  cout << "Auxiliary recurrence solution: " << solution << endl;
	
	// FIXME: The recurrence must not have been rewritten for this to succeed.
	
	// The final solution will be given as a combination of two possibile solutions.
	Expr solution_0 = solution;
	Expr solution_1 = solution;
	// Restore original arity and symbol names.
	for (unsigned int i = 0; i < solution.nops(); ++i) {
	  const Expr& this_term = solution.op(i);
	  if (this_term.is_the_x1_function()) {
	    solution_0 = solution_0.substitute(this_term, 
					       lhs.substitute(real_var_symbol_0, this_term.arg(0)));
	    solution_1 = solution_1.substitute(this_term, 
					       lhs.substitute(real_var_symbol_0, this_term.arg(0)));
	    //	  cout << this_term << " - " << solution << endl;
	  }
	}
	
	insert_initial_conditions(solution_0);
	insert_initial_conditions(solution_1);

	// Replace the substituted symbols back to their place.
	solution_0 = solution_0.substitute(Recurrence::n, real_var_symbol_0);
	solution_1 = solution_1.substitute(Recurrence::n, real_var_symbol_1);
       
	Expr zero = 0;
	// FIXME: Use min if it is allowed.
	// FIXME: Handle the case m = n as well.
	solution = max(real_var_symbol_0 - real_var_symbol_1, zero) / 
	  (real_var_symbol_0 - real_var_symbol_1) * solution_0 +
	  max(real_var_symbol_1 - real_var_symbol_0, zero) / 
	  (real_var_symbol_1 - real_var_symbol_0) * solution_1;

	// If 'n` was substituted, perform the inverse substitution.
	solution = solution.substitute(n_replacement, Recurrence::n);
	lhs = lhs.substitute(n_replacement, Recurrence::n);

      }
    }
    else {
      if (constant_sum) {
#if DEBUG
	cout << "Constant sum." << endl;
#endif
    
	// Solve the recurrence with respect to the decreasing variable.
	for (vector<Expr>::const_iterator i = terms_with_x.begin(); i != terms_with_x.end(); ++i) {
	  const Expr& this_term = (*i);
	  const Expr& real_var_expr_0 = this_term.arg(1).op(decreasing_variable);
	  real_var_symbol_0 = lhs.arg(1).op(decreasing_variable);
	  rhs = rhs.substitute(this_term, x(real_var_expr_0.substitute(real_var_symbol_0, Recurrence::n)));
	}
	
	rec.replace_recurrence(rhs);
	
	outcome = compute_exact_solution_wrapper(rec);
	
	if (outcome == Recurrence::SUCCESS) {
	  rec.exact_solution(solution);
	  
	  if (verbose)
	    cout << "Auxiliary recurrence solution: " << solution << endl;
	  
	  // Restore original arity and symbol names.
	  real_var_symbol_1 = lhs.arg(1).op(increasing_variable);
	  for (unsigned int i = 0; i < solution.nops(); ++i) {
	    const Expr& this_term = solution.op(i);
	    if (this_term.is_the_x1_function()) {
	      solution = solution.substitute(this_term, 
					     lhs.substitute(real_var_symbol_0, this_term.arg(0))
					     .substitute(real_var_symbol_1, real_var_symbol_1 + real_var_symbol_0));
	    }
	  }
	  
	  insert_initial_conditions(solution);
	  
	  // Replace the substituted symbols back to their place.
	  solution = solution.substitute(Recurrence::n, real_var_symbol_0);
	  solution = solution.substitute(n_replacement, Recurrence::n);
	  lhs = lhs.substitute(n_replacement, Recurrence::n);
	}
      }
    }
  }
  else
    error_message("Internal error while finding true parameters");

  return outcome;
}

int
main(int argc, char *argv[]) try {
  program_name = argv[0];

  set_handlers();

  init_symbols();

  process_options(argc, argv);

  if (recs.size() > 1)
    error_message("Systems of recurrences are not yet supported");

  Expr this_rec = recs[0];

  // FIXME: Check that it is an equality and not an inequality.
  if (!this_rec.is_a_relational())
    error_message("Invalid equation: must be in the form `lhs == rhs'");

  Expr lhs = this_rec.op(0);
  Expr rhs = this_rec.op(1);

  if (!lhs.is_the_x2_function())
    error_message("Invalid lhs: must be in the form `x(index, {arg_list})'");

  const Expr& index_expr = lhs.arg(0);
  const Expr& param_list = lhs.arg(1);

  if (!index_expr.is_a_number())
    error_message("Invalid lhs: must be in the form `x(index, {arg_list})'");

  if (!param_list.is_a_Expr_List())
    error_message("Invalid lhs: must be in the form `x(index, {arg_list})'");

  unsigned int index = index_expr.ex_to_number().to_unsigned_int();

  function_args.insert(std::map<index_type, Expr>::value_type(index, param_list));

  std::vector<Expr> terms_with_x;

  if (!find_terms_with_x(rhs, terms_with_x)) {
    std::cerr << "Inconsistent number of arguments for the function x()." << std::endl
	      << std::endl;
    abort();
  }

  Expr exact_solution;
  
  Recurrence::Solver_Status solve;

  solve = multivar_solve(lhs, rhs, terms_with_x, n_replacement, real_var_symbol, exact_solution);
  
#ifdef DEBUG
  cout << "Auxiliary recurrence: " << rhs << endl;
#endif
  
  switch (solve) {
  case Recurrence::SUCCESS: {
    
    switch (output) {
      // FIXME: Prolog must be handled differently.
    case PROLOG:
    case TEXT:
      cout << lhs << " = " << exact_solution << endl;
      break;
    }
    goto finish;
    break;
  }
  case Recurrence::UNSOLVABLE_RECURRENCE:
    cout << endl << "Unsolvable" << endl << endl;
        goto finish;
        break;
  case Recurrence::INDETERMINATE_RECURRENCE:
    cout << endl << "Indeterminate" << endl << endl;
    goto finish;
    break;
  case Recurrence::MALFORMED_RECURRENCE:
    cout << endl << "Malformed" << endl << endl;
    goto finish;
    break;
  case Recurrence::DOMAIN_ERROR:
    cout << endl << "Domain error" << endl << endl;
    goto finish;
    break;
  case Recurrence::TOO_COMPLEX:
    cout << endl << "Too complex" << endl << endl;
  default:
    break;
  }
  
 finish:
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
