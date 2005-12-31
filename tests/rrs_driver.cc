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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
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

#ifdef USE_READLINE
#include "readlinebuf.hh"
#include <memory>
#endif

#define PROFILE_VERIFICATION 0

using namespace std;
using namespace Parma_Recurrence_Relation_Solver;

static struct option long_options[] = {
  {"rhs",               required_argument, 0, 'R'},
  {"initial-condition", required_argument, 0, 'I'},
  {"exact",             no_argument,       0, 'E'},
  {"lower-bound",       no_argument,       0, 'L'},
  {"upper-bound",       no_argument,       0, 'U'},
  {"check-solution",    no_argument,       0, 'C'},
  {"prolog-term",       no_argument,       0, 'P'},
  {"timeout",           required_argument, 0, 'T'},
  {"help",              no_argument,       0, 'h'},
  {"interactive",       no_argument,       0, 'i'},
  // FIXME: Option -l is deprecated.
  {"latex",             no_argument,       0, 'l'},
  // FIXME: Option -m is deprecated.
  {"mathml",            no_argument,       0, 'm'},
  {"output",            required_argument, 0, 'o'},
  {"regress-test",      optional_argument, 0, 'r'},
  {"verbose",           no_argument,       0, 'v'},
  {"version",           no_argument,       0, 'V'},
  {0, 0, 0, 0}
};

enum Formats {TEXT, LATEX, MATHML};

const char* program_name = 0;

void
print_usage() {
  cerr << "Usage: " << program_name << " [OPTION]...\n\n"
    "  -R, --rhs \"<expr>\"       set the right-hand side of the recurrence\n"
    "                           that has to be solved/approximated\n" 
    "  -I, --initial-condition \"<expr>\"\n"
    "                           set an initial condition for the recurrence\n"
    "  -E, --exact              try to solve the recurrence\n"
    "  -L, --lower-bound        try to approximate the solution of the\n"
    "                           recurrence from below\n"
    "  -U, --upper-bound        try to approximate the solution of the\n"
    "                           recurrence from above\n"
    "  -C, --check-solution     try to verify the solution/approximation\n"
    "  -P, --prolog-term        print the solution/approximation as a prolog term\n"
    "  -T, --timeout N          interrupt computation after N seconds\n"
    "  -h, --help               print this help text\n"
    "  -i, --interactive        set interactive mode on\n"
    "  -o, --output=FORMAT      give output in the specified FORMAT:\n"
    "                           text (default), latex, mathml\n"
    // Deprecated: accept them but hide them from help screen.
    //    "  -l, --latex              output LaTeX code\n"
    //    "  -m, --mathml             output MathML code\n"
    "  -r, --regress-test [N]   set regression-testing mode on (N tries)\n"
    "  -v, --verbose            be verbose\n"
    "  -V, --version            show version number and exit"
       << endl;
}

#define OPTION_LETTERS "CEI:LPR:T:Uhilmo:r::vV"

// To avoid mixing incompatible options.
static bool production_mode = false;
static bool test_mode = false;

// When true, the recurrence for the production mode has been specified.
static bool have_recurrence = false;

// When true, the exact solution is required.
static bool exact_solution_required = false;

// When true, an approximation from below of the solution is required.
static bool lower_bound_required = false;

// When true, an approximation from above of the solution is required.
static bool upper_bound_required = false;

// When true, the validation of the solution or the approximation is required.
static bool validation_solution_required = false;

// When true, the output has the form of prolog term.
static bool prolog_term_required = false;

// When greater than zero, gives the timeout threshold.
static long timeout_threshold = 0;

// Interactive mode is on when true.
static bool interactive = false;

// Output format, defaults to text.
static int output = TEXT;

// Regression-testing mode is on when true.
static bool regress_test = false;

// Number of tries for each recurrence. Useful to diagnose weird bugs.
static long tries = 1;

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
         << "test mode options (-i, -r, -v) are mutually exclusive"
         << endl;
    my_exit(1);
  }
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
    if (e.is_the_x_function())
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
       << "must be of the form `x(i)=k'"
       << "for `i' a positive integer\n"
       << "and `k' a number not floating point "
       << "or a symbolic expression\n"
       << "containing the parameters a, b,..., z"
       << " different from `n', `x' and `e'." << endl;
  my_exit(1);
}

static Recurrence* precp = 0;

static void
init_production_recurrence() {
  if (precp == 0)
    precp = new Recurrence();
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

    case 'R':
      {
        production_mode = true;
        do_not_mix_modes();
        if (!optarg) {
          cerr << program_name << ": an expression must follow `-R'"
               << endl;
          my_exit(1);
        }
        // Here `optarg' points to the beginning of the rhs.
        assert(optarg);
        Expr rec;
        if (!parse_expression(optarg, rec))
          invalid_recurrence(optarg);
        have_recurrence = true;
        init_production_recurrence();
        precp->replace_recurrence(rec);
      }
      break;
      
    case 'I':
      {
        production_mode = true;
        do_not_mix_modes();
        // Here `optarg' points to the beginning of the initial condition.
        assert(optarg);
        string cond = optarg;
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
      }
      break;

    case 'E':
      production_mode = true;
      do_not_mix_modes();
      exact_solution_required = true;
      break;

    case 'L':
      production_mode = true;
      do_not_mix_modes();
      lower_bound_required = true;
      break;

    case 'U':
      production_mode = true;
      do_not_mix_modes();
      upper_bound_required = true;
      break;

    case 'C':
      production_mode = true;
      do_not_mix_modes();
      validation_solution_required = true;
      break;

    case 'P':
      production_mode = true;
      do_not_mix_modes();
      prolog_term_required = true;
      break;

    case 'T':
      {
        char* endptr;
        long secs = strtol(optarg, &endptr, 10);
        if (*endptr || secs < 0) {
          cerr << program_name << ": a non-negative integer must follow `-T'"
               << endl;
          my_exit(1);
        }
        else
          timeout_threshold = secs;
      }
      break;

    case '?':
    case 'h':
      print_usage();
      my_exit(0);
      break;

    case 'i':
      interactive = true;
      test_mode = true;
      do_not_mix_modes();
      break;

    // Deprecated: use --output=latex instead.
    case 'l':
      cerr << program_name 
	   << ": WARNING: Option -l is deprecated. Use --output=latex instead." << endl;
      output = LATEX;
      break;

    // Deprecated: use --output=mathml instead.
    case 'm':
      cerr << program_name 
	   << ": WARNING: Option -m is deprecated. Use --output=mathml instead." << endl;
      output = MATHML;
      break;

    case 'o':
      {
	string output_format = optarg;
	if (output_format=="text")
	  output = TEXT;
	else if (output_format=="latex")
	  output = LATEX;
	else if (output_format=="mathml")
	  output= MATHML;
	else {
          cerr << program_name << ": Invalid output format '" << output_format << "'."
               << endl;
          my_exit(1);
        }
      }
      break;

    case 'r':
      {
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
      }
      break;

    case 'v':
      verbose = true;
      test_mode = true;
      do_not_mix_modes();
      break;

    case 'V':
      cerr << "rrs_driver version " << PACKAGE_VERSION << ".\n"
        "This is the Parma Recurrence Relation Solver simple driver.\n"
        "Copyright (C) 2002-2003 Roberto Bagnara <bagnara@cs.unipr.it>\n"
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

void 
errors_summary(const vector<unsigned>& errors_vec, const string& message, bool& failed) {
  if (errors_vec.size() > 0) {
    failed = true;
      cerr << errors_vec.size() << " " << message
	   << " at lines " << dump_vector(errors_vec)
           << endl;
  }
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

#if PROFILE_VERIFICATION
class Profiler {
private:
  time_unit_t accumulator;
  std::string name;

public:
  Profiler(const std::string& n)
    : name(n) {
  }

  void accumulate(time_unit_t i) {
    accumulator += i;
  }

  ~Profiler() {
    std::cerr << "Time spent in " << name
	      << ": " << time_unit_to_usecs(accumulator)/1000 << " ms"
	      << std::endl;
  }
};

Profiler ves_profiler("verify_exact_solution()");
#endif

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
do_production_mode() {
  if (!have_recurrence) {
    cerr << program_name
         << ": must specify a recurrence using the `-R' option" << endl;
    my_exit(1);
  }

  // If nothing has been explicitely requested, we take it
  // as an implicit request for the exact solution.
  if (!exact_solution_required
      && !lower_bound_required
      && !upper_bound_required)
    exact_solution_required = true;

  if (exact_solution_required) {
    switch (compute_exact_solution_wrapper(*precp)) {
    case Recurrence::SUCCESS:
      {
	// OK: get the exact solution and print it.
	Expr exact_solution;
	precp->exact_solution(exact_solution);
	exact_solution
	  = precp->substitute_auxiliary_definitions(exact_solution);
	string validation;
	// Try to verify the solution.
	if (validation_solution_required) {
#if PROFILE_VERIFICATION
	  time_unit_t time_begin = get_time();
#endif
	  Recurrence::Verify_Status status = precp->verify_exact_solution();
	  validation = set_string_validation(status);
#if PROFILE_VERIFICATION
	  time_unit_t time_end = get_time(); 
	  ves_profiler.accumulate(time_end-time_begin);
#endif
	}
	// Get all informations about this recurrence.
	std::vector<string> conditions;
	std::vector<string> blackboard;
	std::map<index_type, Expr> initial_conditions
	  = precp->get_initial_conditions();
	prepare_for_the_output(EXACT, conditions, initial_conditions,
			       blackboard);
	if (prolog_term_required)
	  output_solution_prolog_term(EXACT, exact_solution, conditions,
				      initial_conditions, blackboard,
				      validation);
	else
	  output_solution(EXACT, exact_solution, conditions,
			  initial_conditions, blackboard, validation);
	goto exit;
      }

    case Recurrence::UNSOLVABLE_RECURRENCE:
      cout << "unsolvable." << endl;
      goto exit;

    case Recurrence::INDETERMINATE_RECURRENCE:
      cout << "indeterminate." << endl;
      goto exit;

    case Recurrence::MALFORMED_RECURRENCE:
      cout << "malformed." << endl;
      goto exit;

    case Recurrence::DOMAIN_ERROR:
      cout << "domain error." << endl;
      goto exit;

    case Recurrence::TOO_COMPLEX:
    default:
      cout << "exact(too_complex)." << endl;
      break;
    }
  }

  if (lower_bound_required) {
    switch (precp->compute_lower_bound()) {
    case Recurrence::SUCCESS:
      {
	// OK: get the lower bound and print it.
	Expr lower;
	precp->lower_bound(lower);
	string validation;
	// Try to verify the solution.
	if (validation_solution_required) {
#if PROFILE_VERIFICATION
	  time_unit_t time_begin = get_time();
#endif
	  Recurrence::Verify_Status status = precp->verify_lower_bound();
	  validation = set_string_validation(status);
#if PROFILE_VERIFICATION
	  time_unit_t time_end = get_time(); 
	  ves_profiler.accumulate(time_end-time_begin);
#endif
	}
	// Get all informations about this recurrence.
	std::vector<string> conditions;
	std::map<index_type, Expr> initial_conditions
	  = precp->get_initial_conditions();
	std::vector<string> blackboard;
	prepare_for_the_output(LOWER, conditions, initial_conditions,
			       blackboard);
	// Output must have the form of a prolog term.
	if (prolog_term_required)
	  output_solution_prolog_term(LOWER, lower, conditions,
				      initial_conditions, blackboard,
				      validation);
	else
	  output_solution(LOWER, lower, conditions,
			  initial_conditions, blackboard, validation);
      }
      break;

    case Recurrence::UNSOLVABLE_RECURRENCE:
      cout << "unsolvable." << endl;
      goto exit;

    case Recurrence::INDETERMINATE_RECURRENCE:
      cout << "indeterminate." << endl;
      goto exit;

    case Recurrence::MALFORMED_RECURRENCE:
      cout << "malformed." << endl;
      goto exit;

    case Recurrence::DOMAIN_ERROR:
      cout << "domain error." << endl;
      goto exit;

    case Recurrence::TOO_COMPLEX:
    default:
      cout << "lower_bound(too_complex)." << endl;
      break;
    }
  }

  if (upper_bound_required) {
    switch (precp->compute_upper_bound()) {
    case Recurrence::SUCCESS:
      {
	// OK: get the upper bound and print it.
	Expr upper;
	precp->upper_bound(upper);
	string validation;
	// Try to verify the solution.
	if (validation_solution_required) {
#if PROFILE_VERIFICATION
	  time_unit_t time_begin = get_time();
#endif
	  Recurrence::Verify_Status status = precp->verify_upper_bound();
	  validation = set_string_validation(status);
#if PROFILE_VERIFICATION
	  time_unit_t time_end = get_time(); 
	  ves_profiler.accumulate(time_end-time_begin);
#endif
	}
	// Get all informations about this recurrence.
	std::vector<string> conditions;
	std::map<index_type, Expr> initial_conditions
	  = precp->get_initial_conditions();
	std::vector<string> blackboard;
	prepare_for_the_output(LOWER, conditions, initial_conditions,
			       blackboard);
	// Output must have the form of a prolog term.
	if (prolog_term_required)
	  output_solution_prolog_term(UPPER, upper, conditions,
				      initial_conditions, blackboard,
				      validation);
	else
	  output_solution(UPPER, upper, conditions,
			  initial_conditions, blackboard, validation);
      }
      break;
    case Recurrence::UNSOLVABLE_RECURRENCE:
      cout << "unsolvable." << endl;
      goto exit;

    case Recurrence::INDETERMINATE_RECURRENCE:
      cout << "indeterminate." << endl;
      goto exit;

    case Recurrence::MALFORMED_RECURRENCE:
      cout << "malformed." << endl;
      goto exit;

    case Recurrence::DOMAIN_ERROR:
      cout << "domain error." << endl;
      goto exit;

    case Recurrence::TOO_COMPLEX:
    default:
      cout << "upper_bound(too_complex)." << endl;
      break;
    }
  }
 exit:
  my_exit(0);
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


int
main(int argc, char *argv[]) try {
  program_name = argv[0];

  set_handlers();

  init_symbols();

  process_options(argc, argv);

  if (production_mode)
    do_production_mode();

  vector<unsigned int> unexpected_solution_or_bounds_for_it;
  vector<unsigned int> unexpected_exact_failures;
  vector<unsigned int> unexpected_lower_failures;
  vector<unsigned int> unexpected_upper_failures;
  vector<unsigned int> unexpected_unsolvability_diagnoses;
  vector<unsigned int> unexpected_failures_to_diagnose_unsolvability;
  vector<unsigned int> unexpected_failures_to_diagnose_indeterminably;
  vector<unsigned int> unexpected_failures_to_diagnose_malformation;
  vector<unsigned int> unexpected_failures_to_diagnose_domain_error;
  vector<unsigned int> unexpected_failures_to_verify;
  vector<unsigned int> unexpected_failures_to_partially_verify;
  vector<unsigned int> unexpected_failures_to_disprove;
  vector<unsigned int> unexpected_conclusive_verifications;

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

  if (output==LATEX)
    cout << "\\documentclass{article}\n\\begin{document}" << endl;
  else if (output==MATHML)
    cout << "<?xml version=\"1.0\"?>\n"
      "<?xml-stylesheet type=\"text/xsl\" href=\"mathml.xsl\"?>\n"
      "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
      "<head><title>MathML Output</title></head>\n"
      "<body>\n";

  long try_number = 0;
  string the_line;

  while (input_stream) {

    if (try_number == 0) { 
      ++line_number;
      getline(input_stream, the_line);
    };
      
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
    if (!parse_expression(rhs_string, rhs)) {
      message(parse_error_message);
      continue;
    }
    
    if (verbose) {
      if (output==LATEX)
	cout << "\\bigskip\n\n\\noindent\n\\textbf{Line";
      if (!interactive) {
	cout << line_number;
	if (tries == 1)
	  cout << ": ";
	else
	  cout << "/" << try_number << ": ";
      };
      if (output==LATEX)
	cout << "} $";
      cout << "x(n) = ";
      if (output==LATEX) {
	rhs.latex_print(cout);
	cout << "$";
      }
      else
	cout << rhs;
      cout << endl;
      if (output==LATEX)
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
	      || expect_partial_provably_correct_result
	      || expect_provably_wrong_result
	      || expect_inconclusive_verification) {
#if PROFILE_VERIFICATION
	    time_unit_t time_begin = get_time();
#endif
	    Recurrence::Verify_Status status = rec.verify_exact_solution();
#if PROFILE_VERIFICATION
	    time_unit_t time_end = get_time(); 
	    ves_profiler.accumulate(time_end-time_begin);
#endif
	    result_of_the_verification(0, status,
				       unexpected_failures_to_verify,
				       unexpected_failures_to_partially_verify,
				       unexpected_failures_to_disprove,
				       unexpected_conclusive_verifications);
	  }
	  if (interactive) {
	    if (output==MATHML)
	      cout << "<p>\n";
	    cout << "*** EXACT SOLUTION ***" << endl;
	    if (output==MATHML)
	      cout << "</p>\n";
	    switch (output) {
	    case TEXT:
	      cout << exact_solution;
	      break;
	    case MATHML:
	      exact_solution.mathml_output(cout);
	      break;
	    case LATEX:
	      cout << "\\medskip\\indent\n\n $$ x(n) = ";
	      exact_solution.latex_print(cout);
	      cout << " $$\n\n";
	      break;
	    }
	    cout << endl
		 << (output==MATHML?"<p>\n":"")
		 << "**********************"
		 << (output==MATHML?"</p>\n":"")
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
	  // FIXME: Double LaTeX output in non-interactive mode.
	  if (output==LATEX) {
	    cout << "\\medskip\\indent\n\\(\n x(n) = ";
	    exact_solution.latex_print(cout);
	  }
	}
	else
	  // There is not the exact solution.
	  if (verbose) {
	    cerr << "*** unexpected failure to find exact solution" << endl;
	    unexpected_exact_failures.push_back(line_number);
	  }
      }
      if (expect_lower_bound)
	if (rec.compute_lower_bound() == Recurrence::SUCCESS) {
	  // Get the lower bound.
	  rec.lower_bound(lower);
	  if (expect_provably_correct_result
	      || expect_partial_provably_correct_result
	      || expect_provably_wrong_result
	      || expect_inconclusive_verification) {
	    Recurrence::Verify_Status status = rec.verify_lower_bound();
	    result_of_the_verification(1, status,
				       unexpected_failures_to_verify,
				       unexpected_failures_to_partially_verify,
				       unexpected_failures_to_disprove,
				       unexpected_conclusive_verifications);
	  }
	  if (interactive)
	    cout << "*** LOWER BOUND ***"
		 << endl
		 << lower
		 << endl
		 << "*******************"
		 << endl << endl;
	  if (output==LATEX) {
	    cout << "\\medskip\\indent\n\\(\n x(n) = ";
	    lower.latex_print(cout);
	  }
	}
	else {
	  // There is not the lower bound.
	  if (verbose)
	    cerr << "*** unexpected failure to find lower bound for "
		 << "the solution" << endl;
	  unexpected_lower_failures.push_back(line_number);
	}
      if (expect_upper_bound)
	if (rec.compute_upper_bound() == Recurrence::SUCCESS) {
	  // Get the upper bound.
	  rec.upper_bound(upper);
	  if (expect_provably_correct_result
	      || expect_partial_provably_correct_result
	      || expect_provably_wrong_result
	      || expect_inconclusive_verification) {
	    Recurrence::Verify_Status status = rec.verify_upper_bound();
	    result_of_the_verification(2, status,
				       unexpected_failures_to_verify,
				       unexpected_failures_to_partially_verify,
				       unexpected_failures_to_disprove,
				       unexpected_conclusive_verifications);
	  }
	  if (interactive)
	    cout << "*** UPPER BOUND ***"
		 << endl
		 << upper
		 << endl
		 << "*******************"
		 << endl << endl;
	  if (output==LATEX) {
	    cout << "\\medskip\\indent\n\\(\n x(n) = ";
	    upper.latex_print(cout);
	  }
	}
	else {
	  // There is not the upper bound.
	  if (verbose)
	    cerr << "*** unexpected failure to find upper bound for "
		 << "the solution" << endl;
	  unexpected_upper_failures.push_back(line_number);
	}
      if (expect_not_to_be_solved)
	if (compute_exact_solution_wrapper(rec) == Recurrence::SUCCESS
	    || rec.compute_lower_bound() == Recurrence::SUCCESS
	    || rec.compute_upper_bound() == Recurrence::SUCCESS) {
	  if (verbose)
	    cerr << "*** unexpected solution or bounds for it" << endl;
	  unexpected_solution_or_bounds_for_it.push_back(line_number);
	}
      if (expect_not_diagnose_unsolvable)
	if (compute_exact_solution_wrapper(rec)
	    == Recurrence::UNSOLVABLE_RECURRENCE
	    || rec.compute_lower_bound()
	    == Recurrence::UNSOLVABLE_RECURRENCE
	    || rec.compute_upper_bound()
	    == Recurrence::UNSOLVABLE_RECURRENCE) {
	  if (verbose)
	    cerr << "*** unexpected unsolvability diagnosis" << endl;
	  unexpected_unsolvability_diagnoses.push_back(line_number);
	    if (interactive)
	      cout << "Unsolvable." << endl;
	  }
	if (expect_diagnose_unsolvable)
	  if (compute_exact_solution_wrapper(rec)
	      != Recurrence::UNSOLVABLE_RECURRENCE
	      || rec.compute_lower_bound()
	      != Recurrence::UNSOLVABLE_RECURRENCE
	      || rec.compute_upper_bound()
	      != Recurrence::UNSOLVABLE_RECURRENCE) {
	    if (verbose)
	      cerr << "*** unexpected failure to diagnose unsolvability" << endl;
	    unexpected_failures_to_diagnose_unsolvability.push_back(line_number);
	  }
	if (expect_diagnose_indeterminate)
	  if (compute_exact_solution_wrapper(rec)
	      != Recurrence::INDETERMINATE_RECURRENCE
	      || rec.compute_lower_bound()
	      != Recurrence::INDETERMINATE_RECURRENCE
	      || rec.compute_upper_bound()
	      != Recurrence::INDETERMINATE_RECURRENCE) {
	    if (verbose)
	      cerr << "*** unexpected failure to diagnose indeterminably"
		   << endl;
	    unexpected_failures_to_diagnose_indeterminably.push_back(line_number);
	  }
	if (expect_diagnose_malformed)
	  if (compute_exact_solution_wrapper(rec)
	      != Recurrence::MALFORMED_RECURRENCE
	      || rec.compute_lower_bound()
	      != Recurrence::MALFORMED_RECURRENCE
	      || rec.compute_upper_bound()
	      != Recurrence::MALFORMED_RECURRENCE) {
	    if (verbose)
	      cerr << "*** unexpected failure to diagnose malformation"
		   << endl;
	    unexpected_failures_to_diagnose_malformation.push_back(line_number);
	  }
	if (expect_diagnose_domain_error)
	  if (compute_exact_solution_wrapper(rec)
	      != Recurrence::DOMAIN_ERROR
	      || rec.compute_lower_bound()
	      != Recurrence::DOMAIN_ERROR
	      || rec.compute_upper_bound()
	      != Recurrence::DOMAIN_ERROR) {
	    if (verbose)
	      cerr << "*** unexpected failure to diagnose domain error"
		   << endl;
	    unexpected_failures_to_diagnose_domain_error.push_back(line_number);
	  }
      
      } // *** regress test
      else {
	switch (compute_exact_solution_wrapper(rec)) {
	case Recurrence::SUCCESS:
	  // OK: get the exact solution and print it.
	  rec.exact_solution(exact_solution);
	  if (interactive) {
	    if (output==MATHML)
	      cout << "<p>\n";
            cout << "*** EXACT SOLUTION ***" << endl;
	    if (output==MATHML)
	      cout << "</p>\n";
	    switch (output) {
	    case TEXT:
	      cout << exact_solution;
	      break;
	    case MATHML:
	      exact_solution.mathml_output(cout);
	      break;
	    case LATEX:
	      cout << "\\medskip\\indent\n\n $$ x(n) = ";
	      exact_solution.latex_print(cout);
	      cout << " $$\n\n";
	      break;
	    }
	    cout << endl
		 << (output==MATHML?"<p>\n":"")
		 << "**********************"
		 << (output==MATHML?"</p>\n":"")
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
        goto finish;
        break;
      case Recurrence::UNSOLVABLE_RECURRENCE:
        if (interactive)
          cout << endl << "Unsolvable" << endl << endl;
        goto finish;
        break;
      case Recurrence::INDETERMINATE_RECURRENCE:
        if (interactive)
          cout << endl << "Indeterminate" << endl << endl;
        goto finish;
        break;
      case Recurrence::MALFORMED_RECURRENCE:
        if (interactive)
          cout << endl << "Malformed" << endl << endl;
        goto finish;
        break;
      case Recurrence::DOMAIN_ERROR:
        if (interactive)
          cout << endl << "Domain error" << endl << endl;
        goto finish;
        break;
      case Recurrence::TOO_COMPLEX:
      default:
        break;
      }
      bool too_complex_printed = false;
      switch (rec.compute_lower_bound()) {
      case Recurrence::SUCCESS:
        // OK: get the lower bound and print it.
        rec.lower_bound(lower);
        if (interactive)
          cout << "*** LOWER BOUND ***"
               << endl
               << lower
               << endl
               << "*******************"
               << endl << endl;
        break;
      case Recurrence::UNSOLVABLE_RECURRENCE:
        if (interactive)
          cout << endl << "Unsolvable" << endl << endl;
        goto finish;
        break;
      case Recurrence::INDETERMINATE_RECURRENCE:
        if (interactive)
          cout << endl << "Indeterminate" << endl;
        goto finish;
        break;
      case Recurrence::MALFORMED_RECURRENCE:
        if (interactive)
          cout << endl << "Malformed" << endl;
        goto finish;
        break;
      case Recurrence::DOMAIN_ERROR:
        if (interactive)
          cout << endl << "Domain error" << endl;
        goto finish;
        break;
      case Recurrence::TOO_COMPLEX:
        if (!too_complex_printed && interactive) {
          cout << endl << "Too complex" << endl << endl;
          too_complex_printed = true;
        }
        break;
      default:
        break;
      }
      switch (rec.compute_upper_bound()) {
      case Recurrence::SUCCESS:
        // OK: get the upper bound and print it.
        rec.upper_bound(upper);
        if (interactive)
          cout << "*** UPPER BOUND ***"
               << endl
               << upper
               << endl
               << "*******************"
               << endl << endl;
        break;
      case Recurrence::UNSOLVABLE_RECURRENCE:
        if (interactive)
          cout << endl << "Unsolvable" << endl << endl;
        goto finish;
        break;
      case Recurrence::INDETERMINATE_RECURRENCE:
        if (interactive)
          cout << "Indeterminate" << endl;
        goto finish;
        break;
      case Recurrence::MALFORMED_RECURRENCE:
        if (interactive)
          cout << "Malformed" << endl;
        goto finish;
        break;
      case Recurrence::DOMAIN_ERROR:
        if (interactive)
          cout << "Domain error" << endl;
        goto finish;
        break;
      case Recurrence::TOO_COMPLEX:
        if (!too_complex_printed && interactive)
          cout << endl << "Too complex" << endl << endl;
        break;
      default:
        break;
      }
    }
  finish:
    if (interactive) {
      switch (output) {
      case MATHML:
	cout << "<p>\n";
	break;
      case LATEX:
	cout << "\\begin{verbatim}\n";
	break;
      }
      rec.dump(cout);
      switch (output) {
      case MATHML:
	cout << "</p>\n";
	break;
      case LATEX:
	cout << "\\end{verbatim}\n";
	break;
      }
    }
  
    try_number++;

    // Exit regression testing after <CODE>tries</CODE> attempts.
    if (try_number == tries) {
      try_number = 0;
      continue;
    }

  } // while (input_stream)
  
  if (output==LATEX)
    cout << "\\end{document}" << endl;
  else if (output==MATHML)
    cout << "</body></html>\n";

  if (regress_test) {
    bool failed = false;
    errors_summary(unexpected_failures_to_verify, "unexpected failures to verify solutions", failed);
    errors_summary(unexpected_failures_to_partially_verify, "unexpected failures to partially verify solutions", failed);
    errors_summary(unexpected_failures_to_disprove,  "unexpected failures to disprove", failed);
    errors_summary(unexpected_conclusive_verifications, "unexpected conclusive verifications", failed);
    errors_summary(unexpected_solution_or_bounds_for_it, "unexpected solutions or bounds for it", failed);
    errors_summary(unexpected_exact_failures, "unexpected failures to find exact solutions", failed);
    errors_summary(unexpected_lower_failures, "unexpected failures to find lower bound for solutions", failed);
    errors_summary(unexpected_upper_failures, "unexpected failures to find upper bound for solutions", failed);
    errors_summary(unexpected_unsolvability_diagnoses, "unexpected unsolvability diagnoses", failed);
    errors_summary(unexpected_failures_to_diagnose_unsolvability, "unexpected failures to diagnose unsolvability", failed);
    errors_summary(unexpected_failures_to_diagnose_indeterminably, "unexpected failures to diagnose indeterminably", failed);
    errors_summary(unexpected_failures_to_diagnose_malformation, "unexpected failures to diagnose malformation", failed);
    errors_summary(unexpected_failures_to_diagnose_domain_error, "unexpected failures to diagnose a domain error", failed);
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
