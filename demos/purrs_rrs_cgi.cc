/* Simple CGI program to test the recurrence equation solver.
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

#include <config.h>

#include <csignal>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <new>
#include <vector>
#include <cassert>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include <cgicc/CgiDefs.h>
#include <cgicc/Cgicc.h>
#include <cgicc/HTTPHTMLHeader.h>
#include <cgicc/HTMLClasses.h>

#include "purrs_install.hh"
#include "timings.hh"

#if !defined(HAVE_GETRLIMIT) || !defined(HAVE_SETRLIMIT)
#error "We must have a way of limiting time and space!"
#endif

using std::string;
using std::cout;
using std::endl;
using std::vector;

using namespace Parma_Recurrence_Relation_Solver;
using namespace cgicc;

const int MAX_SECONDS_OF_CPU_TIME = 5;
const int MAX_VIRTUAL_MEMORY = 16*1024*1024;

static void
set_alarm_on_cpu_time(unsigned int seconds, void (*handler)(int)) {
  sigset_t mask;
  sigemptyset(&mask);

  struct sigaction s;
  s.sa_handler = handler;
  s.sa_mask = mask;
  s.sa_flags = SA_ONESHOT;

  int r;
  r = sigaction(SIGXCPU, &s, 0);
  if (r)
    throw("sigaction failed");

  struct rlimit t;
  r = getrlimit(RLIMIT_CPU, &t);
  if (r)
    throw("getrlimit failed");

  if (seconds < t.rlim_cur) {
    t.rlim_cur = seconds;
    r = setrlimit(RLIMIT_CPU, &t);
    if (r)
      throw("setrlimit failed");
  }
}

static void
limit_virtual_memory(unsigned int bytes) {
  struct rlimit t;
  int r = getrlimit(RLIMIT_AS, &t);
  if (r)
    throw("getrlimit failed");

  if (bytes < t.rlim_cur) {
    t.rlim_cur = bytes;
    r = setrlimit(RLIMIT_AS, &t);
    if (r)
      throw("setrlimit failed");
  }
}

#include <cerrno>
#include <cstring>

// Create a new Cgicc object containing all the CGI data.
static Cgicc cgi;

static void
footer() {
  // Try again!
  cout << p() << cgicc::div().set("align", "center");
  cout << a("Try again!").set("href", cgi.getEnvironment().getReferrer()) 
       << endl;

  cout << cgicc::div() << br() << hr().set("class", "half") << endl;
  cout << cgicc::div().set("align", "center") << endl;
  cout << a("PURRS").set("href", "http://www.cs.unipr.it/purrs/")
       << span(" recurrence relation solver", set("class", "red")) << br()
       << " by the "
       << a("PURRS development team")
    .set("href", "mailto:purrs-devel@cs.unipr.it") << "." << br() << br()
       << "A free service brought to you by " << br()
       << a(img()
            .set("src", "http://www.cs.unipr.it/images/cs_at_parma")
            .set("alt", "cs@parma")
            .set("border", "0"))
    .set("href", "http://www.cs.unipr.it/") << br()
       << endl;

  // End of document.
  cout << body() << html() << endl;
}

static void
error(const string& message) {
  // Reset all the HTML elements that might have been used
  // to their initial state, so that we get valid output.
  html::reset();        head::reset();                 body::reset();
  title::reset();       h1::reset();                 h4::reset();
  comment::reset();     td::reset();                 tr::reset(); 
  table::reset();       cgicc::div::reset();        p::reset(); 
  a::reset();           h2::reset();                 colgroup::reset();

  // Output the HTTP headers for an HTML document, and the HTML 4.0 DTD info.
  cout << HTTPHTMLHeader() << HTMLDoctype(HTMLDoctype::eStrict) << endl
       << html().set("lang", "en").set("dir", "ltr") << endl;

  // Set up the page's header and title.
  cout << head() << endl;

  // Output the style sheet portion of the header
  cout << style() << comment() << endl
       << "body { color: black; background-color: white; }" << endl
       << "hr.half { width: 60%; align: center; }" << endl
       << "span.red, STRONG.red { color: red; }" << endl
       << "div.notice { border: solid thin; padding: 1em; margin: 1em 0; "
       << "background: #ddd; }" << endl
       << comment() << style() << endl;

  cout << title("PURRS Demo Error") << endl
       << meta()
    .set("name", "author")
    .set("content", "PURRS development team") << endl
       << head() << endl;
    
  cout << body() << endl
       << h1() << "PURRS Demo " << span("Error", set("class", "red"))
       << h1() << endl
       << cgicc::div().set("align", "center").set("class", "notice") << endl
       << h2(message) << endl
       << cgicc::div() << endl;

  footer();
  exit(0);
}

void
my_timeout(int) {
  error("timeout");
}

typedef void (*PFV)(void);

static char* safety;
static PFV old_new_handler;

static void
my_new_handler () {
  delete safety;
  std::set_new_handler(old_new_handler);
}

void
my_unexpected_exception() {
  error("unexpected exception");
}

void
my_uncaught_exception() {
  error("uncaught exception");
}

#if 0
static void
my_exit(int status) {
  //(void) purrs_finalize();
  exit(status);
}
#endif

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

void
mark_verified_solution() {
  cout << a(img()
	    .set("src", "http://www.cs.unipr.it/purrs/images/verified")
	    .set("alt", "Verified solution")
	    .set("border", "0"))
    .set("href", "http://www.cs.unipr.it/purrs/details#verification")
       << " ";
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
  catch (std::exception& e) {
    std::ostringstream m;
    m << "parse error: " << e.what();
    parse_error_message = m.str();
    return false;
  }
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
invalid_initial_condition(const string& culprit) {
    std::ostringstream m;
    const char* nl = "<br />";
    m << "Invalid initial condition `" << culprit << "';" << nl
      << "Initial conditions must be given in the form `x(i)=k' where:"  << nl
      << "`i' is a non-negative integer" << nl
      << "`k' is a number (not a floating point) " << nl
      << "or a symbolic expression" << nl
      << "containing the parameters a, b,..., z" << nl
      << " different from `n', `x' and `e'." << endl;
    error(m.str());
}

#if 0
static void
invalid_list_initial_conditions(const string& s) {
  cerr << "Invalid list of initial conditions `" << s << "'.\n"
       << "If the initial conditions are more than one\n"
       << "they must be separated by a semicolon (;)." << endl;
  my_exit(1);
}
#endif

void
check_single_initial_condition(const string& cond,
			       std::map<index_type, Expr>&
			       initial_conditions) {
  // Special case: `cond' is a string containing only blanks.
  bool only_blank_spaces = true;
  for (unsigned i = cond.size(); i-- > 0; )
    if (cond[i] != ' ') {
      only_blank_spaces = false;
      break;
    }
  if (only_blank_spaces)
    return;

  string::size_type eq_pos = cond.find('=');
  if (eq_pos == string::npos)
    invalid_initial_condition(cond);
  string lhs(cond, 0, eq_pos);
  string rhs(cond, eq_pos+1);
  Expr l;
  if (!parse_expression(lhs, l))
    invalid_initial_condition(cond);
  if (!l.is_the_x_function())
    invalid_initial_condition(cond);
  Number index;
  if (!l.arg(0).is_a_number(index)
      || !index.is_nonnegative_integer()
      || index > LONG_MAX)
    invalid_initial_condition(cond);
  Expr r;
  if (!parse_expression(rhs, r))
    invalid_initial_condition(cond);
  if (invalid_initial_condition(r))
    invalid_initial_condition(cond);
  // Insert the pair (index, r), which represents the initial
  // condition `x(index) = r', in the map `initial_conditions'.
  initial_conditions.insert(std::map<index_type, Expr>
			    ::value_type(index.to_unsigned_int(), r));
}

void
check_sintax_initial_conditions(const string& list_i_c,
				std::map<index_type, Expr>&
				initial_conditions) {
  string::size_type pos_first = 0;
  for (unsigned int i = 0; i < list_i_c.size(); ++i) {
    string::size_type pos_second = list_i_c.find(";", pos_first);
    string cond;
    if (pos_second == string::npos)
      cond = list_i_c.substr(pos_first);
    else
      cond = list_i_c.substr(pos_first, pos_second-pos_first);
    check_single_initial_condition(cond, initial_conditions);
    if (pos_second == string::npos)
      break;
    pos_first = pos_second+1;
  }
}

bool
can_split_at(const char c) {
  switch (c) {
  case '+':
  case '-':
  case '*':
  case '/':
    return true;
  default:
    return false;
  }
}

void
portray(std::ostream& s, const Expr& e) {
  std::ostringstream b;
  b << e;
  const std::string& t = b.str();
  int l = t.length();
  // Ideal line length.
  const int ill = 80;
  // `x' will be the ideal split point.
  int x;
  // `y' will be the split point found, if any.
  int y;
  // Loop until there are characters after the ideal split point.
  for (x = ill; x < l; x = y+ill) {
    // Look for a split point at or before `x'.
    for (y = x; y > x-ill; --y)
      if (can_split_at(t[y]))
	goto found;
    // Split point not found before `x': try after.
    for (y = x+1; y < l; ++y)
      if (can_split_at(t[y]))
	goto found;
    // If no split point has been found, will print the rest of `t'.
  found:
    s << t.substr(x-ill, y-(x-ill)) << std::endl;
  }
  // Print what is left.
  s << t.substr(x-ill, l-(x-ill));
}

int
main() try {
  // Limit the amount of resources we may consume.
  set_alarm_on_cpu_time(MAX_SECONDS_OF_CPU_TIME, my_timeout);
  limit_virtual_memory(MAX_VIRTUAL_MEMORY);

  // Even if we run low on heap memory,
  // we would still like to run to completion.
  safety = new char[500000];
  old_new_handler = std::set_new_handler(&my_new_handler);

  std::set_unexpected(my_unexpected_exception);
  std::set_terminate(my_uncaught_exception);

#ifdef HAVE_NICE
  // Be very nice to the system.
  if (nice(20) == -1)
    error("cannot be nice");
#endif

  // Get the expression (rhs of the recurrence), if any.
  const_form_iterator expr = cgi.getElement("expr");

  // Get options, if any.
  bool verify = false;
  vector<FormEntry> options;
  cgi.getElement("options", options);
  for(string::size_type i = 0; i < options.size(); ++i) {
    const string& option = options[i].getValue();
    if (option == "verify")
      verify = true;
    else
      error("internal error");
  }

  if(expr == (*cgi).end() || expr->isEmpty())
    error("you did not type an expression for the recurrence!!!");

  init_symbols();

  Expr rhs = Expr(**expr, symbols);
  Recurrence recurrence(rhs);

  // Get the expression (initial conditions for the recurrence), if any.
  const_form_iterator i_c = cgi.getElement("i_c");

  string list_i_c;
  if(i_c == (*cgi).end() || i_c->isEmpty())
    list_i_c = "";
  else
    list_i_c = **i_c;

  std::map<index_type, Expr> initial_conditions;
  if (!list_i_c.empty())
    check_sintax_initial_conditions(list_i_c, initial_conditions);
  
  if (!initial_conditions.empty())
    recurrence.set_initial_conditions(initial_conditions);

  bool have_exact_solution = false;
  bool have_lower_bound = false;
  bool have_upper_bound = false;

  bool have_verified_exact_solution = false;
  bool have_verified_lower_bound = false;
  bool have_verified_upper_bound = false;

  Expr exact_solution;
  Expr lower_bound;
  Expr upper_bound;
  unsigned int first_valid_index_for_solution = 0;
  string Sc_function;

  time_unit_t start_time;

  long solution_time_msecs = 0;
  long verification_time_msecs = 0;

  try {
    start_time = get_time();
    switch (recurrence.compute_exact_solution()) {
    case Recurrence::SUCCESS:
      have_exact_solution = true;
      recurrence.exact_solution(exact_solution);
      exact_solution
        = recurrence.substitute_auxiliary_definitions(exact_solution);
      first_valid_index_for_solution
	= recurrence.first_valid_index_for_solution();
      solution_time_msecs += time_unit_to_msecs(get_time() - start_time);

      start_time = get_time();
      if (verify
	  &&
	  recurrence.verify_exact_solution() == Recurrence::PROVABLY_CORRECT) {
	have_verified_exact_solution = true;
	verification_time_msecs += time_unit_to_msecs(get_time() - start_time);
      }
      goto done;
      break;
    case Recurrence::UNSOLVABLE_RECURRENCE:
      error("this recurrence is unsolvable");
      break;
    case Recurrence::INDETERMINATE_RECURRENCE:
      error("this recurrence has infinitely many solutions");
      break;
    case Recurrence::MALFORMED_RECURRENCE:
      error("this recurrence is malformed");
      break;
    case Recurrence::TOO_COMPLEX:
    default:
      ;
    }
  }
  catch (const char* s) {
    error(s);
  }

  try {
    start_time = get_time();
    switch (recurrence.compute_lower_bound()) {
    case Recurrence::SUCCESS:
      have_lower_bound = true;
      recurrence.lower_bound(lower_bound);
      first_valid_index_for_solution
	= recurrence.first_valid_index_for_solution();
      Sc_function = recurrence.definition_Sc();
      solution_time_msecs += time_unit_to_msecs(get_time() - start_time);

      start_time = get_time();
      if (verify
	  &&
	  recurrence.verify_lower_bound() == Recurrence::PROVABLY_CORRECT) {
	have_verified_lower_bound = true;
	verification_time_msecs += time_unit_to_msecs(get_time() - start_time);
      }
      break;
    case Recurrence::UNSOLVABLE_RECURRENCE:
      error("this recurrence is unsolvable");
      break;
    case Recurrence::INDETERMINATE_RECURRENCE:
      error("this recurrence has infinitely many solutions");
      break;
    case Recurrence::MALFORMED_RECURRENCE:
      error("this recurrence is malformed");
      break;
    case Recurrence::TOO_COMPLEX:
    default:
      ;
    }
  }
  catch (const char* s) {
    error(s);
  }

  try {
    start_time = get_time();
    switch (recurrence.compute_upper_bound()) {
    case Recurrence::SUCCESS:
      have_upper_bound = true;
      recurrence.upper_bound(upper_bound);
      first_valid_index_for_solution
	= recurrence.first_valid_index_for_solution();
      Sc_function = recurrence.definition_Sc();
      solution_time_msecs += time_unit_to_msecs(get_time() - start_time);

      start_time = get_time();
      if (verify
	  &&
	  recurrence.verify_upper_bound() == Recurrence::PROVABLY_CORRECT) {
	have_verified_upper_bound = true;
	verification_time_msecs += time_unit_to_msecs(get_time() - start_time);
      }
      break;
    case Recurrence::UNSOLVABLE_RECURRENCE:
      error("this recurrence is unsolvable");
      break;
    case Recurrence::INDETERMINATE_RECURRENCE:
      error("this recurrence has infinitely many solutions");
      break;
    case Recurrence::MALFORMED_RECURRENCE:
      error("this recurrence is malformed");
      break;
    case Recurrence::TOO_COMPLEX:
    default:
      ;
    }
  }
  catch (const char* s) {
    error(s);
  }

  if (!have_exact_solution
      && !have_lower_bound
      && !have_upper_bound)
    error("sorry, this is too difficult");

 done:

  // Output the HTTP headers for an HTML document, and the HTML 4.0 DTD info.
  cout << HTTPHTMLHeader() << HTMLDoctype(HTMLDoctype::eStrict) << endl
       << html().set("lang", "en").set("dir", "ltr") << endl;

  // Set up the page's header and title.
  cout << head() << endl;

    // Output the style sheet portion of the header
  cout << style() << comment() << endl
       << "body { color: black; background-color: white; }" << endl
       << "hr.half { width: 60%; align: center; }" << endl
       << "span.red, strong.red { color: red; }" << endl
       << "span.green, strong.green { color: green; }" << endl
       << "div.smaller { font-size: small; }" << endl
       << "div.bigger { font-size: large; }" << endl
       << "div.notice { border: solid thin; padding: 1em; margin: 1em 0; "
       << "background: #ddd; }" << endl
       << "span.blue { color: blue; }" << endl
       << "col.title { color: white; background-color: black; "
       << "font-weight: bold; text-align: center; }" << endl
       << "col.data { background-color: #DDD; text-align: left; }" << endl
       << "td.data, tr.data { background-color: #ddd; text-align: left; }"
       << endl
       << "td.grayspecial { background-color: #ddd; text-align: left; }"
       << endl
       << "td.ltgray, tr.ltgray { background-color: #ddd; }" << endl
       << "td.dkgray, tr.dkgray { background-color: #bbb; }" << endl
       << "col.black, td.black, td.title, tr.title { color: white; " 
       << "background-color: black; font-weight: bold; text-align: center; }"
       << endl
       << "col.gray, td.gray { background-color: #ddd; text-align: center; }"
       << endl
       << "table.cgi { left-margin: auto; right-margin: auto; width: 90%; }"
       << endl
       << comment() << style() << endl;

  cout << title() << "PURRS Demo Results" 
       << title() << endl;
  cout << meta()
    .set("name", "author")
    .set("content", "PURRS development team") << endl;

  cout << head() << endl;

  // Start the HTML body.
  cout << body() << endl
       << h1() << "PURRS Demo " << span("Results", set("class", "green"))
       << h1() << endl;

  cout << h2();
  // Print the recurrence.
  if (have_exact_solution) {
    if (have_verified_exact_solution)
      cout << span("Verified", set("class", "red")) << " exact solution";
    else
      cout << "Exact solution";
  }
  else if (have_verified_lower_bound || have_verified_upper_bound)
    cout << span("Verified", set("class", "red")) << " bounds";
  else
    cout << "Bounds";
  cout << " for x(n) = " << rhs
       << br() << endl;

  // Print the possibly initial conditions.
  if ((have_exact_solution || have_lower_bound || have_upper_bound)
      && !initial_conditions.empty()) {
    cout << "for the initial conditions" << br() << endl;
    for (std::map<index_type, Expr>::const_iterator i
	   = initial_conditions.begin(),
	   initial_conditions_end = initial_conditions.end();
	 i != initial_conditions_end; ++i)
      cout << "  x(" << i->first << ")"
	   << " = " << i->second << br() << endl;
  }
  cout << h2() << endl;

  // Print the solution or the bounds.
  if (have_exact_solution) {
    if (have_verified_exact_solution)
      mark_verified_solution();
    cout << "x(n) = ";
    portray(cout, exact_solution);
    cout << endl;
  }
  else {
    if (have_lower_bound) {
      if (have_verified_lower_bound)
	mark_verified_solution();
      cout << "x(n) >= ";
      portray(cout, lower_bound);
      cout << br() << endl;
    }
    if (have_upper_bound) {
      if (have_verified_upper_bound)
	mark_verified_solution();
      cout << "x(n) <= ";
      portray(cout, upper_bound);
      cout << endl;
    }
  }

  if (have_exact_solution || have_lower_bound || have_upper_bound)
    cout << br() << "for each n >= " << first_valid_index_for_solution << endl;

  if (have_lower_bound || have_upper_bound) {
    // In the bound occurs the symbolic initial condition `x(1)'.
    if (Sc_function.empty()) {
      if (initial_conditions.empty())
	cout << ", assuming x(1) >= 0," << endl;
    }
    // In the bound occurs the symbolic initial condition `x(Sc(n, b))'.
    else {
      if (initial_conditions.empty())
	cout << ", assuming x(" << Sc_function.substr(0, 8) << ") >= 0,"
	     << endl;
      cout << br() << "where " << Sc_function << endl;
      if (!initial_conditions.empty())
	cout << br() << "and" << endl;
      for (std::map<index_type, Expr>::const_iterator i
	     = initial_conditions.begin(),
	     initial_conditions_end = initial_conditions.end(), j = i;
	   i != initial_conditions_end; ++i) {
	cout << "  x(" << i->first << ")"
	     << " = " << i->second;
	if (++j != initial_conditions_end)
	  cout << ", " << endl;
	else
	  cout << ". " << endl;
      }
    }
  }
  
  // Get a pointer to the environment.
  const CgiEnvironment& env = cgi.getEnvironment();

  // Timings and thank you.
  cout << br() << br()
       <<cgicc::div().set("align", "center").set("class", "bigger") << endl;
  cout << "Computing the "
       << (have_exact_solution ? "exact solution"
	   : ((have_lower_bound && have_upper_bound)
	      ? "approximations" : "approximation"))
       << " took about " << solution_time_msecs << " ms of CPU time";
  if (have_verified_exact_solution
      || have_verified_lower_bound || have_verified_upper_bound) {
  cout << ";"
       << endl
       << br()
       << "verifying "
       << ((have_lower_bound && have_upper_bound) ? "them" : "it")
       << " took about " << verification_time_msecs << " ms of CPU time.";
  }
  else {
    cout << "." << endl;
  }
  cout << br() << br() << endl;
  string host = env.getRemoteHost();
  if (host.empty())
    host = env.getRemoteAddr();
  cout << "Thanks for using PURRS, "
       << host << "!" << endl
       << cgicc::div() << endl;

  footer();
  return 0;
}
catch (std::bad_alloc&) {
  error("memory limit exceeded");
}
catch (std::exception &e) {
  error(e.what());
}
catch (const char* s) {
  error(s);
}
