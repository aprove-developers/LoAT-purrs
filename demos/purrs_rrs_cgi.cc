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

  // Get the expression, if any.
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
    error("you did not type anything!!!");

  init_symbols();

  Expr rhs = Expr(**expr, symbols);
  Recurrence recurrence(rhs);

  bool have_exact_solution = false;
  bool have_lower_bound = false;
  bool have_upper_bound = false;

  bool have_verified_exact_solution = false;
  bool have_verified_lower_bound = false;
  bool have_verified_upper_bound = false;

  Expr exact_solution;
  Expr lower_bound;
  Expr upper_bound;
  unsigned int first_valid_index_for_solution;

  time_unit_t start_time;

  long solution_time_usecs = 0;
  long verification_time_usecs = 0;

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
      solution_time_usecs += time_unit_to_usecs(get_time() - start_time);

      start_time = get_time();
      if (verify
	  &&
	  recurrence.verify_exact_solution() == Recurrence::PROVABLY_CORRECT) {
	have_verified_exact_solution = true;
	verification_time_usecs += time_unit_to_usecs(get_time() - start_time);
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
      solution_time_usecs += time_unit_to_usecs(get_time() - start_time);

      start_time = get_time();
      if (verify
	  &&
	  recurrence.verify_lower_bound() == Recurrence::PROVABLY_CORRECT) {
	have_verified_lower_bound = true;
	verification_time_usecs += time_unit_to_usecs(get_time() - start_time);
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
      solution_time_usecs += time_unit_to_usecs(get_time() - start_time);

      start_time = get_time();
      if (verify
	  &&
	  recurrence.verify_upper_bound() == Recurrence::PROVABLY_CORRECT) {
	have_verified_upper_bound = true;
	verification_time_usecs += time_unit_to_usecs(get_time() - start_time);
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
  cout<< " for x(n) = " << rhs
      << h2() << endl;

  if (have_exact_solution) {
    if (have_verified_exact_solution)
      mark_verified_solution();
    cout << "x(n) = " << exact_solution << endl;
  }
  else {
    if (have_lower_bound) {
      if (have_verified_lower_bound)
	mark_verified_solution();
      cout << "x(n) >= " << lower_bound << br() << endl;
    }
    if (have_upper_bound) {
      if (have_verified_upper_bound)
	mark_verified_solution();
      cout << "x(n) <= " << upper_bound << endl;
    }
  }
  if (have_exact_solution || have_lower_bound || have_upper_bound)
    cout << br() << "for each n >= " << first_valid_index_for_solution << endl;
  
  // Get a pointer to the environment.
  const CgiEnvironment& env = cgi.getEnvironment();

  // Timings and thank you.
  cout << br() << br()
       <<cgicc::div().set("align", "center").set("class", "bigger") << endl;
  cout << "Computing the "
       << (have_exact_solution ? "exact solution"
	   : ((have_lower_bound && have_upper_bound)
	      ? "approximations" : "approximation"))
       << " took about " << solution_time_usecs << " us of CPU time";
  if (have_verified_exact_solution
      || have_verified_lower_bound || have_verified_upper_bound) {
  cout << ";"
       << endl
       << br()
       << "verifying "
       << ((have_lower_bound && have_upper_bound) ? "them" : "it")
       << " took about " << verification_time_usecs << " us of CPU time.";
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
