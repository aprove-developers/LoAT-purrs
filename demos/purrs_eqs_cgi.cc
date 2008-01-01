/* Simple CGI program to test the algebraic equation solver.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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

#if !defined(HAVE_GETRLIMIT) || !defined(HAVE_SETRLIMIT)
#error "We must have a way of limiting time and space!"
#endif

using std::cout;
using std::endl;
using std::string;

using namespace Parma_Recurrence_Relation_Solver;
using namespace cgicc;

const int MAX_SECONDS_OF_CPU_TIME = 2;
const int MAX_VIRTUAL_MEMORY = 8*1024*1024;

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
  cout << a("PURRS ").set("href", "http://www.cs.unipr.it/purrs/")
       << span("algebraic equation solver", set("class", "red")) << br()
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
  html::reset(); 	head::reset(); 		body::reset();
  title::reset(); 	h1::reset(); 		h4::reset();
  comment::reset(); 	td::reset(); 		tr::reset(); 
  table::reset();	cgicc::div::reset();	p::reset(); 
  a::reset();		h2::reset(); 		colgroup::reset();

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

  if(expr == (*cgi).end() || expr->isEmpty())
    error("you did not type anything!!!");

  Symbol x("x");
  Expr p = Expr(**expr, Expr_List(x));
  if (p.is_zero()
      || !p.is_integer_polynomial(x)
      || p.is_a_number()) {
    std::ostringstream s;
    s << "you call '<TT>" << p << "</TT>' a polynomial in <TT>x</TT> "
      << "with integer coefficients?";
    error(s.str());
  }

#if HAVE_GETRUSAGE
  timeval start;
  rusage rsg;
  if (getrusage(RUSAGE_SELF, &rsg) != 0)
    error("getrusage failed");
  else
    start = rsg.ru_utime;
#endif

  std::vector<Polynomial_Root> roots;
  bool all_distinct;
  if (!find_roots(p, x, roots, all_distinct))
    error("sorry, this is too difficult");

#if HAVE_GETRUSAGE
  timeval end;
  if (getrusage(RUSAGE_SELF, &rsg) != 0)
    error("getrusage failed");
  else
    end = rsg.ru_utime;

  long us_of_cpu_time = ((end.tv_sec - start.tv_sec) * 1000000)
    + (end.tv_usec - start.tv_usec);
#endif

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

  cout << h2() << "Roots of " << p << h2() << endl;

  cout << cgicc::div().set("align", "center") << endl
       << table()
    .set("border", "0").set("rules", "none").set("frame", "void")
    .set("cellspacing", "2").set("cellpadding", "2")
    .set("class", "cgi") << endl
       << colgroup().set("span", "2") << endl
       << col().set("align", "center").set("span", "2") << endl
       << colgroup() << endl
       << tr().set("class", "title") << td("Multiplicity") 
       << td("Value or Approximation") << tr() << endl;
  
  size_t n = roots.size();
  for (size_t i = 0; i < n; ++i) {
    std::ostringstream m;
    m << roots[i].multiplicity();
    std::ostringstream v;
    v << roots[i].value();
    cout << tr().set("class", "data")
	 << td(m.str()) 
	 << td(v.str())
	 << tr() << endl;
  }
  cout << table() << cgicc::div() << endl;

  // Get a pointer to the environment.
  const CgiEnvironment& env = cgi.getEnvironment();

  // Timings and thank you.
  cout << br() << br()
       <<cgicc::div().set("align", "center").set("class", "bigger") << endl;
#if HAVE_GETRUSAGE
  cout << "The computation of roots took about "
       << (double) (us_of_cpu_time/1000000.0) << " s of CPU time."
       << br() << br() << endl;
#endif
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
