/* Simple CGI program to test the algebraic equation solver.
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

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include <cgicc/CgiDefs.h>
#include <cgicc/Cgicc.h>
#include <cgicc/HTTPHeaders.h>
#include <cgicc/HTMLClasses.h>

#include "purrs_install.hh"

#if !defined(HAVE_GETRLIMIT) || !defined(HAVE_SETRLIMIT)
#error "We must have a way of limiting time and space!"
#endif

using std::cout;
using std::endl;
using std::string;

using namespace GiNaC;
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
  cout << a("PURRS").set("href", "http://www.cs.unipr.it/purrs/")
       << span(" algebraic equation solver", set("class", "red")) << br()
       << " by the "
       << a("PURRS development team")
    .set("href", "mailto:purrs-devel@cs.unipr.it") << "." <<br()
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
  cout << HTTPHTMLHeader() << HTMLDoctype(HTMLDoctype::eStrict) << endl;
  cout << html().set("lang", "en").set("dir", "ltr") << endl;

  // Set up the page's header and title.
  cout << head() << endl;

  // Output the style sheet portion of the header
  cout << style() << comment() << endl;
  cout << "body { color: black; background-color: white; }" << endl;
  cout << "hr.half { width: 60%; align: center; }" << endl;
  cout << "span.red, STRONG.red { color: red; }" << endl;
  cout << "div.notice { border: solid thin; padding: 1em; margin: 1em 0; "
       << "background: #ddd; }" << endl;

  cout << comment() << style() << endl;

  cout << title("PURRS Demo Error") << endl;
  cout << meta().set("name", "author")
    .set("content", "PURRS development team") << endl;
  cout << head() << endl;
    
  cout << body() << endl;
    
  cout << h1() << "PURRS Demo " << span("Error", set("class", "red"))
       << h1() << endl; 
  
  cout << cgicc::div().set("align", "center").set("class", "notice") << endl;

  cout << h2(message) << endl;
  cout << cgicc::div() << endl;

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

  // Get the expression, if any.
  const_form_iterator expr = cgi.getElement("expr");

  if(expr == (*cgi).end() || expr->isEmpty())
    error("you did not type anything!!!");

#if HAVE_GETTIMEOFDAY
  timeval start;
  gettimeofday(&start, NULL);
#endif

  GSymbol x("x");
  GExpr p = GExpr(**expr, lst(x));
  if (p == GExpr(0)
      || !p.info(info_flags::integer_polynomial)
      || p.info(info_flags::numeric)) {
    std::string message;
    std::ostringstream s(message);
    s << "you call '<TT>" << p << "</TT>' a polynomial in <TT>x</TT> "
      << "with integer coefficients?";
    error(s.str());
  }

  std::vector<Polynomial_Root> roots;
  bool all_distinct;
  if (!find_roots(p, x, roots, all_distinct))
    error("sorry, this is too difficult.");


  // Output the HTTP headers for an HTML document, and the HTML 4.0 DTD info.
  cout << HTTPHTMLHeader() << HTMLDoctype(HTMLDoctype::eStrict) << endl;
  cout << html().set("lang", "en").set("dir", "ltr") << endl;

  // Set up the page's header and title.
  cout << head() << endl;

  cout << title() << "PURRS Demo Results" 
       << title() << endl;
  cout << meta().set("name", "author").set("content", "PURRS") 
       << endl;

  cout << head() << endl;
    
  // Start the HTML body
  cout << body() << endl;

  cout << h1() << "PURRS Demo Results" << h1() << endl;

  // Get a pointer to the environment.
  const CgiEnvironment& env = cgi.getEnvironment();
    
  // Generic thank you message.
  cout << comment() << "This page has been generated by PURRS for "
       << env.getRemoteHost() << comment() << endl;
  cout << h4() << "Thanks for using PURRS, "
       << env.getRemoteHost()
       << '(' << env.getRemoteAddr() << ")!" << h4() << endl; 

  size_t n = roots.size();
  for (size_t i = 0; i < n; ++i) {
    GExpr value = roots[i].value();
    GNumber multiplicity = roots[i].multiplicity();
    cout << "x_" << i+1 << " = " << value;
    if (multiplicity > 1)
      cout << " (multiplicity " << multiplicity << ")";
    cout << br() << endl;
    if (!is_a<numeric>(value))
      cout << "****  x_" << i+1 << " ~= "
	   << value.evalf() << br() << endl;
  }

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
