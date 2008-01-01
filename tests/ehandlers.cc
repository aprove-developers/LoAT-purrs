/* Implementation of exception handlers useful for debugging purposes.
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

#include <exception>
#include <iostream>
#include <cstdlib>

using namespace std;

static void
my_unexpected_exception() {
  cerr << "unexpected exception thrown" << endl;
  exit(1);
}

static void
my_uncaught_exception() {
  cerr << "uncaught exception" << endl;
  exit(1);
}

void
set_handlers() {
  set_unexpected(my_unexpected_exception);
  set_terminate(my_uncaught_exception);
}
