/* Implementation of the public function for approximating expressions.
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

#include "approximate.hh"
#include "approximate_impl.hh"

namespace PURRS = Parma_Recurrence_Relation_Solver;

namespace {

struct Identity {
  bool approximate(const PURRS::Symbol& s, PURRS::Expr& ae, CInterval&) const {
    ae = s;
    return false;
  }
};

}

PURRS::Expr
PURRS::approximate(const PURRS::Expr& e) {
  return generic_approximate(e, Identity());
}
