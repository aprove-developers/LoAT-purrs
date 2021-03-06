/* Init class declaration.
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

#ifndef PURRS_Init_defs_hh
#define PURRS_Init_defs_hh 1

#include "Init.types.hh"

#ifdef PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS
//! Class for initialization and finalization.
/*!
  <EM>Nifty Counter</EM> initialization class,
  ensuring that the library is initialized only once
  and before its first use.
  A count of the number of translation units using the library
  is maintained. A static object of Init type will be declared
  by each translation unit using the library.  As a result,
  only one of them will initialize and properly finalize
  the library.
*/
#endif // PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS

class Parma_Recurrence_Relation_Solver::Init {
private:
  //! Count the number of objects created.
  static unsigned int count;

public:
  //! Initializes the PURRS.
  Init();

  //! Finalizes the PURRS.
  ~Init();
};

#include "Init.inlines.hh"

#endif // !defined(PURRS_Init_defs_hh)
