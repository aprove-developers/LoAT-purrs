/* To be written.
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

#ifndef PURRS_complint_hh
#define PURRS_complint_hh 1

#include <ginac/ginac.h>
#include "CInterval.types.hh"

namespace GiNaC {

const unsigned TINFO_complint = 0x42420001U;

class complint : public basic {
  GINAC_DECLARE_REGISTERED_CLASS(complint, basic)

public:
  complint(const CInterval& i);

  void complint::print(const print_context& c, unsigned level) const;

private:
  CInterval ci;
};

} // namespace GiNaC

#endif // !defined(PURRS_complint_hh)
