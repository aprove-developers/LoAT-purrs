/* Bool3 class implementation (non-inline functions).
   Copyright (C) 2002 Roberto Bagnara <bagnara@cs.unipr.it>

This file is part of the Parma Interval Library (PIL).

The PIL is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The PIL is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.

For the most up-to-date information see the Parma Interval Library
site: http://www.cs.unipr.it/pil/ . */

#include <config.h>

#include "Bool3.defs.hh"
#if OUTLINE
#include "Bool3.inlines.hh"
#endif

#include <iostream>

ostream&
operator << (ostream& s, Bool3 b) {
  return s << (b == ALWAYS ? "ALWAYS" : b == NEVER ? "NEVER" : "MAYBE");
}
