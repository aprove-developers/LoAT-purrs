/* Boundary class implementation (non-inline functions).
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

#include "Boundary.defs.hh"
#if OUTLINE
#include "Boundary.inlines.hh"
#endif

#include <iostream>

using std::ws;

#if USE_DIRECTED_ROUNDING
fpu_ctrl_reg_type _saved_fpcr;	// For round_save and round_restore
#endif

// Output operator for lower boundary
ostream&
operator << (ostream& s, const LBoundary& b)
{
  s.precision(60);
  return s << (b.is_closed() ? "[" : "(") << b.value;
}

// Input operator for lower boundary
istream&
operator >> (istream& s, LBoundary& b)
{
  char ch;
  char v[256];
  LBoundary::Open_Closed f;
  s >> ws;
  s >> ch;
  if (ch == '(')
    f=LBoundary::OPEN;
  else if (ch == '[')
    f=LBoundary::CLOSED;
  else {
    s.putback(ch);
    //    s.clear(_bad);
    return s;
  }
  s >> ws;
  s >> v;
  b = LBoundary(v, f);
  return s;
}

// Output operator for upper boundary
ostream&
operator << (ostream& s, const UBoundary& b)
{
  s.precision(60);
  return s << b.value << (b.is_closed() ? "]" : ")");
}

// Input operator for upper boundary
istream&
operator >> (istream& s, UBoundary& b)
{
  char ch;
  char v[256];
  UBoundary::Open_Closed f;
  s >> ws;
  s >> v;
  s >> ws;
  s >> ch;
  if (ch == ')')
    f=UBoundary::OPEN;
  else if (ch == ']')
    f=UBoundary::CLOSED;
  else {
    s.putback(ch);
    //    s.clear(_bad);
    return s;
  }
  b = UBoundary(v, f);
  return s;
}
