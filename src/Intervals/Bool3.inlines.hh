/* Bool3 class implementation: inline functions.
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

INLINE Bool3
lor(Bool3 a, Bool3 b) {
  return (a == ALWAYS || b == ALWAYS) ? ALWAYS
    : (a == NEVER && b == NEVER) ? NEVER
      : MAYBE;
}

INLINE Bool3
operator || (Bool3 a, Bool3 b) {
  return lor(a, b);
}

INLINE Bool3
land(Bool3 a, Bool3 b) {
  return
    (a == NEVER || b == NEVER) ? NEVER
    : (a == ALWAYS && b == ALWAYS) ? ALWAYS
    : MAYBE;
}

INLINE Bool3
operator && (Bool3 a, Bool3 b) {
  return land(a, b);
}

INLINE Bool3
lnot(Bool3 a) {
  return
    a == NEVER ? ALWAYS
    : a == ALWAYS ? NEVER
    : MAYBE;
}

INLINE Bool3
operator ! (Bool3 a) {
  return lnot(a);
}

INLINE bool
is_strong(Bool3 a) {
  return a == ALWAYS || a == NEVER;
}

INLINE bool
stronger(Bool3 a, Bool3 b) {
  return b == MAYBE && (a == ALWAYS || a == NEVER);
}

INLINE bool
stronger_eq(Bool3 a, Bool3 b) {
  return a == b || (b == MAYBE && (a == ALWAYS || a == NEVER));
}
