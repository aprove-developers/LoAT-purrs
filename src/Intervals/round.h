/* Definitions to set the rounding direction.
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

/* Machine and OS dependant stuff */

/*
#if __GNUG__
extern "C" {
#endif
*/

#if USE_DIRECTED_ROUNDING
#if defined(__GNUC__) && (defined(i386) || defined(sparc) || defined(m68k) || defined(alpha))
#if defined(i386)
#include "i386.h"
#elif defined(sparc)
#include "sparc.h"
#elif defined(m68k)
#include "m68k.h"
#elif defined(alpha)
#include "alpha.h"
#endif /* defined(alpha) */

#elif defined(HAVE_FENV_H)
#include "c99.h"

#elif defined(HAVE_IEEE_LIBRARY)
#include <ieeefp.h>
static char __ieee_direction[10];
#define round_save() ieee_flags("get", "direction", 0, _ieee_direction)
#define round_restore() ieee_flags("set", "direction", _ieee_direction, 0)
#define round_up() ieee_flags("set", "direction", "positive", 0)
#define round_down() ieee_flags("set", "direction", "negative", 0)
#else
#warning "I don't know how to set rounding mode"
#define round_save()
#define round_restore()
#define round_up()
#define round_down()
#endif

#else

#define round_save()
#define round_restore()
#define round_up()
#define round_down()

#endif
