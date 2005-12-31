/* FPU control for compilers conforming to the C99 standard.
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
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the Parma Interval Library
site: http://www.cs.unipr.it/pil/ . */

#include <cassert>
#include <fenv.h>

typedef int fpu_ctrl_reg_type;

#define FPU_INEXACT_RESULT	0x20
typedef unsigned short int fpu_stat_reg_type;

__inline__ fpu_stat_reg_type
fpu_get_status_reg() {
  fpu_stat_reg_type stat_reg;
  __asm__ __volatile__ ("fnstsw %0" : "=a" (stat_reg));
  return stat_reg;
}

__inline__ void
fpu_reset_status_reg() {
  __asm__ ("fclex" : /* No outputs.  */);
}

extern fpu_ctrl_reg_type _saved_fpcr;

__inline__ void
round_save() {
  _saved_fpcr = fegetround();
  assert(_saved_fpcr == FE_TONEAREST);
}

__inline__ void
round_restore() {
  fesetround(_saved_fpcr);
}

__inline__ void
round_up() {
  fesetround(FE_UPWARD);
}

__inline__ void
round_down() {
  fesetround(FE_DOWNWARD);
}

#if CHECK_INEXACT_RESULT

__inline__ void
reset_inexact() {
  //feclearexcept(FE_INEXACT);
  feclearexcept(FE_ALL_EXCEPT);
}

__inline__ int
is_inexact() {
 return fetestexcept(FE_INEXACT) & FE_INEXACT;
}

#endif /* CHECK_INEXACT_RESULT */
