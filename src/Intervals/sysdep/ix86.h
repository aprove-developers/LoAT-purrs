/* FPU control for the ix86 architecture.
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

#include <cassert>

#define FPU_ROUNDING_MASK	0xc00
#define FPU_ROUNDING_NEAREST	0x0
#define FPU_ROUNDING_DOWN	0x400
#define FPU_ROUNDING_UP		0x800
#define FPU_ROUNDING_ZERO	0xc00

typedef unsigned short int fpu_ctrl_reg_type;

__inline__ fpu_ctrl_reg_type
fpu_get_control_reg() {
  __volatile__ fpu_ctrl_reg_type ctrl_reg; /* Force in memory */
  __asm__ __volatile__ ("fnstcw %0" : "=m" (ctrl_reg));
  return ctrl_reg;
}

__inline__ void
fpu_set_control_reg(__volatile__ fpu_ctrl_reg_type ctrl_reg) {
  __asm__ __volatile__ ("fldcw %0" : : "m" (ctrl_reg));
}

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
  _saved_fpcr = fpu_get_control_reg();
  assert((_saved_fpcr & FPU_ROUNDING_MASK) == FPU_ROUNDING_NEAREST);
}

__inline__ void
round_restore() {
  fpu_set_control_reg(_saved_fpcr);
}

__inline__ void
round_up() {
  fpu_set_control_reg((_saved_fpcr&(~FPU_ROUNDING_MASK))|FPU_ROUNDING_UP);
}

__inline__ void
round_down() {
  fpu_set_control_reg((_saved_fpcr&(~FPU_ROUNDING_MASK))|FPU_ROUNDING_DOWN);
}

#if CHECK_INEXACT_RESULT

#define IE_ASSIGN 1
#define IE_STRTOD 0
#define IE_ADD 1
#define IE_SUB 1
#define IE_MUL 1
#define IE_DIV 1
#define IE_EXP 1
#define IE_LOG 1
#define IE_SQRT 1
#define IE_SIN 1
#define IE_COS 1
#define IE_TAN 1
#define IE_POW 0

__inline__ void
reset_inexact() {
  fpu_reset_status_reg();
}

__inline__ int
is_inexact() {
 return fpu_get_status_reg() & FPU_INEXACT_RESULT;
}

#endif /* CHECK_INEXACT_RESULT */
