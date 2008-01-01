/* FPU control for the M68K architecture.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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

#define FPU_ROUNDING_MASK	3 << 4
#define FPU_ROUNDING_NEAREST	0 << 4
#define FPU_ROUNDING_ZERO	1 << 4
#define FPU_ROUNDING_DOWN	2 << 4
#define FPU_ROUNDING_UP		3 << 4

typedef unsigned int fpu_ctrl_reg_type;

__inline__ fpu_ctrl_reg_type
fpu_get_control_reg() {
  fpu_ctrl_reg_type ctrl_reg;
  __asm__ __volatile__ ("fmove%.l %!,%0" : "=dm" (ctrl_reg));
  return ctrl_reg;
}
__inline__ void
fpu_set_control_reg(fpu_ctrl_reg_type ctrl_reg) {
  __asm__ __volatile__ ("fmove%.l %0,%!" : : "dm" (ctrl_reg));
}

#define FPU_INEXACT_RESULT 1 << 3

typedef unsigned long int fpu_stat_reg_type;

__inline__ fpu_stat_reg_type
fpu_get_status_reg() {
  fpu_stat_reg_type stat_reg;
  __asm__ __volatile__ ("fmove%.l %/fpsr,%0" : "=dm" (stat_reg));
  return stat_reg;
}

__inline__ void
fpu_set_status_reg(fpu_stat_reg_type stat_reg) {
  __asm__ __volatile__ ("fmove%.l %0,%/fpsr" : : "dm" (stat_reg));
}

extern fpu_ctrl_reg_type _saved_fpcr;

__inline__ void
round_save() {
  _saved_fpcr = fpu_get_control_reg();
}

__inline__ void
round_restore() {
  fpu_set_control_reg(_saved_fpcr);
}

__inline__ void
round_up() {
  fpu_set_control_reg((_saved_fpcr & (~FPU_ROUNDING_MASK))
		      | FPU_ROUNDING_UP);
}

__inline__ void
round_down() {
  fpu_set_control_reg((_saved_fpcr & (~FPU_ROUNDING_MASK))
		      | FPU_ROUNDING_DOWN);
}

#if CHECK_INEXACT_RESULT
__inline__ void
reset_inexact() {
  /* We could reset only the flag FPU_INEXACT_RESULT,
     but then we would have to read the status register too.
  */
  fpu_set_status_reg(0);
}

__inline__ bool
is_inexact() {
 return fpu_get_status_reg() & FPU_INEXACT_RESULT;
}

#endif /* CHECK_INEXACT_RESULT */
