/* FPU control for the SPARC architecture.
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

#define FPU_ROUNDING_MASK	(3U << 30)
#define FPU_ROUNDING_NEAREST	(0U << 30)
#define FPU_ROUNDING_ZERO	(1U << 30)
#define FPU_ROUNDING_UP		(2U << 30)
#define FPU_ROUNDING_DOWN	(3U << 30)

#define FPU_INEXACT_RESULT	(1 << 5)

typedef unsigned long int fpu_ctrl_reg_type;

__inline__ fpu_ctrl_reg_type
fpu_get_control_reg() {
  __volatile__ fpu_ctrl_reg_type ctrl_reg;
  __asm__ __volatile__ ("st %%fsr,%0" : "=m" (ctrl_reg));
  return ctrl_reg;
}

__inline__ void
fpu_set_control_reg(__volatile__ fpu_ctrl_reg_type ctrl_reg) {
  __asm__ __volatile__ ("ld %0,%%fsr" : : "m" (ctrl_reg));
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

#if CHECK_INEXACT_RESULT
/* The sparc uses the same register for both the inexact flag
   and the rounding direction.  For efficiency, we reset
   the inexact flag when setting the rounding direction.
*/
__inline__ void
round_up() {
  fpu_set_control_reg((_saved_fpcr & (~FPU_ROUNDING_MASK|FPU_INEXACT_RESULT))
		      | FPU_ROUNDING_UP);
}

__inline__ void
round_down() {
  fpu_set_control_reg((_saved_fpcr & (~FPU_ROUNDING_MASK|FPU_INEXACT_RESULT))
		      | FPU_ROUNDING_DOWN);
}

__inline__ void
reset_inexact() {
}

__inline__ bool
is_inexact() {
  return fpu_get_control_reg() & FPU_INEXACT_RESULT;
}

#else /* !CHECK_INEXACT_RESULT */

__inline__ void
round_up() {
  fpu_set_control_reg((_saved_fpcr & (~FPU_ROUNDING_MASK))|FPU_ROUNDING_UP);
}

__inline__ void
round_down() {
  fpu_set_control_reg((_saved_fpcr & (~FPU_ROUNDING_MASK))|FPU_ROUNDING_DOWN);
}

#endif /* CHECK_INEXACT_RESULT */
