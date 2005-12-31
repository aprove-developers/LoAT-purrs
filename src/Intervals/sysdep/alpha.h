// FPU control for the Alpha architecture.

// Copyright (C) 2002 Roberto Bagnara <bagnara@cs.unipr.it>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

#define FPU_ROUNDING_MASK	(3UL << 58)
#define FPU_ROUNDING_ZERO	(0UL << 58)
#define FPU_ROUNDING_DOWN	(1UL << 58)
#define FPU_ROUNDING_NEAREST	(2UL << 58)
#define FPU_ROUNDING_UP		(3UL << 58)

typedef unsigned long fpu_ctrl_reg_type;

__inline__ fpu_ctrl_reg_type
fpu_get_control_reg() {
  fpu_ctrl_reg_type ctrl_reg;
  __asm__ __volatile__ ("excb; mf_fpcr %0" : "=f" (ctrl_reg));
  return ctrl_reg;
}
__inline__ void
fpu_set_control_reg(fpu_ctrl_reg_type ctrl_reg) {
  __asm__ __volatile__ ("mt_fpcr %0; excb" : : "f" (ctrl_reg));
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
}

__inline__ bool
is_inexact() {
  /* The only safe thing we can do is pretend the result is always inexact. */
  return true;
}

#endif /* CHECK_INEXACT_RESULT */
