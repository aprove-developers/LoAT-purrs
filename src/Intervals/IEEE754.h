/* Definitions for the handling of IEEE 754 floating point numbers.
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

#if defined(__GNUG__)
extern "C" {
#endif

#if __BYTE_ORDER == BIG_ENDIAN
#define _pinf_bytes_d { 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
#define _ninf_bytes_d { 0xff, 0xf0, 0, 0, 0, 0, 0, 0 }
#define _pinf_bytes_f { 0x7f, 0x80, 0, 0 }
#define _ninf_bytes_f { 0xff, 0x80, 0, 0 }
#endif
#if __BYTE_ORDER == LITTLE_ENDIAN
#define _pinf_bytes_d { 0, 0, 0, 0, 0, 0, 0xf0, 0x7f }
#define _ninf_bytes_d { 0, 0, 0, 0, 0, 0, 0xf0, 0xff }
#define _pinf_bytes_f { 0, 0, 0x80, 0x7f }
#define _ninf_bytes_f { 0, 0, 0x80, 0xff }
#endif

#ifdef  __GNUC__
#define __u_d(p) (__extension__ ((union { unsigned char __c[sizeof(double)]; \
					    double __d; })        \
				 { p }).__d)
#define PINF_D __u_d(_pinf_bytes_d)
#define NINF_D __u_d(_ninf_bytes_d)
#define __u_f(p) (__extension__ ((union { unsigned char __c[sizeof(float)]; \
					    float __f; })        \
                  { p }).__f)
#define PINF_F __u_f(_pinf_bytes_f)
#define NINF_F __u_f(_ninf_bytes_f)
#else
static char _pinf_d[] = _pinf_bytes_d;
static char _ninf_d[] = _ninf_bytes_d;
static char _pinf_f[] = _pinf_bytes_f;
static char _ninf_f[] = _ninf_bytes_f;
#define PINF_D (*(double *) _pinf_d)
#define NINF_D (*(double *) _ninf_d)
#define PINF_F (*(float *) _pinf_f)
#define NINF_F (*(float *) _ninf_f)
#endif

#define MAX_ODD_D 9007199254740991.0
#define MAX_ODD_F 16777215.0

#define MAX_FRACT_D 4503599627370495.5
#define MAX_FRACT_F 8388607.5

#define _exponent_max_d 0x7ff
union ieee754_d
{
  double v;

  struct
    {
#if __BYTE_ORDER == BIG_ENDIAN
      unsigned int negative:1;
      unsigned int exponent:11;
      unsigned int mantissa0:20;
      unsigned int mantissa1:32;
#endif				/* Big endian.  */
#if __BYTE_ORDER == LITTLE_ENDIAN
      unsigned int mantissa1:32;
      unsigned int mantissa0:20;
      unsigned int exponent:11;
      unsigned int negative:1;
#endif				/* Little endian.  */
    } ieee;
};

#define _exponent_max_f 0xff
union ieee754_f
{
  float v;

  struct
    {
#if __BYTE_ORDER == BIG_ENDIAN
      unsigned int negative:1;
      unsigned int exponent:8;
      unsigned int mantissa:23;
#endif				/* Big endian.  */
#if __BYTE_ORDER == LITTLE_ENDIAN
      unsigned int mantissa:23;
      unsigned int exponent:8;
      unsigned int negative:1;
#endif				/* Little endian.  */
    } ieee;
};

inline
int is_infinity_d(double value)
{
  union ieee754_d u;
  u.v=value;
  /* An IEEE 754 infinity has an exponent with the
     maximum possible value and a zero mantissa.  */
  if (u.ieee.exponent == _exponent_max_d &&
      u.ieee.mantissa0 == 0 && u.ieee.mantissa1 == 0)
    return u.ieee.negative ? -1 : 1;
  return 0;
}

inline
int is_infinity_f(float value)
{
  union ieee754_f u;
  u.v=value;
  /* An IEEE 754 infinity has an exponent with the
     maximum possible value and a zero mantissa.  */
  if (u.ieee.exponent == _exponent_max_f &&
      u.ieee.mantissa == 0)
    return u.ieee.negative ? -1 : 1;
  return 0;
}

#if defined(__GNUG__)
}
#endif
