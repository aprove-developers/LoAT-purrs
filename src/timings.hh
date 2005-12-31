/* Some machinery to perform CPU time measurements.
   Copyright (C) 2002 Roberto Bagnara <bagnara@cs.unipr.it>

This file is part of the Parma University's Recurrence Relation
Solver (PURRS).

The PURRS is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The PURRS is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#ifndef PURRS_timings_hh
#define PURRS_timings_hh 1

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#if USE_TSC
#include "tsc.hh"
#endif

namespace Parma_Recurrence_Relation_Solver {

#if USE_TSC

typedef tsc_t time_unit_t;

inline time_unit_t
get_time() {
  return read_tsc();
}

inline long
time_unit_to_usecs(time_unit_t t) {
  return tsc_to_usecs(t);
}

inline long
time_unit_to_msecs(time_unit_t t) {
  return tsc_to_msecs(t);
}

#elif HAVE_GETRUSAGE

typedef timeval time_unit_t;

inline time_unit_t
get_time() {
  rusage r;
  if (getrusage(RUSAGE_SELF, &r) != 0)
    r.ru_utime.tv_sec = r.ru_utime.tv_usec = 0;
  return r.ru_utime;
}

inline time_unit_t&
operator+=(time_unit_t& t1, const time_unit_t& t2) {
  t1.tv_sec += t2.tv_sec;
  t1.tv_usec += t2.tv_usec;
  if (t1.tv_usec >= 1000000) {
    t1.tv_sec += 1;
    t1.tv_usec -= 1000000;
  }
  return t1;
}

inline time_unit_t
operator-(const time_unit_t& t1, const time_unit_t& t2) {
  time_unit_t t;
  if (t1.tv_usec < t2.tv_usec) {
    t.tv_sec = t1.tv_sec - t2.tv_sec - 1;
    t.tv_usec = (1000000 - t2.tv_usec) + t1.tv_usec;
  }
  else {
    t.tv_sec = t1.tv_sec - t2.tv_sec;
    t.tv_usec = t1.tv_usec - t2.tv_usec;
  }
  return t;
}

inline long
time_unit_to_usecs(time_unit_t t) {
  return t.tv_sec*1000000 + t.tv_usec;
}

inline long
time_unit_to_msecs(time_unit_t t) {
  return time_unit_to_usecs(t)/1000;
}

#else

#error "No way to measure time!!!"

typedef long time_unit_t;

inline time_unit_t
get_time() {
  return 0;
}

inline long
time_unit_to_usecs(time_unit_t t) {
  return t;
}

inline long
time_unit_to_msecs(time_unit_t t) {
  return t;
}

#endif


} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_timings_hh)
