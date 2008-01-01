/* Performance monitoring using the TSC register.
   Copyright (C) 2001-2008 Roberto Bagnara <bagnara@cs.unipr.it>

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

#ifndef PURRS_tsc_hh
#define PURRS_tsc_hh 1

#include <sys/time.h>
#include <unistd.h>
#include <iostream>

namespace Parma_Recurrence_Relation_Solver {

typedef long long tsc_t;

inline tsc_t
read_tsc() {
  tsc_t tsc;
  asm __volatile__("rdtsc" : "=A" (tsc));
  return tsc;
}

namespace {

double cpu_clock = 0.0;
double tsc_per_msec = 0.0;
double tsc_per_usec = 0.0;

} // namespace

inline tsc_t
msecs_to_tsc(unsigned int msecs) {
  return (tsc_t) (msecs * tsc_per_msec);
}

inline tsc_t
usecs_to_tsc(unsigned int usecs) {
  return (tsc_t) (usecs * tsc_per_usec);
}

inline long
tsc_to_usecs(tsc_t tsc) {
  return (long) (tsc / tsc_per_usec);
}

inline long
tsc_to_msecs(tsc_t tsc) {
  return (long) (tsc / tsc_per_msec);
}

namespace {

#ifdef DJGPP
typedef uclock_t Time_t;
# define gettime(v) v = uclock()
# define time_diff(a, b) ((double) ((a) - (b)) / UCLOCKS_PER_SEC)
#else
typedef struct timeval Time_t;
# define gettime(v) gettimeofday(&v, 0)
# define time_diff(a, b) ((a).tv_sec - (b).tv_sec + ((a).tv_usec - (b).tv_usec) / 1e6)
#endif

inline void
delay(unsigned int msec) {
  usleep(msec * 1000);
}

double
detect_cpu_clock() {
  Time_t tm_begin;
  // Warmup.
  gettime(tm_begin);

  tsc_t tsc_begin = read_tsc();
  gettime(tm_begin);

  delay(1000);

  tsc_t tsc_end = read_tsc(); 
  Time_t tm_end;
  gettime(tm_end);

  return (tsc_end - tsc_begin) / time_diff(tm_end, tm_begin);
}

void
cpu_clock_init() {
  cpu_clock = detect_cpu_clock();
  tsc_per_msec = cpu_clock / 1e3;
  tsc_per_usec = cpu_clock / 1e6;
#if 0
  std::cerr << "Detected CPU clock: " << cpu_clock << " Hz\n" << std::endl;
#endif
}

class Initializer {
public:
  Initializer() {
    cpu_clock_init();
  }
};

Initializer dummy;

} // namespace

} // namespace Parma_Recurrence_Relation_Solver

#endif // !defined(PURRS_tsc_hh)
