/* Interval class declaration.
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

#ifndef _Interval_defs_hh
#define _Interval_defs_hh 1

#ifdef Integer_INTERFACE
#include "Integer.types.hh"
#endif

#include "Boundary.defs.hh"
//#include "Rel.defs.hh"
#include "Bool3.defs.hh"

// Error checking can be turned on or off for that little extra speed lift.

#if !NO_WARNINGS
#define warn_check(cond,message,action) \
  do\
  {\
    if (cond)\
      {\
        cerr << message << endl;\
        action;\
      }\
  } while (0)
#else
#define warn_check(cond,message,action)
#endif


// The Interval Class.
class Interval {
public:
  // This is the type of the flag we use for discriminating
  // between integer and "real" intervals.
  enum Type { REAL=0, INTEGER=1 };

  // Construction Operations
  Interval();
  Interval(Boundary::Value x);
  Interval(int x);
#ifdef Integer_INTERFACE
  Interval(const Integer& x);
#endif
#if !BOUNDARY_TYPE_IS_DOUBLE
  Interval(double x);
#endif

  Interval(const char* x);
  Interval(const LBoundary& lb, const UBoundary& ub, Type flag=REAL);
  void destructure(Boundary::Value& lb_value, bool& lb_closed, 
	      Boundary::Value& ub_value, bool& ub_closed,
	      bool& integer) const;

  Interval& widen_assign (const Interval& old);
  Interval& narrow_assign (const Interval& old);
  Interval& operator |= (const Interval& y);
  Interval& operator &= (const Interval& y);

  Interval& negate ();
  Interval& operator += (const Interval& y);
  Interval& operator -= (const Interval& y);
  Interval& operator *= (const Interval& y);
  Interval& operator /= (const Interval& y);
  Interval& operator += (Boundary::Value y);
  Interval& operator -= (Boundary::Value y);
  Interval& operator *= (Boundary::Value y);
  Interval& operator /= (Boundary::Value y);
  Interval& abs_assign();
  Interval& sqr_assign();
  Interval& sqrt_assign();
  Interval& min_assign(const Interval& y);
  Interval& max_assign(const Interval& y);
  Interval& pow_assign(const Interval& y);
  Interval& div_assign(const Interval& y);
  Interval& mod_assign(const Interval& y);
  Interval& exp_assign();
  Interval& log_assign();
  Interval& floor_assign();
  Interval& ceil_assign();
  Interval& integer_assign();
  Interval& sin_assign();
  Interval& cos_assign();
  Interval& tan_assign();

  bool ok() const;
  bool is_one_point() const;
  bool is_empty() const;
  bool is_full(Type flag=REAL) const;
  bool is_integer() const;
  bool is_lower_min() const;
  bool is_upper_max() const;

  void set_empty();
  void set_lower_min();
  void set_upper_max();
  void set_full(Type flag=REAL);
  bool is_strong() const;
  Boundary::Value value() const;

#if 0
  friend void refine(Interval& x, Rel rel, Interval& y);
  friend void refine1(Interval& x, Rel rel, const Interval& y);
  friend void refine2(const Interval& x, Rel rel, Interval& y);
  friend Rel findrel (const Interval& x, const Interval& y);
  friend void linear_refine(const Interval& k, Rel r, int sz,
			    Boundary::Value* coeffs, Interval** vars);
#endif

  friend bool subset(const Interval& x, const Interval& y);
  friend bool superset(const Interval& x, const Interval& y);
  friend bool subset_eq(const Interval& x, const Interval& y);
  friend bool superset_eq(const Interval& x, const Interval& y);
  friend bool operator == (const Interval& x, const Interval& y);
  friend bool operator != (const Interval& x, const Interval& y);
  friend bool operator && (const Interval& x, const Interval& y);
  friend bool operator && (const Interval& x, Boundary::Value y);
  friend bool lt (const Interval& x, const Interval& y);
  friend bool gt (const Interval& x, const Interval& y);
  friend bool le (const Interval& x, const Interval& y);
  friend bool ge (const Interval& x, const Interval& y);
  friend bool eq (const Interval& x, const Interval& y);
  friend bool ne (const Interval& x, const Interval& y);
  friend bool lt (const Interval& x, Boundary::Value y);
  friend bool gt (const Interval& x, Boundary::Value y);
  friend bool le (const Interval& x, Boundary::Value y);
  friend bool ge (const Interval& x, Boundary::Value y);
  friend bool eq (const Interval& x, Boundary::Value y);
  friend bool ne (const Interval& x, Boundary::Value y);
  friend int lcompare (const Interval& x, const Interval& y);
  friend ostream& operator << (ostream& s, const Interval& x);
  friend istream& operator >> (istream& s, Interval& x);

private:
  LBoundary lower;
  UBoundary upper;
  Type flag;

  void set_maybe_real();
  void set_integer();
  bool contains_multiple(Boundary::Value offset,
			 Boundary::Value period);
};

Bool3 lt3(const Interval& x, const Interval& y);
Bool3 gt3(const Interval& x, const Interval& y);
Bool3 le3(const Interval& x, const Interval& y);
Bool3 ge3(const Interval& x, const Interval& y);
Bool3 eq3(const Interval& x, const Interval& y);
Bool3 ne3(const Interval& x, const Interval& y);
Bool3 lt3(const Interval& x, Boundary::Value y);
Bool3 gt3(const Interval& x, Boundary::Value y);
Bool3 le3(const Interval& x, Boundary::Value y);
Bool3 ge3(const Interval& x, Boundary::Value y);
Bool3 eq3(const Interval& x, Boundary::Value y);
Bool3 ne3(const Interval& x, Boundary::Value y);
#if 0
Bool3 cmp3(const Interval& x, Rel r, const Interval& y);
Bool3 cmp3(const Interval& x, Rel r, Boundary::Value y);
#endif

Interval operator - (const Interval& x);
Interval abs(const Interval& x);
Interval sqr(const Interval& x);
Interval sqrt(const Interval& x);
Interval exp(const Interval& x);
Interval log(const Interval& x);
Interval floor(const Interval& x);
Interval ceil(const Interval& x);
Interval integer(const Interval& x);
Interval sin(const Interval& x);
Interval cos(const Interval& x);
Interval tan(const Interval& x);
Interval widen(const Interval& x, const Interval& y);
Interval narrow(const Interval& x, const Interval& y);
Interval operator & (const Interval& x, const Interval& y);
Interval operator | (const Interval& x, const Interval& y);
Interval operator + (const Interval& x, const Interval& y);
Interval operator - (const Interval& x, const Interval& y);
Interval operator * (const Interval& x, const Interval& y);
Interval operator / (const Interval& x, const Interval& y);
Interval min(const Interval& x, const Interval& y);
Interval max(const Interval& x, const Interval& y);
Interval pow(const Interval& x, const Interval& y);
Interval div(const Interval& x, const Interval& y);
Interval mod(const Interval& x, const Interval& y);

#if !OUTLINE
#include "Interval.inlines.hh"
#endif

#endif
