/*

 File: cinterval.h, 2002/03/21

 CoStLy (COmplex interval STandard functions LibrarY), Version 0.2

 Copyright (C) Markus Neher, markus.neher@math.uni-karlsruhe.de
               Ingo Eble,    ingoeble@web.de

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef CINTERVAL_H_INCLUDED
#define CINTERVAL_H_INCLUDED

#include "Interval.h" //filib++ Header: Macro Version
#include <iostream>
#ifdef HAS_Complex
#include <Complex.h>
#endif

#ifdef FILIB_NAMESPACES
typedef filib::Interval Interval;
#endif

class cinterval
{
private:

  Interval real_part,imag_part;

public:

  cinterval() {}
  cinterval(const Interval& x, const Interval& y) : real_part(x), imag_part(y) {}
  cinterval(const cinterval& z) : real_part(z.real_part), imag_part(z.imag_part) {}

  explicit cinterval(const Interval& x) : real_part(x), imag_part(0.0) {}
  explicit cinterval(const double& d) : real_part(d), imag_part(0.0) {}
#ifdef HAS_Complex
  explicit cinterval(const Complex& c) : real_part(c.real()), imag_part(c.imag()) {}
#endif

  cinterval& operator = (const cinterval& z)
  {
    if( this == &z ) return *this;

    real_part = z.real_part;
    imag_part = z.imag_part;

    return *this;
  }

  Interval& re() { return real_part; }  
  Interval& im() { return imag_part; }
  const Interval& re() const { return real_part; }  
  const Interval& im() const { return imag_part; }

  cinterval& operator += (const cinterval& z)
  {
    real_part += z.real_part;
    imag_part += z.imag_part;
    return *this;
  }
  
  cinterval& operator -= (const cinterval& z)
  {
    real_part -= z.real_part;
    imag_part -= z.imag_part;
    return *this;
  }

  cinterval& operator *= (const cinterval& z)
  {
    Interval tmp( real_part );
    real_part = real_part*z.real_part - imag_part*z.imag_part;
    imag_part = tmp      *z.imag_part + imag_part*z.real_part;
    return *this;
  }

  cinterval& operator /= (const cinterval&);
};

inline Interval& Re(cinterval& z) { return z.re(); }
inline Interval& Im(cinterval& z) { return z.im(); }
inline const Interval& Re(const cinterval& z) { return z.re(); }
inline const Interval& Im(const cinterval& z) { return z.im(); }

Interval abs (const cinterval&);
#ifdef HAS_Complex
Complex  mid (const cinterval&);
Complex  diam(const cinterval&);
#endif
cinterval operator - (const cinterval&);

// cinterval o cinterval, o \in { +,-,* }.

cinterval operator + (const cinterval&, const cinterval&);
cinterval operator - (const cinterval&, const cinterval&);
cinterval operator * (const cinterval&, const cinterval&);
cinterval operator / (const cinterval&, const cinterval&);

// Interval o cinterval, o \in { +,-,* }.

cinterval operator + (const Interval&, const cinterval&);
cinterval operator - (const Interval&, const cinterval&);
cinterval operator * (const Interval&, const cinterval&);
cinterval operator / (const Interval&, const cinterval&);

// cinterval o Interval, o \in { +,-,* }.

cinterval operator + (const cinterval&, const Interval&);
cinterval operator - (const cinterval&, const Interval&);
cinterval operator * (const cinterval&, const Interval&);
cinterval operator / (const cinterval&, const Interval&);

// double o cinterval, o \in { +,-,* }.

cinterval operator + (const double&, const cinterval&);
cinterval operator - (const double&, const cinterval&);
cinterval operator * (const double&, const cinterval&);
cinterval operator / (const double&, const cinterval&);

// cinterval o double, o \in { +,-,* }.

cinterval operator + (const cinterval&, const double&);
cinterval operator - (const cinterval&, const double&);
cinterval operator * (const cinterval&, const double&);
cinterval operator / (const cinterval&, const double&);

// output

std::ostream& operator << (std::ostream&, const cinterval&);
std::string&  operator << (std::string&, const cinterval&);

//

bool operator <= (const double&, const cinterval&);

//

cinterval operator & (const cinterval&, const cinterval&);

#endif

/*

  End of File: cinterval.h

*/
