/*

 File: cimath.h, 2002/12/06

 CoStLy (COmplex interval STandard functions LibrarY), Version 0.3

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

#ifndef _CIMATH_H_INCLUDED
#define _CIMATH_H_INCLUDED

#ifdef FILIB_VERSION

#include "cinterval.h"
typedef Interval interval;

#else //C-XSC Version

#include "cinterval.hpp"
using namespace cxsc;

#endif

#include <list> //STL

cinterval exp(const cinterval&);
cinterval cos(const cinterval&);
cinterval sin(const cinterval&);
cinterval cosh(const cinterval&);
cinterval sinh(const cinterval&);

cinterval tan(const cinterval&);
cinterval cot(const cinterval&);
cinterval tanh(const cinterval&);
cinterval coth(const cinterval&);

interval arg(const cinterval&);

#ifdef FILIB_VERSION
cinterval log(const cinterval&);
#else //C-XSC Version
cinterval ln(const cinterval&);
#endif

cinterval sqr(const cinterval&);

cinterval sqrt(const cinterval&);
cinterval root(const cinterval&,unsigned int);

std::list<cinterval> sqrt_all(const cinterval&); 
std::list<cinterval> root_all(const cinterval&,unsigned int);

cinterval asin(const cinterval&);
cinterval acos(const cinterval&);
cinterval asinh(const cinterval&);
cinterval acosh(const cinterval&);
cinterval atan(const cinterval&);
cinterval acot(const cinterval&);
cinterval atanh(const cinterval&);
cinterval acoth(const cinterval&);

cinterval power(const cinterval&,int);

cinterval pow(const cinterval&,const interval&);
cinterval pow(const cinterval&,const cinterval&);

std::list<cinterval> pow_all(const cinterval&,const interval&);
std::list<cinterval> pow_all(const cinterval&,const cinterval&);

#endif

/*

  End of File: cimath.h

*/
