/*

 File: cimath.h, 2002/03/04

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

#ifndef _CIMATH_H_INCLUDED
#define _CIMATH_H_INCLUDED

#include "cinterval.hpp"

cinterval exp(const cinterval&);
cinterval cos(const cinterval&);
cinterval sin(const cinterval&);
cinterval cosh(const cinterval&);
cinterval sinh(const cinterval&);

cinterval tan(const cinterval&);
cinterval cot(const cinterval&);
cinterval tanh(const cinterval&);
cinterval coth(const cinterval&);

Interval arg(const cinterval&);

cinterval ln(const cinterval&);

cinterval sqr(const cinterval&);
cinterval sqrt(const cinterval&);
cinterval sqrt(const cinterval&,int);
cinterval asin(const cinterval&);
cinterval acos(const cinterval&);
cinterval asinh(const cinterval&);
cinterval acosh(const cinterval&);
cinterval atan(const cinterval&);
cinterval acot(const cinterval&);
cinterval atanh(const cinterval&);
cinterval acoth(const cinterval&);

cinterval power(const cinterval&,int);
cinterval pow(const cinterval&,const Interval&);
cinterval pow(const cinterval&,const cinterval&);

#endif

/*

  End of File: cimath.h

*/
