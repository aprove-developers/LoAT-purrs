/*

 File: error.h, 2001/12/03

 CoStLy (COmplex interval STandard functions LibrarY), Version 0.2

 Copyright (C) Markus Neher, markus.neher@math.uni-karlsruhe.de
               Ingo Eble,    IngoEble@web.de

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

#ifndef _ERROR_H_INCLUDED
#define _ERROR_H_INCLUDED

#include <stdexcept>

class syntax_error : public logic_error 
{
 public:
  syntax_error(const string& what_arg) : logic_error(what_arg) {}
};

class division_by_zero : public logic_error  
{
 public:
  division_by_zero() : logic_error("") {}
};

class function_not_defined : public invalid_argument
{
 public:
  function_not_defined() : invalid_argument("") {}
};

class wrong_dimensions : public length_error
{
 public:
  wrong_dimensions() : length_error("") {}
};

#endif

/*

  End of File: error.h

*/
