/*

 File: error.h, 2002/05/16

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

#ifndef _ERROR_H_INCLUDED
#define _ERROR_H_INCLUDED

#include <stdexcept>

class syntax_error : public std::logic_error 
{
 public:
  syntax_error(const std::string& what_arg) : std::logic_error(what_arg) { abort(); }
};

class division_by_zero : public std::logic_error  
{
 public:
  division_by_zero() : std::logic_error("") { abort(); }
};

class function_not_defined : public std::invalid_argument
{
 public:
  function_not_defined() : std::invalid_argument("") { abort(); }
};

class wrong_dimensions : public std::length_error
{
 public:
  wrong_dimensions() : std::length_error("") { abort(); }
};

#endif

/*

  End of File: error.h

*/
