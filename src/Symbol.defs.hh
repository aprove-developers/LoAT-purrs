/* Symbol class declaration.
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

#ifndef PURRS_Symbol_defs_hh
#define PURRS_Symbol_defs_hh 1

#include "Symbol.types.hh"
#include "Expr.types.hh"

#include <ginac/ginac.h>
#include <string>

namespace Parma_Recurrence_Relation_Solver {

//! The set of symbols for parameters.
/*!
  ...
*/
class Symbol {
public:
  //! Builds a symbol with a unique name.
  Symbol();

  //! Builds a symbol named \p n.
  Symbol(const char* n);

  //! Builds a symbol named \p str.
  explicit Symbol(const std::string& str);

  //! Copy-constructor.
  Symbol(const Symbol& y);

  //! Destructor.
  ~Symbol();

  //! Assignment operator.
  Symbol& operator=(const Symbol& y);

  //! Return the symbol's name.
  std::string get_name() const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a symbol generated
  //! by the system; returns <CODE>false</CODE> otherwise.
  /*!
    A symbol is generated by the system if its name has the shape "symbol"
    followed by a number.
  */
  bool is_system_generated() const;

#ifdef PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS
 //! Binary predicate defining the total ordering on the names of symbols.
#endif // PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS
  struct NameCompare {
    //! \brief
    //! Returns <CODE>true</CODE> if the name of \p x comes before
    //! (lexicographically) the name of \p y.
    bool operator()(const Symbol& x, const Symbol& y) const;
  };

  typedef std::set<Symbol, Symbol::NameCompare> SymbolSet;

private:
  friend class Expr;
  friend class Expr_List;

  friend Expr quo(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr rem(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr prem(const Expr& a, const Expr& b, const Symbol& x);

  friend bool operator==(const Expr& e, const Symbol& s);

  GiNaC::symbol s;

  //! Builds the symbol corresponding to \p gs.
  Symbol(const GiNaC::symbol& gs);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Symbol.inlines.hh"

#endif // !defined(PURRS_Symbol_defs_hh)
