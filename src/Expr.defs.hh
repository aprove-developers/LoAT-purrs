/* Expr class declaration.
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
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.

For the most up-to-date information see the PURRS site:
http://www.cs.unipr.it/purrs/ . */

#ifndef PURRS_Expr_defs_hh
#define PURRS_Expr_defs_hh 1

#include "Expr.types.hh"
#include "Symbol.types.hh"
#include "Expr_List.types.hh"
#include "Number.types.hh"
#include "Constant.types.hh"

#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

//! Output operator.
/*! \relates Expr */
std::ostream& operator<<(std::ostream& os, const Expr& exp);

//! Returns \f$ x \f$.
Expr operator+(const Expr& x);

//! Returns \f$ - x \f$.
Expr operator-(const Expr& x);

//! Returns \f$ x + y \f$.
Expr operator+(const Expr& x, const Expr& y);

//! Returns \f$ x - y \f$.
Expr operator-(const Expr& x, const Expr& y);

//! Returns \f$ x \cdot y \f$.
Expr operator*(const Expr& x, const Expr& y);

//! If \f$ y \neq 0 \f$, returns \f$ x / y \f$.
/*!
  \exception std::runtime_error thrown if <CODE>y.is_zero()</CODE>.
*/
Expr operator/(const Expr& x, const Expr& y);

//! Assigns \f$ x + y \f$ to \f$ x \f$ and returns the result.
Expr& operator+=(Expr& x, const Expr& y);

//! Assigns \f$ x - y \f$ to \f$ x \f$ and returns the result.
Expr& operator-=(Expr& x, const Expr& y);

//! Assigns \f$ x \cdot y \f$ to \f$ x \f$ and returns the result.
Expr& operator*=(Expr& x, const Expr& y);

//! If \f$ y \neq 0 \f$, assigns \f$ x / y \f$ to \f$ x \f$
//! and returns the result.
/*!
  \exception std::runtime_error thrown if \f$ y = 0 \f$.
*/
Expr& operator/=(Expr& x, const Expr& y);

#if 0
bool operator==(const Expr& x, const Expr& y);
bool operator!=(const Expr& x, const Expr& y);
#endif


// FIXME: meglio con argomento di default `unsigned label = 0'?
Expr wild(unsigned label);

Expr power(const Expr& b, const Expr& e);
Expr sqrt(const Expr& e);
Expr sin(const Expr& e);
Expr cos(const Expr& e);
Expr acos(const Expr& e);
Expr tan(const Expr& e);
Expr exp(const Expr& e);
Expr log(const Expr& e);

Expr quo(const Expr& a, const Expr& b, const Symbol& x);
Expr rem(const Expr& a, const Expr& b, const Symbol& x);
Expr prem(const Expr& a, const Expr& b, const Symbol& x);
Expr gcd(const Expr& a, const Expr& b);
Expr lcm(const Expr& a, const Expr& b);
Expr sqrfree(const Expr& e, const Expr_List& lst);

Expr lsolve(const Expr_List& lst1, const Expr_List& lst2);

Expr x(const Expr& e);

class Expr : private GiNaC::ex {
private:
  typedef GiNaC::ex Base;

public:
  //! Default constructor.
  Expr();

  //! Builds the integer expression \p i.
  Expr(int i);

  //! Builds the numeric expression \p n.
  Expr(const Number& n);

  //! Builds the symbolic expression \p s.
  Expr(const Symbol& s);

  //! Builds the constant expression \p k.
  Expr(const Constant& k);

  //! Builds the expression from a string \p st and a list of symbol \p lst.
  Expr(const std::string& st, const Expr_List& lst);

  //! Builds the relational expression \f$ lh == rh \f$.
  // FIXME: temporary
  Expr(const Expr& lh, const Expr& rh);

  //! Copy-constructor.
  Expr(const Expr& exp);

  //! Destructor.
  ~Expr();

  //! Assignment operator.
  Expr& operator=(const Expr& exp);

  //! \brief
  //! Accedes to \f$ i \f$-th term/factor of an addiction/multiplication of
  //! expressions.
  Expr operator[](int i) const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a symbolic
  //! expression.
  bool is_a_symbol() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a numberic
  //! expression.
  bool is_a_number() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a constant
  //! expression.
  bool is_a_constant() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is an addiction
  //! of expressions.
  bool is_a_add() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a multiplication
  //! of expressions.
  bool is_a_mul() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is an expression
  //! in the form of power.
  bool is_a_power() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a function
  //! with like argument an expression.
  bool is_a_function() const;

  //! Returns <CODE>true</CODE> if and only if \p *this is a matrix expression.
  bool is_a_matrix() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a list of
  //! expressions.
  bool is_a_Expr_List() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is an expression
  //! in the form of relational.
  bool is_a_relational() const;

  bool is_the_abs_function() const;
  bool is_the_exp_function() const;
  bool is_the_log_function() const;
  bool is_the_sin_function() const;
  bool is_the_cos_function() const;
  bool is_the_tan_function() const;
  bool is_the_acos_function() const;

  //! Returns the numeric value of \p *this.
  Number ex_to_number() const;

  // FIXME: info, temporary
  bool is_integer_polynomial() const;
  bool is_rational_polynomial() const;
  bool is_relation_equal() const;

  //! If \p *this is an addiction or a multiplication of expressions
  //! returns the terms'number or the factors'number, respectively.
  //! If \p *this is a power returns \f$ 2 \f$.
  //! If \p *this is a function returns \f$ 1 \f$.
  //! Returns \f$ 0 \f$ otherwise.
  unsigned nops() const;

  //! If \p *this is an addiction or a multiplication of expressions
  //! returns \f$ i \f$-th (\f$ i = 0, \dotsc, nops()-1 \f$) term or factor,
  //! respectively.
  //! If \p *this is a power then <CODE>op(0)</CODE> and <CODE>op(1)</CODE>
  //! return base and exponent of the power.
  //! If \p *this is a function then <CODE>op(0)</CODE>
  //! returns the function's argument.
  Expr op(unsigned i) const;

  //! Returns <CODE>true</CODE> if and only if \p *this is sinctatically
  //! equal to \p e.
  bool is_equal(const Expr& e) const;

  //! Returns <CODE>true</CODE> if and only if \p *this is sinctatically
  //! zero.
  bool is_zero() const;

  //! Substitutes in \p *this the occurrences of \p exp1 with \p exp2.
  Expr subs(const Expr& exp1, const Expr& exp2) const;

  //! Allows the substitution at the same time of the occurences in \p *this
  //! of the expressions contained in \p to_replace with the relative
  //! expressions in \p replacements.
  Expr subs(const Expr_List& to_replace, const Expr_List& replacements) const;

  //! Returns <CODE>true</CODE> if and only if \p *this matches \p pattern.
  bool match(const Expr& pattern) const;

  bool match(const Expr& pattern, Expr_List& replacements) const;
  bool has(const Expr& pattern) const;

  Expr expand() const;
  Expr collect(const Expr_List& lst) const;
  int degree(const Symbol& symb) const;
  int ldegree(const Symbol& symb) const;
  Expr coeff(const Symbol& symb, int k) const;
  Expr lcoeff(const Symbol& symb) const;
  Expr tcoeff(const Symbol& symb) const;
  Expr primpart(const Symbol& symb) const;
  Expr numerator() const;
  Expr denominator() const;
  void numerator_denominator(Expr& numer, Expr& denom) const;
  Expr lhs() const;
  Expr rhs() const;

  Expr diff(const Symbol& symb, unsigned nth = 1);

  Expr to_rational(Expr_List& lst);

private:
  friend std::ostream& operator<<(std::ostream& os, const Expr& exp);

  friend Expr operator+(const Expr& x);
  friend Expr operator-(const Expr& x);
  friend Expr operator+(const Expr& x, const Expr& y);
  friend Expr operator-(const Expr& x, const Expr& y);
  friend Expr operator*(const Expr& x, const Expr& y);
  friend Expr operator/(const Expr& x, const Expr& y);
  friend Expr& operator+=(Expr& x, const Expr& y);
  friend Expr& operator-=(Expr& x, const Expr& y);
  friend Expr& operator*=(Expr& x, const Expr& y);
  friend Expr& operator/=(Expr& x, const Expr& y);

#if 0
  friend Expr& operator==(const Expr& x, const Expr& y);
  friend Expr& operator!=(const Expr& x, const Expr& y);
#endif

  friend Expr wild(unsigned label);

  friend Expr power(const Expr& b, const Expr& e);
  friend Expr sqrt(const Expr& e);
  friend Expr sin(const Expr& e);
  friend Expr cos(const Expr& e);
  friend Expr acos(const Expr& e);
  friend Expr tan(const Expr& e);
  friend Expr exp(const Expr& e);
  friend Expr log(const Expr& e);

  friend Expr quo(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr rem(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr prem(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr gcd(const Expr& a, const Expr& b);
  friend Expr lcm(const Expr& a, const Expr& b);
  friend Expr sqrfree(const Expr& e, const Expr_List& lst);
  
  friend Expr lsolve(const Expr_List& lst1, const Expr_List& lst2);
  
  friend Expr x(const Expr& e);

  friend class Number;
  friend class Symbol;
  friend class Constant;
  friend class Expr_List;
  friend class Matrix;

  Expr(const GiNaC::ex& ge);
  Expr(const GiNaC::function& gf);
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Expr.inlines.hh"

#endif // !defined(PURRS_Expr_defs_hh)
