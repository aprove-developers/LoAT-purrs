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
#include "Expr_List.types.hh"
#include "Number.types.hh"
#include "Constant.types.hh"
#include "Complex_Interval.types.hh"
#include "Symbol.defs.hh"

#include <set>
#include <ginac/ginac.h>

namespace Parma_Recurrence_Relation_Solver {

//! Type of a functor object.
typedef unsigned Functor;

//! Output operator.
/*! \relates Expr */
std::ostream& operator<<(std::ostream& s, const Expr& x);

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
  \exception std::logic_error thrown if \f$ y = 0 \f$.
*/
Expr operator/(const Expr& x, const Expr& y);

//! Assigns \f$ x + y \f$ to \f$ x \f$ and returns the result.
Expr& operator+=(Expr& x, const Expr& y);

//! Assigns \f$ x - y \f$ to \f$ x \f$ and returns the result.
Expr& operator-=(Expr& x, const Expr& y);

//! Assigns \f$ x \cdot y \f$ to \f$ x \f$ and returns the result.
Expr& operator*=(Expr& x, const Expr& y);

//! \brief
//! If \f$ y \neq 0 \f$, assigns \f$ x / y \f$ to \f$ x \f$
//! and returns the result.
/*!
  \exception std::logic_error thrown if \f$ y = 0 \f$.
*/
Expr& operator/=(Expr& x, const Expr& y);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! equal to \p s.
bool operator==(const Expr& e, const Symbol& s);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! different from \p s.
bool operator!=(const Expr& e, const Symbol& s);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! equal to \p s.
bool operator==(const Symbol& s, const Expr& e);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! different from \p s.
bool operator!=(const Symbol& s, const Expr& e);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! equal to \p c.
bool operator==(const Expr& e, const Constant& c);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! different from \p c.
bool operator!=(const Expr& e, const Constant& c);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! equal to \p c.
bool operator==(const Constant& c, const Expr& e);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically
//! different from \p c.
bool operator!=(const Constant& c, const Expr& e);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically \p n.
bool operator==(const Expr& e, const Number& n);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically \p n.
bool operator==(const Number& n, const Expr& e);

//! Returns <CODE>true</CODE> if and only if \p e is syntactically \p i.
bool operator==(const Expr& e, long i);

//! \brief
//! Returns <CODE>true</CODE> if and only if \p x and \p y are
//! syntactically equal.
bool operator==(const Expr& x, const Expr& y);

//! \brief
//! Returns <CODE>true</CODE> if and only if \p x and \p y are
//! syntactically different.
bool operator!=(const Expr& x, const Expr& y);

//! Compares two numeric expressions.
/*! \relates Expr
  Compares the numeric expressions \p x and \p y and returns the following
  values:
  - -1, if \p x is smaller than \p y;
  - 0,  if x and y are equal;
  - +1, if \p y is smaller than \p x;
  - +2, if the system is not able to compare \p x and \p y.
*/
int compare(const Expr& x, const Expr& y);

//! \brief
//! Returns the application of functor \f$ f \f$ to operand \f$ x \f$,
//! i.e., the expression \f$ f(x) \f$.
Expr apply(Functor f, const Expr& x);

//! \brief
//! Returns the application of functor \f$ f \f$ to operands \f$ x_1 \f$
//! and \f$ x_2 \f$, i.e., the expression \f$ f(x_1, x_2) \f$.
Expr apply(Functor f, const Expr& x1, const Expr& x2);

//! \brief
//! Returns the application of functor \f$ f \f$ to operands \f$ x_1, x_2 \f$
//! and \f$ x_3 \f$, i.e., the expression \f$ f(x_1, x_2, x_3) \f$.
Expr apply(Functor f, const Expr& x1, const Expr& x2, const Expr& x3);

//! \brief
//! Returns the application of functor \f$ f \f$ to operands
//! \f$ x_1, x_2, x_3 \f$ and \f$ x_4 \f$, i.e., the expression
//! \f$ f(x_1, x_2, x_3, x_4) \f$.
Expr apply(Functor f,
	   const Expr& x1, const Expr& x2, const Expr& x3, const Expr& x4);

//! \brief
//! Returns the application of functor \f$ f \f$ to operands
//! \f$ x_1, \dots, x_k \f$, i.e., the expression \f$ f(x_1, \dots, x_k) \f$.
Expr apply(Functor f, const std::vector<Expr>& x);

//! \brief
//! If \f$ x \f$ and \f$ y \f$ are not zero or \f$ x = 0 \f$ and \f$ y \f$
//! is a positive rational number, returns \f$ x^y \f$.
/*!
  \exception std::logic_error thrown if \f$ x = y = 0 \f$.
  \exception std::logic_error thrown if \f$ x = 0 \f$ and \f$ y \f$
                              is not a positive rational number.
*/
Expr pwr(const Expr& x, const Expr& y);

//! \brief
//! If \f$ x \f$ is an exact number, returns the number \f$ y \f$ such that
//! \f$ y^2 = x \f$; otherwise, only for the output, returns the expression
//! \f$ sqrt(x) \f$ not valued. 
Expr sqrt(const Expr& x);

//! 
Expr sin(const Expr& x);
Expr cos(const Expr& x);
Expr acos(const Expr& x);
Expr tan(const Expr& x);
Expr exp(const Expr& x);
Expr log(const Expr& x);
Expr factorial(const Expr& x);
Expr gamma(const Expr& x);

// FIXME: what is a polynomial in a variable `x'?
// The answer is necessary for the function quo(), rem(), prem(), gcd(),
// lcm(), degree(), ldegree(), coeff(), lcoeff(), tcoeff(), expand(),
// collect(), primpart() and content(). 

//! \brief
//! Computes the quotient \f$ q(x) \f$ of polynomials \f$ a(x) \f$ and
//! \f$ b(x) \f$ in \f$ Q[x] \f$.
/*!
  The quozient \f$ q(x) \f$ satisfies \f$ a(x) = b(x) * q(x) + r(x) \f$,
  where \f$ r(x) \f$ is the remainder.

  \exception std::overflow_error   thrown if \f$ b = 0 \f$.
  \exception std::invalid_argument thrown if \f$ a \f$ or \f$ b \f$ are not
                                   rational polynomials.
*/
Expr quo(const Expr& a, const Expr& b, const Symbol& x);

//! \brief
//! Computes the remainder \f$ r(x) \f$ of polynomials \f$ a(x) \f$ and
//! \f$ b(x) \f$ in \f$ Q[x] \f$.
/*
  The remainder \f$ r(x) \f$ satisfies \f$ a(x) = b(x) * q(x) + r(x) \f$,
  where \f$ q(x) \f$ is the quozient.

  \exception std::overflow_error   thrown if \f$ b = 0 \f$.
  \exception std::invalid_argument thrown if \f$ a \f$ or \f$ b \f$ are not
                                   rational polynomials.
*/
Expr rem(const Expr& a, const Expr& b, const Symbol& x);
  
//! \brief
//! Computes the pseudo-remainder \f$ pseudo\_r(x) \f$ of polynomials
//! \f$ a(x) \f$ and \f$ b(x) \f$ in \f$ Z[x] \f$.
/*!
  The pseudo-remainder satisfies
  \f$ c * a(x) = b(x) * pseudo\_q(x) + pseudo\_r(x) \f$, where
  \f$ c = b\_lcoeff^(deg\_a - deg\_b + 1) \f$ with \f$ b\_lcoeff \f$ the
  leading coefficient of \f$ b(x) \f$ and \f$ deg\_a \f$ and \f$ deg\_b \f$
  the degree of \f$ a(x) \f$ and \f$ b(x) \f$, respectively.

  \exception std::overflow_error   thrown if \f$ b = 0 \f$.
  \exception std::invalid_argument thrown if \f$ a \f$ or \f$ b \f$ are not
                                   rational polynomials.
*/
// FIXME: what is the pseudo-quozient?
Expr prem(const Expr& a, const Expr& b, const Symbol& x);

//! \brief
//! Computes the GCD (greatest common divisor) of multivariate polynomials
//! \f$ a(X) \f$ and \f$ b(X) \f$ in \f$ Z[X] \f$.
/*!
  \exception std::invalid_argument thrown if \f$ a \f$ or \f$ b \f$ are not
                                   rational polynomials.
*/
// FIXME: what is a multivariate polynomial?
Expr gcd(const Expr& a, const Expr& b);

//! \brief
//! Computes the LCM (least common multiple) of multivariate polynomials
//! \f$ a(X) \f$ and \f$ b(X) \f$ in \f$ Z[X] \f$.
/*!
  \exception std::invalid_argument thrown if \f$ a \f$ or \f$ b \f$ are not
                                   rational polynomials.
*/
// FIXME: what is a multivariate polynomial?
Expr lcm(const Expr& a, const Expr& b);

//! \brief
//! Computes a square-free factorization for the polynomial \p x with respect
//! to the variable in the list \p y.
/*
  A polynomial \f$ p(X) \in C[X] \f$ is said <EM>square-free</EM>
  if, whenever any two polynomials \f$ q(X) \f$ and \f$ r(X) \f$
  are such that
  \f[
    p(X) = q(X)^2 r(X),
  \f]
  we have \f$ q(X) \in C \f$.
  This means that \f$ p(X) \f$ has no repeated factors, apart
  eventually from constants.
  Given a polynomial \f$ p(X) \in C[X] \f$, we say that the
  decomposition
  \f[
    p(X) = b \cdot p_1(X)^{a_1} \cdot p_2(X)^{a_2} \cdots p_r(X)^{a_r}
  \f]
  is a <EM>square-free factorization</EM> of \f$ p(X) \f$ if the
  following conditions hold:
  - \f$ b \in C \f$ and \f$ b \neq 0 \f$;
  - \f$ a_i \f$ is a positive integer for \f$ i = 1, \ldots, r \f$;
  - the degree of the polynomial \f$ p_i \f$ is strictly positive
    for \f$ i = 1, \ldots, r \f$;
  - the polynomial \f$ \Pi_{i=1}^r p_i(X) \f$ is square-free.

  Square-free factorizations need not be unique.  For example, if
  \f$ a_i \f$ is even, we could change the polynomial \f$ p_i(X) \f$
  into \f$ -p_i(X) \f$.
  Observe also that the factors \f$ p_i(X) \f$ need not be irreducible
  polynomials.
*/
Expr sqrfree(const Expr& x, const Expr_List& y);

//! Solves a linear system of equations.
/*!
  \p x is a list of of equations in the form of relational expressions and
  \p y is a list of symbols.

  \exception std::invalid_argument thrown if \p x and \p y are not lists and,
                                   in particular, lst of equations or lst of
				   symbols, respectively.
  \exception std::logic_error      thrown if the system is not linear. 
*/
Expr lsolve(const Expr_List& x, const Expr_List& y);

//! Returns the function \f$ x(y) \f$.
Expr x(const Expr& y);

//! Returns the function \f$ sum(index,lower,upper,summand) \f$.
Expr sum(const Expr& index, const Expr& lower, const Expr& upper,
	 const Expr& summand);

//! Returns the function \f$ prod(index,lower,upper,factor) \f$.
Expr prod(const Expr& index, const Expr& lower, const Expr& upper,
	  const Expr& factor);

//! Returns the function \f$ mod(x, y) \f$.
Expr mod(const Expr& x, const Expr& y);

class Expr : private GiNaC::ex {
private:
  typedef GiNaC::ex Base;

public:
  //! Default constructor.
  Expr();

  //! Builds the integer expression \p i.
  Expr(int i);

  //! Builds the integer expression \p i.
  Expr(unsigned int i);

  //! Builds the integer expression \p i.
  Expr(long i);

  //! Builds the integer expression \p i.
  Expr(unsigned long i);

  //! Builds the numeric expression \p y.
  Expr(const Number& y);

  //! Builds the symbolic expression \p y.
  Expr(const Symbol& y);

  //! Builds the constant expression \p y.
  Expr(const Constant& y);

  //! Builds the complex interval expression \p y.
  Expr(const Complex_Interval& y);

  //! Builds the expression from a string \p s and a list of symbol \p y.
  Expr(const std::string& s, const Expr_List& y);

  //! Builds the relational expression \f$ x == y \f$.
  // FIXME: temporary
  Expr(const Expr& x, const Expr& y);

  //! Copy-constructor.
  Expr(const Expr& y);

  //! Destructor.
  ~Expr();

  //! Assignment operator.
  Expr& operator=(const Expr& y);

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a symbolic
  //! expression.
  bool is_a_symbol() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a symbolic
  //! expression, in which case its value is assigned to \p s.
  //! Returns <CODE>false</CODE> leaving \p n untouched otherwise.
  bool is_a_symbol(Symbol& s) const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a numeric
  //! expression.
  bool is_a_number() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a numeric
  //! expression, in which case its value is assigned to \p n.
  //! Returns <CODE>false</CODE> leaving \p n untouched otherwise.
  bool is_a_number(Number& n) const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a constant
  //! expression.
  bool is_a_constant() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a complex
  //! interval expression.
  bool is_a_complex_interval() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is an addition
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
  //! Returns <CODE>true</CODE> if and only if \p *this is a constant
  //! power with respect to the variable \p x.
  bool is_a_constant_power(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a function
  //! with like argument an expression.
  bool is_a_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is a constant
  //! function with respect to the variable \p x.
  bool is_a_constant_function(const Symbol& x) const;

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

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>abs()</CODE>.
  bool is_the_abs_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>exp()</CODE>.
  bool is_the_exp_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>log()</CODE>.
  bool is_the_log_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>sin()</CODE>.
  bool is_the_sin_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>cos()</CODE>.
  bool is_the_cos_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>tan()</CODE>.
  bool is_the_tan_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>acos()</CODE>.
  bool is_the_acos_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>factorial()</CODE>.
  bool is_the_factorial_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>gamma()</CODE>.
  /*!
    Note that the method <CODE>tgamma</CODE> of <CODE>GiNaC</CODE> does
    not work on number, it aborts! When we call this method with a
    <CODE>Number</CODE> then automatically the number is converted
    in <CODE>Expr</CODE> by a constructor of the class <CODE>Expr</CODE>.
  */
  bool is_the_gamma_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>x()</CODE>.
  bool is_the_x_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>sum()</CODE>.
  bool is_the_sum_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>prod()</CODE>.
  bool is_the_prod_function() const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is the function
  //! <CODE>mod()</CODE>.
  bool is_the_mod_function() const;

  //! If \p *this is a number, returns the corresponding value.
  Number ex_to_number() const;

  //! If \p *this is symbolic, returns the corresponding symbol.
  Symbol ex_to_symbol() const;

  //! \brief
  //! If \p *this is an interval expression,
  //! returns the corresponding complex interval.
  Complex_Interval ex_to_complex_interval() const;

  // FIXME: info, temporary
  bool is_relation_equal() const;

  //! Returns the functor of \p *this.
  /*!
    It is assumed that \p *this is a function application;
    the behavior is undefined otherwise.
  */
  Functor functor() const;

  //! Returns the operand's number of \p *this.
  /*!
    If \p *this is an addiction or a multiplication of expressions
    returns the terms'number or the factors'number, respectively.
    If \p *this is a power returns \f$ 2 \f$.
    If \p *this is a function returns \f$ 1 \f$.
    Returns \f$ 0 \f$ otherwise.
  */
  unsigned nops() const;

  //! \brief
  //! Returns the \f$ i \f$-th (\f$ i = 0, \dotsc, nops()-1 \f$) operand
  //! of \p *this.
  /*!
    If \p *this is an addiction or a multiplication of expressions
    returns \f$ i \f$-th (\f$ i = 0, \dotsc, nops()-1 \f$) term or factor,
    respectively.
    If \p *this is a power then <CODE>op(0)</CODE> and <CODE>op(1)</CODE>
    return base and exponent of the power.
    If \p *this is a function then <CODE>op(0)</CODE>
    returns the function's argument.
    
    \exception std::out_of_range thrown if
                                 \f$ i \notin \{0, \dotsc, nops() - 1 \} \f$.
  */
  Expr op(unsigned i) const;

#ifdef UNSAFE_ARG
  const Expr& arg(unsigned i) const;
  Expr& arg(unsigned i);
#else
  Expr arg(unsigned i) const;
#endif

  //! \brief
  //! Returns <CODE>true</CODE> if and only if \p *this is sinctatically
  //! zero.
  bool is_zero() const;

  //! \brief
  //! Substitutes every occurrence of the expression \p s in \p * this
  //! with the expression \p r.
  /*! relates Expr */
  Expr substitute(const Expr& s, const Expr& r) const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if a subexpression of \p *this
  //! matches \p x.
  bool has(const Expr& x) const;

  //! Expands \p *this, i.e., distributes multiplication over addition.
  Expr distribute_mul_over_add() const;

  //! \brief
  //! Expands further \p *this, i.e., distributes multiplication over
  //! addition and executes the powers with exponents positive integers
  //! numbers.
  Expr expand() const;

  //! \brief
  //! Returns <CODE>true</CODE> if the expression \p *this is expanded;
  //! returns <CODE>false</CODE> otherwise.
  bool is_expanded() const;

  //! \brief
  //! Sorts expanded expression \p *this collecting, or, in other words,
  //! putting in evidence, the expressions in \p x.
  Expr collect(const Expr_List& x) const;

  //! \brief
  //! If \p *this is a polynomial or a rational function, returns the degree
  //! of the polynomial or the asymptotic degree of the rational function,
  //! respectively; otherwise the behavior of <CODE>degree()</CODE>
  //! is undefined. 
  //! The behavior of <CODE>degree()</CODE> is undefined also if \p x is not
  //! contained in \p *this.
  /*!
    For degree or asymptotic degree we intend the highest integer exponent
    of the powers in \p *this with base equal to \p x.
    This method only works reliably if the input expression is collected
    in terms of \p x.

    \exception std::out_of_range thrown if there is at least one exponent of
                                 the powers with base equal to \p x not
				 integer.
  */
  unsigned degree(const Symbol& x) const;

  //! \brief
  //! If \p *this is a polynomial or a rational function, returns the degree
  //! of the polynomial or the asymptotic degree of the rational function,
  //! respectively; otherwise the behavior of <CODE>degree()</CODE>
  //! is undefined.
  //! The behavior of <CODE>ldegree()</CODE> is undefined also if \p x is not
  //! contained in \p *this.
  /*!
    For low degree or asymptotic low degree we intend the lowest integer
    exponent of the powers in \p *this with base equal to \p x.
    This method only works reliably if the input expression is collected
    in terms of \p x.

    \exception std::out_of_range thrown if there is at least one exponent of
                                 the powers with base equal to \p x not
				 integer.
  */
  unsigned ldegree(const Symbol& x) const;

  //! \brief
  //! Returns coefficient of the term of degree equal to \p k, i. e.,
  //! the coefficient of the power in the expanded polynomial \p *this
  //! with base equal to \p x and degree equal to \p k.
  //! The behavior of <CODE>coeff()</CODE> is undefined if \p *this
  //! is not a polynomial or \p x is not contained in \p *this or there
  //! is not a term of \p k's degree.
  // FIXME: lavora anche sulle espressioni non polinomiali!!
  Expr coeff(const Symbol& x, int k) const;

  //! \brief
  //! Returns coefficient of the term with highest degree, i. e., the
  //! coefficient of the highest power in the expanded polynomial \p *this
  //! with base equal to \p x.
  //! It is equivalent to <CODE>coeff(x, degree(x))</CODE>.
  //! The behavior of <CODE>lcoeff()</CODE> is undefined if \p *this
  //! is not a polynomial or \p x is not contained in \p *this or there
  //! is not a term of \p k's degree.
  // FIXME: lavora anche sulle espressioni non polinomiali!!
  Expr lcoeff(const Symbol& x) const;

  //! \brief
  //! Returns coefficient of the term with lowest degree, i. e., the
  //! coefficient of the lowest power in the expanded polynomial \p *this
  //! with base equal to \p x.
  //! It is equivalent to <CODE>coeff(x, ldegree(x))</CODE>.
  //! The behavior of <CODE>tcoeff()</CODE> is undefined if \p *this
  //! is not a polynomial or \p x is not contained in \p *this or there
  //! is not a term of \p k's degree.
  // FIXME: lavora anche sulle espressioni non polinomiali!!
  Expr tcoeff(const Symbol& x) const;

  //! \brief
  //! Computes the primitive part of a multivariate polynomial in
  //! \f$ Z[x] \f$ with respect to \p x.
  /*!
    By definition, <EM>primitive part</EM> of polynomial \f$ f \f$ is
    a polynomial \f$ g \f$ whose coefficients are relatively prime such
    that \f$ f = r * g \f$ for some element \f$ r \f$ of the coefficient
    ring.
    The product of content part and primitive part is the polynomial itself.

    \exception std::out_of_range     thrown if there is at least one exponent
                                     of the powers with base equal to \p x not
				     integer.
    \exception std::domain_error     thrown if VEDI TEST primpart.cc   
    \exception std::invalid_argument thrown if VEDI TEST primpart.cc
  */
  Expr primpart(const Symbol& x) const;

  //! \brief
  //! Computes content part of a multivariate polynomial in \f$ Z[x]
  //! \f$ with respect to \p x.
  /*!
    By definition, <EM>content part</EM> of polynomial is the greatest
    common divisor of its coefficients.
    The product of content part and primitive part is the polynomial itself.
    FIXME: rivedere eccezioni.
  */
  Expr content(const Symbol& x) const;

  //! Returns numerator \p x and denominator \p y of \p *this.
  /*!
    If the expression is not of the normal form `numerator/denominator'
    (where numerator and denominator are relatively prime) polynomials,
    it is first converted to this form and then the numerator and
    denominator are returned.
  */
  // FIXME: what is a polynomial?
  void numerator_denominator(Expr& x, Expr& y) const;

  //! Returns left hand side of relational expression \p *this.
  /*!
    \exception std::runtime_error thrown if \p *this is not a relational
                                  expression.
  */
  Expr lhs() const;

  //! Returns right hand side of relational expression \p *this.
  /*!
    \exception std::runtime_error thrown if \p *this is not a relational
                                  expression.
  */
  Expr rhs() const;

  //! \brief
  //! Returns partial derivative of \p nth order of \p *this with respect
  //! to the variable \p x.
  Expr diff(const Symbol& x, unsigned nth = 1);

  void latex_print(std::ostream& s);

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a scalar rapresentation for
  //! \p x; returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of <EM>scalar representation
    for \f$ x \f$</EM>, where \f$ x \f$ is any symbol.
    This is more briefly written <EM>scalar</EM> and defined inductively
    as follows:
    - every number is a scalar;
    - every symbolic constant is a scalar;
    - every parameter different from \f$ x \f$ is a scalar;
    - if \f$ f \f$ is any function and \f$ x_1, \dots, x_k \f$ are scalars,
      then \f$ f(x_1, \dots, x_k) \f$ is a scalar;
    - if \f$ a \f$ and \f$ b \f$ are scalars then
      \f$ a+b \f$, \f$ a*b \f$, and \f$ a^b \f$ are scalars.
  */
  bool is_scalar_representation(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is an integer scalar
  //! rapresentation for \p x; returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of <EM>integer scalar
    representation for \f$ x \f$</EM>, where \f$ x \f$ is any symbol.
    This is defined inductively as follows:
    - every integer number is a integer scalar representation;
    - every parameter different from \f$ x \f$ is an integer scalar
      representation;
    - if \f$ f \f$ is any function and \f$ x_1, \dots, x_k \f$ are parameters
      different from \f$ x \f$, then \f$ f(x_1, \dots, x_k) \f$ is an
      integer scalar representation (note: is not a constant function
      with respect to the variable \p x, because otherwise also \f$ log(2) \f$
      would be an integerscalar representation);
    - if \f$ a \f$ and \f$ b \f$ are integer scalar representations then
      \f$ a+b \f$, \f$ a*b \f$ are integer scalar representations;
    - if \f$ a \f$ is an integer scalar representation in \f$ x \f$ and
      \f$ b \f$ is a positive integer or a parameter different from
      \f$ x \f$, then \f$ a^b \f$ is an integer scalar representation.
  */
  bool is_integer_scalar_representation(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a rational scalar
  //! rapresentation for \p x; returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of <EM>rational scalar
    representation for \f$ x \f$</EM>, where \f$ x \f$ is any symbol.
    This is defined inductively as follows:
    - every rational number is a rational scalar representation;
    - every parameter different from \f$ x \f$ is a rational scalar
      representation;
    - if \f$ f \f$ is any function and \f$ x_1, \dots, x_k \f$ are parameters
      different from \f$ x \f$, then \f$ f(x_1, \dots, x_k) \f$ is a
      rational scalar representation;
    - if \f$ a \f$ and \f$ b \f$ are rational scalar representations then
      \f$ a+b \f$, \f$ a*b \f$ are rational scalar representations;
    - if \f$ a \f$ is a rational scalar representation in \f$ x \f$ and
      \f$ b \f$ is a positive integer or a parameter different from
      \f$ x \f$, then \f$ a^b \f$ is a rational scalar representation.
  */
  bool is_rational_scalar_representation(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a scalar rapresentation;
  //! returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of <EM>scalar representation</EM>,
    defined inductively as follows:
    - every number is a scalar representation;
    - every symbolic constant is a scalar representation;
    - if \f$ f \f$ is any function and \f$ x_1, \dots, x_k \f$ are scalar
      representations, then \f$ f(x_1, \dots, x_k) \f$ is a scalar
      representation;
    - if \f$ a \f$ and \f$ b \f$ are scalar representations, then
      \f$ a+b \f$, \f$ a*b \f$ and \f$ a^b \f$ are scalars.
  */
  bool is_scalar_representation() const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a polynomial in \p x;
  //! returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of <EM>polynomial in \f$ x \f$</EM>,
    where \f$ x \f$ is any symbol.
    This is more briefly written <EM>polynomial</EM> and defined inductively
    as follows:
    - every scalar representation for \f$ x \f$ is a polynomial;
    - \f$ x \f$ is a polynomial;
    - if \f$ a \f$ is a polynomial in \f$ x \f$ and \f$ b \f$ is a positive
      integer, then \f$ a^b \f$ is a polynomial;
    - if \f$ a \f$ and \f$ b \f$ are polynomials then
      \f$ a + b \f$ and \f$ a * b \f$ are polynomials.
  */
  bool is_polynomial(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is an integer polynomial in \p x;
  //! returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of
    <EM>integer polynomial in \f$ x \f$</EM>, where \f$ x \f$ is any symbol.
    This is more briefly written <EM>integer polynomial</EM> and defined
    inductively as follows:
    - every integer scalar representation for \f$ x \f$ is an integer
      polynomial;
    - \f$ x \f$ is an integer polynomial;
    - if \f$ a \f$ is an integer polynomial in \f$ x \f$ and \f$ b \f$
      is a positive integer, then \f$ a^b \f$ is an integer polynomial;
    - if \f$ a \f$ and \f$ b \f$ are integer polynomials then
      \f$ a + b \f$ and \f$ a * b \f$ are integer polynomials.
  */
  bool is_integer_polynomial(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a rational polynomial in \p x;
  //! returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of
    <EM>rational polynomial in \f$ x \f$</EM>, where \f$ x \f$ is any symbol.
    This is more briefly written <EM>rational polynomial</EM> and defined
    inductively as follows:
    - every rational scalar representation for \f$ x \f$ is a rational
      polynomial;
    - \f$ x \f$ is a rational polynomial;
    - if \f$ a \f$ is a rational polynomial in \f$ x \f$ and \f$ b \f$
      is a positive integer, then \f$ a^b \f$ is a rational polynomial;
    - if \f$ a \f$ and \f$ b \f$ are rational polynomials then
      \f$ a + b \f$ and \f$ a * b \f$ are rational polynomials.
  */
  bool is_rational_polynomial(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a multivariate polynomial;
  //! returns <CODE>false</CODE> otherwise.
  /*!
    This function realizes the definition of
    <EM>multivariate polynomial\f$</EM>, defined inductively as follows:
    - every scalar representation is a multivariate polynomial;
    - every parameter is a multivariate polynomial;
    - if \f$ a \f$ is a multivariate polynomial and \f$ b \f$ is a positive
      integer, then \f$ a^b \f$ is a multivariate polynomial;
    - if \f$ a \f$ and \f$ b \f$ are multivariate polynomials then
      \f$ a + b \f$ and \f$ a * b \f$ are multivariate polynomials.
  */
  bool is_multivariate_polynomial() const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a rational function in \p x;
  //! returns <CODE>false</CODE> otherwise.
  /*!
    A quotient of two polynomials \f$ P(x) \f$ and \f$ Q(x) \f$,
    \f[
      R(x) = \frac{P(x)}{Q(x)},
    \f]
    is called a <EM>rational function</EM>.
  */
  bool is_rational_function(const Symbol& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this contains non rational numbers;
  //! returns <CODE>false</CODE> otherwise.
  bool has_non_rational_numbers() const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this contains the function
  //! \f$ x() \f$ with its argument containing the expression \p e;
  //! returns <CODE>false</CODE> otherwise.
  bool has_x_function(const Expr& y) const;

  //! \brief
  //! Returns <CODE>true</CODE> if and only if the functions \f$ x() \f$
  //! of \p *this are all initial conditions; returns <CODE>false</CODE>
  //! otherwise, i.e., if there is at least one function \f$ x() \f$
  //! that is not an initial condition.
  /*!
    A function \f$ x() \f$ is an initial condition in two cases:
    - the argument is a positive integer;
    - the argument is equal to \f$ mod(n, k) + h \f$, with
      \f$ k, h \in \Nset \f$.
  */
  bool has_x_function_only_ic() const;

  //! \brief
  //! Let \p x be a simple not numeric expression, i.e., it is not a sum or
  //! a product or a power, which is assumed to be already expanded.
  //! Returns a new expression obtained from \p *this collecting the
  //! expression \p x and puts in \p coeff_x the new coefficient of the
  //! expression \p x. If \p x it is not in \p *this, then returns \p *this.
  Expr collect_term(const Expr& x, Expr& coeff_x) const;

  //! Returns \p *this evaluated numerically. 
  Expr unsafe_fp_approximation() const;

  //! Collects all symbols contained in \p *this in two sets dividing
  //! the symbols generated by the system (i.e., symbols with the name
  //! of the shape "symbol" followed by a number) from all the other symbols.
  void collect_symbols(Symbol::SymbolSet& system_generated_symbols,
		       Symbol::SymbolSet& new_symbols) const;

private:
  friend class Number;
  friend class Symbol;
  friend class Constant;
  friend class Expr_List;
  friend class Matrix;

  friend std::ostream& operator<<(std::ostream& s, const Expr& x);

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

  friend bool operator==(const Expr& e, const Symbol& s);
  friend bool operator==(const Expr& e, const Constant& c);
  friend bool operator==(const Expr& e, const Number& n);
  friend bool operator==(const Expr& e, long i);

  friend bool operator==(const Expr& x, const Expr& y);

  friend int compare(const Expr& x, const Expr& y);
  friend Expr apply(Functor f, const Expr& x);
  friend Expr apply(Functor f, const Expr& x1, const Expr& x2);
  friend Expr apply(Functor f, const Expr& x1, const Expr& x2, const Expr& x3);
  friend Expr apply(Functor f,
		    const Expr& x1, const Expr& x2,
		    const Expr& x3, const Expr& x4);
  friend Expr apply(Functor f, const std::vector<Expr>& x);
  friend Expr pwr(const Expr& x, const Expr& y);
  friend Expr sqrt(const Expr& x);
  friend Expr sin(const Expr& x);
  friend Expr cos(const Expr& x);
  friend Expr acos(const Expr& x);
  friend Expr tan(const Expr& x);
  friend Expr exp(const Expr& x);
  friend Expr log(const Expr& x);
  friend Expr factorial(const Expr& x);
  friend Expr gamma(const Expr& x);
  friend Expr quo(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr rem(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr prem(const Expr& a, const Expr& b, const Symbol& x);
  friend Expr gcd(const Expr& a, const Expr& b);
  friend Expr lcm(const Expr& a, const Expr& b);
  friend Expr sqrfree(const Expr& x, const Expr_List& y);
  friend Expr lsolve(const Expr_List& x, const Expr_List& y);
  friend Expr x(const Expr& y);
  friend Expr sum(const Expr& index, const Expr& lower, const Expr& upper,
		  const Expr& summand);
  friend Expr prod(const Expr& index, const Expr& lower, const Expr& upper,
		   const Expr& factor);
  friend Expr mod(const Expr& x, const Expr& y);

  //! Builds the expression corresponding to \p ge.
  Expr(const GiNaC::ex& ge);

  //! Builds the function corresponding to \p gf.
  Expr(const GiNaC::function& gf);
};

/*!
  Assuming <CODE>sh.size_norm(const Symbol& s)</CODE> computes
  the norm of <CODE>s</CODE>, returns the size-norm of \p e.
*/
template <typename SymbolHandler>
unsigned
generic_size_norm(const Expr& e, const SymbolHandler& sh);

} // namespace Parma_Recurrence_Relation_Solver

#include "Expr.inlines.hh"

#endif // !defined(PURRS_Expr_defs_hh)
