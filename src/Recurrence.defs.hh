/* This is the main object PURRS operates upon: a recurrence.
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

#ifndef PURRS_Recurrence_defs_hh
#define PURRS_Recurrence_defs_hh 1

#include "Recurrence.types.hh"
#include "Cached_Expr.defs.hh"
#include "Non_Linear_Info.types.hh"
#include "Weighted_Average_Info.types.hh"
#include "Blackboard.defs.hh"
#include "Finite_Order_Info.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Expr.defs.hh"
#include "alg_eq_solver.hh"
#include <iterator>
#include <map>
#include <iosfwd>

namespace Parma_Recurrence_Relation_Solver {

//! \brief
//! An unsigned integral type for representing different types of
//! indices of the recurrence (e.g. the order of the recurrence).
typedef unsigned int index_type;

//! \brief
//! Assuming that \p x and \p y represent sequences of reals,
//! returns <CODE>true</CODE> if it can be proved that all
//! the recurrences represented by \p x are less than all
//! the recurrences represented by \p y.
/*!
  \relates Recurrence
  Let \f$ S_x, S_y \in \wp(\Nset \to \Rset) \f$ be the sets of sequences
  represented by \p x and \p y, respectively.  If the system can prove that
  \f[
    \forall s_x \in S_x \itc \forall s_y \in S_y \itc \forall n \in \Nset
      \itc s_x(n) < s_y(n),
  \f]
  then <CODE>true</CODE> is returned;
  <CODE>false</CODE> is returned otherwise.
*/
bool less_than(const Recurrence& x, const Recurrence& y);

//! The base class for recurrence relations.
/*!
  An object of this class abstracts a (possibly infinite) set of
  sequences of real numbers.  Formally, an object of class
  Recurrence is an element of \f$ \wp(\Nset \to \Rset) \f$.  Elements
  of \f$ \wp(\Nset \to \Rset) \f$ can be specified by imposing a
  recurrence relation (whence the class' name) that the sequences must
  satisy, by imposing some initial conditions for such recurrences,
  and by computing approximations of set union and intersection.

  In this class we will introduce constructors and methods in order to
  create the object Recurrence, to compute the exact solution or the
  approximations, to verify the results and to specify possible
  initial conditions.

  The system works with recurrences in \ref normal_form "normal form"
  \f[
    x(n) = r.
  \f]
  The symbol \f$ n \f$ and the function \f$ x \f$ have a special meaning
  in the context of recurrence relations: \f$ n \f$ is the index of the
  recurrence and \f$ x \f$ indicates the function defined from the same
  recurrence. Consequently \f$ n \f$ and \f$ x \f$ can not be used as
  any other parameter.
  Writing <CODE>Symbol n("n")</CODE>, the symbol created is different
  from the built-in <CODE>Recurrence::n</CODE>, even if the name is the
  same. You can define <CODE>const Symbol& n = Recurrence::n</CODE>
  once for all or you must specify which is the symbol <CODE>n</CODE>
  you want to consider using every time <CODE>Recurrence::n</CODE>.

  We say that a recurrence is
  \anchor syntactically_correct <EM>syntactically correct</EM>
  if its right-hand side is sum, product and powers of the following object:
  - numbers (no floating-point);
  - symbols for parameters;
  - mathematical functions (ex. \f$ \log() \f$, \f$ \dots \f$);
  - unary functions \f$ x() \f$ with argument of the following form:
    \f$ n + k \f$ with \f$ k \in \Zset \f$,
    \f$ n / k \f$ with \f$ k \in \Qset \f$.
  
  The following examples show some different ways in order to build
  objects Recurrence:

  \par Example 1
  \code
    const Symbol& n = Recurrence::n;
    Recurrence rec1;
    rec1 = Recurrence(3/n*x(n-1));
    Recurrence rec2(x(n-1)+4);
    Recurrence rec3(rec1);
  \endcode

  \par Example 2
  The following code builds the parametric recurrence
  \f$ x(n) = a x(n-1)+n \f$:
  \code
    const Symbol& n = Recurrence::n;
    Symbol a("a");
    Recurrence rec(a*x(n-1)+n);
  \endcode

  We wish to remark the possible problem arising from
  an incorrect use of the special symbol \f$ n \f$:

  \par Example 3
  This definition produces a time compilation error because the symbol
  \f$ n \f$ is unknown:
  \code
    Recurrence rec(3*x(n-1)+1);
  \endcode
  To add the symbol's definition
  \code
    Symbol n("n");
  \endcode
  introduces an error more complicated than the previous: the system
  does not consider \f$ n \f$ like the index of the recurrence producing
  an unexpected result.

  To define and to use the symbol \f$ x \f$, as for the symbol \f$ n \f$,
  creates problems because it is hidden to the compiler the meaning that
  this special symbol has in the recurrence.

  If the user inserts recurrences not in normal form then the system
  tries to rewrite them.
  If the rewriting's process has success then in all the successive
  computations the system works with the recurrence rewritten in normal form.

  \par Example 4
  The recurrence of the second order
  \f$ x(n) = 5 x(n+1) - 6 x(n+2) \f$
  will be rewritten in the normal form
  \f$ x(n) = 5/6 x(n-1) - 1/6 x(n-2) \f$.

  \par Example 5
  The weighted-average recurrence
  \f$ x(n) = 2 \sum_{k=3}^n x(k) - 2 n \f$
  will be transformed in the normal form and, moreover,
  the recurrence will be rewritten, if possible, so that the lower limit of
  the sum is \f$ 0 \f$ and the upper limit is \f$ n-1 \f$.
  In general,
  \f[
    x(n) = f(n) \sum_{k=n_0}^n x(k) + g(n),
  \f]
  is rewritable in
  \f[
    x(n) = \frac{f(n)}{1-f(n)} \sum_{k=n_0}^{n-1} x(k) + \frac{g(n)}{1-f(n)}
  \f]
  when \f$ f(n) \neq 1 \f$ or in
  \f[
    x(n) = - \sum_{k=n_0+1}^{n-1} x(k) - g(n+1)
  \f]
  when \f$ f(n) = 1 \f$.
  Hence, in this example, the recurrence becomes
  \f$ x(n) = -2 \sum_{k=3}^{n-1} x(k) + 2 n \f$.
  Now we apply the second transformation: from
  \f[
    x(n) = f(n) \sum_{k=n_0}^{n-1} x(k) + g(n)
  \f]
  in
  \f[
    x(n) = f(n+n_0) \sum_{k=0}^{n-1} x(k) + g(n+n_0).
  \f]
  In conclusion, the system will work with the weighted-average recurrence
  \f$ x(n) = -2 \sum_{k=0}^{n-1} x(k) + 2 n + 6 \f$.
*/
class Recurrence {
public:
  //! Default constructor: builds the singleton satisfying \f$ x(n) = 0 \f$.
  Recurrence();

  //! \brief
  //! Builds the set of sequences satisfying \f$ x(n) = e \f$,
  //! with \f$ e \f$ the right-hand side of the recurrence.
  explicit Recurrence(const Expr& e);

  //! Copy-constructor.
  Recurrence(const Recurrence& y);

  //! Destructor.
  ~Recurrence();

  //! Assignment operator.
  Recurrence& operator=(const Recurrence& y);

  //! \brief
  //! Assigns to \p *this an approximation
  //! of the union of \p *this and \p y.
  void approximate_union_assign(const Recurrence& y);

  //! \brief
  //! Assigns to \p *this an approximation
  //! of the intersection of \p *this and \p y.
  void approximate_intersection_assign(const Recurrence& y);

  //! Sets to \f$ e \f$ the right-hand side of \p *this.
  void replace_recurrence(const Expr& e);

  //! \brief
  //! To use in the case of a system of recurrences: sets to \f$ e \f$
  //! the right-hand side of the recurrence of index \f$ k \f$.
  //! The system of recurrences will then include \f$ x_k(n) = e \f$.
  void replace_recurrence(unsigned int k, const Expr& e);

  //! Returns a new symbol \f$ z \f$ and records the equation \f$ z = e \f$.
  Symbol insert_auxiliary_definition(const Expr& e) const;

  //! \brief
  //! Returns an expression obtained from \p e by substituting the symbols
  //! that are also on the blackboard with their definition.
  Expr substitute_auxiliary_definitions(const Expr& e) const;

   //! \brief
  //! Returns the smallest number of initial conditions so that
  //! \p *this is uniquely identify.
  /*!
    \exception std::logic_error      thrown if \p *this is not classified yet
                                     and the classification's process
				     called by this method fails.
  */
  unsigned int number_initial_conditions() const;

  //! \brief
  //! Returns the smallest index of initial conditions acceptable
  //! for \p *this.
  /*!
    \exception std::logic_error      thrown if \p *this is not classified yet
                                     and the classification's process
				     called by this method fails.
  */
  // FIXME: It now just returns `first_valid_index'.
  index_type first_valid_initial_condition() const;

  //! \brief
  //! Returns
  //! \ref first_valid_index_for_solution "first valid index for solution".
  /*!
    \exception std::logic_error      thrown if this method is called
                                     but no exact solution, no lower bound
				     and no upper bound were computed.
  */
  index_type first_valid_index_for_solution() const;

  //! \brief
  //! In the case of functional equations, if the divisor \f$ b \f$ of
  //! the homogeneous term \f$ a x(n/b) \f$ in \p *this is different
  //! from \f$ 2 \f$, returns the definition of the function \f$ Sc() \f$:
  //! \f[
  //!   Sc(n, b) = \lfloor
  //!                \frac{n}{b^{\lfloor \frac{\log n}{\log b} \rfloor}}
  //!              \rfloor.
  //! \f]
  //! Returns the empty string if \f$ b = 2 \f$.
  /*!
    \exception std::logic_error      thrown if this method is called
                                     but no lower bound and no upper bound
				     were computed.
  */
  std::string definition_Sc() const;

  //! \brief
  //! Defines initial conditions of the recurrence that have been stored
  //! in the map \p initial_conditions.
  /*!
    The map \p initial_conditions represents a set of initial
    conditions \f$ x(i) = e \f$, corresponding to a pair \f$ (i, e) \f$
    in the map.
    This set of initial conditions must be sufficient to uniquely
    determine, together with the recurrence and given an assignment
    of values to the parameters, a sequence of complex numbers.
    In this context, the notion of <EM>sufficient</EM> depends
    on the type of the recurrence.
    In particular, called \f$ j \f$ the
    \ref first_valid_index "first valid index":
    - let \p *this be a
      \ref linear_finite_order "linear finite order recurrence"
      or a
      \ref non_linear "non linear finite order recurrence" of order \f$ k \f$.
      Then the initial conditions in the map must include the \f$ k \f$ pairs
      \f[
           (j-k+p, e_1), (j-k+p+1, e_2), (j-k+p+2, e_3), \ldots, (j+p-1, e_k)
      \f]
      for \f$ p \geq 1 \f$.  If \p initial_conditions contains such a sequence
      for more than one value of \f$ p \f$ the biggest of such values is used.
      Conditions in  \p initial_conditions for indices less than \f$ j-k+p \f$,
      if any, are interpreted as special cases.

    - Let \p *this be a \ref weighted_average "weighted-average recurrence".
      The only initial condition of the form \f$ (j, e) \f$ is sufficient
      to uniquely determine the sequence of complex numbers.
      If \p initial_conditions contains more initial conditions of the form
      \f$ (j+p, e_p) \f$, for \f$ p \geq 0 \f$, the pair with the biggest
      left value is used, while all the other conditions are interpreted as
      special cases.

    - a \ref generalized_recurrence "functional equation":         
                                     to be written.
 
    Calling this method overrides any previously given set of initial
    conditions.

    Initial conditions can be forgotten by means of the method
    <CODE>reset_initial_conditions()</CODE>.

    \par Example 1
    We consider the second order recurrence relation with constant
    coefficients \f$ x_n = 2*x_{n-1} - x_{n-2} \f$.
    In this example we first compute the symbolic exact solution
    \f$ x(n) = x_1*n + x_0 * (1-n) \f$; successively we evaluate it
    basing us on two different set of initial conditions
    (notice that are necessary at least \f$ 2 \f$ initial conditions
    with consecutive indices); at the end we forget the initial conditions
    given thanks to the method <CODE>reset_initial_conditions()</CODE>.
    \code
      Recurrence rec(2*x(n-1)-x(n-2));
      Recurrence::Solver_Status status = rec.compute_exact_solution();
      if (status == Recurrence::SUCCESS) {
        Expr exact_solution;
	rec.exact_solution(exact_solution);
	cout << "x(n) = " << exact_solution << endl;

	std::map<index_type, Expr> initial_conditions;
	initial_conditions[3] = 5;
	initial_conditions[4] = 2;
	initial_conditions[5] = 0;
	initial_conditions[6] = 1;
	rec.set_initial_conditions(initial_conditions);
	rec.exact_solution(exact_solution);
	cout << "x(n) = " << exact_solution << endl;
	
	initial_conditions.clear();
	Symbol a("a");
	Symbol b("b");
	initial_conditions[0] = a;
	initial_conditions[1] = b;
	rec.set_initial_conditions(initial_conditions);
	rec.exact_solution(exact_solution);
	cout << "x(n) = " << exact_solution << endl;
	
	rec.reset_initial_conditions();
	rec.exact_solution(exact_solution);
	cout << "x(n) = " << exact_solution << endl;
      }
    \endcode
    These instructions printed the following rows:
    \code
      x(n) = x(1)*n+x(0)*(1-n)
      x(n) = -5+n
      x(n) = -n*a+a+n*b
      x(n) = x(1)*n+x(0)*(1-n)
    \endcode
    
    \exception std::logic_error      thrown if \p *this is not classified yet
                                     and the classification's process
				     called by this method fails.
    
    \exception std::invalid_argument thrown if the initial conditions
                                     in the map \p initial_conditions
				     are not a sufficient set to uniquely
				     determine a sequence of complex
				     numbers or if these initial
				     conditions are not in agreement with
				     the type of the recurrence.
  */
  void
  set_initial_conditions(const std::map<index_type, Expr>& initial_conditions);

  //! Resets initial conditions, if any.
  void reset_initial_conditions();

  //! Returns the map containing the set of initial conditions.
  const std::map<index_type, Expr>& get_initial_conditions() const;

  //! Returns the exact solution of \p *this evaluated for \f$ n = x \f$.
  /*!
    This method replaces <CODE>Recurrence::n</CODE> with \p x in
    the expression containing the exact solution.
    In general must be \f$ x \geq 0 \f$ and, moreover, \p x must respect
    some rules depends on the type of the recurrence.
    In particular, called \f$ j \f$ the
    \ref first_valid_index "first valid index":
    - let \p *this be a
      \ref linear_finite_order "linear finite order recurrence" of order
      \f$ k \f$.
      If \f$ h = j - k + 1 > 0 \f$ must be \f$ x \geq h \f$ or,
      if \f$ x \in [0, h) \f$, must exist an initial condition with the
      index equal to \p x.

    - Let \p *this be a \ref weighted_average "weighted-average recurrence".
      Then it must be \f$ x > j \f$  or must exist an initial condition
      with the index equal to \p x.

    - a \ref generalized_recurrence "functional equation":         
                                     to be written.

    - a \ref non_linear "non linear finite order recurrence"
                                     to be written.

    Notice that if the initial conditions were specified this method
    returns a numeric expression; otherwise it returns a symbolic expression.

    \exception std::logic_error       thrown if this method is called
                                      but no exact solution was computed.

    \exception std::invalid_argument  thrown if \p x is not in agreement
                                      with the solution of the recurrence.
  */
  Expr evaluate_exact_solution(const Number& x) const;

  //! \brief
  //! Evaluates the exact solution for \f$ n = i \f$,
  //! where \f$ i \f$ assumes all the values of the interval
  //! \f$ [begin, end) \f$. Puts the results in a container marked by \p oi.
  /*!
    This method replaces <CODE>Recurrence::n</CODE> with all numbers
    in the interval \f$ [begin, end) \f$ in the expression containing the
    exact solution.
    \f$ i \in [begin, end) \f$ must be a non-negative integer in agreement with
    the least non-negative integer \f$ j \f$ such that the 
    recurrence is well-defined for \f$ n \geq j \f$ (for more details
    see the method <CODE>evaluate_exact_solution()</CODE> with only
    one argument).

    Notice that if the initial conditions were specified the expressions
    returned by this method can be numerics; otherwise they will be symbolics.

    \exception std::logic_error       thrown if this method is called
                                      but no exact solution was computed.

    \exception std::invalid_argument  thrown if \f$ begin > end \f$ or
                                      \p begin is not in agreement
                                      with the solution of the recurrence.
  */
  template <class OutputIterator>
  void
  evaluate_exact_solution(const Number& begin, const Number& end,
			  OutputIterator oi) const;
  
  //! \brief
  //! Returns the lower bound of \p *this evaluated for \f$ n = x \f$.
  /*!
    This method replaces <CODE>Recurrence::n</CODE> with \p x in
    the expression containing the lower bound.
    In general must be \f$ x \geq 0 \f$ and, moreover, \p x must respect
    some rules depends on the type of the recurrence (for more details
    see the method <CODE>evaluate_exact_solution()</CODE> with only
    one argument).

    Notice that if the initial conditions were specified this method
    returns a numeric expression; otherwise it returns a symbolic expression.

    \exception std::logic_error       thrown if this method is called
                                      but no lower bound was computed.
  
    \exception std::invalid_argument  thrown if \p x is not in agreement
                                      with the lower bound.
  */
  Expr evaluate_lower_bound(const Number& x) const;

  //! \brief
  //! Evaluates the lower bound for \f$ n = i \f$,
  //! where \f$ i \f$ assumes all the values of the interval
  //! \f$ [begin, end) \f$. Puts the results in a container marked by \p oi.
  /*!
    This method replaces <CODE>Recurrence::n</CODE> with \p x in
    the expression containing the upper bound.
    In general must be \f$ x \geq 0 \f$ and, moreover, \p x must respect
    some rules depends on the type of the recurrence (for more details
    see the method <CODE>evaluate_exact_solution()</CODE> with only
    one argument).

    Notice that if the initial conditions were specified the expressions
    returned by this method can be numerics; otherwise they will be symbolics.

    \exception std::logic_error       thrown if this method is called
                                      but no lower bound was computed.

    \exception std::invalid_argument  thrown if \f$ begin > end \f$ or
                                      \p begin is not in agreement
                                      with the lower bound.
  */
  template <class OutputIterator>
  void
  evaluate_lower_bound(const Number& begin, const Number& end,
		       OutputIterator oi) const;

  //! Returns the upper bound of \p *this evaluated for \f$ n = x \f$.
  /*!
    This method replaces <CODE>Recurrence::n</CODE> with \p x in
    the expression containing the upper bound.
    \p x must be a non-negative integer in agreement with
    the least non-negative integer \f$ j \f$ such that the 
    recurrence is well-defined for \f$ n \geq j \f$ (for more details
    see the method <CODE>evaluate_exact_solution()</CODE> with only
    one argument).

    Notice that if the initial conditions were specified this method
    returns a numeric expression; otherwise it returns a symbolic expression.

    \exception std::logic_error       thrown if this method is called
                                      but no upper bound was computed.
  
    \exception std::invalid_argument  thrown if \p x is not in agreement
                                      with the upper bound.
  */
  Expr evaluate_upper_bound(const Number& x) const;

  //! \brief
  //! Evaluates the upper bound for \f$ n = i \f$,
  //! where \f$ i \f$ assumes all the values of the interval
  //! \f$ [begin, end) \f$. Puts the results in a container marked by \p oi.
  /*!
    This method replaces <CODE>Recurrence::n</CODE> with all numbers
    in the interval \f$ [begin, end) \f$ in the expression containing the
    upper bound.
    \f$ i \in [begin, end) \f$ must be a non-negative integer in agreement with
    the least non-negative integer \f$ j \f$ such that the 
    recurrence is well-defined for \f$ n \geq j \f$ (for more details
    see the method <CODE>evaluate_exact_solution()</CODE> with only
    one argument).

    Notice that if the initial conditions were specified the expressions
    returned by this method can be numerics; otherwise they will be symbolics.

    \exception std::logic_error       thrown if this method is called
                                      but no upper bound was computed.

    \exception std::invalid_argument  thrown if \f$ begin > end \f$ or
                                      \p begin is not in agreement
                                      with the upper bound.
  */
  template <class OutputIterator>
  void
  evaluate_upper_bound(const Number& begin, const Number& end,
		       OutputIterator oi) const;

  //! Returns the right-hand side of \p *this evaluated for \f$ n = x \f$.
  /*!
    

    \exception std::logic_error       thrown if \p *this is not classified yet
                                      and the classification's process
				      called by this method fails.

    \exception std::invalid_argument  thrown if \p x is not in agreement
                                      with the recurrence.
  */
  Expr evaluate_rhs(const Number& x) const;

  //! \brief
  //! Evaluates the right-hand side of \p *this for \f$ n = i \f$,
  //! where \f$ i \f$ assumes all the values of the interval
  //! \f$ [begin, end) \f$. Puts the results in a container marked by \p oi.
  /*!
    \exception std::logic_error       thrown if \p *this is not classified yet
                                      and the classification's process
				      called by this method fails.

    \exception std::invalid_argument  thrown if \f$ begin > end \f$  or
                                      \p begin is not in agreement
                                      with the recurrence.

  */
  template <class OutputIterator>
  void
  evaluate_rhs(const Number& begin, const Number& end,
	       OutputIterator oi) const;

#ifdef PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS
  //! Checks if all the invariants are satisfied.
  /*!
    The check is performed so as to intrude as little as possible.
    In case invariants are violated error messages are written on
    <CODE>std::cerr</CODE>.
  */
#endif // PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS
  bool OK() const; 

private:
  //! Kinds of bounds for the solution that is possible to compute.
  enum Bound {
    /*! An upper bound. */
    UPPER,
    
    /*! A lower bound. */
    LOWER
  };

  //! \brief
  //! Returns the right-hand side of the auxiliary equation \f$ z = e \f$,
  //! if such an auxiliary equation exists;
  //! returns the expression \f$ z \f$ otherwise.
  Expr get_auxiliary_definition(const Symbol& z) const;

  //! Sets to \p e the <CODE>inhomogeneous_term</CODE>.
  void set_inhomogeneous_term(const Expr& e) const;

  //! \brief
  //! If in the map <CODE>initial_conditions</CODE> there is the
  //! expression \f$ e \f$ correspondent to \p k then returns \f$ e \f$;
  //! returns \f$ x(k) \f$ otherwise.
  Expr get_initial_condition(unsigned int k) const;

  //! \brief
  //! If the map <CODE>initial_conditions</CODE> is not empty returns
  //! the maximum integer among the first elements of it;
  //! returns \f$ 0 \f$ if the map is empty.
  unsigned int get_max_index_initial_condition() const;

  Expr evaluate(unsigned int kind, const Number& num) const;

  void check_number_for_evaluation(const char* method, const char* name,
				   const Number& num) const;

  //! \brief
  //! Throw a <CODE>std::invalid_argument</CODE> exception
  //! containing the appropriate error message.
  void
  throw_invalid_argument(const char* method, const char* reason) const;

public:
  //! The possible states of the recurrence.
  enum Solver_Status {
    /*!
      Solution, or approximation, was successful.
    */
    SUCCESS,

    /*!
      The recurrence is indeterminate, hence it has infinitely
      many solutions.
    */
    INDETERMINATE_RECURRENCE,

    /*!
      The recurrence is not solvable.
    */
    UNSOLVABLE_RECURRENCE,

    /*!
      The recurrence is not
      \ref syntactically_correct "syntactically correct".
    */
    MALFORMED_RECURRENCE,

    /*!
      The recurrence is not well-defined
      (ex. \f$ x_n = (n-2) \sum_{k=0}^{n-1} x(k) \f$).
    */
    DOMAIN_ERROR,

    /*!
      Catchall: the recurrence is generically too complex for the solver.
    */
    TOO_COMPLEX
  };

  //! \brief
  //! Kinds of responses of validation's process of the solution, or
  //! of the approximation, of the recurrence.
  enum Verify_Status {
    /*!
      The system can prove that the recurrence \p *this has been
      successfully solved or approximated.
    */
    PROVABLY_CORRECT,

    /*!
      The system can not prove whether the solution, or the approximation,
      of \p *this, is correct or not, but it can prove that the solution,
      or the approximation, of another recurrence used to compute the
      solution, or the approximation, of \p *this, is correct.
      This is the case of
      -  linear recurrences of finite order solved with the order reduction
         method. The system tries to verify the solution (expanded,
	 e.g, without using the function \f$ mod() \f$) of the original
	 recurrence and if it fails tries to verify the solution of the
	 reduced recurrence.
      -  non-linear recurrences of finite order. If the verification
         of the solution of the non-linear recurrence fails then it
	 tries to verify the solution of the linear recurrence associated
	 to \p *this.
      -  weighted-average recurrences. If the verification
         of the solution of the weighted-average recurrence fails then it
	 tries to verify the solution of the first order recurrence associated
	 to \p *this.
    */
    PARTIAL_PROVABLY_CORRECT,

    /*!
      The system can prove that the solution, or the approximation, of
      \p *this, is wrong.
    */
    PROVABLY_INCORRECT,

    /*!
      The system can not prove whether the solution, or the approximation,
      of \p *this, is correct or not. 
    */
    INCONCLUSIVE_VERIFICATION
  };

  //! \brief
  //! Tries to solve \p *this exactly and returns <CODE>SUCCESS</CODE>
  //! if the system finds the exact solution.
  /*!
    The system tries to compute the exact solution of every type
    of recurrence, linear and non-linear of finite order, weighted-average
    and also with functional equations, possibly calling the methods
    that compute lower bound and upper bound and verifying if they are
    the same.

    If the method returns <CODE>SUCCESS</CODE>, you can get the
    exact solution calling the method <CODE>exact_solution()</CODE>;
    otherwise you are in one of the cases explained in the enumerate
    type <CODE>Solver_Status</CODE>.
  */
  Solver_Status compute_exact_solution() const;

  //! Gets the exact solution and puts it in \p e.
  /*!
    We assume that the exact solution has already been computed:
    in this case the method <CODE>compute_exact_solution()</CODE>
    returns <CODE>SUCCESS</CODE>.

    \par Example 1
    A correct use of this method can be the following:
    \code
        Recurrence rec(2*x(n-1)+1);
	if (rec.compute_exact_solution() == Recurrence::SUCCESS) {
	  Expr exact_solution;
	  rec.exact_solution(exact_solution);
	}
    \endcode
    The solution of the recurrence \f$ x_n = 2 x_{n-1} + 1 \f$ is
    \f[
      x(n) = x_0 2^n + 2^n - 1
    \f]
    and its right-hand side is contained in the variable
    <CODE>exact_solution</CODE>.

    We remark that in the exact solution will appear symbolically
    the initial conditions until the user will not define it: this
    can be performed with the method
    <CODE>set_initial_conditions()</CODE>.

    \exception std::logic_error thrown if this method is called
                                but no exact solution was computed.
  */
  void exact_solution(Expr& e) const;

  //! \brief
  //! Tries to get lower bound for \p *this and returns <CODE>SUCCESS</CODE>
  //! if the system finds the lower bound.
  /*!
    The system tries to compute the lower bound for the solution of
    every type of recurrence, not only with functional equations.
    In the case of linear and non-linear recurrences of finite order
    and weighted-average recurrences, it tries to compute the
    exact solution and, if it was successful, the solution found obviously
    represent also a lower bound for the solution.

    If the method returns <CODE>SUCCESS</CODE>, you can get the
    lower bound calling the method <CODE>lower_bound()</CODE>;
    otherwise you are in one of the cases explained in the enumerate
    type <CODE>Solver_Status</CODE>.    
  */
  Solver_Status compute_lower_bound() const;

  //! Gets the lower bound for the solution and puts it in \p e.
  /*!
    We assume that the lower bound has already been computed:
    in this case the method <CODE>compute_lower_bound()</CODE>
    returns <CODE>SUCCESS</CODE>.

    \par Example 1
    A correct use of this method can be the following:
    \code
        Recurrence rec(7*x(n/5)+10);
	if (rec.compute_lower_bound() == Recurrence::SUCCESS) {
	  Expr lower_bound;
	  rec.lower_bound(lower_bound);
	}
    \endcode
    The lower bound of the solution of the functional equation
    \f$ x(n) = 7 x(n/5)+10 \f$ is
    \f[
      x(n) \geq \frac{1}{7} x(1) n^{\log 7 / \log 5} - \frac{5}{3}
                + \frac{5}{21} n^{\log 7 / \log 5}
    \f]
    and its right-hand side is contained in the variable
    <CODE>lower_bound</CODE>.

    We remark that in the lower bound will appear symbolically
    the initial conditions until the user will not define it: this
    can be performed with the method
    <CODE>set_initial_conditions()</CODE>.

    \exception std::logic_error thrown if this method is called
                                but no lower bounds was computed.
  */
  void lower_bound(Expr& e) const;

  //! \brief
  //! Tries to get upper bound for \p *this and returns <CODE>SUCCESS</CODE>
  //! if the system finds the upper bound.
  /*!
    The system tries to compute the upper bound for the solution of
    every type of recurrence, not only with functional equations.
    In the case of linear and non-linear recurrences of finite order
    and weighted-average recurrences, it tries to compute the
    exact solution and, if it was successful, the solution found obviously
    represent also an upper bound for the solution.

    If the method returns <CODE>SUCCESS</CODE>, you can get the
    upper bound calling the method <CODE>upper_bound()</CODE>;
    otherwise you are in one of the cases explained in the enumerate
    type <CODE>Solver_Status</CODE>.    
  */
  Solver_Status compute_upper_bound() const;

  //! Gets the upper bound for the solution and puts it in \p e.
  /*!
    We assume that the upper bound has already been computed:
    in this case the method <CODE>compute_upper_bound()</CODE>
    returns <CODE>SUCCESS</CODE>.

    \par Example 1
    A correct use of this method can be the following:
    \code
        Recurrence rec(7*x(n/5)+10);
	if (rec.compute_upper_bound() == Recurrence::SUCCESS) {
	  Expr upper_bound;
	  rec.upper_bound(upper_bound);
	}
    \endcode
    The upper bound of the solution of the functional equation
    \f$ x(n) = 7 x(n/5)+10 \f$ is
    \f[
      x(n) \leq x(1) n^{\log 7 / \log 5} - \frac{5}{3}
                + \frac{5}{3} n^{\log 7 / \log 5}
    \f]
    and its right-hand side is contained in the variable
    <CODE>upper_bound</CODE>.

    We remark that in the upper bound will appear symbolically
    the initial conditions until the user will not define it: this
    can be performed with the method
    <CODE>set_initial_conditions()</CODE>.

    \exception std::logic_error thrown if this method is called
                                but no upper bounds was computed.
  */
  void upper_bound(Expr& e) const;

  //! \brief
  //! Tries to verify the exact solution of \p *this and returns
  //! <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  /*!
    We assume that the exact solution has already been computed
    (for a description about the computation of the exact solution
    see the documentation of the method <CODE>compute_exact_solution()</CODE>).

    There are three different answers when the system does not succeed
    to verify the solution:
    -  returns <CODE>PROVABLY_INCORRECT</CODE> if the system proved
       that the solution of the recurrence \p *this is wrong;
    -  returns <CODE>PARTIAL_PROVABLY_CORRECT</CODE> in the following
       cases:
       -  when the system fails in the verification of the solution (expanded,
          e.g, without using the function \f$ mod() \f$) of linear finite
	  order recurrences solved with the order reduction method but
	  it has success in the validation of the solution of the associated
	  reduced recurrences.
       -  when the system fails in the verification of the solution of the
          non-linear recurrence but it has success in the validation
	  of the solution of the associated linear recurrence.
       -  when the system fails in the verification of the solution of the
          weighted-average recurrence but it has success in the validation
	  of the solution of the associated first order recurrence.
    -  returns <CODE>INCONCLUSIVE_VERIFICATION</CODE> when the system is
       not able to prove if the solution of the recurrence \p *this
       is correct or not.

    \exception std::logic_error thrown if this method is called
                                when no exact solution was computed.
  */
  Verify_Status verify_exact_solution(bool partial_verification = false) const;
 
  //! \brief
  //! Tries to verify the lower bound of the solution of \p *this and
  //! returns <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  /*!
    We assume that the lower bound of the solution has already been computed
    (for a description about the computation of the approximation
    from below see the documentation of the method
    <CODE>compute_lower_bound()</CODE>).

    There are two different answers when the system does not succeed
    to verify the solution:
    -  returns <CODE>PROVABLY_INCORRECT</CODE> if the system proved
       that the lower bound of the solution of \p *this is wrong;
    -  returns <CODE>INCONCLUSIVE_VERIFICATION</CODE> when the system is
       not able to prove if the lower bound of the solution of \p *this
       is correct or not.

    \exception std::logic_error thrown if this method is called
                                when no lower bound was computed.
  */
  Verify_Status verify_lower_bound() const;

  //! \brief
  //! Tries to verify the upper bound of the solution of \p *this and
  //! returns <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  /*!
    We assume that the lower bound of the solution has already been computed
    (for a description about the computation of the approximation
    from over see the documentation of the method
    <CODE>compute_upper_bound()</CODE>).

    There are two different answers when the system does not succeed
    to verify the solution:
    -  returns <CODE>PROVABLY_INCORRECT</CODE> if the system proved
       that the upper bound of the solution of \p *this is wrong;
    -  returns <CODE>INCONCLUSIVE_VERIFICATION</CODE> when the system is
       not able to prove if the upper bound of the solution of \p *this
       is correct or not.

    \exception std::logic_error thrown if this method is called
                                when no upper bound was computed.
  */
  Verify_Status verify_upper_bound() const;

  //! The index of the recurrence.
  /*!
    The symbol \p n is special in the context of recurrence relations:
    in order to avoid ambiguities is better do not use this symbol as
    any other parameter.
  */
  static const Symbol& n;

  //! Dumps all the data members of \p *this onto \p s.
  void dump(std::ostream& s) const;
  
  // FIXME: this public section is temporary.
public:
  //private:
  //! The possible classification of the recurrence.
  enum Classifier_Status {
    /*!
      The classification's process was successful: the type of the
      recurrence is known and the solver is able to work with it
    */
    CL_SUCCESS,

    /*!
      The recurrence is indeterminate, hence it has infinitely many
      solutions.
    */
    CL_INDETERMINATE_RECURRENCE,

    /*!
      The recurrence is not solvable.
    */
    CL_UNSOLVABLE_RECURRENCE,

    /*!
      The recurrence is not
      \ref syntactically_correct "syntactically correct".
    */
    CL_MALFORMED_RECURRENCE,

    /*!
      The recurrence is not well-defined
      (ex. \f$ x_n = (n-2) \sum_{k=0}^{n-1} x(k) \f$).
    */
    CL_DOMAIN_ERROR,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is not an integer.
    */
    CL_HAS_NON_INTEGER_DECREMENT,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is too big to be handled
      by the standard solution techniques.
    */
    CL_HAS_HUGE_DECREMENT,

    /*!
      Catchall: the recurrence is generically too complex for the solver.
    */
    CL_TOO_COMPLEX,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is a negative integer.
    */
    CL_HAS_NEGATIVE_DECREMENT,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n)</CODE>.
    */
    CL_HAS_NULL_DECREMENT,

    /*!
      The recurrence is not yet classified.
    */
    NOT_CLASSIFIED_YET
  };

  mutable Classifier_Status classifier_status_;

private:
  //! \brief
  //! Returns the status of <CODE>Solver_Status</CODE> associated
  //! to \p classifier_status.
  static Solver_Status map_status(Classifier_Status classifier_status);

  // FIXME: this public section is temporary.
public:
  //! \brief
  //! Classifies the recurrence \p *this calling the method
  //! <CODE>classify()</CODE>. 
  Classifier_Status classify_and_catch_special_cases() const;

private:
  //! \brief
  //! Analyzes the \f$ i \f$-th addend of the right hand side \p rhs
  //! of the recurrence \p *this.
  Classifier_Status classification_summand(const Expr& r, Expr& rhs,
					   Expr& e, index_type& order,
					   std::vector<Expr>& coefficients,
					   int& gcd_among_decrements,
					   int num_term,
					   std::map<Number, Expr>&
					   homogeneous_terms) const;

  //! \brief
  //! Computes the exact solution of \p *this, where \p *this is a
  //! linear recurrence of finite order.
  Solver_Status compute_exact_solution_finite_order() const;

  //! \brief
  //! Computes the exact solution of \p *this, where \p *this is a
  //! functional equation.
  Solver_Status compute_exact_solution_functional_equation() const;

  //! \brief
  //! Computes the exact solution of \p *this, where \p *this is a
  //! non linear recurrence of finite order.
  Solver_Status compute_exact_solution_non_linear() const;
  
  //! \brief
  //! Computes the exact solution of \p *this, where \p *this is a
  //! weighted-average recurrences.
  Solver_Status compute_exact_solution_weighted_average() const;

  //! \brief
  //! Computes the exact solution of \p *this, where \p *this is a
  //! recurrencs defined as the maximum of two quantities.
  Solver_Status compute_exact_solution_max() const;

  //! \brief
  //! Computes the requested bound (lower or upper, specified by
  //! \p kind_of_bound) of the functional equation \p *this.
  Solver_Status
  compute_bound_functional_equation(Bound kind_of_bound) const;

  //! \brief
  //! Computes the requested bound (lower or upper, specified by
  //! \p kind_of_bound) of the non-linear recurrence \p *this.
  Solver_Status
  compute_bound_non_linear(Bound kind_of_bound) const;

  //! \brief
  //! Computes the solution of \p *this applying the order reduction method:
  //! builds a new object <CODE>Recurrence</CODE> with the reduced recurrence;
  //! solves the reduced recurrence; from the solution of the reduced
  //! recurrence finds the solution of the original recurrence \p *this.
  Solver_Status apply_order_reduction() const;

  //! \brief
  //! Computes the solution of \p *this: transforms \p *this in a linear
  //! recurrence and stores it in a new object <CODE>Recurrence</CODE>;
  //! solves the linear recurrence; from the solution of the linear
  //! recurrence finds the solution of the original recurrence \p *this.
  Solver_Status compute_non_linear_recurrence(Expr& solution_or_bound,
					      unsigned int type) const;

  //! \brief
  //! Solves the weighted-average recurrence \p *this transforming it
  //! in a linear recurrence of finite order stored in a new object
  //! <CODE>Recurrence</CODE>.
  Solver_Status compute_weighted_average_recurrence(Expr& solution) const;

  Expr
  compute_solution_linear_finite_order_on_i_c(const Expr& solution) const;

  Expr
  compute_solution_non_linear_finite_order_on_i_c(const Expr& solution) const;

  Expr
  compute_bound_functional_equation_on_i_c(const Expr& bound) const;

  Expr
  compute_solution_weighted_average_on_i_c(const Expr& solution) const;

  //! \brief
  //! Let \p solution_or_bound be the expression that represent the
  //! solution or the bound computed for the recurrence \p *this.
  //! This function substitutes possible initial conditions specified
  //! by the user shifting the solution or the bound if necessary.
  Expr compute_solution_or_bound_on_i_c(const Expr& solution_or_bound) const;

  //! Classifies the recurrence \p *this sliding it recursively.
  Classifier_Status classify() const;

  //! \brief
  //! Solves the linear recurrence of finite order.
  //! The recurrence can have constant or variable coefficients.
  Solver_Status solve_linear_finite_order() const;

  //! Solves the linear recurrence of first order with constant coefficient.
  Expr solve_constant_coeff_order_1(const std::vector<Polynomial_Root>&
				    roots) const;

  //! Solves the linear recurrence of first order with variable coefficient. 
  Solver_Status
  solve_variable_coeff_order_1(const std::vector<Expr>& coefficients,
			       Expr& solution) const;

  //! Solves the linear recurrence of second order with constant coefficients. 
  Expr solve_constant_coeff_order_2(Expr& g_n, bool all_distinct,
				    const std::vector<Number>& coefficients,
				    const std::vector<Polynomial_Root>&
				    roots) const;

  //! \brief
  //! Solves the linear recurrence of finite order \f$ k \f$, \f$ k > 2 \f$,
  //! with constant coefficients.
  Expr solve_constant_coeff_order_k(Expr& g_n, bool all_distinct,
				    const std::vector<Number>& coefficients,
				    const std::vector<Polynomial_Root>&
				    roots) const;

  //! \brief
  //! Performs the second step of the validation's process of exact solutions
  //! of the recurrence relation of finite order \p *this verifying
  //! \ref initial_conditions "symbolic initial conditions".
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  Verify_Status
  validate_initial_conditions(bool symbolic_solution, index_type order,
			      const std::vector<Expr>& coefficients_i_c,
			      const Expr& summands_with_i_c,
			      const Expr& summands_without_i_c) const;

  Verify_Status
  traditional_step_3(index_type order_rec,
		     const Expr& summands_with_i_c) const;

  //! \brief
  //! Apply the method explained in the paper
  //! <EM>Checking and Confining the Solutions of Recurrence Realtions</EM>
  //! in order to validate solution of linear finite order with constant
  // coefficients.
  bool verify_new_method_const_coeff(index_type order_rec, const Expr& e,
				     const std::vector<Polynomial_Root> roots,
				     bool inhomogeneous_part) const;

  //! \brief
  //! Apply the method explained in the paper
  //! <EM>Checking and Confining the Solutions of Recurrence Realtions</EM>
  //! in order to validate solution of linear finite order with variable
  // coefficients.
  bool verify_new_method_var_coeff(index_type order_rec,
				   const Expr& summands_without_i_c) const;

  //! \brief
  //! Verifies the exact solution of the finite order recurrence \p *this.
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  Verify_Status verify_finite_order(bool partial_verification) const;

  //! \brief
  //! Verifies the exact solution of the weighted-average recurrence \p *this.
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  Verify_Status verify_weighted_average() const;

  //! \brief
  //! Performs the first step of the validation's process of the lower or
  //! the upper bound of the functional equation \p *this verifying
  //! \ref initial_conditions "symbolic initial conditions".
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  Verify_Status validate_initial_conditions(Bound kind_of_bound,
					    const Expr& bound,
					    unsigned int index) const;

  //! \brief
  //! Verifies the lower or the upper bound of the functional equation
  //! \p *this.
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system is successful in
  //! the verification.
  Verify_Status verify_bound(Bound kind_of_bound) const;

  //! \brief
  //! If \p kind_of_bound is <CODE>LOWER</CODE> computes the lower bound
  //! for the functional equation \p *this; otherwise, if \p kind_of_bound
  //! is <CODE>UPPER</CODE> computes the upper bound.
  Solver_Status approximate_functional_equation(Bound kind_of_bound) const;

  //! \brief
  //! Holds the right-hand side of the global recurrence to be solved.
  //! This may have been set directly by the constructor or it may be the
  //! result of transforming a system into a single recurrence.
  //! The global recurrence is thus of the form
  //! <CODE>x(n) = recurrence_rhs</CODE>.
  Expr recurrence_rhs;

  //! \brief
  //! It is <CODE>true</CODE> if the recurrence has been rewritten, i. e.,
  //! in the case there are null or negative decrements or if has been applied
  //! the order reduction or the recurrence is non-linear (and then, in
  //! order to solve it we have rewritten the recurrence in linear);
  //! it is <CODE>false</CODE> in all the other cases.
  bool recurrence_rewritten;

  //! \brief
  //! Stores the inhomogeneous part of \p *this, i. e., those terms
  //! that not involve \f$ x \f$ functions like \f$ x(n-k) \f$ or
  //! \f$ x(n/k) \f$, where \f$ k \f$ is a positive integer.
  mutable Expr inhomogeneous_term;

  //! \brief
  //! Holds the right-hand sides of a system of  recurrence equations.
  //! If <CODE>i == system_rhs.find(k)</CODE> then
  //! <CODE>x(k,n) = (*i).second()</CODE>
  //! is one of the equations of the system.
  std::map<unsigned int, Expr> system_rhs;

  // FIXME: this public section is temporary.
public:
  //! The recurrence type.
  enum Type {
    /*!
      Special recurrence of the form \f$ x(n) = rhs \f$, where \f$ rhs \f$
      contains only functions of  \f$ n \f$, parameters and
      \f$ x(k_1), \dots, x(k_m) \f$ where \f$ m >= 0 \f$ and
      \f$ k_1, \dots, k_m \f$ are non-negative integers.
    */
    ORDER_ZERO,

    /*!
      Linear recurrence of finite order with constant coefficient.
    */
    LINEAR_FINITE_ORDER_CONST_COEFF,
    
    /*!
      Linear recurrence of finite order with variable coefficient.
    */
    LINEAR_FINITE_ORDER_VAR_COEFF,
    
    /*!
      Non-linear recurrence of finite order.
    */
    NON_LINEAR_FINITE_ORDER,
    
    /*!
      Weighted-average recurrence. 
    */
    WEIGHTED_AVERAGE,
    
    /*!
      Special recurrence of the form \f$ x(n) = a x(n / b) + d n^e \f$ . 
    */
    FUNCTIONAL_EQUATION,

    /*!
      Special recurrence of the form 
      \f$ x(n) = max(f(x(0),..x(n)), g(x(0),..x(n))) \f$ . 
    */
    MAX_FUNCTION,

    /*!
      The recurrence is not yet classified, then the type of the
      recurrence is unknown.
    */
    UNKNOWN
  };

  mutable Type type_;

private:
  mutable Finite_Order_Info* finite_order_p;
  mutable Functional_Equation_Info* functional_eq_p;
  mutable Non_Linear_Info* non_linear_p;
  mutable Weighted_Average_Info* weighted_average_p;

  //! Returns <CODE>type_</CODE>.
  const Type& type() const;

  //! Returns <CODE>type_</CODE>.
  Type& type();

  //! Sets <CODE>type_</CODE> with \p t.
  void set_type(const Type& t) const;

  //! Sets <CODE>type_recurrence = ORDER_ZERO</CODE>.
  void set_order_zero() const;

  //! \brief
  //! Returns <CODE>true</CODE> if the recurrence is linear
  //! of finite order with constant coefficient;
  //! returns <CODE>false</CODE> otherwise.
  bool is_linear_finite_order_const_coeff() const;

  //! Sets <CODE>type_recurrence = LINEAR_FINITE_ORDER_CONST_COEFF</CODE>.
  void set_linear_finite_order_const_coeff() const;

  //! \brief
  //! Returns <CODE>true</CODE> if the recurrence is linear
  //! of finite order with variable coefficient;
  //! returns <CODE>false</CODE> otherwise.
  bool is_linear_finite_order_var_coeff() const;

  //! Sets <CODE>type_recurrence = LINEAR_FINITE_ORDER_VAR_COEFF</CODE>.
  void set_linear_finite_order_var_coeff() const;

  //! \brief
  //! Returns <CODE>true</CODE> if the recurrence is linear
  //! of finite order; returns <CODE>false</CODE> otherwise.
  bool is_linear_finite_order() const;

  //! \brief
  //! Returns <CODE>true</CODE> if the recurrence is non linear
  //! of finite order; returns <CODE>false</CODE> otherwise.
  bool is_non_linear_finite_order() const;

  //! Sets <CODE>type_recurrence = NON_LINEAR_FINITE_ORDER</CODE>.
  void set_non_linear_finite_order() const;

  //! \brief
  //! Returns <CODE>true</CODE> if is the case of functional equation;
  //! returns <CODE>false</CODE> otherwise.
  bool is_functional_equation() const;

  //! Sets <CODE>type_recurrence = FUNCTIONAL_EQUATION</CODE>.
  void set_functional_equation() const;

  //! \brief
  //! Returns <CODE>true</CODE> if \p *this is a weighted-average
  //! recurrence; returns <CODE>false</CODE> otherwise.
  bool is_weighted_average() const;

  //! Sets <CODE>type_recurrence = WEIGHTED_AVERAGE</CODE>.
  void set_weighted_average() const;


  // Methods to access to private data of `Finite_Order_Info'.

  //! Returns the order of the finite order recurrence.
  index_type order() const;

  //! Returns the coefficients of the linear finite order recurrence.
  const std::vector<Expr>& coefficients() const;

  //! Returns the coefficients of the linear finite order recurrence.
  std::vector<Expr>& coefficients();

  // FIXME: this public section is temporary.
public:
  //! \brief
  //! Returns the greatest common divisor among the decrements \f$ k \f$
  //! of the terms \f$ x(n-k) \f$ of a linear finite order recurrence.
  //! Returns \f$ 0 \f$ if the order of the recurrence is \f$ 0 \f$.
  unsigned int gcd_among_decrements() const;

private:
  //! \brief
  //! Returns the expression \f$ \prod_{i}^n a(k)\f$,
  //! where \f$ i \f$ is the \ref first_valid_index "first valid index"
  //! and \f$ a(n) \f$ is the coefficient of the recurrence of first order
  //! with variable coefficient.
  const Expr& product_factor() const;

  //! \brief
  //! Returns the expression \f$ \prod_{i}^n a(k)\f$,
  //! where \f$ i \f$ is the \ref first_valid_index "first valid index"
  //! and \f$ a(n) \f$ is the coefficient of the recurrence of first order
  //! with variable coefficient.
  Expr& product_factor();

  //! \brief
  //! Sets to \p x the expression \f$ \prod_{i}^n a(k)\f$,
  //! where \f$ i \f$ is the \ref first_valid_index "first valid index"
  //! and \f$ a(n) \f$ is the coefficient of the recurrence of first order
  //! with variable coefficient.
  void set_product_factor(const Expr& x) const;

  //! \brief
  //! Returns <CODE>true</CODE> if the order reduction method has been
  //! applied; returns <CODE>false</CODE> otherwise.
  bool applied_order_reduction() const;

  //! \brief
  //! Sets the flag that indicates if the order reduction method has been
  //! applied to <CODE>true</CODE>.
  void set_order_reduction() const;

  //! \brief
  //! Sets the flag that indicates if the order reduction method has been
  //! applied to <CODE>false</CODE>.
  void unset_order_reduction() const;
  

  // Methods to access to private data of `Functional_Equation_Info'.

  //! \brief
  //! Returns the positive integer starting from which the inhomogeneous term
  //! of a functional equation is a non negative, non decreasing function.
  index_type applicability_condition() const;

  //! \brief
  //! \p c is the positive integer starting from which the inhomogeneous term
  //! of a functional equation is a non negative, non decreasing function.
  void set_applicability_condition(index_type c) const;

  //! \brief
  //! Returns the rank of the functional equation, i. e., the number of terms
  //! of the form \f$ a x(n/b) \f$ where \f$ b \f$ is a rational number
  //! larger than one.
  index_type rank() const;

  //! \brief
  //! If the divisor \f$ b \f$ of the homogeneous term \f$ a x(n/b) \f$ in
  //! \p *this is different from \f$ 2 \f$, sets the function \f$ Sc() \f$
  //! in the following way:
  //! \f[
  //!   Sc(n, b) = \lfloor
  //!                \frac{n}{b^{\lfloor \frac{\log n}{\log b} \rfloor}}
  //!              \rfloor.
  //! \f]
  //! Set the definition of the function \f$ Sc() \f$ to the empty string
  //! if \f$ b = 2 \f$.
  void set_definition_Sc() const;

  // Method to access to private data of `Non_Linear_Info'.

  //! \brief
  //! In the case which the system is able to rewrite the non-linear
  //! recurrence \p *this in linear, this method returns the linear
  //! recurrence computed (in order to know the cases of rewritable
  //! non-linear recurrences see the function
  //! <CODE>rewrite_non_linear_recurrence()</CODE>).
  const Recurrence& associated_linear_rec() const;

  //! \brief
  //! In the case which the system is able to rewrite the non-linear
  //! recurrence \p *this in linear, this method returns the linear
  //! recurrence computed (in order to know the cases of rewritable
  //! non-linear recurrences see the function
  //! <CODE>rewrite_non_linear_recurrence()</CODE>).
  Recurrence& associated_linear_rec();

  //! \brief
  //! In the case of simple non-linear recurrence of the form
  //! \f$ c x(n-1)^{\alpha} \f$ contains \f$ c \f$; in all the
  //! other cases contains \f$ 0 \f$.
  const Number& coeff_simple_non_linear_rec() const;

  //! \brief
  //! In the case of simple non-linear recurrence of the form
  //! \f$ c x(n-1)^{\alpha} \f$ contains \f$ c \f$; in all the
  //! other cases contains \f$ 0 \f$.
  Number& coeff_simple_non_linear_rec();

  //! \brief
  //! Contains the value that will be the logarithm's base or the
  //! exponential's base used in the rewriting of the non-linear recurrence
  //! in the correspondent linear recurrence.
  const Expr& base_exp_log() const;

  //! \brief
  //! Contains the value that will be the logarithm's base or the
  //! exponential's base used in the rewriting of the non-linear recurrence
  //! in the correspondent linear recurrence.
  Expr& base_exp_log();

  //! \brief
  //! Stores the symbols associated to the possible negative numbers
  //! that will be the arguments of the logarithms.
  const std::vector<Symbol>& auxiliary_symbols() const;

  //! \brief
  //! Stores the symbols associated to the possible negative numbers
  //! that will be the arguments of the logarithms.
  std::vector<Symbol>& auxiliary_symbols();


  // Method to access to private data of `Weighted_Average_Info'.

  //! \brief
  //! In the case which the system is able to rewrite the weighted-average
  //! recurrence \p *this in first order recurrence, this method returns
  //! the first order recurrence computed (in order to know the cases of
  //! rewritable weighted-average recurrences see the function
  //! <CODE>rewrite_weighted_average_recurrence()</CODE>).
  const Recurrence& associated_first_order_rec() const;

  //! \brief
  //! In the case which the system is able to rewrite the weighted-average
  //! recurrence \p *this in first order recurrence, this method returns
  //! the first order recurrence computed (in order to know the cases of
  //! rewritable weighted-average recurrence see the function
  //! <CODE>rewrite_weighted_average_recurrence()</CODE>).
  Recurrence& associated_first_order_rec();

  //! \brief
  //! Returns the factor \f$ f(n) \f$ of the weighted-average recurrence
  //! \f[
  //!   T(n) = f(n) \sum_{k=0}^{n-1} T(k) + g(n).
  //! \f]
  const Expr& weight() const;

  //! \brief
  //! Returns the factor \f$ f(n) \f$ of the weighted-average recurrence
  //! \f[
  //!   T(n) = f(n) \sum_{k=0}^{n-1} T(k) + g(n).
  //! \f]
  Expr& weight();

  //! \brief
  //! Returns the factor \f$ f(n) \f$ of the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ n_0 \in \Nset \cup \{ 0 \} \f$ and
  //! \f$ u(n) \in \{ n-1, n \} \f$, which has been rewritten
  //! in the form of weighted-average recurrence.
  const Expr& original_weight() const;

  //! \brief
  //! Returns the factor \f$ f(n) \f$ of the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ n_0 \in \Nset \cup \{ 0 \} \f$ and
  //! \f$ u(n) \in \{ n-1, n \} \f$, which has been rewritten
  //! in the form of weighted-average recurrence.
  Expr& original_weight();

  //! \brief
  //! Returns the inhomogeneous term \f$ g(n) \f$ of the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ n_0 \in \Nset \cup \{ 0 \} \f$ and
  //! \f$ u(n) \in \{ n-1, n \} \f$, which has been rewritten
  //! in the form of weighted-average recurrence.
  const Expr& original_inhomogeneous() const;

  //! \brief
  //! Returns the inhomogeneous term \f$ g(n) \f$ of the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ n_0 \in \Nset \cup \{ 0 \} \f$ and
  //! \f$ u(n) \in \{ n-1, n \} \f$, which has been rewritten
  //! in the form of weighted-average recurrence.
  Expr& original_inhomogeneous();

  //! \brief
  //! Returns the lower limit \f$ n_0 \in \Nset \cup \{ 0 \} \f$
  //! of the sum in the recurrence
  //! \f[
  //!   x(n) = f(n) \sum_{k=n_0}^{u(n)} x(k) + g(n),
  //! \f]
  //! where \f$ u(n) \in \{ n-1, n \} \f$, which has been rewritten
  //! in the form of weighted-average recurrence.
  unsigned int lower_limit() const;

  //! \brief
  //! Sets with \p weight, \p inhomogeneous, \p lower and \p upper
  //! the elements of <CODE>Components_Original_Rec</CODE>.
  void set_original_rhs(const Expr& weight, const Expr& inhomogeneous,
			unsigned int lower, const Expr& upper) const;
  
  //! Contains the \ref first_valid_index "first valid index".
  mutable index_type first_valid_index;

  //! \brief
  //! Contains the
  //! \ref first_valid_index_for_solution "first valid index for solution".
  mutable index_type first_valid_index_for_solution_;

  //! Sets to \p i <CODE>first_valid_index</CODE>.
  void set_first_valid_index(index_type i) const;

  //! Sets <CODE>first_valid_index_for_solution_</CODE>.
  /*!
    Computes
    \ref first_valid_index_for_solution "first_valid_index_for_solution"
    knowing <CODE>first_valid_index</CODE> and the type of \p *this.
    Once it has been computed, assign it to
    <CODE>first_valid_index_for_solution_</CODE>.
  */
  void set_first_valid_index_for_solution() const;

  //! \brief
  //! Stores the value, if it exists, of the exact solution
  //! (with symbolic initial conditions) and other informations about it.
  mutable Cached_Expr exact_solution_;

  //! \brief
  //! Stores the value, if it exists, of the lower bound
  //! (with symbolic initial conditions) and other informations about it.
  mutable Cached_Expr lower_bound_;

  //! \brief
  //! Stores the value, if it exists, of the upper bound
  //! (with symbolic initial conditions) and other informations about it.
  mutable Cached_Expr upper_bound_;

  //! \brief
  //! Stores the value, if it exists, of the exact solution
  //! evaluated on the initial conditions stores in the map
  //! <CODE>initial _conditions</CODE> and other informations about it.
  mutable Cached_Expr evaluated_exact_solution_;

  //! \brief
  //! Stores the value, if it exists, of the lower bound
  //! evaluated on the initial conditions stores in the map
  //! <CODE>initial _conditions</CODE> and other informations about it.
  mutable Cached_Expr evaluated_lower_bound_;

  //! \brief
  //! Stores the value, if it exists, of the upper bound
  //! evaluated on the initial conditions stores in the map
  //! <CODE>initial _conditions</CODE> and other informations about it.
  mutable Cached_Expr evaluated_upper_bound_;

  //! \brief
  //! If \p tried_to_compute_exact_solution is true then the system has
  //! already tried to compute the exact solution.
  mutable bool tried_to_compute_exact_solution;

  mutable Blackboard blackboard;

  std::map<index_type, Expr> initial_conditions_;

private:
  static Classifier_Status
  compute_order(const Number& decrement, index_type& order,
		unsigned long& index, unsigned long max_size);
  static Expr
  write_expanded_solution(const Recurrence& rec,
			  unsigned int gcd_among_decrements);

  //! \brief
  //! This function must have access to the private data
  //! <CODE>blackboard</CODE> of the class <CODE>Recurrence</CODE>.
  friend Expr Cached_Expr::
  replace_system_generated_symbols(const Recurrence& rec) const;
};

} // namespace Parma_Recurrence_Relation_Solver

//#include "Recurrence.inlines.hh"

#endif // !defined(PURRS_Recurrence_defs_hh)

