/* This is the main object PURRS operates upon: a recurrence.
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

#ifndef PURRS_Recurrence_defs_hh
#define PURRS_Recurrence_defs_hh 1

#include "Recurrence.types.hh"
#include "Blackboard.defs.hh"
#include "Finite_Order_Info.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Expr.defs.hh"
#include "alg_eq_solver.hh"
#include <map>
#include <iosfwd>

namespace Parma_Recurrence_Relation_Solver {

//! \brief
//! Assuming that \p x and \p y represent sequences of reals,
//! returns <CODE>true</CODE> if it can be proved that all
//! the recurrences represented by \p x are less than all
//! the recurrences represented by \p y.
/*!
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


/*!
  An object of this class abstracts a (possibly infinite) set of
  sequences of real numbers.  Formally, an object of class
  Recurrence is an element of \f$ \wp(\Nset \to \Rset) \f$.  Elements
  of \f$ \wp(\Nset \to \Rset) \f$ can be specified by imposing a
  recurrence relation (whence the class' name) that the sequences must
  satisy, by imposing some initial conditions for such recurrences,
  and by computing approximations of set union and intersection.
*/
class Recurrence {
public:
  //! Builds the singleton satisfying \f$ x(n) = 0 \f$.
  Recurrence();

  //! Builds the set of sequences satisfying \f$ x(n) = e \f$.
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

  //! WRITEME
  void replace_recurrence(const Expr& e);

  //! \brief
  //! Sets to \f$ e \f$ the right-hand side of the recurrence
  //! of index \f$ k \f$.  The system of recurrences will then
  //! include \f$ x_k(n) = e \f$.
  void replace_recurrence(unsigned k, const Expr& e);

  //! Returns a new symbol \f$ z \f$ and records the equation \f$ z = e \f$.
  Symbol insert_auxiliary_definition(const Expr& e) const;

  //! \brief
  //! Substitutes eventually auxiliary definitions contained in
  //! \p e with their original values stored in the blackboard.
  Expr substitute_auxiliary_definitions(const Expr& e) const;

  //! Checks if all the invariants are satisfied.
  /*!
    The check is performed so as to intrude as little as possible.
    In case invariants are violated error messages are written on
    <CODE>std::cerr</CODE>.
  */
  bool OK() const; 

private:
  //! \brief
  //! Returns the right-hand side of the auxiliary equation \f$ z = e \f$,
  //! if such an auxiliary equation exists;
  //! returns the expression \f$ z \f$ otherwise.
  Expr get_auxiliary_definition(const Symbol& z) const;

  //! Sets to \p e the <CODE>inhomogeneous_term</CODE>.
  void set_inhomogeneous_term(const Expr& e) const;

public:
  enum Solver_Status {
    /*!
      Solution was successful.
    */
    SUCCESS,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is not an integer.
    */
    HAS_NON_INTEGER_DECREMENT,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is a negative integer.
    */
    HAS_NEGATIVE_DECREMENT,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n-k)</CODE> where <CODE>k</CODE> is too big to be handled
      by the standard solution techniques.
    */
    HAS_HUGE_DECREMENT,

    /*!
      The right-hand side of the recurrence contains at least an occurrence
      of <CODE>x(n)</CODE>.
    */
    HAS_NULL_DECREMENT,

    /*!
      The recurrence is not linear.
    */
    NON_LINEAR_RECURRENCE,

    /*!
      The recurrence is not solvable.
    */
    UNSOLVABLE_RECURRENCE,

    /*!
      Catchall: the recurrence is generically too complex for the solver.
    */
    TOO_COMPLEX
  };


  Solver_Status solve() const;
  Expr exact_solution() const;
  Expr approximated_solution() const;

  enum Verify_Status {
    /*!
      The system can prove that the recurrence has been successfully
      solved.
    */
    PROVABLY_CORRECT,

    /*!
      The system can prove that the solution of \p *this is wrong.
    */
    PROVABLY_INCORRECT,

    /*!
      The system can not prove if the solution of \p *this is correct
      or not. 
    */
    INCONCLUSIVE_VERIFICATION
  };

  //! \brief
  //! Returns <CODE>true</CODE> if the <CODE>solution</CODE> is right.
  //! Returns <CODE>false</CODE> if the solution is wrong or the
  //! simplifications are failed, i.e., the solution is right but
  //! this function does not realize.
  Verify_Status verify_solution() const;

  //! The index of the recurrence.
  static const Symbol& n;

  //! Dumps all the data members of \p *this onto \p s.
  void dump(std::ostream& s) const;

private:
  Solver_Status solve_easy_cases() const;
  Solver_Status solve_try_hard() const;
  Solver_Status classification_recurrence(const Expr& rhs,
					  int& gcd_among_decrements) const;
  Solver_Status classification_summand(const Expr& r, Expr& e,
				       std::vector<Expr>& coefficients,
				       unsigned int& order,
				       int& gcd_among_decrements,
				       int num_term,
				       Expr& coefficient,
				       unsigned& divisor_arg) const;
  void add_initial_conditions(const Expr& g_n,
			      const std::vector<Number>& coefficients) const;
  Solver_Status solve_linear_finite_order(int gcd_among_decrements) const;
  Solver_Status
  solve_constant_coeff_order_1(const std::vector<Polynomial_Root>&
			       roots) const;
  Solver_Status solve_variable_coeff_order_1(const Expr& coefficient) const;
  Expr
  solve_constant_coeff_order_2(Expr& g_n, bool all_distinct,
			       const std::vector<Number>& coefficients,
			       const std::vector<Polynomial_Root>& roots) const;
  Expr
  solve_constant_coeff_order_k(Expr& g_n, bool all_distinct,
			       const std::vector<Number>& coefficients,
			       const std::vector<Polynomial_Root>& roots) const;
  Solver_Status approximate_functional_equation() const;

  //! \brief
  //! Holds the right-hand side of the global recurrence to be solved.
  //! This may have been set directly by the constructor or it may be the
  //! result of transforming a system into a single recurrence.
  //! The global recurrence is thus of the form
  //! <CODE>x(n) = recurrence_rhs</CODE>.
  mutable Expr recurrence_rhs;

  //! \brief
  //! When is applied the order reduction stores the value of the
  //! <CODE>recurrence_rhs</CODE> before to solve the reduced recurrence;
  //! stores 0 otherwise.
  mutable Expr old_recurrence_rhs;

  //! \brief
  //! When is applied the order reduction stores the greatest common
  //! divisor among the decrements <CODE>d</CODE> of the terms
  //! <CODE>x(n-d)</CODE> present in the right hand side of the recurrence
  //! before to solve the reduced recurrence;
  //! stores 0 otherwise.
  mutable unsigned gcd_decrements_old_rhs;

  //! \brief
  //! When is applied the order reduction stores the solution of the
  //! reduced order recurrence; stores 0 otherwise.
  mutable Expr solution_order_reduced;

  //! \brief
  //! It is <CODE>true</CODE> if the recurrence has been rewritten, i. e.,
  //! in the case there are null or negative decrements or if has been applied
  //! the order reduction; it is <CODE>false</CODE> in all the other cases.
  mutable bool recurrence_rhs_rewritten;

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
  std::map<unsigned, Expr> system_rhs;

  enum Type {
    /*!
      The type of the recurrence is unknown.
    */
    UNKNOWN,

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
      Linear recurrence of infinite order. 
    */
    LINEAR_INFINITE_ORDER,
    
    /*!
      Special recurrence of the form \f$ x(n) = a x(n / b) + d n^e \f$ . 
    */
    FUNCTIONAL_EQUATION
  };

  mutable Type type;
  mutable Finite_Order_Info* finite_order_p;
  mutable Functional_Equation_Info* functional_eq_p;

  //! \brief
  //! Returns <CODE>true</CODE> if the recurrence's type is unknown;
  //! returns <CODE>false</CODE> otherwise.
  bool is_unknown() const;

  //! \brief
  //! Returns <CODE>true</CODE> if the recurrence is a special case,
  //! i. e., is of order zero; returns <CODE>false</CODE> otherwise.
  bool is_order_zero() const;

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

  //! Returns the order of the finite order recurrence.
  unsigned int order() const;

  //! Returns the order of the finite order recurrence.
  unsigned int& order();

  //! \brief
  //! Returns the smallest positive integer for which the finite order
  //! recurrence is well-defined: the initial conditions will start from it.
  unsigned first_initial_condition() const;

  //! \brief
  //! Returns the smallest positive integer for which the finite order
  //! recurrence is well-defined: the initial conditions will start from it.
  unsigned& first_initial_condition();

  //! \brief
  //! \p i_c is the smallest positive integer for which the finite order
  //! recurrence is well-defined: the initial conditions will start from it. 
  void set_first_initial_condition(unsigned i_c) const;

  //! Returns the coefficients of the finite order recurrence.
  const std::vector<Expr>& coefficients() const;

  //! Returns the coefficients of the finite order recurrence.
  std::vector<Expr>& coefficients();

  //! \brief
  //! Returns the coefficient \f$ a \f$ of the functional equation
  //! \f$ x_n = a x_{n/b} + p(n).
  Expr coefficient() const;

  //! \brief
  //! Returns the coefficient \f$ a \f$ of the functional equation
  //! \f$ x_n = a x_{n/b} + p(n).
  Expr& coefficient();

  //! \brief
  //! Returns \f$ b \f$, the divisor of \f$ n \f$ in the functional equation
  //! \f$ x_n = a x_{n/b} + p(n)..
  unsigned divisor_arg() const;

  //! \brief
  //! Returns \f$ b \f$, the divisor of \f$ n \f$ in the functional equation
  //! \f$ x_n = a x_{n/b} + p(n)..
  unsigned& divisor_arg();

  mutable bool solved;

  mutable Expr solution;

  mutable Blackboard blackboard;

private:
  static Solver_Status
  check_powers_and_functions(const Expr& e);
  static Solver_Status
  find_non_linear_recurrence(const Expr& e);
  static Solver_Status
  compute_order(const Number& decrement, unsigned int& order,
		unsigned long& index, unsigned long max_size);

#if 0
    Expr poly_char;
    Expr symb_solution;

public:
  insert_exact_solution(const Expr& e);
  insert_side_condition(int n, const Expr& e);
  // Restituisce le condizioni iniziali non assegnate.
  std::set<unsigned int> undefined_side_conditions();
  complex_interval approximate(const std::vector<Expr> side_conditions,
			       int n);
//    Expr big_o(...);
//    Expr big_omega(...);
#endif
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Recurrence.inlines.hh"

#endif // !defined(PURRS_Recurrence_defs_hh)

