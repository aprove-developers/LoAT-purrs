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
#include "Cached_Expr.defs.hh"
#include "Non_Linear_Info.defs.hh"
#include "Infinite_Order_Info.defs.hh"
#include "Blackboard.defs.hh"
#include "Finite_Order_Info.defs.hh"
#include "Functional_Equation_Info.defs.hh"
#include "Expr.defs.hh"
#include "alg_eq_solver.hh"
#include <map>
#include <iosfwd>

namespace Parma_Recurrence_Relation_Solver {

//! \brief
//! An unsigned integral type for representing different types of
//! indexes of the recurrence (e.g. the order of the recurrence).
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
  void replace_recurrence(unsigned int k, const Expr& e);

  //! Returns a new symbol \f$ z \f$ and records the equation \f$ z = e \f$.
  Symbol insert_auxiliary_definition(const Expr& e) const;

  //! \brief
  //! Returns an expression obtained from \p e by substituting the symbols
  //! that are also on the blackboard with their definition.
  Expr substitute_auxiliary_definitions(const Expr& e) const;

  //! \brief
  //! Replaces the values in the \f$ k \f$-th position of the map
  //! <CODE>initial_conditions</CODE> with the expression \p e. 
  void replace_initial_condition(unsigned int k, const Expr& e);

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

  //! \brief
  //! If in the map <CODE>initial_conditions</CODE> there is the
  //! expression \f$ e \f$ correspondent to \p k then returns \f$ e \f$;
  //! returns \f$ x(k) \f$ otherwise.
  Expr get_initial_condition(unsigned int k) const;

public:
#ifdef PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS
  //! Kinds of bounds for the solution that is possible to compute.
#endif // PURRS_DOXYGEN_INCLUDE_IMPLEMENTATION_DETAILS
  enum Bound {
    /*! An upper bound. */
    UPPER,
    
    /*! A lower bound. */
    LOWER
  };

  //! The possible states of the recurrence.
  enum Solver_Status {
    /*!
      Solution was successful.
    */
    SUCCESS,

    /*!
      The recurrence is indeterminate, hence it has infinite solutions.
    */
    INDETERMINATE_RECURRENCE,

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
      The recurrence is not solvable.
    */
    UNSOLVABLE_RECURRENCE,

    /*!
      The recurrence does not have any sense.
    */
    MALFORMED_RECURRENCE,

    /*!
      The recurrence is not well-defined.
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
      The system can prove that the recurrence has been successfully
      solved or approximated.
    */
    PROVABLY_CORRECT,

    /*!
      The system can prove that the solution, or the approximation, of
      \p *this, is wrong.
    */
    PROVABLY_INCORRECT,

    /*!
      The system can not prove if the solution, or the approximation, of
      \p *this, is correct or not. 
    */
    INCONCLUSIVE_VERIFICATION
  };

  //! \brief
  //! Tries to solve \p *this exactly and returns <CODE>SUCCESS</CODE>
  //! if the system finds the exact solution.
  Solver_Status compute_exact_solution() const;

  //! Gets the exact solution and puts it in \p e
  void exact_solution(Expr& e) const;

  //! \brief
  //! Tries to get lower bound for \p *this and returns <CODE>SUCCESS</CODE>
  //! if the system finds the lower bound.
  Solver_Status compute_lower_bound() const;

  //! Gets the lower bound for the solution and puts it in \p e.
  void lower_bound(Expr& e) const;

  //! \brief
  //! Tries to get upper bound for \p *this and returns <CODE>SUCCESS</CODE>
  //! if the system finds the upper bound.
  Solver_Status compute_upper_bound() const;

  //! Gets the upper bound for the solution and puts it in \p e.
  void upper_bound(Expr& e) const;

  Expr approximated_solution() const;

  //! \brief
  //! Verifies the exact solution of \p *this; note that the system
  //! supposes the exact solution already has been computed.
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system proved
  //! that the recurrence \p *this has been successfully solved.
  //! Returns <CODE>PROVABLY_INCORRECT</CODE> if the system proved
  //! that the solution of the recurrence \p *this is wrong.
  //! Returns <CODE>INCONCLUSIVE_VERIFICATION</CODE> when the system is
  //! not able to prove if the solution of the recurrence \p *this
  //! is correct or not.
  /*!
    Since the system supposes the exact solution already has been
    computed, the user must call this method after to have checked
    the system has successfully computed the exact solution:
    this is possibile invoking the method
    <CODE>compute_exact_solution</CODE> and verifying it has returned
    <CODE>SUCCESS</CODE>.

    \exception std::logic_error thrown if this method is called
                                when no exact solution was computed.
  */
  Verify_Status verify_exact_solution() const;
 
  // @@@
  //! \brief
  //! Verifies the lower bound for \p *this; note that the system
  //! supposes the lower bound already has been computed.
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system proved
  //! that the recurrence \p *this has been successfully approximated
  //! from below.
  //! Returns <CODE>PROVABLY_INCORRECT</CODE> if the system proved
  //! that the lower bound for the solution of the recurrence \p *this
  //! is wrong.
  //! Returns <CODE>INCONCLUSIVE_VERIFICATION</CODE> when the system is
  //! not able to prove if the lower bound for the solution of the
  //! recurrence \p *this is correct or not.
  /*!
    \exception std::logic_error thrown if this method is called
                                when no lower bound was computed.
  */
  Verify_Status verify_lower_bound() const;

  // @@@
  //! \brief
  //! Verifies the upper bound for \p *this; note that the system
  //! supposes the upper bound already has been computed.
  //! Returns <CODE>PROVABLY_CORRECT</CODE> if the system proved
  //! that the recurrence \p *this has been successfully approximated
  //! from above.
  //! Returns <CODE>PROVABLY_INCORRECT</CODE> if the system proved
  //! that the upper bound for the solution of the recurrence \p *this
  //! is wrong.
  //! Returns <CODE>INCONCLUSIVE_VERIFICATION</CODE> when the system is
  //! not able to prove if the upper bound for the solution of the
  //! recurrence \p *this is correct or not.
  /*!
    \exception std::logic_error thrown if this method is called
                                when no upper bound was computed.
  */
  Verify_Status verify_upper_bound() const;

  //! \brief
  //! Returns <CODE>false</CODE> if there are no undefined initial conditions;
  //! otherwise returns <CODE>true</CODE> and adds the corresponding
  //! indexes to \p undefined.
  bool undefined_initial_conditions(std::set<unsigned int>& undefined) const;

  //! The index of the recurrence.
  static const Symbol& n;

  //! Dumps all the data members of \p *this onto \p s.
  void dump(std::ostream& s) const;

private:
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
  //! Returns <CODE>SUCCESS</CODE> if the system is able to solve the
  //! recurrence of infinite order
  //! \f$ x(n) = f(n) \sum_{k=n_0}^n-1 x(k) + g(n) \f$, where
  //! \f$ f(n) \f$ is stored in \p weight and  \f$ g(n) \f$ is stored
  //! in \p inhomogeneous: in this case the solution is returned in
  //! \p solution.
  Solver_Status
  solve_new_infinite_order_rec(const Expr& weight, const Expr& inhomogeneous,
			       index_type first_valid_index,
			       Expr& solution) const;

  //! \brief
  //! Solves the recurrence of infinite order \p *this transforming it
  //! in a linear recurrence of finite order stored in a new object
  //! <CODE>Recurrence</CODE>.
  Solver_Status compute_infinite_order_recurrence(Expr& solution) const;

  //! \brief
  //! Let \p solution_or_bound be the expression that represent the
  //! solution or the bound computed for the recurrence \p *this.
  //! This function substitutes eventual initial conditions specified
  //! by the user shifting the solution or the bound if necessary.
  Expr substitute_i_c_shifting(const Expr& solution_or_bound) const;

  //! Classifies the recurrence \p *this sliding it recursively.
  Solver_Status classify() const;

  //! \brief
  //! Classifies the recurrence \p *this calling the method
  //! <CODE>classify()</CODE>. 
  Solver_Status classify_and_catch_special_cases() const;

  //! \brief
  //! Analyzes the \f$ i \f$-th addend of the right hand side \p rhs
  //! of the recurrence \p *this.
  Solver_Status classification_summand(const Expr& rhs, const Expr& r,
				       Expr& e, index_type& order,
				       std::vector<Expr>& coefficients,
				       int& gcd_among_decrements,
				       int num_term,
				       std::map<Number, Expr>&
				       homogeneous_terms) const;

  //! \brief
  //! Solves the linear recurrence of finite order.
  //! The recurrence can have constant or variable coefficients.
  Solver_Status solve_linear_finite_order() const;

  //! Solves the linear recurrence of first order with constant coefficient.
  Expr solve_constant_coeff_order_1(const std::vector<Polynomial_Root>&
				    roots) const;

  //! Solves the linear recurrence of first order with variable coefficient. 
  Solver_Status
  solve_variable_coeff_order_1(const std::vector<Expr>& coefficients) const;

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

  //! @@@
  Verify_Status verify_bound(Bound kind_of_bound) const;

  //! Computes the lower bound for the functional equation. 
  Solver_Status approximate_functional_equation_lower() const;

  //! Computes the upper bound for the functional equation. 
  Solver_Status approximate_functional_equation_upper() const;

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
  //! It is <CODE>true</CODE> when the recurrence has been deduced from
  //! a non-linear recurence; it is <CODE>false</CODE> in all the other cases.
  bool come_from_non_linear_rec;

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
      Linear recurrence of infinite order. 
    */
    LINEAR_INFINITE_ORDER,
    
    /*!
      Special recurrence of the form \f$ x(n) = a x(n / b) + d n^e \f$ . 
    */
    FUNCTIONAL_EQUATION
  };

  mutable Type type_;

  //! \brief
  //! It is <CODE>true</CODE> if the recurrence already it has been
  //! classified; it is <CODE>false</CODE> otherwise.
  bool is_classified;

  mutable Finite_Order_Info* finite_order_p;
  mutable Functional_Equation_Info* functional_eq_p;
  mutable Non_Linear_Info* non_linear_p;
  mutable Infinite_Order_Info* infinite_order_p;

  //! Returns <CODE>type_</CODE>.
  Type type() const;

  //! Returns <CODE>type_</CODE>.
  Type& type();

  //! Sets <CODE>type_</CODE> with \p t.
  void set_type(const Type& t) const;

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

  //! \brief
  //! Returns <CODE>true</CODE> if the recurrence is linear
  //! of infinite order; returns <CODE>false</CODE> otherwise.
  bool is_linear_infinite_order() const;

  //! Sets <CODE>type_recurrence =  LINEAR_INFINITE_ORDER</CODE>.
  void set_linear_infinite_order() const;


  // Methods to access to private data of `Finite_Order_Info'.

  //! Returns the order of the finite order recurrence.
  index_type order() const;

  //! \brief
  //! Returns the smallest positive integer for which the finite order
  //! recurrence is well-defined: the initial conditions will start from it.
  index_type first_valid_index() const;

  //! \brief
  //! Sets to \p i_c is the smallest positive integer for which the finite
  //! order recurrence is well-defined: the initial conditions will start
  //! from it. 
  void set_first_valid_index(index_type i_c) const;

  //! Returns the coefficients of the linear finite order recurrence.
  const std::vector<Expr>& coefficients() const;

  //! Returns the coefficients of the linear finite order recurrence.
  std::vector<Expr>& coefficients();

  //! \brief
  //! Returns the greatest common divisor among the decrements \f$ k \f$
  //! of the terms \f$ x(n-k) \f$ of a linear finite order recurrence.
  //! Returns \f$ 0 \f$ if the order of the recurrence is \f$ 0 \f$.
  unsigned int gcd_among_decrements() const;

  //! \brief
  //! Returns the expression \f$ \prod_{i}^n a(k)\f$,
  //! where \f$ i \f$ is \p first_valid_index() and \f$ a(n) \f$
  //! is the coefficient of the recurrence of first order with variable
  //! coefficient.
  Expr& product_factor();

  //! \brief
  //! Returns the expression \f$ \prod_{i}^n a(k)\f$,
  //! where \f$ i \f$ is \p first_valid_index() and \f$ a(n) \f$
  //! is the coefficient of the recurrence of first order with variable
  //! coefficient.
  Expr product_factor() const;

  //! \brief
  //! Sets to \p x the expression \f$ \prod_{i}^n a(k)\f$,
  //! where \f$ i \f$ is \p first_valid_index() and \f$ a(n) \f$
  //! is the coefficient of the recurrence of first order with variable
  //! coefficient.
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


  // Method to access to private data of `Non_Linear_Info'.

  //! \brief
  //! If the non-linear recurrence is rewritable in a linear recurrence
  //! then this data contains the right hand side of the linear recurrence
  //! associated to the original non-linear recurrence.
  Expr rhs_transformed_in_linear() const;

  //! \brief
  //! If the non-linear recurrence is rewritable in a linear recurrence
  //! then this data contains the right hand side of the linear recurrence
  //! associated to the original non-linear recurrence. 
  Expr& rhs_transformed_in_linear();

  //! \brief
  //! In the case of simple non-linear recurrence of the form
  //! \f$ c x(n-1)^{\alpha} \f$ contains \f$ c \f$; in all the
  //! other cases contains \f$ 0 \f$.
  Number coeff_simple_non_linear_rec() const;

  //! \brief
  //! In the case of simple non-linear recurrence of the form
  //! \f$ c x(n-1)^{\alpha} \f$ contains \f$ c \f$; in all the
  //! other cases contains \f$ 0 \f$.
  Number& coeff_simple_non_linear_rec();

  //! \brief
  //! Contains the value that will be the logarithm's base or the
  //! exponential's base used in the rewriting of the non-linear recurrence
  //! in the correspondent linear recurrence.
  Expr base_exp_log() const;

  //! \brief
  //! Contains the value that will be the logarithm's base or the
  //! exponential's base used in the rewriting of the non-linear recurrence
  //! in the correspondent linear recurrence.
  Expr& base_exp_log();

  //! \brief
  //! Stores the symbols associated to the eventual negative numbers
  //! that will be the arguments of the logarithms.
  const std::vector<Symbol>& auxiliary_symbols() const;

  //! \brief
  //! Stores the symbols associated to the eventual negative numbers
  //! that will be the arguments of the logarithms.
  std::vector<Symbol>& auxiliary_symbols();

  //! \brief
  //! Returns the order of the finite order recurrence associated to
  //! the non-linear recurrence .
  index_type order_if_linear() const;

  //! \brief
  //! Sets to \p x the order of the linear recurrence associated to
  //! the non linear. 
  void set_order_if_linear(index_type x) const;

  //! \brief
  //! When the non-linear recurrence is rewritable in a linear recurrence
  //! of finite order then this method returns the smallest positive
  //! integer for which the finite order recurrence is well-defined:
  //! the initial conditions will start from it.
  index_type first_valid_index_if_linear() const;

  //! \brief
  //! When the non-linear recurrence is rewritable in a linear recurrence
  //! of finite order then this method sets to \p i_c is the smallest
  //! positive integer for which the finite order recurrence is well-defined:
  //! the initial conditions will start from it. 
  void set_first_valid_index_if_linear(index_type i_c) const;


  // Method to access to private data of `Infinite_Order_Info'.

  //! \brief
  //! If the infinite order recurrence is rewritable in a first order linear
  //! recurrence then this last recurrence is returned from this method.
  Expr rhs_transformed_in_first_order() const;

  //! \brief
  //! If the infinite order recurrence is rewritable in a first order linear
  //! recurrence then this last recurrence is returned from this method.
  Expr& rhs_transformed_in_first_order();

  //! \brief
  //! If the infinite order recurrence is rewritable in a first order linear
  //! recurrence then the coefficient of this last recurrence is returned
  //! from this method.
  Expr coeff_first_order() const;

  //! \brief
  //! If the infinite order recurrence is rewritable in a first order linear
  //! recurrence then the coefficient of this last recurrence is returned
  //! from this method.
  Expr& coeff_first_order();

  //! \brief
  //! If the infinite order recurrence is rewritable in a first order linear
  //! recurrence then the non-homogeneous part of this last recurrence is
  //! returned from this method.
  Expr inhomog_first_order() const;

  //! \brief
  //! If the infinite order recurrence is rewritable in a first order linear
  //! recurrence then the non-homogeneous part of this last recurrence is
  //! returned from this method.
  Expr& inhomog_first_order();

  //! \brief
  //! Returns the factor \f$ f(n) \f$ of the infinite order recurrence
  //! \f[
  //!   T(n) = f(n) \sum_{k=n_0}^{u(n)} T(k) + g(n).
  //! \f]
  Expr weight_inf_order() const;

  //! \brief
  //! Returns the factor \f$ f(n) \f$ of the infinite order recurrence
  //! \f[
  //!   T(n) = f(n) \sum_{k=n_0}^{u(n)} T(k) + g(n).
  //! \f]
  Expr& weight_inf_order();

  //! \brief
  //! Returns the lower bound \f$ n_0 \f$ of the infinite order recurrence
  //! \f[
  //!   T(n) = f(n) \sum_{k=n_0}^{u(n)} T(k) + g(n).
  //! \f]
  unsigned int lower_bound_sum() const;

  //! \brief
  //! Returns the upper bound \f$ u(n) \f$ of the infinite order recurrence
  //! \f[
  //!   T(n) = f(n) \sum_{k=n_0}^{u(n)} T(k) + g(n).
  //! \f]
  Expr upper_bound_sum() const;

  //! \brief
  //! Returns the upper bound \f$ u(n) \f$ of the infinite order recurrence
  //! \f[
  //!   T(n) = f(n) \sum_{k=n_0}^{u(n)} T(k) + g(n).
  //! \f]
  Expr& upper_bound_sum();

  //! \brief
  //! Stores the smallest positive integer for which the infinite
  //! order recurrence is well-defined: the initial conditions will
  //! start from it.
  index_type first_valid_index_inf_order() const;

  //! \brief
  //! Stores the smallest positive integer for which the infinite
  //! order recurrence is well-defined: the initial conditions will
  //! start from it.
  void set_first_valid_index_inf_order(index_type i_c) const;



  mutable Cached_Expr exact_solution_;
  mutable Cached_Expr lower_bound_;
  mutable Cached_Expr upper_bound_;

  //! \brief
  //! If \p tested_exact_solution is true then the system have already
  //! try to compute the exact solution.
  bool tested_exact_solution;

  mutable Blackboard blackboard;

  std::map<unsigned int, Expr> initial_conditions;

private:
  static Solver_Status
  compute_order(const Number& decrement, index_type& order,
		unsigned long& index, unsigned long max_size);
  static Expr
  write_expanded_solution(const Recurrence& rec,
			  unsigned int gcd_among_decrements);

  //! \brief
  //! This function must have access to the private data
  //! <CODE>blackboard</CODE>of the class <CODE>Recurrence</CODE>.
  friend Expr Cached_Expr::
  replace_system_generated_symbols(const Recurrence& rec) const;
};

} // namespace Parma_Recurrence_Relation_Solver

#include "Recurrence.inlines.hh"

#endif // !defined(PURRS_Recurrence_defs_hh)

