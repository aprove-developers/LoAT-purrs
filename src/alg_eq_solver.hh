
//! Let \p p be a polynomial with integer coefficients in \p x and
//! \p roots be a (possibly non-empty) vector.
//! This function appends to \p roots some or all the roots of \p p.
//! Let \f$n\f$ be the degree of \p p and let \f$h\f$ and \f$k\f$
//! be the value of <CODE>roots.size()</CODE> on entry and on exit
//! to and from find_roots(), respectively.
//! If \f$n = k-h\f$, then the positions \f$h, h+1, \ldots, k-1\f$
//! of \p roots contain <EM>all</EM> the (possibly complex) roots of \p p
//! and the function returns <CODE>true</CODE>.
//! If \f$n \neq k-h\f$, then the positions \f$h, h+1, \ldots, k-1\f$
//! of \p roots contain <EM>some</EM> of the roots of \p p and the function
//! returns <CODE>true</CODE> if it was possible to prove that <EM>all</EM>
//! the roots of \p p of maximal modulus are among those inserted
//! into \p roots;
//! the function returns <CODE>false</CODE> otherwise.
//! The parameter \p all_distinct is set to <CODE>true</CODE> if it was
//! possible to prove that all the roots of \p p are distinct;
//! \p all_distinct is set to <CODE>false</CODE> otherwise.
bool
find_roots(const GExpr& p, const GSymbol& x,
	   vector<GExpr>& roots, bool& all_distinct);
