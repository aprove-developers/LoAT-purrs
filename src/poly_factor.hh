#ifndef PURRS_poly_factor_hh
#define PURRS_poly_factor_hh 1

#include "Expr.types.hh"
#include "Symbol.types.hh"
#include <vector>

namespace Parma_Recurrence_Relation_Solver {

int
poly_factor(const Expr& p, const Symbol& x, std::vector<Expr>& factors);

} // namespace Parma_Recurrence_Relation_Solver

#endif
