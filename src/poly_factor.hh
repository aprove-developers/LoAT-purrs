#ifndef PURRS_poly_factor_hh
#define PURRS_poly_factor_hh 1

#include "globals.hh"
#include <vector>

int
poly_factor(const GExpr& p, const GSymbol& x, std::vector<GExpr>& factors);

#endif
