
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cassert>

#include "purrs_install.hh"

using namespace GiNaC;

void
solve_equation_4(const GNumber& a1, const GNumber& a2,
		 const GNumber& a3, const GNumber& a4,
		 GExpr& x1, GExpr& x2, GExpr& x3, GExpr& x4);

void foo() {
  GExpr x1, x2, x3, x4;
  solve_equation_4(0, -2, 0, 1, x1, x2, x3, x4);
  cout << "x1 = " << x1 << endl;
  cout << "****** " << x1.evalf() << endl;
  cout << "x2 = " << x2 << endl;
  cout << "****** " << x2.evalf() << endl;
  cout << "x3 = " << x3 << endl;
  cout << "****** " << x3.evalf() << endl;
  cout << "x4 = " << x4 << endl;
  cout << "****** " << x4.evalf() << endl;
}

int
main() try {
  GSymbol x("x");
  GList l(x);

  foo();

  GExpr p;
  while (cin) {
    string s;
    getline(cin, s);
    if (!cin) {
      cout << endl;
      break;
    }
    p = GExpr(s, l);
    std::vector<GExpr> roots;
    bool all_distinct;
    if (!find_roots(p, x, roots, all_distinct))
      cout << "Sorry, this is too difficult." << endl;
    else {
      size_t n = roots.size();
      for (size_t i = 0; i < n; ++i) {
	cout << "x_" << i << " = " << roots[i] << endl;
	if (!is_a<numeric>(roots[i]))
	  cout << "  x_" << i << " = " << roots[i].evalf() << endl;
      }
    }
  }
  return 0;
} catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
}
