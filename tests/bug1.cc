#include <iostream>
#include <stdexcept>

#include "purrs_install.hh"

using namespace GiNaC;
using namespace std;

void
solve_equation_4(const GNumber& a1, const GNumber& a2,
		 const GNumber& a3, const GNumber& a4,
		 GExpr& x1, GExpr& x2, GExpr& x3, GExpr& x4);

int
main() try {
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
  return 0;
} catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
}
