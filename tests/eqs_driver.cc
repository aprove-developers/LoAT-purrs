
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cassert>

#include "purrs_install.hh"

using namespace GiNaC;

int
main() try {
  GSymbol x("x");
  GList l(x);
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
	  cout << "****  x_" << i+1 << " ~= " << roots[i].evalf() << endl;
      }
    }
  }
  return 0;
} catch (exception &p) {
  cerr << "Exception caught: " << p.what() << endl;
}
