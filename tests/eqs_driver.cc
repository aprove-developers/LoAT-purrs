#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cassert>

#include "purrs_install.hh"

using namespace GiNaC;

#define NOISY 1

int
main() try {
  GSymbol x("x");
  GList symbols(x);
  GExpr p;
  while (cin) {
    string s;
    getline(cin, s);

    if (!cin)
      return 0;

    // Skip comments.
    if (s.find("%") == 0)
      continue;

    if (s.find("symbols ") == 0) {
      for (int i = 7, l = s.length(); i < l; ++i) {
	char c = s[i];
	if (c >= 'a' && c <= 'z') {
	  char name[2];
	  name[0] = c;
	  name[1] = '\0';
	  GSymbol new_symbol(name);
	  symbols.append(new_symbol);
	}
      }
      continue;
    }

    p = GExpr(s, symbols);
    if (p == GExpr(0))
      continue;

#if NOISY
    cout << "Trying to solve " << p << " = 0" << endl;
#endif

    std::vector<GExpr> roots;
    bool all_distinct;
    if (!find_roots(p, x, roots, all_distinct)) {
#if NOISY
      cout << "Sorry, this is too difficult." << endl;
#endif
    }
#if NOISY
    else {
      size_t n = roots.size();
      for (size_t i = 0; i < n; ++i) {
	cout << "x_" << i+1 << " = " << roots[i] << endl;
	if (!is_a<numeric>(roots[i]))
	  cout << "****  x_" << i+1 << " ~= " << roots[i].evalf() << endl;
      }
    }
#endif
  }
  return 0;
} catch (exception &p) {
#if NOISY
  cerr << "Exception caught: " << p.what() << endl;
#endif
  return 1;
}
