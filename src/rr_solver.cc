
void
clear(GList& l) {
  for (unsigned n = l.nops(); n-- > 0; )
    l.remove_first();
  assert(l.nops() == 0);
}

static GExpr
get_binding(const GList& l, unsigned wild_index) {
  assert(wild_index < l.nops());
  assert(l.op(wild_index).info(GiNaC::info_flags::relation_equal));
  assert(l.op(wild_index).lhs() == GiNaC::wild(wild_index));
  return l.op(wild_index).rhs();
}

static bool
get_linear_decrement(const GExpr& e, const GSymbol& n, GNumber& decrement) {
  cout << "get_l_d(" << e << ", " << n << ", " << decrement << ")" << endl;
  static GExpr n_plus_d = n + GiNaC::wild(0);
  GList substitution;
  if (match(e, n_plus_d, substitution)) {
    GExpr d = get_binding(substitution, 0);
    if (GiNaC::is_a<GiNaC::numeric>(d)) {
      decrement = -GiNaC::ex_to<GiNaC::numeric>(d);
      return true;
    }
  }
  return false;
}

bool
solve(const GExpr& rhs, const GSymbol& n) {
  cout << "Trying to solve a(n) = " << rhs << endl;

  static GExpr x_i = x(GiNaC::wild(0));
  static GExpr x_i_plus_r = x_i + GiNaC::wild(1);
  static GExpr a_times_x_i = GiNaC::wild(1)*x_i;
  static GExpr a_times_x_i_plus_r = a_times_x_i + GiNaC::wild(2);

  static GList substitution;

  int degree = -1;
  vector<GNumber> coefficients;
  GExpr e = rhs;
  bool failed = false;
  do {
    cout << "e = " << e << endl;
    GExpr i;
    GExpr a;
    // The following matches are attempted starting from the most common,
    // then the second most common and so forth.
    if (clear(substitution), match(e, x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = get_binding(substitution, 1);
    }
    else if (clear(substitution), match(e, a_times_x_i_plus_r, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = get_binding(substitution, 2);
    }
    else if (clear(substitution), match(e, a_times_x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = get_binding(substitution, 1);
      e = 0;
    }
    else if (clear(substitution), match(e, x_i, substitution)) {
      i = get_binding(substitution, 0);
      a = 1;
      e = 0;
    }
    else
      break;

    GNumber decrement;
    if (!get_linear_decrement(i, n, decrement)) {
      failed = true;
      break;
    }
    cout << "decrement = " << decrement << endl;
    if (!decrement.is_integer()
	|| decrement < 0
	|| decrement >= coefficients.max_size()) {
      failed = true;
      break;
    }
    if (!GiNaC::is_a<GiNaC::numeric>(a)) {
      failed = true;
      break;
    }
    GNumber coefficient = GiNaC::ex_to<GiNaC::numeric>(a);
    unsigned long index = decrement.to_long();
    if (degree < 0 || index > unsigned(degree))
      degree = index;
    if (index > coefficients.size())
      coefficients.insert(coefficients.end(),
			  index - coefficients.size(),
			  GNumber(0));
    if (index == coefficients.size())
      coefficients.insert(coefficients.end(), coefficient);
    else
      coefficients[index] += coefficient;
  } while (e != 0);

  if (e.has(x_i))
    failed = true;

  if (failed)
    return false;
  else {
    cout << "Degree = " << degree << endl;
    cout << "Coefficients = ";
    for (int i = 1; i <= degree; ++i)
      cout << coefficients[i] << " ";
    cout << endl;
    cout << "Inhomogeneous term = " << e << endl;
    switch (degree) {
    case 1:
      cout << "1 root: " << coefficients[1] << endl;
      break;
    case 2:
      {
	GExpr x1;
	GExpr x2;
	bool d = solve_equation_2(-coefficients[1], -coefficients[2], x1, x2);
	cout << "2 " << (d ? "distinct" : "non-distinct") << " roots:" << endl
	     << x1 << endl << x2 << endl;
      }
      break;
    case 3:
      {
	bool all_real;
	GExpr x1;
	GExpr x2;
	GExpr x3;
	bool d = solve_equation_3(-coefficients[1], -coefficients[2],
				  -coefficients[3],
				  x1, x2, x3, all_real);
	cout << "3 " << (d ? "distinct" : "non-distinct") << " roots:" << endl
	     << x1 << endl << x2 << endl << x3 << endl
	     << x1.evalf() << endl
	     << x2.evalf() << endl
	     << x3.evalf() << endl
	     << x1.expand() << endl << x2 << endl << x3 << endl;
      }
      break;
    case 4:
      {
	GExpr x1;
	GExpr x2;
	GExpr x3;
	GExpr x4;
	bool d = solve_equation_4(-coefficients[1], -coefficients[2],
				  -coefficients[3], -coefficients[4],
				  x1, x2, x3, x4);
	cout << "4 " << (d ? "distinct" : "non-distinct") << " roots:" << endl
	     << x1 << endl << x2 << endl << x3 << endl << x4 << endl;
      }
      break;
    default:
      abort();
      break;
    }
    return true;
  }
}
