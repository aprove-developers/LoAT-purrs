
#include <config.h>
#include "complint.hh"

int
lexicographic_compare(const double& x, const double& y) {
  if (x < y)
    return -1;
  else if (x > y)
    return 1;
  else
    return 0;
}

int
lexicographic_compare(const Interval& x, const Interval& y) {
  if (x.isEmpty())
    return y.isEmpty() ? 0 : -1;
  else if (y.isEmpty())
    return 1;
  else if (int r = lexicographic_compare(x.inf(), y.inf()))
    return r;
  else if (int r = lexicographic_compare(x.sup(), y.sup()))
    return r;
  else
    return 0;
}

int
lexicographic_compare(const CInterval& x, const CInterval& y) {
  if (int r = lexicographic_compare(x.re(), y.re()))
    return r;
  else if (int r = lexicographic_compare(x.im(), y.im()))
    return r;
  else
    return 0;
}

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS(complint, basic)

complint::complint()
  : inherited(TINFO_complint) {
}

complint::complint(const CInterval& i)
  : inherited(TINFO_complint), ci(i) {
}

void
complint::destroy(bool call_parent) {
  if (call_parent)
    inherited::destroy(call_parent);
}

void
complint::copy(const complint& other) {
  inherited::copy(other);
  ci = other.ci;
}

int
complint::compare_same_type(const basic& other) const {
  const complint& o = static_cast<const complint&>(other);
  return lexicographic_compare(ci, o.ci);
}

void
complint::print(const print_context& c, unsigned int /* level */) const {
  c.s << "complint(" << ci << ')';
}

void
complint::archive(archive_node&) const {
  abort();
}

ex
complint::unarchive(const archive_node&, lst&) {
  abort();
}

} // namespace GiNaC
