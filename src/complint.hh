
#include <ginac/ginac.h>
#include <Interval.h>
#include <cinterval.h>

typedef filib::Interval Interval;
typedef cinterval CInterval;

namespace GiNaC {

const unsigned TINFO_complint = 0x42420001U;

class complint : public basic {
  GINAC_DECLARE_REGISTERED_CLASS(complint, basic)

public:
  complint(const CInterval& i);

  void complint::print(const print_context& c, unsigned level) const;

private:
  CInterval ci;
};

} // namespace GiNaC
