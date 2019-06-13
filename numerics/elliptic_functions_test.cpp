
#include "numerics/elliptic_functions.hpp"

#include "gtest/gtest.h"

namespace principia {
namespace numerics {

class EllipticFunctionsTest : public ::testing::Test {
 protected:
  struct Arguments {
    double n;
    double m;
    double phi_over_PI;
  };
  struct Values {
    double elb;
    double eld;
    double elj;
  };
};

TEST_F(EllipticFunctionsTest, Xgscd) {
  // These expected values come from the output given by Fukushima at the end of
  // his xgscd.txt.  For our purpose, an unit test with assertions is more
  // useful than eyeballing the output.
  std::vector<std::pair<Arguments, Values>> const xgscd_expected_ = {};
  constexpr double PI = 3.1415926535897932384626433;
  constexpr double PIHALF = PI * 0.5;
  double dmc, mc, m, du, u, s, c, d;
  int jend, iend;
  jend = 10;
  iend = 8;
  dmc = 1.0 / static_cast<double>(jend);
  std::printf("%10s%10s%25s%25s%25s\n", "m", "u", "s", "c", "d");
  for (int j = 1; j <= jend; ++j) {
    mc = static_cast<double>(j) * dmc;
    m = 1.0 - mc;
    du = Elk(mc) / static_cast<double>(iend);
    for (int i = 0; i <= iend * 8; ++i) {
      u = du * static_cast<double>(i);
      Gscd(u, mc, s, c, d);
      std::printf("%10.5f%10.5f%25.15e%25.15e%25.15e\n", m, u, s, c, d);
    }
    std::printf("\n");
  }
}

}  // namespace numerics
}  // namespace principia
