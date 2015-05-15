
#include "numerics/чебышёв_series.hpp"

#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {

using geometry::Instant;
using si::Second;

namespace numerics {

class ЧебышёвSeriesTest : public ::testing::Test {
 protected:
  ЧебышёвSeriesTest()
      : t_min_(-1 * Second),
        t_max_(3 * Second) {}

  Instant t_min_;
  Instant t_max_;
};

TEST_F(ЧебышёвSeriesTest, T2) {
  ЧебышёвSeries<double> t2({0, 0, 1}, t_min_, t_max_);
  EXPECT_EQ(-1, t2.Evaluate(Instant(1 * Second)));
  EXPECT_EQ(1, t2.Evaluate(Instant(3 * Second)));
}

}  // namespace numerics
}  // namespace principia
