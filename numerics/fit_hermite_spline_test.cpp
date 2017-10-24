
#include "numerics/fit_hermite_spline.hpp"

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {

using geometry::Instant;
using quantities::Length;
using quantities::Speed;
using quantities::si::Metre;

namespace numerics {

class FitHermiteSplineTest : public ::testing::Test {
 protected:
  struct Sample {
    Instant t;
    Length x;
    Speed v;
  };
};

TEST_F(FitHermiteSplineTest, JustMakeItCompile) {
  std::vector<Sample> samples;
  FitHermiteSpline(samples,
                   [](auto&& sample) -> auto&& { return sample.t; },
                   [](auto&& sample) -> auto&& { return sample.x; },
                   [](auto&& sample) -> auto&& { return sample.v; },
                   1 * Metre);
}

}  // namespace numerics
}  // namespace principia