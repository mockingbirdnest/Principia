
#include "numerics/quadrature.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {
namespace quadrature {

using quantities::Angle;
using quantities::Cos;
using quantities::Sin;
using quantities::si::Radian;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_⑴;

class QuadratureTest : public ::testing::Test {};

TEST_F(QuadratureTest, Sin) {
  EXPECT_THAT(
      GaussLegendre<5>(
          [](Angle const x) { return Sin(x); }, -2.0 * Radian, 5.0 * Radian),
      RelativeErrorFrom((Cos(2.0 * Radian) - Cos(5.0 * Radian)) * Radian,
                        IsNear(8.6e-4_⑴)));
  EXPECT_THAT(
      GaussLegendre<10>(
          [](Angle const x) { return Sin(x); }, -2.0 * Radian, 5.0 * Radian),
      AlmostEquals((Cos(2.0 * Radian) - Cos(5.0 * Radian)) * Radian, 2498));
  EXPECT_THAT(
      GaussLegendre<20>(
          [](Angle const x) { return Sin(x); }, -2.0 * Radian, 5.0 * Radian),
      AlmostEquals((Cos(2.0 * Radian) - Cos(5.0 * Radian)) * Radian, 3));
}

}  // namespace quadrature
}  // namespace numerics
}  // namespace principia