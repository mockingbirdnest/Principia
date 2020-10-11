
#include "numerics/quadrature.hpp"

#include <limits>

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
using ::testing::Eq;

class QuadratureTest : public ::testing::Test {};

TEST_F(QuadratureTest, Sin) {
  int evaluations = 0;
  auto const f = [&evaluations](Angle const x) {
    ++evaluations;
    return Sin(x);
  };
  auto const ʃf = (Cos(2.0 * Radian) - Cos(5.0 * Radian)) * Radian;
  EXPECT_THAT(GaussLegendre<1>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(10_⑴)));
  EXPECT_THAT(GaussLegendre<2>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.3_⑴)));
  EXPECT_THAT(GaussLegendre<3>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.0e-1_⑴)));
  EXPECT_THAT(GaussLegendre<4>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.4e-2_⑴)));
  EXPECT_THAT(GaussLegendre<5>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(8.6e-4_⑴)));
  EXPECT_THAT(GaussLegendre<6>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.1e-5_⑴)));
  EXPECT_THAT(GaussLegendre<7>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.6e-7_⑴)));
  EXPECT_THAT(GaussLegendre<8>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.7e-9_⑴)));
  EXPECT_THAT(GaussLegendre<9>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.8e-11_⑴)));
  EXPECT_THAT(GaussLegendre<10>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 2495, 2498));
  EXPECT_THAT(GaussLegendre<11>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 20));
  EXPECT_THAT(GaussLegendre<12>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 6, 7));
  EXPECT_THAT(GaussLegendre<13>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 0, 1));
  EXPECT_THAT(GaussLegendre<14>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 1, 2));
  EXPECT_THAT(GaussLegendre<15>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 3, 4));

  evaluations = 0;
  EXPECT_THAT(AutomaticClenshawCurtis(nullptr, f,
                                      -2.0 * Radian,
                                      5.0 * Radian,
                                      std::numeric_limits<double>::epsilon()),
              AlmostEquals(ʃf, 3, 4));
  EXPECT_THAT(evaluations, Eq(134));
}

TEST_F(QuadratureTest, Sin2) {
  auto const f = [](Angle const x) { return Sin(2 * x); };
  auto const ʃf = (Cos(4 * Radian) - Cos(10 * Radian)) / 2 * Radian;
  EXPECT_THAT(GaussLegendre<1>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(9.7_⑴)));
  EXPECT_THAT(GaussLegendre<2>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(7.6_⑴)));
  EXPECT_THAT(GaussLegendre<3>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(7.6_⑴)));
  EXPECT_THAT(GaussLegendre<4>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.4_⑴)));
  EXPECT_THAT(GaussLegendre<5>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.2e-1_⑴)));
  EXPECT_THAT(GaussLegendre<6>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(4.6e-2_⑴)));
  EXPECT_THAT(GaussLegendre<7>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.5e-3_⑴)));
  EXPECT_THAT(GaussLegendre<8>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.0e-4_⑴)));
  EXPECT_THAT(GaussLegendre<9>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(8.5e-6_⑴)));
  EXPECT_THAT(GaussLegendre<10>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(2.9e-7_⑴)));;
  EXPECT_THAT(GaussLegendre<11>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(8.1e-9_⑴)));
  EXPECT_THAT(GaussLegendre<12>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(1.9e-10_⑴)));
  EXPECT_THAT(GaussLegendre<13>(f, -2.0 * Radian, 5.0 * Radian),
              RelativeErrorFrom(ʃf, IsNear(3.7e-12_⑴)));
  EXPECT_THAT(GaussLegendre<14>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 422, 425));
  EXPECT_THAT(GaussLegendre<15>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 36, 50));
  EXPECT_THAT(GaussLegendre<16>(f, -2.0 * Radian, 5.0 * Radian),
              AlmostEquals(ʃf, 152, 159));
}

TEST_F(QuadratureTest, Sin10) {
  int evaluations = 0;
  auto const f = [&evaluations](Angle const x) {
    ++evaluations;
    return Sin(10 * x);
  };
  auto const ʃf = (Cos(20 * Radian) - Cos(50 * Radian)) / 10 * Radian;
  EXPECT_THAT(AutomaticClenshawCurtis(nullptr, f,
                                      -2.0 * Radian,
                                      5.0 * Radian,
                                      std::numeric_limits<double>::epsilon()),
              AlmostEquals(ʃf, 144, 277));
  EXPECT_THAT(evaluations, Eq(263));
}

}  // namespace quadrature
}  // namespace numerics
}  // namespace principia
