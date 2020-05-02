
#include "numerics/pid.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {

using geometry::Frame;
using geometry::Velocity;
using quantities::Inverse;
using quantities::Time;
using quantities::si::Second;
using testing_utilities::IsNear;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::operator""_⑴;

class PIDTest : public ::testing::Test {
 protected:
  using F = Frame<enum class FTag>;

  static constexpr double kp = 0.5;
  static constexpr Inverse<Time> ki = 0.3 / Second;
  static constexpr Time kd = 0.2 * Second;
};

TEST_F(PIDTest, Double) {
  PID<double, /*horizon=*/10, /*finite_difference_order=*/3> pid(
      kp, ki, kd);
  auto actual = [](int const step) { return 1.1 * step * step; };
  auto apparent = [](int const step) { return step; };
  auto actual_primitive = [](int const step) {
    return 1.1 * step * step * step / 3.0;
  };
  auto apparent_primitive = [](int const step) { return step * step / 2.0; };
  auto actual_derivative = [](int const step) { return 2.2 * step; };
  auto apparent_derivative = [](int const step) { return 1.0; };

  // Prime the PID.
  for (int i = 0; i < 20; ++i) {
    pid.ComputeValue(apparent(i), actual(i), 1 * Second);
  }

  // Compute the various elements for i = 20.
  double const proportional = apparent(20) - actual(20);
  Time const integral = ((apparent_primitive(20) - apparent_primitive(11)) -
                         (actual_primitive(20) - apparent_primitive(11))) *
                        Second;
  Inverse<Time> const derivative =
      (apparent_derivative(20) - actual_derivative(20)) / Second;

  // TODO(phl): The large relative error is because our integration is lame.
  // Fix it, one day.
  EXPECT_THAT(pid.ComputeValue(apparent(20), actual(20), 1 * Second),
              RelativeErrorFrom(actual(20) + kp * proportional + ki * integral +
                                    kd * derivative,
                                IsNear(0.08_⑴)));
}

TEST_F(PIDTest, Geometry) {
  PID<Velocity<F>, /*horizon=*/10, /*finite_difference_order=*/3> pid(
    kp, ki, kd);
  // Merely a compilation test.
  pid.ComputeValue(Velocity<F>{}, Velocity<F>{}, 1 * Second);
}

}  // namespace numerics
}  // namespace principia
