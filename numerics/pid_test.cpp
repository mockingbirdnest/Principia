
#include "numerics/pid.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using geometry::Frame;
using geometry::Velocity;
using quantities::Inverse;
using quantities::Time;
using quantities::si::Second;

class PIDTest : public ::testing::Test {
 protected:
  using F = Frame<enum class FTag>;

  static constexpr double kp = 0.5;
  static constexpr Inverse<Time> ki = 0.3 / Second;
  static constexpr Time kd = 0.2 * Second;
};

TEST_F(PIDTest, Double) {
  PID<double, /*past_horizon=*/10, /*finite_difference_order=*/3> pid(
      kp, ki, kd);
}

TEST_F(PIDTest, Geometry) {
  PID<Velocity<F>, /*past_horizon=*/10, /*finite_difference_order=*/3> pid(
    kp, ki, kd);
}

}  // namespace numerics
}  // namespace principia