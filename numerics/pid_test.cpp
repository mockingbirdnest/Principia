
#include "numerics/pid.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
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
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using quantities::Current;
using quantities::Inverse;
using quantities::Time;
using quantities::si::Ampere;
using quantities::si::Coulomb;
using quantities::si::Metre;
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
  PID<double, double, /*horizon=*/10, /*finite_difference_order=*/3> pid(
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
    pid.ComputeControlVariable(apparent(i), actual(i), Instant() + i * Second);
  }

  // Compute the various elements for i = 20.
  double const proportional = actual(20) - apparent(20);
  Time const integral = ((actual_primitive(20) - actual_primitive(11)) -
                         (apparent_primitive(20) - apparent_primitive(11))) *
                        Second;
  Inverse<Time> const derivative =
      (actual_derivative(20) - apparent_derivative(20)) / Second;

  // TODO(phl): The large relative error is because our integration is lame.
  // Fix it, one day.
  EXPECT_THAT(
      pid.ComputeControlVariable(
          apparent(20), actual(20), Instant() + 20 * Second),
      RelativeErrorFrom(kp * proportional + ki * integral + kd * derivative,
                        IsNear(0.09_⑴)));
}

TEST_F(PIDTest, Geometry) {
  using FunkyPID = PID<Position<F>,
                       Vector<Current, F>,
                       /*horizon=*/10,
                       /*finite_difference_order=*/3>;
  FunkyPID::Kp const kp = 1 * Ampere / Metre;
  FunkyPID::Ki const ki = 1 * Ampere / Metre / Second;
  FunkyPID::Kd const kd = 1 * Coulomb / Metre;

  // Merely a compilation test.
  FunkyPID pid(kp, ki, kd);
  auto i = pid.ComputeControlVariable(F::origin, F::origin, Instant());
}

}  // namespace numerics
}  // namespace principia
