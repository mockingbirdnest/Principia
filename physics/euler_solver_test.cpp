
#include "physics/euler_solver.hpp"

#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace physics {

using geometry::Instant;
using geometry::R3Element;
using quantities::AngularMomentum;
using quantities::MomentOfInertia;
using quantities::SIUnit;
using testing_utilities::AlmostEquals;

class EulerSolverTest : public ::testing::Test {
};

TEST_F(EulerSolverTest, InitialState) {
  R3Element<MomentOfInertia> const& moments_of_inertia{
      2 * SIUnit<MomentOfInertia>(),
      3 * SIUnit<MomentOfInertia>(),
      5 * SIUnit<MomentOfInertia>()};
  EulerSolver::AngularMomentumBivector initial_angular_momentum(
      {6 * SIUnit<AngularMomentum>(),
       7 * SIUnit<AngularMomentum>(),
       8 * SIUnit<AngularMomentum>()});
  EulerSolver const solver(
      moments_of_inertia, initial_angular_momentum, Instant());

  auto const computed_initial_angular_momentum =
      solver.AngularMomentumAt(Instant());

  EXPECT_THAT(computed_initial_angular_momentum,
              AlmostEquals(initial_angular_momentum, 0));
}

}  // namespace physics
}  // namespace principia
