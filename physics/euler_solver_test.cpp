
#include "physics/euler_solver.hpp"

#include <algorithm>
#include <array>
#include <random>

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

// Check that we are able to retrieve the initial state for random choices of
// the moments of inertia and the angular momentum.
TEST_F(EulerSolverTest, InitialState) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> moment_of_inertia_distribution(0.0, 10.0);
  std::uniform_real_distribution<> angular_momentum_distribution(-10.0, 10.0);
  for (int i = 0; i < 1000; ++i) {
    // Make sure that the moments of inertia are properly ordered.
    std::array<double, 3> randoms{moment_of_inertia_distribution(random),
                                  moment_of_inertia_distribution(random),
                                  moment_of_inertia_distribution(random)};
    std::sort(randoms.begin(), randoms.end());
    R3Element<MomentOfInertia> const& moments_of_inertia{
        randoms[0] * SIUnit<MomentOfInertia>(),
        randoms[1] * SIUnit<MomentOfInertia>(),
        randoms[2] * SIUnit<MomentOfInertia>()};

    EulerSolver::AngularMomentumBivector initial_angular_momentum(
        {angular_momentum_distribution(random) * SIUnit<AngularMomentum>(),
         angular_momentum_distribution(random) * SIUnit<AngularMomentum>(),
         angular_momentum_distribution(random) * SIUnit<AngularMomentum>()});

    EulerSolver const solver(
        moments_of_inertia, initial_angular_momentum, Instant());
    auto const computed_initial_angular_momentum =
        solver.AngularMomentumAt(Instant());

    EXPECT_THAT(computed_initial_angular_momentum,
                AlmostEquals(initial_angular_momentum, 0, 20))
        << moments_of_inertia << " " << initial_angular_momentum;
  }
}

}  // namespace physics
}  // namespace principia
