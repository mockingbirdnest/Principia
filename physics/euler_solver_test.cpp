
#include "physics/euler_solver.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <set>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using geometry::Bivector;
using geometry::Frame;
using geometry::Instant;
using geometry::R3Element;
using quantities::Abs;
using quantities::AngularFrequency;
using quantities::AngularMomentum;
using quantities::Cos;
using quantities::Energy;
using quantities::MomentOfInertia;
using quantities::Sin;
using quantities::SIUnit;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::VanishesBefore;
using ::testing::Lt;
using ::testing::Matcher;

class EulerSolverTest : public ::testing::Test {
 protected:
  using ICRSFrame = Frame<serialization::Frame::SolarSystemTag,
                          serialization::Frame::ICRS,
                          /*frame_is_inertial*/ true>;
  using PrincipalAxesFrame = Frame<serialization::Frame::PhysicsTag,
                                   serialization::Frame::PRINCIPAL_AXES,
                                   /*frame_is_inertial*/ false>;

  using Solver = EulerSolver<ICRSFrame, PrincipalAxesFrame>;

  EulerSolverTest()
      : identity_attitude_(
            EulerSolver<ICRSFrame,
                        PrincipalAxesFrame>::AttitudeRotation::Identity()) {}

  // Checks that the angular momentum transformed by the attitude to the
  // inertial frame satisfies the given matcher.
  void CheckAngularMomentumConservation(
      std::vector<Solver::AngularMomentumBivector> const& angular_momenta,
      std::vector<Solver::AttitudeRotation> const& attitudes,
      Matcher<Bivector<AngularMomentum, ICRSFrame>> const& matcher) {
    CHECK_EQ(angular_momenta.size(), attitudes.size());
    for (int i = 0; i < angular_momenta.size(); ++i) {
      Bivector<AngularMomentum, ICRSFrame> const angular_momentum_in_inertial =
          attitudes[i](angular_momenta[i]);
      EXPECT_THAT(angular_momentum_in_inertial, matcher);
    }
  }

  // Checks that the kinetic energy computed using the attitude has the right
  // value.
  void CheckPoinsotConstruction(
      R3Element<MomentOfInertia> const& moments_of_inertia,
      std::vector<Solver::AngularMomentumBivector> const& angular_momenta,
      std::vector<Solver::AttitudeRotation> const& attitudes,
      int const ulps) {
    CHECK_EQ(angular_momenta.size(), attitudes.size());
    Energy maximum_kinetic_energy;
    Energy minimum_kinetic_energy = quantities::Infinity<Energy>();
    for (int i = 0; i < angular_momenta.size(); ++i) {
      Bivector<AngularMomentum, PrincipalAxesFrame> const angular_momentum =
          angular_momenta[i];
      Bivector<AngularFrequency, PrincipalAxesFrame> const angular_velocity(
          {angular_momentum.coordinates().x / moments_of_inertia.x,
           angular_momentum.coordinates().y / moments_of_inertia.y,
           angular_momentum.coordinates().z / moments_of_inertia.z});
      Bivector<AngularMomentum, ICRSFrame> const angular_momentum_in_inertial =
          attitudes[i](angular_momentum);
      Bivector<AngularFrequency, ICRSFrame> const
          angular_velocity_in_inertial = attitudes[i](angular_velocity);
      Energy const kinetic_energy = 0.5 *
                                    InnerProduct(angular_momentum_in_inertial,
                                                 angular_velocity_in_inertial) /
                                    Radian / Radian;
      maximum_kinetic_energy = std::max(maximum_kinetic_energy, kinetic_energy);
      minimum_kinetic_energy = std::min(minimum_kinetic_energy, kinetic_energy);
    }
    EXPECT_THAT(minimum_kinetic_energy,
                AlmostEquals(maximum_kinetic_energy, ulps));
  }

  // TODO(phl): Test with different initial attitudes.
  Solver::AttitudeRotation identity_attitude_;
};

// Check that we are able to retrieve the initial state for random choices of
// the moments of inertia and the angular momentum.
TEST_F(EulerSolverTest, InitialStateRandom) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> moment_of_inertia_distribution(0.0, 10.0);
  std::uniform_real_distribution<> angular_momentum_distribution(-10.0, 10.0);
  for (int i = 0; i < 1000; ++i) {
    // Make sure that the moments of inertia are properly ordered.
    std::array<double, 3> randoms{moment_of_inertia_distribution(random),
                                  moment_of_inertia_distribution(random),
                                  moment_of_inertia_distribution(random)};
    std::sort(randoms.begin(), randoms.end());
    R3Element<MomentOfInertia> const moments_of_inertia{
        randoms[0] * SIUnit<MomentOfInertia>(),
        randoms[1] * SIUnit<MomentOfInertia>(),
        randoms[2] * SIUnit<MomentOfInertia>()};

    Solver::AngularMomentumBivector
        initial_angular_momentum(
            {angular_momentum_distribution(random) * SIUnit<AngularMomentum>(),
             angular_momentum_distribution(random) * SIUnit<AngularMomentum>(),
             angular_momentum_distribution(random) *
                 SIUnit<AngularMomentum>()});

    Solver const solver(moments_of_inertia,
                        initial_angular_momentum,
                        identity_attitude_,
                        Instant());
    auto const computed_initial_angular_momentum =
        solver.AngularMomentumAt(Instant());

    EXPECT_THAT(computed_initial_angular_momentum,
                AlmostEquals(initial_angular_momentum, 0, 209))
        << moments_of_inertia << " " << initial_angular_momentum;
  }
}

// Same as above, but exercises the symmetrical cases where at least two moments
// of inertia are equal.
TEST_F(EulerSolverTest, InitialStateSymmetrical) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> angular_momentum_distribution(-10.0, 10.0);

  R3Element<MomentOfInertia> const moments_of_inertia1{
      2 * SIUnit<MomentOfInertia>(),
      2 * SIUnit<MomentOfInertia>(),
      3 * SIUnit<MomentOfInertia>()};
  R3Element<MomentOfInertia> const moments_of_inertia2{
      2 * SIUnit<MomentOfInertia>(),
      3 * SIUnit<MomentOfInertia>(),
      3 * SIUnit<MomentOfInertia>()};
  R3Element<MomentOfInertia> const moments_of_inertia3{
      3 * SIUnit<MomentOfInertia>(),
      3 * SIUnit<MomentOfInertia>(),
      3 * SIUnit<MomentOfInertia>()};

  for (int i = 0; i < 100; ++i) {
    Solver::AngularMomentumBivector
        initial_angular_momentum(
            {angular_momentum_distribution(random) * SIUnit<AngularMomentum>(),
             angular_momentum_distribution(random) * SIUnit<AngularMomentum>(),
             angular_momentum_distribution(random) *
                 SIUnit<AngularMomentum>()});
    {
      Solver const solver1(moments_of_inertia1,
                           initial_angular_momentum,
                           identity_attitude_,
                           Instant());
      auto const computed_initial_angular_momentum1 =
          solver1.AngularMomentumAt(Instant());

      EXPECT_THAT(computed_initial_angular_momentum1,
                  AlmostEquals(initial_angular_momentum, 0, 87))
          << moments_of_inertia1 << " " << initial_angular_momentum;
    }
    {
      Solver const solver2(moments_of_inertia2,
                           initial_angular_momentum,
                           identity_attitude_,
                           Instant());
      auto const computed_initial_angular_momentum2 =
          solver2.AngularMomentumAt(Instant());

      EXPECT_THAT(computed_initial_angular_momentum2,
                  AlmostEquals(initial_angular_momentum, 0, 50))
          << moments_of_inertia2 << " " << initial_angular_momentum;
    }
    {
      Solver const solver3(moments_of_inertia3,
                           initial_angular_momentum,
                           identity_attitude_,
                           Instant());
      auto const computed_initial_angular_momentum3 =
          solver3.AngularMomentumAt(Instant());

      EXPECT_THAT(computed_initial_angular_momentum3,
                  AlmostEquals(initial_angular_momentum, 0, 0))
          << moments_of_inertia3 << " " << initial_angular_momentum;
      // TODO(phl): Test the attitude in this case.
    }
  }
}

// Same as above, but exercises all the formulæ.  We compute an angular
// momentum by fixing its first coordinate and picking the third coordinate so
// that it falls in the right interval.  (The second coordinate turns out to be
// irrelevant.)
TEST_F(EulerSolverTest, InitialStateFormulæ) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> moment_of_inertia_distribution(0.0, 10.0);
  std::uniform_real_distribution<> angular_momentum_distribution(-10.0, 10.0);
  for (int i = 0; i < 1000; ++i) {
    // Make sure that the moments of inertia are properly ordered.
    std::array<double, 3> randoms{moment_of_inertia_distribution(random),
                                  moment_of_inertia_distribution(random),
                                  moment_of_inertia_distribution(random)};
    std::sort(randoms.begin(), randoms.end());
    auto const I₁ = randoms[0] * SIUnit<MomentOfInertia>();
    auto const I₂ = randoms[1] * SIUnit<MomentOfInertia>();
    auto const I₃ = randoms[2] * SIUnit<MomentOfInertia>();
    R3Element<MomentOfInertia> const moments_of_inertia{I₁, I₂, I₃};

    // G² = T * (I₁ + I₂)
    {
      auto const mx =
          angular_momentum_distribution(random) * SIUnit<AngularMomentum>();
      auto mz = mx * Sqrt(((I₂ - I₁) * I₃) / ((2.0 * I₃ - I₂ - I₁) * I₁));
      if (i % 2 == 0) {
        mz = -mz;
      }
      Solver::AngularMomentumBivector
          initial_angular_momentum({mx, SIUnit<AngularMomentum>(), mz});
      Solver const solver(moments_of_inertia,
                          initial_angular_momentum,
                          identity_attitude_,
                          Instant());

      auto const computed_initial_angular_momentum =
          solver.AngularMomentumAt(Instant());
      EXPECT_THAT(computed_initial_angular_momentum,
                  AlmostEquals(initial_angular_momentum, 0, 356))
          << moments_of_inertia << " " << initial_angular_momentum;
    }

    // G² = 2 * T * I₂
    {
      auto const mx =
          angular_momentum_distribution(random) * SIUnit<AngularMomentum>();
      auto mz = mx * Sqrt(((I₂ - I₁) * I₃) / ((I₃ - I₂) * I₁));
      if (i % 2 == 0) {
        mz = -mz;
      }
      Solver::AngularMomentumBivector
          initial_angular_momentum({mx, SIUnit<AngularMomentum>(), mz});
      Solver const solver(moments_of_inertia,
                          initial_angular_momentum,
                          identity_attitude_,
                          Instant());

      auto const computed_initial_angular_momentum =
          solver.AngularMomentumAt(Instant());
      // NOTE(phl): The largest error happens to actually go through
      // Formula::ii and is on the z component (x and y are fine).  That's
      // probably related to the fact that Δ₂ is very small.
      EXPECT_THAT(computed_initial_angular_momentum,
                  AlmostEquals(initial_angular_momentum, 0, 11126))
          << moments_of_inertia << " " << initial_angular_momentum;
    }

    // G² = T * (I₂ + I₃)
    {
      auto const mx =
          angular_momentum_distribution(random) * SIUnit<AngularMomentum>();
      auto mz = mx * Sqrt(((I₂ + I₃ - 2.0 * I₁) * I₃) / ((I₃ - I₂) * I₁));
      if (i % 2 == 0) {
        mz = -mz;
      }
      Solver::AngularMomentumBivector
          initial_angular_momentum({mx, SIUnit<AngularMomentum>(), mz});
      Solver const solver(moments_of_inertia,
                          initial_angular_momentum,
                          identity_attitude_,
                          Instant());

      auto const computed_initial_angular_momentum =
          solver.AngularMomentumAt(Instant());
      EXPECT_THAT(computed_initial_angular_momentum,
                  AlmostEquals(initial_angular_momentum, 0, 2711))
          << moments_of_inertia << " " << initial_angular_momentum;
    }
  }
}

// This test and the next come from
// http://n.ethz.ch/~stiegerc/HS09/Mechanik/Unterlagen/Lecture19.pdf.
TEST_F(EulerSolverTest, ShortFatSymmetricTopPrecession) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      3.0 * SIUnit<MomentOfInertia>(),
      3.0 * SIUnit<MomentOfInertia>(),
      9.0 * SIUnit<MomentOfInertia>()};

  Solver::AngularMomentumBivector const
      initial_angular_momentum({0.0 * SIUnit<AngularMomentum>(),
                                5.0 * SIUnit<AngularMomentum>(),
                                7.0 * SIUnit<AngularMomentum>()});
  Solver::AttitudeRotation const initial_attitude = identity_attitude_;

  // Correspondence with the referential of lecture 19: x = e1, y = e2, z = e3.
  AngularFrequency Ω = initial_angular_momentum.coordinates().z *
                       (moments_of_inertia[0] - moments_of_inertia[2]) /
                       (moments_of_inertia[0] * moments_of_inertia[2]);

  Solver const solver(moments_of_inertia,
                      initial_angular_momentum,
                      initial_attitude,
                      Instant());

  std::vector<Solver::AngularMomentumBivector> angular_momenta;
  std::vector<Solver::AttitudeRotation> attitudes;
  for (Time t = 0 * Second; t < 5.0 * Second; t += 0.1 * Second) {
    auto const angular_momentum_at_t = solver.AngularMomentumAt(Instant() + t);
    EXPECT_THAT(
        angular_momentum_at_t,
        AlmostEquals(
            Solver::AngularMomentumBivector(
                {5.0 * Sin(Ω * t) * SIUnit<AngularMomentum>(),
                 5.0 * Cos(Ω * t) * SIUnit<AngularMomentum>(),
                 7.0 * SIUnit<AngularMomentum>()}),
            0,
            102))
        << t;
    angular_momenta.push_back(angular_momentum_at_t);
    attitudes.push_back(solver.AttitudeAt(angular_momentum_at_t,
                                          Instant() + t));
  }

  Bivector<AngularMomentum, ICRSFrame> const reference_momentum =
      initial_attitude(initial_angular_momentum);
  CheckAngularMomentumConservation(
      angular_momenta,
      attitudes,
      Componentwise(VanishesBefore(1 * SIUnit<AngularMomentum>(), 0, 32),
                    AlmostEquals(reference_momentum.coordinates().y, 0, 4),
                    AlmostEquals(reference_momentum.coordinates().z, 0, 2)));
  CheckPoinsotConstruction(
      moments_of_inertia, angular_momenta, attitudes, /*ulps=*/13);
}

TEST_F(EulerSolverTest, TallSkinnySymmetricTopPrecession) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      3.0 * SIUnit<MomentOfInertia>(),
      9.0 * SIUnit<MomentOfInertia>(),
      9.0 * SIUnit<MomentOfInertia>()};

  Solver::AngularMomentumBivector const
      initial_angular_momentum({7.0 * SIUnit<AngularMomentum>(),
                                0.0 * SIUnit<AngularMomentum>(),
                                5.0 * SIUnit<AngularMomentum>()});
  Solver::AttitudeRotation const initial_attitude = identity_attitude_;

  // Correspondence with the referential of lecture 19:  x = e3, y = e1, z = e2.
  AngularFrequency Ω = initial_angular_momentum.coordinates().x *
                       (moments_of_inertia[1] - moments_of_inertia[0]) /
                       (moments_of_inertia[1] * moments_of_inertia[0]);

  Solver const solver(moments_of_inertia,
                      initial_angular_momentum,
                      initial_attitude,
                      Instant());

  std::vector<Solver::AngularMomentumBivector> angular_momenta;
  std::vector<Solver::AttitudeRotation> attitudes;
  for (Time t = 0 * Second; t < 5.0 * Second; t += 0.1 * Second) {
    auto const angular_momentum_at_t = solver.AngularMomentumAt(Instant() + t);
    EXPECT_THAT(
        angular_momentum_at_t,
        AlmostEquals(
            Solver::AngularMomentumBivector(
                {7.0 * SIUnit<AngularMomentum>(),
                 5.0 * Sin(Ω * t) * SIUnit<AngularMomentum>(),
                 5.0 * Cos(Ω * t) * SIUnit<AngularMomentum>()}),
            0,
            34))
        << t;
    angular_momenta.push_back(angular_momentum_at_t);
    attitudes.push_back(solver.AttitudeAt(angular_momentum_at_t,
                                          Instant() + t));
  }

  Bivector<AngularMomentum, ICRSFrame> const reference_momentum =
      initial_attitude(initial_angular_momentum);
  CheckAngularMomentumConservation(
      angular_momenta,
      attitudes,
      Componentwise(AlmostEquals(reference_momentum.coordinates().x, 0, 3),
                    VanishesBefore(1 * SIUnit<AngularMomentum>(), 0, 24),
                    AlmostEquals(reference_momentum.coordinates().z, 0, 6)));
  CheckPoinsotConstruction(
      moments_of_inertia, angular_momenta, attitudes, /*ulps=*/4);
}

// This test demonstrates the Джанибеков effect, also known as tennis racket
// theorem: the rotation of an object around its second principal axis is not
// stable.  Here we choose the initial angular momentum to be mostly in the y
// direction with a small component in the z direction.  This causes the object
// to periodically flip, rotating along y or along -y.
TEST_F(EulerSolverTest, ДжанибековEffect) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      3.0 * SIUnit<MomentOfInertia>(),
      5.0 * SIUnit<MomentOfInertia>(),
      9.0 * SIUnit<MomentOfInertia>()};

  Solver::AngularMomentumBivector const
      initial_angular_momentum({0.0 * SIUnit<AngularMomentum>(),
                                2.0 * SIUnit<AngularMomentum>(),
                                0.01 * SIUnit<AngularMomentum>()});
  Solver::AttitudeRotation const initial_attitude = identity_attitude_;

  Solver const solver(moments_of_inertia,
                      initial_angular_momentum,
                      initial_attitude,
                      Instant());

  // Find the maxima, minima and zeroes of the y coordinate of the angular
  // momentum.
  AngularMomentum previous_my = initial_angular_momentum.coordinates().y;
  Instant previous_t;
  std::vector<Instant> maxima;
  std::vector<Instant> minima;
  std::vector<Instant> zeroes;
  bool is_abs_decreasing = false;
  bool is_decreasing = false;
  bool is_increasing = true;
  std::vector<Solver::AngularMomentumBivector> angular_momenta;
  std::vector<Solver::AttitudeRotation> attitudes;
  for (Instant t; t < Instant() + 100.0 * Second; t += 0.1 * Second) {
    auto const angular_momentum_at_t = solver.AngularMomentumAt(t);
    auto const my = angular_momentum_at_t.coordinates().y;
    if (is_increasing && my > 1.99 * SIUnit<AngularMomentum>()) {
      if (my < previous_my) {
        maxima.push_back(previous_t);
        is_abs_decreasing = true;
        is_decreasing = true;
        is_increasing = false;
      }
    }
    if (is_decreasing && my < -1.99 * SIUnit<AngularMomentum>()) {
      if (my > previous_my) {
        minima.push_back(previous_t);
        is_abs_decreasing = true;
        is_decreasing = false;
        is_increasing = true;
      }
    }
    if (is_abs_decreasing && Abs(my) < 0.1 * SIUnit<AngularMomentum>()) {
      if (Abs(my) > Abs(previous_my)) {
        zeroes.push_back(previous_t);
        is_abs_decreasing = false;
      }
    }
    previous_my = my;
    previous_t = t;
    angular_momenta.push_back(angular_momentum_at_t);
    attitudes.push_back(solver.AttitudeAt(angular_momentum_at_t,
                                          t));
  }

  // Check that the maxima, minima and zeroes properly alternate and are
  // roughly equidistant.
  std::set<Instant> all;
  all.insert(maxima.begin(), maxima.end());
  all.insert(minima.begin(), minima.end());
  all.insert(zeroes.begin(), zeroes.end());
  EXPECT_EQ(maxima.size() + minima.size() + zeroes.size(), all.size());
  Time const quarter_period = (*all.rbegin() - *all.begin()) / (all.size() - 1);
  for (auto it = all.begin(); it != all.end(); ++it) {
    auto const t = *it;
    int const i = std::distance(all.begin(), it);
    if (i % 4 == 0) {
      EXPECT_EQ(maxima[i / 4], t);
    }
    if (i % 4 == 2) {
      EXPECT_EQ(minima[i / 4], t);
    }
    if (i % 4 == 1 || i % 4 == 3) {
      EXPECT_EQ(zeroes[i / 2], t);
    }
    if (it != all.begin()) {
      EXPECT_THAT(RelativeError(quarter_period, t - *std::prev(it)),
                  Lt(0.0023));
    }
  }

  Bivector<AngularMomentum, ICRSFrame> const reference_momentum =
      initial_attitude(initial_angular_momentum);
  CheckAngularMomentumConservation(
      angular_momenta,
      attitudes,
      Componentwise(VanishesBefore(1 * SIUnit<AngularMomentum>(), 0, 10),
                    AlmostEquals(reference_momentum.coordinates().y, 0, 13),
                    AlmostEquals(reference_momentum.coordinates().z, 0, 965)));
  CheckPoinsotConstruction(
      moments_of_inertia, angular_momenta, attitudes, /*ulps=*/33);
}

}  // namespace physics
}  // namespace principia
