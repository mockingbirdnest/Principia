
#include "physics/euler_solver.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <set>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/time_scales.hpp"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/permutation.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using astronomy::ICRS;
using astronomy::operator""_UTC;
using geometry::AngleBetween;
using geometry::AngularVelocity;
using geometry::Arbitrary;
using geometry::Bivector;
using geometry::DefinesFrame;
using geometry::EulerAngles;
using geometry::EvenPermutation;
using geometry::Frame;
using geometry::Handedness;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Permutation;
using geometry::R3Element;
using geometry::RadiusLatitudeLongitude;
using geometry::Rotation;
using quantities::Abs;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::AngularMomentum;
using quantities::Cos;
using quantities::Energy;
using quantities::MomentOfInertia;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteErrorFrom;
using testing_utilities::AlmostEquals;
using testing_utilities::ApproximateQuantity;
using testing_utilities::Componentwise;
using testing_utilities::EqualsProto;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::RelativeErrorFrom;
using testing_utilities::VanishesBefore;
using testing_utilities::operator""_⑴;
using ::testing::Le;
using ::testing::Lt;
using ::testing::Matcher;
namespace si = quantities::si;

class EulerSolverTest : public ::testing::Test {
 protected:
  using PrincipalAxes = Frame<serialization::Frame::TestTag,
                              Arbitrary,
                              Handedness::Right,
                              serialization::Frame::TEST>;

  using Solver = EulerSolver<ICRS, PrincipalAxes>;

  EulerSolverTest()
      : identity_attitude_(
            EulerSolver<ICRS, PrincipalAxes>::AttitudeRotation::Identity()),
        e1_({1, 0, 0}),
        e2_({0, 1, 0}),
        e3_({0, 0, 1}) {}

  // Checks that the angular momentum transformed by the attitude to the
  // inertial frame satisfies the given matcher.
  void CheckAngularMomentumConservation(
      std::vector<Bivector<AngularMomentum, PrincipalAxes>> const&
          angular_momenta,
      std::vector<Solver::AttitudeRotation> const& attitudes,
      Matcher<Bivector<AngularMomentum, ICRS>> const& matcher) {
    CHECK_EQ(angular_momenta.size(), attitudes.size());
    for (int i = 0; i < angular_momenta.size(); ++i) {
      Bivector<AngularMomentum, ICRS> const angular_momentum_in_inertial =
          attitudes[i](angular_momenta[i]);
      EXPECT_THAT(angular_momentum_in_inertial, matcher);
    }
  }

  // Checks that the kinetic energy computed using the attitude has the right
  // value.
  void CheckPoinsotConstruction(
      Solver const& solver,
      std::vector<Bivector<AngularMomentum, PrincipalAxes>> const&
          angular_momenta,
      std::vector<Solver::AttitudeRotation> const& attitudes,
      int const min_ulps,
      int const max_ulps) {
    CHECK_EQ(angular_momenta.size(), attitudes.size());
    Energy maximum_kinetic_energy;
    Energy minimum_kinetic_energy = quantities::Infinity<Energy>();
    for (int i = 0; i < angular_momenta.size(); ++i) {
      Bivector<AngularMomentum, PrincipalAxes> const angular_momentum =
          angular_momenta[i];
      AngularVelocity<PrincipalAxes> const angular_velocity =
          solver.AngularVelocityFor(angular_momentum);
      Bivector<AngularMomentum, ICRS> const angular_momentum_in_inertial =
          attitudes[i](angular_momentum);
      AngularVelocity<ICRS> const
          angular_velocity_in_inertial = attitudes[i](angular_velocity);
      Energy const kinetic_energy = 0.5 *
                                    InnerProduct(angular_momentum_in_inertial,
                                                 angular_velocity_in_inertial) /
                                    Radian / Radian;
      maximum_kinetic_energy = std::max(maximum_kinetic_energy, kinetic_energy);
      minimum_kinetic_energy = std::min(minimum_kinetic_energy, kinetic_energy);
    }
    EXPECT_THAT(minimum_kinetic_energy,
                AlmostEquals(maximum_kinetic_energy, min_ulps, max_ulps));
  }

  Solver::AttitudeRotation identity_attitude_;
  Bivector<double, PrincipalAxes> const e1_;
  Bivector<double, PrincipalAxes> const e2_;
  Bivector<double, PrincipalAxes> const e3_;
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
        randoms[0] * si::Unit<MomentOfInertia>,
        randoms[1] * si::Unit<MomentOfInertia>,
        randoms[2] * si::Unit<MomentOfInertia>};

    Bivector<AngularMomentum, PrincipalAxes>
        initial_angular_momentum(
            {angular_momentum_distribution(random) * si::Unit<AngularMomentum>,
             angular_momentum_distribution(random) * si::Unit<AngularMomentum>,
             angular_momentum_distribution(random) *
                 si::Unit<AngularMomentum>});

    Solver const solver(moments_of_inertia,
                        identity_attitude_(initial_angular_momentum),
                        identity_attitude_,
                        Instant());
    auto const computed_initial_angular_momentum =
        solver.AngularMomentumAt(Instant());

    EXPECT_THAT(computed_initial_angular_momentum,
                AlmostEquals(initial_angular_momentum, 0, 2136))
        << moments_of_inertia << " " << initial_angular_momentum;
  }
}

// Same as above, but exercises the symmetrical cases where at least two moments
// of inertia are equal.
TEST_F(EulerSolverTest, InitialStateSymmetrical) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> angular_momentum_distribution(-10.0, 10.0);

  R3Element<MomentOfInertia> const moments_of_inertia1{
      2 * si::Unit<MomentOfInertia>,
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>};
  R3Element<MomentOfInertia> const moments_of_inertia2{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>};
  R3Element<MomentOfInertia> const moments_of_inertia3{
      3 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>};

  for (int i = 0; i < 100; ++i) {
    Bivector<AngularMomentum, PrincipalAxes>
        initial_angular_momentum(
            {angular_momentum_distribution(random) * si::Unit<AngularMomentum>,
             angular_momentum_distribution(random) * si::Unit<AngularMomentum>,
             angular_momentum_distribution(random) *
                 si::Unit<AngularMomentum>});
    {
      Solver const solver1(moments_of_inertia1,
                           identity_attitude_(initial_angular_momentum),
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
                           identity_attitude_(initial_angular_momentum),
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
                           identity_attitude_(initial_angular_momentum),
                           identity_attitude_,
                           Instant());
      auto const computed_initial_angular_momentum3 =
          solver3.AngularMomentumAt(Instant());

      EXPECT_THAT(computed_initial_angular_momentum3,
                  AlmostEquals(initial_angular_momentum, 0, 0))
          << moments_of_inertia3 << " " << initial_angular_momentum;
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
    auto const I₁ = randoms[0] * si::Unit<MomentOfInertia>;
    auto const I₂ = randoms[1] * si::Unit<MomentOfInertia>;
    auto const I₃ = randoms[2] * si::Unit<MomentOfInertia>;
    R3Element<MomentOfInertia> const moments_of_inertia{I₁, I₂, I₃};

    // G² = T * (I₁ + I₂)
    {
      auto const mx =
          angular_momentum_distribution(random) * si::Unit<AngularMomentum>;
      auto mz = mx * Sqrt(((I₂ - I₁) * I₃) / ((2.0 * I₃ - I₂ - I₁) * I₁));
      if (i % 2 == 0) {
        mz = -mz;
      }
      Bivector<AngularMomentum, PrincipalAxes>
          initial_angular_momentum({mx, si::Unit<AngularMomentum>, mz});
      Solver const solver(moments_of_inertia,
                          identity_attitude_(initial_angular_momentum),
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
          angular_momentum_distribution(random) * si::Unit<AngularMomentum>;
      auto mz = mx * Sqrt(((I₂ - I₁) * I₃) / ((I₃ - I₂) * I₁));
      if (i % 2 == 0) {
        mz = -mz;
      }
      Bivector<AngularMomentum, PrincipalAxes>
          initial_angular_momentum({mx, si::Unit<AngularMomentum>, mz});
      Solver const solver(moments_of_inertia,
                          identity_attitude_(initial_angular_momentum),
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
          angular_momentum_distribution(random) * si::Unit<AngularMomentum>;
      auto mz = mx * Sqrt(((I₂ + I₃ - 2.0 * I₁) * I₃) / ((I₃ - I₂) * I₁));
      if (i % 2 == 0) {
        mz = -mz;
      }
      Bivector<AngularMomentum, PrincipalAxes>
          initial_angular_momentum({mx, si::Unit<AngularMomentum>, mz});
      Solver const solver(moments_of_inertia,
                          identity_attitude_(initial_angular_momentum),
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
      3.0 * si::Unit<MomentOfInertia>,
      3.0 * si::Unit<MomentOfInertia>,
      9.0 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const
      initial_angular_momentum({0.0 * si::Unit<AngularMomentum>,
                                5.0 * si::Unit<AngularMomentum>,
                                7.0 * si::Unit<AngularMomentum>});
  Solver::AttitudeRotation const initial_attitude = identity_attitude_;

  // Correspondence with the referential of lecture 19: x = e1, y = e2, z = e3.
  AngularFrequency Ω = initial_angular_momentum.coordinates().z *
                       (moments_of_inertia[0] - moments_of_inertia[2]) /
                       (moments_of_inertia[0] * moments_of_inertia[2]);

  Solver const solver(moments_of_inertia,
                      initial_attitude(initial_angular_momentum),
                      initial_attitude,
                      Instant());

  std::vector<Bivector<AngularMomentum, PrincipalAxes>> angular_momenta;
  std::vector<Solver::AttitudeRotation> attitudes;
  for (Time t = 0 * Second; t < 5.0 * Second; t += 0.1 * Second) {
    auto const angular_momentum_at_t = solver.AngularMomentumAt(Instant() + t);
    EXPECT_THAT(
        angular_momentum_at_t,
        AlmostEquals(
            Bivector<AngularMomentum, PrincipalAxes>(
                {5.0 * Sin(Ω * t) * si::Unit<AngularMomentum>,
                 5.0 * Cos(Ω * t) * si::Unit<AngularMomentum>,
                 7.0 * si::Unit<AngularMomentum>}),
            0,
            102))
        << t;
    angular_momenta.push_back(angular_momentum_at_t);
    attitudes.push_back(solver.AttitudeAt(angular_momentum_at_t,
                                          Instant() + t));
  }

  Bivector<AngularMomentum, ICRS> const reference_momentum =
      initial_attitude(initial_angular_momentum);
  CheckAngularMomentumConservation(
      angular_momenta,
      attitudes,
      Componentwise(VanishesBefore(1 * si::Unit<AngularMomentum>, 0, 32),
                    AlmostEquals(reference_momentum.coordinates().y, 0, 4),
                    AlmostEquals(reference_momentum.coordinates().z, 0, 2)));
  CheckPoinsotConstruction(solver,
                           angular_momenta,
                           attitudes,
                           /*min_ulps=*/8,
                           /*max_ulps=*/8);
}

TEST_F(EulerSolverTest, TallSkinnySymmetricTopPrecession) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      3.0 * si::Unit<MomentOfInertia>,
      9.0 * si::Unit<MomentOfInertia>,
      9.0 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const
      initial_angular_momentum({7.0 * si::Unit<AngularMomentum>,
                                0.0 * si::Unit<AngularMomentum>,
                                5.0 * si::Unit<AngularMomentum>});
  Solver::AttitudeRotation const initial_attitude = identity_attitude_;

  // Correspondence with the referential of lecture 19:  x = e3, y = e1, z = e2.
  AngularFrequency Ω = initial_angular_momentum.coordinates().x *
                       (moments_of_inertia[1] - moments_of_inertia[0]) /
                       (moments_of_inertia[1] * moments_of_inertia[0]);

  Solver const solver(moments_of_inertia,
                      initial_attitude(initial_angular_momentum),
                      initial_attitude,
                      Instant());

  std::vector<Bivector<AngularMomentum, PrincipalAxes>> angular_momenta;
  std::vector<Solver::AttitudeRotation> attitudes;
  for (Time t = 0 * Second; t < 5.0 * Second; t += 0.1 * Second) {
    auto const angular_momentum_at_t = solver.AngularMomentumAt(Instant() + t);
    EXPECT_THAT(
        angular_momentum_at_t,
        AlmostEquals(
            Bivector<AngularMomentum, PrincipalAxes>(
                {7.0 * si::Unit<AngularMomentum>,
                 5.0 * Sin(Ω * t) * si::Unit<AngularMomentum>,
                 5.0 * Cos(Ω * t) * si::Unit<AngularMomentum>}),
            0,
            34))
        << t;
    angular_momenta.push_back(angular_momentum_at_t);
    attitudes.push_back(solver.AttitudeAt(angular_momentum_at_t,
                                          Instant() + t));
  }

  Bivector<AngularMomentum, ICRS> const reference_momentum =
      initial_attitude(initial_angular_momentum);
  CheckAngularMomentumConservation(
      angular_momenta,
      attitudes,
      Componentwise(AlmostEquals(reference_momentum.coordinates().x, 0, 3),
                    VanishesBefore(1 * si::Unit<AngularMomentum>, 0, 24),
                    AlmostEquals(reference_momentum.coordinates().z, 0, 6)));
  CheckPoinsotConstruction(solver,
                           angular_momenta,
                           attitudes,
                           /*min_ulps=*/2,
                           /*max_ulps=*/2);
}

// This test demonstrates the Джанибеков effect, also known as tennis racket
// theorem: the rotation of an object around its second principal axis is not
// stable.  Here we choose the initial angular momentum to be mostly in the y
// direction with a small component in the z direction.  This causes the object
// to periodically flip, rotating along y or along -y.
TEST_F(EulerSolverTest, ДжанибековEffect) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      3.0 * si::Unit<MomentOfInertia>,
      5.0 * si::Unit<MomentOfInertia>,
      9.0 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const
      initial_angular_momentum({0.0 * si::Unit<AngularMomentum>,
                                2.0 * si::Unit<AngularMomentum>,
                                0.01 * si::Unit<AngularMomentum>});
  Solver::AttitudeRotation const initial_attitude = identity_attitude_;

  Solver const solver(moments_of_inertia,
                      initial_attitude(initial_angular_momentum),
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
  std::vector<Bivector<AngularMomentum, PrincipalAxes>> angular_momenta;
  std::vector<Solver::AttitudeRotation> attitudes;
  for (Instant t; t < Instant() + 100.0 * Second; t += 0.1 * Second) {
    auto const angular_momentum_at_t = solver.AngularMomentumAt(t);
    auto const my = angular_momentum_at_t.coordinates().y;
    if (is_increasing && my > 1.99 * si::Unit<AngularMomentum>) {
      if (my < previous_my) {
        maxima.push_back(previous_t);
        is_abs_decreasing = true;
        is_decreasing = true;
        is_increasing = false;
      }
    }
    if (is_decreasing && my < -1.99 * si::Unit<AngularMomentum>) {
      if (my > previous_my) {
        minima.push_back(previous_t);
        is_abs_decreasing = true;
        is_decreasing = false;
        is_increasing = true;
      }
    }
    if (is_abs_decreasing && Abs(my) < 0.1 * si::Unit<AngularMomentum>) {
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

  Bivector<AngularMomentum, ICRS> const reference_momentum =
      initial_attitude(initial_angular_momentum);
  CheckAngularMomentumConservation(
      angular_momenta,
      attitudes,
      Componentwise(VanishesBefore(1 * si::Unit<AngularMomentum>, 0, 10),
                    AlmostEquals(reference_momentum.coordinates().y, 0, 16),
                    AlmostEquals(reference_momentum.coordinates().z, 0, 965)));
  CheckPoinsotConstruction(solver,
                           angular_momenta,
                           attitudes,
                           /*min_ulps=*/39,
                           /*max_ulps=*/43);
}

// A general body that doesn't rotate.
TEST_F(EulerSolverTest, GeneralBodyNoRotation) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      5 * si::Unit<MomentOfInertia>};

  Solver const solver(moments_of_inertia,
                      Bivector<AngularMomentum, ICRS>(),
                      identity_attitude_,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  EXPECT_THAT(actual_angular_momentum,
              AlmostEquals(Bivector<AngularMomentum, PrincipalAxes>(), 0));
  EXPECT_THAT(actual_angular_velocity,
              AlmostEquals(AngularVelocity<PrincipalAxes>(), 0));
  EXPECT_THAT(actual_attitude(e1_), AlmostEquals(identity_attitude_(e1_), 0));
  EXPECT_THAT(actual_attitude(e2_), AlmostEquals(identity_attitude_(e2_), 0));
  EXPECT_THAT(actual_attitude(e3_), AlmostEquals(identity_attitude_(e3_), 0));
}

// A general body that rotates along the first principal axis.
TEST_F(EulerSolverTest, GeneralBodyRotationAlongFirstAxis) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      5 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {5.0 * si::Unit<AngularMomentum>,
       0.0 * si::Unit<AngularMomentum>,
       0.0 * si::Unit<AngularMomentum>});
  Solver const solver(moments_of_inertia,
                      identity_attitude_(initial_angular_momentum),
                      identity_attitude_,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  auto const expected_angular_frequency = actual_angular_velocity.Norm();
  auto const expected_attitude =
      identity_attitude_ *
      Rotation<PrincipalAxes, PrincipalAxes>(
          expected_angular_frequency * 10 * Second, initial_angular_momentum);

  EXPECT_THAT(actual_angular_momentum,
              AlmostEquals(initial_angular_momentum, 0));
  EXPECT_THAT(
      actual_angular_velocity,
      AlmostEquals(solver.AngularVelocityFor(initial_angular_momentum), 0));
  EXPECT_THAT(actual_attitude(e1_), AlmostEquals(expected_attitude(e1_), 0));
  EXPECT_THAT(actual_attitude(e2_), AlmostEquals(expected_attitude(e2_), 17));
  EXPECT_THAT(actual_attitude(e3_), AlmostEquals(expected_attitude(e3_), 17));
}

// A general body that rotates along the second principal axis.
TEST_F(EulerSolverTest, GeneralBodyRotationAlongSecondAxis) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      5 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {0.0 * si::Unit<AngularMomentum>,
       6.0 * si::Unit<AngularMomentum>,
       0.0 * si::Unit<AngularMomentum>});
  Solver const solver(moments_of_inertia,
                      identity_attitude_(initial_angular_momentum),
                      identity_attitude_,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  auto const expected_angular_frequency = actual_angular_velocity.Norm();
  auto const expected_attitude =
      identity_attitude_ *
      Rotation<PrincipalAxes, PrincipalAxes>(
          expected_angular_frequency * 10 * Second, initial_angular_momentum);

  EXPECT_THAT(actual_angular_momentum,
              AlmostEquals(initial_angular_momentum, 0));
  EXPECT_THAT(
      actual_angular_velocity,
      AlmostEquals(solver.AngularVelocityFor(initial_angular_momentum), 0));
  EXPECT_THAT(
      actual_attitude(e1_),
      Componentwise(AlmostEquals(expected_attitude(e1_).coordinates().x, 0),
                    VanishesBefore(1, 0),
                    AlmostEquals(expected_attitude(e1_).coordinates().z, 1)));
  EXPECT_THAT(
      actual_attitude(e2_),
      Componentwise(VanishesBefore(1, 1),
                    AlmostEquals(expected_attitude(e2_).coordinates().y, 0),
                    VanishesBefore(1, 1)));
  EXPECT_THAT(
      actual_attitude(e3_),
      Componentwise(AlmostEquals(expected_attitude(e3_).coordinates().x, 1),
                    VanishesBefore(1, 2),
                    AlmostEquals(expected_attitude(e3_).coordinates().z, 0)));
}

// A general body that rotates along the third principal axis.
TEST_F(EulerSolverTest, GeneralBodyRotationAlongThirdAxis) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      5 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {0.0 * si::Unit<AngularMomentum>,
       0.0 * si::Unit<AngularMomentum>,
       7.0 * si::Unit<AngularMomentum>});
  Solver const solver(moments_of_inertia,
                      identity_attitude_(initial_angular_momentum),
                      identity_attitude_,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  auto const expected_angular_frequency = actual_angular_velocity.Norm();
  auto const expected_attitude =
      identity_attitude_ *
      Rotation<PrincipalAxes, PrincipalAxes>(
          expected_angular_frequency * 10 * Second, initial_angular_momentum);

  EXPECT_THAT(actual_angular_momentum,
              AlmostEquals(initial_angular_momentum, 0));
  EXPECT_THAT(
      actual_angular_velocity,
      AlmostEquals(solver.AngularVelocityFor(initial_angular_momentum), 0));
  EXPECT_THAT(actual_attitude(e1_), AlmostEquals(expected_attitude(e1_), 136));
  EXPECT_THAT(actual_attitude(e2_), AlmostEquals(expected_attitude(e2_), 136));
  EXPECT_THAT(actual_attitude(e3_), AlmostEquals(expected_attitude(e3_), 0));
}

// A general body that has an angular momentum close to the third principal
// axis.
TEST_F(EulerSolverTest, GeneralBodyRotationCloseToThirdAxis) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      5 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {std::numeric_limits<double>::epsilon() * si::Unit<AngularMomentum>,
       0.0 * si::Unit<AngularMomentum>,
       1.0 * si::Unit<AngularMomentum>});
  Solver const solver(moments_of_inertia,
                      identity_attitude_(initial_angular_momentum),
                      identity_attitude_,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  // The expected attitude ignores any precession and is just a rotation
  // around z.
  auto const expected_angular_frequency = actual_angular_velocity.Norm();
  auto const expected_attitude =
      identity_attitude_ * Rotation<PrincipalAxes, PrincipalAxes>(
                               expected_angular_frequency * 10 * Second,
                               Bivector<AngularMomentum, PrincipalAxes>(
                                   {0.0 * si::Unit<AngularMomentum>,
                                    0.0 * si::Unit<AngularMomentum>,
                                    1.0 * si::Unit<AngularMomentum>}));

  EXPECT_THAT(
      actual_attitude(e1_),
      Componentwise(
          AlmostEquals(expected_attitude(e1_).coordinates().x, 8, 12),
          AlmostEquals(expected_attitude(e1_).coordinates().y, 0, 2),
          VanishesBefore(std::numeric_limits<double>::epsilon(), 4, 8)));
  EXPECT_THAT(
      actual_attitude(e2_),
      Componentwise(AlmostEquals(expected_attitude(e2_).coordinates().x, 0, 2),
                    AlmostEquals(expected_attitude(e2_).coordinates().y, 8, 12),
                    VanishesBefore(1, 2)));
  EXPECT_THAT(
      actual_attitude(e3_),
      Componentwise(VanishesBefore(1, 2),
                    VanishesBefore(1, 1),
                    AlmostEquals(expected_attitude(e3_).coordinates().z, 0)));
}

// A general body that has an angular momentum very close to the third principal
// axis.
TEST_F(EulerSolverTest, GeneralBodyRotationVeryCloseToThirdAxis) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      5 * si::Unit<MomentOfInertia>};

  // There is a small region from 3.4 to 4.0 in the expression below where the
  // check fails.
  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {3.5 * Sqrt(std::numeric_limits<double>::denorm_min()) *
           si::Unit<AngularMomentum>,
       0.0 * si::Unit<AngularMomentum>,
       1.0 * si::Unit<AngularMomentum>});
  Solver const solver(moments_of_inertia,
                      identity_attitude_(initial_angular_momentum),
                      identity_attitude_,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  // The expected attitude ignores any precession and is just a rotation
  // around z.
  auto const expected_angular_frequency = actual_angular_velocity.Norm();
  auto const expected_attitude =
      identity_attitude_ * Rotation<PrincipalAxes, PrincipalAxes>(
                               expected_angular_frequency * 10 * Second,
                               Bivector<AngularMomentum, PrincipalAxes>(
                                   {0.0 * si::Unit<AngularMomentum>,
                                    0.0 * si::Unit<AngularMomentum>,
                                    1.0 * si::Unit<AngularMomentum>}));
  EXPECT_THAT(
      actual_attitude(e1_),
      Componentwise(AlmostEquals(expected_attitude(e1_).coordinates().x, 0, 12),
                    AlmostEquals(expected_attitude(e1_).coordinates().y, 0, 2),
                    VanishesBefore(1, 0)));
  EXPECT_THAT(
      actual_attitude(e2_),
      Componentwise(AlmostEquals(expected_attitude(e2_).coordinates().x, 0, 2),
                    AlmostEquals(expected_attitude(e2_).coordinates().y, 0, 12),
                    VanishesBefore(1, 0)));
  EXPECT_THAT(actual_attitude(e3_),
              Componentwise(
                  VanishesBefore(1, 0),
                  VanishesBefore(1, 0),
                  AlmostEquals(expected_attitude(e3_).coordinates().z, 0, 0)));
}

// A sphere that rotates around an axis with a random orientation.  Also, a
// non-trivial initial attitude.
TEST_F(EulerSolverTest, SphereRotationAlongRandomAxis) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      3 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {1.0 * si::Unit<AngularMomentum>,
       2.0 * si::Unit<AngularMomentum>,
       4.0 * si::Unit<AngularMomentum>});
  Solver::AttitudeRotation const initial_attitude(
      0.5 * Radian,
      1.5 * Radian,
      -2.5 * Radian,
      EulerAngles::YZY,
      DefinesFrame<PrincipalAxes>{});
  Solver const solver(moments_of_inertia,
                      initial_attitude(initial_angular_momentum),
                      initial_attitude,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  auto const expected_angular_velocity =
      AngularVelocity<PrincipalAxes>({1.0 / 3.0 * Radian / Second,
                                      2.0 / 3.0 * Radian / Second,
                                      4.0 / 3.0 * Radian / Second});
  auto const expected_angular_frequency = expected_angular_velocity.Norm();
  auto const expected_attitude =
      initial_attitude *
      Rotation<PrincipalAxes, PrincipalAxes>(
          expected_angular_frequency * 10 * Second, initial_angular_momentum);

  EXPECT_THAT(actual_angular_momentum,
              AlmostEquals(initial_angular_momentum, 8, 16));
  EXPECT_THAT(actual_angular_velocity,
              AlmostEquals(expected_angular_velocity, 5, 10));
  EXPECT_THAT(actual_attitude(e1_),
              AlmostEquals(expected_attitude(e1_), 31, 60));
  EXPECT_THAT(actual_attitude(e2_),
              AlmostEquals(expected_attitude(e2_), 44, 114));
  EXPECT_THAT(actual_attitude(e3_),
              AlmostEquals(expected_attitude(e3_), 9, 20));
}

// Rotation on the separatrix with a constant momentum.
TEST_F(EulerSolverTest, SeparatrixConstantMomentum) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      2 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>,
      3 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {0.0 * si::Unit<AngularMomentum>,
       4.0 * si::Unit<AngularMomentum>,
       5.0 * si::Unit<AngularMomentum>});
  Solver const solver(moments_of_inertia,
                      identity_attitude_(initial_angular_momentum),
                      identity_attitude_,
                      Instant());

  auto const t = Instant() + 10 * Second;
  auto const actual_angular_momentum = solver.AngularMomentumAt(t);
  auto const actual_angular_velocity =
      solver.AngularVelocityFor(actual_angular_momentum);
  auto const actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

  auto const expected_angular_velocity =
      AngularVelocity<PrincipalAxes>({0.0 * Radian / Second,
                                      4.0 / 3.0 * Radian / Second,
                                      5.0 / 3.0 * Radian / Second});
  auto const expected_angular_frequency = expected_angular_velocity.Norm();
  auto const expected_attitude =
      identity_attitude_ *
      Rotation<PrincipalAxes, PrincipalAxes>(
          expected_angular_frequency * 10 * Second, initial_angular_momentum);

  EXPECT_THAT(actual_angular_momentum,
              AlmostEquals(initial_angular_momentum, 1));
  EXPECT_THAT(actual_angular_velocity,
              AlmostEquals(expected_angular_velocity, 2));
  EXPECT_THAT(actual_attitude(e1_), AlmostEquals(expected_attitude(e1_), 8));
  EXPECT_THAT(actual_attitude(e2_), AlmostEquals(expected_attitude(e2_), 4, 6));
  EXPECT_THAT(actual_attitude(e3_), AlmostEquals(expected_attitude(e3_), 6, 7));
}

// The data in this test are from Takahashi, Busch and Scheeres, Spin state and
// moment of inertia characterization of 4179 Toutatis, 2013 [TBS13].
TEST_F(EulerSolverTest, Toutatis) {
  Instant const epoch = "1992-11-09T17:49:47"_UTC;

  // [TBS13] adopt a bizarre convention where their x axis is our I₂, their y
  // axis is our I₃ and their z axis is our I₁, see Table 2.  This appears to
  // contradict Figure 1, but it is consistent with ω₁, ω₂, ω₃ being along their
  // x, y, z axes respectively.
  using TakahashiPrincipalAxes = Frame<enum class TakahashiPrincipalAxesTag>;
  using TakahashiAttitudeRotation = Rotation<TakahashiPrincipalAxes, ICRS>;
  using TakahashiPermutation = Permutation<TakahashiPrincipalAxes,
                                           PrincipalAxes>;
  TakahashiPermutation const takahashi_to_vanilla(EvenPermutation::ZXY);

  R3Element<MomentOfInertia> const takahashi_moments_of_inertia{
      3.0836 * si::Unit<MomentOfInertia>,
      3.235 * si::Unit<MomentOfInertia>,
      1 * si::Unit<MomentOfInertia>};

  TakahashiAttitudeRotation const takahashi_initial_attitude(
      /*α=*/145.498 * Degree,
      /*β=*/65.865 * Degree,
      /*γ=*/241.524 * Degree,
      EulerAngles::ZXZ,
      DefinesFrame<TakahashiPrincipalAxes>{});

  AngularVelocity<TakahashiPrincipalAxes> const
      takahashi_initial_angular_velocity({14.51 * Degree / Day,
                                          33.529 * Degree / Day,
                                          -98.709 * Degree / Day});

  Bivector<AngularMomentum, TakahashiPrincipalAxes> const
      takahashi_initial_angular_momentum(
          {takahashi_initial_angular_velocity.coordinates().x *
               takahashi_moments_of_inertia.x,
           takahashi_initial_angular_velocity.coordinates().y *
               takahashi_moments_of_inertia.y,
           takahashi_initial_angular_velocity.coordinates().z *
               takahashi_moments_of_inertia.z});

  // From Zhao et al., Orientation and rotational parameters of asteroid 4179
  // Toutatis: new insights from Change'e-2's close flyby [Zha+15].
  auto const angular_momentum_orientation_in_inertial = Bivector<double, ICRS>(
      RadiusLatitudeLongitude(1.0, -54.75 * Degree, 180.2 * Degree)
          .ToCartesian());

  Solver::AttitudeRotation const initial_attitude =
      takahashi_initial_attitude *
      takahashi_to_vanilla.Inverse().Forget<Rotation>();
  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum =
      takahashi_to_vanilla(takahashi_initial_angular_momentum);

  // Check that the angular momentum in the body frame is consistent with that
  // in the inertial frame transformed by the initial attitude.  The data come
  // from different sources, so this is not a trivial check.
  Bivector<double, PrincipalAxes> const angular_momentum_orientation_in_body =
      initial_attitude.Inverse()(angular_momentum_orientation_in_inertial);
  EXPECT_THAT(
      Normalize(initial_angular_momentum),
      Componentwise(RelativeErrorFrom(
                        angular_momentum_orientation_in_body.coordinates().x,
                        IsNear(0.016_⑴)),
                    RelativeErrorFrom(
                        angular_momentum_orientation_in_body.coordinates().y,
                        IsNear(0.061_⑴)),
                    RelativeErrorFrom(
                        angular_momentum_orientation_in_body.coordinates().z,
                        IsNear(0.001_⑴))));

  // Same check as above, but in the inertial frame.  The y coordinate is small
  // because the longitute is close to 180°, so we only check the absolute error
  // for it.
  EXPECT_THAT(Normalize(initial_attitude(initial_angular_momentum)),
              Componentwise(
                  RelativeErrorFrom(
                      angular_momentum_orientation_in_inertial.coordinates().x,
                      IsNear(0.028_⑴)),
                  AbsoluteErrorFrom(
                      angular_momentum_orientation_in_inertial.coordinates().y,
                      IsNear(0.008_⑴)),
                  RelativeErrorFrom(
                      angular_momentum_orientation_in_inertial.coordinates().z,
                      IsNear(0.014_⑴))));

  Solver const solver(takahashi_to_vanilla(takahashi_moments_of_inertia),
                      initial_attitude(initial_angular_momentum),
                      initial_attitude,
                      epoch);

  struct Observation {
    Instant t;
    Angle α;
    Angle β;
    Angle γ;
    AngularFrequency ω₁;
    AngularFrequency ω₂;
    AngularFrequency ω₃;
    ApproximateQuantity<double> angular_velocity_norm_error;
    ApproximateQuantity<Angle> angular_velocity_direction_error;
    ApproximateQuantity<Angle> e1_direction_error;
    ApproximateQuantity<Angle> e2_direction_error;
    ApproximateQuantity<Angle> e3_direction_error;
  };

  // The errors for 1992 seem consistent with the residuals from [TBS13] (3 σ
  // for the attitude on 1992-12-04, 1.5 σ on 1992-12-07, for instance, figures
  // 7 and 8) and they are generally smaller for the angular velocity than for
  // the attitude, in agreement with that paper.
  // The errors for 2008 are caused by the fact that we ignore the torque
  // exerted by the Earth and the Sun.  As shown in figure 6 of [TBS13], this
  // results in a completely different orientation, estimated in the text to be
  // around 100°, which is again generally consistent with our results.
  std::vector<Observation> const observations = {
      {"1992-12-02T21:40:00"_UTC,
       122.2 * Degree, 86.5 * Degree, 107.0 * Degree,
       -35.6 * Degree / Day, 7.2 * Degree / Day, -97.0 * Degree / Day,

       0.011_⑴, 4.5_⑴ * Degree,
       7.9_⑴ * Degree, 5.4_⑴ * Degree, 5.9_⑴ * Degree},
      {"1992-12-03T19:30:00"_UTC,
       86.3 * Degree, 81.8 * Degree, 24.5 * Degree,
       -16.4 * Degree / Day, -29.1 * Degree / Day, -91.9 * Degree / Day,

       0.076_⑴, 0.5_⑴ * Degree,
       6.6_⑴ * Degree, 7.7_⑴ * Degree, 4.7_⑴ * Degree},
      {"1992-12-04T18:10:00"_UTC,
       47.8 * Degree, 60.7 * Degree, 284.0 * Degree,
       29.1 * Degree / Day, -23.2 * Degree / Day, -97.8 * Degree / Day,

       0.005_⑴, 5.0_⑴ * Degree,
       13_⑴ * Degree, 15_⑴ * Degree, 9.4_⑴ * Degree},
      {"1992-12-05T18:50:00"_UTC,
       14.6 * Degree, 39.4 * Degree, 207.1 * Degree,
       33.3 * Degree / Day, 8.2 * Degree / Day, -92.2 * Degree / Day,

       0.064_⑴, 1.1_⑴ * Degree,
       9.4_⑴ * Degree, 4.3_⑴ * Degree, 8.4_⑴ * Degree},
      {"1992-12-06T17:30:00"_UTC,
       331.3 * Degree, 23.7 * Degree, 151.6 * Degree,
       6.6 * Degree / Day, 34.5 * Degree / Day, -95.8 * Degree / Day,

       0.032_⑴, 0.8_⑴ * Degree,
       2.0_⑴ * Degree, 0.6_⑴ * Degree, 2.1_⑴ * Degree},
      {"1992-12-07T17:20:00"_UTC,
       222.5 * Degree, 25.4 * Degree, 143.9 * Degree,
       12.8 * Degree / Day, 25.4 * Degree / Day, -104.1 * Degree / Day,

       0.028_⑴, 24_⑴ * Degree,
       5.7_⑴ * Degree, 15_⑴ * Degree, 15_⑴ * Degree},
      {"1992-12-08T16:40:00"_UTC,
       169.8 * Degree, 45.5 * Degree, 106.9 * Degree,
       -31.1 * Degree / Day, -21.9 * Degree / Day, -97.7 * Degree / Day,

       0.0001_⑴, 2.7_⑴ * Degree,
       3.0_⑴ * Degree, 2.6_⑴ * Degree, 2.6_⑴ * Degree},
      {"1992-12-09T17:50:00"_UTC,
       137.3 * Degree, 71.3 * Degree, 22.3 * Degree,
       11.8 * Degree / Day, -36.9 * Degree / Day, -94.9 * Degree / Day,

       0.028_⑴, 3.4_⑴ * Degree,
       3.4_⑴ * Degree, 3.4_⑴ * Degree, 1.1_⑴ * Degree},
      {"1992-12-10T17:20:00"_UTC,
       103.1 * Degree, 85.2 * Degree, 292.6 * Degree,
       35.8 * Degree / Day, -8.9 * Degree / Day, -97.9 * Degree / Day,

       0.0009_⑴, 0.8_⑴ * Degree,
       3.1_⑴ * Degree, 1.4_⑴ * Degree, 3.4_⑴ * Degree},
      {"1992-12-11T09:40:00"_UTC,
       77.0 * Degree, 85.7 * Degree, 225.5 * Degree,
       31.0 * Degree / Day, 17.0 * Degree / Day, -96.3 * Degree / Day,

       0.022_⑴, 1.0_⑴ * Degree,
       2.8_⑴ * Degree, 1.6_⑴ * Degree, 2.7_⑴ * Degree},
      {"1992-12-12T09:20:00"_UTC,
       42.8 * Degree, 70.2 * Degree, 133.2 * Degree,
       -1.3 * Degree / Day, 37.0 * Degree / Day, -95.9 * Degree / Day,

       0.025_⑴, 2.2_⑴ * Degree,
       2.0_⑴ * Degree, 3.2_⑴ * Degree, 2.6_⑴ * Degree},
      {"1992-12-13T08:10:00"_UTC,
       13.7 * Degree, 44.4 * Degree, 51.9 * Degree,
       -38.3 * Degree / Day, 17.9 * Degree / Day, -97.3 * Degree / Day,

       0.013_⑴, 3.4_⑴ * Degree,
       2.5_⑴ * Degree, 5.3_⑴ * Degree, 5.3_⑴ * Degree},
      {"1992-12-14T07:50:00"_UTC,
       323.7 * Degree, 14.0 * Degree, 0.0 * Degree,
       -70.5 * Degree / Day, -30.6 * Degree / Day, -91.1 * Degree / Day,

       0.119_⑴, 22_⑴ * Degree,
       4.1_⑴ * Degree, 3.8_⑴ * Degree, 4.5_⑴ * Degree},
      {"1992-12-15T07:50:00"_UTC,
       193.2 * Degree, 24.4 * Degree, 21.4 * Degree,
       22.1 * Degree / Day, -26.6 * Degree / Day, -96.6 * Degree / Day,

       0.026_⑴, 4.9_⑴ * Degree,
       3.2_⑴ * Degree, 8.2_⑴ * Degree, 8.0_⑴ * Degree},
      {"1992-12-16T07:10:00"_UTC,
       165.1 * Degree, 46.4 * Degree, 310.6 * Degree,
       33.4 * Degree / Day, -3.4 * Degree / Day, -93.7 * Degree / Day,

       0.052_⑴, 2.6_⑴ * Degree,
       7.3_⑴ * Degree, 5.8_⑴ * Degree, 8.0_⑴ * Degree},
      {"1992-12-17T06:49:00"_UTC,
       130.6 * Degree, 76.1 * Degree, 234.9 * Degree,
       12.6 * Degree / Day, 33.9 * Degree / Day, -94.0 * Degree / Day,

       0.045_⑴, 2.1_⑴ * Degree,
       0.8_⑴ * Degree, 1.1_⑴ * Degree, 1.2_⑴ * Degree},
      {"1992-12-18T07:09:00"_UTC,
       91.6 * Degree, 81.6 * Degree, 142.4 * Degree,
       -24.3 * Degree / Day, 29.6 * Degree / Day, -102.0 * Degree / Day,

       0.036_⑴, 2.1_⑴ * Degree,
       5.6_⑴ * Degree, 6.3_⑴ * Degree, 4.1_⑴ * Degree},
      {"2008-11-23T10:45:00"_UTC,
       86.2 * Degree, 85.0 * Degree, 0.3 * Degree,
       -0.4 * Degree / Day, -36.2 * Degree / Day, -98.9 * Degree / Day,

       0.0005_⑴, 40_⑴ * Degree,
       29_⑴ * Degree, 122_⑴ * Degree, 131_⑴ * Degree}};

  for (auto const& observation : observations) {
    Instant const& t = observation.t;
    TakahashiAttitudeRotation const takahashi_expected_attitude(
        /*α=*/observation.α,
        /*β=*/observation.β,
        /*γ=*/observation.γ,
        EulerAngles::ZXZ,
        DefinesFrame<TakahashiPrincipalAxes>{});
    AngularVelocity<TakahashiPrincipalAxes> const
        takahashi_expected_angular_velocity(
            {observation.ω₁, observation.ω₂, observation.ω₃});

    Solver::AttitudeRotation const expected_attitude =
        takahashi_expected_attitude *
        takahashi_to_vanilla.Inverse().Forget<Rotation>();
    AngularVelocity<PrincipalAxes> const expected_angular_velocity =
        takahashi_to_vanilla(takahashi_expected_angular_velocity);

    auto actual_angular_momentum = solver.AngularMomentumAt(t);
    auto actual_angular_velocity =
        solver.AngularVelocityFor(actual_angular_momentum);
    auto actual_attitude = solver.AttitudeAt(actual_angular_momentum, t);

    EXPECT_THAT(
        actual_angular_velocity.Norm(),
        RelativeErrorFrom(expected_angular_velocity.Norm(),
                          IsNear(observation.angular_velocity_norm_error)));
    EXPECT_THAT(
        AngleBetween(actual_angular_velocity, expected_angular_velocity),
        IsNear(observation.angular_velocity_direction_error));

    EXPECT_THAT(AngleBetween(actual_attitude(e1_), expected_attitude(e1_)),
                IsNear(observation.e1_direction_error));
    EXPECT_THAT(AngleBetween(actual_attitude(e2_), expected_attitude(e2_)),
                IsNear(observation.e2_direction_error));
    EXPECT_THAT(AngleBetween(actual_attitude(e3_), expected_attitude(e3_)),
                IsNear(observation.e3_direction_error));
  }
}

TEST_F(EulerSolverTest, Serialization) {
  R3Element<MomentOfInertia> const moments_of_inertia{
      3.0 * si::Unit<MomentOfInertia>,
      5.0 * si::Unit<MomentOfInertia>,
      9.0 * si::Unit<MomentOfInertia>};

  Bivector<AngularMomentum, PrincipalAxes> const initial_angular_momentum(
      {0.0 * si::Unit<AngularMomentum>,
       2.0 * si::Unit<AngularMomentum>,
       0.01 * si::Unit<AngularMomentum>});
  Solver::AttitudeRotation const initial_attitude = identity_attitude_;

  Solver const solver1(moments_of_inertia,
                       initial_attitude(initial_angular_momentum),
                       initial_attitude,
                       Instant() + 3 * Second);

  serialization::EulerSolver message1;
  solver1.WriteToMessage(&message1);

  auto const solver2 = Solver::ReadFromMessage(message1);

  EXPECT_EQ(solver1.AngularMomentumAt(Instant() + 5 * Second),
            solver2.AngularMomentumAt(Instant() + 5 * Second));

  serialization::EulerSolver message2;
  solver2.WriteToMessage(&message2);

  EXPECT_THAT(message2, EqualsProto(message1));
}

}  // namespace physics
}  // namespace principia
