#include "physics/barycentric_rotating_dynamic_frame.hpp"

#include <memory>

#include "astronomy/frames.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::Bivector;
using geometry::Instant;
using geometry::Rotation;
using geometry::Vector;
using quantities::Time;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using ::testing::Lt;

namespace physics {

namespace {

constexpr char kBig[] = "Big";
constexpr char kSmall[] = "Small";

}  // namespace

class BarycentricRotatingDynamicFrameTest : public ::testing::Test {
 protected:
  // The rotating frame centred on the barycentre of the two bodies.
  using BigSmall = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, false /*inertial*/>;

  BarycentricRotatingDynamicFrameTest()
      : period_(10 * π * sqrt(5.0 / 7.0) * Second),
        centre_of_mass_initial_state_(Position<ICRFJ2000Equator>(),
                                      Velocity<ICRFJ2000Equator>()),
        big_initial_state_(Position<ICRFJ2000Equator>(),
                           Velocity<ICRFJ2000Equator>()),
        small_initial_state_(Position<ICRFJ2000Equator>(),
                             Velocity<ICRFJ2000Equator>()) {
    solar_system_.Initialize(
        SOLUTION_DIR / "astronomy" / "gravity_model_two_bodies_test.proto.txt",
        SOLUTION_DIR / "astronomy" / "initial_state_two_bodies_test.proto.txt");
    t0_ = solar_system_.epoch();
    ephemeris_ = solar_system_.MakeEphemeris(
                     integrators::McLachlanAtela1992Order4Optimal<
                         Position<ICRFJ2000Equator>>(),
                     10 * Milli(Second),
                     1 * Milli(Metre));
    ephemeris_->Prolong(t0_ + 2 * period_);
    big_initial_state_ = solar_system_.initial_state(kBig);
    big_gravitational_parameter_ = solar_system_.gravitational_parameter(kBig);
    small_initial_state_ = solar_system_.initial_state(kSmall);
    small_gravitational_parameter_ =
        solar_system_.gravitational_parameter(kSmall);
    centre_of_mass_initial_state_ =
        Barycentre<ICRFJ2000Equator, GravitationalParameter>(
            {big_initial_state_, small_initial_state_},
            {big_gravitational_parameter_, small_gravitational_parameter_});
    big_small_frame_ =
        std::make_unique<
            BarycentricRotatingDynamicFrame<ICRFJ2000Equator, BigSmall>>(
                ephemeris_.get(),
                solar_system_.massive_body(*ephemeris_, kBig),
                solar_system_.massive_body(*ephemeris_, kSmall));
  }

  Time const period_;
  Instant t0_;
  DegreesOfFreedom<ICRFJ2000Equator> centre_of_mass_initial_state_;
  DegreesOfFreedom<ICRFJ2000Equator> big_initial_state_;
  DegreesOfFreedom<ICRFJ2000Equator> small_initial_state_;
  GravitationalParameter big_gravitational_parameter_;
  GravitationalParameter small_gravitational_parameter_;
  std::unique_ptr<BarycentricRotatingDynamicFrame<ICRFJ2000Equator, BigSmall>>
      big_small_frame_;
  std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  SolarSystem<ICRFJ2000Equator> solar_system_;
};


TEST_F(BarycentricRotatingDynamicFrameTest, ToBigSmallFrameAtTime) {
  int const kSteps = 100;

  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / kSteps) {
    auto const to_big_small_frame_at_t = big_small_frame_->ToThisFrameAtTime(t);
    DegreesOfFreedom<BigSmall> const centre_of_mass_in_big_small_at_t =
        to_big_small_frame_at_t(centre_of_mass_initial_state_);
    EXPECT_THAT(AbsoluteError(centre_of_mass_in_big_small_at_t.position() -
                                  BigSmall::origin,
                              Displacement<BigSmall>()),
                Lt(1.0E-11 * Metre));
    EXPECT_THAT(AbsoluteError(centre_of_mass_in_big_small_at_t.velocity(),
                              Velocity<BigSmall>()),
                Lt(0 * Milli(Metre) / Second));
  }
}

}  // namespace physics
}  // namespace principia
