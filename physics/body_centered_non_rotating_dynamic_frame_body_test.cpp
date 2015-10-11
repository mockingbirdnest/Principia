#include "physics/body_centered_non_rotating_dynamic_frame.hpp"

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
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::ICRFJ2000Equator;
using geometry::Bivector;
using geometry::Instant;
using geometry::Rotation;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using ::testing::Lt;

namespace physics {

namespace {

constexpr char kBig[] = "Big";
constexpr char kSmall[] = "Small";

}  // namespace

class BodyCentredNonRotatingDynamicFrameTest : public ::testing::Test {
 protected:
  // The non-rotating frame centred on the big body.
  using Big = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, false /*inertial*/>;
  using Small = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, false /*inertial*/>;

  BodyCentredNonRotatingDynamicFrameTest()
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
    big_frame_ = std::make_unique<
                     BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Big>>(
                         ephemeris_.get(),
                         solar_system_.massive_body(*ephemeris_, kBig));
    small_initial_state_ = solar_system_.initial_state(kSmall);
    small_gravitational_parameter_ =
        solar_system_.gravitational_parameter(kSmall);
    small_frame_ = std::make_unique<
                     BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator,
                                                        Small>>(
                         ephemeris_.get(),
                         solar_system_.massive_body(*ephemeris_, kSmall));
    centre_of_mass_initial_state_ =
        Barycentre<ICRFJ2000Equator, GravitationalParameter>(
            {big_initial_state_, small_initial_state_},
            {big_gravitational_parameter_, small_gravitational_parameter_});
  }

  Time const period_;
  Instant t0_;
  DegreesOfFreedom<ICRFJ2000Equator> centre_of_mass_initial_state_;
  DegreesOfFreedom<ICRFJ2000Equator> big_initial_state_;
  DegreesOfFreedom<ICRFJ2000Equator> small_initial_state_;
  GravitationalParameter big_gravitational_parameter_;
  GravitationalParameter small_gravitational_parameter_;
  std::unique_ptr<
      BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Big>> big_frame_;
  std::unique_ptr<
      BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Small>> small_frame_;
  std::unique_ptr<Ephemeris<ICRFJ2000Equator>> ephemeris_;
  SolarSystem<ICRFJ2000Equator> solar_system_;
};

TEST_F(BodyCentredNonRotatingDynamicFrameTest, ToBigFrame) {
  int const kSteps = 100;
  RelativeDegreesOfFreedom<ICRFJ2000Equator> const big_to_centre =
      centre_of_mass_initial_state_ - big_initial_state_;

  Bivector<double, ICRFJ2000Equator> const axis({0, 0, 1});
  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / kSteps) {
    auto const angle_at_t = 2 * π * (t - t0_) * Radian / period_;
    auto const to_big_frame = big_frame_->ToThisFrameAtTime(t);
    Displacement<Big> const displacement_at_t =
        Rotation<ICRFJ2000Equator, Big>(-angle_at_t, axis)(
            big_to_centre.displacement());
    Velocity<Big> const velocity_at_t =
        Rotation<ICRFJ2000Equator, Big>(-angle_at_t, axis)(
            big_to_centre.velocity());
    EXPECT_THAT(AbsoluteError(
                    to_big_frame(centre_of_mass_initial_state_).position() -
                        Big::origin,
                    displacement_at_t),
                Lt(0.2 * Milli(Metre)));
    EXPECT_THAT(AbsoluteError(
                    to_big_frame(centre_of_mass_initial_state_).velocity(),
                    velocity_at_t),
                Lt(2 * Milli(Metre) / Second));
  }
}

}  // namespace physics
}  // namespace principia
