#include "physics/body_centred_non_rotating_reference_frame.hpp"

#include <memory>

#include "astronomy/frames.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/rigid_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using ::testing::IsNull;
using ::testing::Lt;
using ::testing::Not;
using namespace principia::astronomy::_frames;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_vanishes_before;

namespace {

constexpr char big[] = "Big";
constexpr char small[] = "Small";

}  // namespace

class BodyCentredNonRotatingReferenceFrameTest : public ::testing::Test {
 protected:
  // The non-rotating frame centred on the big body.
  using Big = Frame<serialization::Frame::TestTag,
                    Arbitrary,
                    Handedness::Right,
                    serialization::Frame::TEST>;
  using Small = Frame<serialization::Frame::TestTag,
                      Arbitrary,
                      Handedness::Right,
                      serialization::Frame::TEST1>;

  BodyCentredNonRotatingReferenceFrameTest()
      : period_(10 * π * sqrt(5.0 / 7.0) * Second),
        solar_system_(SOLUTION_DIR / "astronomy" /
                          "test_gravity_model_two_bodies.proto.txt",
                      SOLUTION_DIR / "astronomy" /
                          "test_initial_state_two_bodies_circular.proto.txt"),
        t0_(solar_system_.epoch()),
        ephemeris_(solar_system_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymplecticRungeKuttaNyströmIntegrator<
                    McLachlanAtela1992Order4Optimal,
                    Ephemeris<ICRS>::NewtonianMotionEquation>(),
                /*step=*/10 * Milli(Second)))),
        big_initial_state_(solar_system_.degrees_of_freedom(big)),
        big_gravitational_parameter_(
            solar_system_.gravitational_parameter(big)),
        small_initial_state_(solar_system_.degrees_of_freedom(small)),
        small_gravitational_parameter_(
            solar_system_.gravitational_parameter(small)) {
    EXPECT_OK(ephemeris_->Prolong(t0_ + 2 * period_));
    big_frame_ =
        std::make_unique<BodyCentredNonRotatingReferenceFrame<ICRS, Big>>(
            ephemeris_.get(), solar_system_.massive_body(*ephemeris_, big));
    small_frame_ =
        std::make_unique<BodyCentredNonRotatingReferenceFrame<ICRS, Small>>(
            ephemeris_.get(), solar_system_.massive_body(*ephemeris_, small));
  }

  Time const period_;
  SolarSystem<ICRS> solar_system_;
  Instant const t0_;
  std::unique_ptr<Ephemeris<ICRS>> const ephemeris_;
  DegreesOfFreedom<ICRS> const big_initial_state_;
  GravitationalParameter const big_gravitational_parameter_;
  DegreesOfFreedom<ICRS> const small_initial_state_;
  GravitationalParameter const small_gravitational_parameter_;
  std::unique_ptr<BodyCentredNonRotatingReferenceFrame<ICRS, Big>> big_frame_;
  std::unique_ptr<BodyCentredNonRotatingReferenceFrame<ICRS, Small>>
      small_frame_;
};


// Check that the small body has the right degrees of freedom in the frame of
// the big body.
TEST_F(BodyCentredNonRotatingReferenceFrameTest, SmallBodyInBigFrame) {
  int const steps = 100;
  Bivector<double, ICRS> const axis({0, 0, 1});

  RelativeDegreesOfFreedom<ICRS> const initial_big_to_small =
      small_initial_state_ - big_initial_state_;
  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    DegreesOfFreedom<ICRS> const small_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, small).
            EvaluateDegreesOfFreedom(t);

    auto const rotation_in_big_frame_at_t = Rotation<ICRS, Big>(
        2 * π * (t - t0_) * Radian / period_, axis, DefinesFrame<Big>{});
    DegreesOfFreedom<Big> const small_in_big_frame_at_t(
        rotation_in_big_frame_at_t(initial_big_to_small.displacement()) +
            Big::origin,
        rotation_in_big_frame_at_t(initial_big_to_small.velocity()));

    auto const to_big_frame_at_t = big_frame_->ToThisFrameAtTime(t);
    EXPECT_THAT(AbsoluteError(
                    to_big_frame_at_t(small_in_inertial_frame_at_t).position(),
                    small_in_big_frame_at_t.position()),
                Lt(0.3 * Milli(Metre)));
    EXPECT_THAT(AbsoluteError(
                    to_big_frame_at_t(small_in_inertial_frame_at_t).velocity(),
                    small_in_big_frame_at_t.velocity()),
                Lt(4 * Milli(Metre) / Second));
  }
}

TEST_F(BodyCentredNonRotatingReferenceFrameTest, Inverse) {
  int const steps = 100;
  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    auto const from_big_frame_at_t = big_frame_->FromThisFrameAtTime(t);
    auto const to_big_frame_at_t = big_frame_->ToThisFrameAtTime(t);
    auto const small_initial_state_transformed_and_back =
        from_big_frame_at_t(to_big_frame_at_t(small_initial_state_));
    auto const position_coordinates =
        (small_initial_state_.position() - ICRS::origin).coordinates();
    auto const velocity_coordinates =
        small_initial_state_.velocity().coordinates();
    EXPECT_THAT(
        small_initial_state_transformed_and_back.position() - ICRS::origin,
        Componentwise(VanishesBefore(position_coordinates.y, 0),
                      AlmostEquals(position_coordinates.y, 0),
                      VanishesBefore(position_coordinates.y, 0)));
    EXPECT_THAT(small_initial_state_transformed_and_back.velocity(),
                Componentwise(AlmostEquals(velocity_coordinates.x, 0, 1),
                              VanishesBefore(velocity_coordinates.x, 0),
                              VanishesBefore(velocity_coordinates.x, 0)));
  }
}

TEST_F(BodyCentredNonRotatingReferenceFrameTest, GeometricAcceleration) {
  int const steps = 10;
  RelativeDegreesOfFreedom<ICRS> const initial_big_to_small =
      small_initial_state_ - big_initial_state_;
  Length const big_to_small = initial_big_to_small.displacement().Norm();
  Acceleration const small_on_big =
      small_gravitational_parameter_ / (big_to_small * big_to_small);
  for (Length y = big_to_small / steps;
       y < big_to_small;
       y += big_to_small / steps) {
    Position<Big> const position(Big::origin +
                                     Displacement<Big>({0 * Kilo(Metre),
                                                        y,
                                                        0 * Kilo(Metre)}));
    Acceleration const big_on_position =
        -big_gravitational_parameter_ / (y * y);
    Acceleration const small_on_position =
        small_gravitational_parameter_ /
            ((big_to_small - y) * (big_to_small - y));
    Vector<Acceleration, Big> const expected_acceleration(
                  {0 * si::Unit<Acceleration>,
                   small_on_position + big_on_position - small_on_big,
                   0 * si::Unit<Acceleration>});
    EXPECT_THAT(AbsoluteError(
                    big_frame_->GeometricAcceleration(
                        t0_,
                        DegreesOfFreedom<Big>(position, Big::unmoving)),
                    expected_acceleration),
                Lt(1e-10 * si::Unit<Acceleration>));
  }
}

TEST_F(BodyCentredNonRotatingReferenceFrameTest, Serialization) {
  serialization::ReferenceFrame message;
  small_frame_->WriteToMessage(&message);

  EXPECT_TRUE(message.HasExtension(
      serialization::BodyCentredNonRotatingReferenceFrame::extension));
  auto const extension = message.GetExtension(
      serialization::BodyCentredNonRotatingReferenceFrame::extension);
  EXPECT_TRUE(extension.has_centre());
  EXPECT_EQ(1, extension.centre());

  auto const read_small_frame =
      RigidReferenceFrame<ICRS, Small>::ReadFromMessage(message,
                                                        ephemeris_.get());
  EXPECT_THAT(read_small_frame, Not(IsNull()));

  Instant const t = t0_ + period_;
  DegreesOfFreedom<Small> const point_dof =
      {Displacement<Small>({10 * Metre, 20 * Metre, 30 * Metre}) +
           Small::origin,
       Velocity<Small>({3 * Metre / Second,
                        2 * Metre / Second,
                        1 * Metre / Second})};
  EXPECT_EQ(small_frame_->GeometricAcceleration(t, point_dof),
            read_small_frame->GeometricAcceleration(t, point_dof));
}

}  // namespace physics
}  // namespace principia
