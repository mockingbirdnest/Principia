
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"

#include <memory>

#include "astronomy/frames.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
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
namespace physics {
namespace internal_body_centred_non_rotating_dynamic_frame {

using astronomy::ICRFJ2000Equator;
using geometry::Barycentre;
using geometry::Bivector;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::Frame;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::McLachlanAtela1992Order4Optimal;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::SIUnit;
using quantities::Time;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using ::testing::IsNull;
using ::testing::Lt;
using ::testing::Not;

namespace {

constexpr char big[] = "Big";
constexpr char small[] = "Small";

}  // namespace

class BodyCentredNonRotatingDynamicFrameTest : public ::testing::Test {
 protected:
  // The non-rotating frame centred on the big body.
  using Big = Frame<serialization::Frame::TestTag,
                    serialization::Frame::TEST, /*inertial=*/false>;
  using Small = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, /*inertial=*/false>;

  BodyCentredNonRotatingDynamicFrameTest()
      : period_(10 * π * sqrt(5.0 / 7.0) * Second),
        solar_system_(SOLUTION_DIR / "astronomy" /
                          "test_gravity_model_two_bodies.proto.txt",
                      SOLUTION_DIR / "astronomy" /
                          "test_initial_state_two_bodies_circular.proto.txt"),
        t0_(solar_system_.epoch()),
        ephemeris_(solar_system_.MakeEphemeris(
            /*fitting_tolerance=*/1 * Milli(Metre),
            Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
                SymplecticRungeKuttaNyströmIntegrator<
                    McLachlanAtela1992Order4Optimal,
                    Position<ICRFJ2000Equator>>(),
                /*step=*/10 * Milli(Second)))),
        big_initial_state_(solar_system_.degrees_of_freedom(big)),
        big_gravitational_parameter_(
            solar_system_.gravitational_parameter(big)),
        small_initial_state_(solar_system_.degrees_of_freedom(small)),
        small_gravitational_parameter_(
            solar_system_.gravitational_parameter(small)) {
    ephemeris_->Prolong(t0_ + 2 * period_);
    big_frame_ = std::make_unique<
                     BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Big>>(
                         ephemeris_.get(),
                         solar_system_.massive_body(*ephemeris_, big));
    small_frame_ = std::make_unique<
                     BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator,
                                                        Small>>(
                         ephemeris_.get(),
                         solar_system_.massive_body(*ephemeris_, small));
  }

  Time const period_;
  SolarSystem<ICRFJ2000Equator> solar_system_;
  Instant const t0_;
  std::unique_ptr<Ephemeris<ICRFJ2000Equator>> const ephemeris_;
  DegreesOfFreedom<ICRFJ2000Equator> const big_initial_state_;
  GravitationalParameter const big_gravitational_parameter_;
  DegreesOfFreedom<ICRFJ2000Equator> const small_initial_state_;
  GravitationalParameter const small_gravitational_parameter_;
  std::unique_ptr<
      BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Big>> big_frame_;
  std::unique_ptr<
      BodyCentredNonRotatingDynamicFrame<ICRFJ2000Equator, Small>> small_frame_;
};


// Check that the small body has the right degrees of freedom in the frame of
// the big body.
TEST_F(BodyCentredNonRotatingDynamicFrameTest, SmallBodyInBigFrame) {
  int const steps = 100;
  Bivector<double, ICRFJ2000Equator> const axis({0, 0, 1});

  RelativeDegreesOfFreedom<ICRFJ2000Equator> const initial_big_to_small =
      small_initial_state_ - big_initial_state_;
  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    DegreesOfFreedom<ICRFJ2000Equator> const small_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, small).
            EvaluateDegreesOfFreedom(t);

    auto const rotation_in_big_frame_at_t =
        Rotation<ICRFJ2000Equator, Big>(2 * π * (t - t0_) * Radian / period_,
                                        axis,
                                        DefinesFrame<Big>{});
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

TEST_F(BodyCentredNonRotatingDynamicFrameTest, Inverse) {
  int const steps = 100;
  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    auto const from_big_frame_at_t = big_frame_->FromThisFrameAtTime(t);
    auto const to_big_frame_at_t = big_frame_->ToThisFrameAtTime(t);
    auto const small_initial_state_transformed_and_back =
        from_big_frame_at_t(to_big_frame_at_t(small_initial_state_));
    EXPECT_THAT(small_initial_state_transformed_and_back.position(),
                AlmostEquals(small_initial_state_.position(), 0, 1));
    EXPECT_THAT(small_initial_state_transformed_and_back.velocity(),
                AlmostEquals(small_initial_state_.velocity(), 0, 1));
  }
}

TEST_F(BodyCentredNonRotatingDynamicFrameTest, GeometricAcceleration) {
  int const steps = 10;
  RelativeDegreesOfFreedom<ICRFJ2000Equator> const initial_big_to_small =
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
                  {0 * SIUnit<Acceleration>(),
                   small_on_position + big_on_position - small_on_big,
                   0 * SIUnit<Acceleration>()});
    EXPECT_THAT(AbsoluteError(
                    big_frame_->GeometricAcceleration(
                        t0_,
                        DegreesOfFreedom<Big>(position, Velocity<Big>())),
                    expected_acceleration),
                Lt(1e-10 * SIUnit<Acceleration>()));
  }
}

TEST_F(BodyCentredNonRotatingDynamicFrameTest, Serialization) {
  serialization::DynamicFrame message;
  small_frame_->WriteToMessage(&message);

  EXPECT_TRUE(message.HasExtension(
      serialization::BodyCentredNonRotatingDynamicFrame::extension));
  auto const extension = message.GetExtension(
      serialization::BodyCentredNonRotatingDynamicFrame::extension);
  EXPECT_TRUE(extension.has_centre());
  EXPECT_EQ(1, extension.centre());

  auto const read_small_frame =
      DynamicFrame<ICRFJ2000Equator, Small>::ReadFromMessage(
          message, ephemeris_.get());
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

}  // namespace internal_body_centred_non_rotating_dynamic_frame
}  // namespace physics
}  // namespace principia
