
#include "physics/barycentric_rotating_dynamic_frame.hpp"

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
#include "physics/mock_continuous_trajectory.hpp"
#include "physics/mock_ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace physics {
namespace internal_barycentric_rotating_dynamic_frame {

using astronomy::ICRS;
using base::check_not_null;
using geometry::Barycentre;
using geometry::Bivector;
using geometry::Frame;
using geometry::Instant;
using geometry::Rotation;
using geometry::Vector;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::McLachlanAtela1992Order4Optimal;
using quantities::Time;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using ::testing::Eq;
using ::testing::InSequence;
using ::testing::IsNull;
using ::testing::Lt;
using ::testing::Not;
using ::testing::Return;
using ::testing::StrictMock;
using ::testing::_;

namespace {

char constexpr big[] = "Big";
char constexpr small[] = "Small";

}  // namespace

class BarycentricRotatingDynamicFrameTest : public ::testing::Test {
 protected:
  // The rotating frame centred on the barycentre of the two bodies.
  using BigSmallFrame = Frame<serialization::Frame::TestTag,
                              serialization::Frame::TEST, /*inertial=*/false>;
  using MockFrame = Frame<serialization::Frame::TestTag,
                          serialization::Frame::TEST1, /*inertial=*/false>;

  BarycentricRotatingDynamicFrameTest()
      : period_(10 * π * sqrt(5.0 / 7.0) * Second),
        solar_system_(SOLUTION_DIR / "astronomy" /
                          "test_gravity_model_two_bodies.proto.txt",
                      SOLUTION_DIR / "astronomy" /
                          "test_initial_state_two_bodies_circular.proto.txt"),
        t0_(solar_system_.epoch()),
        ephemeris_(solar_system_.MakeEphemeris(
            /*fitting_tolerance=*/1 * Milli(Metre),
            Ephemeris<ICRS>::FixedStepParameters(
                SymplecticRungeKuttaNyströmIntegrator<
                    McLachlanAtela1992Order4Optimal,
                    Position<ICRS>>(),
                /*step=*/10 * Milli(Second)))),
        big_(solar_system_.massive_body(*ephemeris_, big)),
        big_initial_state_(solar_system_.degrees_of_freedom(big)),
        big_gravitational_parameter_(
            solar_system_.gravitational_parameter(big)),
        small_(solar_system_.massive_body(*ephemeris_, small)),
        small_initial_state_(solar_system_.degrees_of_freedom(small)),
        small_gravitational_parameter_(
            solar_system_.gravitational_parameter(small)),
        centre_of_mass_initial_state_(
            Barycentre<DegreesOfFreedom<ICRS>, GravitationalParameter>(
                {big_initial_state_, small_initial_state_},
                {big_gravitational_parameter_,
                 small_gravitational_parameter_})) {
    EXPECT_CALL(mock_ephemeris_,
                trajectory(solar_system_.massive_body(*ephemeris_, big)))
        .WillOnce(Return(&mock_big_trajectory_));
    EXPECT_CALL(mock_ephemeris_,
                trajectory(solar_system_.massive_body(*ephemeris_, small)))
        .WillOnce(Return(&mock_small_trajectory_));
    mock_frame_ =
        std::make_unique<BarycentricRotatingDynamicFrame<ICRS, MockFrame>>(
            &mock_ephemeris_, big_, small_);

    ephemeris_->Prolong(t0_ + 2 * period_);
    big_small_frame_ =
        std::make_unique<BarycentricRotatingDynamicFrame<ICRS, BigSmallFrame>>(
            ephemeris_.get(), big_, small_);
  }

  Time const period_;
  SolarSystem<ICRS> solar_system_;
  Instant const t0_;
  std::unique_ptr<Ephemeris<ICRS>> const ephemeris_;
  MassiveBody const* const big_;
  DegreesOfFreedom<ICRS> const big_initial_state_;
  GravitationalParameter const big_gravitational_parameter_;
  MassiveBody const* const small_;
  DegreesOfFreedom<ICRS> const small_initial_state_;
  GravitationalParameter const small_gravitational_parameter_;
  DegreesOfFreedom<ICRS> const centre_of_mass_initial_state_;
  StrictMock<MockEphemeris<ICRS>> mock_ephemeris_;

  std::unique_ptr<BarycentricRotatingDynamicFrame<ICRS, MockFrame>> mock_frame_;
  std::unique_ptr<BarycentricRotatingDynamicFrame<ICRS, BigSmallFrame>>
      big_small_frame_;
  StrictMock<MockContinuousTrajectory<ICRS>> mock_big_trajectory_;
  StrictMock<MockContinuousTrajectory<ICRS>> mock_small_trajectory_;
};


TEST_F(BarycentricRotatingDynamicFrameTest, ToBigSmallFrameAtTime) {
  int const steps = 100;

  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    auto const to_big_small_frame_at_t = big_small_frame_->ToThisFrameAtTime(t);

    // Check that the centre of mass is at the origin and doesn't move.
    DegreesOfFreedom<BigSmallFrame> const centre_of_mass_in_big_small_at_t =
        to_big_small_frame_at_t(centre_of_mass_initial_state_);
    EXPECT_THAT(AbsoluteError(centre_of_mass_in_big_small_at_t.position(),
                              BigSmallFrame::origin),
                Lt(1.0e-11 * Metre));
    EXPECT_THAT(AbsoluteError(centre_of_mass_in_big_small_at_t.velocity(),
                              Velocity<BigSmallFrame>()),
                Lt(1.1e-11 * Metre / Second));

    // Check that the bodies don't move and are at the right locations.
    DegreesOfFreedom<ICRS> const big_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, big).EvaluateDegreesOfFreedom(t);
    DegreesOfFreedom<ICRS> const small_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, small)
            .EvaluateDegreesOfFreedom(t);

    DegreesOfFreedom<BigSmallFrame> const big_in_big_small_at_t =
        to_big_small_frame_at_t(big_in_inertial_frame_at_t);
    DegreesOfFreedom<BigSmallFrame> const small_in_big_small_at_t =
        to_big_small_frame_at_t(small_in_inertial_frame_at_t);
    EXPECT_THAT(AbsoluteError(big_in_big_small_at_t.position(),
                              Displacement<BigSmallFrame>({
                                  -10.0 / 7.0 * Kilo(Metre),
                                  0 * Kilo(Metre),
                                  0 * Kilo(Metre)}) + BigSmallFrame::origin),
                Lt(1.0e-6 * Metre));
    EXPECT_THAT(AbsoluteError(big_in_big_small_at_t.velocity(),
                              Velocity<BigSmallFrame>()),
                Lt(1.0e-4 * Metre / Second));
    EXPECT_THAT(AbsoluteError(small_in_big_small_at_t.position(),
                              Displacement<BigSmallFrame>({
                                  25.0 / 7.0 * Kilo(Metre),
                                  0 * Kilo(Metre),
                                  0 * Kilo(Metre)}) + BigSmallFrame::origin),
                Lt(1.0e-5 * Metre));
    EXPECT_THAT(AbsoluteError(small_in_big_small_at_t.velocity(),
                              Velocity<BigSmallFrame>()),
                Lt(1.0e-4 * Metre / Second));
  }
}

TEST_F(BarycentricRotatingDynamicFrameTest, Inverse) {
  int const steps = 100;
  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    auto const from_big_small_frame_at_t =
        big_small_frame_->FromThisFrameAtTime(t);
    auto const to_big_small_frame_at_t = big_small_frame_->ToThisFrameAtTime(t);
    auto const small_initial_state_transformed_and_back =
        from_big_small_frame_at_t(to_big_small_frame_at_t(
            small_initial_state_));
    EXPECT_THAT(
        AbsoluteError(small_initial_state_transformed_and_back.position(),
                      small_initial_state_.position()),
        Lt(1.0e-11 * Metre));
    EXPECT_THAT(
        AbsoluteError(small_initial_state_transformed_and_back.velocity(),
                      small_initial_state_.velocity()),
        Lt(1.0e-11 * Metre / Second));
  }
}

// Two bodies in rotation with their barycentre at rest.  The test point is at
// the origin and in motion.  The acceleration is purely due to Coriolis.
TEST_F(BarycentricRotatingDynamicFrameTest, CoriolisAcceleration) {
  Instant const t = t0_ + 0 * Second;
  // The velocity is opposed to the motion and away from the centre.
  DegreesOfFreedom<MockFrame> const point_dof =
      {Displacement<MockFrame>({0 * Metre, 0 * Metre, 0 * Metre}) +
           MockFrame::origin,
       Velocity<MockFrame>({(80 - 30) * Metre / Second,
                            (-60 - 40) * Metre / Second,
                            0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const big_dof =
      {Displacement<ICRS>({0.8 * Metre, -0.6 * Metre, 0 * Metre}) +
       ICRS::origin,
       Velocity<ICRS>(
           {-16 * Metre / Second, 12 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const small_dof =
      {Displacement<ICRS>({5 * Metre, 5 * Metre, 0 * Metre}) + ICRS::origin,
       Velocity<ICRS>(
           {40 * Metre / Second, -30 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const barycentre_dof =
      Barycentre<DegreesOfFreedom<ICRS>, GravitationalParameter>(
          {big_dof, small_dof},
          {big_gravitational_parameter_, small_gravitational_parameter_});
  EXPECT_THAT(barycentre_dof.position() - ICRS::origin,
              Eq(Displacement<ICRS>({2 * Metre, 1 * Metre, 0 * Metre})));
  EXPECT_THAT(barycentre_dof.velocity(), Eq(Velocity<ICRS>()));

  EXPECT_CALL(mock_big_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(big_dof));
  EXPECT_CALL(mock_small_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(small_dof));
  {
    InSequence s;
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(big_), t))
        .WillOnce(Return(Vector<Acceleration, ICRS>({
                             120 * Metre / Pow<2>(Second),
                             160 * Metre / Pow<2>(Second),
                             0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(small_), t))
        .WillOnce(Return(Vector<Acceleration, ICRS>({
                             -300 * Metre / Pow<2>(Second),
                             -400 * Metre / Pow<2>(Second),
                             0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMasslessBody(_, t))
        .WillOnce(Return(Vector<Acceleration, ICRS>()));
  }

  // The Coriolis acceleration is towards the centre and opposed to the motion.
  EXPECT_THAT(mock_frame_->GeometricAcceleration(t, point_dof),
              AlmostEquals(Vector<Acceleration, MockFrame>({
                               (-1200 - 800) * Metre / Pow<2>(Second),
                               (-1600 + 600) * Metre / Pow<2>(Second),
                               0 * Metre / Pow<2>(Second)}), 0));
}

// Two bodies in rotation with their barycentre at rest.  The test point doesn't
// move so the acceleration is purely centrifugal.
TEST_F(BarycentricRotatingDynamicFrameTest, CentrifugalAcceleration) {
  Instant const t = t0_ + 0 * Second;
  DegreesOfFreedom<MockFrame> const point_dof =
      {Displacement<MockFrame>({10 * Metre, 20 * Metre, 30 * Metre}) +
           MockFrame::origin,
       Velocity<MockFrame>({0 * Metre / Second,
                            0 * Metre / Second,
                            0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const big_dof =
      {Displacement<ICRS>({0.8 * Metre, -0.6 * Metre, 0 * Metre}) +
       ICRS::origin,
       Velocity<ICRS>(
           {-16 * Metre / Second, 12 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const small_dof =
      {Displacement<ICRS>({5 * Metre, 5 * Metre, 0 * Metre}) + ICRS::origin,
       Velocity<ICRS>(
           {40 * Metre / Second, -30 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const barycentre_dof =
      Barycentre<DegreesOfFreedom<ICRS>, GravitationalParameter>(
          {big_dof, small_dof},
          {big_gravitational_parameter_, small_gravitational_parameter_});
  EXPECT_THAT(barycentre_dof.position() - ICRS::origin,
              Eq(Displacement<ICRS>({2 * Metre, 1 * Metre, 0 * Metre})));
  EXPECT_THAT(barycentre_dof.velocity(), Eq(Velocity<ICRS>()));

  EXPECT_CALL(mock_big_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(big_dof));
  EXPECT_CALL(mock_small_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(small_dof));
  {
    InSequence s;
    EXPECT_CALL(
        mock_ephemeris_,
        ComputeGravitationalAccelerationOnMassiveBody(check_not_null(big_), t))
        .WillOnce(
            Return(Vector<Acceleration, ICRS>({120 * Metre / Pow<2>(Second),
                                               160 * Metre / Pow<2>(Second),
                                               0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(small_), t))
        .WillOnce(
            Return(Vector<Acceleration, ICRS>({-300 * Metre / Pow<2>(Second),
                                               -400 * Metre / Pow<2>(Second),
                                               0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMasslessBody(_, t))
        .WillOnce(Return(Vector<Acceleration, ICRS>()));
  }

  EXPECT_THAT(mock_frame_->GeometricAcceleration(t, point_dof),
              AlmostEquals(Vector<Acceleration, MockFrame>({
                               1e3 * Metre / Pow<2>(Second),
                               2e3 * Metre / Pow<2>(Second),
                               0 * Metre / Pow<2>(Second)}), 2));
}

// Two bodies in rotation with their barycentre at rest, with a tangential
// acceleration that increases their rotational speed.  The test point doesn't
// move.  The resulting acceleration combines centrifugal and Euler.
TEST_F(BarycentricRotatingDynamicFrameTest, EulerAcceleration) {
  Instant const t = t0_ + 0 * Second;
  DegreesOfFreedom<MockFrame> const point_dof =
      {Displacement<MockFrame>({10 * Metre, 20 * Metre, 30 * Metre}) +
           MockFrame::origin,
       Velocity<MockFrame>({0 * Metre / Second,
                            0 * Metre / Second,
                            0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const big_dof =
      {Displacement<ICRS>({0.8 * Metre, -0.6 * Metre, 0 * Metre}) +
       ICRS::origin,
       Velocity<ICRS>(
           {-16 * Metre / Second, 12 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const small_dof =
      {Displacement<ICRS>({5 * Metre, 5 * Metre, 0 * Metre}) + ICRS::origin,
       Velocity<ICRS>(
           {40 * Metre / Second, -30 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const barycentre_dof =
      Barycentre<DegreesOfFreedom<ICRS>, GravitationalParameter>(
          {big_dof, small_dof},
          {big_gravitational_parameter_, small_gravitational_parameter_});
  EXPECT_THAT(barycentre_dof.position() - ICRS::origin,
              Eq(Displacement<ICRS>({2 * Metre, 1 * Metre, 0 * Metre})));
  EXPECT_THAT(barycentre_dof.velocity(), Eq(Velocity<ICRS>()));

  EXPECT_CALL(mock_big_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(big_dof));
  EXPECT_CALL(mock_small_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(small_dof));
  {
    // The acceleration is centripetal + tangential.
    InSequence s;
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(big_), t))
        .WillOnce(Return(Vector<Acceleration, ICRS>({
                             (120 - 160) * Metre / Pow<2>(Second),
                             (160 + 120) * Metre / Pow<2>(Second),
                             0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(small_), t))
        .WillOnce(Return(
            Vector<Acceleration, ICRS>({(-300 + 400) * Metre / Pow<2>(Second),
                                        (-400 - 300) * Metre / Pow<2>(Second),
                                        0 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMasslessBody(_, t))
        .WillOnce(Return(Vector<Acceleration, ICRS>()));
  }

  // The acceleration is centrifugal + Euler.
  EXPECT_THAT(mock_frame_->GeometricAcceleration(t, point_dof),
              AlmostEquals(Vector<Acceleration, MockFrame>({
                               (1e3 + 2e3) * Metre / Pow<2>(Second),
                               (2e3 - 1e3) * Metre / Pow<2>(Second),
                               0 * Metre / Pow<2>(Second)}), 1));
}

// Two bodies in rotation with their barycentre at rest, with a linear
// acceleration identical for both bodies.  The test point doesn't move.  The
// resulting acceleration combines centrifugal and linear.
TEST_F(BarycentricRotatingDynamicFrameTest, LinearAcceleration) {
  Instant const t = t0_ + 0 * Second;
  DegreesOfFreedom<MockFrame> const point_dof =
      {Displacement<MockFrame>({10 * Metre, 20 * Metre, 30 * Metre}) +
           MockFrame::origin,
       Velocity<MockFrame>({0 * Metre / Second,
                            0 * Metre / Second,
                            0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const big_dof =
      {Displacement<ICRS>({0.8 * Metre, -0.6 * Metre, 0 * Metre}) +
       ICRS::origin,
       Velocity<ICRS>(
           {-16 * Metre / Second, 12 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const small_dof =
      {Displacement<ICRS>({5 * Metre, 5 * Metre, 0 * Metre}) + ICRS::origin,
       Velocity<ICRS>(
           {40 * Metre / Second, -30 * Metre / Second, 0 * Metre / Second})};
  DegreesOfFreedom<ICRS> const barycentre_dof =
      Barycentre<DegreesOfFreedom<ICRS>, GravitationalParameter>(
          {big_dof, small_dof},
          {big_gravitational_parameter_, small_gravitational_parameter_});
  EXPECT_THAT(barycentre_dof.position() - ICRS::origin,
              Eq(Displacement<ICRS>({2 * Metre, 1 * Metre, 0 * Metre})));
  EXPECT_THAT(barycentre_dof.velocity(), Eq(Velocity<ICRS>()));

  EXPECT_CALL(mock_big_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(big_dof));
  EXPECT_CALL(mock_small_trajectory_, EvaluateDegreesOfFreedom(t))
      .Times(2)
      .WillRepeatedly(Return(small_dof));
  {
    // The acceleration is linear + centripetal.
    InSequence s;
    EXPECT_CALL(
        mock_ephemeris_,
        ComputeGravitationalAccelerationOnMassiveBody(check_not_null(big_), t))
        .WillOnce(Return(
            Vector<Acceleration, ICRS>({(-160 + 120) * Metre / Pow<2>(Second),
                                        (120 + 160) * Metre / Pow<2>(Second),
                                        300 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMassiveBody(
                    check_not_null(small_), t))
        .WillOnce(Return(
            Vector<Acceleration, ICRS>({(-160 - 300) * Metre / Pow<2>(Second),
                                        (120 - 400) * Metre / Pow<2>(Second),
                                        300 * Metre / Pow<2>(Second)})));
    EXPECT_CALL(mock_ephemeris_,
                ComputeGravitationalAccelerationOnMasslessBody(_, t))
        .WillOnce(Return(Vector<Acceleration, ICRS>()));
  }

  // The acceleration is linear + centrifugal.
  EXPECT_THAT(mock_frame_->GeometricAcceleration(t, point_dof),
              AlmostEquals(Vector<Acceleration, MockFrame>({
                               1e3 * Metre / Pow<2>(Second),
                               (200 + 2e3) * Metre / Pow<2>(Second),
                               300 * Metre / Pow<2>(Second)}), 2));
}

TEST_F(BarycentricRotatingDynamicFrameTest, GeometricAcceleration) {
  Instant const t = t0_ + period_;
  DegreesOfFreedom<BigSmallFrame> const point_dof =
      {Displacement<BigSmallFrame>({10 * Metre, 20 * Metre, 30 * Metre}) +
           BigSmallFrame::origin,
       Velocity<BigSmallFrame>({3 * Metre / Second,
                                2 * Metre / Second,
                                1 * Metre / Second})};
  // We trust the functions to compute the values correctly, but this test
  // ensures that we don't get NaNs.
  EXPECT_THAT(big_small_frame_->GeometricAcceleration(t, point_dof),
              AlmostEquals(Vector<Acceleration, BigSmallFrame>({
                  -2.25461792868978819e3 * Metre / Pow<2>(Second),
                  -3.47432622325080658e1 * Metre / Pow<2>(Second),
                  -5.16651053897896801e1 * Metre / Pow<2>(Second)}), 1));
}

TEST_F(BarycentricRotatingDynamicFrameTest, Serialization) {
  serialization::DynamicFrame message;
  big_small_frame_->WriteToMessage(&message);

  EXPECT_TRUE(message.HasExtension(
      serialization::BarycentricRotatingDynamicFrame::extension));
  auto const extension = message.GetExtension(
      serialization::BarycentricRotatingDynamicFrame::extension);
  EXPECT_TRUE(extension.has_primary());
  EXPECT_TRUE(extension.has_secondary());
  EXPECT_EQ(0, extension.primary());
  EXPECT_EQ(1, extension.secondary());

  auto const read_big_small_frame =
      DynamicFrame<ICRS, BigSmallFrame>::ReadFromMessage(message,
                                                         ephemeris_.get());
  EXPECT_THAT(read_big_small_frame, Not(IsNull()));

  Instant const t = t0_ + period_;
  DegreesOfFreedom<BigSmallFrame> const point_dof =
      {Displacement<BigSmallFrame>({10 * Metre, 20 * Metre, 30 * Metre}) +
           BigSmallFrame::origin,
       Velocity<BigSmallFrame>({3 * Metre / Second,
                                2 * Metre / Second,
                                1 * Metre / Second})};
  EXPECT_EQ(big_small_frame_->GeometricAcceleration(t, point_dof),
            read_big_small_frame->GeometricAcceleration(t, point_dof));
}

}  // namespace internal_barycentric_rotating_dynamic_frame
}  // namespace physics
}  // namespace principia
