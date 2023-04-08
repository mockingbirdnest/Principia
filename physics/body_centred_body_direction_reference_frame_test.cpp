#include "physics/body_centred_body_direction_reference_frame.hpp"

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
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using ::testing::IsNull;
using ::testing::Lt;
using ::testing::Not;
using ::testing::Return;
using ::testing::_;
using namespace principia::astronomy::_frames;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symplectic_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_vanishes_before;

namespace {

char constexpr big[] = "Big";
char constexpr small[] = "Small";

}  // namespace

class BodyCentredBodyDirectionReferenceFrameTest : public ::testing::Test {
 protected:
  // The rotating frame centred on the big body and directed to the small one.
  using BigSmallFrame = Frame<serialization::Frame::TestTag,
                              Arbitrary,
                              Handedness::Right,
                              serialization::Frame::TEST>;

  BodyCentredBodyDirectionReferenceFrameTest()
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
        big_(solar_system_.massive_body(*ephemeris_, big)),
        big_initial_state_(solar_system_.degrees_of_freedom(big)),
        big_gravitational_parameter_(
            solar_system_.gravitational_parameter(big)),
        small_(solar_system_.massive_body(*ephemeris_, small)),
        small_initial_state_(solar_system_.degrees_of_freedom(small)),
        small_gravitational_parameter_(
            solar_system_.gravitational_parameter(small)) {
    EXPECT_OK(ephemeris_->Prolong(t0_ + 2 * period_));
    big_small_frame_ = std::make_unique<
        BodyCentredBodyDirectionReferenceFrame<ICRS, BigSmallFrame>>(
        ephemeris_.get(), big_, small_);
  }

  Time const period_;
  SolarSystem<ICRS> solar_system_;
  Instant const t0_;
  std::unique_ptr<Ephemeris<ICRS>> const ephemeris_;
  MassiveBody const* big_;
  DegreesOfFreedom<ICRS> big_initial_state_;
  GravitationalParameter big_gravitational_parameter_;
  MassiveBody const* small_;
  DegreesOfFreedom<ICRS> small_initial_state_;
  GravitationalParameter small_gravitational_parameter_;

  std::unique_ptr<BodyCentredBodyDirectionReferenceFrame<ICRS, BigSmallFrame>>
      big_small_frame_;
};


TEST_F(BodyCentredBodyDirectionReferenceFrameTest, ToBigSmallFrameAtTime) {
  int const steps = 100;

  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    auto const to_big_small_frame_at_t = big_small_frame_->ToThisFrameAtTime(t);

    // Check that the bodies don't move and are at the right locations.
    DegreesOfFreedom<ICRS> const big_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, big).
            EvaluateDegreesOfFreedom(t);
    DegreesOfFreedom<ICRS> const small_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, small).
            EvaluateDegreesOfFreedom(t);

    DegreesOfFreedom<BigSmallFrame> const big_in_big_small_at_t =
        to_big_small_frame_at_t(big_in_inertial_frame_at_t);
    DegreesOfFreedom<BigSmallFrame> const small_in_big_small_at_t =
        to_big_small_frame_at_t(small_in_inertial_frame_at_t);
    EXPECT_THAT(AbsoluteError(big_in_big_small_at_t.position(),
                              BigSmallFrame::origin),
                Lt(1.0e-6 * Metre));
    EXPECT_THAT(AbsoluteError(big_in_big_small_at_t.velocity(),
                              BigSmallFrame::unmoving),
                Lt(1.0e-4 * Metre / Second));
    EXPECT_THAT(AbsoluteError(small_in_big_small_at_t.position(),
                              Displacement<BigSmallFrame>({
                                  5.0 * Kilo(Metre),
                                  0 * Kilo(Metre),
                                  0 * Kilo(Metre)}) + BigSmallFrame::origin),
                Lt(1.0e-5 * Metre));
    EXPECT_THAT(AbsoluteError(small_in_big_small_at_t.velocity(),
                              BigSmallFrame::unmoving),
                Lt(1.0e-4 * Metre / Second));
  }
}

TEST_F(BodyCentredBodyDirectionReferenceFrameTest, Inverse) {
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

TEST_F(BodyCentredBodyDirectionReferenceFrameTest, GeometricAcceleration) {
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
                  -9.54502614154908457e5 * Metre / Pow<2>(Second),
                  -1.90900949256416853e6 * Metre / Pow<2>(Second),
                  -2.86351378905829182e6 * Metre / Pow<2>(Second)}), 0));
}

TEST_F(BodyCentredBodyDirectionReferenceFrameTest, Serialization) {
  serialization::RigidReferenceFrame message;
  big_small_frame_->WriteToMessage(&message);

  EXPECT_TRUE(message.HasExtension(
      serialization::BodyCentredBodyDirectionReferenceFrame::extension));
  auto const extension = message.GetExtension(
      serialization::BodyCentredBodyDirectionReferenceFrame::extension);
  EXPECT_TRUE(extension.has_primary());
  EXPECT_TRUE(extension.has_secondary());
  EXPECT_EQ(0, extension.primary());
  EXPECT_EQ(1, extension.secondary());

  auto const read_big_small_frame =
      RigidReferenceFrame<ICRS, BigSmallFrame>::ReadFromMessage(
          message, ephemeris_.get());
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

TEST_F(BodyCentredBodyDirectionReferenceFrameTest, ConstructFromOneBody) {
  // A discrete trajectory that remains motionless at the barycentre.  Since
  // both bodies don't have the same mass, this means it has an intrinsic
  // acceleration.
  DiscreteTrajectory<ICRS> barycentre_trajectory;
  for (Time t; t <= period_; t += period_ / 16) {
    auto const big_dof =
        ephemeris_->trajectory(big_)->EvaluateDegreesOfFreedom(t0_ + t);
    auto const small_dof =
        ephemeris_->trajectory(small_)->EvaluateDegreesOfFreedom(t0_ + t);
    auto const barycentre =
        Barycentre<DegreesOfFreedom<ICRS>, GravitationalParameter>(
            {big_dof, small_dof},
            {big_->gravitational_parameter(),
             small_->gravitational_parameter()});
    EXPECT_THAT(barycentre.velocity().Norm(),
                VanishesBefore(1 * Kilo(Metre) / Second, 0, 50));
    EXPECT_OK(barycentre_trajectory.Append(t0_ + t, barycentre));
  }
  BodyCentredBodyDirectionReferenceFrame<ICRS, BigSmallFrame>
      barycentric_from_discrete{
          ephemeris_.get(),
          [&t = barycentre_trajectory]() -> auto& { return t; },
          small_};
  BarycentricRotatingReferenceFrame<ICRS, BigSmallFrame>
      barycentric_from_both_bodies{ephemeris_.get(), big_, small_};
  for (Time t = period_ / 32; t <= period_ / 2; t += period_ / 32) {
    auto const dof_from_discrete =
        barycentric_from_discrete.ToThisFrameAtTime(t0_ + t)(
            {ICRS::origin, ICRS::unmoving});
    auto const dof_from_both_bodies =
        barycentric_from_both_bodies.ToThisFrameAtTime(t0_ + t)(
            {ICRS::origin, ICRS::unmoving});
    EXPECT_THAT(
        (dof_from_discrete.position() - dof_from_both_bodies.position()).Norm(),
        VanishesBefore(1 * Kilo(Metre), 0, 15));
    EXPECT_THAT(
        (dof_from_discrete.velocity() - dof_from_both_bodies.velocity()).Norm(),
        VanishesBefore(1 * Kilo(Metre) / Second, 0, 93));
    // For the moment, the |BodyCentredBodyDirectionReferenceFrame| assumes that
    // its reference trajectories are free-falling, and gives us the wrong
    // geometric acceleration when this is not the case.
    auto const intrinsic_acceleration =
        ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(
            ICRS::origin + Displacement<ICRS>({0 * Kilo(Metre),
                                               10.0 / 7.0 * Kilo(Metre),
                                               0 * Kilo(Metre)}),
            t0_ + t);
    EXPECT_THAT(
        (barycentric_from_discrete.GeometricAcceleration(t0_ + t,
                                                         dof_from_discrete) -
         barycentric_from_both_bodies.GeometricAcceleration(
             t0_ + t,
             dof_from_both_bodies)).Norm(),
         AlmostEquals(intrinsic_acceleration.Norm(), 0, 142));
  }
}

}  // namespace physics
}  // namespace principia
