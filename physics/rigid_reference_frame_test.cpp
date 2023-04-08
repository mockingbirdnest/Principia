#include "physics/rigid_reference_frame.hpp"

#include "geometry/frame.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/mock_rigid_reference_frame.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {

using ::testing::Gt;
using ::testing::Invoke;
using ::testing::Lt;
using ::testing::Return;
using ::testing::StrictMock;
using ::testing::_;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;
using namespace principia::testing_utilities::_vanishes_before;

class ReferenceFrameTest : public testing::Test {
 protected:
  using Inertial = Frame<struct InertialFrameTag, geometry::Inertial>;
  using Rotating = Frame<struct RotatingFrameTag, Arbitrary>;
  using Translating = Rotating;  // A better name for linear acceleration.

  // Computes a first-order approximation of the acceleration using the
  // potential returned by |mock_frame_|.  Useful for checking that the
  // potential and the accelaration are consistent.
  Vector<Acceleration, Rotating> FirstOrderAccelerationFromPotential(
      Instant const& t,
      Position<Rotating> const& position,
      Length const& Δl) {
    auto const potential = mock_frame_.GeometricPotential(t, position);
    auto const potential_Δx = mock_frame_.GeometricPotential(
        t, position + Displacement<Rotating>({Δl, 0 * Metre, 0 * Metre}));
    auto const potential_Δy = mock_frame_.GeometricPotential(
        t, position + Displacement<Rotating>({0 * Metre, Δl, 0 * Metre}));
    auto const potential_Δz = mock_frame_.GeometricPotential(
        t, position + Displacement<Rotating>({0 * Metre, 0 * Metre, Δl}));

    return -Vector<Acceleration, Rotating>({(potential_Δx - potential) / Δl,
                                            (potential_Δy - potential) / Δl,
                                            (potential_Δz - potential) / Δl});
  }

  StrictMock<MockRigidReferenceFrame<Inertial, Rotating>> mock_frame_;
  Instant const t0_;

  // General values used to check that the acceleration does not depend on some
  // factors.
  Position<Inertial> const general_position_ =
      Inertial::origin +
      Displacement<Inertial>({11 * Metre, -13 * Metre, 17 * Metre});
  Displacement<Rotating> const general_displacement_ =
      Displacement<Rotating>({-23 * Metre, 29 * Metre, -31 * Metre});
  Displacement<Rotating> const general_displacement_z_ =
      Displacement<Rotating>({0 * Metre, 0 * Metre, 37 * Metre});
  Velocity<Rotating> const general_velocity_ = Velocity<Rotating>(
      {-41 * Metre / Second, 43 * Metre / Second, -47 * Metre / Second});
  Velocity<Rotating> const general_velocity_z_ = Velocity<Rotating>(
      {0 * Metre / Second, 0 * Metre / Second, 53 * Metre / Second});
};

// A frame uniformly accelerated along the z axis.  The test point is in general
// motion.  The acceleration is purely linear.
TEST_F(ReferenceFrameTest, LinearAcceleration) {
  EXPECT_CALL(mock_frame_, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([this](Instant const& t) {
        Acceleration const γ = 10 * Metre / Second / Second;
        auto const acceleration_of_to_frame_origin =
            Vector<Acceleration, Inertial>({0 * Metre / Second / Second,
                                            0 * Metre / Second / Second,
                                            γ});
        auto const velocity_of_to_frame_origin =
            Velocity<Inertial>({0 * Metre / Second,
                                0 * Metre / Second,
                                γ * (t - t0_)});
        auto const position_of_to_frame_origin =
            general_position_ +
            Displacement<Inertial>({0 * Metre,
                                    0 * Metre,
                                    γ * Pow<2>(t - t0_) / 2});
        RigidTransformation<Inertial, Translating> const
            rigid_transformation(
                /*from_origin=*/position_of_to_frame_origin,
                /*to_origin=*/Translating::origin,
                OrthogonalMap<Inertial, Translating>::Identity());
        AngularVelocity<Inertial> const angular_velocity_of_to_frame;
        RigidMotion<Inertial, Translating> const rigid_motion(
            rigid_transformation,
            angular_velocity_of_to_frame,
            velocity_of_to_frame_origin);
        Bivector<AngularAcceleration, Inertial> const
            angular_acceleration_of_to_frame;
        return AcceleratedRigidMotion<Inertial, Translating>(
            rigid_motion,
            angular_acceleration_of_to_frame,
            acceleration_of_to_frame_origin);
      }));

  // No gravity.
  Vector<Acceleration, Inertial> const gravitational_acceleration;
  SpecificEnergy const gravitational_potential;
  EXPECT_CALL(mock_frame_, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(gravitational_acceleration));
  EXPECT_CALL(mock_frame_, GravitationalPotential(_, _))
      .WillRepeatedly(Return(gravitational_potential));

  // The test point is in general position and velocity.
  DegreesOfFreedom<Translating> const initial_state_in_translating_frame = {
      Translating::origin + general_displacement_,
      general_velocity_};
  DegreesOfFreedom<Inertial> const initial_state_in_inertial_frame =
      mock_frame_.MotionOfThisFrame(t0_).rigid_motion().Inverse()(
          initial_state_in_translating_frame);

  // The time interval for evaluating the first order effect.
  Time const Δt = 1 * Milli(Second);

  Position<Inertial> const final_position_in_inertial_frame =
      initial_state_in_inertial_frame.position() +
      initial_state_in_inertial_frame.velocity() * Δt;

  Position<Translating> const final_position_in_translating_frame =
      mock_frame_.MotionOfThisFrame(t0_ + Δt)
          .rigid_motion()
          .rigid_transformation()(final_position_in_inertial_frame);

  Position<Translating> const first_order_final_position_in_translating_frame =
      initial_state_in_translating_frame.position() +
      initial_state_in_translating_frame.velocity() * Δt;

  Displacement<Translating> const higher_order_effect =
      final_position_in_translating_frame -
      first_order_final_position_in_translating_frame;

  // The second order effect is the linear acceleration.  This computation only
  // depends on the stub motion defined above.
  EXPECT_THAT(higher_order_effect,
              Componentwise(AlmostEquals(0 * Metre, 0),
                            AlmostEquals(0 * Metre, 0),
                            IsNear(-5.0_(1) * Micro(Metre))));

  // The linear acceleration matches that computed based on the motion to the
  // second order.
  EXPECT_THAT(mock_frame_.GeometricAcceleration(
                  t0_, initial_state_in_translating_frame) *
                  Pow<2>(Δt) / 2,
              AlmostEquals(Displacement<Translating>({0 * Metre,
                                                      0 * Metre,
                                                      -5 * Micro(Metre)}),
                           0));

  // The linear acceleration shows up for a point at rest.
  EXPECT_THAT(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                  t0_, initial_state_in_translating_frame.position()) *
                  Pow<2>(Δt) / 2,
              AlmostEquals(Displacement<Translating>({0 * Metre,
                                                      0 * Metre,
                                                      -5 * Micro(Metre)}),
                           0));

  // The linear acceleration derives from a potential.
  EXPECT_THAT(
      FirstOrderAccelerationFromPotential(
          t0_,
          initial_state_in_translating_frame.position(),
          /*Δl=*/1 * Micro(Metre)),
      RelativeErrorFrom(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                            t0_, initial_state_in_translating_frame.position()),
                        IsNear(3.1e-9_(1))));
}

// A frame in uniform rotation around the origin.  The test point is at the
// origin and in motion along the x axis.  The acceleration is purely due to
// Coriolis.
TEST_F(ReferenceFrameTest, CoriolisAcceleration) {
  EXPECT_CALL(mock_frame_, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([this](Instant const& t) {
        AngularFrequency const ω = 10 * Radian / Second;
        AngularVelocity<Inertial> const angular_velocity_of_to_frame(
            {0 * Radian / Second, 0 * Radian / Second, ω});
        Rotation<Inertial, Rotating> const rotation(
            ω * (t - t0_),
            angular_velocity_of_to_frame,
            DefinesFrame<Rotating>{});
        RigidTransformation<Inertial, Rotating> const
            rigid_transformation(
                /*from_origin=*/general_position_,
                /*to_origin=*/Rotating::origin,
                rotation.Forget<OrthogonalMap>());
        Velocity<Inertial> const velocity_of_to_frame_origin;
        RigidMotion<Inertial, Rotating> const rigid_motion(
            rigid_transformation,
            angular_velocity_of_to_frame,
            velocity_of_to_frame_origin);
        Bivector<AngularAcceleration, Inertial> const
            angular_acceleration_of_to_frame;
        Vector<Acceleration, Inertial> const
            acceleration_of_to_frame_origin;
        return AcceleratedRigidMotion<Inertial, Rotating>(
            rigid_motion,
            angular_acceleration_of_to_frame,
            acceleration_of_to_frame_origin);
      }));

  // No gravity.
  Vector<Acceleration, Inertial> const gravitational_acceleration;
  SpecificEnergy const gravitational_potential;
  EXPECT_CALL(mock_frame_, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(gravitational_acceleration));
  EXPECT_CALL(mock_frame_, GravitationalPotential(_, _))
      .WillRepeatedly(Return(gravitational_potential));

  // The velocity is along the x axis.
  DegreesOfFreedom<Rotating> const initial_state_in_rotating_frame =
      {Rotating::origin,
       general_velocity_z_ + Velocity<Rotating>({100 * Metre / Second,
                                                 0 * Metre / Second,
                                                 0 * Metre / Second})};
  DegreesOfFreedom<Inertial> const initial_state_in_inertial_frame =
      mock_frame_.MotionOfThisFrame(t0_).rigid_motion().Inverse()(
          initial_state_in_rotating_frame);

  // The time interval for evaluating the first order effect.
  Time const Δt = 1 * Milli(Second);

  Position<Inertial> const final_position_in_inertial_frame =
      initial_state_in_inertial_frame.position() +
      initial_state_in_inertial_frame.velocity() * Δt;

  Position<Rotating> const final_position_in_rotating_frame =
      mock_frame_.MotionOfThisFrame(t0_ + Δt)
          .rigid_motion()
          .rigid_transformation()(final_position_in_inertial_frame);

  Position<Rotating> const first_order_final_position_in_rotating_frame =
      initial_state_in_rotating_frame.position() +
      initial_state_in_rotating_frame.velocity() * Δt;

  Displacement<Rotating> const higher_order_effect =
      final_position_in_rotating_frame -
      first_order_final_position_in_rotating_frame;

  // The second order effect is the Coriolis acceleration, the higher order
  // effects are irrelevant.  This computation only depends on the stub motion
  // defined above.
  EXPECT_THAT(higher_order_effect,
              Componentwise(IsNear(-5.0_(1) * Micro(Metre)),
                            IsNear(-1.0_(1) * Milli(Metre)),
                            VanishesBefore(1 * Metre, 4)));

  // The Coriolis acceleration matches that computed based on the motion to the
  // second order.
  EXPECT_THAT(
      mock_frame_.GeometricAcceleration(t0_, initial_state_in_rotating_frame) *
          Pow<2>(Δt) / 2,
      AlmostEquals(Displacement<Rotating>({0 * Metre,
                                           -1 * Milli(Metre),
                                           0 * Metre}),
                   0));

  // No Coriolis acceleration when at rest.
  EXPECT_THAT(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                  t0_, initial_state_in_rotating_frame.position()),
              AlmostEquals(Vector<Acceleration, Rotating>(), 0));

  // No (scalar) Coriolis potential.  There is a bit of noise in x and y because
  // moving the points for the first-order computation introduces a tiny
  // centrifugal force.
  EXPECT_THAT(FirstOrderAccelerationFromPotential(
                  t0_,
                  initial_state_in_rotating_frame.position(),
                  /*Δl=*/1 * Micro(Metre)),
              Componentwise(IsNear(50.0_(1) * Micro(Metre) / Second / Second),
                            IsNear(50.0_(1) * Micro(Metre) / Second / Second),
                            AlmostEquals(0 * Metre / Second / Second, 0)));
}

// A frame in uniform rotation around the origin.  The test point is on the x
// axis.  The acceleration is purely due to centrifugal effects.
TEST_F(ReferenceFrameTest, CentrifugalAcceleration) {
  EXPECT_CALL(mock_frame_, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([this](Instant const& t) {
        AngularFrequency const ω = 10 * Radian / Second;
        AngularVelocity<Inertial> const angular_velocity_of_to_frame(
            {0 * Radian / Second, 0 * Radian / Second, ω});
        Rotation<Inertial, Rotating> const rotation(
            ω * (t - t0_),
            angular_velocity_of_to_frame,
            DefinesFrame<Rotating>{});
        RigidTransformation<Inertial, Rotating> const
            rigid_transformation(
                /*from_origin=*/general_position_,
                /*to_origin=*/Rotating::origin,
                rotation.Forget<OrthogonalMap>());
        Velocity<Inertial> const velocity_of_to_frame_origin;
        RigidMotion<Inertial, Rotating> const rigid_motion(
            rigid_transformation,
            angular_velocity_of_to_frame,
            velocity_of_to_frame_origin);
        Bivector<AngularAcceleration, Inertial> const
            angular_acceleration_of_to_frame;
        Vector<Acceleration, Inertial> const
            acceleration_of_to_frame_origin;
        return AcceleratedRigidMotion<Inertial, Rotating>(
            rigid_motion,
            angular_acceleration_of_to_frame,
            acceleration_of_to_frame_origin);
      }));

  // No gravity.
  Vector<Acceleration, Inertial> const gravitational_acceleration;
  SpecificEnergy const gravitational_potential;
  EXPECT_CALL(mock_frame_, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(gravitational_acceleration));
  EXPECT_CALL(mock_frame_, GravitationalPotential(_, _))
      .WillRepeatedly(Return(gravitational_potential));

  // The test point is on the x axis.
  DegreesOfFreedom<Rotating> const initial_state_in_rotating_frame = {
      Rotating::origin +
          general_displacement_z_ + Displacement<Rotating>({100 * Metre,
                                                            0 * Metre,
                                                            0 * Metre}),
      general_velocity_z_};
  DegreesOfFreedom<Inertial> const initial_state_in_inertial_frame =
      mock_frame_.MotionOfThisFrame(t0_).rigid_motion().Inverse()(
          initial_state_in_rotating_frame);

  // The time interval for evaluating the first order effect.
  Time const Δt = 1 * Milli(Second);

  Position<Inertial> const final_position_in_inertial_frame =
      initial_state_in_inertial_frame.position() +
      initial_state_in_inertial_frame.velocity() * Δt;

  Position<Rotating> const final_position_in_rotating_frame =
      mock_frame_.MotionOfThisFrame(t0_ + Δt)
          .rigid_motion()
          .rigid_transformation()(final_position_in_inertial_frame);

  Position<Rotating> const first_order_final_position_in_rotating_frame =
      initial_state_in_rotating_frame.position() +
      initial_state_in_rotating_frame.velocity() * Δt;

  Displacement<Rotating> const higher_order_effect =
      final_position_in_rotating_frame -
      first_order_final_position_in_rotating_frame;

  // The second order effect is the centrifugal acceleration, the higher order
  // effects are irrelevant.  This computation only depends on the stub motion
  // defined above.
  EXPECT_THAT(higher_order_effect,
              Componentwise(IsNear(5.0_(1) * Milli(Metre)),
                            IsNear(-33.3_(1) * Micro(Metre)),
                            AlmostEquals(0 * Metre, 0)));

  // The centrifugal acceleration matches that computed based on the motion to
  // the second order.
  EXPECT_THAT(
      mock_frame_.GeometricAcceleration(t0_, initial_state_in_rotating_frame) *
          Pow<2>(Δt) / 2,
      AlmostEquals(Displacement<Rotating>({5 * Milli(Metre),
                                           0 * Metre,
                                           0 * Metre}),
                   0));

  // The centrifugal acceleration shows up for a point at rest.
  EXPECT_THAT(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                  t0_, initial_state_in_rotating_frame.position()) *
                  Pow<2>(Δt) / 2,
              AlmostEquals(Displacement<Rotating>({5 * Milli(Metre),
                                                   0 * Metre,
                                                   0 * Metre}),
                           0));

  // The centrifugal acceleration derives from a potential.
  EXPECT_THAT(
      FirstOrderAccelerationFromPotential(
          t0_,
          initial_state_in_rotating_frame.position(),
          /*Δl=*/1 * Micro(Metre)),
      RelativeErrorFrom(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                            t0_, initial_state_in_rotating_frame.position()),
                        IsNear(5.9e-9_(1))));
}

// A frame initially nonrotating and in uniformly accelerated rotation around
// the origin.  The test point is on the x axis.  The acceleration is purely due
// to Euler.
TEST_F(ReferenceFrameTest, EulerAcceleration) {
  EXPECT_CALL(mock_frame_, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([this](Instant const& t) {
        AngularAcceleration const ωʹ = 10 * Radian / Second / Second;
        Bivector<AngularAcceleration, Inertial> const
            angular_acceleration_of_to_frame({0 * Radian / Second / Second,
                                              0 * Radian / Second / Second,
                                              ωʹ});
        Rotation<Inertial, Rotating> const rotation(
            ωʹ * Pow<2>(t - t0_) / 2,
            angular_acceleration_of_to_frame,
            DefinesFrame<Rotating>{});
        RigidTransformation<Inertial, Rotating> const
            rigid_transformation(
                /*from_origin=*/general_position_,
                /*to_origin=*/Rotating::origin,
                rotation.Forget<OrthogonalMap>());
        AngularVelocity<Inertial> const angular_velocity_of_to_frame;
        Velocity<Inertial> const velocity_of_to_frame_origin;
        RigidMotion<Inertial, Rotating> const rigid_motion(
            rigid_transformation,
            angular_velocity_of_to_frame,
            velocity_of_to_frame_origin);
        Vector<Acceleration, Inertial> const
            acceleration_of_to_frame_origin;
        return AcceleratedRigidMotion<Inertial, Rotating>(
            rigid_motion,
            angular_acceleration_of_to_frame,
            acceleration_of_to_frame_origin);
      }));

  // No gravity.
  Vector<Acceleration, Inertial> const gravitational_acceleration;
  SpecificEnergy const gravitational_potential;
  EXPECT_CALL(mock_frame_, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(gravitational_acceleration));
  EXPECT_CALL(mock_frame_, GravitationalPotential(_, _))
      .WillRepeatedly(Return(gravitational_potential));

  // The test point is on the x axis.
  DegreesOfFreedom<Rotating> const initial_state_in_rotating_frame = {
      Rotating::origin +
          general_displacement_z_ + Displacement<Rotating>({100 * Metre,
                                                            0 * Metre,
                                                            0 * Metre}),
      general_velocity_};
  DegreesOfFreedom<Inertial> const initial_state_in_inertial_frame =
      mock_frame_.MotionOfThisFrame(t0_).rigid_motion().Inverse()(
          initial_state_in_rotating_frame);

  // The time interval for evaluating the first order effect.
  Time const Δt = 1 * Milli(Second);

  Position<Inertial> const final_position_in_inertial_frame =
      initial_state_in_inertial_frame.position() +
      initial_state_in_inertial_frame.velocity() * Δt;

  Position<Rotating> const final_position_in_rotating_frame =
      mock_frame_.MotionOfThisFrame(t0_ + Δt)
          .rigid_motion()
          .rigid_transformation()(final_position_in_inertial_frame);

  Position<Rotating> const first_order_final_position_in_rotating_frame =
      initial_state_in_rotating_frame.position() +
      initial_state_in_rotating_frame.velocity() * Δt;

  Displacement<Rotating> const higher_order_effect =
      final_position_in_rotating_frame -
      first_order_final_position_in_rotating_frame;

  // The second order effect is the Euler acceleration, the higher order effects
  // are irrelevant.  This computation only depends on the stub motion defined
  // above.
  EXPECT_THAT(higher_order_effect,
              Componentwise(IsNear(214_(1) * Nano(Metre)),
                            IsNear(-0.5_(1) * Milli(Metre)),
                            AlmostEquals(0 * Metre, 0)));

  // The centrifugal acceleration matches that computed based on the motion to
  // the second order.
  EXPECT_THAT(
      mock_frame_.GeometricAcceleration(t0_, initial_state_in_rotating_frame) *
          Pow<2>(Δt) / 2,
      AlmostEquals(Displacement<Rotating>({0 * Metre,
                                           -0.5 * Milli(Metre),
                                           0 * Metre}),
                   0));

  // No Euler acceleration when at rest.
  EXPECT_THAT(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                  t0_, initial_state_in_rotating_frame.position()),
              AlmostEquals(Vector<Acceleration, Rotating>(), 0));

  // No (scalar) Euler potential.
  EXPECT_THAT(FirstOrderAccelerationFromPotential(
                  t0_,
                  initial_state_in_rotating_frame.position(),
                  /*Δl=*/1 * Micro(Metre)),
              AlmostEquals(Vector<Acceleration, Rotating>(), 0));
}

// A bench vice with a right-hand screw.  The x axis is along the screw,
// pointing toward the operator, and the z axis is vertical pointing up.  The
// vice is being untightened (the handle is on the top and turns to the left).
// In the frame of the handle, a point on top of the moving jaw has a circular
// trajectory, while a point on top of the fixed jaw has an helical trajectory.
// This in turns determines two distinct Frenet frames which share the same
// normal (directed along -z) but different tangents: for the fixed jaw, the
// tangent is directed along the lead angle of the screw.
TEST_F(ReferenceFrameTest, FrenetFrame) {
  using Handle = Frame<struct HangleTag, geometry::Inertial>;
  using FixedJaw = Frame<struct FixedJawTag, geometry::Inertial>;
  using MovingJaw = Frame<struct MovingJawTag, geometry::Inertial>;

  StrictMock<MockRigidReferenceFrame<FixedJaw, MovingJaw>> jaw_to_jaw_frame;
  StrictMock<MockRigidReferenceFrame<MovingJaw, Handle>> jaw_to_handle_frame;

  Speed const v = 1 * Centi(Metre) / Second;
  AngularFrequency const ω = 1 * Radian / Second;

  EXPECT_CALL(jaw_to_jaw_frame, ToThisFrameAtTime(_))
      .WillRepeatedly(Invoke([this, v](Instant const& t) {
        auto const velocity_of_to_frame_origin =
            Velocity<FixedJaw>({v,
                                0 * Metre / Second,
                                0 * Metre / Second});
        auto const position_of_to_frame_origin =
            FixedJaw::origin +
            Displacement<FixedJaw>({v * (t - t0_),
                                    0 * Metre,
                                    0 * Metre});
        RigidTransformation<FixedJaw, MovingJaw> const
            rigid_transformation(
                /*from_origin=*/position_of_to_frame_origin,
                /*to_origin=*/MovingJaw::origin,
                OrthogonalMap<FixedJaw, MovingJaw>::Identity());
        AngularVelocity<FixedJaw> const angular_velocity_of_to_frame;
        return RigidMotion<FixedJaw, MovingJaw>(
            rigid_transformation,
            angular_velocity_of_to_frame,
            velocity_of_to_frame_origin);
      }));

  EXPECT_CALL(jaw_to_jaw_frame, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([&jaw_to_jaw_frame](Instant const& t) {
        Bivector<AngularAcceleration, FixedJaw> const
            angular_acceleration_of_to_frame;
        Vector<Acceleration, FixedJaw> const acceleration_of_to_frame_origin;
        return AcceleratedRigidMotion<FixedJaw, MovingJaw>(
            jaw_to_jaw_frame.ToThisFrameAtTime(t),
            angular_acceleration_of_to_frame,
            acceleration_of_to_frame_origin);
      }));

  EXPECT_CALL(jaw_to_handle_frame, ToThisFrameAtTime(_))
      .WillRepeatedly(Invoke([this, ω](Instant const& t) {
        AngularVelocity<MovingJaw> const angular_velocity_of_to_frame(
            {ω, 0 * Radian / Second, 0 * Radian / Second});
        Rotation<MovingJaw, Handle> const rotation(
            ω * (t - t0_),
            angular_velocity_of_to_frame,
            DefinesFrame<Handle>{});
        RigidTransformation<MovingJaw, Handle> const
            rigid_transformation(
                /*from_origin=*/MovingJaw::origin,
                /*to_origin=*/Handle::origin,
                rotation.Forget<OrthogonalMap>());
        Velocity<MovingJaw> const velocity_of_to_frame_origin;
        return RigidMotion<MovingJaw, Handle>(
            rigid_transformation,
            angular_velocity_of_to_frame,
            velocity_of_to_frame_origin);
      }));

  EXPECT_CALL(jaw_to_handle_frame, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([&jaw_to_handle_frame](Instant const& t) {
        Bivector<AngularAcceleration, MovingJaw> const
            angular_acceleration_of_to_frame;
        Vector<Acceleration, MovingJaw> const
            acceleration_of_to_frame_origin;
        return AcceleratedRigidMotion<MovingJaw, Handle>(
            jaw_to_handle_frame.ToThisFrameAtTime(t),
            angular_acceleration_of_to_frame,
            acceleration_of_to_frame_origin);
      }));

  // No gravity (or put it another way, the gravity is compensated by a reaction
  // force).
  EXPECT_CALL(jaw_to_jaw_frame, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(Vector<Acceleration, FixedJaw>()));
  EXPECT_CALL(jaw_to_handle_frame, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(Vector<Acceleration, MovingJaw>()));

  // Fixed points on top of each jaw.
  Length const r = 0.2 * Metre;
  DegreesOfFreedom<FixedJaw> fixed_jaw_degrees_of_freedom(
      FixedJaw::origin +
          Displacement<FixedJaw>({0 * Metre, 0 * Metre, r}),
      FixedJaw::unmoving);
  DegreesOfFreedom<MovingJaw> moving_jaw_degrees_of_freedom(
      MovingJaw::origin +
          Displacement<MovingJaw>({0 * Metre, 0 * Metre, r}),
      MovingJaw::unmoving);

  auto const fixed_jaw_frenet_frame = jaw_to_handle_frame.FrenetFrame(
      t0_,
      jaw_to_handle_frame.ToThisFrameAtTime(t0_)(
          jaw_to_jaw_frame.ToThisFrameAtTime(t0_)(
              fixed_jaw_degrees_of_freedom)));
  auto const moving_jaw_frenet_frame = jaw_to_handle_frame.FrenetFrame(
      t0_,
      jaw_to_handle_frame.ToThisFrameAtTime(t0_)(
          moving_jaw_degrees_of_freedom));

  Vector<double, Handle> fixed_jaw_tangent =
      fixed_jaw_frenet_frame(Vector<double, Frenet<FixedJaw>>({1, 0, 0}));
  Vector<double, Handle> fixed_jaw_normal =
      fixed_jaw_frenet_frame(Vector<double, Frenet<FixedJaw>>({0, 1, 0}));

  Vector<double, Handle> moving_jaw_tangent =
      moving_jaw_frenet_frame(Vector<double, Frenet<FixedJaw>>({1, 0, 0}));
  Vector<double, Handle> moving_jaw_normal =
      moving_jaw_frenet_frame(Vector<double, Frenet<FixedJaw>>({0, 1, 0}));

  Length const pitch = v * (2 * π * Radian / ω);
  Angle const lead_angle = ArcTan(pitch, 2 * π * r);
  EXPECT_THAT(fixed_jaw_normal,
              Componentwise(VanishesBefore(1, 1),
                            VanishesBefore(1, 0),
                            AlmostEquals(-1, 2)));
  EXPECT_THAT(fixed_jaw_tangent,
              Componentwise(Lt(0), Gt(0), VanishesBefore(1, 1)));
  EXPECT_THAT(AngleBetween(fixed_jaw_tangent, moving_jaw_tangent),
              AlmostEquals(lead_angle, 13));

  EXPECT_THAT(moving_jaw_normal, Componentwise(0, 0, -1));
  EXPECT_THAT(moving_jaw_tangent, Componentwise(0, 1, 0));
}

}  // namespace physics
}  // namespace principia
