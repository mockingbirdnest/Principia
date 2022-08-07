
#include "physics/dynamic_frame.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace physics {
namespace internal_dynamic_frame {

using geometry::AngularVelocity;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::Frame;
using geometry::InfiniteFuture;
using geometry::InfinitePast;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::RigidTransformation;
using geometry::Velocity;
using quantities::AngularAcceleration;
using quantities::AngularFrequency;
using quantities::GravitationalParameter;
using quantities::Sqrt;
using quantities::Time;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::IsNear;
using testing_utilities::operator""_;
using ::testing::Invoke;
using ::testing::Return;
using ::testing::StrictMock;
using ::testing::_;
namespace si = quantities::si;

namespace {

using Circular = Frame<enum class CircularTag, geometry::Inertial>;
using Helical = Frame<enum class HelicalTag, geometry::Inertial>;

Vector<Acceleration, Circular> Gravity(Instant const& t,
                                       Position<Circular> const& q) {
  Displacement<Circular> const r = q - Circular::origin;
  auto const r² = r.Norm²();
  return -si::Unit<GravitationalParameter> * r / (Sqrt(r²) * r²);
}

// An inertial frame.
template<typename OtherFrame, typename ThisFrame>
class InertialFrame : public DynamicFrame<OtherFrame, ThisFrame> {
 public:
  InertialFrame(
      DegreesOfFreedom<OtherFrame> const& origin_degrees_of_freedom_at_epoch,
      Instant const& epoch,
      OrthogonalMap<OtherFrame, ThisFrame> const& orthogonal_map,
      std::function<Vector<Acceleration, OtherFrame>(
          Instant const& t,
          Position<OtherFrame> const& q)> gravity);

  Instant t_min() const override {
    return InfinitePast;
  }

  Instant t_max() const override {
    return InfiniteFuture;
  }

  RigidMotion<OtherFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::DynamicFrame*> message) const override;

 private:
  Vector<Acceleration, OtherFrame> GravitationalAcceleration(
      Instant const& t,
      Position<OtherFrame> const& q) const override;
  AcceleratedRigidMotion<OtherFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const override;

  DegreesOfFreedom<OtherFrame> const origin_degrees_of_freedom_at_epoch_;
  Instant const epoch_;
  OrthogonalMap<OtherFrame, ThisFrame> const orthogonal_map_;
  std::function<Vector<Acceleration, OtherFrame>(
      Instant const& t,
      Position<OtherFrame> const& q)> gravity_;
};

template<typename OtherFrame, typename ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::InertialFrame(
    DegreesOfFreedom<OtherFrame> const& origin_degrees_of_freedom_at_epoch,
    Instant const& epoch,
    OrthogonalMap<OtherFrame, ThisFrame> const& orthogonal_map,
    std::function<Vector<Acceleration, OtherFrame>(
        Instant const& t,
        Position<OtherFrame> const& q)> gravity)
    : origin_degrees_of_freedom_at_epoch_(origin_degrees_of_freedom_at_epoch),
      epoch_(epoch),
      orthogonal_map_(orthogonal_map),
      gravity_(std::move(gravity)) {}

template<typename OtherFrame, typename ThisFrame>
RigidMotion<OtherFrame, ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  return RigidMotion<OtherFrame, ThisFrame>(
             RigidTransformation<OtherFrame, ThisFrame>(
                 origin_degrees_of_freedom_at_epoch_.position() +
                     (t - epoch_) *
                         origin_degrees_of_freedom_at_epoch_.velocity(),
                 ThisFrame::origin,
                 orthogonal_map_),
             OtherFrame::nonrotating,
             origin_degrees_of_freedom_at_epoch_.velocity());
}

template<typename OtherFrame, typename ThisFrame>
void InertialFrame<OtherFrame, ThisFrame>::WriteToMessage(
    not_null<serialization::DynamicFrame*> message) const {}

template<typename OtherFrame, typename ThisFrame>
Vector<Acceleration, OtherFrame>
InertialFrame<OtherFrame, ThisFrame>::GravitationalAcceleration(
    Instant const& t,
    Position<OtherFrame> const& q) const {
  return gravity_(t, q);
}

template<typename OtherFrame, typename ThisFrame>
AcceleratedRigidMotion<OtherFrame, ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::MotionOfThisFrame(
    Instant const& t) const {
  return AcceleratedRigidMotion<OtherFrame, ThisFrame>(
      ToThisFrameAtTime(t),
      /*angular_acceleration_of_to_frame=*/{},
      /*acceleration_of_to_frame_origin=*/{});
}

template<typename InertialFrame, typename ThisFrame>
class MockDynamicFrame : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  MockDynamicFrame() = default;

  MOCK_METHOD((RigidMotion<InertialFrame, ThisFrame>),
              ToThisFrameAtTime,
              (Instant const& t),
              (const, override));
  MOCK_METHOD((RigidMotion<ThisFrame, InertialFrame>),
              FromThisFrameAtTime,
              (Instant const& t),
              (const, override));

  MOCK_METHOD(Instant, t_min, (), (const, override));
  MOCK_METHOD(Instant, t_max, (), (const, override));

  MOCK_METHOD(void,
              WriteToMessage,
              (not_null<serialization::DynamicFrame*> message),
              (const, override));

  MOCK_METHOD((Vector<Acceleration, InertialFrame>),
              GravitationalAcceleration,
              (Instant const& t, Position<InertialFrame> const& q),
              (const, override));
  MOCK_METHOD((AcceleratedRigidMotion<InertialFrame, ThisFrame>),
              MotionOfThisFrame,
              (Instant const& t),
              (const, override));
};

}  // namespace

class DynamicFrameTest : public testing::Test {
 protected:
  using Inertial = Frame<enum class InertialFrameTag, geometry::Inertial>;
  using Rotating = Frame<enum class TestFrameTag, Arbitrary>;

  StrictMock<MockDynamicFrame<Inertial, Rotating>> mock_frame_;

  DegreesOfFreedom<Circular> circular_degrees_of_freedom_ = {
      Circular::origin +
          Displacement<Circular>({1 * Metre, 0 * Metre, 0 * Metre}),
      Velocity<Circular>(
          {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second})};

  InertialFrame<Circular, Helical> helix_frame_ =
      InertialFrame<Circular, Helical>(
          {Circular::origin,
           Velocity<Circular>(
               {0 * Metre / Second, 0 * Metre / Second, 1 * Metre / Second})},
          /*epoch=*/Instant() ,
          OrthogonalMap<Circular, Helical>::Identity(),
          &Gravity);
};

// A frame in uniform rotation around the origin.  The test point is at the
// origin and in motion along the x axis.  The acceleration is purely due to
// Coriolis.  The motion elements that don't have specific values have no effect
// on the acceleration.
TEST_F(DynamicFrameTest, CoriolisAcceleration) {
  Instant const t0;

  EXPECT_CALL(mock_frame_, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([t0](Instant const& t) {
        AngularFrequency const ω = 10 * Radian / Second;
        AngularVelocity<Inertial> const angular_velocity_of_to_frame(
            {0 * Radian / Second, 0 * Radian / Second, ω});
        Rotation<Inertial, Rotating> const rotation(
            ω * (t - t0),
            angular_velocity_of_to_frame,
            DefinesFrame<Rotating>{});
        RigidTransformation<Inertial, Rotating> const
            rigid_transformation(
                /*from_origin=*/Inertial::origin,
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
  EXPECT_CALL(mock_frame_, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(gravitational_acceleration));

  // The velocity is along the x axis.
  DegreesOfFreedom<Rotating> const initial_state_in_rotating_frame =
      {Rotating::origin,
       Velocity<Rotating>({100 * Metre / Second,
                           0 * Metre / Second,
                           0 * Metre / Second})};
  DegreesOfFreedom<Inertial> const initial_state_in_inertial_frame =
      mock_frame_.MotionOfThisFrame(t0).rigid_motion().Inverse()(
          initial_state_in_rotating_frame);

  // The time interval for evaluating the first order effect.
  Time const Δt = 1 * Milli(Second);

  Position<Inertial> const final_position_in_inertial_frame =
      initial_state_in_inertial_frame.position() +
      initial_state_in_inertial_frame.velocity() * Δt;

  Position<Rotating> const final_position_in_rotating_frame =
      mock_frame_.MotionOfThisFrame(t0 + Δt)
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
                            AlmostEquals(0 * Metre, 0)));

  // The Coriolis acceleration matches that computed based on the motion to the
  // second order.  This validates that we don't have sign errors in the actual
  // frame implementation.
  EXPECT_THAT(
      mock_frame_.GeometricAcceleration(t0, initial_state_in_rotating_frame) *
          Pow<2>(Δt) / 2,
      AlmostEquals(Displacement<Rotating>({0 * Metre,
                                           -1 * Milli(Metre),
                                           0 * Metre}),
                   0));

  // No Coriolis acceleration when at rest.
  EXPECT_THAT(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                  t0, initial_state_in_rotating_frame.position()),
              AlmostEquals(Vector<Acceleration, Rotating>(), 0));
}

// A frame in uniform rotation around the origin.  The test point is on the x
// axis.  The acceleration is purely due to centrifugal effects.  The motion
// elements that don't have specific values have no effect on the acceleration.
TEST_F(DynamicFrameTest, CentrifugalAcceleration) {
  Instant const t0;

  EXPECT_CALL(mock_frame_, MotionOfThisFrame(_))
      .WillRepeatedly(Invoke([t0](Instant const& t) {
        AngularFrequency const ω = 10 * Radian / Second;
        AngularVelocity<Inertial> const angular_velocity_of_to_frame(
            {0 * Radian / Second, 0 * Radian / Second, ω});
        Rotation<Inertial, Rotating> const rotation(
            ω * (t - t0),
            angular_velocity_of_to_frame,
            DefinesFrame<Rotating>{});
        RigidTransformation<Inertial, Rotating> const
            rigid_transformation(
                /*from_origin=*/Inertial::origin,
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
  EXPECT_CALL(mock_frame_, GravitationalAcceleration(_, _))
      .WillRepeatedly(Return(gravitational_acceleration));

  // The test point is on the x axis.
  DegreesOfFreedom<Rotating> const initial_state_in_rotating_frame = {
      Rotating::origin + Displacement<Rotating>({100 * Metre,
                                                           0 * Metre,
                                                           0 * Metre}),
      Rotating::unmoving};
  DegreesOfFreedom<Inertial> const initial_state_in_inertial_frame =
      mock_frame_.MotionOfThisFrame(t0).rigid_motion().Inverse()(
          initial_state_in_rotating_frame);

  // The time interval for evaluating the first order effect.
  Time const Δt = 1 * Milli(Second);

  Position<Inertial> const final_position_in_inertial_frame =
      initial_state_in_inertial_frame.position() +
      initial_state_in_inertial_frame.velocity() * Δt;

  Position<Rotating> const final_position_in_rotating_frame =
      mock_frame_.MotionOfThisFrame(t0 + Δt)
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
  // the second order.  This validates that we don't have sign errors in the
  // actual frame implementation.
  EXPECT_THAT(
      mock_frame_.GeometricAcceleration(t0, initial_state_in_rotating_frame) *
          Pow<2>(Δt) / 2,
      AlmostEquals(Displacement<Rotating>({5 * Milli(Metre),
                                           0 * Metre,
                                           0 * Metre}),
                   0));

  // The centrifugal acceleration shows up for a point at rest.
  EXPECT_THAT(mock_frame_.RotationFreeGeometricAccelerationAtRest(
                  t0, initial_state_in_rotating_frame.position()) *
                  Pow<2>(Δt) / 2,
              AlmostEquals(Displacement<Rotating>({5 * Milli(Metre),
                                                   0 * Metre,
                                                   0 * Metre}),
                           0));
}

TEST_F(DynamicFrameTest, Helix) {
  auto helix_frenet_frame = helix_frame_.FrenetFrame(
      Instant(),
      helix_frame_.ToThisFrameAtTime(Instant())(circular_degrees_of_freedom_));
  Vector<double, Helical> tangent =
      helix_frenet_frame(Vector<double, Frenet<Helical>>({1, 0, 0}));
  Vector<double, Helical> normal =
      helix_frenet_frame(Vector<double, Frenet<Helical>>({0, 1, 0}));
  EXPECT_THAT(normal, Componentwise(-1, 0, 0));
  EXPECT_THAT(tangent,
              Componentwise(0,
                            AlmostEquals(Sqrt(0.5), 1),
                            AlmostEquals(-Sqrt(0.5), 1)));
}

}  // namespace internal_dynamic_frame
}  // namespace physics
}  // namespace principia
