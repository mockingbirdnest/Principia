
#include "physics/dynamic_frame.hpp"

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"

namespace principia {
namespace physics {
namespace internal_dynamic_frame {

using geometry::AngularVelocity;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::Frame;
using geometry::Inertial;
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
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using ::testing::Invoke;
using ::testing::Return;
using ::testing::StrictMock;
using ::testing::_;
namespace si = quantities::si;

namespace {

using Circular = Frame<enum class CircularTag, Inertial>;
using Helical = Frame<enum class HelicalTag, Inertial>;

Vector<Acceleration, Circular> Gravity(Instant const& t,
                                       Position<Circular> const& q) {
  Displacement<Circular> const r = q - Circular::origin;
  auto const r² = r.Norm²();
  return -si::Unit<GravitationalParameter> * r / (Sqrt(r²) * r²);
}

// An inertial frame.
template<typename OtherFrame, typename ThisFrame>
class FakeDynamicFrame : public DynamicFrame<OtherFrame, ThisFrame> {
 public:
  FakeDynamicFrame(
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
FakeDynamicFrame<OtherFrame, ThisFrame>::FakeDynamicFrame(
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
FakeDynamicFrame<OtherFrame, ThisFrame>::ToThisFrameAtTime(
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
void FakeDynamicFrame<OtherFrame, ThisFrame>::WriteToMessage(
    not_null<serialization::DynamicFrame*> message) const {}

template<typename OtherFrame, typename ThisFrame>
Vector<Acceleration, OtherFrame>
FakeDynamicFrame<OtherFrame, ThisFrame>::GravitationalAcceleration(
    Instant const& t,
    Position<OtherFrame> const& q) const {
  return gravity_(t, q);
}

template<typename OtherFrame, typename ThisFrame>
AcceleratedRigidMotion<OtherFrame, ThisFrame>
FakeDynamicFrame<OtherFrame, ThisFrame>::MotionOfThisFrame(
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
  using InertialFrame = Frame<enum class InertialFrameTag,
                              Inertial>;
  using TestFrame = Frame<enum class TestFrameTag,
                          Arbitrary>;

  StrictMock<MockDynamicFrame<InertialFrame, TestFrame>> mock_frame_;

  DegreesOfFreedom<Circular> circular_degrees_of_freedom_ = {
      Circular::origin +
          Displacement<Circular>({1 * Metre, 0 * Metre, 0 * Metre}),
      Velocity<Circular>(
          {0 * Metre / Second, 1 * Metre / Second, 0 * Metre / Second})};

  FakeDynamicFrame<Circular, Helical> helix_frame_ =
      FakeDynamicFrame<Circular, Helical>(
          {Circular::origin,
           Velocity<Circular>(
               {0 * Metre / Second, 0 * Metre / Second, 1 * Metre / Second})},
          /*epoch=*/Instant() ,
          OrthogonalMap<Circular, Helical>::Identity(),
          &Gravity);
};

// A frame in uniform rotation around |from_origin|.  The test point is at the
// origin and in motion.  The acceleration is purely due to Coriolis.  The
// motion elements that don't have specific values have no effect on the
// acceleration.
TEST_F(DynamicFrameTest, CoriolisAcceleration) {
  Instant const t0;

  // The velocity is opposed to the motion and away from the centre.
  DegreesOfFreedom<TestFrame> const point_dof =
      {TestFrame::origin,
       Velocity<TestFrame>({50 * Metre / Second,
                            -100 * Metre / Second,
                            0 * Metre / Second})};

  EXPECT_CALL(mock_frame_, MotionOfThisFrame(_))
      .WillOnce(Invoke([&point_dof, t0](Instant const& t) {
        AngularFrequency const ω = 10 * Radian / Second;
        AngularVelocity<InertialFrame> const angular_velocity_of_to_frame(
            {0 * Radian / Second, 0 * Radian / Second, ω});
        Rotation<InertialFrame, TestFrame> const rotation(
            ω * (t - t0),
            angular_velocity_of_to_frame,
            DefinesFrame<TestFrame>{});
        Position<InertialFrame> const from_origin =
            InertialFrame::origin +
            Displacement<InertialFrame>({2 * Metre, 1 * Metre, 0 * Metre});
        Position<TestFrame> const to_origin = point_dof.position();
        RigidTransformation<InertialFrame, TestFrame> const
            rigid_transformation(
                from_origin, to_origin, rotation.Forget<OrthogonalMap>());
        Velocity<InertialFrame> const velocity_of_to_frame_origin;
        RigidMotion<InertialFrame, TestFrame> const rigid_motion(
            rigid_transformation,
            angular_velocity_of_to_frame,
            velocity_of_to_frame_origin);
        Bivector<AngularAcceleration, InertialFrame> const
            angular_acceleration_of_to_frame;
        Vector<Acceleration, InertialFrame> const
            acceleration_of_to_frame_origin;
        return AcceleratedRigidMotion<InertialFrame, TestFrame>(
            rigid_motion,
            angular_acceleration_of_to_frame,
            acceleration_of_to_frame_origin);
      }));

  Vector<Acceleration, InertialFrame> const gravitational_acceleration;
  EXPECT_CALL(mock_frame_, GravitationalAcceleration(_, _))
      .WillOnce(Return(gravitational_acceleration));

  // The Coriolis acceleration is towards the centre and opposed to the motion.
  EXPECT_THAT(mock_frame_.GeometricAcceleration(t0, point_dof),
              AlmostEquals(Vector<Acceleration, TestFrame>(
                               {-2000 * Metre / Pow<2>(Second),
                                -1000 * Metre / Pow<2>(Second),
                                0 * Metre / Pow<2>(Second)}),
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
