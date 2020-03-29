
#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "physics/dynamic_frame.hpp"
#include "physics/rigid_motion.hpp"

namespace principia {
namespace physics {
namespace internal_dynamic_frame {

template<typename InertialFrame, typename ThisFrame>
class MockDynamicFrame : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  MOCK_CONST_METHOD1_T(ToThisFrameAtTime,
                       RigidMotion<InertialFrame, ThisFrame>(Instant const& t));
  MOCK_CONST_METHOD1_T(FromThisFrameAtTime,
                       RigidMotion<ThisFrame, InertialFrame>(Instant const& t));

  MOCK_CONST_METHOD0(t_min, Instant());
  MOCK_CONST_METHOD0(t_max, Instant());

  MOCK_CONST_METHOD2_T(
      GeometricAcceleration,
      Vector<Acceleration, ThisFrame>(
          Instant const& t,
          DegreesOfFreedom<ThisFrame> const& degrees_of_freedom));

  using Rot = Rotation<Frenet<ThisFrame>, ThisFrame>;

  MOCK_CONST_METHOD2_T(
      FrenetFrame,
      Rot(Instant const& t,
          DegreesOfFreedom<ThisFrame> const& degrees_of_freedom));

  MOCK_CONST_METHOD1_T(
      WriteToMessage,
      void(not_null<serialization::DynamicFrame*> message));

 private:
  MOCK_CONST_METHOD2_T(
      GravitationalAcceleration,
      Vector<Acceleration, InertialFrame>(Instant const& t,
                                          Position<InertialFrame> const& q));
  MOCK_CONST_METHOD1_T(
      MotionOfThisFrame,
      AcceleratedRigidMotion<InertialFrame, ThisFrame>(Instant const& t));
};

}  // namespace internal_dynamic_frame

using internal_dynamic_frame::MockDynamicFrame;

}  // namespace physics
}  // namespace principia
