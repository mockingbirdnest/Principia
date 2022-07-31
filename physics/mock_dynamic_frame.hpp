
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

  MOCK_METHOD((Vector<Acceleration, ThisFrame>),
              GeometricAcceleration,
              (Instant const& t,
               DegreesOfFreedom<ThisFrame> const& degrees_of_freedom),
              (const, override));

  using Rot = Rotation<Frenet<ThisFrame>, ThisFrame>;

  MOCK_METHOD(Rot,
              FrenetFrame,
              (Instant const& t,
               DegreesOfFreedom<ThisFrame> const& degrees_of_freedom),
              (const, override));

  MOCK_METHOD(void,
              WriteToMessage,
              (not_null<serialization::DynamicFrame*> message),
              (const, override));

 private:
  MOCK_METHOD((Vector<Acceleration, InertialFrame>),
              GravitationalAcceleration,
              (Instant const& t, Position<InertialFrame> const& q),
              (const, override));
  MOCK_METHOD(SpecificEnergy,
              GravitationalPotential,
              (Instant const& t, Position<InertialFrame> const& q),
              (const override));

  MOCK_METHOD((AcceleratedRigidMotion<InertialFrame, ThisFrame>),
              MotionOfThisFrame,
              (Instant const& t),
              (const, override));
};

}  // namespace internal_dynamic_frame

using internal_dynamic_frame::MockDynamicFrame;

}  // namespace physics
}  // namespace principia
