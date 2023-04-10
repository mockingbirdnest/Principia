#pragma once

#include "physics/rigid_reference_frame.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "physics/rigid_motion.hpp"

namespace principia {
namespace physics {
namespace _rigid_reference_frame {
namespace internal {

using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;

template<typename InertialFrame, typename ThisFrame>
class MockRigidReferenceFrame : public RigidReferenceFrame<InertialFrame,
                                                           ThisFrame> {
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

  MOCK_METHOD(void,
              WriteToMessage,
              (not_null<serialization::ReferenceFrame*> message),
              (const, override));

  MOCK_METHOD((Vector<Acceleration, InertialFrame>),
              GravitationalAcceleration,
              (Instant const& t, Position<InertialFrame> const& q),
              (const, override));
  MOCK_METHOD(SpecificEnergy,
              GravitationalPotential,
              (Instant const& t, Position<InertialFrame> const& q),
              (const, override));

  MOCK_METHOD((AcceleratedRigidMotion<InertialFrame, ThisFrame>),
              MotionOfThisFrame,
              (Instant const& t),
              (const, override));
};

}  // namespace internal

using internal::MockRigidReferenceFrame;

}  // namespace _rigid_reference_frame
}  // namespace physics
}  // namespace principia
