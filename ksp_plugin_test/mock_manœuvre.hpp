
#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/manœuvre.hpp"

namespace principia {
namespace ksp_plugin {

template <typename InertialFrame, typename Frame>
class MockManœuvre : public Manœuvre<InertialFrame, Frame>{
 public:
  MockManœuvre(
      Force const& thrust,
      Mass const& initial_mass,
      SpecificImpulse const& specific_impulse,
      Vector<double, Frenet<Frame>> const& direction,
      not_null<std::unique_ptr<DynamicFrame<InertialFrame, Frame> const>> frame)
      : Manœuvre(thrust,
                 initial_mass,
                 specific_impulse,
                 direction,
                 std::move(frame)) {}

  MOCK_CONST_METHOD0_T(inertial_direction, Vector<double, InertialFrame>());
};

}  // namespace ksp_plugin
}  // namespace principia
