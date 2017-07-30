
#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/manœuvre.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_manœuvre {

template<typename InertialFrame, typename Frame>
class MockManœuvre : public Manœuvre<InertialFrame, Frame>{
 public:
  MockManœuvre(
      Force const& thrust,
      Mass const& initial_mass,
      SpecificImpulse const& specific_impulse,
      Vector<double, Frenet<Frame>> const& direction,
      not_null<std::unique_ptr<DynamicFrame<InertialFrame, Frame> const>> frame,
      bool const is_inertially_fixed)
      : Manœuvre<InertialFrame, Frame>(thrust,
                                       initial_mass,
                                       specific_impulse,
                                       direction,
                                       std::move(frame),
                                       is_inertially_fixed) {}

  MOCK_CONST_METHOD0_T(InertialDirection, Vector<double, InertialFrame>());

  MOCK_CONST_METHOD0_T(FrenetFrame,
                       OrthogonalMap<Frenet<Frame>, InertialFrame>());
};

}  // namespace internal_manœuvre

using internal_manœuvre::MockManœuvre;

}  // namespace ksp_plugin
}  // namespace principia
