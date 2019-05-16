
#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/manœuvre.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_manœuvre {

template<typename InertialFrame, typename Frame>
class MockManœuvre : public Manœuvre<InertialFrame, Frame>{
 public:
  using Intensity = typename Manœuvre<InertialFrame, Frame>::Intensity;
  using Timing = typename Manœuvre<InertialFrame, Frame>::Timing;

  MockManœuvre(
      Force const& thrust,
      Mass const& initial_mass,
      SpecificImpulse const& specific_impulse,
      Intensity const& intensity,
      Timing const& timing,
      not_null<std::shared_ptr<DynamicFrame<InertialFrame, Frame> const>> frame,
      bool const is_inertially_fixed)
    : Manœuvre<InertialFrame, Frame>(thrust,
                                      initial_mass,
                                      specific_impulse,
                                      intensity,
                                      timing,
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
