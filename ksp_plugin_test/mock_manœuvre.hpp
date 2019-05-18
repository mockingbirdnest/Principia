
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
  using Burn = typename Manœuvre<InertialFrame, Frame>::Burn;

  MockManœuvre(Mass const& initial_mass,
               Burn const& burn)
    : Manœuvre<InertialFrame, Frame>(initial_mass, burn) {}

  MOCK_CONST_METHOD0_T(InertialDirection, Vector<double, InertialFrame>());

  MOCK_CONST_METHOD0_T(FrenetFrame,
                       OrthogonalMap<Frenet<Frame>, InertialFrame>());
};

}  // namespace internal_manœuvre

using internal_manœuvre::MockManœuvre;

}  // namespace ksp_plugin
}  // namespace principia
