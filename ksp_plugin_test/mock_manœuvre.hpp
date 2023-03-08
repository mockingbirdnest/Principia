#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/manœuvre.hpp"

namespace principia {
namespace ksp_plugin {
namespace _manœuvre {
namespace internal {

template<typename InertialFrame, typename Frame>
class MockManœuvre : public Manœuvre<InertialFrame, Frame>{
 public:
  using Intensity = typename Manœuvre<InertialFrame, Frame>::Intensity;
  using Timing = typename Manœuvre<InertialFrame, Frame>::Timing;
  using Burn = typename Manœuvre<InertialFrame, Frame>::Burn;

  MockManœuvre(Mass const& initial_mass,
               Burn const& burn)
    : Manœuvre<InertialFrame, Frame>(initial_mass, burn) {}

  MOCK_METHOD((Vector<double, InertialFrame>),
              InertialDirection,
              (),
              (const, override));

  MOCK_METHOD((OrthogonalMap<Frenet<Frame>, InertialFrame>),
              FrenetFrame,
              (),
              (const, override));
};

}  // namespace internal

using internal::MockManœuvre;

}  // namespace _manœuvre
}  // namespace ksp_plugin
}  // namespace principia
