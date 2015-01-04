
#pragma once

#include "physics/n_body_system.hpp"

#include "gmock/gmock.h"

namespace principia {
namespace physics {

template<typename InertialFrame>
class MockNBodySystem : public NBodySystem<InertialFrame> {
 public:
  MockNBodySystem() = default;

  MOCK_CONST_METHOD6_T(
      Integrate,
      void(SymplecticIntegrator<Length, Speed> const& integrator,
           Instant const& tmax,
           Time const& Δt,
           int const sampling_period,
           bool const tmax_is_exact,
           typename NBodySystem<InertialFrame>::Trajectories const&
               trajectories));
};

}  // namespace physics
}  // namespace principia
