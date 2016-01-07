#pragma once

#include "kepler_orbit.hpp"

namespace principia {
namespace physics {

template<typename Frame>
KeplerOrbit<Frame>::KeplerOrbit(
    TwoBodySystem const& system,
    Instant const& epoch,
    KeplerianElements<Frame> const& elements_at_epoch)
    : system_gravitational_parameter_(
          system.primary->gravitational_parameter() +
          system.secondary->is_massless
              ? GravitationalParameter{}
              : CHECK_NOTNULL(
                    dynamic_cast<MassiveBody const*>(&*system.secondary))->
                        gravitational_parameter()),
      epoch_(epoch),
      elements_at_epoch_(element_at_epoch) {}

template<typename Frame>
KeplerOrbit<Frame>::Evaluate(Instant const& t) const {}


}  // namespace physics
}  // namespace principia
