#pragma once

#include "physics/jacobi_coordinates.hpp"

namespace principia {

namespace physics {

template<typename Frame>
JacobiCoordinates<Frame>::JacobiCoordinates(MassiveBody const& primary) {
  DegreesOfFreedom<PrimocentricFrame> motionless_origin = {
      PrimocentricFrame::origin, Velocity<PrimocentricFrame>()};
  primocentric_dof_.push_back(motionless_origin);
  system_barycentre_.Add(primocentric_dof_.back(),
                         primary.gravitational_parameter());
}

template<typename Frame>
void JacobiCoordinates<Frame>::Add(
    MassiveBody const& body,
    RelativeDegreesOfFreedom<Frame> const& dof_wrt_system) {
  primocentric_dof_.emplace_back(system_barycentre_.Get() +
                                 id_fp_(dof_wrt_system));
  system_barycentre_.Add(primocentric_dof_.back(),
                         body.gravitational_parameter());
}

template<typename Frame>
void JacobiCoordinates<Frame>::Add(
    MassiveBody const& body,
    KeplerianElements<Frame> const& osculating_elements_wrt_system) {
  Instant const epoch;
  Add(body,
      KeplerOrbit<Frame>(/*primary=*/System(),
                         /*secondary=*/body,
                         osculating_elements_wrt_system,
                         epoch).StateVectors(epoch));
}

template<typename Frame>
MassiveBody JacobiCoordinates<Frame>::System() const {
  return MassiveBody(system_barycentre_.weight());
}

template<typename Frame>
std::vector<RelativeDegreesOfFreedom<Frame>>
JacobiCoordinates<Frame>::BarycentricCoordinates() const {
  std::vector<RelativeDegreesOfFreedom<Frame>> result;
  for (auto const& dof : primocentric_dof_) {
    result.emplace_back(id_pf_(dof - system_barycentre));
  }
  return result;
}

}  // namespace physics
}  // namespace principia
