
#pragma once

#include "physics/jacobi_coordinates.hpp"

#include <vector>

namespace principia {
namespace physics {
namespace internal_jacobi_coordinates {

using geometry::Instant;
using geometry::Velocity;

template<typename Frame>
JacobiCoordinates<Frame>::JacobiCoordinates(MassiveBody const& primary) {
  static DegreesOfFreedom<PrimocentricFrame> const motionless_origin = {
      PrimocentricFrame::origin, Velocity<PrimocentricFrame>()};
  primocentric_dof_.emplace_back(motionless_origin);
  system_barycentre_.Add(primocentric_dof_.back(),
                         primary.gravitational_parameter());
}

template<typename Frame>
void JacobiCoordinates<Frame>::Add(
    MassiveBody const& body,
    RelativeDegreesOfFreedom<Frame> const& dof_relative_to_system) {
  primocentric_dof_.emplace_back(system_barycentre_.Get() +
                                 id_fp_(dof_relative_to_system));
  system_barycentre_.Add(primocentric_dof_.back(),
                         body.gravitational_parameter());
}

template<typename Frame>
void JacobiCoordinates<Frame>::Add(
    MassiveBody const& body,
    KeplerianElements<Frame> const& osculating_elements_relative_to_system) {
  Instant const epoch;
  Add(body,
      KeplerOrbit<Frame>(/*primary=*/System(),
                         /*secondary=*/body,
                         osculating_elements_relative_to_system,
                         epoch).StateVectors(epoch));
}

template<typename Frame>
MassiveBody JacobiCoordinates<Frame>::System() const {
  // A point mass.
  return MassiveBody(MassiveBody::Parameters(system_barycentre_.weight()));
}

template<typename Frame>
std::vector<RelativeDegreesOfFreedom<Frame>>
JacobiCoordinates<Frame>::BarycentricDegreesOfFreedom() const {
  std::vector<RelativeDegreesOfFreedom<Frame>> result;
  DegreesOfFreedom<PrimocentricFrame> system_barycentre =
      system_barycentre_.Get();
  for (auto const& dof : primocentric_dof_) {
    result.emplace_back(id_pf_(dof - system_barycentre));
  }
  return result;
}

template<typename Frame>
Identity<typename JacobiCoordinates<Frame>::PrimocentricFrame, Frame> const
    JacobiCoordinates<Frame>::id_pf_;
template<typename Frame>
Identity<Frame, typename JacobiCoordinates<Frame>::PrimocentricFrame> const
    JacobiCoordinates<Frame>::id_fp_;

}  // namespace internal_jacobi_coordinates
}  // namespace physics
}  // namespace principia
