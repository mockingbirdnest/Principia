
#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/identity.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"

namespace principia {
namespace physics {
namespace internal_jacobi_coordinates {

using base::not_null;
using geometry::BarycentreCalculator;
using geometry::Identity;
using quantities::GravitationalParameter;

// An utility for converting a linearly ordered system of massive bodies given
// in Jacobi coordinates to barycentric coordinates.
template<typename Frame>
class JacobiCoordinates final {
 public:
  explicit JacobiCoordinates(MassiveBody const& primary);

  // Adds |body| with the given |DegreesOfFreedom| with respect to the
  // barycentre of the existing bodies.
  void Add(MassiveBody const& body,
           RelativeDegreesOfFreedom<Frame> const& dof_relative_to_system);

  // Adds |body| with the |RelativeDegreesOfFreedom| of a |KeplerOrbit| with the
  // given |KeplerianElements| around the barycentre of the existing bodies.
  // |osculating_elements_relative_to_system| must be a valid argument to the
  // constructor of |KeplerOrbit|.
  void Add(
      MassiveBody const& body,
      KeplerianElements<Frame> const& osculating_elements_relative_to_system);

  // A body with the total mass of the existing bodies.
  MassiveBody System() const;

  // Returns the degrees of freedom of the bodies with respect to their
  // barycentre, in the order in which they were added (starting with the
  // primary).
  std::vector<RelativeDegreesOfFreedom<Frame>> BarycentricDegreesOfFreedom()
      const;

 private:
  // A reference frame parallel to |Frame|, in which the primary is motionless
  // at the origin.
  using PrimocentricFrame = geometry::Frame<enum class PrimocentricFrameTag>;
  static Identity<PrimocentricFrame, Frame> const id_pf_;
  static Identity<Frame, PrimocentricFrame> const id_fp_;

  // The degrees of freedom of the bodies with respect to the primary,
  // in the order in which they were added.
  std::vector<DegreesOfFreedom<PrimocentricFrame>> primocentric_dof_;

  // The barycentre of the system, weighted by its total gravitational
  // parameter.
  BarycentreCalculator<DegreesOfFreedom<PrimocentricFrame>,
                       GravitationalParameter> system_barycentre_;
};

}  // namespace internal_jacobi_coordinates

using internal_jacobi_coordinates::JacobiCoordinates;

}  // namespace physics
}  // namespace principia

#include "physics/jacobi_coordinates_body.hpp"
