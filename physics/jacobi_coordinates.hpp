#pragma once

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/identity.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"

namespace principia {

using base::not_null;
using geometry::Identity;

namespace physics {

// An utility for converting a linearly ordered system of massive bodies given
// in Jacobi coordinates to barycentric coordinates.
template<typename Frame>
class JacobiCoordinates {
 public:
  JacobiCoordinates(MassiveBody const& primary);

  // Adds |body| with the given |DegreesOfFreedom| with respect to the
  // barycentre of the existing bodies.
  void Add(MassiveBody const& body,
           RelativeDegreesOfFreedom<Frame> const& dof_wrt_system);

  // Adds |body| with the |RelativeDegreesOfFreedom| of a |KeplerOrbit| with the
  // given |KeplerianElements| around the barycentre of the existing bodies.
  // |osculating_elements_wrt_system| must be a valid argument to the
  // constructor of |KeplerOrbit|.
  // Equivalent to
  //   Add(body,
  //       KeplerOrbit<Frame>(/*primary=*/System(),
  //                          /*secondary=*/body,
  //                          osculating_elements_wrt_system,
  //                          epoch).StateVectors(epoch));
  // for any |epoch|.
  void Add(
      MassiveBody const& body,
      KeplerianElements<Frame> const& osculating_elements_wrt_system);

  // A body with the total mass of the existing bodies.
  MassiveBody System() const;

  // Returns the degrees of freedom of the bodies with respect to their
  // barycentre, in the order in which they were added (starting with the
  // primary).
  std::vector<RelativeDegreesOfFreedom<Frame>> BarycentricCoordinates() const;

 private:
  // A reference frame parallel to |Frame|, in which the primary is motionless
  // at the origin.
  enum class PrivateFrameTag { kFrameTag };
  using PrimocentricFrame = geometry::Frame<PrivateFrameTag,
                                            PrivateFrameTag::kFrameTag,
                                            /*frame_is_inertial=*/false>;
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

template<typename Frame>
class HierarchicalSystem {
 public:
  struct BarycentricSystem {
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    std::vector<DegreesOfFreedom<Frame>> barycentric_degrees_of_freedom;
  };

  HierarchicalSystem(not_null<std::unique_ptr<MassiveBody const>> primary);

  // Adds the given |body| with the given |parent|.  |parent| must already have
  // been inserted.  |jacobi_osculating_elements| must be a valid argument to
  // the constructor of |KeplerOrbit|.
  void Add(not_null<std::unique_ptr<MassiveBody const>> body,
           not_null<MassiveBody const*> const parent,
           KeplerianElements<Frame> const& jacobi_osculating_elements);

  // Puts the barycentre of the system at the motionless origin of |Frame|;
  // |*this| is invalid after a call to |Get()|.
  BarycentricSystem Get();

 private:
  struct Subsystem;
  struct System {
    System(not_null<std::unique_ptr<MassiveBody const>> p)
        : primary(std::move(p)){};
    virtual ~System() = default;
    not_null<std::unique_ptr<MassiveBody const>> primary;
    std::vector<not_null<std::unique_ptr<Subsystem>>> satellites;
  };
  struct Subsystem : public System {
    Subsystem(not_null<std::unique_ptr<MassiveBody const>> p)
        : System(std::move(p)){};
    KeplerianElements<Frame> jacobi_osculating_elements;
  };

  System system_;
  // Invariant: |subsystems_[p]->primary.get() == p|.
  // None of these pointers should be null, but I want to use operator[].
  std::map<not_null<MassiveBody const*>, System*> subsystems_;
};

}  // namespace physics
}  // namespace principia

#include "physics/jacobi_coordinates_body.hpp"
