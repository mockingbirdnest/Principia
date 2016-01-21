#pragma once

#include <map>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/identity.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/jacobi_coordinates.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"

namespace principia {

using base::not_null;

namespace physics {

template<typename Frame>
class HierarchicalSystem {
 public:
  struct BarycentricSystem {
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
    std::vector<DegreesOfFreedom<Frame>> degrees_of_freedom;
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
  std::map<not_null<MassiveBody const*>, System*> systems_;
};

}  // namespace physics
}  // namespace principia

#include "physics/hierarchical_system_body.hpp"
