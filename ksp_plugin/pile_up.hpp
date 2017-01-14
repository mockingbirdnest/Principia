
#pragma once

#include <list>
#include <map>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "ksp_plugin/frames.hpp"

namespace principia {
namespace ksp_plugin {

FORWARD_DECLARE_FROM(vessel, class, Vessel);

namespace internal_pile_up {

using base::not_null;
using geometry::Frame;
using geometry::Instant;
using geometry::Vector;
using physics::DiscreteTrajectory;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::MasslessBody;
using physics::RelativeDegreesOfFreedom;
using quantities::Force;
using quantities::Mass;

// A |PileUp| handles a connected component of the graph of |Vessels| under
// physical contact.  It advances the history and prolongation of its component
// |Vessels|, modeling them as a massless body at their centre of mass.
class PileUp final {
 public:
  explicit PileUp(std::list<not_null<Vessel*>>&& vessels);

  void set_mass(Mass const& mass);
  void set_intrinsic_force(Vector<Force, Barycentric> const& intrinsic_force);

  std::list<not_null<Vessel*>> const& vessels() const;

  // Flows the history authoritatively as far as possible up to |t|, advances
  // the histories of the vessels.  After this call, the histories of |*this|
  // and of its vessels have a (possibly ahistorical) final point exactly at
  // |t|.
  void AdvanceTime(
      Ephemeris<Barycentric>& ephemeris,
      Instant const& t,
      Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_step_parameters);

 private:
  std::list<not_null<Vessel*>> vessels_;
  Mass mass_;
  Vector<Force, Barycentric> intrinsic_force_;
  DiscreteTrajectory<Barycentric> psychohistory_;
  // True if the last point of the |psychohistory_| should be flowed from;
  // otherwise, the last point should be removed before flowing the trajectory
  // (in that case, there is a penultimate point, and it is historical).
  bool psychohistory_is_authoritative_;
  // The |PileUp| is seen as a (currently non-rotating) rigid body; the degrees
  // of freedom of the vessels in the frame of that body can be set, however
  // their motion is not integrated; this is simply applied as an offset from
  // the rigid body motion of the |PileUp|.
  // The origin of |RigidPileUp| is the centre of mass of the pile up.
  // its axes are those of Barycentric for now; eventually we will probably want
  // to use the inertia ellipsoid.
  using RigidPileUp = Frame<serialization::Frame::PluginTag,
                            serialization::Frame::RIGID_PILE_UP,
                            /*frame_is_inertial=*/false>;
  std::map<not_null<Vessel*>, DegreesOfFreedom<RigidPileUp>>
      vessel_degrees_of_freedom_;
};

}  // namespace internal_pile_up

using internal_pile_up::PileUp;

}  // namespace ksp_plugin
}  // namespace principia
