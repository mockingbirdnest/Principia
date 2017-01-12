
#pragma once

#include <list>

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
using geometry::Instant;
using geometry::Vector;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using quantities::Force;
using quantities::Mass;

// A |PileUp| handles a connected component of the graph of |Vessels| under
// physical contact.  It advances the history and prolongation of its component
// |Vessels|, modeling them as a massless body at their centre of mass.
class PileUp final {
 public:
  explicit PileUp(std::list<not_null<Vessel*>>&& vessels);

  void set_mass_and_intrinsic_force(
      Mass const& mass,
      Vector<Force, Barycentric> const& intrinsic_force);

  std::list<not_null<Vessel*>> const& vessels() const;

  // Flows the history authoritatively as far as possible up to |t|, advances
  // the histories of the vessels.  After this call, the histories of |*this|
  // and of its vessels have a (possibly non-authoritative) final point exactly
  // at |t|.
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
  DiscreteTrajectory<Barycentric> history_;
  // True if the last point should be flowed from; otherwise, the last point
  // should be removed before flowing the trajectory (in that case, there is a
  // penultimate point, and it is authoritative).
  bool last_point_is_authoritative_;
};

}  // namespace internal_pile_up

using internal_pile_up::PileUp;

}  // namespace ksp_plugin
}  // namespace principia
