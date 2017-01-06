
#pragma once

#include <list>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"

namespace principia {
namespace ksp_plugin {

FORWARD_DECLARE_FROM(vessel, class, Vessel);

namespace internal_pile_up {

using base::not_null;
using geometry::Vector;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using quantities::Force;

// A |PileUp| handles a connected component of the graph of |Vessels| under
// physical contact.  It advances the history and prolongation of its component
// |Vessels|, modeling them as a massless body at their centre of mass.
class PileUp final {
 public:
  explicit PileUp(std::list<not_null<Vessel*>>&& vessels);

  std::list<not_null<Vessel*>> const& vessels() const;

 private:
  std::list<not_null<Vessel*>> vessels_;
};

}  // namespace internal_pile_up

using internal_pile_up::PileUp;

}  // namespace ksp_plugin
}  // namespace principia
