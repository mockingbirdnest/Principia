#pragma once

#include "ksp_plugin/vessel.hpp"

#include <list>

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using geometry::Vector;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using quantities::Force;

// A |PileUp| handles a connected component of the graph of |Vessels| under
// physical contact.  It advances the history and prolongation of its component
// |Vessels|, modeling them as a massless body at their centre of mass.
class PileUp {
 public:
  PileUp(std::list<not_null<Vessel*>> vessels);

  /*
  void AdvanceTimeInFreeFall(Instant const& t);
  void AdvanceTimeUnderThrust(Instant const& t,
                              Vector<Force, Barycentric> net_force);
                              */
  std::list<not_null<Vessel*>> const& vessels() const;

 private:
  std::list<not_null<Vessel*>> vessels_;

  /*
  MasslessBody const body_;
  not_null<Ephemeris<Barycentric>*> const ephemeris_;
  DiscreteTrajectory<Barycentric> history_;
  not_null<DiscreteTrajectory<Barycentric>*> prolongation_;
  */
};


}  // namespace internal_pile_up

using internal_pile_up::PileUp;

}  // namespace ksp_plugin
}  // namespace principia
