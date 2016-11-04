#pragma once

#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_coupling {

using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;

// A |Coupling| handles a connected component of the graph of |Vessels| under
// physical contact.  It advances the history and prolongation of its component
// |Vessels|, modeling them as a massless body at their centre of mass.
class Coupling {
 public:
 private:
  MasslessBody const body_;
  not_null<Ephemeris<Barycentric>*> const ephemeris_;
  std::unique_ptr<DiscreteTrajectory<Barycentric>> history_;
  DiscreteTrajectory<Barycentric>* prolongation_;
};

}  // namespace internal_coupling
}  // namespace ksp_plugin
}  // namespace principia
