
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel.hpp"

#include <list>

namespace principia {
namespace ksp_plugin {

using geometry::BarycentreCalculator;
using geometry::Position;
using physics::DegreesOfFreedom;

namespace internal_pile_up {

PileUp::PileUp(std::list<not_null<Vessel*>>&& vessels)
    : vessels_(std::move(vessels)) {
  BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> barycentre;
  for (not_null<Vessel*> vessel : vessels_) {
    barycentre.Add(vessel->prolongation().last().degrees_of_freedom(),
                   vessel->mass());
  }
  history_.Append(vessels_.front()->prolongation().last().time(),
                  barycentre.Get());
  last_point_is_authoritative_ = true;
}

void PileUp::set_mass_and_intrinsic_force(
    Mass const& mass,
    Vector<Force, Barycentric> const& intrinsic_force) {
  mass_ = mass;
  intrinsic_force_ = intrinsic_force;
}

std::list<not_null<Vessel*>> const& PileUp::vessels() const {
  return vessels_;
}

void PileUp::AdvanceTime(Ephemeris<Barycentric>& ephemeris,
                         Instant const& t) {
  if (!last_point_is_authoritative_) {
    auto const penultimate = --history_.last();
    history_.ForgetAfter(penultimate.time());
    last_point_is_authoritative_ = true;
  }

  if (intrinsic_force_ == Vector<Force, Barycentric>{}) {
    ephemeris.FlowWithFixedStep(
        {&history_},
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        t,
        fixed_step_parameters_);
    if (history_.last().time() < t) {
      // TODO(egg): this is clumsy, we need an option for FlowWithAdaptiveStep
      // to only add the last point.  Amusingly, this is the unwanted behaviour
      // of FlowWithFixedStep (#228).
      DiscreteTrajectory<Barycentric> prolongation;
      prolongation.Append(history_.last().time(),
                          history_.last().degrees_of_freedom());
      ephemeris.FlowWithAdaptiveStep(
          &history_,
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          t,
          adaptive_step_parameters_,
          Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
      history_.Append(prolongation.last().time(),
                      prolongation.last().degrees_of_freedom());
      last_point_is_authoritative_ = false;
    }
  } else {
    auto intrinsic_acceleration = [this](Instant const& t) {
      return intrinsic_force_ / mass_;
    };
    ephemeris.FlowWithAdaptiveStep(
        &history_,
        intrinsic_acceleration,
        t,
        adaptive_step_parameters_,
        Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
  }
  // TODO(egg): somehow append to the histories of the vessels.
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
