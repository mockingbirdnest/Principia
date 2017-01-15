
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel.hpp"

#include <list>

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using base::FindOrDie;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::OrthogonalMap;
using geometry::Position;
using physics::DegreesOfFreedom;
using physics::RigidMotion;
using physics::RigidTransformation;

PileUp::PileUp(std::list<not_null<Vessel*>>&& vessels)
    : vessels_(std::move(vessels)) {
  BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> barycentre;
  Vector<Force, Barycentric> total_intrinsic_force;
  for (not_null<Vessel*> const vessel : vessels_) {
    CHECK(vessel->psychohistory_is_authoritative());
    total_intrinsic_force += vessel->intrinsic_force();
    barycentre.Add(vessel->psychohistory().last().degrees_of_freedom(),
                   vessel->mass());
  }
  mass_ = barycentre.weight();
  intrinsic_force_ = total_intrinsic_force;

  psychohistory_.Append(vessels_.front()->psychohistory().last().time(),
                        barycentre.Get());
  psychohistory_is_authoritative_ = true;
}

void PileUp::set_mass(Mass const& mass) {
  mass_ = mass;
}

void PileUp::set_intrinsic_force(
    Vector<Force, Barycentric> const& intrinsic_force) {
  intrinsic_force_ = intrinsic_force;
}

std::list<not_null<Vessel*>> const& PileUp::vessels() const {
  return vessels_;
}

void PileUp::AdvanceTime(
    Ephemeris<Barycentric>& ephemeris,
    Instant const& t,
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  if (!psychohistory_is_authoritative_) {
    CHECK_GT(psychohistory_.Size(), 1);
    auto const penultimate = --psychohistory_.last();
    psychohistory_.ForgetAfter(penultimate.time());
    psychohistory_is_authoritative_ = true;
  }
  auto const last_preexisting_authoritative_point = psychohistory_.last();

  if (intrinsic_force_ == Vector<Force, Barycentric>{}) {
    ephemeris.FlowWithFixedStep(
        {&psychohistory_},
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        t,
        fixed_step_parameters);
    if (psychohistory_.last().time() < t) {
      Ephemeris<Barycentric>::AdaptiveStepParameters
          psychohistory_adaptive_step_parameters = adaptive_step_parameters;
      psychohistory_adaptive_step_parameters.set_last_point_only(true);
      ephemeris.FlowWithAdaptiveStep(
          &psychohistory_,
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          t,
          psychohistory_adaptive_step_parameters,
          Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
      psychohistory_is_authoritative_ = false;
    }
  } else {
    auto const a = intrinsic_force_ / mass_;
    auto const intrinsic_acceleration = [a](Instant const& t) { return a; };
    ephemeris.FlowWithAdaptiveStep(
        &psychohistory_,
        intrinsic_acceleration,
        t,
        adaptive_step_parameters,
        Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
  }
  auto it = last_preexisting_authoritative_point;
  ++it;
  for (; it != psychohistory_.End(); ++it) {
    auto const& pile_up_dof = it.degrees_of_freedom();
    RigidMotion<Barycentric, RigidPileUp> const from_barycentric(
        RigidTransformation<Barycentric, RigidPileUp>(
            pile_up_dof.position(),
            RigidPileUp::origin,
            OrthogonalMap<Barycentric, RigidPileUp>::Identity()),
        AngularVelocity<Barycentric>{},
        pile_up_dof.velocity());
    auto const to_barycentric = from_barycentric.Inverse();
    bool const authoritative =
        psychohistory_is_authoritative_ || it != psychohistory_.last();
    for (not_null<Vessel*> const vessel : vessels_)  {
      vessel->AppendToPsychohistory(
          it.time(),
          to_barycentric(FindOrDie(vessel_degrees_of_freedom_, vessel)),
          authoritative);
    }
  }
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
