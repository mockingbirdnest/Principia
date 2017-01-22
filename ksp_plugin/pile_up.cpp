
#include "ksp_plugin/pile_up.hpp"

#include <list>

#include "geometry/identity.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/rigid_motion.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

using base::FindOrDie;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Identity;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Velocity;
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

void PileUp::SetVesselApparentDegreesOfFreedom(
    not_null<Vessel*> const vessel,
    DegreesOfFreedom<Barycentric> const degrees_of_freedom) {
  static RigidMotion<Barycentric, ApparentBarycentric> const
      barycentric_to_apparent_barycentric(
          RigidTransformation<Barycentric, ApparentBarycentric>::Identity(),
          AngularVelocity<Barycentric>(),
          Velocity<Barycentric>());

  std::map<not_null<Vessel*>,
           DegreesOfFreedom<ApparentBarycentric>>::iterator it;
  bool inserted;
  std::tie(it, inserted) = apparent_vessel_degrees_of_freedom_.emplace(
      vessel, barycentric_to_apparent_barycentric(degrees_of_freedom));
  CHECK(inserted) << "Duplicate vessel " << vessel << " at "
                  << degrees_of_freedom;
}

void PileUp::UpdateVesselsInPileUp() {
  // A consistency check that |SetVesselApparentDegreesOfFreedom| was called for
  // all the vessels.
  CHECK_EQ(vessels_.size(), apparent_vessel_degrees_of_freedom_.size());
  for (auto it = vessels_.cbegin(); it != vessels_.cend(); ++it) {
    CHECK(apparent_vessel_degrees_of_freedom_.find(*it) !=
          apparent_vessel_degrees_of_freedom_.cend());
  }

  // Compute the apparent centre of mass of the vessels.
  BarycentreCalculator<DegreesOfFreedom<ApparentBarycentric>, Mass> calculator;
  for (auto it = apparent_vessel_degrees_of_freedom_.cbegin();
       it != apparent_vessel_degrees_of_freedom_.cend();
       ++it) {
    auto const vessel = it->first;
    auto const apparent_vessel_degrees_of_freedom = it->second;
    calculator.Add(apparent_vessel_degrees_of_freedom, vessel->mass());
  }
  auto const apparent_centre_of_mass = calculator.Get();

  // A motion that maps the apparent centre of mass of the vessels to the actual
  // centre of mass of the pile-up.
  RigidTransformation<ApparentBarycentric, RigidPileUp> const
      apparent_barycentric_to_pile_up_transformation(
          apparent_centre_of_mass.position(),
          RigidPileUp::origin,
          Identity<ApparentBarycentric, RigidPileUp>().Forget());
  RigidMotion<ApparentBarycentric, RigidPileUp> const
      apparent_barycentric_to_pile_up_motion(
          apparent_barycentric_to_pile_up_transformation,
          AngularVelocity<ApparentBarycentric>(),
          apparent_centre_of_mass.velocity());

  // Now update the positions of the vessels in the pile-up frame.
  vessel_degrees_of_freedom_.clear();
  for (auto it = apparent_vessel_degrees_of_freedom_.cbegin();
       it != apparent_vessel_degrees_of_freedom_.cend();
       ++it) {
    auto const vessel = it->first;
    auto const apparent_vessel_degrees_of_freedom = it->second;
    vessel_degrees_of_freedom_.emplace(vessel,
                                       apparent_barycentric_to_pile_up_motion(
                                           apparent_vessel_degrees_of_freedom));
  }
  apparent_vessel_degrees_of_freedom_.clear();
}

DegreesOfFreedom<Barycentric> PileUp::GetVesselActualDegreesOfFreedom(
    not_null<Vessel*> const vessel) const {
  auto const actual_centre_of_mass =
      psychohistory_.last().degrees_of_freedom();
  RigidTransformation<Barycentric, RigidPileUp> const
      barycentric_to_pile_up_transformation(
          actual_centre_of_mass.position(),
          RigidPileUp::origin,
          Identity<Barycentric, RigidPileUp>().Forget());
  RigidMotion<Barycentric, RigidPileUp> const barycentric_to_pile_up_motion(
      barycentric_to_pile_up_transformation,
      AngularVelocity<Barycentric>(),
      actual_centre_of_mass.velocity());

  auto const it = vessel_degrees_of_freedom_.find(vessel);
  CHECK(it != vessel_degrees_of_freedom_.cend());
  return barycentric_to_pile_up_motion.Inverse()(it->second);
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
      ephemeris.FlowWithAdaptiveStep(
          &psychohistory_,
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          t,
          adaptive_step_parameters,
          Ephemeris<Barycentric>::unlimited_max_ephemeris_steps,
          /*last_point_only=*/true);
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
        Ephemeris<Barycentric>::unlimited_max_ephemeris_steps,
        /*last_point_only=*/false);
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
