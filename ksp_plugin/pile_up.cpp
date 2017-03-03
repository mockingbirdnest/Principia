
#include "ksp_plugin/pile_up.hpp"

#include <list>
#include <map>

#include "geometry/identity.hpp"
#include "ksp_plugin/part.hpp"
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

PileUp::PileUp(std::list<not_null<Part*>>&& parts, Instant const& t)
    : parts_(std::move(parts)) {
  BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
  Vector<Force, Barycentric> total_intrinsic_force;
  for (not_null<Part*> const part : parts_) {
    total_intrinsic_force += part->intrinsic_force();
    calculator.Add(part->degrees_of_freedom(), part->mass());
  }
  mass_ = calculator.weight();
  intrinsic_force_ = total_intrinsic_force;
  psychohistory_.Append(t, calculator.Get());

  RigidMotion<Barycentric, RigidPileUp> const barycentric_to_pile_up{
      RigidTransformation<Barycentric, RigidPileUp>{
          calculator.Get().position(),
          RigidPileUp::origin,
          OrthogonalMap<Barycentric, RigidPileUp>::Identity()},
      AngularVelocity<Barycentric>{},
      calculator.Get().velocity()};
  for (not_null<Part*> const part : parts_) {
    actual_part_degrees_of_freedom_.emplace(
        part,
        barycentric_to_pile_up(part->degrees_of_freedom()));
  }
}

void PileUp::set_mass(Mass const& mass) {
  mass_ = mass;
}

void PileUp::set_intrinsic_force(
    Vector<Force, Barycentric> const& intrinsic_force) {
  intrinsic_force_ = intrinsic_force;
}

std::list<not_null<Part*>> const& PileUp::parts() const {
  return parts_;
}

void PileUp::SetPartApparentDegreesOfFreedom(
    not_null<Part*> const part,
    DegreesOfFreedom<ApparentBubble> const& degrees_of_freedom) {
  std::map<not_null<Part*>, DegreesOfFreedom<ApparentBubble>>::iterator it;
  bool inserted;
  std::tie(it, inserted) =
      apparent_part_degrees_of_freedom_.emplace(part, degrees_of_freedom);
  CHECK(inserted) << "Duplicate part " << part << " at "
                  << degrees_of_freedom;
}

void PileUp::DeformPileUpIfNeeded() {
  if (apparent_part_degrees_of_freedom_.empty()) {
    return;
  }
  // A consistency check that |SetPartApparentDegreesOfFreedom| was called for
  // all the parts.
  CHECK_EQ(parts_.size(), apparent_part_degrees_of_freedom_.size());
  for (not_null<Part*> const part : parts_) {
    CHECK_GT(apparent_part_degrees_of_freedom_.count(part), 0);
  }

  // Compute the apparent centre of mass of the parts.
  BarycentreCalculator<DegreesOfFreedom<ApparentBubble>, Mass> calculator;
  for (auto const& pair : apparent_part_degrees_of_freedom_) {
    auto const part = pair.first;
    auto const& apparent_part_degrees_of_freedom = pair.second;
    calculator.Add(apparent_part_degrees_of_freedom, part->mass());
  }
  auto const apparent_centre_of_mass = calculator.Get();

  // A motion that maps the apparent centre of mass of the parts to the actual
  // centre of mass of the pile-up.
  RigidTransformation<ApparentBubble, RigidPileUp> const
      apparent_bubble_to_pile_up_transformation(
          apparent_centre_of_mass.position(),
          RigidPileUp::origin,
          Identity<ApparentBubble, RigidPileUp>().Forget());
  RigidMotion<ApparentBubble, RigidPileUp> const
      apparent_bubble_to_pile_up_motion(
          apparent_bubble_to_pile_up_transformation,
          AngularVelocity<ApparentBubble>(),
          apparent_centre_of_mass.velocity());

  // Now update the positions of the parts in the pile-up frame.
  actual_part_degrees_of_freedom_.clear();
  for (auto const& pair : apparent_part_degrees_of_freedom_) {
    auto const part = pair.first;
    auto const& apparent_part_degrees_of_freedom = pair.second;
    actual_part_degrees_of_freedom_.emplace(
        part,
        apparent_bubble_to_pile_up_motion(apparent_part_degrees_of_freedom));
  }
  apparent_part_degrees_of_freedom_.clear();
}

void PileUp::AdvanceTime(
    Ephemeris<Barycentric>& ephemeris,
    Instant const& t,
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  CHECK_EQ(psychohistory_.Size(), 1);

  DiscreteTrajectory<Barycentric> prolongation;

  if (intrinsic_force_ == Vector<Force, Barycentric>{}) {
    ephemeris.FlowWithFixedStep(
        {&psychohistory_},
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        t,
        fixed_step_parameters);
    if (psychohistory_.last().time() < t) {
      prolongation.Append(prolongation.last().time(),
                          prolongation.last().degrees_of_freedom());
      ephemeris.FlowWithAdaptiveStep(
          &prolongation,
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          t,
          adaptive_step_parameters,
          Ephemeris<Barycentric>::unlimited_max_ephemeris_steps,
          /*last_point_only=*/true);
      CHECK_EQ(prolongation.Size(), 2);
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
  auto it = psychohistory_.Begin();
  ++it;
  for (; it != psychohistory_.End(); ++it) {
    AppendToPartTails(it, /*authoritative=*/true);
  }
  // TODO(phl): you have reinvented the map.  Now I want an empty :-p
  if (prolongation.Size() != 0) {
    CHECK_EQ(prolongation.Size(), 2);
    AppendToPartTails(prolongation.last(), /*authoritative=*/false);
  }
  psychohistory_.ForgetBefore(psychohistory_.last().time());
}

void PileUp::NudgeParts() const {
  auto const actual_centre_of_mass =
      psychohistory_.last().degrees_of_freedom();

  RigidMotion<Barycentric, RigidPileUp> const barycentric_to_pile_up{
      RigidTransformation<Barycentric, RigidPileUp>{
          actual_centre_of_mass.position(),
          RigidPileUp::origin,
          Identity<Barycentric, RigidPileUp>().Forget()},
      AngularVelocity<Barycentric>(),
      actual_centre_of_mass.velocity()};
  auto const pile_up_to_barycentric = barycentric_to_pile_up.Inverse();
  for (not_null<Part*> const part : parts_) {
    part->set_degrees_of_freedom(pile_up_to_barycentric(
        FindOrDie(actual_part_degrees_of_freedom_, part)));
  }
}

void PileUp::AppendToPartTails(
    DiscreteTrajectory<Barycentric>::Iterator const it,
    bool const authoritative) const {
  auto const& pile_up_dof = it.degrees_of_freedom();
  RigidMotion<Barycentric, RigidPileUp> const from_barycentric(
      RigidTransformation<Barycentric, RigidPileUp>(
          pile_up_dof.position(),
          RigidPileUp::origin,
          OrthogonalMap<Barycentric, RigidPileUp>::Identity()),
      AngularVelocity<Barycentric>{},
      pile_up_dof.velocity());
  auto const to_barycentric = from_barycentric.Inverse();
  for (not_null<Part*> const part : parts_) {
    part->tail().Append(
        it.time(),
        to_barycentric(FindOrDie(actual_part_degrees_of_freedom_, part)));
    part->set_tail_is_authoritative(authoritative);
  }
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
