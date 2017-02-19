
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

PileUp::PileUp(std::list<not_null<Part*>>&& parts,
               DegreesOfFreedom<Barycentric> const& bubble_barycentre,
               Instant const& t)
    : parts_(std::move(parts)) {
  BarycentreCalculator<DegreesOfFreedom<Bubble>, Mass> barycentre;
  Vector<Force, Barycentric> total_intrinsic_force;
  for (not_null<Part*> const part : parts_) {
    total_intrinsic_force += part->intrinsic_force();
    barycentre.Add(*part->degrees_of_freedom(), part->mass());
  }
  mass_ = barycentre.weight();
  intrinsic_force_ = total_intrinsic_force;

  RigidMotion<Bubble, Barycentric> const bubble_to_barycentric =
      RigidMotion<Barycentric, Bubble>{
          RigidTransformation<Barycentric, Bubble>{
              bubble_barycentre.position(),
              Bubble::origin,
              OrthogonalMap<Barycentric, Bubble>::Identity()},
          AngularVelocity<Barycentric>{},
          bubble_barycentre.velocity()}.Inverse();
  psychohistory_.Append(t, bubble_to_barycentric(barycentre.Get()));
  psychohistory_is_authoritative_ = true;

  RigidMotion<Bubble, RigidPileUp> bubble_to_pile_up{
      RigidTransformation<Bubble, RigidPileUp>{
          barycentre.Get().position(),
          RigidPileUp::origin,
          OrthogonalMap<Bubble, RigidPileUp>::Identity()},
      AngularVelocity<Bubble>{},
      barycentre.Get().velocity()};
  for (not_null<Part*> const part : parts_) {
    actual_part_degrees_of_freedom_.emplace(
        part,
        bubble_to_pile_up(*part->degrees_of_freedom()));
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
  for (auto it = parts_.cbegin(); it != parts_.cend(); ++it) {
    CHECK(apparent_part_degrees_of_freedom_.find(*it) !=
          apparent_part_degrees_of_freedom_.cend());
  }

  // Compute the apparent centre of mass of the parts.
  BarycentreCalculator<DegreesOfFreedom<ApparentBubble>, Mass> calculator;
  for (auto it = apparent_part_degrees_of_freedom_.cbegin();
       it != apparent_part_degrees_of_freedom_.cend();
       ++it) {
    auto const part = it->first;
    auto const apparent_part_degrees_of_freedom = it->second;
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
  for (auto it = apparent_part_degrees_of_freedom_.cbegin();
       it != apparent_part_degrees_of_freedom_.cend();
       ++it) {
    auto const part = it->first;
    auto const apparent_part_degrees_of_freedom = it->second;
    actual_part_degrees_of_freedom_.emplace(part,
                                       apparent_bubble_to_pile_up_motion(
                                           apparent_part_degrees_of_freedom));
  }
  apparent_part_degrees_of_freedom_.clear();
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
    for (not_null<Part*> const part : parts_)  {
      part->AppendToPsychohistory(
          it.time(),
          to_barycentric(FindOrDie(actual_part_degrees_of_freedom_, part)),
          authoritative);
    }
  }
}

void PileUp::NudgeParts(
    DegreesOfFreedom<Barycentric> const& bubble_barycentre) const {
  auto const actual_centre_of_mass =
      psychohistory_.last().degrees_of_freedom();

  RigidMotion<Barycentric, RigidPileUp> const barycentric_to_pile_up{
      RigidTransformation<Barycentric, RigidPileUp>{
          actual_centre_of_mass.position(),
          RigidPileUp::origin,
          Identity<Barycentric, RigidPileUp>().Forget()},
      AngularVelocity<Barycentric>(),
      actual_centre_of_mass.velocity()};
  RigidMotion<Barycentric, Bubble> const barycentric_to_bubble{
      RigidTransformation<Barycentric, Bubble>{
          bubble_barycentre.position(),
          Bubble::origin,
          Identity<Barycentric, Bubble>().Forget()},
      AngularVelocity<Barycentric>(),
      bubble_barycentre.velocity()};

  auto const it = actual_part_degrees_of_freedom_.find(part);
  CHECK(it != actual_part_degrees_of_freedom_.cend());
  return (barycentric_to_bubble * barycentric_to_pile_up.Inverse())(it->second);
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
