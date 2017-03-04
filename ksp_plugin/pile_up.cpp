
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
using base::make_not_null_unique;
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
    : parts_(std::move(parts)),
      psychohistory_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {
  BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
  Vector<Force, Barycentric> total_intrinsic_force;
  for (not_null<Part*> const part : parts_) {
    total_intrinsic_force += part->intrinsic_force();
    calculator.Add(part->degrees_of_freedom(), part->mass());
  }
  mass_ = calculator.weight();
  intrinsic_force_ = total_intrinsic_force;
  DegreesOfFreedom<Barycentric> const barycentre = calculator.Get();
  psychohistory_->Append(t, barycentre);

  RigidMotion<Barycentric, RigidPileUp> const barycentric_to_pile_up{
      RigidTransformation<Barycentric, RigidPileUp>{
          barycentre.position(),
          RigidPileUp::origin,
          Identity<Barycentric, RigidPileUp>().Forget()},
      AngularVelocity<Barycentric>{},
      barycentre.velocity()};
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
  CHECK_EQ(psychohistory_->Size(), 1);

  DiscreteTrajectory<Barycentric> prolongation;

  if (intrinsic_force_ == Vector<Force, Barycentric>{}) {
    ephemeris.FlowWithFixedStep(
        {psychohistory_.get()},
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        t,
        fixed_step_parameters);
    if (psychohistory_->last().time() < t) {
      // TODO(phl): Consider using forks.
      prolongation.Append(psychohistory_->last().time(),
                          psychohistory_->last().degrees_of_freedom());
      CHECK(ephemeris.FlowWithAdaptiveStep(
                &prolongation,
                Ephemeris<Barycentric>::NoIntrinsicAcceleration,
                t,
                adaptive_step_parameters,
                Ephemeris<Barycentric>::unlimited_max_ephemeris_steps,
                /*last_point_only=*/true));
      CHECK_EQ(prolongation.Size(), 2);
    }
  } else {
    auto const a = intrinsic_force_ / mass_;
    auto const intrinsic_acceleration = [a](Instant const& t) { return a; };
    CHECK(ephemeris.FlowWithAdaptiveStep(
              psychohistory_.get(),
              intrinsic_acceleration,
              t,
              adaptive_step_parameters,
              Ephemeris<Barycentric>::unlimited_max_ephemeris_steps,
              /*last_point_only=*/false));
  }
  auto it = psychohistory_->Begin();
  ++it;
  for (; it != psychohistory_->End(); ++it) {
    AppendToPartTails(it, /*authoritative=*/true);
  }
  // TODO(phl): you have reinvented the map.  Now I want an empty :-p
  if (prolongation.Size() != 0) {
    CHECK_EQ(prolongation.Size(), 2);
    AppendToPartTails(prolongation.last(), /*authoritative=*/false);
  }
  psychohistory_->ForgetBefore(psychohistory_->last().time());
}

void PileUp::NudgeParts() const {
  auto const actual_centre_of_mass =
      psychohistory_->last().degrees_of_freedom();

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

void PileUp::WriteToMessage(not_null<serialization::PileUp*> message) const {
  for (auto const part : parts_) {
    message->add_part_id(part->part_id());
  }
  mass_.WriteToMessage(message->mutable_mass());
  intrinsic_force_.WriteToMessage(message->mutable_intrinsic_force());
  psychohistory_->WriteToMessage(message->mutable_psychohistory(),
                                 /*forks=*/{});
  for (auto const& pair : actual_part_degrees_of_freedom_) {
    auto const part = pair.first;
    auto const& degrees_of_freedom = pair.second;
    degrees_of_freedom.WriteToMessage(&(
        (*message->mutable_actual_part_degrees_of_freedom())[part->part_id()]));
  }
  for (auto const& pair : apparent_part_degrees_of_freedom_) {
    auto const part = pair.first;
    auto const& degrees_of_freedom = pair.second;
    degrees_of_freedom.WriteToMessage(&(
        (*message
              ->mutable_apparent_part_degrees_of_freedom())[part->part_id()]));
  }
}

not_null<std::unique_ptr<PileUp>> PileUp::ReadFromMessage(
    serialization::PileUp const& message,
    std::function<not_null<Part*>(PartId)> const& part_id_to_part) {
  std::list<not_null<Part*>> parts;
  for (auto const part_id : message.part_id()) {
    parts.push_back(part_id_to_part(part_id));
  }
  not_null<std::unique_ptr<PileUp>> pile_up =
      std::unique_ptr<PileUp>(new PileUp(std::move(parts)));
  pile_up->mass_ = Mass::ReadFromMessage(message.mass());
  pile_up->intrinsic_force_ =
      Vector<Force, Barycentric>::ReadFromMessage(message.intrinsic_force());
  pile_up->psychohistory_ =
      DiscreteTrajectory<Barycentric>::ReadFromMessage(message.psychohistory(),
                                                       /*forks=*/{});
  for (auto const& pair : message.actual_part_degrees_of_freedom()) {
    std::uint32_t const part_id = pair.first;
    serialization::Pair const& degrees_of_freedom = pair.second;
    pile_up->actual_part_degrees_of_freedom_.emplace(
        part_id_to_part(part_id),
        DegreesOfFreedom<RigidPileUp>::ReadFromMessage(degrees_of_freedom));
  }
  for (auto const& pair : message.apparent_part_degrees_of_freedom()) {
    std::uint32_t const part_id = pair.first;
    serialization::Pair const& degrees_of_freedom = pair.second;
    pile_up->apparent_part_degrees_of_freedom_.emplace(
        part_id_to_part(part_id),
        DegreesOfFreedom<ApparentBubble>::ReadFromMessage(degrees_of_freedom));
  }
  return pile_up;
}

PileUp::PileUp(std::list<not_null<Part*>>&& parts)
    : parts_(std::move(parts)),
      psychohistory_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {}

void PileUp::AppendToPartTails(
    DiscreteTrajectory<Barycentric>::Iterator const it,
    bool const authoritative) const {
  auto const& pile_up_dof = it.degrees_of_freedom();
  RigidMotion<Barycentric, RigidPileUp> const barycentric_to_pile_up(
      RigidTransformation<Barycentric, RigidPileUp>(
          pile_up_dof.position(),
          RigidPileUp::origin,
          Identity<Barycentric, RigidPileUp>().Forget()),
      AngularVelocity<Barycentric>{},
      pile_up_dof.velocity());
  auto const pile_up_to_barycentric = barycentric_to_pile_up.Inverse();
  for (not_null<Part*> const part : parts_) {
    part->tail().Append(it.time(),
                        pile_up_to_barycentric(
                            FindOrDie(actual_part_degrees_of_freedom_, part)));
    part->set_tail_is_authoritative(authoritative);
  }
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
