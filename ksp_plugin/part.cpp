
#pragma once

#include "part.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

Part::Part(PartId part_id, Mass const& mass)
    : part_id_(part_id),
      mass_(mass) {}

PartId Part::part_id() const {
  return part_id_;
}

void Part::set_mass(Mass const& mass) {
  mass_ = mass;
}

Mass const& Part::mass() const {
  return mass_;
}

void Part::clear_intrinsic_force() {
  intrinsic_force_ = Vector<Force, Barycentric>{};
}

void Part::increment_intrinsic_force(
    Vector<Force, Barycentric> const& intrinsic_force) {
  intrinsic_force_ += intrinsic_force;
}

Vector<Force, Barycentric> const& Part::intrinsic_force() const {
  return intrinsic_force_;
}

void Part::clear_degrees_of_freedom() {
  degrees_of_freedom_ = std::experimental::nullopt;
}

void Part::set_degrees_of_freedom(
    DegreesOfFreedom<Bubble> const& degrees_of_freedom) {
  degrees_of_freedom_ = degrees_of_freedom;
}

std::experimental::optional<DegreesOfFreedom<Bubble>> const&
Part::degrees_of_freedom() {
  return degrees_of_freedom_;
}

DiscreteTrajectory<Barycentric>& Part::trajectory() {
  return trajectory_;
}

DiscreteTrajectory<Barycentric> const& Part::trajectory() const {
  return trajectory_;
}

void Part::set_containing_pile_up(IteratorOn<std::list<PileUp>> pile_up) {
  CHECK(!is_piled_up());
  containing_pile_up_ = pile_up;
}

std::experimental::optional<IteratorOn<std::list<PileUp>>>
Part::containing_pile_up() const {
  return containing_pile_up_;
}

bool Part::is_piled_up() const {
  // TODO(egg): |has_value()| once we have a standard |optional|.
  return static_cast<bool>(containing_pile_up_);
}

void Part::clear_pile_up() {
  if (is_piled_up()) {
    IteratorOn<std::list<PileUp>> pile_up = *containing_pile_up_;
    for (not_null<Part*> const part : pile_up.iterator()->parts()) {
      part->containing_pile_up_ = std::experimental::nullopt;
      part->psychohistory_is_authoritative_ = true;
    }
    CHECK(!is_piled_up());
    pile_up.Erase();
  }
}

void Part::WriteToMessage(not_null<serialization::Part*> const message) const {
  //TODO(phl): fix.
  //degrees_of_freedom_.WriteToMessage(message->mutable_degrees_of_freedom());
  //mass_.WriteToMessage(message->mutable_mass());
  //gravitational_acceleration_to_be_applied_by_ksp_.WriteToMessage(
  //    message->mutable_gravitational_acceleration_to_be_applied_by_ksp());
}

Part Part::ReadFromMessage(serialization::Part const& message) {
  //TODO(phl): fix.
  //return Part(DegreesOfFreedom<Frame>::ReadFromMessage(
  //                message.degrees_of_freedom()),
  //            Mass::ReadFromMessage(message.mass()),
  //            Vector<Acceleration, Frame>::ReadFromMessage(
  //                message.gravitational_acceleration_to_be_applied_by_ksp()));
}

std::ostream& operator<<(std::ostream& out, Part const& part) {
  return out << "{"
             << part.part_id() << ", "
             << part.mass() << "}";
}

}  // namespace internal_part
}  // namespace ksp_plugin
}  // namespace principia
