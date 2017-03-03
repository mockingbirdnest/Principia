
#pragma once

#include "ksp_plugin/part.hpp"

#include "base/not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using base::make_not_null_unique;

Part::Part(
    PartId const part_id,
    Mass const& mass,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
    std::function<void()> deletion_callback)
    : part_id_(part_id),
      mass_(mass),
      degrees_of_freedom_(degrees_of_freedom),
      tail_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      subset_node_(make_not_null_unique<Subset<Part>::Node>()),
      deletion_callback_(std::move(deletion_callback)) {}

Part::~Part() {
  clear_pile_up();
  if (deletion_callback_ != nullptr) {
    deletion_callback_();
  }
}

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

void Part::set_degrees_of_freedom(
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  degrees_of_freedom_ = degrees_of_freedom;
}

DegreesOfFreedom<Barycentric> const&
Part::degrees_of_freedom() const {
  return degrees_of_freedom_;
}

DiscreteTrajectory<Barycentric>& Part::tail() {
  return *tail_;
}

DiscreteTrajectory<Barycentric> const& Part::tail() const {
  return *tail_;
}

bool Part::tail_is_authoritative() const {
  return tail_is_authoritative_;
}

void Part::set_tail_is_authoritative(bool const tail_is_authoritative) {
  tail_is_authoritative_ = tail_is_authoritative;
}

void Part::set_containing_pile_up(IteratorOn<std::list<PileUp>> const pile_up) {
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
    }
    CHECK(!is_piled_up());
    pile_up.Erase();
  }
}

void Part::WriteToMessage(not_null<serialization::Part*> const message) const {
  message->set_part_id(part_id_);
  mass_.WriteToMessage(message->mutable_mass());
  intrinsic_force_.WriteToMessage(message->mutable_intrinsic_force());
  if (containing_pile_up_) {
    // TODO(phl): Implement.
  }
  degrees_of_freedom_.WriteToMessage(message->mutable_degrees_of_freedom());
  tail_->WriteToMessage(message->mutable_tail(), /*forks=*/{});
  message->set_tail_is_authoritative(tail_is_authoritative_);
}

not_null<std::unique_ptr<Part>> Part::ReadFromMessage(
    serialization::Part const& message,
    std::function<void()> deletion_callback) {
  not_null<std::unique_ptr<Part>> part =
      make_not_null_unique<Part>(message.part_id(),
                                 DegreesOfFreedom<Barycentric>::ReadFromMessage(
                                     message.degrees_of_freedom()),
                                 Mass::ReadFromMessage(message.mass()),
                                 std::move(deletion_callback));
  part->increment_intrinsic_force(
      Vector<Force, Barycentric>::ReadFromMessage(message.intrinsic_force()));
  if (message.has_containing_pile_up()) {
    // TODO(phl): Implement.
  }
  part->tail_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(message.tail(),
                                                                 /*forks=*/{});
  part->set_tail_is_authoritative(message.tail_is_authoritative());
  return part;
}

std::ostream& operator<<(std::ostream& out, Part const& part) {
  return out << "{"
             << part.part_id() << ", "
             << part.mass() << "}";
}

}  // namespace internal_part
}  // namespace ksp_plugin
}  // namespace principia
