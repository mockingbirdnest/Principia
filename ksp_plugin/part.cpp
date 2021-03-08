﻿
#include "ksp_plugin/part.hpp"

#include <list>
#include <string>
#include <utility>

#include "base/allocated_by.hpp"
#include "base/allocator_new.hpp"
#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "base/not_null.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using base::AllocateWith;
using base::Array;
using base::AssertAllocatedBy;
using base::HexadecimalEncoder;
using base::make_not_null_unique;
using base::UniqueArray;
using geometry::OrthogonalMap;
using geometry::R3x3Matrix;
using physics::RigidTransformation;
using quantities::Cbrt;
using quantities::Density;
using quantities::MomentOfInertia;
using quantities::Pow;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;

constexpr Mass untruthful_part_mass = 1 * Kilogram;

Part::Part(PartId const part_id,
           std::string const& name,
           Mass const& mass,
           Position<EccentricPart> const& centre_of_mass,
           InertiaTensor<RigidPart> const& inertia_tensor,
           RigidMotion<EccentricPart, Barycentric> const& rigid_motion,
           std::function<void()> deletion_callback)
    : Part(part_id,
           name,
           /*truthful=*/true,
           mass,
           inertia_tensor,
           rigid_motion * MakeRigidToEccentricMotion(centre_of_mass),
           std::move(deletion_callback)) {}

Part::Part(PartId part_id,
           std::string const& name,
           DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
           std::function<void()> deletion_callback)
    : Part(part_id,
           name,
           /*truthful=*/false,
           untruthful_part_mass,
           MakeWaterSphereInertiaTensor(untruthful_part_mass),
           RigidMotion<RigidPart, Barycentric>::MakeNonRotatingMotion(
               degrees_of_freedom),
           std::move(deletion_callback)) {}

Part::~Part() {
  LOG(INFO) << "Destroying part " << ShortDebugString();
  if (deletion_callback_ != nullptr) {
    deletion_callback_();
  }
}

PartId Part::part_id() const {
  return part_id_;
}

bool Part::truthful() const {
  return truthful_;
}

void Part::make_truthful() {
  truthful_ = true;
}

void Part::set_mass(Mass const& mass) {
  mass_change_ = mass - mass_;
  mass_ = mass;
}

Mass const& Part::mass() const {
  return mass_;
}

void Part::set_centre_of_mass(Position<EccentricPart> const& centre_of_mass) {
  centre_of_mass_ = centre_of_mass;
}

RigidMotion<RigidPart, EccentricPart> Part::MakeRigidToEccentricMotion() const {
  return MakeRigidToEccentricMotion(centre_of_mass_);
}

void Part::set_inertia_tensor(InertiaTensor<RigidPart> const& inertia_tensor) {
  inertia_tensor_ = inertia_tensor;
}

InertiaTensor<RigidPart> const& Part::inertia_tensor() const {
  return inertia_tensor_;
}

void Part::set_is_solid_rocket_motor(bool const is_solid_rocket_motor) {
  is_solid_rocket_motor_ = is_solid_rocket_motor;
}

bool Part::is_solid_rocket_motor() const {
  return is_solid_rocket_motor_;
}

Mass const& Part::mass_change() const {
  return mass_change_;
}

void Part::clear_intrinsic_force() {
  intrinsic_force_ = Vector<Force, Barycentric>{};
}

void Part::apply_intrinsic_force(
    Vector<Force, Barycentric> const& intrinsic_force) {
  intrinsic_force_ += intrinsic_force;
}

Vector<Force, Barycentric> const& Part::intrinsic_force() const {
  return intrinsic_force_;
}

void Part::clear_intrinsic_torque() {
  intrinsic_torque_ = Bivector<Torque, Barycentric>{};
}

void Part::apply_intrinsic_torque(
    Bivector<Torque, Barycentric> const& intrinsic_torque) {
  intrinsic_torque_ += intrinsic_torque;
}

Bivector<Torque, Barycentric> const& Part::intrinsic_torque() const {
  return intrinsic_torque_;
}

void Part::ApplyIntrinsicForceWithLeverArm(
    Vector<Force, Barycentric> const& intrinsic_force,
    Displacement<Barycentric> const& lever_arm) {
  apply_intrinsic_force(intrinsic_force);
  apply_intrinsic_torque(Wedge(lever_arm, intrinsic_force) * Radian);
}

void Part::set_rigid_motion(
    RigidMotion<RigidPart, Barycentric> const& rigid_motion) {
  rigid_motion_ = rigid_motion;
}

RigidMotion<RigidPart, Barycentric> const& Part::rigid_motion() const {
  return rigid_motion_;
}

DiscreteTrajectory<Barycentric>::Iterator Part::history_begin() {
  // Make sure that we skip the point of the prehistory.
  auto it = history_->Fork();
  return ++it;
}

DiscreteTrajectory<Barycentric>::Iterator Part::history_end() {
  return history_->end();
}

DiscreteTrajectory<Barycentric>::Iterator Part::psychohistory_begin() {
  if (psychohistory_ == nullptr) {
    psychohistory_ = history_->NewForkAtLast();
  }
  // Make sure that we skip the fork, which may be the point of the prehistory.
  auto it = psychohistory_->Fork();
  return ++it;
}

DiscreteTrajectory<Barycentric>::Iterator Part::psychohistory_end() {
  if (psychohistory_ == nullptr) {
    psychohistory_ = history_->NewForkAtLast();
  }
  return psychohistory_->end();
}

void Part::AppendToHistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  if (psychohistory_ != nullptr) {
    history_->DeleteFork(psychohistory_);
  }
  history_->Append(time, degrees_of_freedom);
}

void Part::AppendToPsychohistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  if (psychohistory_ == nullptr) {
    psychohistory_ = history_->NewForkAtLast();
  }
  psychohistory_->Append(time, degrees_of_freedom);
}

void Part::ClearHistory() {
  if (psychohistory_ != nullptr) {
    history_->DeleteFork(psychohistory_);
  }
  prehistory_->DeleteFork(history_);
  history_ = prehistory_->NewForkAtLast();
}

void Part::set_containing_pile_up(
    not_null<std::shared_ptr<PileUp>> const& pile_up) {
  CHECK(!is_piled_up());
  LOG(INFO) << "Adding part " << ShortDebugString() << " to the pile up at "
            << pile_up;
  containing_pile_up_ = static_cast<std::shared_ptr<PileUp> const&>(pile_up);
}

PileUp* Part::containing_pile_up() const {
  return containing_pile_up_.get();
}

bool Part::is_piled_up() const {
  return containing_pile_up_ != nullptr;
}

void Part::reset_containing_pile_up() {
  LOG_IF(INFO, containing_pile_up_ != nullptr)
      << "Removing part " << ShortDebugString() << " from its pile up at "
      << containing_pile_up_;
  containing_pile_up_.reset();
}

void Part::WriteToMessage(not_null<serialization::Part*> const message,
                          PileUp::SerializationIndexForPileUp const&
                              serialization_index_for_pile_up) const {
  message->set_part_id(part_id_);
  message->set_name(name_);
  message->set_truthful(truthful_);
  mass_.WriteToMessage(message->mutable_mass());
  centre_of_mass_.WriteToMessage(message->mutable_centre_of_mass());
  inertia_tensor_.WriteToMessage(message->mutable_inertia_tensor());
  intrinsic_force_.WriteToMessage(message->mutable_intrinsic_force());
  intrinsic_torque_.WriteToMessage(message->mutable_intrinsic_torque());
  if (containing_pile_up_) {
    message->set_containing_pile_up(
        serialization_index_for_pile_up(containing_pile_up_.get()));
  }
  rigid_motion_.WriteToMessage(message->mutable_rigid_motion());
  prehistory_->WriteToMessage(message->mutable_prehistory(),
                              /*forks=*/{history_, psychohistory_});
}

not_null<std::unique_ptr<Part>> Part::ReadFromMessage(
    serialization::Part const& message,
    std::function<void()> deletion_callback) {
  bool const is_pre_cesàro = message.has_tail_is_authoritative();
  bool const is_pre_fréchet = message.has_mass() &&
                              message.has_degrees_of_freedom();
  bool const is_pre_frenet =
      is_pre_fréchet || (message.has_pre_frenet_inertia_tensor() &&
                         !message.has_intrinsic_torque());
  bool const is_pre_galileo = !message.has_centre_of_mass();

  std::unique_ptr<Part> part;
  if (is_pre_fréchet) {
    auto const degrees_of_freedom =
        DegreesOfFreedom<Barycentric>::ReadFromMessage(
            message.degrees_of_freedom());
    part = std::unique_ptr<Part>(AssertAllocatedBy<std::allocator<Part>>(
        new (AllocateWith<std::allocator<Part>>{}) Part(
            message.part_id(),
            message.name(),
            message.truthful(),
            Mass::ReadFromMessage(message.mass()),
            MakeWaterSphereInertiaTensor(Mass::ReadFromMessage(message.mass())),
            RigidMotion<RigidPart, Barycentric>::MakeNonRotatingMotion(
                degrees_of_freedom),
            std::move(deletion_callback))));
  } else if (is_pre_frenet) {
    part = std::unique_ptr<Part>(AssertAllocatedBy<std::allocator<Part>>(
        new (AllocateWith<std::allocator<Part>>{}) Part(
            message.part_id(),
            message.name(),
            message.truthful(),
            Mass::ReadFromMessage(message.pre_frenet_inertia_tensor().mass()),
            InertiaTensor<RigidPart>::ReadFromMessage(
                message.pre_frenet_inertia_tensor().form()),
            RigidMotion<RigidPart, Barycentric>::ReadFromMessage(
                message.rigid_motion()),
            std::move(deletion_callback))));
  } else {
    part = std::unique_ptr<Part>(AssertAllocatedBy<std::allocator<Part>>(
        new (AllocateWith<std::allocator<Part>>{}) Part(
            message.part_id(),
            message.name(),
            message.truthful(),
            Mass::ReadFromMessage(message.mass()),
            InertiaTensor<RigidPart>::ReadFromMessage(message.inertia_tensor()),
            RigidMotion<RigidPart, Barycentric>::ReadFromMessage(
                message.rigid_motion()),
            std::move(deletion_callback))));
  }

  if (!is_pre_galileo) {
    part->set_centre_of_mass(
        Position<EccentricPart>::ReadFromMessage(message.centre_of_mass()));
  }

  part->apply_intrinsic_force(
      Vector<Force, Barycentric>::ReadFromMessage(message.intrinsic_force()));
  if (!is_pre_frenet) {
    part->apply_intrinsic_torque(
        Bivector<Torque, Barycentric>::ReadFromMessage(
            message.intrinsic_torque()));
  }

  if (is_pre_cesàro) {
    auto tail = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.prehistory(),
        /*forks=*/{});
    // The |prehistory_| and |history_| have been created by the constructor
    // above.  Construct the various trajectories from the |tail|.
    for (auto it = tail->begin(); it != tail->end();) {
      auto const& [time, degrees_of_freedom] = *it;
      ++it;
      if (it == tail->end() && !message.tail_is_authoritative()) {
        part->AppendToPsychohistory(time, degrees_of_freedom);
      } else {
        part->AppendToHistory(time, degrees_of_freedom);
      }
    }
  } else {
    part->history_ = nullptr;
    part->prehistory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.prehistory(),
        /*forks=*/{&part->history_, &part->psychohistory_});
  }
  return std::move(part);
}

void Part::FillContainingPileUpFromMessage(
    serialization::Part const& message,
    PileUp::PileUpForSerializationIndex const&
        pile_up_for_serialization_index) {
  if (message.has_containing_pile_up()) {
    containing_pile_up_ =
        pile_up_for_serialization_index(message.containing_pile_up());
  }
}

std::string Part::ShortDebugString() const {
  Array<std::uint8_t const> id_bytes(
      reinterpret_cast<std::uint8_t const*>(&part_id_), sizeof(part_id_));
  HexadecimalEncoder</*null_terminated=*/true> encoder;
  auto const hex_id = encoder.Encode(id_bytes);
  return name_ + " (" + hex_id.data.get() + ")";
}

Part::Part(PartId const part_id,
           std::string name,
           bool const truthful,
           Mass const& mass,
           InertiaTensor<RigidPart> const& inertia_tensor,
           RigidMotion<RigidPart, Barycentric> rigid_motion,
           std::function<void()> deletion_callback)
    : part_id_(part_id),
      name_(std::move(name)),
      truthful_(truthful),
      mass_(mass),
      inertia_tensor_(inertia_tensor),
      rigid_motion_(std::move(rigid_motion)),
      prehistory_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      subset_node_(make_not_null_unique<Subset<Part>::Node>()),
      deletion_callback_(std::move(deletion_callback)) {
  prehistory_->Append(astronomy::InfinitePast,
                      {Barycentric::origin, Barycentric::unmoving});
  history_ = prehistory_->NewForkAtLast();
}

RigidMotion<RigidPart, EccentricPart> Part::MakeRigidToEccentricMotion(
    Position<EccentricPart> const& centre_of_mass) {
  return RigidMotion<RigidPart, EccentricPart>(
      RigidTransformation<RigidPart, EccentricPart>(
          RigidPart::origin,
          centre_of_mass,
          OrthogonalMap<RigidPart, EccentricPart>::Identity()),
      RigidPart::nonrotating,
      RigidPart::unmoving);
}

InertiaTensor<RigidPart> MakeWaterSphereInertiaTensor(Mass const& mass) {
  static constexpr MomentOfInertia zero;
  static constexpr Density ρ_of_water = 1000 * Kilogram / Pow<3>(Metre);
  MomentOfInertia const I =
      Cbrt(9 * Pow<5>(mass) / (250 * Pow<2>(π * ρ_of_water)));
  return InertiaTensor<RigidPart>(R3x3Matrix<MomentOfInertia>({I, zero, zero},
                                                              {zero, I, zero},
                                                              {zero, zero, I}));
}

std::ostream& operator<<(std::ostream& out, Part const& part) {
  return out << "{" << part.part_id() << ", " << part.mass() << "}";
}

}  // namespace internal_part
}  // namespace ksp_plugin
}  // namespace principia
