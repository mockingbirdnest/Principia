#include "ksp_plugin/part.hpp"

#include <memory>
#include <string>
#include <utility>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/space_transformations.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace _part {
namespace internal {

using namespace principia::base::_array;
using namespace principia::base::_hexadecimal;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_space_transformations;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_si;

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

DiscreteTrajectory<Barycentric>::iterator Part::history_begin() {
  return history_->begin();
}

DiscreteTrajectory<Barycentric>::iterator Part::history_end() {
  return history_->end();
}

DiscreteTrajectory<Barycentric>::iterator Part::psychohistory_begin() {
  if (psychohistory_ == trajectory_.segments().end()) {
    psychohistory_ = trajectory_.NewSegment();
  }
  // TODO(phl): This used to say:
  // Make sure that we skip the fork, which may be the point of the prehistory.
  return psychohistory_->begin();
}

DiscreteTrajectory<Barycentric>::iterator Part::psychohistory_end() {
  if (psychohistory_ == trajectory_.segments().end()) {
    psychohistory_ = trajectory_.NewSegment();
  }
  return psychohistory_->end();
}

void Part::AppendToHistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  if (psychohistory_ != trajectory_.segments().end()) {
    trajectory_.DeleteSegments(psychohistory_);
  }
  trajectory_.Append(time, degrees_of_freedom).IgnoreError();
}

void Part::AppendToPsychohistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  if (psychohistory_ == trajectory_.segments().end()) {
    psychohistory_ = trajectory_.NewSegment();
  }
  trajectory_.Append(time, degrees_of_freedom).IgnoreError();
}

void Part::ClearHistory() {
  trajectory_.clear();
  psychohistory_ = trajectory_.segments().end();
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
  trajectory_.WriteToMessage(message->mutable_prehistory(),
                             /*tracked=*/{history_, psychohistory_},
                             /*exact=*/{});
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
  bool const is_pre_hamilton = message.prehistory().segment_size() == 0;
  LOG_IF(WARNING, is_pre_hamilton)
      << "Reading pre-"
      << (is_pre_cesàro    ? "Cesàro"
          : is_pre_fréchet ? "Fréchet"
          : is_pre_frenet  ? "Frenet"
          : is_pre_galileo ? "Galileo"
                           : "Hamilton") << " Part";

  std::unique_ptr<Part> part;
  if (is_pre_fréchet) {
    auto const degrees_of_freedom =
        DegreesOfFreedom<Barycentric>::ReadFromMessage(
            message.degrees_of_freedom());
    part = std::unique_ptr<Part>(new Part(
        message.part_id(),
        message.name(),
        message.truthful(),
        Mass::ReadFromMessage(message.mass()),
        MakeWaterSphereInertiaTensor(Mass::ReadFromMessage(message.mass())),
        RigidMotion<RigidPart, Barycentric>::MakeNonRotatingMotion(
            degrees_of_freedom),
        std::move(deletion_callback)));
  } else if (is_pre_frenet) {
    part = std::unique_ptr<Part>(new Part(
        message.part_id(),
        message.name(),
        message.truthful(),
        Mass::ReadFromMessage(message.pre_frenet_inertia_tensor().mass()),
        InertiaTensor<RigidPart>::ReadFromMessage(
            message.pre_frenet_inertia_tensor().form()),
        RigidMotion<RigidPart, Barycentric>::ReadFromMessage(
            message.rigid_motion()),
        std::move(deletion_callback)));
  } else {
    part = std::unique_ptr<Part>(new Part(
        message.part_id(),
        message.name(),
        message.truthful(),
        Mass::ReadFromMessage(message.mass()),
        InertiaTensor<RigidPart>::ReadFromMessage(message.inertia_tensor()),
        RigidMotion<RigidPart, Barycentric>::ReadFromMessage(
            message.rigid_motion()),
        std::move(deletion_callback)));
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
        /*tracked=*/{});
    // The `history_` has been created by the constructor above.  Construct the
    // various trajectories from the `tail`.
    for (auto it = tail.begin(); it != tail.end();) {
      auto const& [time, degrees_of_freedom] = *it;
      ++it;
      if (it == tail.end() && !message.tail_is_authoritative()) {
        part->AppendToPsychohistory(time, degrees_of_freedom);
      } else {
        part->AppendToHistory(time, degrees_of_freedom);
      }
    }
  } else if (is_pre_hamilton) {
    DiscreteTrajectorySegmentIterator<Barycentric> history;
    DiscreteTrajectorySegmentIterator<Barycentric> psychohistory;
    auto prehistory = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.prehistory(),
        /*tracked=*/{&history, &psychohistory});
    // We want to get rid of the prehistory segment, so the easiest is to detach
    // the history and psychohistory.  We must take care of the segment
    // iterators as they get invalidated in this case.
    part->trajectory_ = prehistory.DetachSegments(history);
    part->history_ = part->trajectory_.segments().begin();
    part->psychohistory_ = std::next(part->history_);
  } else {
    part->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.prehistory(),
        /*tracked=*/{&part->history_, &part->psychohistory_});
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
      history_(trajectory_.segments().begin()),
      psychohistory_(trajectory_.segments().end()),
      subset_node_(make_not_null_unique<Subset<Part>::Node>()),
      deletion_callback_(std::move(deletion_callback)) {}

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

}  // namespace internal
}  // namespace _part
}  // namespace ksp_plugin
}  // namespace principia
