#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <string>

#include "base/disjoint_sets.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/part_subsets.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace _part {
namespace internal {

using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::DiscreteTrajectorySegmentIterator;
using physics::RigidMotion;
using namespace principia::base::_disjoint_sets;
using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_named_quantities;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Represents a KSP part.
class Part final {
 public:
  // A truthful part.
  Part(PartId part_id,
       std::string const& name,
       Mass const& mass,
       Position<EccentricPart> const& centre_of_mass,
       InertiaTensor<RigidPart> const& inertia_tensor,
       RigidMotion<EccentricPart, Barycentric> const& rigid_motion,
       std::function<void()> deletion_callback);

  // An untruthful part.
  Part(PartId part_id,
       std::string const& name,
       DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
       std::function<void()> deletion_callback);

  // Calls the deletion callback passed at construction, if any.  This part must
  // not be piled up.
  ~Part();

  PartId part_id() const;

  // When a part is not truthful, all its properties except for name and part_id
  // are lies and should not be propagated to the game.
  bool truthful() const;
  void make_truthful();

  // Sets or returns the mass and inertia tensor.  Even though a part is
  // massless in the sense that it doesn't exert gravity, it has a mass and an
  // inertia used to determine its intrinsic acceleration and rotational
  // properties.
  void set_mass(Mass const& mass);
  Mass const& mass() const;
  void set_centre_of_mass(Position<EccentricPart> const& centre_of_mass);
  // Returns the transformation from the part-centre-of-mass-centred to the
  // part-position-centred frame; the offset is set by |set_centre_of_mass|.
  RigidMotion<RigidPart, EccentricPart> MakeRigidToEccentricMotion() const;
  void set_inertia_tensor(InertiaTensor<RigidPart> const& inertia_tensor);
  InertiaTensor<RigidPart> const& inertia_tensor() const;
  // Whether this part is a solid rocket motor, whose lost mass is expelled with
  // its angular momentum.
  void set_is_solid_rocket_motor(bool is_solid_rocket_motor);
  bool is_solid_rocket_motor() const;

  // The difference between successive values passed to |set_mass()|.
  Mass const& mass_change() const;

  // Clears, increments or returns the intrinsic force exerted on the part by
  // its engines (or a tractor beam).
  void clear_intrinsic_force();
  void apply_intrinsic_force(
      Vector<Force, Barycentric> const& intrinsic_force);
  Vector<Force, Barycentric> const& intrinsic_force() const;

  void clear_intrinsic_torque();
  void apply_intrinsic_torque(
      Bivector<Torque, Barycentric> const& intrinsic_torque);
  Bivector<Torque, Barycentric> const& intrinsic_torque() const;

  void ApplyIntrinsicForceWithLeverArm(
      Vector<Force, Barycentric> const& intrinsic_force,
      Displacement<Barycentric> const& lever_arm);

  // Sets or returns the rigid motion of the part.
  void set_rigid_motion(
      RigidMotion<RigidPart, Barycentric> const& rigid_motion);
  RigidMotion<RigidPart, Barycentric> const& rigid_motion() const;

  // Return iterators to the beginning and end of the history and psychohistory
  // of the part, respectively.  Either trajectory may be empty, but they are
  // not both empty.
  DiscreteTrajectory<Barycentric>::iterator history_begin();
  DiscreteTrajectory<Barycentric>::iterator history_end();
  DiscreteTrajectory<Barycentric>::iterator psychohistory_begin();
  DiscreteTrajectory<Barycentric>::iterator psychohistory_end();

  // Appends a point to the history or psychohistory of this part.  These
  // temporarily hold the trajectory of the part and are constructed by
  // |PileUp::AdvanceTime|.  They are consumed by |Vessel::AdvanceTime| for the
  // containing |Vessel|.
  // Note that |AppendToHistory| clears the psychohistory so the order of the
  // calls matter.
  void AppendToHistory(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);
  void AppendToPsychohistory(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Clears the history and psychohistory.
  void ClearHistory();

  // Requires |!is_piled_up()|.  The part assumes co-ownership of the |pile_up|.
  void set_containing_pile_up(not_null<std::shared_ptr<PileUp>> const& pile_up);

  // A pointer to the containing pile up, if any.
  PileUp* containing_pile_up() const;

  // Whether this part is in a |PileUp|, i.e., has a non-null
  // |containing_pile_up|.
  bool is_piled_up() const;

  // Remove this part from its pile-up, if any.  This may cause the pile-up to
  // be destroyed if this was the last part owning the pile-up.
  void reset_containing_pile_up();

  void WriteToMessage(not_null<serialization::Part*> message,
                      PileUp::SerializationIndexForPileUp const&
                          serialization_index_for_pile_up) const;
  static not_null<std::unique_ptr<Part>> ReadFromMessage(
      serialization::Part const& message,
      std::function<void()> deletion_callback);
  void FillContainingPileUpFromMessage(
      serialization::Part const& message,
      PileUp::PileUpForSerializationIndex const&
          pile_up_for_serialization_index);

  // Returns "part name (part ID)".
  std::string ShortDebugString() const;

 private:
  Part(PartId part_id,
       std::string name,
       bool truthful,
       Mass const& mass,
       InertiaTensor<RigidPart> const& inertia_tensor,
       RigidMotion<RigidPart, Barycentric> rigid_motion,
       std::function<void()> deletion_callback);

  static RigidMotion<RigidPart, EccentricPart> MakeRigidToEccentricMotion(
      Position<EccentricPart> const& centre_of_mass);

  PartId const part_id_;
  std::string const name_;
  bool truthful_;
  Mass mass_;
  Position<EccentricPart> centre_of_mass_ = EccentricPart::origin;
  // NOTE(eggrobin): |mass_change_| and |is_solid_rocket_motor_| are set by
  // |InsertOrKeepLoadedPart|, and used by |PileUp::RecomputeFromParts|.
  // Ultimately, both are called in the adapter in |WaitedForFixedUpdate|.
  // They therefore do not need to be serialized.
  Mass mass_change_;
  bool is_solid_rocket_motor_ = false;
  InertiaTensor<RigidPart> inertia_tensor_;
  Vector<Force, Barycentric> intrinsic_force_;
  Bivector<Torque, Barycentric> intrinsic_torque_;

  std::shared_ptr<PileUp> containing_pile_up_;

  RigidMotion<RigidPart, Barycentric> rigid_motion_;

  // See the comments in pile_up.hpp for an explanation of the terminology.

  // The trajectory of the part, composed of (at most) two segments, the
  // history and the psychohistory.
  DiscreteTrajectory<Barycentric> trajectory_;

  // The |history_| is nearly always present, except in some transient
  // situations.
  DiscreteTrajectorySegmentIterator<Barycentric> history_;

  // The |psychohistory_| is destroyed by |AppendToHistory| and is recreated
  // as needed by |AppendToPsychohistory|.
  DiscreteTrajectorySegmentIterator<Barycentric> psychohistory_;

  // We will use union-find algorithms on |Part|s.
  not_null<std::unique_ptr<Subset<Part>::Node>> const subset_node_;
  friend class Subset<Part>::Node;

  // Called in the destructor.
  std::function<void()> deletion_callback_;
};

// A factory that creates an inertia tensor for a solid sphere of water having
// the given mass.  Useful, e.g., for save compatibility.
InertiaTensor<RigidPart> MakeWaterSphereInertiaTensor(Mass const& mass);

std::ostream& operator<<(std::ostream& out, Part const& part);

}  // namespace internal

using internal::Part;
using internal::MakeWaterSphereInertiaTensor;

}  // namespace _part
}  // namespace ksp_plugin

namespace base {

template<>
inline not_null<Subset<ksp_plugin::_part::Part>::Node*>
Subset<ksp_plugin::_part::Part>::Node::Get(ksp_plugin::_part::Part& element) {
  return element.subset_node_.get();
}

}  // namespace base
}  // namespace principia

namespace principia::ksp_plugin {
using namespace principia::ksp_plugin::_part;
}  // namespace principia::ksp_plugin
