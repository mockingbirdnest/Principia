
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
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using base::not_null;
using base::Subset;
using geometry::Bivector;
using geometry::InertiaTensor;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::RigidMotion;
using quantities::Force;
using quantities::Mass;
using quantities::Torque;

// Represents a KSP part.
class Part final {
 public:
  Part(PartId part_id,
       std::string const& name,
       Mass const& mass,
       InertiaTensor<RigidPart> const& inertia_tensor,
       RigidMotion<RigidPart, Barycentric> const& rigid_motion,
       std::function<void()> deletion_callback);

  // Calls the deletion callback passed at construction, if any.  This part must
  // not be piled up.
  ~Part();

  PartId part_id() const;

  // Sets or returns the mass and inertia tensor.  Even though a part is
  // massless in the sense that it doesn't exert gravity, it has a mass and an
  // inertia used to determine its intrinsic acceleration and rotational
  // properties.
  void set_mass(Mass const& mass);
  Mass const& mass() const;
  void set_inertia_tensor(InertiaTensor<RigidPart> const& inertia_tensor);
  InertiaTensor<RigidPart> const& inertia_tensor() const;

  // Clears, increments or returns the intrinsic force exerted on the part by
  // its engines (or a tractor beam).
  // TODO(phl): Keep track of the point where the force is applied.
  void clear_intrinsic_force();
  void increment_intrinsic_force(
      Vector<Force, Barycentric> const& intrinsic_force);
  Vector<Force, Barycentric> const& intrinsic_force() const;

  void clear_intrinsic_torque();
  void increment_intrinsic_torque(
      Bivector<Torque, Barycentric> const& intrinsic_torque);
  Bivector<Torque, Barycentric> const& intrinsic_torque() const;

  // Sets or returns the rigid motion of the part.
  void set_rigid_motion(
      RigidMotion<RigidPart, Barycentric> const& rigid_motion);
  RigidMotion<RigidPart, Barycentric> const& rigid_motion() const;

  // A convenience selector.
  // TODO(phl): Should probably be eliminated at some point.
  DegreesOfFreedom<Barycentric> degrees_of_freedom() const;

  // Return iterators to the beginning and end of the history and psychohistory
  // of the part, respectively.  Either trajectory may be empty, but they are
  // not both empty.
  DiscreteTrajectory<Barycentric>::Iterator history_begin();
  DiscreteTrajectory<Barycentric>::Iterator history_end();
  DiscreteTrajectory<Barycentric>::Iterator psychohistory_begin();
  DiscreteTrajectory<Barycentric>::Iterator psychohistory_end();

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
  PartId const part_id_;
  std::string const name_;
  Mass mass_;
  InertiaTensor<RigidPart> inertia_tensor_;
  Vector<Force, Barycentric> intrinsic_force_;
  Bivector<Torque, Barycentric> intrinsic_torque_;

  std::shared_ptr<PileUp> containing_pile_up_;

  RigidMotion<RigidPart, Barycentric> rigid_motion_;

  // See the comments in pile_up.hpp for an explanation of the terminology.

  // The |prehistory_| always has a single point at time -∞.  It sole purpose is
  // to make it convenient to hook the |psychohistory_| even if there is no
  // point in the |history_| (it's not possible to fork-at-last an empty root
  // trajectory, but it works for a non-root).
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> prehistory_;
  // The |history_| is nearly always not null, except in some transient
  // situations.  It's a fork of the |prehistory_|.
  DiscreteTrajectory<Barycentric>* history_ = nullptr;
  // The |psychohistory_| is destroyed by |AppendToHistory| and is recreated
  // as needed by |AppendToPsychohistory| or by |tail|.  That's because
  // |NewForkAtLast| is relatively expensive so we only call it when necessary.
  DiscreteTrajectory<Barycentric>* psychohistory_ = nullptr;

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

}  // namespace internal_part

using internal_part::Part;
using internal_part::MakeWaterSphereInertiaTensor;

}  // namespace ksp_plugin

namespace base {

template<>
inline not_null<Subset<ksp_plugin::Part>::Node*>
Subset<ksp_plugin::Part>::Node::Get(ksp_plugin::Part& element) {
  return element.subset_node_.get();
}

}  // namespace base
}  // namespace principia
