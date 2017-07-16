
#pragma once

#include <experimental/optional>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <string>

#include "base/container_iterator.hpp"
#include "base/disjoint_sets.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/identification.hpp"
#include "ksp_plugin/part_subsets.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
namespace internal_part {

using base::IteratorOn;
using base::not_null;
using base::Subset;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using quantities::Force;
using quantities::Mass;

// Represents a KSP part.
class Part final {
 public:
  Part(PartId part_id,
       std::string const& name,
       Mass const& mass,
       DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
       std::function<void()> deletion_callback);

  // Calls the deletion callback passed at construction, if any.  This part must
  // not be piled up.
  ~Part();

  PartId part_id() const;

  // Sets or returns the mass.  Event though a part is massless in the sense
  // that it doesn't exert gravity, it has a mass used to determine its
  // intrinsic acceleration.
  void set_mass(Mass const& mass);
  Mass const& mass() const;

  // Clears, increments or returns the intrinsic force exerted on the part by
  // its engines (or a tractor beam).
  // TODO(phl): Keep track of the point where the force is applied.
  void clear_intrinsic_force();
  void increment_intrinsic_force(
      Vector<Force, Barycentric> const& intrinsic_force);
  Vector<Force, Barycentric> const& intrinsic_force() const;

  // Sets or returns the degrees of freedom of the part.
  void set_degrees_of_freedom(
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);
  DegreesOfFreedom<Barycentric> const& degrees_of_freedom() const;

  // This temporarily holds the trajectory followed by the part during the call
  // to |PileUp::AdvanceTime| for the containing |PileUp|.  It is read and
  // cleared by |Vessel::AdvanceTime| for the containing |Vessel|.
  DiscreteTrajectory<Barycentric>& tail();
  DiscreteTrajectory<Barycentric> const& tail() const;

  // True if and only if the last point of the tail is authoritative, i.e.,
  // corresponds to a point in the psychohistory of the enclosing Part.
  bool tail_is_authoritative() const;
  void set_tail_is_authoritative(bool tail_is_authoritative);

  // Requires |!is_piled_up()|.
  void set_containing_pile_up(IteratorOn<std::list<PileUp>> pile_up);

  // An iterator to the containing |PileUp|, if any.  Do not |Erase| this
  // iterator, use |clear_pile_up| instead, which will take care of letting all
  // parts know that their |PileUp| is gone.
  std::experimental::optional<IteratorOn<std::list<PileUp>>>
  containing_pile_up() const;

  // Whether this part is in a |PileUp|.  Equivalent to |containing_pile_up()|.
  bool is_piled_up() const;

  // If this part is in a |PileUp|, erases that |PileUp|.  After this call, all
  // parts in that |PileUp| are no longer piled up.
  void clear_pile_up();

  void WriteToMessage(not_null<serialization::Part*> message) const;
  static not_null<std::unique_ptr<Part>> ReadFromMessage(
      serialization::Part const& message,
      std::function<void()> deletion_callback);
  void FillContainingPileUpFromMessage(
      serialization::Part const& message,
      not_null<std::list<PileUp>*> const pile_ups);

  // Returns "part name (part ID)".
  std::string ShortDebugString() const;

 private:
  PartId const part_id_;
  std::string const name_;
  Mass mass_;
  Vector<Force, Barycentric> intrinsic_force_;

  std::experimental::optional<IteratorOn<std::list<PileUp>>>
      containing_pile_up_;

  DegreesOfFreedom<Barycentric> degrees_of_freedom_;
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> tail_;
  bool tail_is_authoritative_ = false;

  // TODO(egg): we may want to keep track of the moment of inertia, angular
  // momentum, etc.

  // We will use union-find algorithms on |Part|s.
  not_null<std::unique_ptr<Subset<Part>::Node>> const subset_node_;
  friend class Subset<Part>::Node;

  // Called in the destructor.
  std::function<void()> deletion_callback_;
};

std::ostream& operator<<(std::ostream& out, Part const& part);

}  // namespace internal_part

using internal_part::Part;

}  // namespace ksp_plugin

namespace base {

template<>
inline not_null<Subset<ksp_plugin::Part>::Node*>
Subset<ksp_plugin::Part>::Node::Get(ksp_plugin::Part& element) {
  return element.subset_node_.get();
}

}  // namespace base
}  // namespace principia
