
#pragma once

#include <experimental/optional>
#include <list>
#include <map>
#include <memory>

#include "base/container_iterator.hpp"
#include "base/disjoint_sets.hpp"
#include "ksp_plugin/frames.hpp"
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

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
using PartId = std::uint32_t;

// Represents a KSP part.
class Part final {
 public:
  Part(PartId part_id, Mass const& mass);

  virtual PartId part_id() const;

  // Sets or returns the mass.  Event though a part is massless in the sense
  // that it doesn't exert gravity, it has a mass used to determine its
  // intrinsic acceleration.
  virtual void set_mass(Mass const& mass);
  virtual Mass const& mass() const;

  // Clears, increments or returns the intrinsic force exerted on the part by
  // its engines (or a tractor beam).
  // TODO(phl): Keep track of the point where the force is applied.
  virtual void clear_intrinsic_force();
  virtual void increment_intrinsic_force(
      Vector<Force, Barycentric> const& intrinsic_force);
  virtual Vector<Force, Barycentric> const& intrinsic_force() const;

  // Clears, sets or returns the degrees of freedom of the part.  A part that is
  // not in the bubble (i.e., that belongs to a pile-up that is not in the
  // bubble) doesn't have degrees of freedom.
  virtual void clear_degrees_of_freedom();
  virtual void set_degrees_of_freedom(
      DegreesOfFreedom<Bubble> const& degrees_of_freedom);
  virtual std::experimental::optional<DegreesOfFreedom<Bubble>> const&
  degrees_of_freedom() const;

  // This temporarily holds the trajectory followed by the part during the call
  // to |PileUp::AdvanceTime| for the containing |PileUp|.  It read and cleared
  // by |Vessel::AdvanceTime| for the containing |Vessel|.
  virtual DiscreteTrajectory<Barycentric>& tail();
  virtual DiscreteTrajectory<Barycentric> const& tail() const;

  // True if and only if the last point of the tail is authoritative, i.e.,
  // corresponds to a point in the psychohistory of the enclosing vessel.
  virtual bool tail_is_authoritative() const;
  virtual void set_tail_is_authoritative(bool tail_is_authoritative);

  // Requires |!is_piled_up()|.
  virtual void set_containing_pile_up(IteratorOn<std::list<PileUp>> pile_up);

  // An iterator to the containing |PileUp|, if any.  Do not |Erase| this
  // iterator, use |clear_pile_up| instead, which will take care of letting all
  // parts know that their |PileUp| is gone.
  virtual std::experimental::optional<IteratorOn<std::list<PileUp>>>
  containing_pile_up() const;

  // Whether this part is in a |PileUp|.  Equivalent to |containing_pile_up()|.
  virtual bool is_piled_up() const;

  // If this part is in a |PileUp|, erases that |PileUp|.  After this call, all
  // parts in that |PileUp| are no longer piled up.
  virtual void clear_pile_up();

  void WriteToMessage(not_null<serialization::Part*> message) const;
  static not_null<std::unique_ptr<Part>> ReadFromMessage(
      serialization::Part const& message);

 private:
  PartId const part_id_;
  Mass mass_;
  Vector<Force, Barycentric> intrinsic_force_;

  std::experimental::optional<IteratorOn<std::list<PileUp>>>
      containing_pile_up_;

  std::experimental::optional<DegreesOfFreedom<Bubble>> degrees_of_freedom_;
  not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> tail_;
  bool tail_is_authoritative_ = false;

  // TODO(egg): we may want to keep track of the moment of inertia, angular
  // momentum, etc.

  // We will use union-find algorithms on |Part|s.
  not_null<std::unique_ptr<Subset<Part>::Node>> const subset_node_;
  friend class Subset<Part>::Node;
};

std::ostream& operator<<(std::ostream& out, Part const& part);

using PartIdToOwnedPart = std::map<PartId, not_null<std::unique_ptr<Part>>>;
using IdAndOwnedPart = PartIdToOwnedPart::value_type;

}  // namespace internal_part

using internal_part::IdAndOwnedPart;
using internal_part::Part;
using internal_part::PartId;
using internal_part::PartIdToOwnedPart;

}  // namespace ksp_plugin

namespace base {

template<>
inline not_null<Subset<ksp_plugin::Part>::Node*>
Subset<ksp_plugin::Part>::Node::Get(ksp_plugin::Part& element) {
  return element.subset_node_.get();
}

}  // namespace base
}  // namespace principia
