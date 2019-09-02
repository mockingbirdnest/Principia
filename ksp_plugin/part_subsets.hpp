
#pragma once

#include <list>

#include "base/disjoint_sets.hpp"

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"

// This ksp_plugin file is in namespace |base| to specialize templates declared
// therein.

namespace principia {

namespace ksp_plugin {
FORWARD_DECLARE_FROM(part, class, Part);
}  // namespace ksp_plugin

namespace base {

// Within an union-find on |Part|s, we maintain lists of the elements in the
// disjoint sets.  Moreover, we keep track of the inclusion relations of those
// sets to the sets of |Part|s in existing |PileUp|s, destroying existing
// |PileUp|s as we learn that they will not appear in the new arrangement.
// The |Collect| operation finalizes this, destroying existing |PileUp| which
// are strict supersets of the new sets, and creating the new |PileUp|s.
template<>
class Subset<ksp_plugin::Part>::Properties final {
  using PileUps = std::list<ksp_plugin::PileUp*>;

 public:
  // |*part| must outlive the constructed object.
  explicit Properties(not_null<ksp_plugin::Part*> part);

  // If |*this| and |other| are subsets of different |PileUp|s, or one is a
  // subset and not the other, the relevant |PileUp|s are erased.
  // Otherwise, |this->subset_of_existing_pile_up_| keeps track of the number of
  // missing parts.
  // Maintains |parts_| by joining the lists.
  void MergeWith(Properties& other);

  // “What’s this thing suddenly coming towards me very fast? Very very fast.
  // So big and flat and round, it needs a big wide sounding name like … ow …
  // ound … round … ground! That’s it! That’s a good name – ground!  I wonder if
  // it will be friends with me?”
  void Ground();
  bool grounded() const;

  // If |collected_|, performs no action.
  // Otherwise, sets |collected_|, and:
  // - if |EqualsExistingPileUp()|, performs no action.
  // - if |StrictSubsetOfExistingPileUp()|, erases the existing |PileUp| inserts
  //   a new |PileUp| into |pile_ups| with the parts in |parts_|.
  // - if |!subset_of_existing_pile_up_|, inserts a new |PileUp| into |pile_ups|
  //   with the parts in |parts_|.  The new |PileUp| is created using the given
  //   parameters.
  void Collect(
      PileUps& pile_ups,
      geometry::Instant const& t,
      physics::Ephemeris<ksp_plugin::Barycentric>::AdaptiveStepParameters const&
          adaptive_step_parameters,
      physics::Ephemeris<ksp_plugin::Barycentric>::FixedStepParameters const&
          fixed_step_parameters,
      not_null<physics::Ephemeris<ksp_plugin::Barycentric>*> ephemeris);

 private:
  // Whether |left| and |right| are both subsets of the same existing |PileUp|.
  // Implies |left.SubsetOfExistingPileUp()| and
  // |right.SubsetOfExistingPileUp()|.
  static bool SubsetsOfSamePileUp(Properties const& left,
                                  Properties const& right);
  // Whether the set of |Part|s in |parts_| is equal to the set of |Part|s
  // in an existing |PileUp|.  Implies |SubsetOfExistingPileUp()|.
  bool EqualsExistingPileUp() const;
  // Whether the set of |Part|s in |parts_| is a subset of the set of
  // |Part|s in an existing |PileUp|.  In that case
  // |parts_.front()->containing_pile_up()| is that |PileUp|.
  bool SubsetOfExistingPileUp() const;
  // Whether the set of |Part|s in |parts_| is a strict subset of the set of
  // |Part|s in an existing |PileUp|.  Implies |SubsetOfExistingPileUp()|.
  bool StrictSubsetOfExistingPileUp() const;

  // Whether |Collect| has been called.
  bool collected_ = false;
  // if |SubsetOfExistingPileUp()|, |missing_| is the number of parts in that
  // |PileUp| that are not in this subset.
  int missing_;
  // The list of parts in this subset.
  std::list<not_null<ksp_plugin::Part*>> parts_;
  // The sum of the masses of the |parts_|.
  quantities::Mass total_mass_;
  // The sum of the |intrinsic_force|s on the |parts_|.
  geometry::Vector<quantities::Force, ksp_plugin::Barycentric>
      total_intrinsic_force_;
  // Whether the subset touches the ground.
  bool grounded_ = false;
};

}  // namespace base
}  // namespace principia
