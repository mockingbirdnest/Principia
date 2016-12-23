#pragma once

#include <list>

#include "base/disjoint_sets.hpp"

#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel.hpp"

// This ksp_plugin file is in namespace |base| to specialize templates declared
// therein.

namespace principia {
namespace base {

// Within an union-find on |Vessel|s, we maintain lists of the elements in the
// disjoint sets.  Moreover, we keep track of the inclusion relations of those
// sets to the sets of |Vessel|s in existing |PileUp|s, destroying existing
// |PileUp|s as we learn that they will not appear in the new arrangement.
// The |Collect| operation finalizes this, destroying existing |PileUp| which
// are strict supersets of the new sets, and creating the new |PileUp|s.
template<>
class Subset<ksp_plugin::Vessel>::Properties final {
  using PileUps = std::list<ksp_plugin::PileUp>;

 public:
  Properties(not_null<ksp_plugin::Vessel*> vessel);

  // If |*this| and |other| are subsets of different |PileUp|s, or one is a
  // subset and not the other, the relevant |PileUp|s are erased.
  // Otherwise, |this->subset_of_existing_pile_up_| keeps track of the number of
  // missing vessels.
  // Maintains |vessels_| by joining the lists.
  void MergeWith(Properties& other);

  // If |collected_|, performs no action.
  // Otherwise, sets |collected_|, and:
  // - if |EqualsExistingPileUp()|, performs no action.
  // - if |StrictSubsetOfExistingPileUp()|, erases the existing |PileUp| inserts
  //   a new |PileUp| into |pile_ups| with the vessels in |vessels_|.
  // - if |!subset_of_existing_pile_up_|, inserts a new |PileUp| into |pile_ups|
  //   with the vessels in |vessels_|.
  void Collect(not_null<PileUps*> pile_ups);

 private:
  // Whether |left| and |right| are both subsets of the same existing |PileUp|.
  static bool SubsetsOfSamePileUp(Properties const& left,
                                  Properties const& right);
  // Whether the set of |Vessel|s in |vessels_| is equal to the set of |Vessel|s
  // in an existing |PileUp|.  In that case
  // |subset_of_existing_pile_up_->pile_up_| is that |PileUp|.
  bool EqualsExistingPileUp() const;
  // Whether the set of |Vessel|s in |vessels_| is a strict subset of the set of
  // |Vessel|s in an existing |PileUp|.  In that case
  // |subset_of_existing_pile_up_->pile_up_| is that |PileUp|.
  bool StrictSubsetOfExistingPileUp() const;

  // Information about a subset of the set of |Vessel|s in a |PileUp|.  Keeps
  // track of the difference of cardinalities.
  class SubsetOfExistingPileUp final {
   public:
    explicit SubsetOfExistingPileUp(IteratorOn<PileUps> pile_up);

    // The existing |PileUp|.
    IteratorOn<PileUps> const pile_up_;
    // The difference between the cardinalities; never negative.
    int missing_;
  };

  // Whether |Collect| has been called.
  bool collected_ = false;
  // Has a value if and only if the set of elements of |vessels_| is a subset
  // of some pre-existing |PileUp|.
  std::experimental::optional<SubsetOfExistingPileUp>
      subset_of_existing_pile_up_;
  // The list of vessels in this subset.
  std::list<not_null<ksp_plugin::Vessel*>> vessels_;
};

// TODO(egg): figure out why the compiler complains if I put the implementation
// in the cpp.
template<> not_null<Subset<ksp_plugin::Vessel>::Node*>
Subset<ksp_plugin::Vessel>::Node::Get(ksp_plugin::Vessel& element) {
  return element.subset_node_.get();
}

}  // namespace base
}  // namespace principia
