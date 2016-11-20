#pragma once

#include <list>

#include "base/disjoint_sets.hpp"

#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel.hpp"

// This ksp_plugin file is in namespace |base| to specialize templates declared
// therein.

namespace principia {
namespace base {

template<>
class Subset<ksp_plugin::Vessel>::Properties {
  using PileUps = std::list<ksp_plugin::PileUp>;
 public:

  Properties(not_null<ksp_plugin::Vessel*> vessel);

  // If |*this| and |other| have |subset_of_existing_pile_up_| corresponding
  // to different |PileUp|s (or one is a subset and not the other), the relevant
  // |PileUp|s are erased. Otherwise, |this->subset_of_existing_pile_up_| keeps
  // track of the number of missing vessels.
  // Maintains |vessels_| by joining the lists.
  void MergeWith(Properties& other);

  // If |collected_|, performs no action.
  // Otherwise, sets |collected_|, and:
  // - if |subset_of_existing_pile_up_| and
  //   |subset_of_existing_pile_up_->missing == 0|, performs no action.
  // - if |subset_of_existing_pile_up_| but
  //   |subset_of_existing_pile_up_->missing > 0|, erases
  //   |subset_of_existing_pile_up_| inserts a new |PileUp| into |pile_ups| with
  //   the vessels in |vessels_|.
  // - if |!subset_of_existing_pile_up_|, inserts a new |PileUp| into |pile_ups|
  //   with the vessels in |vessels_|.
  void Collect(not_null<PileUps*> pile_ups);

 private:
  class SubsetOfExistingPileUp {
   public:
    explicit SubsetOfExistingPileUp(ContainerIterator<PileUps> pile_up);

    ContainerIterator<PileUps> const pile_up_;
    int missing_;
  };

  bool collected_ = false;
  // Has a value if and only if the set of elements of |vessels_| is a subset
  // of some pre-existing |PileUp|. |subset_of_existing_pile_up_->pile_up_|
  // refers to that |PileUp|, and |subset_of_existing_pile_up_->missing_| is the
  // difference of the cardinalities.
  std::experimental::optional<SubsetOfExistingPileUp>
      subset_of_existing_pile_up_;
  // The list of vessels in this subset.
  std::list<not_null<ksp_plugin::Vessel*>> vessels_;
};

template<>
not_null<Subset<ksp_plugin::Vessel>::Node*>
Subset<ksp_plugin::Vessel>::Node::Get(ksp_plugin::Vessel& element);

}  // namespace base
}  // namespace principia
