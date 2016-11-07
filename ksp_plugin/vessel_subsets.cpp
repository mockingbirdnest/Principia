#include "ksp_plugin/vessel_subsets.hpp"

#include <list>

#include "ksp_plugin/pile_up.hpp"

namespace principia {

using ksp_plugin::Vessel;
using ksp_plugin::PileUp;

namespace base {

Subset<Vessel>::Properties::SubsetOfExistingPileUp::SubsetOfExistingPileUp(
    not_null<std::list<ksp_plugin::PileUp>*> pile_ups,
    std::list<ksp_plugin::PileUp>::iterator pile_up)
    : pile_ups_(pile_ups),
      pile_up_(pile_up) {
  missing_ = pile_up_->vessels().size() - 1;
}

Subset<Vessel>::Properties::Properties(
    not_null<Vessel*> vessel,
    std::experimental::optional<SubsetOfExistingPileUp>
        subset_of_existing_pile_up)
    : subset_of_existing_pile_up_(subset_of_existing_pile_up) {
  vessels_.emplace_back(vessel);
}

void Subset<Vessel>::Properties::MergeWith(Properties& other) {
  if (subset_of_existing_pile_up_ && other.subset_of_existing_pile_up_ &&
      subset_of_existing_pile_up_->pile_up_ ==
          other.subset_of_existing_pile_up_->pile_up_) {
    CHECK_EQ(subset_of_existing_pile_up_->missing_ - other.vessels_.size(),
             other.subset_of_existing_pile_up_->missing_ - vessels_.size());
    subset_of_existing_pile_up_->missing_ -= other.vessels_.size();
    CHECK_GE(subset_of_existing_pile_up_->missing_, 0);
  } else {
    if (subset_of_existing_pile_up_) {
      subset_of_existing_pile_up_->pile_ups_->erase(
          subset_of_existing_pile_up_->pile_up_);
    }
    if (other.subset_of_existing_pile_up_) {
      other.subset_of_existing_pile_up_->pile_ups_->erase(
          other.subset_of_existing_pile_up_->pile_up_);
    }
    subset_of_existing_pile_up_ = std::experimental::nullopt;
  }
  vessels_.splice(vessels_.end(), other.vessels_);
}

template<>
not_null<Subset<ksp_plugin::Vessel>::Node*>
Subset<ksp_plugin::Vessel>::Node::Get(ksp_plugin::Vessel& element) {
  return element.subset_node_.get();
}

}  // namespace base
}  // namespace principia
