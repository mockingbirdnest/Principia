#include "ksp_plugin/vessel_subsets.hpp"

#include <list>

#include "ksp_plugin/pile_up.hpp"

namespace principia {

using ksp_plugin::Vessel;
using ksp_plugin::PileUp;

namespace base {

Subset<Vessel>::Properties::SubsetOfExistingPileUp::SubsetOfExistingPileUp(
    ContainerIterator<PileUps> pile_up)
    : pile_up_(pile_up) {
  missing_ = pile_up_.iterator->vessels().size() - 1;
}

Subset<Vessel>::Properties::Properties(not_null<ksp_plugin::Vessel*> vessel) {
  if (vessel->piled_up()) {
    subset_of_existing_pile_up_.emplace(
        SubsetOfExistingPileUp(*vessel->containing_pile_up()));
  }
  vessels_.emplace_back(vessel);
}

void Subset<ksp_plugin::Vessel>::Properties::Collect(
    not_null<PileUps*> pile_ups) {
  if (!collected_ && !(EqualsExistingPileUp())) {
    collected_ = true;
    pile_ups->emplace_front(std::move(vessels_));
    auto const it = pile_ups->begin();
    for (not_null<Vessel*> const vessel : it->vessels()) {
      vessel->set_containing_pile_up(ContainerIterator<PileUps>(pile_ups, it));
    }
  }
}

bool Subset<ksp_plugin::Vessel>::Properties::SubsetsOfSamePileUp(
    Properties const& left,
    Properties const& right) {
  return left.subset_of_existing_pile_up_ &&
         right.subset_of_existing_pile_up_ &&
         left.subset_of_existing_pile_up_->pile_up_.iterator ==
             right.subset_of_existing_pile_up_->pile_up_.iterator;
}

bool Subset<ksp_plugin::Vessel>::Properties::EqualsExistingPileUp() const {
  return subset_of_existing_pile_up_ &&
         subset_of_existing_pile_up_->missing_ == 0;
}

bool Subset<ksp_plugin::Vessel>::Properties::StrictSubsetOfExistingPileUp()
    const {
  return subset_of_existing_pile_up_ &&
         subset_of_existing_pile_up_->missing_ > 0;
}

void Subset<Vessel>::Properties::MergeWith(Properties& other) {
  if (SubsetsOfSamePileUp(*this, other)) {
    // The subsets |*this| and |other| are disjoint.
    CHECK_EQ(subset_of_existing_pile_up_->missing_ - other.vessels_.size(),
             other.subset_of_existing_pile_up_->missing_ - vessels_.size());
    subset_of_existing_pile_up_->missing_ -= other.vessels_.size();
    CHECK_GE(subset_of_existing_pile_up_->missing_, 0);
  } else {
    if (subset_of_existing_pile_up_) {
      subset_of_existing_pile_up_->pile_up_.Erase();
    }
    if (other.subset_of_existing_pile_up_) {
      other.subset_of_existing_pile_up_->pile_up_.Erase();
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
