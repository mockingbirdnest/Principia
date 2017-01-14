
#include "ksp_plugin/vessel_subsets.hpp"

#include <list>

#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel.hpp"

namespace principia {

using ksp_plugin::Vessel;
using ksp_plugin::PileUp;

namespace base {

Subset<Vessel>::Properties::Properties(not_null<ksp_plugin::Vessel*> vessel) {
  if (vessel->is_piled_up()) {
    missing_ = vessel->containing_pile_up()->iterator()->vessels().size() - 1;
  }
  vessels_.emplace_back(vessel);
  total_mass_ = vessel->mass();
  total_intrinsic_force_ = vessel->intrinsic_force();
}

void Subset<ksp_plugin::Vessel>::Properties::Collect(
    not_null<PileUps*> const pile_ups) {
  if (collected_) {
    return;
  }
  collected_ = true;
  if (!EqualsExistingPileUp()) {
    if (StrictSubsetOfExistingPileUp()) {
      vessels_.front()->clear_pile_up();
    }
    pile_ups->emplace_front(std::move(vessels_));
    auto const it = pile_ups->begin();
    for (not_null<Vessel*> const vessel : it->vessels()) {
      vessel->set_containing_pile_up(IteratorOn<PileUps>(pile_ups, it));
    }
    pile_up = it;
  } else {
    PileUp& pile_up = *vessels_.front()->containing_pile_up()->iterator();
    pile_up.set_mass(total_mass_);
    pile_up.set_intrinsic_force(total_intrinsic_force_);
  }
}

bool Subset<ksp_plugin::Vessel>::Properties::SubsetsOfSamePileUp(
    Properties const& left,
    Properties const& right) {
  return left.SubsetOfExistingPileUp() && right.SubsetOfExistingPileUp() &&
         left.vessels_.front()->containing_pile_up()->iterator() ==
             right.vessels_.front()->containing_pile_up()->iterator();
}

bool Subset<ksp_plugin::Vessel>::Properties::EqualsExistingPileUp() const {
  return SubsetOfExistingPileUp() && missing_ == 0;
}

bool Subset<ksp_plugin::Vessel>::Properties::SubsetOfExistingPileUp() const {
  return vessels_.front()->is_piled_up();
}

bool Subset<ksp_plugin::Vessel>::Properties::StrictSubsetOfExistingPileUp()
    const {
  return SubsetOfExistingPileUp() && missing_ > 0;
}

void Subset<Vessel>::Properties::MergeWith(Properties& other) {
  if (SubsetsOfSamePileUp(*this, other)) {
    // The subsets |*this| and |other| are disjoint.
    CHECK_EQ(missing_ - other.vessels_.size(),
             other.missing_ - vessels_.size());
    missing_ -= other.vessels_.size();
    CHECK_GE(missing_, 0);
  } else {
    vessels_.front()->clear_pile_up();
    other.vessels_.front()->clear_pile_up();
  }
  vessels_.splice(vessels_.end(), other.vessels_);
  total_mass_ += other.total_mass_;
  total_intrinsic_force_ += other.total_intrinsic_force_;
}

}  // namespace base
}  // namespace principia
