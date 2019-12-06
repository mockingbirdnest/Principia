
#include "ksp_plugin/part_subsets.hpp"

#include <list>
#include <memory>

#include "base/not_null.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/pile_up.hpp"

namespace principia {

using base::make_not_null_shared;
using geometry::Instant;
using ksp_plugin::Barycentric;
using ksp_plugin::Part;
using ksp_plugin::PileUp;
using physics::DegreesOfFreedom;
using physics::Ephemeris;

namespace base {

Subset<Part>::Properties::Properties(not_null<ksp_plugin::Part*> const part)
    : total_mass_(part->inertia_tensor().mass()),
      total_intrinsic_force_(part->intrinsic_force()) {
  if (part->is_piled_up()) {
    missing_ = part->containing_pile_up()->parts().size() - 1;
  }
  parts_.emplace_back(part);
}

void Subset<Part>::Properties::MergeWith(Properties& other) {
  if (SubsetsOfSamePileUp(*this, other)) {
    // The subsets |*this| and |other| are disjoint.
    CHECK_EQ(missing_ - other.parts_.size(),
             other.missing_ - parts_.size());
    missing_ -= other.parts_.size();
    CHECK_GE(missing_, 0);
  } else {
    for (auto const part : parts_) {
      part->reset_containing_pile_up();
    }
    for (auto const part : other.parts_) {
      part->reset_containing_pile_up();
    }
  }
  parts_.splice(parts_.end(), other.parts_);
  total_mass_ += other.total_mass_;
  total_intrinsic_force_ += other.total_intrinsic_force_;
  grounded_ |= other.grounded_;
}

void Subset<Part>::Properties::Ground() {
  grounded_ = true;
}

bool Subset<Part>::Properties::grounded() const {
  return grounded_;
}

void Subset<Part>::Properties::Collect(
    PileUps& pile_ups,
    Instant const& t,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    Ephemeris<Barycentric>::FixedStepParameters const& fixed_step_parameters,
    not_null<Ephemeris<Barycentric>*> ephemeris) {
  if (collected_) {
    return;
  }
  collected_ = true;
  if (EqualsExistingPileUp()) {
    PileUp& pile_up = *parts_.front()->containing_pile_up();
    pile_up.set_mass(total_mass_);
    pile_up.set_intrinsic_force(total_intrinsic_force_);
  } else {
    if (StrictSubsetOfExistingPileUp()) {
      for (auto const part : parts_) {
        part->reset_containing_pile_up();
      }
    }
    CHECK(!parts_.empty());

    // First push a nullptr to be able to capture an iterator to the new
    // location in the list in the deletion callback.
    pile_ups.push_front(nullptr);
    auto deletion_callback = [it = pile_ups.begin(), &pile_ups]() {
      pile_ups.erase(it);
    };
    auto const pile_up =
        make_not_null_shared<PileUp>(std::move(parts_),
                                     t,
                                     adaptive_step_parameters,
                                     fixed_step_parameters,
                                     ephemeris,
                                     std::move(deletion_callback));
    *pile_ups.begin() = pile_up.get();

    // The pile-up is now co-owned by all its parts.
    for (not_null<Part*> const part : pile_up->parts()) {
      part->set_containing_pile_up(pile_up);
    }
  }
}

bool Subset<Part>::Properties::SubsetsOfSamePileUp(
    Properties const& left,
    Properties const& right) {
  return left.SubsetOfExistingPileUp() && right.SubsetOfExistingPileUp() &&
         left.parts_.front()->containing_pile_up() ==
             right.parts_.front()->containing_pile_up();
}

bool Subset<Part>::Properties::EqualsExistingPileUp() const {
  return SubsetOfExistingPileUp() && missing_ == 0;
}

bool Subset<Part>::Properties::SubsetOfExistingPileUp() const {
  return parts_.front()->is_piled_up();
}

bool Subset<Part>::Properties::StrictSubsetOfExistingPileUp()
    const {
  return SubsetOfExistingPileUp() && missing_ > 0;
}

}  // namespace base
}  // namespace principia
