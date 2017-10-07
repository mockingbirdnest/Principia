
#include "ksp_plugin/part_subsets.hpp"

#include <list>

#include "ksp_plugin/part.hpp"
#include "ksp_plugin/pile_up.hpp"

namespace principia {

using geometry::Instant;
using ksp_plugin::Barycentric;
using ksp_plugin::Part;
using ksp_plugin::PileUp;
using physics::DegreesOfFreedom;
using physics::Ephemeris;

namespace base {

Subset<Part>::Properties::Properties(not_null<ksp_plugin::Part*> const part,
                                     bool const grounded)
    : total_mass_(part->mass()),
      total_intrinsic_force_(part->intrinsic_force()),
      grounded_(grounded) {
  if (part->is_piled_up()) {
    missing_ = part->containing_pile_up()->iterator()->parts().size() - 1;
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
    parts_.front()->clear_pile_up();
    other.parts_.front()->clear_pile_up();
  }
  parts_.splice(parts_.end(), other.parts_);
  total_mass_ += other.total_mass_;
  total_intrinsic_force_ += other.total_intrinsic_force_;
  grounded_ |= other.grounded_;
}

void Subset<ksp_plugin::Part>::Properties::Ground() {
  grounded_ = true;
}

void Subset<ksp_plugin::Part>::Properties::Collect(
    not_null<PileUps*> const pile_ups,
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
    PileUp& pile_up = *parts_.front()->containing_pile_up()->iterator();
    pile_up.set_mass(total_mass_);
    pile_up.set_intrinsic_force(total_intrinsic_force_);
  } else {
    if (StrictSubsetOfExistingPileUp()) {
      parts_.front()->clear_pile_up();
    }
    pile_ups->emplace_front(std::move(parts_),
                            t,
                            adaptive_step_parameters,
                            fixed_step_parameters,
                            ephemeris);
    auto const it = pile_ups->begin();
    for (not_null<Part*> const part : it->parts()) {
      part->set_containing_pile_up(IteratorOn<PileUps>(pile_ups, it));
    }
  }
}

bool Subset<ksp_plugin::Part>::Properties::SubsetsOfSamePileUp(
    Properties const& left,
    Properties const& right) {
  return left.SubsetOfExistingPileUp() && right.SubsetOfExistingPileUp() &&
         left.parts_.front()->containing_pile_up()->iterator() ==
             right.parts_.front()->containing_pile_up()->iterator();
}

bool Subset<ksp_plugin::Part>::Properties::EqualsExistingPileUp() const {
  return SubsetOfExistingPileUp() && missing_ == 0;
}

bool Subset<ksp_plugin::Part>::Properties::SubsetOfExistingPileUp() const {
  return parts_.front()->is_piled_up();
}

bool Subset<ksp_plugin::Part>::Properties::StrictSubsetOfExistingPileUp()
    const {
  return SubsetOfExistingPileUp() && missing_ > 0;
}

}  // namespace base
}  // namespace principia
