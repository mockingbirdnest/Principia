#pragma once

#include "base/disjoint_sets.hpp"

#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel.hpp"

// This ksp_plugin file is in namespace |base| to specialize templates declared
// therein.

namespace principia {
namespace base {

template<>
class Subset<ksp_plugin::Vessel>::Properties {
 public:
  class SubsetOfExistingPileUp {
   public:
    SubsetOfExistingPileUp(not_null<std::list<ksp_plugin::PileUp>*> pile_ups,
                           std::list<ksp_plugin::PileUp>::iterator pile_up);

   private:
    not_null<std::list<ksp_plugin::PileUp>*> const pile_ups_;
    std::list<ksp_plugin::PileUp>::iterator const pile_up_;
    int missing_;
    friend class Subset<ksp_plugin::Vessel>::Properties;
  };

  Properties(not_null<ksp_plugin::Vessel*> vessel,
             std::experimental::optional<SubsetOfExistingPileUp>
                 subset_of_existing_pile_up);

  void MergeWith(Properties& other);

 private:
  std::experimental::optional<SubsetOfExistingPileUp>
      subset_of_existing_pile_up_;
  std::list<not_null<ksp_plugin::Vessel*>> vessels_;
};

template<>
not_null<Subset<ksp_plugin::Vessel>::Node*>
Subset<ksp_plugin::Vessel>::Node::Get(ksp_plugin::Vessel& element);

}  // namespace base
}  // namespace principia
