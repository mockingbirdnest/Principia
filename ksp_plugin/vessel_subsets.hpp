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
  class SubsetOfExistingPileUp {
   public:
    explicit SubsetOfExistingPileUp(ContainerIterator<PileUps> pile_up);

   private:
    ContainerIterator<PileUps> const pile_up_;
    int missing_;
    friend class Subset<ksp_plugin::Vessel>::Properties;
  };

  Properties(not_null<ksp_plugin::Vessel*> vessel);

  void MergeWith(Properties& other);

  void Collect(not_null<PileUps*> pile_ups);

 private:
  bool collected_ = false;
  std::experimental::optional<SubsetOfExistingPileUp>
      subset_of_existing_pile_up_ ;
  std::list<not_null<ksp_plugin::Vessel*>> vessels_;
};

template<>
not_null<Subset<ksp_plugin::Vessel>::Node*>
Subset<ksp_plugin::Vessel>::Node::Get(ksp_plugin::Vessel& element);

}  // namespace base
}  // namespace principia
