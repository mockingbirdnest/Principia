
#include "ksp_plugin/pile_up.hpp"

#include <list>

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

PileUp::PileUp(std::list<not_null<Vessel*>>&& vessels)
    : vessels_(std::move(vessels)) {
  LOG(ERROR) << "Constructing PileUp " << this
             << " with the following vessels:";
  for (auto const vessel : vessels_) {
    LOG(ERROR) << vessel;
  }
}

PileUp::~PileUp() {
  LOG(ERROR) << "Destroying PileUp " << this;
}

std::list<not_null<Vessel*>> const& PileUp::vessels() const {
  return vessels_;
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
