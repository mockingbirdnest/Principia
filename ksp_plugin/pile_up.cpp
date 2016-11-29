
#include "ksp_plugin/pile_up.hpp"

#include <list>

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

PileUp::PileUp(std::list<not_null<Vessel*>> vessels)
    : vessels_(std::move(vessels)) {}

std::list<not_null<Vessel*>> const& PileUp::vessels() const {
  return vessels_;
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
