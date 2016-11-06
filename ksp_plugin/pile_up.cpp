#include "ksp_plugin/pile_up.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_pile_up {

std::list<not_null<Vessel*>> const& PileUp::vessels() const {
  return vessels_;
}

}  // namespace internal_pile_up
}  // namespace ksp_plugin
}  // namespace principia
