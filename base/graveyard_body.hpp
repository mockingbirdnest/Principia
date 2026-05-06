#pragma once

#include "base/graveyard.hpp"

#include <utility>

namespace principia {
namespace base {
namespace _graveyard {
namespace internal {

inline Graveyard::Graveyard(std::int64_t const number_of_threads)
    : gravedigger_(number_of_threads) {}

template<std::movable T>
void Graveyard::Bury(T t) {
  // TODO(egg): Investigate the possibility of a std::move_only_function in the
  // ThreadPool instead of std::function.
  gravedigger_.Add([coffin = new T(std::move(t))]() { delete coffin; });
}

}  // namespace internal
}  // namespace _graveyard
}  // namespace base
}  // namespace principia
