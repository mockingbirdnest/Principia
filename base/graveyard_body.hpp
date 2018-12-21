
#pragma once

#include "base/graveyard.hpp"

namespace principia {
namespace base {

inline Graveyard::Graveyard(std::int64_t const number_of_threads)
    : gravedigger_(number_of_threads) {}

template<typename T>
void Graveyard::Bury(std::unique_ptr<T> t) {
  // TODO(egg): Investigate the possibility of a mutable lambda with
  // std::packaged_task in the ThreadPool instead of std::function.
  gravedigger_.Add([coffin = t.release()]() {
    delete coffin;
  });
}

}  // namespace base
}  // namespace principia
