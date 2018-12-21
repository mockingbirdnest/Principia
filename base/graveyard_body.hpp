
#pragma once

#include "base/graveyard.hpp"

namespace principia {
namespace base {

template<typename T>
Graveyard<T>::Graveyard(std::int64_t const number_of_threads)
    : gravedigger_(number_of_threads) {}

template<typename T>
void Graveyard<T>::Bury(std::unique_ptr<T> t) {
  gravedigger_.Add([coffin = t.release()]() {
    delete coffin;
  });
}

}  // namespace base
}  // namespace principia
