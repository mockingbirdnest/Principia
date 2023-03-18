#pragma once

#include "base/for_all_of.hpp"

#include <tuple>
#include <utility>

namespace principia {
namespace base {
namespace _for_all_of {
namespace internal {

template<typename... Tuple>
constexpr Iteration<Tuple...>::Iteration(Tuple&&... tuple)
    : all_the_tuples_(std::forward_as_tuple(tuple...)) {}

template<typename... Tuple>
template<std::size_t i, typename F>
constexpr void Iteration<Tuple...>::loop(F const& f) {
  if constexpr (i < size) {
    std::apply([&f](Tuple&&... tuple) { f(std::get<i>(tuple)...); },
               all_the_tuples_);
    loop<i + 1, F>(f);
  }
}

template<typename... Tuple>
constexpr Iteration<Tuple...> for_all_of(Tuple&&... tuple) {
  return Iteration<Tuple...>(std::forward<Tuple>(tuple)...);
}

}  // namespace internal
}  // namespace _for_all_of
}  // namespace base
}  // namespace principia
