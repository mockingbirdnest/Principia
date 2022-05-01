#pragma once

#include <tuple>
#include <utility>

namespace principia {
namespace base {

namespace internal_for_all_of {
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

}  // namespace internal_for_all_of

template<typename... Tuple>
constexpr internal_for_all_of::Iteration<Tuple...> for_all_of(
    Tuple&&... tuple) {
  return internal_for_all_of::Iteration<Tuple...>(
      std::forward<Tuple>(tuple)...);
}

}  // namespace base
}  // namespace principia
