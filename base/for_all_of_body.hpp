#pragma once

#include "base/for_all_of.hpp"

#include <tuple>
#include <utility>

#include "base/macros.hpp"  // 🧙 For FORCE_INLINE.

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
    std::apply(
        [&f](Tuple&&... tuple) {
          using namespace std;
          f(get<i>(tuple)...);
        },
        all_the_tuples_);
    loop<i + 1, F>(f);
  }
}

template<typename... Tuple>
template<std::size_t i, typename F>
constexpr void Iteration<Tuple...>::loop_indexed(F const& f) {
  if constexpr (i < size) {
    std::apply(
        [&f](Tuple&&... tuple) {
          using namespace std;
          f.template operator()<i>(get<i>(tuple)...);
        },
        all_the_tuples_);
    loop_indexed<i + 1, F>(f);
  }
}

template<typename... Tuple>
constexpr Iteration<Tuple...> for_all_of(Tuple&&... tuple) {
  return Iteration<Tuple...>(std::forward<Tuple>(tuple)...);
}

template<std::int64_t begin, std::int64_t end>
template<std::int64_t i, typename F>
FORCE_INLINE constexpr void for_integer_range<begin, end>::loop(F const& f) {
  if constexpr (i != end) {
    f.template operator()<i>();
    loop<i + 1, F>(f);
  }
}

}  // namespace internal
}  // namespace _for_all_of
}  // namespace base
}  // namespace principia
