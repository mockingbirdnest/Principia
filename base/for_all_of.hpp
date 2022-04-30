#pragma once
#include <tuple>
#include <utility>

namespace principia {
namespace base {

template<std::size_t i = 0, typename F, typename Tuple>
constexpr void for_all_of(Tuple& tuple, F f) {
  if constexpr (i < std::tuple_size_v<Tuple>) {
    f(std::get<i>(tuple));
    for_each<i + 1, F, Tuple>(tuple, f);
  }
}

}  // namespace base
}  // namespace principia
