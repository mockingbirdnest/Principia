#pragma once

#include <tuple>
#include <utility>

namespace principia {
namespace base {
namespace _for_all_of {
namespace internal {

template<typename... Tuple>
class Iteration {
  static constexpr size_t size =
      (std::tuple_size_v<std::remove_cvref_t<Tuple>>, ...);
  static_assert(((std::tuple_size_v<std::remove_cvref_t<Tuple>> == size) &&
                 ...),
                "Parallel iteration must apply to containers of equal size");

 public:
  constexpr Iteration(Tuple&&... tuple);

  template<std::size_t i = 0, typename F>
  constexpr void loop(F const& f);

 private:
  std::tuple<Tuple&&...> all_the_tuples_;
};

// Iterates over all the tuples in parallel.  |F| must be a functor taking
// elements at corresponding positions in each of the tuples.
// Example:
//   std::tuple const t{"a", 2.5, 3};
//   std::array const a{4, 5, 6};
//   for_all_of(t, a).loop([](auto const tuple_element, int const i) {
//     std::cout << tuple_element << " " << i << "\n";
//   });
template<typename... Tuple>
constexpr Iteration<Tuple...> for_all_of(Tuple&&... tuple);

}  // namespace internal

using internal::for_all_of;

}  // namespace _for_all_of
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_for_all_of;
}  // namespace principia::base

#include "base/for_all_of_body.hpp"
