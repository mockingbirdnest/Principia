#pragma once

#include <tuple>
#include <utility>

namespace principia {
namespace base {
namespace internal_for_all_of {

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

}  // namespace internal_for_all_of

// Iterates over all the tuples in parallel.  |F| must be a functor taking
// elements at corresponding positions in each of the tuples.
// Example:
//   std::tuple const t{1, 2.5, 3};
//   std::array const a{4, 5, 6};
//   for_all_of(t, a).loop([](auto const tuple_element, int const i) {
//     std::cout << tuple_element << " " << i << "\n";
//   });
template<typename... Tuple>
constexpr internal_for_all_of::Iteration<Tuple...> for_all_of(Tuple&&... tuple);

}  // namespace base
}  // namespace principia

#include "base/for_all_of_body.hpp"
