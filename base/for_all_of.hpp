#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

namespace principia {
namespace base {
namespace _for_all_of {
namespace internal {

template<typename... Tuple>
class Iteration {
  static constexpr std::size_t size =
      (std::tuple_size_v<std::remove_cvref_t<Tuple>>, ...);
  static_assert(((std::tuple_size_v<std::remove_cvref_t<Tuple>> == size) &&
                 ...),
                "Parallel iteration must apply to containers of equal size");

 public:
  constexpr Iteration(Tuple&&... tuple);  // NOLINT(runtime/explicit)

  template<std::size_t i = 0, typename F>
  constexpr void loop(F const& f);

  template<std::size_t i = 0, typename F>
  constexpr void loop_indexed(F const& f);

 private:
  std::tuple<Tuple&&...> all_the_tuples_;
};

// Iterates over all the tuples in parallel.  `F` must be a functor taking
// elements at corresponding positions in each of the tuples.
// Example:
//   std::tuple const t{"a", 2.5, 3};
//   std::array const a{4, 5, 6};
//   for_all_of(t, a).loop([](auto const tuple_element,
//                            int const array_element) {
//     std::cout << tuple_element << " " << array_element << "\n";
//   });
// For `loop_indexed`, `F::operator()` must also take an index as a
// template parameter:
//   for_all_of(t, a).loop_indexed([]<int i>(auto const tuple_element,
//                                           int const array_element) {
//     std::cout << i << " " << tuple_element << " " << array_element << "\n";
//   });
template<typename... Tuple>
constexpr Iteration<Tuple...> for_all_of(Tuple&&... tuple);

// Iterates over the integers in [begin, end[.  `F::operator()` must be parameterless and
// take an index as a template parameter:
//   std::tuple t{std::string("a"), 2.5, 3};
//   for_integer_range<0, 3>().loop([&]<int i>() {
//     if constexpr (i == 0) {
//       get<i>(t) += std::to_string(i);
//     } else {
//       get<i>(t) += i;
//     }
//   });
template<int begin, int end>
class for_integer_range {
 public:
  constexpr for_integer_range() = default;

  template<int i = begin, typename F>
  constexpr void loop(F const& f);
};

}  // namespace internal

using internal::for_all_of;
using internal::for_integer_range;

}  // namespace _for_all_of
}  // namespace base
}  // namespace principia

#include "base/for_all_of_body.hpp"
