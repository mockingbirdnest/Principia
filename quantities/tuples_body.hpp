
#pragma once

#include "quantities/tuples.hpp"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <utility>

#include "base/not_constructible.hpp"

namespace principia {
namespace quantities {
namespace internal_tuples {

using base::not_constructible;

// A helper for getting the type of an element even for an index that is not in
// range.  In that case, the member Type is defined to be void.  Doing this with
// conditional_t is awkward.
template<typename Tuple, std::size_t index,
         bool in_range = (index < std::tuple_size_v<Tuple>)>
struct TupleElementOrVoid;

template<typename Tuple, std::size_t index>
struct TupleElementOrVoid<Tuple, index, false> {
  using Type = void;
};

template<typename Tuple, std::size_t index>
struct TupleElementOrVoid<Tuple, index, true> {
  using Type = std::tuple_element_t<index, Tuple>;
};

template<template<typename...> typename Transform,
         typename... Tuples>
struct ApplyGenerator {
  template<std::size_t index>
  using ResultElement =
      Transform<typename TupleElementOrVoid<Tuples, index>::Type...>;

  template<typename index_sequence>
  struct Result;

  template<std::size_t... indices>
  struct Result<std::index_sequence<indices...>> {
    using Type = std::tuple<ResultElement<indices>...>;
  };

  static constexpr std::size_t result_size = std::max({std::tuple_size_v<Tuples>...});

  using Type = typename Result<std::make_index_sequence<result_size>>::Type;
};

template<typename Value, typename Argument, int n, std::size_t... orders>
struct DerivativesGenerator<Value, Argument, n,
                            std::index_sequence<orders...>>
    : not_constructible {
  using Type = std::tuple<Derivative<Value, Argument, orders>...>;
};

}  // namespace internal_tuples
}  // namespace quantities
}  // namespace principia
