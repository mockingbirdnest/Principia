
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

template<template<typename> class Transform,
         typename Tuple, std::size_t... indices>
struct ApplyGenerator<Transform, Tuple, std::index_sequence<indices...>>
    : not_constructible {
  using Type = std::tuple<Transform<std::tuple_element_t<indices, Tuple>>...>;
};

template<template<typename, typename> class Transform,
         typename LTuple, typename RTuple, std::size_t... indices>
struct Apply2Generator<Transform, LTuple, RTuple,
                       std::index_sequence<indices...>>
    : not_constructible {
  using Type = std::tuple<
      Transform<typename TupleElementOrVoid<LTuple, indices>::Type,
                typename TupleElementOrVoid<RTuple, indices>::Type>...>;
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
