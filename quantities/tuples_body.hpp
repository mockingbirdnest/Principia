
#pragma once

#include "quantities/tuples.hpp"

#include <algorithm>
#include <tuple>
#include <utility>

#include "base/not_constructible.hpp"

namespace principia {
namespace quantities {
namespace internal_tuples {

using base::not_constructible;

template<template<typename> class Transform, typename Tuple, int... indices>
struct ApplyGenerator<Transform, Tuple, std::integer_sequence<int, indices...>>
    : not_constructible {
  using Type = std::tuple<Transform<std::tuple_element_t<indices, Tuple>>...>;
};

template<template<typename, typename> class Transform,
         typename LTuple, typename RTuple, std::size_t... indices>
struct Apply2Generator<Transform, LTuple, RTuple,
                       std::index_sequence<indices...>>
    : not_constructible {
  //TODO(phl): no conditional_t, ?:0 helper
  using Type = std::tuple<
      Transform<std::conditional_t<
                    (indices < std::tuple_size_v<LTuple>),
                    std::tuple_element_t<
                        (indices < std::tuple_size_v<LTuple> ? indices : 0),
                        LTuple>,
                    void>,
                std::conditional_t<
                    (indices < std::tuple_size_v<RTuple>),
                    std::tuple_element_t<
                        (indices < std::tuple_size_v<RTuple> ? indices : 0),
                        RTuple>,
                    void>>...>;
};

template<typename Value, typename Argument, int n, int... orders>
struct DerivativesGenerator<Value, Argument, n,
                            std::integer_sequence<int, orders...>>
    : not_constructible {
  using Type = std::tuple<Derivative<Value, Argument, orders>...>;
};

}  // namespace internal_tuples
}  // namespace quantities
}  // namespace principia
