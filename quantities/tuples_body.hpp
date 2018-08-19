
#pragma once

#include "quantities/tuples.hpp"

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
  using Type = std::tuple<
      typename Transform<std::tuple_element_t<indices, Tuple>>...>;
};

template<typename Value, typename Argument, int n, int... orders>
struct DerivativesGenerator<Value,
                            Argument,
                            n,
                            std::integer_sequence<int, orders...>>
    : not_constructible {
  using Type = std::tuple<Derivative<Value, Argument, orders>...>;
};

}  // namespace internal_tuples
}  // namespace quantities
}  // namespace principia
