
#pragma once

#include "quantities/generators.hpp"

#include "base/not_constructible.hpp"
#include "quantities/dimensions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {
namespace internal_generators {

using base::not_constructible;
using internal_dimensions::Dimensions;
using internal_dimensions::DimensionsExponentiationGenerator;
using internal_dimensions::DimensionsNthRootGenerator;
using internal_dimensions::DimensionsProductGenerator;
using internal_dimensions::DimensionsQuotientGenerator;
using internal_dimensions::NoDimensions;

template<typename Q>
struct Collapse : not_constructible {
  using Type = Q;
};

template<>
struct Collapse<Quantity<NoDimensions>> : not_constructible {
  using Type = double;
};

template<typename Q, int n>
struct ExponentiationGenerator : not_constructible {
  using Type =
      typename Collapse<
          Quantity<typename DimensionsExponentiationGenerator<
                                typename Q::Dimensions, n>::Type>>::Type;
};

template<int n>
struct ExponentiationGenerator<double, n> : not_constructible {
  using Type = double;
};

template<typename Q, int n, typename>
struct NthRootGenerator : not_constructible {
  using Type =
      typename Collapse<
          Quantity<typename DimensionsNthRootGenerator<
                                typename Q::Dimensions, n>::Type>>::Type;
};

// NOTE(phl): We use |is_arithmetic| here, not |double|, to make it possible to
// write something like |Sqrt(2)|.  We could use |is_arithmetic| in more places
// but it would make the template magic even harder to follow, so let's not do
// that until we have a good reason.
template<typename Q, int n>
struct NthRootGenerator<Q, n, std::enable_if_t<std::is_arithmetic<Q>::value>>
    : not_constructible {
  using Type = double;
};

template<typename Left, typename Right>
struct ProductGenerator : not_constructible {
  using Type =
      typename Collapse<
          Quantity<typename DimensionsProductGenerator<
                                typename Left::Dimensions,
                                typename Right::Dimensions>::Type>>::Type;
};

template<typename Left>
struct ProductGenerator<Left, double> : not_constructible {
  using Type = Left;
};

template<typename Right>
struct ProductGenerator<double, Right> : not_constructible {
  using Type = Right;
};

template<>
struct ProductGenerator<double, double> : not_constructible {
  using Type = double;
};

template<typename Left, typename Right>
struct QuotientGenerator : not_constructible {
  using Type =
      typename Collapse<
          Quantity<typename DimensionsQuotientGenerator<
                                typename Left::Dimensions,
                                typename Right::Dimensions>::Type>>::Type;
};

template<typename Left>
struct QuotientGenerator<Left, double> : not_constructible {
  using Type = Left;
};

template<typename Right>
struct QuotientGenerator<double, Right> : not_constructible {
  using Type =
      typename Collapse<
          Quantity<typename DimensionsQuotientGenerator<
                                NoDimensions,
                                typename Right::Dimensions>::Type>>::Type;
};

template<>
struct QuotientGenerator<double, double> : not_constructible {
  using Type = double;
};

template<template<typename> class Transform, typename Qs, int... indices>
struct ApplyGenerator<Transform, Qs, std::integer_sequence<int, indices...>>
    : not_constructible {
  using Type = std::tuple<
      typename Transform<std::tuple_element_t<indices, Qs>>...>;
};

}  // namespace internal_generators
}  // namespace quantities
}  // namespace principia
