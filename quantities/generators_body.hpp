
#pragma once

#include "quantities/generators.hpp"

#include <tuple>

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

// The template template parameter |Quantity| on specializations lifts a
// circular dependency.

template<typename Q>
struct Collapse : not_constructible {
  using Type = Q;
};

template<template<typename> typename Quantity>
struct Collapse<Quantity<NoDimensions>> : not_constructible {
  using Type = double;
};

template<template<typename> typename Quantity, typename D, int n>
struct ExponentiationGenerator<Quantity<D>, n> : not_constructible {
  using Type = typename Collapse<
      Quantity<typename DimensionsExponentiationGenerator<D, n>::Type>>::Type;
};

template<int n>
struct ExponentiationGenerator<double, n> : not_constructible {
  using Type = double;
};

template<template<typename> typename Quantity, typename D, int n>
struct NthRootGenerator<Quantity<D>, n, void> : not_constructible {
  using Type = typename Collapse<
      Quantity<typename DimensionsNthRootGenerator<D, n>::Type>>::Type;
};

// NOTE(phl): We use |is_arithmetic| here, not |double|, to make it possible to
// write something like |Sqrt(2)|.  We could use |is_arithmetic| in more places
// but it would make the template magic even harder to follow, so let's not do
// that until we have a good reason.
template<typename Q, int n>
struct NthRootGenerator<Q, n, std::enable_if_t<std::is_arithmetic_v<Q>>>
    : not_constructible {
  using Type = double;
};

template<template<typename> typename Quantity, typename Left, typename Right>
struct ProductGenerator<Quantity<Left>, Quantity<Right>> : not_constructible {
  using Type = typename Collapse<Quantity<
      typename DimensionsProductGenerator<Left, Right>::Type>>::Type;
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

template<template<typename> typename Quantity, typename Left, typename Right>
struct QuotientGenerator<Quantity<Left>, Quantity<Right>> : not_constructible {
  using Type = typename Collapse<Quantity<
      typename DimensionsQuotientGenerator<Left, Right>::Type>>::Type;
};

template<typename Left>
struct QuotientGenerator<Left, double> : not_constructible {
  using Type = Left;
};

template<template<typename> typename Quantity, typename Right>
struct QuotientGenerator<double, Quantity<Right>> : not_constructible {
  using Type = typename Collapse<Quantity<
      typename DimensionsQuotientGenerator<NoDimensions, Right>::Type>>::Type;
};

template<>
struct QuotientGenerator<double, double> : not_constructible {
  using Type = double;
};

}  // namespace internal_generators
}  // namespace quantities
}  // namespace principia
