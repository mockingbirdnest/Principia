#pragma once

#include "quantities/generators.hpp"

#include <tuple>

#include "base/not_constructible.hpp"
#include "boost/multiprecision/number.hpp"
#include "quantities/dimensions.hpp"

namespace principia {
namespace quantities {
namespace _generators {
namespace internal {

using namespace boost::multiprecision;
using namespace principia::base::_not_constructible;
using namespace principia::quantities::_dimensions;

template<typename T>
struct is_dimension : std::false_type {};

template<std::int64_t... exponents>
struct is_dimension<Dimensions<exponents...>> : std::true_type {};

template<typename T>
inline constexpr bool is_dimension_v = is_dimension<T>::value;

template<typename T>
concept dimension = is_dimension_v<T>;

template<typename Q>
concept dimensionful = requires {
  requires dimension<typename Q::Dimensions>;
};

template<typename Q>
concept dimensionless =
    std::floating_point<Q> || std::integral<Q> || is_number<Q>::value;

// The template template parameter `Quantity` on specializations lifts a
// circular dependency.

template<typename Q>
struct Collapse : not_constructible {
  using Type = Q;
};

template<template<typename> typename Quantity>
struct Collapse<Quantity<NoDimensions>> : not_constructible {
  using Type = double;
};

template<template<typename> typename Quantity, dimension D, int n>
struct ExponentiationGenerator<Quantity<D>, n> : not_constructible {
  using Type = typename Collapse<
      Quantity<typename DimensionsExponentiationGenerator<D, n>::Type>>::Type;
};

template<typename Q, int n>
  requires dimensionless<Q>
struct ExponentiationGenerator<Q, n> : not_constructible {
  using Type = Q;
};


template<template<typename> typename Quantity, dimension D, int n>
struct NthRootGenerator<Quantity<D>, n> : not_constructible {
  using Type = typename Collapse<
      Quantity<typename DimensionsNthRootGenerator<D, n>::Type>>::Type;
};

// NOTE(phl): This is designed so that we can write something like `Sqrt(2)`.
template<typename Q, int n>
  requires dimensionless<Q>
struct NthRootGenerator<Q, n> : not_constructible {
  using Type = double;
};


template<template<typename> typename Quantity, dimension Left, dimension Right>
struct ProductGenerator<Quantity<Left>, Quantity<Right>> : not_constructible {
  using Type = typename Collapse<Quantity<
      typename DimensionsProductGenerator<Left, Right>::Type>>::Type;
};

template<dimensionful Q1, dimensionless Q2>
struct ProductGenerator<Q1, Q2> : not_constructible {
  using Type = Q1;
};

template<dimensionless Q1, dimensionful Q2>
struct ProductGenerator<Q1, Q2> : not_constructible {
  using Type = Q2;
};

template<dimensionless Q>
struct ProductGenerator<Q, Q> : not_constructible {
  using Type = Q;
};


template<template<typename> typename Quantity, dimension Left, dimension Right>
struct QuotientGenerator<Quantity<Left>, Quantity<Right>> : not_constructible {
  using Type = typename Collapse<Quantity<
      typename DimensionsQuotientGenerator<Left, Right>::Type>>::Type;
};

template<dimensionful Q1, dimensionless Q2>
struct QuotientGenerator<Q1, Q2> : not_constructible {
  using Type = Q1;
};

template<dimensionless Q1,
         template<typename> typename Quantity,
         dimension Right>
struct QuotientGenerator<Q1, Quantity<Right>> : not_constructible {
  using Type = typename Collapse<Quantity<
      typename DimensionsQuotientGenerator<NoDimensions, Right>::Type>>::Type;
};

template<dimensionless Q>
struct QuotientGenerator<Q, Q> : not_constructible {
  using Type = Q;
};

}  // namespace internal
}  // namespace _generators
}  // namespace quantities
}  // namespace principia
