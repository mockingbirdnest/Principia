#pragma once

#include <cstdint>

#include "quantities/cantor.hpp"

namespace principia {
namespace quantities {
namespace _dimensions {
namespace internal {

using namespace principia::quantities::_cantor;

// Dimensionality of physical quantities.  Note that we strongly type angles.
template<std::int64_t LengthExponent,
         std::int64_t MassExponent,
         std::int64_t TimeExponent,
         std::int64_t CurrentExponent,
         std::int64_t TemperatureExponent,
         std::int64_t AmountExponent,
         std::int64_t LuminousIntensityExponent,
         std::int64_t AngleExponent>
struct Dimensions;

// A double by any other name...
using NoDimensions = Dimensions<0, 0, 0, 0, 0, 0, 0, 0>;

// Instantiating this struct asserts at compile time that the template parameter
// can be serialized.
template<typename Dimensions>
struct DimensionsAreSerializable;

// These structs have a `Type` member that is a `Dimensions` suitable for
// the result of the operation applied to argument(s) having the `Dimensions`
// given as template parameter(s).

template<typename Dimensions, int n>
struct DimensionsExponentiationGenerator;

// Only legal if `n` divides the dimensions.
template<typename Dimensions, int n>
struct DimensionsNthRootGenerator;

template<typename LDimensions, typename RDimensions>
struct DimensionsProductGenerator;

template<typename LDimensions, typename RDimensions>
struct DimensionsQuotientGenerator;

// Metaprogramming.

template<typename T>
struct is_dimension;

template<typename T>
inline constexpr bool is_dimension_v = is_dimension<T>::value;

template<typename T>
concept dimension = is_dimension_v<T>;

template<typename Q>
concept dimensionful = requires {
  requires dimension<typename Q::Dimensions>;
};

template<typename Q>
concept dimensionless = countable<Q> || continuum<Q>;

}  // namespace internal

using internal::dimension;
using internal::dimensionful;
using internal::dimensionless;
using internal::Dimensions;
using internal::DimensionsAreSerializable;
using internal::DimensionsExponentiationGenerator;
using internal::DimensionsNthRootGenerator;
using internal::DimensionsProductGenerator;
using internal::DimensionsQuotientGenerator;
using internal::NoDimensions;

}  // namespace _dimensions
}  // namespace quantities
}  // namespace principia

#include "quantities/dimensions_body.hpp"
